#include <algorithm>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// Constants and enumerations
enum class Side {TOP, BOTTOM, LEFT, RIGHT};
constexpr uint8_t UL_ABOVE = 8;
constexpr uint8_t UR_ABOVE = 4;
constexpr uint8_t LR_ABOVE = 2;
constexpr uint8_t LL_ABOVE = 1;

// Structs
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    Point start;
    Point end;
    Side start_side;
    Side end_side;
    bool visited;
} Segment;

typedef struct {
    std::unordered_multimap<double, Segment> segments;
    int row;
    int col;
} Square;

typedef struct {
    std::vector<Point> line_string;
    double level;
    bool is_closed;
} Contour;

// Globals
int nrows, ncols;
float *input_array;
Square *squares;


double interpolate(double low, double high, double level) {
    // Return the linear fraction of level between low and high
    assert(level > low && level < high);
    return (level - low) / (high - low);
};

void createSquare(Square &square, int row, int col, const std::vector<double>& levels) {
    assert(row < nrows - 1 && col < ncols - 1);
    square.row = row;
    square.col = col;

    // Get pixel values at corners around this square
    double val_ul = input_array[row * ncols + col];
    double val_ur = input_array[row * ncols + col + 1];
    double val_lr = input_array[(row + 1) * ncols + col + 1];
    double val_ll = input_array[(row + 1) * ncols + col];
    double pixel_min = std::min(val_ul, std::min(val_ur, std::min(val_lr, val_ll)));
    double pixel_max = std::max(val_ul, std::max(val_ur, std::max(val_lr, val_ll)));

    // Iterate over all levels, skipping those that are out of range for these pixels
    for (auto& level : levels) {
        if (level < pixel_min || level > pixel_max) {
            continue;
        }

        // TODO: figure out if this is right thing to do regarding tolerance?
        uint8_t key = 0;
        key |= (val_ul > level) ? UL_ABOVE : 0;
        key |= (val_ur > level) ? UR_ABOVE : 0;
        key |= (val_lr > level) ? LR_ABOVE : 0;
        key |= (val_ll > level) ? LL_ABOVE : 0;

        // We define the upper-left pixel to be the point at coordinates (row,col)
        // No need to add 0.5 to each, see https://gis.stackexchange.com/a/122687/4669
        double top = (double) row;
        double left = (double) col;
        double right = left + 1.0;
        double bottom = top + 1.0;

        Segment segment;
        segment.visited = false;
        switch(key) {
            // Lower-left pixel above, others below. Segment traverses left side to bottom side
            case 1:
                segment.start = {left, top + interpolate(val_ul, val_ll, level)};
                segment.end = {right - interpolate(val_lr, val_ll, level), bottom};
                segment.start_side = Side::LEFT;
                segment.end_side = Side::BOTTOM;
                square.segments.insert({level, segment});
                break;
            // Lower-right pixel above, others below. Segment traverses bottom side to right side
            case 2:
                segment.start = {left + interpolate(val_ll, val_lr, level), bottom};
                segment.end = {right, top + interpolate(val_ur, val_lr, level)};
                segment.start_side = Side::BOTTOM;
                segment.end_side = Side::RIGHT;
                square.segments.insert({level, segment});
                break;

            // Remember to disambiguate the saddles.

            default:
                break;
                // Raise an error.
        }
    }
}

// Find (unvisited) segment at given level starting from given side in given square
bool lookupSegmentInSquare(Segment **segment, double level, Square *square, Side start_side,
        bool unvisited_only) {
    int num_segments_found = 0;
    auto range = square->segments.equal_range(level);
    for (auto iter = range.first; iter != range.second; iter++) {
        if (iter->second.start_side == start_side &&
                (!unvisited_only || iter->second.visited == false)) {
            *segment = &(iter->second);
            num_segments_found++;
        }
    }
    assert(num_segments_found <= 1); // Should never have 2+ segments with same square/level/side
    return num_segments_found == 1;
}

// Go from segment at level in square to next segment, updating points in-place
bool traverseToNextSegment(Square **square, Segment **segment, double level) {

    // Traverse left
    if ((*segment)->end_side == Side::LEFT && (*square)->col > 0) {
        *square = &squares[(*square)->row * ncols + (*square)->col - 1];
        return lookupSegmentInSquare(segment, level, *square, Side::RIGHT, true);
    // Traverse right
    } else if ((*segment)->end_side == Side::RIGHT && (*square)->col < ncols - 1) {
        *square = &squares[(*square)->row * ncols + (*square)->col + 1];
        return lookupSegmentInSquare(segment, level, *square, Side::LEFT, true);
    // Traverse up
    } else if ((*segment)->end_side == Side::TOP && (*square)->row > 0) {
        *square = &squares[((*square)->row - 1) * ncols + (*square)->col];
        return lookupSegmentInSquare(segment, level, *square, Side::BOTTOM, true);
    // Traverse down
    } else if ((*segment)->end_side == Side::BOTTOM && (*square)->row < nrows - 1) {
        *square = &squares[((*square)->row + 1) * ncols + (*square)->col];
        return lookupSegmentInSquare(segment, level, *square, Side::TOP, true);
    // Hit the edge of the raster, can't go any farther!
    } else {
        return false;
    }

}

Contour* traverseContour(double level, Square *starting_square, Segment *starting_segment) {

    // Initialize the contour with the first two vertices
    Contour *contour = new Contour();
    contour->level = level;
    contour->line_string.push_back(starting_segment->start);
    contour->line_string.push_back(starting_segment->end);
    starting_segment->visited = true;

    // Traverse through all connected segments
    Square *current_square = starting_square;
    Segment *current_segment = starting_segment;
    while (traverseToNextSegment(&current_square, &current_segment, level)) {
        contour->line_string.push_back(current_segment->end);
        current_segment->visited = true;
    }

    return contour;
}

void traverseNonClosedContours(double level, std::vector<Contour*> *contours) {

    Segment *segment;
    int row = 0, col = 0;
    // Check squares in top row, exclusive of square in rightmost column
    for (; col < ncols - 1; col++) {
        Square *square = &squares[row * ncols + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::TOP, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in right column, exclusive of square in bottom row
    for (; row < nrows - 1; row++) {
        Square *square = &squares[row * ncols + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::RIGHT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in bottom row, exclusive of square in leftmost column
    for (; col > 0; col--) {
        Square *square = &squares[row * ncols + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::BOTTOM, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in left column, exclusive of square in top row
    for (; row > 0; row--) {
        Square *square = &squares[row * ncols + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::LEFT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }

}

void traverseClosedContours(double level, std::vector<Contour*> *contours) {

}

int main(int argc, char **argv) {

    // Initialize data structures that are proportional to the raster size
    scanf("%d %d", &nrows, &ncols);
    input_array = new float[nrows * ncols];
    squares = new Square[nrows * ncols];
    printf("rows: %d, cols: %d\n", nrows, ncols);

    // Populate the input data array
    float val;
    int index = 0;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            scanf("%f", &val);
            input_array[index] = val;
            index += 1;
        }
    }

    // Phase 1: generate segments in each square (in parallel)
    std::vector<double> levels = {0.0, 5.0, 10.0};
    for (int i = 0; i < nrows - 1; i++) {
        for (int j = 0; j < ncols - 1; j++) {
            createSquare(squares[i * ncols + j], i, j, levels);
        }
    }

    // Phase 2: join adjacent segments into contours
    std::vector<Contour *> contours;
    for (auto& level : levels) {
        traverseNonClosedContours(level, &contours);
        traverseClosedContours(level, &contours);
    }

    return 0;
}
