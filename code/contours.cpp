#include <algorithm>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// Constants and enumerations
enum class Side {TOP, BOTTOM, LEFT, RIGHT, ANY};
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

void printContour(const Contour& contour) {
    printf("level: %.2f, closed: %d\n", contour.level, contour.is_closed);
    int i;
    printf("LINESTRING (");
    for (i = 0; i < contour.line_string.size() - 1; i++) {
        printf("%.2f %.2f, ", contour.line_string[i].x, contour.line_string[i].y);
    }
    printf("%.2f %.2f)\n", contour.line_string[i].x, contour.line_string[i].y);
}

inline double interpolate(double left_or_top, double right_or_bottom, double level) {
    // Return the coordinate at level linearly proportional between top/left and bottom/right
    // TODO(maybe): do we need to be careful about tolerance?
    double ret = (level - left_or_top) / (right_or_bottom - left_or_top);
    assert(ret >= 0);
    return ret;
};

inline Point interpolatePoint(Side side, double level, double top, double left,
        double val_ll, double val_lr, double val_ur, double val_ul) {
    switch(side) {
        case Side::LEFT:
            return {left,  top + interpolate(val_ul, val_ll, level)};
        case Side::RIGHT:
            return {left + 1.0, top + interpolate(val_ur, val_lr, level)};
        case Side::TOP:
            return {left + interpolate(val_ul, val_ur, level), top};
        case Side::BOTTOM:
            return {left + interpolate(val_ll, val_lr, level), top + 1.0};
    }
    assert(false);
}

inline Segment buildSegment(Side side_start, Side side_end, double level, double top, double left,
        double val_ll, double val_lr, double val_ur, double val_ul) {
    Segment segment;
    segment.visited = false;
    segment.start_side = side_start;
    segment.end_side = side_end;
    segment.start = interpolatePoint(side_start, level, top, left, val_ll, val_lr, val_ur, val_ul);
    segment.end = interpolatePoint(side_end, level, top, left, val_ll, val_lr, val_ur, val_ul);
    return segment;
}

void createSquare(Square *square, int row, int col, const std::vector<double>& levels) {
    assert(row < nrows - 1 && col < ncols - 1);
    square->row = row;
    square->col = col;

    // Get pixel values at corners around this square
    double val_ll = input_array[(row + 1) * ncols + col];
    double val_lr = input_array[(row + 1) * ncols + col + 1];
    double val_ur = input_array[row * ncols + col + 1];
    double val_ul = input_array[row * ncols + col];
    double pixel_min = std::min(val_ll, std::min(val_lr, std::min(val_ur, val_ul)));
    double pixel_max = std::max(val_ll, std::max(val_lr, std::max(val_ur, val_ul)));

    // Iterate over all levels, skipping those that are out of range for these pixels
    for (auto& level : levels) {
        if (level < pixel_min || level > pixel_max) {
            continue;
        }

        // TODO: figure out if this is right thing to do regarding tolerance?
        uint8_t key = 0;
        key |= (val_ll > level) ? LL_ABOVE : 0;
        key |= (val_lr > level) ? LR_ABOVE : 0;
        key |= (val_ur > level) ? UR_ABOVE : 0;
        key |= (val_ul > level) ? UL_ABOVE : 0;

        // We define the upper-left pixel to be the point at coordinates (row,col)
        // No need to add 0.5 to each, see https://gis.stackexchange.com/a/122687/4669
        double top = (double) row;
        double left = (double) col;

        switch(key) {
            // Cases where 1 pixel is above the level
            case (LL_ABOVE): // Case 1 = 0001
                square->segments.insert({level, buildSegment(Side::LEFT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (LR_ABOVE): // Case 2 = 0010
                square->segments.insert({level, buildSegment(Side::BOTTOM, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UR_ABOVE): // Case  4 = 0100
                square->segments.insert({level, buildSegment(Side::RIGHT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE): // Case  8 = 1000
                square->segments.insert({level, buildSegment(Side::TOP, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // Cases where 3 pixels are above the level
            case (UR_ABOVE | LR_ABOVE | LL_ABOVE): // Case  7 = 0111
                square->segments.insert({level, buildSegment(Side::LEFT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LL_ABOVE | LR_ABOVE): // Case 11 = 1011
                square->segments.insert({level, buildSegment(Side::TOP, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LL_ABOVE): // Case 13 = 1101
                square->segments.insert({level, buildSegment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LR_ABOVE): // Case 14 = 1110
                square->segments.insert({level, buildSegment(Side::BOTTOM, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // Cases where 2 adjacent pixels are above the level
            case (LL_ABOVE | LR_ABOVE): // Case  3 = 0011
                square->segments.insert({level, buildSegment(Side::LEFT, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UR_ABOVE | LR_ABOVE): // Case  6 = 0110
                square->segments.insert({level, buildSegment(Side::BOTTOM, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LL_ABOVE): // Case  9 = 1001
                square->segments.insert({level, buildSegment(Side::TOP, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE): // Case 12 = 1100
                square->segments.insert({level, buildSegment(Side::RIGHT, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // TODO(maybe): disambiguate saddle points.
            // Cases where 2 non-adjacent pixels are above the level, i.e. a saddle
            case (LL_ABOVE | UR_ABOVE): // Case  5 = 0101
                square->segments.insert({level, buildSegment(Side::LEFT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                square->segments.insert({level, buildSegment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LR_ABOVE): // Case 10 = 1010
                square->segments.insert({level, buildSegment(Side::TOP, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                square->segments.insert({level, buildSegment(Side::BOTTOM, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // Should never reach here
            default:
                assert(false);
        }
    }
}

// Find (unvisited) segment at given level starting from given side in given square
bool lookupSegmentInSquare(Segment **segment, double level, Square *square, Side start_side,
        bool unvisited_only) {
    int num_segments_found = 0;
    auto range = square->segments.equal_range(level);
    for (auto iter = range.first; iter != range.second; iter++) {
        if ((start_side == Side::ANY || iter->second.start_side == start_side) &&
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
    // Iterate over all squares, begin traversal at any unvisited segments
    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            Square *square = &squares[row * ncols + col];
            Segment *segment;
            if (lookupSegmentInSquare(&segment, level, square, Side::ANY, true)) {
                Contour *contour = traverseContour(level, square, segment);
                contour->is_closed = false;
                contours->push_back(contour);
            }
        }
    }
}

void countSegments(int *total_segments = nullptr, int *unvisited_segments = nullptr) {
    int total = 0, unvisited = 0;
    Square *square;
    Segment *segment;
    for (int row = 0; row < nrows; row++) {
        for (int col = 0; col < ncols; col++) {
            Square *square = &squares[row * ncols + col];
            for (auto iter = square->segments.begin(); iter != square->segments.end(); iter++) {
                total++;
                if (iter->second.visited == false) {
                    unvisited++;
                }
            }
        }
    }
    if (total_segments != nullptr && unvisited_segments != nullptr) {
        *total_segments = total;
        *unvisited_segments = unvisited;
    } else {
        float percent = ((float) total) / ((float) unvisited);
        printf("Found %i segments, %i (%.2f%%) are unvisited.\n", total, unvisited, percent);
    }
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
    std::vector<double> levels = {180.0, 185.0, 190.0, 195.0};
    for (int i = 0; i < nrows - 1; i++) {
        for (int j = 0; j < ncols - 1; j++) {
            createSquare(&squares[i * ncols + j], i, j, levels);
        }
    }
    countSegments();

    // Phase 2: join adjacent segments into contours
    std::vector<Contour *> contours;
    for (auto& level : levels) {
        traverseNonClosedContours(level, &contours);
        traverseClosedContours(level, &contours);
    }
    countSegments();
    printf("%d contours found\n", (int) contours.size());
    for (const auto& contour : contours) {
        printContour(*contour);
    }

    return 0;
}
