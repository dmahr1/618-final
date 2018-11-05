#include <algorithm>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

int nrows, ncols;
float *input_array;

typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    int end_side_id;
    Point end;
    bool visited;

    // TODO: Consider how necesary this is in the future.
    int start_side_id;
    Point start;
} Segment;

typedef struct {
    int row;
    int col;
    std::unordered_map<double, std::vector<Segment>> segments;
} Square;

typedef struct {
    double level;
    bool is_closed;
    std::vector<Point> line_string;
} Contour;

constexpr uint8_t NW_ABOVE = 8;
constexpr uint8_t NE_ABOVE = 4;
constexpr uint8_t SE_ABOVE = 2;
constexpr uint8_t SW_ABOVE = 1;

enum class Sides {TOP, BOTTOM, LEFT, RIGHT};

double interpolate(double low, double high, double level) {
    return 0.0;
};

// TODO: implement this.
int computeSideId(const int row, const int col, Sides side) {
    return 0;
}

void processSquare(const int row, const int col, const std::vector<double>& levels,
        std::vector<std::unordered_map<int, Segment>> *segment_maps,
        std::vector<std::vector<Segment>> *boundary_segments) {
    assert(row < nrows - 1);
    assert(col < ncols - 1);

    // nw = northwest, ...
    double nw, ne, se, sw;

    nw = input_array[row * ncols + col];
    ne = input_array[row * ncols + col + 1];
    se = input_array[(row + 1) * ncols + col + 1];
    sw = input_array[(row + 1) * ncols + col];

    double pixel_min, pixel_max;
    pixel_min = std::min(nw, std::min(ne, std::min(se, sw)));
    pixel_max = std::max(nw, std::max(ne, std::max(se, sw)));

    Square square = {row, col, {}};

    for (int level_index = 0; level_index < levels.size(); level_index++) {
        double level = levels[level_index];
        // Skip levels that are totally below or totally above this pixel.
        if (level < pixel_min || level > pixel_max) {
            continue;
        }

        uint8_t key = 0;

        // TODO: figure out if this is right thing to do regarding tolerance?
        if (nw > level) {
            key |= NW_ABOVE;
        }
        if (ne > level) {
            key |= NE_ABOVE;
        }
        if (se > level) {
            key |= SE_ABOVE;
        }
        if (sw > level) {
            key |= SW_ABOVE;
        }

        Segment segment0, segment1;
        segment0.visited = false;
        segment1.visited = false;
        double row_f = (double) row;
        double col_f = (double) col;
        // TODO(long-term): Think about how to refactor / clean up this.
        switch(key) {
            // Need to add 0.5 since the points we interpolate between are actually centers of squares.
            case 1:
                segment0.start = {col_f + 0.5, row_f + 0.5 + interpolate(nw, sw, level)};
                segment0.start_side_id = computeSideId(row, col, Sides::LEFT);
                // If left side of this square is a border, then this is a start segment.
                if (col == 0) {
                    (*boundary_segments)[level_index].push_back(segment0);
                }

                segment0.end = {col_f + 0.5 + interpolate(se, sw, level), row_f + 1 + 0.5};
                segment0.end_side_id = computeSideId(row, col, Sides::BOTTOM);
                (*segment_maps)[level_index].insert({segment0.start_side_id, segment0});
                break;
            case 2:
                segment0.start = {col_f + interpolate(se, sw, level), row_f + 1};
                segment0.end = {col_f + 1, row_f + interpolate(ne, se, level)};
                break;

            // Remember to disambiguate the saddles.

            default:
                break;
                // Raise an error.
        }
    }
}

// TODO: check that compiler will do return value optimization, or else pass in
// the vector so that it may be populated in place.
std::vector<Point> traverseSegments(std::unordered_map<int, Segment>* segment_map,
        Segment start_segment) {
    Segment current_segment = start_segment;
    std::vector<Point> line_string;

    line_string.push_back(current_segment.start);
    line_string.push_back(current_segment.end);

    // Don't think we can directly modify start_segment since it's a different copy.
    // So we have to this weird thing where use a segment to look up "itself".
    // Idea is that segment_map should be the single source of truth of whether
    // a segment has been visited or not.
    segment_map->at(current_segment.start_side_id).visited = true;

    while(segment_map->count(current_segment.end_side_id) > 0 &&
            segment_map->at(current_segment.end_side_id).visited == false) {
        // Go to the next segment (the one that starts from the current end side).
        current_segment = segment_map->at(current_segment.end_side_id);
        current_segment.visited = true;
        line_string.push_back(current_segment.end);
    }
    return line_string;
}

void traverseNonClosedContours(std::unordered_map<int, Segment>* segment_map, 
        std::vector<Segment> start_segments, double level, std::vector<Contour> *contours) {
    for (auto& start_segment : start_segments) {
        Contour contour;
        contour.line_string = traverseSegments(segment_map, start_segment);
        contour.is_closed = false;
        contour.level = level;
        contours->push_back(contour);
    }
    return;
}

void traverseClosedContours(Square *squares, double level,
    std::vector<Contour> *contours) {
    return;
}

int main(int argc, char **argv) {
    scanf("%d %d", &nrows, &ncols);

    input_array = (float *) malloc(nrows * ncols * sizeof(float));

    printf("rows: %d, cols: %d\n", nrows, ncols);

    float val;
    int index = 0;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            scanf("%f", &val);
            input_array[index] = val;
            index += 1;
        }
    }

    std::vector<double> levels = {0.0, 5.0, 10.0};
    // Global segment map. The ith map is for isovalue levels[i].
    std::vector<std::unordered_map<int, Segment>> segment_maps(levels.size());
    // We can tell at processing time whether a segment should be the starting point for a traversal
    // of a non-closed segment. Save them here.
    std::vector<std::vector<Segment>> boundary_segments;
    for (int i = 0; i < nrows - 1; i++) {
        for (int j = 0; j < ncols - 1; j++) {
            processSquare(i, j, levels, &segment_maps, &boundary_segments);
        }
    }

    std::vector<Contour> contours;
    for (int i = 0; i < levels.size(); i++) {
        traverseNonClosedContours(&segment_maps[i], boundary_segments[i], levels[i], &contours);
        //traverseClosedContours(squares, level, &contours);
    }


    return 0;
}
