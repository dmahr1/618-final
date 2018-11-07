#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
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

// Structs and typedefs
typedef double val_t;
typedef double coord_t;

typedef struct {
    coord_t x;
    coord_t y;
} Point;

inline val_t interpolate(const val_t left_or_top, const val_t right_or_bottom, const val_t level) {
    // Return the coordinate at level linearly proportional between top/left and bottom/right
    // TODO(maybe): do we need to be careful about tolerance?
    val_t ret = (level - left_or_top) / (right_or_bottom - left_or_top);
    assert(ret >= 0);
    return ret;
};

inline Point interpolatePoint(const Side side, const val_t level, const coord_t top, 
        const coord_t left, const val_t val_ll, const val_t val_lr, const val_t val_ur, 
        const val_t val_ul) {
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

struct Segment {
    Segment(const Side side_start, const Side side_end, const val_t level, 
        const coord_t top, const coord_t left, const val_t val_ll, 
        const val_t val_lr, const val_t val_ur, const val_t val_ul) :
    
        visited(false), start_side(side_start), end_side(side_end) {
        start = interpolatePoint(side_start, level, top, left, val_ll, val_lr, val_ur, val_ul);
        end = interpolatePoint(side_end, level, top, left, val_ll, val_lr, val_ur, val_ul);
    }  

    Point start;
    Point end;
    Side start_side;
    Side end_side;
    bool visited;
};

typedef struct Segment Segment;

typedef struct {
    std::unordered_multimap<val_t, Segment> segments;
    int row;
    int col;
} Square;

typedef struct {
    std::vector<Point> line_string;
    val_t level;
    bool is_closed;
} Contour;

using Time = std::chrono::high_resolution_clock;

// Globals
int nrows, ncols;
val_t *input_array;
Square *squares;
Point raster_ll;
coord_t pixel_size;

static int _argc;
static char **_argv;

const char *get_option_string(const char *option_name,
                                      const char *default_value) {
        for (int i = _argc - 2; i >= 0; i -= 2)
                    if (strcmp(_argv[i], option_name) == 0)
                                    return _argv[i + 1];
            return default_value;
}

int get_option_int(const char *option_name, int default_value) {
        for (int i = _argc - 2; i >= 0; i -= 2)
                    if (strcmp(_argv[i], option_name) == 0)
                                    return atoi(_argv[i + 1]);
            return default_value;
}

float get_option_float(const char *option_name, float default_value) {
        for (int i = _argc - 2; i >= 0; i -= 2)
                    if (strcmp(_argv[i], option_name) == 0)
                                    return (float)atof(_argv[i + 1]);
            return default_value;
}

// Populate square with top-left pixel at (row, col).
void populateSquare(Square *const square, const int row, const int col, 
        const std::vector<val_t>& levels) {
    assert(row < nrows - 1 && col < ncols - 1);
    square->row = row;
    square->col = col;

    // Get pixel values at corners around this square
    val_t val_ll = input_array[(row + 1) * ncols + col];
    val_t val_lr = input_array[(row + 1) * ncols + col + 1];
    val_t val_ur = input_array[row * ncols + col + 1];
    val_t val_ul = input_array[row * ncols + col];
    val_t pixel_min = std::min(val_ll, std::min(val_lr, std::min(val_ur, val_ul)));
    val_t pixel_max = std::max(val_ll, std::max(val_lr, std::max(val_ur, val_ul)));

    // Iterate over all levels, skipping those that are out of range for these pixels
    for (const auto& level : levels) {
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
        // No need to add 0.5 to each, see https://gis.stackexchange.com/a/122687/4669.
        // TODO: ^^ are we sure about this?
        coord_t top = (coord_t) row;
        coord_t left = (coord_t) col;

        switch(key) {
            // Cases where 1 pixel is above the level
            case (LL_ABOVE): // Case 1 = 0001
                square->segments.insert({level, Segment(Side::LEFT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (LR_ABOVE): // Case 2 = 0010
                square->segments.insert({level, Segment(Side::BOTTOM, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UR_ABOVE): // Case  4 = 0100
                square->segments.insert({level, Segment(Side::RIGHT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE): // Case  8 = 1000
                square->segments.insert({level, Segment(Side::TOP, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // Cases where 3 pixels are above the level
            case (UR_ABOVE | LR_ABOVE | LL_ABOVE): // Case  7 = 0111
                square->segments.insert({level, Segment(Side::LEFT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LL_ABOVE | LR_ABOVE): // Case 11 = 1011
                square->segments.insert({level, Segment(Side::TOP, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LL_ABOVE): // Case 13 = 1101
                square->segments.insert({level, Segment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LR_ABOVE): // Case 14 = 1110
                square->segments.insert({level, Segment(Side::BOTTOM, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // Cases where 2 adjacent pixels are above the level
            case (LL_ABOVE | LR_ABOVE): // Case  3 = 0011
                square->segments.insert({level, Segment(Side::LEFT, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UR_ABOVE | LR_ABOVE): // Case  6 = 0110
                square->segments.insert({level, Segment(Side::BOTTOM, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LL_ABOVE): // Case  9 = 1001
                square->segments.insert({level, Segment(Side::TOP, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | UR_ABOVE): // Case 12 = 1100
                square->segments.insert({level, Segment(Side::RIGHT, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            // TODO(maybe): disambiguate saddle points.
            // Cases where 2 non-adjacent pixels are above the level, i.e. a saddle
            case (LL_ABOVE | UR_ABOVE): // Case  5 = 0101
                square->segments.insert({level, Segment(Side::LEFT, Side::TOP,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                square->segments.insert({level, Segment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
            case (UL_ABOVE | LR_ABOVE): // Case 10 = 1010
                square->segments.insert({level, Segment(Side::TOP, Side::RIGHT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                square->segments.insert({level, Segment(Side::BOTTOM, Side::LEFT,
                        level, top, left, val_ll, val_lr, val_ur, val_ul)});
                break;
        }
    }
}

// Find (unvisited) segment at given level starting from given side in given square
bool lookupSegmentInSquare(Segment **segment, val_t level, Square *square, Side start_side,
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
bool traverseToNextSegment(Square **square, Segment **segment, val_t level) {

    // Traverse left
    if ((*segment)->end_side == Side::LEFT && (*square)->col > 0) {
        *square = *square - 1;
        return lookupSegmentInSquare(segment, level, *square, Side::RIGHT, true);
    // Traverse right
    } else if ((*segment)->end_side == Side::RIGHT && (*square)->col < ncols - 2) {
        *square = *square + 1;
        return lookupSegmentInSquare(segment, level, *square, Side::LEFT, true);
    // Traverse up
    } else if ((*segment)->end_side == Side::TOP && (*square)->row > 0) {
        *square = *square - (ncols - 1);
        return lookupSegmentInSquare(segment, level, *square, Side::BOTTOM, true);
    // Traverse down
    } else if ((*segment)->end_side == Side::BOTTOM && (*square)->row < nrows - 2) {
        *square = *square + (ncols - 1);
        return lookupSegmentInSquare(segment, level, *square, Side::TOP, true);
    // Hit the edge of the raster, can't go any farther!
    } else {
        return false;
    }

}

Contour* traverseContour(val_t level, Square *starting_square, Segment *starting_segment) {

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

void traverseNonClosedContours(const val_t level, std::vector<Contour*> *const contours) {
    Segment *segment;
    int row = 0, col = 0;
    // Check squares in top row, exclusive of square in rightmost column
    for (; col < ncols - 2; col++) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::TOP, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in right column, exclusive of square in bottom row
    for (; row < nrows - 2; row++) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::RIGHT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in bottom row, exclusive of square in leftmost column
    for (; col > 0; col--) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::BOTTOM, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
    // Check squares in left column, exclusive of square in top row
    for (; row > 0; row--) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::LEFT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            contours->push_back(contour);
        }
    }
}

void traverseClosedContours(const val_t level, std::vector<Contour*> *const contours) {
    // Iterate over all squares, begin traversal at any unvisited segments
    for (int idx = 0; idx < (nrows - 1) * (ncols - 1); idx++) {
        Square *square = &squares[idx];
        Segment *segment;
        if (lookupSegmentInSquare(&segment, level, square, Side::ANY, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = true;
            contours->push_back(contour);
        }
    }
}

void countSegments(int *total_segments = nullptr, int *unvisited_segments = nullptr) {
    int total = 0, unvisited = 0;
    Square *square;
    Segment *segment;
    for (int idx = 0; idx < (nrows - 1) * (ncols - 1); idx++) {
        Square *square = &squares[idx];
        for (auto iter = square->segments.begin(); iter != square->segments.end(); iter++) {
            total++;
            if (iter->second.visited == false) {
                unvisited++;
            }
        }
    }
    // If int variables were passed in, write to them.
    if (total_segments != nullptr && unvisited_segments != nullptr) {
        *total_segments = total;
        *unvisited_segments = unvisited;
    // Else just print the results
    } else {
        float percent = ((float) total) / ((float) unvisited);
        printf("Found %i segments, %i (%.2f%%) are unvisited.\n", total, unvisited, percent);
    }
}

inline Point transformPoint(Point point) {
    point.x =  point.x * pixel_size + raster_ll.x;
    point.y = (-point.y + nrows) * pixel_size + raster_ll.y;
    return point;
}

void printGeoJSON(FILE *output, const std::vector<Contour *> contours) {
    fprintf(output,  "{ \"type\":\"FeatureCollection\", \"features\": [\n");
    int i;
    for (i = 0; i < contours.size(); i++) {
        Contour *contour = contours[i];
        fprintf(output, "{ \"type\":\"Feature\", ");
        fprintf(output, "\"properties\": {\"level\":%.2f, \"is_closed\":%s}, ",
                contour->level, (contour->is_closed) ? "true" : "false");
        fprintf(output, "\"geometry\":{ \"type\":\"LineString\", \"coordinates\": [");
        int j;
        for (j = 0; j < contour->line_string.size() - 1; j++) {
            Point pt = transformPoint(contour->line_string[j]);
            fprintf(output, "[%.8lf,%.8lf],", pt.x, pt.y);
        }
        Point pt = transformPoint(contour->line_string[j]);
        fprintf(output, "[%.8lf,%.8lf] ] } }", pt.x, pt.y);
        fprintf(output, "%s\n", (i < contours.size() - 1) ? "," : "");
    }
    fprintf(output, "] } \n");
}

// Read header of ASCII grid format (https://en.wikipedia.org/wiki/Esri_grid)
// Arbitary input raster files can be converted to this format with this GDAL command:
//    gdal_translate -of AAIGrid in.tif out.asc -co DECIMAL_PRECISION=3
void readHeader(FILE *input) {
    char str[100];
    std::fscanf(input, "%s %d", str, &ncols);
    assert(std::strcmp(str, "ncols") == 0);
    std::fscanf(input, "%s %d", str, &nrows);
    assert(std::strcmp(str, "nrows") == 0);
    std::fscanf(input, "%s %lf", str, &raster_ll.x);
    assert(std::strcmp(str, "xllcorner") == 0);
    std::fscanf(input, "%s %lf", str, &raster_ll.y);
    assert(std::strcmp(str, "yllcorner") == 0);
    std::fscanf(input, "%s %lf", str, &pixel_size);
    assert(std::strcmp(str, "cellsize") == 0);
}

std::vector<val_t> determineLevels(val_t interval, val_t val_min, val_t val_max) {
    // By default, generate contours for 10 levels, equally spaced between val_min and val_max.
    if (interval < 0) {
        interval = (val_max - val_min) / 10.0;
    }

    // Populate levels every <interval> starting just below val_min to just above val_max
    std::vector<val_t> levels;
    val_t val;
    for (val = val_min - fmod(val_min, interval); val < val_max; val += interval) {
        levels.push_back(val);
    }
    levels.push_back(val);

    return levels;

}

int main(int argc, char **argv) {
    // This directly copies the flag parsing pattern used by assignment 3.
    _argc = argc - 1;
    _argv = argv + 1;
       
    const char *input_filename = get_option_string("-f", nullptr);
    FILE *input = fopen(input_filename, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", input_filename);
        return -1;
    }

    val_t interval = get_option_float("-i", -1.0);

    auto init_start = Time::now();

    // Read entire header, which is in ASCII format, to global variables.
    readHeader(input);

    // Initialize data structures that are proportional to the raster size. Since squares are gaps
    // between pixels, Square should have (nrows-1) * (ncols-1) elements.
    input_array = new val_t[nrows * ncols];
    squares = new Square[(nrows - 1) * (ncols - 1)];

    // We will compute the min and max elevations in the input to determine what levels to generate
    // contours at.
    val_t val_min = std::numeric_limits<double>::max();
    val_t val_max = std::numeric_limits<double>::lowest();

    // Populate the input data array.
    int index = 0;
    val_t val;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            std::fscanf(input, "%lf", &val);
            input_array[index] = val;
            index += 1;
            val_min = std::min(val_min, val);
            val_max = std::max(val_max, val);
        }
    }
    // Done with the input file.
    fclose(input);

    // Based on input args, min and max value of the input data, compute the levels.
    std::vector<val_t> levels = determineLevels(interval, val_min, val_max);

    std::cout << "Initialization time: " << std::chrono::duration_cast<std::chrono::microseconds>(Time::now() - init_start).count() 
        << " microseconds\n";

    auto p1_start = Time::now();
    // Phase 1: generate segments in each square (in parallel)
    for (int i = 0; i < nrows - 1; i++) {
        for (int j = 0; j < ncols - 1; j++) {
            populateSquare(&squares[i * (ncols - 1) + j], i, j, levels);
        }
    }
    std::cout << "Phase 1 time: " << std::chrono::duration_cast<std::chrono::microseconds>(Time::now() - p1_start).count() << " microseconds\n";

    // Phase 2: join adjacent segments into contours
    std::vector<Contour *> contours;

    for (const auto& level : levels) {
        auto p2_start = Time::now();
        traverseNonClosedContours(level, &contours);
        std::cout << "Phase 2.1, level "<< level << " time: " << std::chrono::duration_cast<std::chrono::microseconds>(Time::now() - p2_start).count() << " microseconds\n";
        p2_start = Time::now();
        traverseClosedContours(level, &contours);
        std::cout << "Phase 2.2, level "<< level << " time: " << std::chrono::duration_cast<std::chrono::microseconds>(Time::now() - p2_start).count() << " microseconds\n";
    }


    char output_filename[100];
    sprintf(output_filename, "output_%s.txt", basename(input_filename));
    FILE *output_file = fopen(output_filename, "w");
    printGeoJSON(output_file, contours);
    fclose(output_file);

    return 0;
}
