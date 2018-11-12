#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <omp.h>
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
        const coord_t left, const val_t ll, const val_t lr, const val_t ur, const val_t ul) {
    switch(side) {
        case Side::LEFT:
            return {left,  top + interpolate(ul, ll, level)};
        case Side::RIGHT:
            return {left + 1.0, top + interpolate(ur, lr, level)};
        case Side::TOP:
            return {left + interpolate(ul, ur, level), top};
        case Side::BOTTOM:
            return {left + interpolate(ll, lr, level), top + 1.0};
        default:
            assert(false);
    }
    assert(false);
}

struct Segment {
    Segment(const Side side_start, const Side side_end, const val_t level,
            const coord_t top, const coord_t left, const val_t ll,
            const val_t lr, const val_t ur, const val_t ul) :

        start_side(side_start), end_side(side_end), visited(false) {
        start = interpolatePoint(side_start, level, top, left, ll, lr, ur, ul);
        end = interpolatePoint(side_end, level, top, left, ll, lr, ur, ul);
    }

    Side start_side;
    Side end_side;
    bool visited;
    Point start;
    Point end;
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

typedef struct {
    int first_row;
    int first_col;
    int num_rows;
    int num_cols;
    // std::unordered_multimap<SegmentKey, Segment> segments; // For later
} Block;

using Time = std::chrono::high_resolution_clock;

// Globals from command line arguments
static int _argc;
static char **_argv;
const char *input_filename;
val_t interval = -1.0;  // Negative means subdivide into 10
int num_threads = 1;
int block_dim = 32;

// Globals from input file's header
int nrows, ncols;
Point raster_ll;
coord_t pixel_size;

// Globals from other places
FILE *input_file;
val_t *input_array;
Square *squares;
int nblocksh, nblocksv;
Block *blocks;
auto prev_time = Time::now();


const char *get_option_string(const char *option_name, const char *default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2) {
        if (strcmp(_argv[i], option_name) == 0) {
                return _argv[i + 1];
        }
    }
    return default_value;
}

int get_option_int(const char *option_name, int default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2) {
        if (strcmp(_argv[i], option_name) == 0) {
            return atoi(_argv[i + 1]);
        }
    }
    return default_value;
}

float get_option_float(const char *option_name, float default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2) {
        if (strcmp(_argv[i], option_name) == 0) {
            return (float)atof(_argv[i + 1]);
        }
    }
    return default_value;
}

// Populate square with top-left pixel at (row, col).
void populateSquare(Square *const square, const int row, const int col,
        const std::vector<val_t>& levels) {
    assert(row < nrows - 1 && col < ncols - 1);
    square->row = row;
    square->col = col;

    // Get pixel values at corners around this square
    val_t ll = input_array[(row + 1) * ncols + col];
    val_t lr = input_array[(row + 1) * ncols + col + 1];
    val_t ur = input_array[row * ncols + col + 1];
    val_t ul = input_array[row * ncols + col];
    val_t pixel_min = std::min(ll, std::min(lr, std::min(ur, ul)));
    val_t pixel_max = std::max(ll, std::max(lr, std::max(ur, ul)));

    // Iterate over all levels, skipping those that are out of range for these pixels
    for (const auto& level : levels) {
        if (level < pixel_min || level > pixel_max) {
            continue;
        }
        // TODO: figure out if this is right thing to do regarding tolerance?
        uint8_t key = 0;
        key |= (ll > level) ? LL_ABOVE : 0;
        key |= (lr > level) ? LR_ABOVE : 0;
        key |= (ur > level) ? UR_ABOVE : 0;
        key |= (ul > level) ? UL_ABOVE : 0;

        // We define the upper-left pixel to be the point at coordinates (row,col)
        // No need to add 0.5 to each, see https://gis.stackexchange.com/a/122687/4669.
        // TODO: ^^ are we sure about this?
        coord_t top = (coord_t) row;
        coord_t left = (coord_t) col;

        switch(key) {
            // Cases where 1 pixel is above the level
            case (LL_ABOVE): // Case 1 = 0001
                square->segments.insert({level, Segment(Side::LEFT, Side::BOTTOM,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (LR_ABOVE): // Case 2 = 0010
                square->segments.insert({level, Segment(Side::BOTTOM, Side::RIGHT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UR_ABOVE): // Case  4 = 0100
                square->segments.insert({level, Segment(Side::RIGHT, Side::TOP,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE): // Case  8 = 1000
                square->segments.insert({level, Segment(Side::TOP, Side::LEFT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            // Cases where 3 pixels are above the level
            case (UR_ABOVE | LR_ABOVE | LL_ABOVE): // Case  7 = 0111
                square->segments.insert({level, Segment(Side::LEFT, Side::TOP,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | LL_ABOVE | LR_ABOVE): // Case 11 = 1011
                square->segments.insert({level, Segment(Side::TOP, Side::RIGHT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LL_ABOVE): // Case 13 = 1101
                square->segments.insert({level, Segment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | UR_ABOVE | LR_ABOVE): // Case 14 = 1110
                square->segments.insert({level, Segment(Side::BOTTOM, Side::LEFT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            // Cases where 2 adjacent pixels are above the level
            case (LL_ABOVE | LR_ABOVE): // Case  3 = 0011
                square->segments.insert({level, Segment(Side::LEFT, Side::RIGHT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UR_ABOVE | LR_ABOVE): // Case  6 = 0110
                square->segments.insert({level, Segment(Side::BOTTOM, Side::TOP,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | LL_ABOVE): // Case  9 = 1001
                square->segments.insert({level, Segment(Side::TOP, Side::BOTTOM,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | UR_ABOVE): // Case 12 = 1100
                square->segments.insert({level, Segment(Side::RIGHT, Side::LEFT,
                        level, top, left, ll, lr, ur, ul)});
                break;
            // TODO(maybe): disambiguate saddle points.
            // Cases where 2 non-adjacent pixels are above the level, i.e. a saddle
            case (LL_ABOVE | UR_ABOVE): // Case  5 = 0101
                square->segments.insert({level, Segment(Side::LEFT, Side::TOP,
                        level, top, left, ll, lr, ur, ul)});
                square->segments.insert({level, Segment(Side::RIGHT, Side::BOTTOM,
                        level, top, left, ll, lr, ur, ul)});
                break;
            case (UL_ABOVE | LR_ABOVE): // Case 10 = 1010
                square->segments.insert({level, Segment(Side::TOP, Side::RIGHT,
                        level, top, left, ll, lr, ur, ul)});
                square->segments.insert({level, Segment(Side::BOTTOM, Side::LEFT,
                        level, top, left, ll, lr, ur, ul)});
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
    Segment *segment = NULL;
    int row = 0, col = 0;
    // Check squares in top row, exclusive of square in rightmost column
    for (; col < ncols - 2; col++) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::TOP, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            # pragma omp critical
            {
                contours->push_back(contour);
            }
        }
    }
    // Check squares in right column, exclusive of square in bottom row
    for (; row < nrows - 2; row++) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::RIGHT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            # pragma omp critical
            {
                contours->push_back(contour);
            }
        }
    }
    // Check squares in bottom row, exclusive of square in leftmost column
    for (; col > 0; col--) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::BOTTOM, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            # pragma omp critical
            {
                contours->push_back(contour);
            }
        }
    }
    // Check squares in left column, exclusive of square in top row
    for (; row > 0; row--) {
        Square *square = &squares[row * (ncols - 1) + col];
        if (lookupSegmentInSquare(&segment, level, square, Side::LEFT, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = false;
            # pragma omp critical
            {
                contours->push_back(contour);
            }
        }
    }
}

void traverseClosedContours(const val_t level, std::vector<Contour*> *const contours) {
    // Iterate over all squares, begin traversal at any unvisited segments
    for (int idx = 0; idx < (nrows - 1) * (ncols - 1); idx++) {
        Square *square = &squares[idx];
        Segment *segment = NULL;
        if (lookupSegmentInSquare(&segment, level, square, Side::ANY, true)) {
            Contour *contour = traverseContour(level, square, segment);
            contour->is_closed = true;
            # pragma omp critical
            {
                contours->push_back(contour);
            }
        }
    }
}

void countSegments(int *total_segments = nullptr, int *unvisited_segments = nullptr) {
    int total = 0, unvisited = 0;
    Square *square;
    for (int idx = 0; idx < (nrows - 1) * (ncols - 1); idx++) {
        square = &squares[idx];
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
    size_t i;
    for (i = 0; i < contours.size(); i++) {
        Contour *contour = contours[i];
        fprintf(output, "{ \"type\":\"Feature\", ");
        fprintf(output, "\"properties\": {\"level\":%.2f, \"is_closed\":%s}, ",
                contour->level, (contour->is_closed) ? "true" : "false");
        fprintf(output, "\"geometry\":{ \"type\":\"LineString\", \"coordinates\": [");
        size_t j;
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

void readArguments(int argc, char **argv) {

    // This directly copies the flag parsing pattern used by assignment 3.
    _argc = argc - 1;
    _argv = argv + 1;

    // Open input file
    input_filename = get_option_string("-f", nullptr);
    if (!(input_file = fopen(input_filename, "r"))) {
        printf("Unable to open file: %s.\n", input_filename);
        exit(-1);
    }

    // Read interval size
    interval = get_option_float("-i", interval);

    // Read number of threads and set the omp variable accordingly
    num_threads = get_option_int("-n", num_threads);
    omp_set_num_threads(num_threads);

    // Read block dimension
    block_dim = get_option_int("-b", block_dim);

}

// Read header of ASCII grid format (https://en.wikipedia.org/wiki/Esri_grid)
// Arbitary input raster files can be converted to this format with this GDAL command:
//    gdal_translate -of AAIGrid in.tif out.asc -co DECIMAL_PRECISION=3
void readHeader(FILE *input) {
    char str[100];
    int ret;
    ret = std::fscanf(input, "%s %d", str, &ncols);
    assert(ret > 0 && std::strcmp(str, "ncols") == 0);
    ret = std::fscanf(input, "%s %d", str, &nrows);
    assert(ret > 0 && std::strcmp(str, "nrows") == 0);
    ret = std::fscanf(input, "%s %lf", str, &raster_ll.x);
    assert(ret > 0 && std::strcmp(str, "xllcorner") == 0);
    ret = std::fscanf(input, "%s %lf", str, &raster_ll.y);
    assert(ret > 0 && std::strcmp(str, "yllcorner") == 0);
    ret = std::fscanf(input, "%s %lf", str, &pixel_size);
    assert(ret > 0 && std::strcmp(str, "cellsize") == 0);
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

void profileTime(std::string message) {
    auto new_time = Time::now();
    auto duration = new_time - prev_time;
    auto duration_microsecs = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    std::cout << message << " = " << duration_microsecs << " microseconds\n";
    prev_time = new_time;
}

int main(int argc, char **argv) {

    // Parse command line arguments
    readArguments(argc, argv);

    // Read entire header, which is in ASCII format, to global variables.
    readHeader(input_file);

    // Initialize data structures that are proportional to the raster size. Since squares are gaps
    // between pixels, Square should have (nrows-1) * (ncols-1) elements.
    input_array = new val_t[nrows * ncols];
    squares = new Square[(nrows - 1) * (ncols - 1)];

    // We will compute the min and max elevations in the input to determine what levels to generate
    // contours at.
    val_t val_min = std::numeric_limits<double>::max();
    val_t val_max = std::numeric_limits<double>::lowest();

    // Populate the input data array.
    int index = 0, ret;
    val_t val;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            ret = std::fscanf(input_file, "%lf", &val);
            assert(ret > 0);
            input_array[index] = val;
            index += 1;
            val_min = std::min(val_min, val);
            val_max = std::max(val_max, val);
        }
    }
    // Done with the input file.
    fclose(input_file);

    // Based on input args, min and max value of the input data, compute the levels.
    std::vector<val_t> levels = determineLevels(interval, val_min, val_max);

    // Perform spatial decomposition of squares into blocks, each of which will be handled by 1 thread
    nblocksh = (ncols - 1 + block_dim - 1) / block_dim;
    nblocksv = (nrows - 1 + block_dim - 1) / block_dim;
    printf("Input file = %s, num threads = %d, block size = %d\n",
            input_filename, num_threads, block_dim);
    printf("Dimensions = %d x %d pixels = %d x %d squares = %d x %d blocks, each %d x %d\n",
            nrows, ncols, nrows - 1, ncols - 1, nblocksv, nblocksh, block_dim, block_dim);
    blocks = new Block[nblocksh * nblocksv];
    for (int ii = 0; ii < nblocksv; ii++) {
        for (int jj = 0; jj < nblocksh; jj++) {
            Block *block = &blocks[ii * nblocksh + jj];
            block->first_row = ii * block_dim;
            block->first_col = jj * block_dim;
            block->num_rows = std::min(block_dim, nrows - 1 - block->first_row);
            block->num_cols = std::min(block_dim, ncols - 1 - block->first_col);
        }
    }
    profileTime("Initialization");

    // Phase 1: generate segments in each square (in parallel)
    int block_num;
    # pragma omp parallel default(shared) private(block_num)
    {
        # pragma omp for schedule(static) nowait
        for (block_num = 0; block_num < nblocksh * nblocksv; block_num++) {
            Block *block = &blocks[block_num];
            for (int i = block->first_row; i < block->first_row + block->num_rows; i++) {
                for (int j = block->first_col; j < block->first_col + block->num_cols; j++) {
                    populateSquare(&squares[i * (ncols - 1) + j], i, j, levels);
                }
            }
        }
    }
    profileTime("Phase 1");

    // Phase 2: join adjacent segments into contours
    std::vector<Contour *> contours;

    size_t i;
    # pragma omp parallel default (shared) private(i)
    {
        # pragma omp for schedule(static) nowait
        for (i = 0; i < levels.size(); i++) {
            const auto& level = levels[i];
            traverseNonClosedContours(level, &contours);
            traverseClosedContours(level, &contours);
        }
    }
    profileTime("Phase 2");


    char output_filename[100];
    sprintf(output_filename, "output_%s.txt", basename(input_filename));
    FILE *output_file = fopen(output_filename, "w");
    printGeoJSON(output_file, contours);
    fclose(output_file);

    return 0;
}
