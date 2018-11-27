#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

using Time = std::chrono::high_resolution_clock;

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

// Declaring these two global variables here, as they are used by
// computeSideIndex.
int nrows, ncols;

int computeSideIndex(int square_row, int square_col, Side side) {
    // Assume horizontal sides are numbered first, then vertical.
    switch(side) {
        case (Side::TOP):
            return square_row * (ncols - 1) + square_col;
        case (Side::BOTTOM):
            return (square_row + 1) * (ncols - 1) + square_col;
        case (Side::LEFT):
            return nrows * (ncols - 1) + square_row * ncols + square_col;
        case (Side::RIGHT):
            return nrows * (ncols - 1) + square_row * ncols + square_col + 1;
        case (Side::ANY):
            printf("side argument to computeSideIndex cannot be Side::ANY\n");
            exit(-1);
        default:
            printf("invalid Side enum\n");
            exit(-1);
    }
}

struct Segment {
    Segment(const Side side_start, const Side side_end, 
            const int side_start_idx, const int end_side_index, const val_t level,
            const coord_t top, const coord_t left, const val_t ll, 
            const val_t lr, const val_t ur, const val_t ul) :
        side_start_idx(side_start_idx), end_side_index(end_side_index),
        visited(false) {
        start = interpolatePoint(side_start, level, top, left, ll, lr, ur, ul);
        end = interpolatePoint(side_end, level, top, left, ll, lr, ur, ul);
    }

    int side_start_idx;
    int end_side_index;
    bool visited;
    Point start;
    Point end;
};

typedef struct Segment Segment;

struct Contour {
    Contour(const std::shared_ptr<std::vector<Point>>& line_string,
            const val_t& level, int start_side_idx, int end_side_idx) :
        line_string(line_string), level(level), start_side_idx(start_side_idx),
        end_side_idx(end_side_idx), visited(false), is_closed(false) {}

    std::shared_ptr<std::vector<Point>> line_string;
    val_t level;
    int start_side_idx;
    int end_side_idx;
    bool visited;
    // TODO: Consider whether this is still necessary?
    bool is_closed;
};

typedef struct Contour Contour;

struct pair_hash {
    // Using pair hash from https://stackoverflow.com/a/32685618/1254992.
    std::size_t operator () (const std::pair<val_t, int> &p) const {
        auto h1 = std::hash<val_t>{}(p.first);
        auto h2 = std::hash<int>{}(p.second);
        return h1 ^ h2;
    }
};

typedef std::pair<val_t, int> SegmentKey;

typedef struct {
    int first_row;
    int first_col;
    int num_rows;
    int num_cols;
    std::unordered_multimap<SegmentKey, Segment, pair_hash> interior_segments;
    std::unordered_multimap<SegmentKey, Segment, pair_hash> inbound_segments;
    std::unordered_multimap<SegmentKey, Contour, pair_hash> interior_contours;
    std::unordered_multimap<SegmentKey, Contour, pair_hash> inbound_contours;
} Block;

// Globals from command line arguments
static int _argc;
static char **_argv;
const char *input_filename;
val_t interval = -1.0;  // Negative means subdivide into 10
int num_threads = 1;
int block_dim = 32;

// Globals from input file's header
Point raster_ll;
coord_t pixel_size;

// Globals from other places
FILE *input_file;
val_t *input_array;
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

bool segmentIsInboundToBlock(const Block *const block, const int row,
        const int col, const Side side) {
    switch (side) {
        case (Side::TOP):
            return row == block->first_row;
        case (Side::BOTTOM):
            return row == block->first_row + block->num_rows - 1;
        case (Side::LEFT):
            return col == block->first_col;
        case (Side::RIGHT):
            return col == block->first_col + block->num_rows - 1;
        case (Side::ANY):
            printf("Undefined behavior: checking whether segment is inbound "
                   "without specifying starting side (passed Side::ANY)\n");
            exit(-1);
        default:
            printf("Undefined behavior: checking whether segment is inbound "
                   "with invalid Side\n");
            exit(-1);
    }
}

// Populate square with top-left pixel at (row, col).
void processSquare(Block *const block, const int row, const int col,
        const std::vector<val_t>& levels) {
    assert(row < nrows - 1 && col < ncols - 1);

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


        std::vector<std::pair<Side, Side>> start_end_sides;

        switch(key) {
            // Cases where 1 pixel is above the level
            case (LL_ABOVE): // Case 1 = 0001
                start_end_sides.push_back({Side::LEFT, Side::BOTTOM});
                break;
            case (LR_ABOVE): // Case 2 = 0010
                start_end_sides.push_back({Side::BOTTOM, Side::RIGHT});
                break;
            case (UR_ABOVE): // Case  4 = 0100
                start_end_sides.push_back({Side::RIGHT, Side::TOP});
                break;
            case (UL_ABOVE): // Case  8 = 1000
                start_end_sides.push_back({Side::TOP, Side::LEFT});
                break;
            // Cases where 3 pixels are above the level
            case (UR_ABOVE | LR_ABOVE | LL_ABOVE): // Case  7 = 0111
                start_end_sides.push_back({Side::LEFT, Side::TOP});
                break;
            case (UL_ABOVE | LL_ABOVE | LR_ABOVE): // Case 11 = 1011
                start_end_sides.push_back({Side::TOP, Side::RIGHT});
                break;
            case (UL_ABOVE | UR_ABOVE | LL_ABOVE): // Case 13 = 1101
                start_end_sides.push_back({Side::RIGHT, Side::BOTTOM});
                break;
            case (UL_ABOVE | UR_ABOVE | LR_ABOVE): // Case 14 = 1110
                start_end_sides.push_back({Side::BOTTOM, Side::LEFT});
                break;
            // Cases where 2 adjacent pixels are above the level
            case (LL_ABOVE | LR_ABOVE): // Case  3 = 0011
                start_end_sides.push_back({Side::LEFT, Side::RIGHT});
                break;
            case (UR_ABOVE | LR_ABOVE): // Case  6 = 0110
                start_end_sides.push_back({Side::BOTTOM, Side::TOP});
                break;
            case (UL_ABOVE | LL_ABOVE): // Case  9 = 1001
                start_end_sides.push_back({Side::TOP, Side::BOTTOM});
                break;
            case (UL_ABOVE | UR_ABOVE): // Case 12 = 1100
                start_end_sides.push_back({Side::RIGHT, Side::LEFT});
                break;
            // TODO(maybe): disambiguate saddle points.
            // Cases where 2 non-adjacent pixels are above the level, i.e. a saddle
            case (LL_ABOVE | UR_ABOVE): // Case  5 = 0101
                start_end_sides.push_back({Side::LEFT, Side::TOP});
                start_end_sides.push_back({Side::RIGHT, Side::BOTTOM});
                break;
            case (UL_ABOVE | LR_ABOVE): // Case 10 = 1010
                start_end_sides.push_back({Side::TOP, Side::RIGHT});
                start_end_sides.push_back({Side::BOTTOM, Side::LEFT});
                break;
        }

        // start_end_sides will always have either 1 or 2 elements.
        for (const auto& start_end_side : start_end_sides) {
            // If a segment is inbound to the local block, then we add it to the
            // block's inbound_segments map. Else, we add it to the block's
            // interior_segments map.
            std::unordered_multimap<SegmentKey, Segment, pair_hash> *destination;
            if (segmentIsInboundToBlock(block, row, col, start_end_side.first)) { // inbound segment
                destination = &(block->inbound_segments);
            } else {
                destination = &(block->interior_segments);
            }

            int start_side_index = computeSideIndex(row, col, 
                    start_end_side.first);
            int end_side_index = computeSideIndex(row, col, 
                    start_end_side.second);
            destination->insert({{level, start_side_index}, 
                    Segment(start_end_side.first, start_end_side.second, 
                        start_side_index, end_side_index, level, top, 
                        left, ll, lr, ur, ul)});
        }
   }
}

Segment * getNextSegment(Block *const block, const val_t& level,
        const int next_start_side_index) {
    auto iter = block->interior_segments.find({level, next_start_side_index});
    return iter == block->interior_segments.end() ? nullptr : &(iter->second);
}

void traverseContourInBlock(const val_t& level,
        const int contour_start_side_index, Block *const block, 
        Segment *current_segment, std::unordered_multimap<SegmentKey, Contour, 
        pair_hash> *dest) {
    auto line_string = std::make_shared<std::vector<Point>>();
    int contour_end_side_index = -1;
    line_string->push_back(current_segment->start);
    line_string->push_back(current_segment->end);

    Segment *next_segment;
    while (true) {
        next_segment = getNextSegment(block, level, 
            current_segment->end_side_index);
        if (next_segment == nullptr || next_segment->visited) {
            contour_end_side_index = current_segment->end_side_index;
            break;
        }
        current_segment->visited = true;
        current_segment = next_segment;
        line_string->push_back(current_segment->end);
    }

    // Kinda ugly code-wise, but constructing the Contour in place feels
    // like the right thing to do.
    // TODO: Consider whether it can/should be cleaner.
    dest->emplace(std::make_pair<SegmentKey, Contour>(
                {level, contour_start_side_index},
                {line_string, level, contour_start_side_index, 
                contour_end_side_index}));
}

void traverseNonClosedContours(Block *const block, const val_t level) {
    Segment *current_segment;

    // Iterate over all inbound segments, i.e. segments that start on the edge 
    // of this block.
    for (auto& kv_pair : block->inbound_segments) {
        val_t level = kv_pair.first.first;
        int contour_start_side_index = kv_pair.first.second;
        current_segment = &(kv_pair.second);
        traverseContourInBlock(level, contour_start_side_index, block, 
                current_segment, &(block->inbound_contours));
    }
}

void traverseClosedContours(Block *const block, const val_t level) {
    Segment *current_segment;

    // Iterate over all segments, and do a traversal from any that remain
    // unvisited. These are the segments which must lie on some contour interior
    // to this block.
    for (auto& kv_pair : block->interior_segments) {
        current_segment = &(kv_pair.second);
        if (current_segment->visited) {
            continue;
        }

        val_t level = kv_pair.first.first;
        int contour_start_side_index = kv_pair.first.second;
        traverseContourInBlock(level, contour_start_side_index, block, 
                current_segment, &(block->interior_contours));
    }
}

void traverseNonClosedContourFragments(const val_t level, std::vector<std::list<Contour>> *output_contours) {
    return;
}

void traverseClosedContourFragments(const val_t level, std::vector<std::list<Contour>> *output_contours) {
    return;
}

// TODO: Fix after everything else has been restructured.
void countSegments(int *total_segments = nullptr, int *unvisited_segments = nullptr) {
//    int total = 0, unvisited = 0;
//    for (int idx = 0; idx < (nrows - 1) * (ncols - 1); idx++) {
//        for (auto iter = square->segments.begin(); iter != square->segments.end(); iter++) {
//            total++;
//            if (iter->second.visited == false) {
//                unvisited++;
//            }
//        }
//    }
//    // If int variables were passed in, write to them.
//    if (total_segments != nullptr && unvisited_segments != nullptr) {
//        *total_segments = total;
//        *unvisited_segments = unvisited;
//    // Else just print the results
//    } else {
//        float percent = ((float) total) / ((float) unvisited);
//        printf("Found %i segments, %i (%.2f%%) are unvisited.\n", total, unvisited, percent);
//    }
}

inline Point transformPoint(Point point) {
    point.x =  point.x * pixel_size + raster_ll.x;
    point.y = (-point.y + nrows) * pixel_size + raster_ll.y;
    return point;
}

void printGeoJSON(FILE *output, const std::vector<std::list<Contour>> output_contours) {
    fprintf(output,  "{ \"type\":\"FeatureCollection\", \"features\": [\n");
    size_t i;
    for (i = 0; i < output_contours.size(); i++) {
        const Contour *contour = &(output_contours[i].front());
        fprintf(output, "{ \"type\":\"Feature\", ");
        fprintf(output, "\"properties\": {\"level\":%.2f, \"is_closed\":%s}, ",
                contour->level, (contour->is_closed) ? "true" : "false");
        fprintf(output, "\"geometry\":{ \"type\":\"LineString\", \"coordinates\": [");
        size_t j;
        for (j = 0; j < contour->line_string->size() - 1; j++) {
            Point pt = transformPoint((*(contour->line_string))[j]);
            fprintf(output, "[%.8lf,%.8lf],", pt.x, pt.y);
        }
        Point pt = transformPoint((*(contour->line_string))[j]);
        fprintf(output, "[%.8lf,%.8lf] ] } }", pt.x, pt.y);
        fprintf(output, "%s\n", (i < output_contours.size() - 1) ? "," : "");
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
    // Note: There are (nrow - 1) x (ncols - 1) squares, i.e. one less than the respective pixel dimensions.
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

            // Phase 1
            Block *block = &blocks[block_num];
            for (int i = block->first_row; i < block->first_row + block->num_rows; i++) {
                for (int j = block->first_col; j < block->first_col + block->num_cols; j++) {
                    processSquare(block, i, j, levels);
                }
            }

            // Phase 2
            for (size_t i = 0; i < levels.size(); i++) {
                const auto& level = levels[i];
                traverseNonClosedContours(block, level);
                traverseClosedContours(block, level);
            }

        }
    }

    // Phase 3
    std::vector<std::list<Contour>> output_contours;
    size_t i;
    # pragma omp parallel default (shared) private(i)
    {
        # pragma omp for schedule(static) nowait
        for (i = 0; i < levels.size(); i++) {
            const auto& level = levels[i];
            traverseNonClosedContourFragments(level, &output_contours);
            traverseClosedContourFragments(level, &output_contours);
        }
    }

    char output_filename[100];
    sprintf(output_filename, "output_%s.txt", basename(input_filename));
    FILE *output_file = fopen(output_filename, "w");
    printGeoJSON(output_file, output_contours);
    fclose(output_file);

    return 0;
}
