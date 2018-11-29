#include "contours.h"

#include <sstream>

using Time = std::chrono::high_resolution_clock;


//======================================================
//  GLOBALS
//======================================================

// Globals set by command line arguments
static int _argc;
static char **_argv;
const char *input_filename;
val_t interval = -1.0;  // Negative means subdivide into 10
int num_threads = 1;
int block_dim = 32;
bool skip_writing_output;

// Globals from input file's header
int nrows, ncols;
Point raster_ll;
coord_t pixel_size;

// Globals from other places
FILE *input_file;
val_t *input_array;
int nblocksh, nblocksv;
Block *blocks;
auto prev_time = Time::now();
bool DEBUG = false;


//======================================================
//  CONSTRUCTORS
//======================================================

Segment::Segment(const Side start_side, const Side end_side,
        const int square_row, const int square_col,
        const val_t level, const bool is_inbound_to_raster,
        const coord_t top, const coord_t left, const val_t ll,
        const val_t lr, const val_t ur, const val_t ul) :
        square_row(square_row), square_col(square_col), end_side(end_side),
        visited(false), is_inbound_to_raster(is_inbound_to_raster) {
    start = interpolatePoint(start_side, level, top, left, ll, lr, ur, ul);
    end = interpolatePoint(end_side, level, top, left, ll, lr, ur, ul);
}

Contour::Contour(const std::shared_ptr<std::vector<Point>>& line_string,
        const val_t& level, int start_side_idx, int end_side_idx,
        Side end_side, bool is_closed) :
        line_string(line_string), level(level), start_side_idx(start_side_idx),
        end_side_idx(end_side_idx), end_side(end_side), visited(false), is_closed(is_closed) {
}


//======================================================
//  INPUT HANDLING
//======================================================

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

bool get_option_bool(const char *option_name) {
    for (int i = _argc - 2; i >= 0; i -= 2) {
        if (strcmp(_argv[i], option_name) == 0) {
            return strcmp(_argv[i + 1], "true") == 0;
        }
    }
    return false;
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

    // Read benchmarking mode
    skip_writing_output = get_option_bool("-w");

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

std::vector<val_t> determineLevels(val_t &interval, val_t val_min, val_t val_max) {
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


//======================================================
//  PHASE 1: MARCHING SQUARES PROCESSING
//======================================================

// Process block with with first row and column at block->first_row, block->first_col.
void processBlock(Block *const block, const std::vector<val_t>& levels) {
    // Iterate over all squares in the block.
    for (int row = block->first_row; row < block->first_row + block->num_rows; row++) {
        for (int col = block->first_col; col < block->first_col + block->num_cols; col++) {
            // Process square whose top-left pixel is  at (row, col).
            // Get pixel values at corners around this square
            val_t ll = input_array[(row + 1) * ncols + col];
            val_t lr = input_array[(row + 1) * ncols + col + 1];
            val_t ur = input_array[row * ncols + col + 1];
            val_t ul = input_array[row * ncols + col];
            val_t pixel_min = std::min(ll, std::min(lr, std::min(ur, ul)));
            val_t pixel_max = std::max(ll, std::max(lr, std::max(ur, ul)));

            // Temporary buffer of segments for a single level at a single square (up to 2).
            std::vector<std::pair<Side, Side>> start_end_sides(2);

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

                start_end_sides.clear(); // Always empty sides from the previous square

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

                // Iterate over the 1 or 2 segments in start_end_sides
                for (const auto& start_end_side : start_end_sides) {

                    // Determine which map this segment should be added to, either the
                    // block's inbound_segments map or the block's interior_segments map.
                    std::unordered_map<SegmentKey, Segment, pair_hash> *destination;
                    if (segmentIsInboundToBlock(block, row, col, start_end_side.first)) { // inbound segment
                        destination = &(block->inbound_segments);
                    } else {
                        destination = &(block->interior_segments);
                    }

                    // Construct the segment and add it to the appropriate map
                    int start_side_index = computeSideIndex(row, col, start_end_side.first);
                    bool is_inbound_to_raster = segmentIsInboundToRaster(row, col, start_end_side.first);
                    destination->insert({{level, start_side_index},
                            Segment(start_end_side.first, start_end_side.second,
                                row, col, level, is_inbound_to_raster,
                                top, left, ll, lr, ur, ul)});
                }
           }
        }
    }
}

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

inline bool segmentIsInboundToBlock(const Block *const block, const int row,
        const int col, const Side side) {
    switch (side) {
        case (Side::TOP):
            return row == block->first_row;
        case (Side::BOTTOM):
            return row == block->first_row + block->num_rows - 1;
        case (Side::LEFT):
            return col == block->first_col;
        case (Side::RIGHT):
            return col == block->first_col + block->num_cols - 1;
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

inline bool segmentIsInboundToRaster(const int row, const int col, const Side side) {

    switch (side) {
        case (Side::TOP):
            return row == 0;            // Index of the first row of squares
        case (Side::BOTTOM):
            return row == nrows - 2;    // Index of the last row of squares
        case (Side::LEFT):
            return col == 0;            // Index of the first column of squares
        case (Side::RIGHT):
            return col == ncols - 2;    // Index of the last column of squares
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

inline int computeSideIndex(const int square_row, const int square_col, const Side side) {
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


//======================================================
//  PHASE 2: JOINING SEGMENTS INTO CONTOUR FRAGMENTS
//======================================================

// Traverse all contour fragments in a block at a level that are non-closed
void traverseNonClosedContourFragments(Block *const block, const val_t level) {
    Segment *current_segment;

    // Iterate over all inbound segments, i.e. segments that start on the edge
    // of this block.
    for (auto& kv_pair : block->inbound_segments) {
        val_t level = kv_pair.first.first;
        int contour_start_side_index = kv_pair.first.second;
        current_segment = &(kv_pair.second);
        traverseContourFragment(block, level, contour_start_side_index,
                current_segment);
    }
}

// Traverse all contour fragments in a block at a level that are closed
void traverseClosedContourFragments(Block *const block, const val_t level) {
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
        traverseContourFragment(block, level, contour_start_side_index,
                current_segment);
    }
}

// Traverse a contour in a block at a level starting at a particular side_idx
//   until it closes on itself or exits the block, and then add contour
//   fragment to the specified map.
void traverseContourFragment(Block *const block, const val_t& level,
        const int contour_start_side_index, Segment *current_segment) {

    // Determine which map output contour fragment should be added to
    std::unordered_map<SegmentKey, ContourFragment, pair_hash> *dest;
    if (current_segment->is_inbound_to_raster) {
        dest = &(block->inbound_fragments);
    } else {
        dest = &(block->interior_fragments);
    }

    auto line_string = std::make_shared<std::vector<Point>>();
    int contour_end_side_index = -1;
    Side contour_end_side = Side::ANY;
    line_string->push_back(current_segment->start);

    Segment *next_segment;
    bool is_closed = false;
    while (true) {
        if (DEBUG) printSegment(level, current_segment);
        int end_side_index = computeSideIndex(current_segment->square_row,
                current_segment->square_col, current_segment->end_side);
        next_segment = getNextSegment(block, level, end_side_index);
        current_segment->visited = true;
        line_string->push_back(current_segment->end);
        if (next_segment == nullptr) {
            contour_end_side_index = end_side_index;
            contour_end_side = current_segment->end_side;
            break;
        } else if (next_segment->visited) {
            is_closed = true;
            contour_end_side_index = end_side_index;
            contour_end_side = current_segment->end_side;
            break;
        }
        current_segment = next_segment;
    }

    // Kinda ugly code-wise, but constructing the ContourFragment in place feels
    // like the right thing to do.
    // TODO: Consider whether it can/should be cleaner.
    assert(contour_end_side_index != -1);
    dest->insert({{level, contour_start_side_index},
            {line_string, level, contour_start_side_index,
            contour_end_side_index, contour_end_side, is_closed}});
}

// Retrieve segment in a block at a level starting at particular side_idx
Segment * getNextSegment(Block *const block, const val_t& level,
        const int next_start_side_index) {
    auto iter = block->interior_segments.find({level, next_start_side_index});
    return iter == block->interior_segments.end() ? nullptr : &(iter->second);
}


//======================================================
//  PHASE 3: JOINING CONTOUR FRAGMENTS INTO CONTOURS
//======================================================

// Traverse contours inbound to the raster at a level by merging contour fragments
void traverseInboundContours(const val_t level,
        std::vector<std::list<ContourFragment *>> *output_contours) {

    // Iterate over all blocks
    for (int block_row = 0; block_row < nblocksv; block_row++) {
        for (int block_col = 0; block_col < nblocksh; block_col++) {

            // Skip non-perimeter blocks
            if (block_row > 0 && block_row < nblocksv - 1 &&
                    block_col > 0 && block_col < nblocksh - 1) {
                continue;
            }

            // Iterate over all inboud contour fragments, except those from other levels
            for (auto& kv_pair : blocks[block_row*nblocksh + block_col].inbound_fragments) {
                if (kv_pair.second.level != level) {
                    continue;
                }
                assert(kv_pair.second.visited == false);
                std::list<ContourFragment *> contour;
                contour = traverseContour(&(kv_pair.second), block_row, block_col);
                // Atomatically add contour to global vector
                #pragma omp critical
                {
                    output_contours->push_back(contour);
                }
            }
        }
    }
}

// Traverse contours not inbound to the raster at a level by merging contour fragments
void traverseInteriorContours(const val_t level,
        std::vector<std::list<ContourFragment *>> *output_contours) {
    // Iterate over all blocks
    for (int block_row = 0; block_row < nblocksv; block_row++) {
        for (int block_col = 0; block_col < nblocksh; block_col++) {

            // Iterate over all inboud contour fragments, except those from other levels
            for (auto& kv_pair : blocks[block_row*nblocksh + block_col].interior_fragments) {
                if (kv_pair.second.level != level || kv_pair.second.visited) {
                    continue;
                }
                // If contour fragment is closed within this block, add it as 1-element list
                std::list<ContourFragment *> contour;
                if (kv_pair.second.is_closed) {
                    contour = {&(kv_pair.second)};
                // Else perform full traversal
                } else {
                    contour = traverseContour(&(kv_pair.second), block_row, block_col);
                }
                // Atomatically add contour to global vector
                #pragma omp critical
                {
                    output_contours->push_back(contour);
                }

            }
        }
    }
}

// Retrieve ContourFragment in neighboring block
ContourFragment * getNextContour(const ContourFragment *const current_contour,
        int *block_row, int *block_col) {

    // Lookup next block based on block and end_side of current contour
    int next_block_row = *block_row;
    int next_block_col = *block_col;
    switch(current_contour->end_side) {
        case Side::LEFT:
            next_block_col--;
            break;
        case Side::RIGHT:
            next_block_col++;
            break;
        case Side::TOP:
            next_block_row--;
            break;
        case Side::BOTTOM:
            next_block_row++;
            break;
        default:
            assert(false);
    }
    // Next block is out of bounds
    if (next_block_row < 0 || next_block_row > nblocksv - 1 ||
            next_block_col < 0 || next_block_col > nblocksh - 1) {
        return nullptr;
    }

    // Lookup next contour in next block
    Block *next_block = &blocks[next_block_row * nblocksh + next_block_col];
    if (DEBUG) {
        printf("Cur block = %d,%d, cur contour at %.1f ends at side %d on %d, next block = %d,%d\n",
                *block_row, *block_col, current_contour->level,
                (int) current_contour->end_side, current_contour->end_side_idx,
                next_block_row, next_block_col);
        printContourFragments(&blocks[*block_row * nblocksh + *block_col]);
        printContourFragments(next_block);
    }

    *block_row = next_block_row;
    *block_col = next_block_col;

    return &(next_block->interior_fragments.at(
            {current_contour->level, current_contour->end_side_idx}));
}

// Traverse a contour at a level starting at a particular contour fragment
//   and return it as a linked list of contour fragments
std::list<ContourFragment *> traverseContour(ContourFragment* current_contour,
        int block_row, int block_col) {

    std::list<ContourFragment *> contour_list;
    ContourFragment *next_contour;
    while (true) {
        next_contour = getNextContour(current_contour, &block_row, &block_col);
        current_contour->visited = true;
        contour_list.push_back(current_contour);
        if (next_contour == nullptr || next_contour->visited) {
            break;
        }
        current_contour = next_contour;
    }
    return contour_list;
}


//======================================================
//  OUTPUT GENERATION
//======================================================

void printGeoJSON(FILE *output, const std::vector<std::list<ContourFragment *>> output_contours) {

    // Iterate over contours (i.e. linked lists of contour fragments)
    fprintf(output,  "{ \"type\":\"FeatureCollection\", \"features\": [\n");
    for (size_t i = 0; i < output_contours.size(); i++) {
        fprintf(output, "{ \"type\":\"Feature\", ");
        fprintf(output, "\"properties\": {\"level\":%.2f}, ", output_contours[i].front()->level);
        fprintf(output, "\"geometry\":{ \"type\":\"LineString\", \"coordinates\": [");

        // Iterate over contour fragments in contour
        size_t num_fragments = output_contours[i].size();
        size_t fragment_num = 0;
        for (ContourFragment *contour : output_contours[i]) {

            // Iterate over vertices in contour fragment, skipping 1st vertex if not 1st fragment
            size_t num_vertices = contour->line_string->size();
            for (size_t j = (fragment_num == 0) ? 0 : 1; j < num_vertices; j++) {
                Point pt = transformPoint((*(contour->line_string))[j]);
                bool is_last_point = (j == num_vertices - 1) && (fragment_num == num_fragments - 1);
                fprintf(output, "[%.8lf,%.8lf]%s", pt.x, pt.y, is_last_point ? "" : ",");
            }
            fragment_num++;

        }
        fprintf(output, "] } }");
        fprintf(output, "%s\n", (i < output_contours.size() - 1) ? "," : "");
    }
    fprintf(output, "] } \n");
}

inline Point transformPoint(Point point) {
    point.x =  point.x * pixel_size + raster_ll.x;
    point.y = (-point.y + nrows) * pixel_size + raster_ll.y;
    return point;
}

//======================================================
//  PROFILING AND DEBUGGING
//======================================================

void profileTime(std::string message) {
    auto new_time = Time::now();
    auto duration = new_time - prev_time;
    auto duration_microsecs = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    std::cout << message << " = " << duration_microsecs << " microseconds\n";
    prev_time = new_time;
}

std::string sideToString(Side side) {
    switch(side) {
        case Side::TOP:     return "TOP";
        case Side::BOTTOM:  return "BOT";
        case Side::LEFT:    return "LEFT";
        case Side::RIGHT:   return "RIGHT";
        case Side::ANY:     return "ANY";
        default:            return "OTHER";
    }
}

void printSegment(val_t level, Segment *cs) {
    printf("Segment @ level %.1f @ square (%d, %d), start (%.1f, %.1f), end (%.1f, %.1f)\n",
            level, cs->square_row, cs->square_col, cs->start.x, cs->start.y, cs->end.x, cs->end.y);
}


void printContourFragment(ContourFragment *cf) {
    printf("Contour @ level %.1f has %zu vertices from side %3d to %3d (%s): (%.1f,%.1f)->(%.1f,%.1f)\n",
        cf->level, cf->line_string->size(), cf->start_side_idx, cf->end_side_idx, sideToString(cf->end_side).c_str(),
        (*(cf->line_string))[0].x, (*(cf->line_string))[0].y,
        (*(cf->line_string))[cf->line_string->size() - 1].x,
        (*(cf->line_string))[cf->line_string->size() - 1].y);
}


void printContourFragments(Block *block) {
    printf("Contents of block with UL corner at %d,%d:\n", block->first_row, block->first_col);
    printf("  Inbound fragments (n = %zu)\n", block->inbound_fragments.size());
    for (auto& kv_pair : block->inbound_fragments) {
        printf("    ");
        printContourFragment(&(kv_pair.second));
    }
    printf("  Interior fragments (n = %zu)\n", block->interior_fragments.size());
    for (auto& kv_pair : block->interior_fragments) {
        printf("    ");
        printContourFragment(&(kv_pair.second));
    }
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
    printf("Levels = %zu intervals of size %.1f from %.1f through %.1f\n",
            levels.size(), interval, levels.front(), levels.back());
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

    int block_num;
    # pragma omp parallel default(shared) private(block_num)
    {
        # pragma omp for schedule(dynamic) nowait
        for (block_num = 0; block_num < nblocksh * nblocksv; block_num++) {
            //printf("Thread %d does phase 1 of block %d\n", omp_get_thread_num(), block_num);
            // Phase 1
            Block *block = &blocks[block_num];
            processBlock(block, levels);

            // Phase 2
            for (const auto& level : levels) {
                traverseNonClosedContourFragments(block, level);
                traverseClosedContourFragments(block, level);
            }
            
            if (DEBUG) printContourFragments(block);
        }
    }
    profileTime("Contour fragments generation");

    // Phase 3
    std::vector<std::list<ContourFragment *>> output_contours;
    size_t i;
    # pragma omp parallel default (shared) private(i)
    {
        # pragma omp for schedule(dynamic) nowait
        for (i = 0; i < levels.size(); i++) {
            const auto& level = levels[i];
            if (DEBUG) printf("PROCESSING INBOUND CONTOURS FOR LEVEL = %.1f\n", level);
            traverseInboundContours(level, &output_contours);
            if (DEBUG) printf("PROCESSING INTERIOR CONTOURS FOR LEVEL = %.1f\n", level);
            traverseInteriorContours(level, &output_contours);
        }
    }
    profileTime("Contour joining");

    if (skip_writing_output) {
        printf("Skipping writing of output file\n");
    } else {
        char output_filename[100];
        sprintf(output_filename, "output_%s.txt", basename(input_filename));
        FILE *output_file = fopen(output_filename, "w");
        printGeoJSON(output_file, output_contours);
        fclose(output_file);
        profileTime("File writing");
    }

    return 0;
}
