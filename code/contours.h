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
#include <map>
#include <vector>

#include <assert.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

#include "gdal_priv.h"

enum class Side {TOP, BOTTOM, LEFT, RIGHT, ANY};
constexpr uint8_t UL_ABOVE = 8;
constexpr uint8_t UR_ABOVE = 4;
constexpr uint8_t LR_ABOVE = 2;
constexpr uint8_t LL_ABOVE = 1;

typedef double val_t;
typedef double coord_t;

typedef struct {
    coord_t x;
    coord_t y;
} Point;

typedef struct Segment {
    Segment(const Side start_side, const Side end_side,
            const int square_row, const int square_col,
            const val_t level, const bool is_inbound_to_raster,
            const coord_t top, const coord_t left, const val_t ll,
            const val_t lr, const val_t ur, const val_t ul);
    int square_row, square_col;
    Side end_side;
    bool visited, is_inbound_to_raster;
    Point start, end;
} Segment;

typedef std::pair<val_t, int> SegmentKey;

typedef struct Contour {
    Contour(const std::shared_ptr<std::vector<Point>>& line_string,
            const val_t& level, int start_side_idx, int end_side_idx,
            Side end_side, bool is_closed);

    std::shared_ptr<std::vector<Point>> line_string;
    val_t level;
    int start_side_idx;
    int end_side_idx;
    Side end_side;
    bool visited;
    // TODO: Consider whether this is still necessary?
    bool is_closed;
} ContourFragment;

typedef struct {
    // Using pair hash from https://stackoverflow.com/a/32685618/1254992.
    std::size_t operator () (const std::pair<val_t, int> &p) const {
        auto h1 = std::hash<val_t>{}(p.first);
        auto h2 = std::hash<int>{}(p.second);
        return h1 ^ h2;
    }
} pair_hash;

typedef struct {
    int first_row;
    int first_col;
    int num_rows;
    int num_cols;
    std::unordered_map<SegmentKey, Segment, pair_hash> interior_segments;
    std::unordered_map<SegmentKey, Segment, pair_hash> inbound_segments;
    std::unordered_map<SegmentKey, ContourFragment, pair_hash> inbound_fragments;
    std::unordered_map<SegmentKey, ContourFragment, pair_hash> interior_fragments;
} Block;

// Input handling
int get_option_int(const char *option_name, int default_value);
const char *get_option_string(const char *option_name, const char *default_value);
float get_option_float(const char *option_name, float default_value);
bool get_option_bool(const char *option_name);
void readArguments(int argc, char **argv);
void readHeader(FILE *input);
std::vector<val_t> determineLevels(val_t &interval, val_t val_min, val_t val_max);

// Phase 1: Marching squares
void processSquare(Block *const block, const int row, const int col,
        const std::vector<val_t>& levels);
val_t interpolate(const val_t left_or_top, const val_t right_or_bottom, const val_t level);
Point interpolatePoint(const Side side, const val_t level, const coord_t top,
        const coord_t left, const val_t ll, const val_t lr, const val_t ur, const val_t ul);
bool segmentIsInboundToBlock(const Block *const block, const int row,
        const int col, const Side side);
bool segmentIsInboundToRaster(const int row, const int col, const Side side);
int computeSideIndex(int square_row, int square_col, Side side);

// Phase 2: joining segments into contour fragments
void traverseNonClosedContourFragments(Block *const block, const val_t level);
void traverseClosedContourFragments(Block *const block, const val_t level);
void traverseContourFragment(Block *const block, const val_t& level,
        const int contour_start_side_index, Segment *current_segment);
Segment * getNextSegment(Block *const block, const val_t& level,
        const int next_start_side_index);
void smoothLineString(std::shared_ptr<std::vector<Point>> *line_string, bool is_closed);
void simplifyLineString(std::shared_ptr<std::vector<Point>> *line_string_ptr, bool is_closed);

// Phase 3: joining contour fragments into contours
ContourFragment * getNextContour(const ContourFragment *const current_contour,
        int *block_row, int *block_col);
std::list<ContourFragment *> traverseContour(ContourFragment* current_contour,
        int block_row, int block_col);
void traverseInboundContours(const val_t level,
        std::vector<std::list<ContourFragment *>> *output_contours);
void traverseInteriorContours(const val_t level,
        std::vector<std::list<ContourFragment *>> *output_contours);

// Output generation
void printGeoJSON(FILE *output,
        const std::vector<std::list<ContourFragment *>> output_contours);
Point transformPoint(Point point);

// Profiling and debugging
void profileTime(std::string message,
        std::chrono::time_point<std::chrono::high_resolution_clock> &prev_time);
void printSegment(val_t level, Segment *segment);
void printContourFragment(ContourFragment *cf);
void printContourFragments(Block *block);
std::string sideToString(Side side);
