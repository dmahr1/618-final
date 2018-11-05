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
    Point start;
    Point end;
    bool visited;
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

double interpolate(double low, double high, double level) {
    return 0.0;
};

Square createSquare(int row, int col, const std::vector<double>& levels) {
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

    for (auto& level : levels) {
        // Skip levels that are totally below or totally above this pixel.
        if (level < pixel_min || level > pixel_max) {
            continue;
        }
    
        uint8_t key = 0;
        
        // TODO: figure out right thing to do regarding tolerance
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
        switch(key) {
            case 1:
                segment0.start = {col_f, row_f + interpolate(nw, sw, level)};
                segment0.end = {col_f + interpolate(se, sw, level), row_f + 1};
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

Contour traverseSquares(Square *squares, double level, Square *square,
        Segment starting_segment) {
    Contour contour;
    contour.level = level;

    Segment current_segment = starting_segment;
    contour.line_string.push_back(current_segment.start);
    contour.line_string.push_back(current_segment.end);

    // Mark the actual original segment as visited.
    // While (true)
    // | Figure out (row, column) of next square implied by current segment.
    // | If square is out of bounds, return the existing contour.
    // | Else, determine next square's segment implied by current segment.
    // | If segment unvisited, add end point to line_string.
    // | Mark segent as visited.
    return contour;
}

void traverseNonClosedContours(Square *squares, double level, 
        std::vector<Contour> *contours) {
    // Top
    for (int col = 0; col < ncols - 1; col++) {
        Square square = squares[col];
        auto iter = square.segments.find(level);
        if (iter != square.end()) {
            continue;
        }

        for (const auto& segment : iter->second) {
            if (segment.start.y == 0.0) {
                Contour contour = traverseSquares(squares, level, segment);
                contour.is_closed = false;
            }
        }
    }

    for (int row = 0; row < nrows - 1; row++) {
        Square square = squares[row * ncols + ncols - 1];
        auto iter = square.segments.find(level);
        if (iter != square.end()) {
            continue;
        }

        for (const auto& segment : iter->second) {
            if (segment.start.x == ncols - 1) {
                // Start traversing
            }
        }
    }
    //
    // Right
    //
    // Bottom
    //
    // Left

}

void traverseClosedContours(Square *squares, double level) {
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

    Square *squares = (Square *) malloc(nrows * ncols * sizeof(Square));
    std::vector<double> levels = {0.0, 5.0, 10.0};
    for (int i = 0; i < nrows - 1; i++) {
        for (int j = 0; j < ncols - 1; j++) {
            squares[i * ncols + j] = createSquare(i, j, levels);
        }
    }

    std::vector<Contour> contours;
    for (auto& level : levels) {
        traverseNonClosedContours(squares, levels, &contours);
        traverseClosedContours(squares, levels, &contours);
    }


    return 0;
}
