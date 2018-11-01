#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) {
    int nrows, ncols;

    scanf("%d %d", &nrows, &ncols);
    
    float *input_array = (float *) malloc(nrows * ncols * sizeof(float));

    int val;
    int index = 0;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            scanf("%d", &val);
            input_array[index] = val;
            index += 1;
        }
    }

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            printf("%.2f ", input_array[i * ncols + j]);
        }
        printf("\n");
    }

    return 0;
}
