#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) {
    int nrows, ncols;

    scanf("%d %d", &nrows, &ncols);
    
    float *input_array = (float *) malloc(nrows * ncols * sizeof(float));

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

    int i;
    #pragma omp parallel private(i)
    {
        #pragma omp for
        for (i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                input_array[i * ncols + j] *= 2;
            }
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
