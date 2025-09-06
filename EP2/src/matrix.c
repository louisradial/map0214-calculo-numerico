#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

// Zero initialize a matrix, used for testing only
void zero_matrix(Matrix* m, int rows, int cols) {
    m->rows = rows;
    m->cols = cols;
    m->matrix = (double**)malloc(sizeof(double) * rows);
    for (int i = 0; i < rows; ++i) {
        m->matrix[i] = calloc(cols, sizeof(double));
    }
    return;
}

// Initialize a matrix from a given array, defined at compile time
void initialize_matrix(Matrix* m, double* data, int rows, int cols) {
    m->rows = rows;
    m->cols = cols;
    m->matrix = (double**)malloc(sizeof(double) * rows);
    for (int i = 0; i < rows; ++i) {
        m->matrix[i] = (double*)malloc(sizeof(double) * cols);
        for (int j = 0; j < cols; ++j) {
            m->matrix[i][j] = data[i * cols + j];
        }
    }
    return;
}

// Free memory after use
void destroy_matrix(Matrix* m) {
    for (int i = 0; i < m->rows; ++i) {
        free(m->matrix[i]);
    }
    free(m->matrix);
    m->rows = 0;
    m->cols = 0;
    return;
}

void print_matrix(Matrix* m) {
    for (int i = 0; i < m->rows; ++i) {
        for (int j = 0; j < m->cols; ++j) {
            printf("%+.5lf ", m->matrix[i][j]);
        }
        printf("\n");
    }
    return;
}

void print_column(Matrix *m, int col) {
    for (int i = 0; i < m->rows; ++i) {
        printf("%+.5lf ", m->matrix[i][col]);
    }
    printf("\n");
}

// As everything is on the heap, I can simply swap pointers to swap rows
// assert(row1 < m->rows && row2 < m->rows && row1 >= 0 && row2 >= 0);
void swap_rows(Matrix* m, int row1, int row2) {
    double* row = m->matrix[row1];
    m->matrix[row1] = m->matrix[row2];
    m->matrix[row2] = row;
    return;
}
