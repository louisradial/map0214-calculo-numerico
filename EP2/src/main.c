#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"

#define MAX_ITERATIONS 100

// Partial pivotting
// assert(col < m->rows);
void pivot(Matrix* m, int col) {
    int row = col;
    int max = fabs(m->matrix[row][col]);
    int cur;
    // find pivot by finding the maximum entry
    for (int i = row; i < m->rows; ++i) {
        cur = fabs(m->matrix[i][col]);
        if (max < cur) {
            max = cur;
            row = i;
        }
    }
    swap_rows(m, col, row);
    return;
}

// assert(m->rows < m->cols)
void gaussian_elimination(Matrix* m) {
    double mult = 1;
    // forward elimination
    printf("forward elimination\n");
    for (int row = 0; row < m->rows; ++row) {
        print_matrix(m);
        printf("row %d - partial pivotting\n", row + 1);
        pivot(m, row);
        print_matrix(m);
        printf("forward elimination\n");
        for (int i = row + 1; i < m->rows; ++i) {
            mult = -m->matrix[i][row]/m->matrix[row][row];
            for (int j = 0; j < m->cols; ++j) {
                m->matrix[i][j] += mult*m->matrix[row][j];
            }
        }
    }
    printf("row echelon form\n");
    print_matrix(m);
    // backward substitution
    for (int row = m->rows -1; row >= 0; --row) {
        for (int i = 0; i < row; ++i) {
            mult = -m->matrix[i][row]/m->matrix[row][row];
            for (int j = 0; j < m->cols; ++j) {
                m->matrix[i][j] += mult*m->matrix[row][j];
            }
        }
        mult = 1/m->matrix[row][row];
        for (int j = 0; j < m-> cols; ++j) {
            m->matrix[row][j] *= mult;
        }
    }
    printf("backward substitution\n");
    print_matrix(m);
    return;
}

void jacobi(Matrix* m, double* x, double precision) {
    // solve Ax = b system with the b vector given by m[][rows]
    assert(m->rows >= m->cols - 1);
    for (int i = 0; i < m->rows; ++i) {
        assert(m->matrix[i][i] != 0); // make sure diagonal is invertible
    }
    int n = 0; // iteration counter
    double diff, maxdiff; // halt condition control
    double xnext; // temp variable
    double* xprev = (double*)malloc(sizeof(double) * m->rows);
    // Jacobi method
    do {
        maxdiff = 0;
        for (int i = 0; i < m->rows; ++i) {
            xprev[i] = x[i]; // update xprev
        }
        printf("%2d ", n++);
        for (int i = 0; i < m->rows; ++i) {
            printf("%+.5lf ", xprev[i]);
            xnext = m->matrix[i][m->rows]; // b_i
            for (int j = 0; j < i; ++j) {
                xnext -= m->matrix[i][j]*xprev[j]; // -aij xj
            }
            for (int j = i+1; j < m->rows; ++j) {
                xnext -= m->matrix[i][j]*xprev[j]; // -aij xj
            }
            xnext /= m->matrix[i][i];
            diff = fabs(xnext - xprev[i]);
            if (diff > maxdiff) {
                maxdiff = diff; // update halt condition
            }
            x[i] = xnext; // update x
        }
        printf("%.5lf\n", maxdiff);
    } while (n < MAX_ITERATIONS && maxdiff > precision);
    printf("%2d ", n);
    for (int i = 0; i < m->rows; ++i) {
        printf("%+.5lf ", x[i]);
    }
    printf("\n");
    free(xprev);
    return;
}

void gauss_seidel(Matrix* m, double* x, double precision) {
    // solve Ax = b system with the b vector given by m[][rows]
    assert(m->rows >= m->cols - 1);
    for (int i = 0; i < m->rows; ++i) {
        assert(m->matrix[i][i] != 0); // make sure diagonal is invertible
    }
    int n = 0; // iteration counter
    double diff, maxdiff; // halt condition control
    double xnext; // temp variable
    double* xprev = (double*)malloc(sizeof(double) * m->rows);
    // Gauss-Seidel method
    do {
        maxdiff = 0;
        for (int i = 0; i < m->rows; ++i) {
            xprev[i] = x[i]; // update xprev
        }
        printf("%2d ", n++);
        for (int i = 0; i < m->rows; ++i) {
            printf("%+.5lf ", xprev[i]);
            xnext = m->matrix[i][m->rows]; // b_i
            for (int j = 0; j < i; ++j) {
                // next x[j] has already been computed
                xnext -= m->matrix[i][j]*x[j]; // -aij xj
            }
            for (int j = i+1; j < m->rows; ++j) {
                xnext -= m->matrix[i][j]*xprev[j]; // -aij xj
            }
            xnext /= m->matrix[i][i];
            diff = fabs(xnext - xprev[i]);
            if (diff > maxdiff) {
                maxdiff = diff; // update halt condition
            }
            x[i] = xnext; // update x
        }
        printf("%.5lf\n", maxdiff);
    } while (n < MAX_ITERATIONS && maxdiff > precision);
    printf("%2d ", n);
    for (int i = 0; i < m->rows; ++i) {
        printf("%+.5lf ", x[i]);
    }
    printf("\n");
    free(xprev);
    return;
}

int main() {
    int rows = 3;
    int cols = rows + 1;
    double a[] = {
        0, 5.3, -1.9, 3.1,
        11.9, 0, 1.9, 15.2,
        1, -1, -1, 0
    };
    // item b - gaussian elimination
    Matrix item_b;
    initialize_matrix(&item_b, a, rows, cols);
    print_matrix(&item_b);
    gaussian_elimination(&item_b);
    print_column(&item_b, cols - 1);
    destroy_matrix(&item_b);
    // item c - jacobi method
    Matrix item_c;
    double x[] = {1, 1, 1};
    initialize_matrix(&item_c, a, rows,cols);
    swap_rows(&item_c, 0,1); // swap first two rows
    jacobi(&item_c, x, 1e-3);
    destroy_matrix(&item_c);
    // item d - jacobi method
    Matrix item_d;
    double y[] = {1, 1, 1};
    initialize_matrix(&item_d, a, rows,cols);
    swap_rows(&item_d, 0,1); // swap first two rows
    gauss_seidel(&item_d, y, 1e-3);
    destroy_matrix(&item_d);
    return 0;
}
