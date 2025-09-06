typedef struct {
    int rows;
    int cols; 
    double **matrix;
} Matrix;

// Zero initialize a matrix
void zero_matrix(Matrix* m, int rows, int cols);

// Initialize a matrix from a given array, defined at compile time
void initialize_matrix(Matrix* m, double* data, int rows, int cols);

// Free memory after use
void destroy_matrix(Matrix* m); 

void print_matrix(Matrix* m);

void print_column(Matrix *m, int col);

// Swaps rows row1 and row2
// assert(row1 < m->rows && row2 < m->rows && row1 >= 0 && row2 >= 0);
void swap_rows(Matrix* m, int row1, int row2);
