#ifndef INC_MATRIX_H
#define INC_MATRIX_H

#include <math.h> //NAN

#define MATRIX_NO_ERR 0
#define MATRIX_NULL_POINTER 1
#define MATRIX_MEMORY_ERR 2
#define MATRIX_NULL_ARG 3
#define MATRIX_DATA_NOT_SQUARE 4
#define MATRIX_UNEQUAL_ROW_COUNTS 5
#define MATRIX_INVALID_COLUMN_COUNT 6
#define MATRIX_INVALID_ROW_COUNT 7
#define MATRIX_SINGULAR 8
#define MATRIX_INCOMPATIBLE 9

typedef struct {
  int nrows;
  int ncols;
  double data[];

} Matrix;



Matrix* matrix_new(int nrows, int ncols);

Matrix* matrix_copy(const Matrix* src);

int matrix_get_determinant(const Matrix* p, double *det);

void matrix_free(Matrix **pp);

int matrix_cross_product(const Matrix* a, const Matrix* b, Matrix* out);

int matrix_get_minor_matrix(const Matrix *in, int mrow, int mcol, Matrix *out);

int matrix_get_minors(const Matrix *in, Matrix *out);

double matrix_get_value(const Matrix *p, int row, int col);

void matrix_load_by_column(Matrix *m, double in[]);

void matrix_load_by_row(Matrix *m, double in[]);

int matrix_multiply(const Matrix *a, const Matrix *b, Matrix *c);

void matrix_set_value(Matrix *p, int row, int col, double value);

int matrix_solve_simeq(const Matrix *a, Matrix* x, const Matrix *b);

int matrix_test_singular(const Matrix *m);

#define MATRIX_LOAD_BY_COLUMN(mtrx, ...) matrix_load_by_column(mtrx, (double[]) {__VA_ARGS__, NAN})

#define MATRIX_LOAD_BY_ROW(mtrx, ...) matrix_load_by_row(mtrx, (double[]) {__VA_ARGS__, NAN})

#endif


