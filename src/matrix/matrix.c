#include <math.h> //NAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

// local functions;
int fixColumnIndex(const Matrix *m, int idx);
void getColumnContent(Matrix *dest, const Matrix *src, int col);
void getRowContent(Matrix *dest, const Matrix *src, int row);
int getDatumIndex(const Matrix *m, int row, int col);
void printMatrix(const Matrix *m);
void setColumnContent(Matrix *dest, const Matrix *src, int col);


// Matrix constructor:
Matrix* matrix_new(int nrows, int ncols) {
  Matrix  *out = NULL;
  if ((nrows * ncols) > 0) {

    // Allocate memory for the matrix:
    size_t data_byte_count = nrows * ncols * sizeof(double);
    out = (Matrix*) malloc(sizeof(Matrix) + data_byte_count);

    if (out) {
      *out = (Matrix) {.nrows = nrows, .ncols = ncols};
      memset((void*) out->data, 0, data_byte_count);
    }
  }

  return out;
}


// Create a copy of a matrix:
Matrix* matrix_copy(const Matrix* src)
{
  Matrix *dest = matrix_new(src->nrows, src->ncols);
  if (dest) {
    memcpy((void *) dest->data, (void * ) src->data, src->nrows * src->ncols * sizeof(src->data[0]));
  }

  return dest;
}


// Calculate cross-product of two vectors, represented as 1 x 3 matrices:
int matrix_cross_product(const Matrix* a, const Matrix* b, Matrix* out) {
  int failure = MATRIX_NO_ERR;
  double value = 0.0;

  // Test matrix dimensions:
  if ((a->ncols != 3) || (b->ncols != 3)) {
    failure = MATRIX_INVALID_COLUMN_COUNT;
  } else if ((a->nrows != 1) || (b->nrows != 1)) { 
    failure = MATRIX_INVALID_ROW_COUNT;
  } else {
    value = matrix_get_value(a, 0, 1) * matrix_get_value(b, 0, 2) - 
            matrix_get_value(a, 0, 2) * matrix_get_value(b, 0, 1);
    matrix_set_value(out, 0, 0, value);

    value = matrix_get_value(a, 0, 2) * matrix_get_value(b, 0, 0) -
            matrix_get_value(a, 0, 0) * matrix_get_value(b, 0, 2);
    matrix_set_value(out, 0, 1, value);

    value = matrix_get_value(a, 0, 0) * matrix_get_value(b, 0, 1) -
            matrix_get_value(a, 0, 1) * matrix_get_value(b, 0, 0);
    matrix_set_value(out, 0, 2, value);
  }

  return failure;
}


// Matrix destructor:
void matrix_free(Matrix **pp) {
  Matrix *pm = *pp;
  if (pm) {   
    // Free struct storage:
    free((void*) *pp);

    // Set pointer to NULL:
    *pp = NULL;
  }
}


// Calculate the determinant of a matrix:
int matrix_get_determinant(const Matrix* m, double *det) {
  int failure = MATRIX_NO_ERR;
  double diag, value = 0;
  int col, row, fixed = 0;

  // Test for null pointers in arguments:
  if ((NULL == m) || (NULL == det)) {
    failure = MATRIX_NULL_POINTER;
  } else if ((m->ncols == 0) || (m->nrows == 0) || (m->ncols != m->nrows)) {
    failure = MATRIX_DATA_NOT_SQUARE;
  } else {

    // Continue if the matrix has data and is square:
    *det = 0.0;

    // Special treatment for 2x2 matrix:
    if (m->nrows == 2) {
      *det = (matrix_get_value(m, 0, 0) * matrix_get_value(m, 1, 1) -
          matrix_get_value(m, 0, 1) * matrix_get_value(m, 1, 0));

    } else {

      // Apply the sum of the left->right diagonals:
      for (col = 0; col < m->ncols; col++) {
        diag = 1.0;
        for (row = 0; row < m->nrows; row++) {
          fixed = fixColumnIndex(m, col + row);
          value = matrix_get_value(m, row, fixed);
          diag *= value;
        }
        *det += diag;
      }

      // Apply the sum of the right->left diagonals:
      for (col = m->ncols - 1; col >= 0; col--) {
        diag = 1.0;
        for (row = 0; row < m->nrows; row++) {
          fixed = fixColumnIndex(m, col - row);
          value = matrix_get_value(m, row, fixed);
          diag *= value;
        }
        *det -= diag;
      }
    }
  } 
  
  return failure;
}

// Calculate the determinant of a matrix:
int matrix_get_det(const Matrix* m, double *det) {
  int failure = MATRIX_NO_ERR;
  double diag, value = 0;
  int col, row, fixed = 0;
  // Test for null pointers in arguments:
  if ((NULL == m) || (NULL == det)) {
    failure = MATRIX_NULL_POINTER;
  } else if ((m->nrows != m->ncols) || (m->nrows == 0) || (m->ncols == 0)){
    failure = MATRIX_DATA_NOT_SQUARE;
  } else {
    // Create an augmented matrix:
    int naug_cols = ((2 * m->ncols) - 1);
    Matrix* aug = matrix_new(m->nrows, naug_cols);
    Matrix *content = matrix_new(m->nrows, 1);
    int src_col = 0;
    
    // Copy the original matrix into the augmented matrix:
    for (row = 0; row < m->nrows; row++) {
      for (col = 0; col < m->ncols; col++) {
        value = matrix_get_value(m, row, col);
        matrix_set_value(aug, row, col, value);
      }
    }

    // Load original matrix columns into the augmented matrix columns:
    for (col = m->ncols; col < naug_cols; col++) {
      getColumnContent(content, m, src_col);
      setColumnContent(aug, content, col);
      src_col++;
    }

    // Calculate sum of left-to-right diagonals in the augmented matrix:
    *det = 0.0;
    diag = 1.0;
    for (col = 0; col < m->ncols; col++) {
      diag = 1.0;
      for (row = 0; row < m->nrows; row++) {
        value = matrix_get_value(aug, row, col + row);
        diag *= value;
      }
      *det += diag;
    }

    // Calculate sum of the right-to-left diagonals in the augmented matrix:
    for (col = (m->ncols - 1); col < naug_cols; col++){
      diag = 1.0;
      for (row = 0; row < m->nrows; row++) {
        value = matrix_get_value(aug, row, col - row);
        diag *= value;
      }
      *det -= diag;
    }

    matrix_free(&aug);
    matrix_free(&content);
  }

  return failure;
}


// Get the minor matrix for a given row, col:
int matrix_get_minor_matrix(const Matrix *in, int mrow, int mcol, Matrix *out) {
  int failure = MATRIX_NO_ERR;

  // Test for null pointers in arguments:
  if ((NULL == in) || (NULL == out)) {
    failure = MATRIX_NULL_POINTER;
  }
  // Test if in and out matrix are the correct size:
  else if (((in->nrows - out->nrows) != 1) || ((in->ncols - out->ncols) != 1)) {
    failure = MATRIX_INCOMPATIBLE;
  } 
  // Test if in and out are square: 
  else if ((in->ncols == 0) || (in->nrows == 0) || 
           (in->ncols != in->nrows) || 
           (out->ncols != out->nrows)) {
    failure = MATRIX_DATA_NOT_SQUARE;
  } else {
      int i = 0;
      double data[out->nrows * out->ncols];

      // Collect the minor matrix for (mrow, mcol):
      for (int row = 0; row < in->nrows; row++) {
        for (int col = 0; col < in->nrows; col++) {
          if ((row != mrow) && (col != mcol)) {
            data[i] = matrix_get_value(in, row, col);
            i++;
          }
        }
      }

      // Load the data into the out matrix:
      matrix_load_by_row(out, data);
  }

  return failure;
}


// Get the matrix of minors of a square matrix:
int matrix_get_minors(const Matrix *in, Matrix *out) {
  int failure = MATRIX_NO_ERR;

  // Test for null pointers in arguments:
  if ((NULL == in) || (NULL == out)) {
    failure = MATRIX_NULL_POINTER;
  }
  // Test if in and out matrix are the same size:
  else if ((in->nrows != out->nrows) || (in->ncols != out->ncols)) {
    failure = MATRIX_INCOMPATIBLE;
  } 
  // Test if in and out are square: 
  else if ((in->ncols == 0) || (in->ncols != in->nrows)) {
    failure = MATRIX_DATA_NOT_SQUARE;
  } 
  else {

    // Create a minor matrix:
    Matrix *minor = matrix_new((in->nrows - 1), (in->ncols - 1));
    double det = 0.0;
    int mrow, mcol;
    
    // Calculate a minor value for each position in the in matrix:
    for (mrow = 0; mrow < in->nrows; mrow++) {
      for (mcol = 0; mcol < in->ncols && !failure; mcol++) {

        failure = matrix_get_minor_matrix(in, mrow, mcol, minor);
        if (!failure) {
          // Get the determinant of the minor:
          failure = matrix_get_determinant(minor, &det);
        }

        if (!failure) {
          matrix_set_value(out, mrow, mcol, det);
        }
      }
    }

    // Free the minor matrix:
    matrix_free(&minor);

  }
  return failure;
}


// Get the value of a location in a matrix:
double matrix_get_value(const Matrix *m, int row, int col) {
  int idx = getDatumIndex(m, row, col);

  double value = m->data[idx];

  return value;
}


// Load the data into a matrix by columns:
void matrix_load_by_column(Matrix *m, double in[]) {
  int idx = 0;
  for (int col = 0; col < m->ncols && !isnan(in[idx]); col++) {
    for (int row = 0; row < m->nrows; row++) {
      idx = getDatumIndex(m, row, col);
      m->data[idx] = in[idx];
      idx++;
    }
  }
}


// Load the data into a matrix by rows:
void matrix_load_by_row(Matrix *m, double in[]) {
  int idx = 0;
  for (int row = 0; row < m->nrows && !isnan(in[idx]); row++) {
    for (int col = 0; col < m->ncols; col++) {
      idx = getDatumIndex(m, row, col);
      m->data[idx] = in[idx];
      idx++;
    }
  }
}


// Set the value of a location in a matrix:
void matrix_set_value(Matrix *p, int row, int col, double value) {
  int idx = getDatumIndex(p, row, col);

  p->data[idx] = value;
}


// Multiply matrix A by matrix B, result goes into matrix C.
// Number of columns in A must equal number of rows in B.
// C has the number of rows in A and the number of columns in B.
// Notation:  https://en.wikipedia.org/wiki/Matrix_multiplication
// Example: https://www.mathsisfun.com/algebra/matrix-multiplying.html
int matrix_multiply(const Matrix *a, const Matrix *b, Matrix *c) {
  int failure = MATRIX_NO_ERR;
  int row, col, i;
  double a_value, b_value, value = 0.0;
  
  // Validate incoming matrices:
  if (a->ncols != b->nrows) {
    failure = MATRIX_INCOMPATIBLE;
  } else if (c->nrows != a->nrows) {
    failure = MATRIX_INVALID_ROW_COUNT;
  } else if (c->ncols != b->ncols) {
    failure = MATRIX_INVALID_COLUMN_COUNT;
  } else {

    // Iterate rows of matrix A:
    for (row = 0; row < a->nrows; row++) {
      
      // Iterate columns of matrix B:
      for (col = 0; col < b->ncols; col++) {

        // Calculate value for (row, col) in matrix C:
        value = 0.0;
        for (i = 0; i < a->ncols; i++) {
          a_value = matrix_get_value(a, row, i);
          b_value = matrix_get_value(b, i, col);

          value += (a_value * b_value);
        }
        matrix_set_value(c, row, col, value);
      }
    }
  }

  return failure;
}


// Solve a set of simultaneous equations of the form AX = B:
int matrix_solve_simeq(const Matrix *a, Matrix *x, const Matrix *b) {
  int failure = MATRIX_NO_ERR;

  // Validate incoming matrices:
  if (a->nrows != a->ncols) {
    failure = MATRIX_DATA_NOT_SQUARE;
  } else if (a->nrows != x->nrows) {
    failure = MATRIX_UNEQUAL_ROW_COUNTS;
  } else if (a->nrows != b->nrows) {
    failure = MATRIX_UNEQUAL_ROW_COUNTS;
  } else if (matrix_test_singular(a)) {
    failure = MATRIX_SINGULAR;
  }
  else {
    // Set the determinant of matrix a as the denominator:
    double denom = 0.0;
    failure = matrix_get_determinant(a, &denom);
    if (!failure) {
      double det = 0.0;

      Matrix *aprime = matrix_copy(a);
      Matrix *content = matrix_new(a->nrows, 1);

      // Iterate through columns in a matrix.
      // Substitute the b matrix for the column col.
      // Take the determinant of the resultant matrix.
      // Implement Cramer's Rule.
      for (int col = 0; col < a->ncols && !failure; col++) {
        // Preserve matrix a column content:
        getColumnContent(content, aprime, col);

        // Substitute b matrix for the column in aprime:
        setColumnContent(aprime, b, col);

        // Take determinant of resultant matrix:
        failure = matrix_get_determinant(aprime, &det);
        if (!failure) {

          // Record result in x matrix:
          matrix_set_value(x, col, 0, det/denom);

          // Restore column in matrix a':
          setColumnContent(aprime, content, col);
        }
      }

      matrix_free(&aprime);
      matrix_free(&content);
    }
  }

  return failure;
}


// Test if the matrix is singular:
int matrix_test_singular(const Matrix *m) {
  int singular = 0;

  double det = 0.0; 
  (void) matrix_get_determinant(m, &det);
  if (det == 0.0) {
    singular = 1;
  }

  return singular;
}


// Fix the column value to be within 0 <= col < ncols:
int fixColumnIndex(const Matrix *m, int idx) {
  int fixed = idx;
  while (fixed >= m->ncols) {
    fixed -= m->ncols;
  }
  while (fixed < 0) {
    fixed += m->ncols;
  }

  return fixed;
}


// Get the content of a column in a matrix.
// This is returned as an nrows x 1 matrix.
void getColumnContent(Matrix *dest, const Matrix *src, int col) {

  for (int row = 0; row < src->nrows; row++) {
    double value = matrix_get_value(src, row, col);
    matrix_set_value(dest, row, 0, value);
  }

}


// Get the content of a row in a matrix.
// This is return as an ncols X 1 matrix.
void getRowContent(Matrix *dest, const Matrix *src, int row) {

  for (int col = 0; col < src->ncols; col++) {
    double value = matrix_get_value(src, row, col);
    matrix_set_value(dest, 0, col, value);
  }

}


// Get the index of a specified row/column in the matrix data field:
int getDatumIndex(const Matrix *m, int row, int col) {
  int idx = row * m->ncols + col;

  return idx;
}


// Set the content of a column in a matrix:
void setColumnContent(Matrix *dest, const Matrix *src, int col) {

  for (int row = 0; row < dest->nrows; row++) {
    double value = matrix_get_value(src, row, 0);
    matrix_set_value(dest, row, col, value);
  }
}


void printMatrix(const Matrix *m) {
  for (int row = 0; row < m->nrows; row++) {
    for (int col = 0; col < m->ncols; col++) {
      double value = matrix_get_value(m, row, col);
      printf("%10.2g", value);
    }
    printf("\n");
  }
}
