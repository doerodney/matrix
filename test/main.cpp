#include <gtest.h>
#include <stdio.h>
extern "C" {
#include "matrix.h"
}


TEST(MatrixTest, DetectSingularity) {
  int isSingular = 0;
  const int nrows = 3, ncols = 3;
  double data[] = {0, 0, 0, 1, 1, 1, 2, 2, 2};

  Matrix *a = matrix_new(nrows, ncols);
  for (int row = 0, i = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      matrix_set_value(a, row, col, data[i++]);
    }
  }

  isSingular = matrix_test_singular(a);
  EXPECT_TRUE(isSingular);

  matrix_free(&a);
  EXPECT_EQ(a, nullptr);
}


TEST(MatrixTest, DetectNonSingularity) {
  int isSingular = 0;
  const int nrows = 3, ncols = 3;
  const double data[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};

  Matrix *a = matrix_new(nrows, ncols);
  for (int row = 0, i = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      matrix_set_value(a, row, col, data[i++]);
    }
  }
  isSingular = matrix_test_singular(a);
  EXPECT_FALSE(isSingular);

  matrix_free(&a);
  EXPECT_EQ(a, nullptr);
}


TEST(MatrixTest, DeterminantOfNonSquareMatrix) {
  const int nrows = 2, ncols = 3;
  Matrix *a = matrix_new(nrows, ncols);
  const double data[] = {6, 1, 1, 4, -2, 5};

  for (int row = 0, i = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      matrix_set_value(a, row, col, data[i++]);
    }
  }
  double det = 0.0;
  int failure = matrix_get_determinant(a, &det);
  EXPECT_EQ(failure, MATRIX_DATA_NOT_SQUARE);

  matrix_free(&a);
  EXPECT_EQ(a, nullptr);
}


TEST(MatrixTest, DeterminantOfEmptyMatrix) {
  const int nrows = 0, ncols = 0;
  double det = 0.0;
  Matrix *a = matrix_new(nrows, ncols);
  int failure = matrix_get_determinant(a, &det);
  EXPECT_EQ(failure, MATRIX_NULL_POINTER);

  matrix_free(&a);
  EXPECT_EQ(a, nullptr);
}


TEST(MatrixTest, DeterminantOfValidMatrix) {
  double det = 0.0;
  const int nrows = 3, ncols = 3;
  const double data[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
  double target = -306.0;

  Matrix *a = matrix_new(nrows, ncols);
  for (int row = 0, i = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      matrix_set_value(a, row, col, data[i++]);
    }
  }

  int failure = matrix_get_determinant(a, &det);
  EXPECT_EQ(failure, MATRIX_NO_ERR);
  EXPECT_DOUBLE_EQ(det, target);

  matrix_free(&a);
  EXPECT_EQ(a, nullptr);
}

TEST(MatrixTest, GetMinorMatrix_happy_path) {
  const int nrows = 3, ncols = 3;
  Matrix *in = matrix_new(nrows, ncols);
  Matrix *out = matrix_new(nrows - 1, ncols - 1);
  double data[] = {3, 0, 2, 2, 0, -2, 0, 1, 1};
  int i = 0;

  matrix_load_by_row(in, data);
  int failure = matrix_get_minor_matrix(in, 0, 0, out);
  double value = 0.0;
  
  EXPECT_EQ(failure, MATRIX_NO_ERR);
  double nw_results[] = {0, -2, 1, 1};
  for (int row = 0; row < out->nrows; row++) {
    for (int col = 0; col < out->ncols; col++) {
      value = matrix_get_value(out, row, col);
      EXPECT_DOUBLE_EQ(value, nw_results[i++]); 
    }
  }
}


/*---
TEST(MatrixTest, GetMinors_happy_path) {
  const int nrows = 3, ncols = 3;
  Matrix *in = matrix_new(nrows, ncols);
  Matrix *out = matrix_new(nrows, ncols);

  double data[] = {3, 0, 2, 2, 0, -2, 0, 1, 1};
  double minors[] = {2, 2, 2, -2, 3, 3, 0, -10, 0};
  int i = 0;

  matrix_load_by_row(in, data);

  int failure = matrix_get_minors(in, out);
  EXPECT_EQ(failure, MATRIX_NO_ERR);

  for (int row = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      double observed = matrix_get_value(out, row, col);
      double expected = minors[i];
      i++;
      printf("observed, expected = %g, %g\n", observed, expected);
      EXPECT_DOUBLE_EQ(expected, observed);
    }
  }

  matrix_free(&in);
  matrix_free(&out);
}
---*/


TEST(MatrixTest, SolveSimultaneousEquations) {
  const int nrows = 3, ncols = 3;
  double solution[] = {1.0, 5.0, 10.0};
  Matrix *a = matrix_new(nrows, ncols);
  Matrix *b = matrix_new(nrows, 1);
  Matrix *x = matrix_new(nrows, 1);

  const double a_data[] = {6, 1, 1, 4, -2, 5, 2, 8, 7};
  const double b_data[] = {21, 44, 112};

  for (int row = 0, i = 0; row < nrows; row++) {
    for (int col = 0; col < ncols; col++) {
      matrix_set_value(a, row, col, a_data[i++]);
    }
  }

  for (int row = 0, col = 0, i = 0; row < nrows; row++) {
    matrix_set_value(b, row, col, b_data[i++]);
  }

  int failure = matrix_solve_simeq(a, x, b);
  EXPECT_EQ(failure, MATRIX_NO_ERR);

  for (int row = 0, col = 0; row < nrows; row++) {
    double value = matrix_get_value(x, row, col);
    EXPECT_DOUBLE_EQ(value, solution[row]);
  }

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&x);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(x, nullptr);
}


TEST(MatrixTest, MatrixCrossProduct_NegativeUnitVectorResult) {
  int nrows = 1, ncols = 3;

  Matrix* a = matrix_new(nrows, ncols);
  Matrix* b = matrix_new(nrows, ncols);
  Matrix* out = matrix_new(nrows, ncols);

  matrix_set_value(a, 0, 0, 1);
  matrix_set_value(a, 0, 1, 0);
  matrix_set_value(a, 0, 2, 0);

  matrix_set_value(b, 0, 0, 0);
  matrix_set_value(b, 0, 1, -1);
  matrix_set_value(b, 0, 2, 0);

  double solution[] = {0, 0, -1};

  int failure = matrix_cross_product(a, b, out);
  EXPECT_EQ(failure, MATRIX_NO_ERR);

  // Assert results in out matrix:
  int i = 0;
  for (int row = 0, col = 0; col < ncols; col++) {
    double value = matrix_get_value(out, row, col);
    EXPECT_DOUBLE_EQ(value, solution[i]);
    i++;
  }

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&out);
}


TEST(MatrixTest, MatrixCrossProduct_NilResult) {
  int nrows = 1, ncols = 3;

  Matrix* a = matrix_new(nrows, ncols);
  Matrix* b = matrix_new(nrows, ncols);
  Matrix* out = matrix_new(nrows, ncols);

  matrix_set_value(a, 0, 0, 1);
  matrix_set_value(a, 0, 1, 0);
  matrix_set_value(a, 0, 2, 0);

  matrix_set_value(b, 0, 0, -1);
  matrix_set_value(b, 0, 1, 0);
  matrix_set_value(b, 0, 2, 0);

  double solution[] = {0, 0, 0};

  int failure = matrix_cross_product(a, b, out);
  EXPECT_EQ(failure, MATRIX_NO_ERR);

  // Assert results in out matrix:
  int i = 0;
  for (int row = 0, col = 0; col < ncols; col++) {
    double value = matrix_get_value(out, row, col);
    EXPECT_DOUBLE_EQ(value, solution[i]);
    i++;
  }

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&out);
}


TEST(MatrixTest, MatrixCrossProduct_MATRIX_INVALID_COLUMN_COUNT) {
  int nrows = 1, ncols = 3;

  Matrix* a = matrix_new(nrows, ncols);
  Matrix* b = matrix_new(nrows, (ncols - 1));
  Matrix* out = matrix_new(nrows, ncols);

  int failure = matrix_cross_product(a, b, out);
  EXPECT_EQ(failure, MATRIX_INVALID_COLUMN_COUNT);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&out);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(out, nullptr);
}


TEST(MatrixTest, MatrixCrossProduct_MATRIX_INVALID_ROW_COUNT) {
  int nrows = 1, ncols = 3;

  Matrix* a = matrix_new(nrows, ncols);
  Matrix* b = matrix_new((nrows + 1), ncols);
  Matrix* out = matrix_new(nrows, ncols);

  int failure = matrix_cross_product(a, b, out);
  EXPECT_EQ(failure, MATRIX_INVALID_ROW_COUNT);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&out);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(out, nullptr);
}


TEST(MatrixTest, MatrixMultiply) {
  const int a_rows = 2, a_cols = 3;
  const int b_rows = 3, b_cols = 2;
  const int c_rows = 2, c_cols = 2;

  Matrix *a = matrix_new(a_rows, a_cols);
  Matrix *b = matrix_new(b_rows, b_cols);
  Matrix *c = matrix_new(c_rows, c_cols);
  double solution[] = {58.0, 64.0, 139.0, 154.0};

  const double a_data[] = {1, 2, 3, 4, 5, 6};
  const double b_data[] = {7, 8, 9, 10, 11, 12};

  // Load A matrix:
  for (int row = 0, i = 0; row < a_rows; row++) {
    for (int col = 0; col < a_cols; col++) {
      matrix_set_value(a, row, col, a_data[i++]);
    }
  }

  // Load B matrix:
  for (int row = 0, i = 0; row < b_rows; row++) {
    for (int col = 0; col < b_cols; col++) {
      matrix_set_value(b, row, col, b_data[i++]);
    }
  }

  int failure = matrix_multiply(a, b, c);
  EXPECT_EQ(failure, MATRIX_NO_ERR);

  // Assert results in C matrix:
  int i = 0;
  for (int row = 0; row < c_rows; row++) {
    for (int col = 0; col < c_cols; col++) {
      double value = matrix_get_value(c, row, col);
      EXPECT_DOUBLE_EQ(value, solution[i]);
      i++;
    }
  }

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&c);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(c, nullptr);
}


TEST(MatrixTest, MatrixMultiply_MATRIX_INCOMPATIBLE) {
  // Test MATRIX_INCOMPATIBLE assertion:
  Matrix* a = matrix_new(4, 3);
  Matrix* b = matrix_new(4, 3);
  Matrix* c = matrix_new(2, 2);

  int failure = matrix_multiply(a, b, c);
  EXPECT_EQ(failure, MATRIX_INCOMPATIBLE);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&c);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(c, nullptr);
}


TEST(MatrixTest, MatrixMultiply_MATRIX_INVALID_ROW_COUNT) {
  // Test MATRIX_INVALID_ROW_COUNT assertion:
  Matrix* a = matrix_new(4, 3);
  Matrix* b = matrix_new(3, 4);
  Matrix* c = matrix_new(2, 2);
 
  int failure = matrix_multiply(a, b, c);
  EXPECT_EQ(failure, MATRIX_INVALID_ROW_COUNT);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&c);

  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(c, nullptr);
}



TEST(MatrixTest, MatrixMultiply_MATRIX_INVALID_COLUMN_COUNT) {
  Matrix* a = matrix_new(4, 3);
  Matrix* b = matrix_new(3, 4);
  Matrix* c = matrix_new(4, 2);

  int failure = matrix_multiply(a, b, c);
  EXPECT_EQ(failure, MATRIX_INVALID_COLUMN_COUNT);

  matrix_free(&a);
  matrix_free(&b);
  matrix_free(&c);
  
  EXPECT_EQ(a, nullptr);
  EXPECT_EQ(b, nullptr);
  EXPECT_EQ(c, nullptr);
}
 

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
