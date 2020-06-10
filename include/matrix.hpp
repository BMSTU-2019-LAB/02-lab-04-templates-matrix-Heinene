// Copyright 2018 Your Name <your_email>

#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_
#include <math.h>
#include <limits>
#include <type_traits>
template <class z>
class Matrix {
 private:
  z **mat;
  int row, col;
 public:
  ~Matrix();
  int gRow() const;
  int gCol() const;
  Matrix(int rows_m, int cols_m);
  Matrix(const Matrix &copy);
  Matrix &operator=(const Matrix<z> &rhs);
  Matrix<z> operator+(Matrix<z> &m2);
  Matrix operator-(Matrix &m2);
  Matrix operator*(Matrix &m2);
  double det(Matrix m);
  Matrix deleterows_cols(Matrix m, int nrows, int ncols);
  Matrix Inversion();
  template <class Vecc>
  friend bool operator==(const Matrix<Vecc> &m1, const Matrix<Vecc> &m2);
  template <class Vec1>
  friend bool operator!=(const Matrix<Vec1> &m1, const Matrix<Vec1> &m2);
  z *operator[](size_t i) const;
};
template <class z>
Matrix<z>::~Matrix() {
  for (int i = 0; i < row; i++) {
    delete[] mat[i];
  }
  delete[] mat;
}
template <class z>
int Matrix<z>::gRow() const {
  return row;
}
template <class z>
int Matrix<z>::gCol() const {
  return col;
}
template <class z>
Matrix<z>::Matrix(int rows_m, int cols_m) {
  row = rows_m;
  col = cols_m;
  mat = new z *[rows_m];
  for (int i = 0; i < rows_m; i++) mat[i] = new z[cols_m];
  for (int i = 0; i < row; i++)
    for (int j = 0; j < cols_m; j++) mat[i][j] = 0;
}
template <class z>
Matrix<z>::Matrix(const Matrix &cp) {
  row = cp.row;
  col = cp.col;
  mat = new z *[row];
  for (int i = 0; i < row; i++) mat[i] = new z[col];
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++) mat[i][j] = cp.mat[i][j];
}
template <class z>
z *Matrix<z>::operator[](size_t i) const {
  return mat[i];
}
template <class z>
Matrix<z> &Matrix<z>::operator=(const Matrix<z> &cp) {
  if (&cp != this) {
    row = cp.row;
    col = cp.col;
    for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++) mat[i][j] = cp[i][j];
  }
  return *this;
}
template <class z>
Matrix<z> Matrix<z>::operator+(Matrix<z> &m2) {
  Matrix<z> m(0, 0);
  if (row == m2.row && col == m2.col) {
    Matrix<z> m(row, col);
    for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++) m[i][j] = mat[i][j] + m2[i][j];
    return m;
  }
  return m;
}
template <class z>
Matrix<z> Matrix<z>::operator-(Matrix &m2) {
  if (this->gCol() != m2.gCol() || this->gRow() != m2.gRow()) {
    Matrix<z> res(0, 0);
    return res;
  }
  Matrix<z> res(row, col);
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      res[i][j] = (*this)[i][j] - m2[i][j];
    }
  }
  return res;
}
template <class z>
Matrix<z> Matrix<z>::operator*(Matrix &m2) {
  if (this->gCol() != m2.gRow()) {
    Matrix<z> res(0, 0);
    return res;
  }
  Matrix<z> resul(row, m2.gCol());
  for (int i = 0; i < resul.gRow(); i++) {
    for (int j = 0; j < resul.gCol(); j++) {
      for (int k = 0; k < col; k++) {
        resul[i][j] += (*this)[i][k] * m2[k][j];
      }
    }
  }
  return resul;
}
template <class z>
bool operator==(const Matrix<z> &m1, const Matrix<z> &m2) {
  if (std::is_floating_point<z>::value) {
    for (int i = 0; i < m1.gRow(); i++) {
      for (int j = 0; j < m1.gCol(); j++) {
        if (abs(m1[i][j] - m2[i][j]) > std::numeric_limits<double>::epsilon()) {
          return false;
        }
      }
    }
	return true;
  } else {
    for (int i = 0; i < m1.gRow(); i++) {
      for (int j = 0; j < m1.gCol(); j++) {
        if (m1[i][j] != m2[i][j]) {
          return false;
        }
      }
    }
    return true;
  }
}
template <>
bool operator==(const Matrix<float> &m1, const Matrix<float> &m2) {
  for (int i = 0; i < m1.gRow(); i++) {
    for (int j = 0; j < m1.gCol(); j++) {
      if (abs(m1[i][j] - m2[i][j]) > std::numeric_limits<float>::epsilon()) {
        return false;
      }
    }
  }
  return true;
}
template <class z>
bool operator!=(const Matrix<z> &m1, const Matrix<z> &m2) {
  return !(m1 == m2);
}
template <class z>
Matrix<z> Matrix<z>::deleterows_cols(Matrix<z> mat, int nrows, int ncols) {
  Matrix<z> res(mat.gRow() - 1, mat.gCol() - 1);
  int numrow = 0;
  int numcol = 0;
  for (int i = 0; i < mat.gRow(); i++) {
    if (i != nrows) {
      for (int j = 0; j < mat.gCol(); j++) {
        if (j != ncols) {
          res[numrow][numcol] = mat[i][j];
          numcol += 1;
        } else {
          continue;
        }
      }
      numrow += 1;
      numcol = 0;
    } else {
      continue;
    }
  }
  return res;
}
template <class z>
double Matrix<z>::det(Matrix<z> matr) {
  double db = 0;
  if (matr.gRow() > 2) {
    for (int i = 0; i < matr.gRow(); i++) {
      db += pow(-1, i) * matr[0][i] * det(deleterows_cols(matr, 0, i));
    }
  } else {
    db = matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0];
  }
  return db;
}
template <class z>
Matrix<z> Matrix<z>::Inversion() {
  Matrix<z> inv(this->row, this->col);
  double det_m = det(*this);
  for (int i = 0; i < inv.gRow(); i++) {
    for (int j = 0; j < inv.gCol(); j++) {
      inv[i][j] = pow(-1, i + j) * det(deleterows_cols(*this, i, j));
    }
  }
  Matrix<z> invT(inv.gRow(), inv.gCol());
  for (int i = 0; i < invT.gRow(); i++) {
    for (int j = 0; j < invT.gCol(); j++) {
      z temp = inv[j][i];
      z det_temp = temp / det_m;
      invT[i][j] = det_temp;
    }
  }
  return invT;
}
#endif  // INCLUDE_MATRIX_HPP_