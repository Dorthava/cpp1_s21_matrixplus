#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(3), cols_(3), matrix_(nullptr) { AllocMatrix(); }

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_(nullptr) {
  AllocMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other) : S21Matrix() { *this = other; }

S21Matrix::S21Matrix(S21Matrix&& other) : S21Matrix() {
  *this = std::move(other);
}

S21Matrix::~S21Matrix() { DeallocMatrix(); }

void S21Matrix::AllocMatrix() {
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument(
        "cols_ and rows_ variables cannot be less than or equal to zero");
  }
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::DeallocMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::SetSize(int rows, int cols) {
  S21Matrix trash = *this;
  DeallocMatrix();
  if (rows != rows_) {
    rows_ = rows;
  } else if (cols != cols_) {
    cols_ = cols;
  }
  AllocMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] =
          i >= trash.rows_ || j >= trash.cols_ ? 0. : trash.matrix_[i][j];
    }
  }
}

void S21Matrix::SetRows(int rows) { SetSize(rows, cols_); }

void S21Matrix::SetCols(int cols) { SetSize(rows_, cols); }

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

bool S21Matrix::EqMatrix(const S21Matrix& other) noexcept {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return FAILURE;
  }
  bool error = SUCCESS;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        error = FAILURE;
      }
    }
  }
  return error;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "the number of columns of the first matrix is not equal to the "
        "number "
        "of rows of the second matrix");
  }

  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      for (int k = 0; k < other.rows_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(result);
}

void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

S21Matrix S21Matrix::Transpose() noexcept {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) result.matrix_[j][i] = matrix_[i][j];
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument("the matrix is not square");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = 1;
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        result.matrix_[i][j] = std::pow(-1, i + j) * MinorElement(i, j);
      }
    }
  }
  return result;
}

double S21Matrix::MinorElement(const int i, const int j) {
  S21Matrix minor(rows_ - 1, cols_ - 1);
  FindDeterminant(i, j, minor);
  return minor.CalcDeterminant();
}

void S21Matrix::FindDeterminant(const int i, const int j, S21Matrix& minor) {
  int resRow = 0, resCol = 0;
  for (int k = 0; k < rows_; k++) {
    for (int m = 0; m < cols_; m++) {
      if (k != i && m != j) {
        minor.matrix_[resRow][resCol++] = matrix_[k][m];
      }
      if (resCol == minor.cols_) {
        resCol = 0;
        resRow += 1;
      }
    }
  }
}

double S21Matrix::CalcDeterminant() {
  double result = 0.;
  if (rows_ == 1) {
    result = matrix_[0][0];
  } else {
    S21Matrix new_matrix(rows_ - 1, cols_ - 1);
    for (int i = 0; i < cols_; i++) {
      FindDeterminant(0, i, new_matrix);
      result = i % 2 ? result - matrix_[0][i] * new_matrix.CalcDeterminant()
                     : result + matrix_[0][i] * new_matrix.CalcDeterminant();
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument("the matrix is not square");
  }
  double result = 0.;
  if (cols_ == 1) {
    result = matrix_[0][0];
  } else {
    result = CalcDeterminant();
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determ_result = Determinant();
  if (fabs(determ_result) < EPS) {
    throw std::invalid_argument("the matrix is not square");
  }
  S21Matrix new_matrix = CalcComplements().Transpose();
  S21Matrix result = std::move(new_matrix * (1 / determ_result));
  return result;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SubMatrix(other);
  return result;
}
S21Matrix S21Matrix::operator*(const double number) {
  S21Matrix result = *this;
  result.MulNumber(number);
  return result;
}
S21Matrix operator*(const double number, const S21Matrix& matrix) {
  S21Matrix result = matrix;
  result.MulNumber(number);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result = *this;
  result.MulMatrix(other);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    DeallocMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    AllocMatrix();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) matrix_[i][j] = other.matrix_[i][j];
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    DeallocMatrix();
    rows_ = other.rows_;
    other.rows_ = 0;
    cols_ = other.cols_;
    other.cols_ = 0;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double number) {
  this->MulNumber(number);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

double& S21Matrix::operator()(const int i, const int j) {
  if (!(i < rows_ && i >= 0 && j < cols_ && j >= 0)) {
    throw std::invalid_argument("index is outside the matrix");
  }
  return matrix_[i][j];
}