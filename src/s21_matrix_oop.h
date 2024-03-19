#include <cmath>
#include <exception>
#include <iostream>

#define EPS 1e-6
#define SUCCESS 1
#define FAILURE 0

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& other) noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() noexcept;
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  void SetSize(int rows, int cols);
  void SetRows(int rows);
  void SetCols(int cols);
  int GetRows() const noexcept;
  int GetCols() const noexcept;

  bool operator==(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const double number);
  friend S21Matrix operator*(const double number, const S21Matrix& matrix);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double number);
  S21Matrix& operator*=(const S21Matrix& other);
  double& operator()(const int i, const int j);

 private:
  int rows_, cols_;
  double** matrix_;
  void AllocMatrix();
  void DeallocMatrix();
  double MinorElement(const int i, const int j);
  void FindDeterminant(const int i, const int j, S21Matrix& minor);
  double CalcDeterminant();
};
