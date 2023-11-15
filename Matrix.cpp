#include <iostream>
#include <windows.h> // Russian fonts for console
#include <random> // rand(), srand() functions
#include <time.h> // clock(), time() functions
#define EPS 1e-15
#define MIN -9
#define MAX 9

using namespace std;

class Matrix {
public:
	double* matrix_pointer; ///???
	int row, col;
	double** matrix;

	Matrix() {
		row = 1;
		col = 1;
		matrix_pointer = new double;
		*matrix_pointer = 0;
		matrix = new double* [row];
		for (int i = 0; i < row; i++) {
			matrix[i] = new double[col];
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = 0;
			}
		}
	}
	Matrix(int N) {
		row = N;
		col = N;
		matrix_pointer = new double;
		*matrix_pointer = 0;
		matrix = new double* [row];
		for (int i = 0; i < row; i++) {
			matrix[i] = new double[col];
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = 0;
			}
		}
	}
	Matrix(int M, int N) {
		row = M;
		col = N;
		matrix_pointer = new double;
		*matrix_pointer = 0;
		matrix = new double* [row];
		for (int i = 0; i < row; i++) {
			matrix[i] = new double[col];
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = 0;
			}
		}
	}
	Matrix(int M, int N, double** mtrx) {
		row = M;
		col = N;
		matrix_pointer = new double;
		*matrix_pointer = 0;
		matrix = mtrx;
		/*
		matrix = new double* [row];
		for (int i = 0; i < row; i++) {
			matrix[i] = new double[col];
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = mtrx[i][j];
			}
		}
		*/
	}
	~Matrix() {}

	void create_randox_matrix() {
		srand(time(NULL) + rand()); // set seed for random as current date
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = MIN + (double)(rand()) / ((double)(RAND_MAX / (MAX - MIN))); // inputs random double from MIN to MAX into matrix
			}
		}
	}

	void print() {
		cout.precision(2);
		cout << "Матрица " << row << "x" << col << endl;
		for (int i = 0; i < row; i++) {
			cout << fixed;
			for (int j = 0; j < col; j++) {
				cout << matrix[i][j] << "\t";
			}
			cout << endl;
		}
	}

	double return_element(int i, int j) {
		if ((i < row) && (j < col)) {
			return matrix[i][j];
		}
		cout << "Индекс элемента вышел за пределы матрицы" << endl;
		return -1;
	}

	Matrix addition(Matrix B) {
		if ((row != B.row) || (col != B.col)){
			Matrix C;
			cout << "Количество строк или столбцов в матрицах не совпадает" << endl;
			return C;
		}
		Matrix C(row, col);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				C.matrix[i][j] = matrix[i][j] + B.matrix[i][j];
			}
		}
		return C;
	}

	Matrix subtraction(Matrix B) {
		if ((row != B.row) || (col != B.col)) {
			Matrix C;
			cout << "Количество строк или столбцов в матрицах не совпадает" << endl;
			return C;
		}
		Matrix C(row, col);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				C.matrix[i][j] = matrix[i][j] - B.matrix[i][j];
			}
		}
		return C;
	}

	Matrix multiplication(Matrix B) {
		if (col != B.row) {
			Matrix C;
			cout << "Количество столбцов первой матрицы не совпадает с количеством строк второй матрицы" << endl;
			return C;
		}
		Matrix C(row, B.col);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < B.col; ++j) {
				for (int k = 0; k < col; ++k) {
					C.matrix[i][j] += matrix[i][k] * B.matrix[k][j];
				}
			}
		}
		return C;
	}

	void getCofactor(Matrix A, Matrix temp, int p, int q, int n) {
		int i = 0, j = 0;
		// Looping for each element of the matrix
		for (int r = 0; r < n; r++) {
			for (int c = 0; c < n; c++) {
				//  Copying into temporary matrix only those element
				//  which are not in given row and column
				if (r != p && c != q) {
					temp.matrix[i][j++] = A.matrix[r][c];
					// Row is filled, so increase row index and
					// reset col index
					if (j == n - 1) {
						j = 0;
						i++;
					}
				}
			}
		}
	}

	double determinant(Matrix A, int n) {
		if (A.row != A.col) {
			cout << "матрица не квадратная, невозможно вычислить определитель" << endl;
			return 0;
		}
		double D = 0; // Initialize result
		//  Base case : if matrix contains single element
		if (n == 1)
			return A.matrix[0][0];

		Matrix temp(A.row); // To store cofactors
		int sign = 1;  // To store sign multiplier
		 // Iterate for each element of first row
		for (int f = 0; f < n; f++) {
			// Getting Cofactor of A[0][f]
			getCofactor(A, temp, 0, f, n);
			D += double(sign * A.matrix[0][f] * determinant(temp, n - 1));
			// terms are to be added with alternate sign
			sign = -sign;
		}
		return D;
	}

	void adjoint(Matrix A, Matrix adj) {
		if ((A.row != A.col) || (adj.row != adj.col)) {
			cout << "Матрица не квадратная, невозможно построить матрицу алгебраических дополнений" << endl;
			return;
		}
		int N = A.row;
		if (N == 1)
		{
			adj.matrix[0][0] = 1;
			return;
		}
		// temp is used to store cofactors of A[][]
		int sign = 1;
		Matrix temp(A.row);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				// Get cofactor of A[i][j]
				getCofactor(A, temp, i, j, N);
				// sign of adj[j][i] positive if sum of row
				// and column indexes is even.
				sign = ((i + j) % 2 == 0) ? 1 : -1;
				// Interchanging rows and columns to get the
				// transpose of the cofactor matrix
				adj.matrix[j][i] = (sign) * (determinant(temp, N - 1));
			}
		}
	}

	bool inverse(Matrix A, Matrix inverse) {
		// Find determinant of A[][]
		if ((A.row != A.col) || (inverse.row != inverse.col)) {
			cout << "Матрица не квадратная, невозможно найти обратную матрицу" << endl;
			return false;
		}
		int N = A.row;
		double det = determinant(A, N);
		if (fabs(det) < EPS) {
			cout << "Матрица вырожденная, невозможно найти обратную матрицу" << endl;
			return false;
		}
		// Find adjoint
		Matrix adj(N);
		adjoint(A, adj);
		// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				inverse.matrix[i][j] = adj.matrix[i][j] / double(det);
		return true;
	}

	Matrix division(Matrix B) {
		Matrix inverse(B.row);
		//Matrix C(row, B.col);
		Matrix C = multiplication(inverse);
		return C;
	}
	void multiplication_by_number(double num) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = matrix[i][j] * num;
			}
		}
	}
	void addition_by_number(double num) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				matrix[i][j] = matrix[i][j] + num;
			}
		}
	}
};

int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	Matrix M;
	M.print();

	cout << "Создание трех матриц разной размерности, заполненных случайными числами:\n" << endl;
	Matrix five_by_two(5, 2);
	five_by_two.create_randox_matrix();
	five_by_two.print();
	Matrix three_by_five(3, 5);
	three_by_five.create_randox_matrix();
	three_by_five.print();
	Matrix new_three_by_two(3, 2);
	new_three_by_two.create_randox_matrix();
	new_three_by_two.print();

	cout << "\nУмножение двух матриц разной (но подходящей) размерности:\n3х5 * 5х2 =" << endl;
	Matrix result = three_by_five.multiplication(five_by_two);
	result.print();
	cout << "\nУмножение полученной матрицы на число:\n3x2 * 0.5 =" << endl;
	result.multiplication_by_number(0.5);
	result.print();
	cout << "\nСложение полученной матрицы с числом:\n3x2 + 100 =" << endl;
	result.addition_by_number(100);
	result.print();
	cout << "\nСложение двух матриц:\n3x2 + 3x2 =" << endl;
	Matrix add = result.addition(new_three_by_two);
	add.print();
	cout << "\nСложение двух матриц:\n3x2 - 3x2 =" << endl;
	Matrix sub = add.subtraction(result);
	sub.print();

	system("Pause");
	return 0;
}