/*Код, который я написал выше, является программой на языке C++, которая использует линейную алгебру для вычисления интегралов, двойных интегралов и интегралов по контуру. Рассмотрим его подробнее.

Первым шагом в программе я определяю структуру Vector для хранения вектора и функции для создания и освобождения памяти, занятой вектором. Вектор представляет собой одномерный массив элементов типа double, а его размерность задается переменной dimension. Функция createVector принимает размерность и массив элементов и возвращает созданный вектор, а функция freeVector освобождает память, занятую вектором.

Затем я определяю функции для вычисления скалярного произведения двух векторов (dot) и нормы вектора (norm). Скалярное произведение двух векторов определяется как сумма произведений соответствующих элементов каждого вектора, а норма вектора вычисляется как квадратный корень из суммы квадратов элементов вектора.

Далее я определяю структуру Matrix для хранения матрицы и функции для создания и освобождения памяти, занятой матрицей, а также для умножения матрицы на вектор. Матрица представляет собой двумерный массив элементов типа double, а ее размерность задается переменными rowDimension (количество строк) и columnDimension (количество столбцов). Функция createMatrix принимает размерности и двумерный массив элементов и возвращает созданную матрицу, а функция freeMatrix освобождает память, занятую матрицей. Функция multiply умножает матрицу на вектор, используя стандартный алгоритм умножения матриц.

После этого я определяю функции для вычисления определенного интеграла методом прямоугольников (definiteIntegral), двойного интеграла методом прямоугольников (doubleIntegral) и интеграла по контуру методом прямоугольников (contourIntegral). В каждой из этих функций используются функции из предыдущего кода для работы с векторами и матрицами.

Функция definiteIntegral вычисляет определенный интеграл функции f на интервале [a, b] методом прямоугольников. Для этого интервал разбивается на n равных частей, и на каждой части вычисляется значение функции. Значения функции на каждой части умножаются на ширину части (то есть на dx = (b - a) / n) и складываются, чтобы получить приближенное значение интеграла.

Функция doubleIntegral вычисляет двойной интеграл функции f на прямоугольной области [a1, b1] x [a2, b2] методом прямоугольников. Для этого область разбивается на сетку из n1 строк и n2 столбцов, и на каждой ячейке сетки вычисляется значение функции. Значения функции на каждой ячейке умножаются на площадь ячейки (то есть на dx = (b1 - a1) / n1 * dy = (b2 - a2) / n2) и складываются, чтобы получить приближенное значение интеграла.

*/
#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdio>

using namespace std;

// Structure for storing a vector
struct Vector {
    int dimension;
    double* elements;
};

// Function to create a vector
Vector* createVector(int dimension, double* elements) {
    Vector* result = new Vector;
    result->dimension = dimension;
    result->elements = new double[dimension];
    for (int i = 0; i < dimension; i++) {
        result->elements[i] = elements[i];
    }
    return result;
}

// Function to free the memory occupied by a vector
void freeVector(Vector* v) {
    delete[] v->elements;
    delete v;
}

// Function to calculate the dot product of two vectors
double dot(Vector* v1, Vector* v2) {
    assert(v1->dimension == v2->dimension);

    double result = 0;
    for (int i = 0; i < v1->dimension; i++) {
        result += v1->elements[i] * v2->elements[i];
    }

    return result;
}

// Function to calculate the norm of a vector
double norm(Vector* v) {
    double result = 0;
    for (int i = 0; i < v->dimension; i++) {
        result += v->elements[i] * v->elements[i];
    }

    return sqrt(result);
}

// Structure for storing a matrix
struct Matrix {
    int rowDimension;
    int columnDimension;
    double** elements;
};

// Function to create a matrix
Matrix* createMatrix(int rowDimension, int columnDimension, double** elements) {
    Matrix* result = new Matrix;
    result->rowDimension = rowDimension;
    result->columnDimension = columnDimension;
    result->elements = new double*[rowDimension];
    for (int i = 0; i < rowDimension; i++) {
        result->elements[i] = new double[columnDimension];
        for (int j = 0; j < columnDimension; j++) {
            result->elements[i][j] = elements[i][j];
        }
    }
    return result;
}

// Function to free the memory occupied by a matrix
void freeMatrix(Matrix* m) {
    for (int i = 0; i < m->rowDimension; i++) {
        delete[] m->elements[i];
    }
    delete[] m->elements;
    delete m;
}

// Function to multiply a matrix by a vector
Vector* multiply(Matrix* m, Vector* v) {
    assert(m->columnDimension == v->dimension);

    Vector* result = new Vector;
    result->dimension = m->rowDimension;
    result->elements = new double[m->rowDimension];

    for (int i = 0; i < m->rowDimension; i++) {
        result->elements[i] = 0;
        for (int j = 0; j < m->columnDimension; j++) {
            result->elements[i] += m->elements[i][j] * v->elements[j];
        }
    }

    return result;
}

// Function to calculate the definite integral using the rectangle method
double definiteIntegral(double a, double b, int n, double (*f)(double)) {
    if (n < 1) {
        return -1;
    }

    if (n > 1000000) {
        return -2;
    }

    double dx = (b - a) / n;
    double result = 0;
    for (int i = 0; i < n; i++) {
        double x = a + i * dx;
        result += f(x) * dx;
    }
    return result;
}

// Function to calculate the double integral using the rectangle method
double doubleIntegral(double a1, double b1, int n1, double a2, double b2, int n2, double (*f)(double, double)) {
    if (n1 < 1 || n2 < 1) {
        return -1;
    }

    if (n1 > 1000 || n2 > 1000) {
        return -2;
    }

    double dx = (b1 - a1) / n1;
    double dy = (b2 - a2) / n2;
    double result = 0;
    for (int i = 0; i < n1; i++) {
        double x = a1 + i * dx;
        for (int j = 0; j < n2; j++) {
            double y = a2 + j * dy;
            result += f(x, y) * dx * dy;
        }
    }
    return result;
}

// Function to calculate the contour integral using the rectangle method
double contourIntegral(Vector* contour, int n, double (*f)(double)) {
    if (n < 1) {
        return -1;
    }

    if (n > 10000) {
        return -2;
    }

    double a = contour->elements[0];
    double b = contour->elements[1];
    double dx = (b - a) / n;
    double result = 0;
    for (int i = 0; i < n; i++) {
        double x = a + i * dx;
        double y = f(x);
        Vector* v = createVector(2, new double[2] {x, y});
        double length = norm(v);
        freeVector(v);
        result += length * dx;
    }
    return result;
}

// Function to calculate the indefinite integral of a function
double indefiniteIntegral(double x, double (*f)(double)) {
    // Perform numerical integration with a step size of 0.0001
    double dx = 0.0001;
    double result = 0;
    for (double t = 0; t <= x; t += dx) {
        result += f(t) * dx;
    }
    return result;
}

// Function to calculate the improper integral of a function
double improperIntegral(double a, double (*f)(double)) {
    // Upper limit is infinity, using a large value as approximation
    double infinity = 1000000;
    double dx = 0.001;
    double result = 0;
    for (double t = a; t < infinity; t += dx) {
        result += f(t) * dx;
    }
    return result;
}

// Testing the functions
int main() {
    int choice;

    do {
        cout << "Menu:" << endl;
        cout << "1. Calculate definite integral" << endl;
        cout << "2. Multiply matrix by vector" << endl;
        cout << "3. Calculate contour integral" << endl;
        cout << "4. Calculate indefinite and improper integrals" << endl;
        cout << "5. Exit" << endl;
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case 1: {
                double a, b;
                int n;

                cout << "Enter the lower limit of the integral: ";
                cin >> a;
                cout << "Enter the upper limit of the integral: ";
                cin >> b;
                cout << "Enter the number of subdivisions (n): ";
                cin >> n;

                double result = definiteIntegral(a, b, n, sin);
                cout << "Definite integral: " << result << endl;
                break;
            }
            case 2: {
                int rowDimension, columnDimension;

                cout << "Enter the row dimension of the matrix: ";
                cin >> rowDimension;
                cout << "Enter the column dimension of the matrix: ";
                cin >> columnDimension;

                double** matrixElements = new double*[rowDimension];
                for (int i = 0; i < rowDimension; i++) {
                    matrixElements[i] = new double[columnDimension];
                    for (int j = 0; j < columnDimension; j++) {
                        cout << "Enter the element at position (" << i << ", " << j << "): ";
                        cin >> matrixElements[i][j];
                    }
                }

                double* vectorElements = new double[columnDimension];
                for (int i = 0; i < columnDimension; i++) {
                    cout << "Enter the element of the vector at position " << i << ": ";
                    cin >> vectorElements[i];
                }

                Matrix* m = createMatrix(rowDimension, columnDimension, matrixElements);
                Vector* v = createVector(columnDimension, vectorElements);
                Vector* result = multiply(m, v);

                cout << "Result: (";
                for (int i = 0; i < result->dimension; i++) {
                    cout << result->elements[i];
                    if (i < result->dimension - 1) {
                        cout << ", ";
                    }
                }
                cout << ")" << endl;

                freeMatrix(m);
                freeVector(v);
                freeVector(result);

                delete[] matrixElements;
                delete[] vectorElements;

                break;
            }
            case 3: {
                double start, end;
                int n;

                cout << "Enter the starting point of the contour: ";
                cin >> start;
                cout << "Enter the ending point of the contour: ";
                cin >> end;
                cout << "Enter the number of subdivisions (n): ";
                cin >> n;

                Vector* contour = createVector(2, new double[2]{start, end});
                double result = contourIntegral(contour, n, sin);
                cout << "Contour integral: " << result << endl;

                freeVector(contour);
                break;
            }
            case 4: {
                double x;

                cout << "Enter the value of x for the indefinite integral: ";
                cin >> x;

                double resultIndefinite = indefiniteIntegral(x, sin);
                cout << "Indefinite integral: " << resultIndefinite << endl;

                double resultImproper = improperIntegral(0, sin);
                cout << "Improper integral: " << resultImproper << endl;

                break;
            }
            case 5: {
                cout << "Exiting..." << endl;
                break;
            }
            default: {
                cout << "Invalid choice! Please try again." << endl;
                break;
            }
        }
        cout << endl;

    } while (choice != 5);

    return 0;
}
