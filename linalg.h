#ifndef _LINALG_H_
#define _LINALG_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>



// define "USING_STATIC_ARRAY" when you use static array instead of STACK_ALLOC or variable-length array(VLA)
#define USE_STATIC_ARRAY





#ifdef USE_STATIC_ARRAY 

// set static buffer size
#define LINALG_MATRIX_BUFFER_SIZE (25)
#define LINALG_VECTOR_BUFFER_SIZE (5)

#else

#include <malloc.h>
// define the function name that allocates stack memory.
#define STACK_ALLOC _alloca

#endif // USING_STATIC_ARRAY


#define LINALG_MAX_ITER (50)
#define LINALG_MAXOFFDIAG_EPS (1e-6f)

// Define Matrix_t structure
typedef struct {
    float* data;
    unsigned int rows;
    unsigned int cols;

} Matrix_t;

#define Mat_at(matrix, rows_, cols_) ((matrix).data[(rows_) * (matrix).cols + (cols_)])

/**
 * @brief Create new Matrix object
 * 
 * @param rows number of rows
 * @param cols number of columns
 * @return Matrix_t new Matrix object
 */
Matrix_t CreateMatrix(unsigned int rows, unsigned int cols);


/**
 * @brief Copy a matrix object with a new Matrix object
 * 
 * @param source matrix object to input
 * @return Matrix_t new Matrix object
 */
Matrix_t CopyMatrix(const Matrix_t* source);


/**
 * @brief Create new Matrix object and it will be initilized with input value
 * 
 * @param rows number of rows
 * @param cols number of columns
 * @param init_value input value to initilize the Matrix
 * @return Matrix_t new Matrix object 
 */
Matrix_t CreateMatrixWithInitValue(unsigned int rows, unsigned int cols, float init_value);


/**
 * @brief Create a Identity Matrix object
 * 
 * @param size matrix size
 * @return Matrix_t new Matrix object 
 */
Matrix_t CreateIdentityMatrix(unsigned int size);


/**
 * @brief Create an Matrix array
 * 
 * @param pp_mat a pointer of pointer to Matrix 
 * @param rows number of rows of array element Matrix
 * @param cols number of columns of array element Matrix
 * @param array_size length of Matrix array
 * @return int error code : returns 0 if executed successfully
 */
int CreateMatrixArray(Matrix_t** pp_mat, unsigned int rows, unsigned int cols, size_t array_size);

/**
 * @brief Assign values to matrix elements with float array elements
 * 
 * @param matrix  pointer to Matrix to be assigned
 * @param array an array to assign
 * @param size length of the array
 */
void AssignMatrix(Matrix_t* matrix, float* _array, unsigned int size);


/**
 * @brief Assign a value to all matrix elements
 * 
 * @param matrix pointer to Matrix to be assigned
 * @param value an value to assign
 */
void AssignMatrixWithValue(Matrix_t* matrix, float value);

// Function to assign a float values to a matrix.
/**
 * @brief Assign identity to matrix
 * 
 * @param matrix pointer to Matrix
 */
void AssignIdentityMatrix(Matrix_t* matrix);

/**
 * @brief Get the Row Vector
 * 
 * @param source pointer to Matrix to be extracted
 * @param destination pointer to Matrix object to be assigned row vector elements
 * @param n_cols number of extracting columns
 * @return int error code : returns 0 if executed successfully
 */
int GetRowVector(const Matrix_t* source, Matrix_t* destination, unsigned int n_cols);


/**
 * @brief Get the Column Vector object
 * 
 * @param source pointer to Matrix object to be extracted
 * @param destination pointer to Matrix object to be assigned column vector elements
 * @param n_rows number of extracting rows
 * @return int error code : returns 0 if executed successfully
 */
int GetColVector(const Matrix_t* source, Matrix_t* destination, unsigned int n_rows);


/**
 * @brief Free allocated heap memory held Matrix object
 * 
 * @param matrix pointer to Matrix object to free
 */
void FreeMatrix(Matrix_t* matrix);

/**
 * @brief Free allocated heap memory held Matrix array. after that, the array object will be free and assigned its pointer.
 * 
 * @param array pointer to Matrix array
 * @param size length of Matrix array
 */
void FreeMatrixArray(Matrix_t** mat_array, size_t size);

/**
 * @brief Operate Matrix addition on Matrix a and Matrix b. (result = a + b)
 * 
 * @param a pointer to  input Matrix object (left side)
 * @param b pointer to another input Matrix object (right side)
 * @param result pointer to output Matrix object
 * @return int error code : returns 0 if executed successfully
 */
int AddMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result);

/**
 * @brief Operate Matrix substraction on Matrix a and Matrix b. (resukt = a - b)
 * 
 * @param a pointer to input Matrix object that will be substracted (left side)
 * @param b pointer to another input Matrix object that will substract Matrix a (right side)
 * @param result pointer to output Matrix object
 * @return int error code : returns 0 if executed successfully
 */
int SubtractMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result);

/**
 * @brief Operate weighted appending with input Matrix and input float value (target += weight * input)
 * 
 * @param target pointer to output Matrix object
 * @param input pointer to  input Matrix
 * @param weight a weighted value
 * @return int error code : returns 0 if executed successfully
 */
int AppendMatrixWithWeight(Matrix_t* target, const Matrix_t* input, float weight);

/**
 * @brief Operate Matrix Transposition
 * 
 * @param input pointer to input Matrix object
 * @param output pointer to output Matrix object
 * @return int error code : returns 0 if executed successfully
 */
int TransposeMatrix(const Matrix_t* input, Matrix_t* output);

/**
 * @brief Operate Matrix multiplication
 * 
 * @param a pointer to input Matrix object (left side)
 * @param b pointer to another input Matrix object (right side)
 * @param result pointer to output Matrix object
 * @return int error code : returns 0 if executed successfully
 */
int MultiplyMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result);

/**
 * @brief Operate Matrix multiplication of on a matrix and scalar (output = scalar * input)
 * 
 * @param input pointer to input Matrix object
 * @param scalar  an input float value
 * @param output pointer to output Matrix object
 * @return int 
 */
int MultiplyMatrixByScalar(const Matrix_t* input, float scalar, Matrix_t* output);

/**
 * @brief Operate Matrix invertion
 * 
 * @param input pointer to input Matrix object
 * @param output pointer to output Matrix object
 * @return int error code : returns 0 if executed successfully
 */
int InvertMatrix(const Matrix_t* input, Matrix_t* output);


/**
 * @brief Operate the Cholesky decomposition
 * 
 * @param source pointer to input Matrix object
 * @param destination pointer to output Matrix object: It will be assigned lower triangular matrix
 * @return int error code : returns 0 if executed successfully
 */
int CholeskyDecomposition(const Matrix_t* source, Matrix_t* destination);

// function to decomposition from source matrix to LU-marged matrix and Pivot matrix
/**
 * @brief Operate LU decomposition
 * 
 * @param source pointer to input Matrix object
 * @param LU pointer to output Matrix object: It will be assigned marged-LU matrix
 * @param P pointer to output Matrix object: It will be assigned pivot selection matrix
 * @return int error code : returns 0 if executed successfully
 * @sa SolveFromPLU
 */
int LU_Decomposition(const Matrix_t* source, Matrix_t* LU, Matrix_t* P);


/**
 * @brief Solve linear equation with LU decomposition components
 * 
 * @param P pointer to input Matrix object: Set a pivot selection matrix of P-LU decomposition
 * @param LU pointer to input Matrix object: Set a marged-LU matrix of P-LU decomposition
 * @param v pointer to input Matrix object: Set a right-hand side vector of simultaneous linear equations
 * @param destination pointer to output Matrix object: It will be assigned the solution vector
 * @return int error code : returns 0 if executed successfully
 * @sa LU_decomposition
 */
int SolveFromPLU(const Matrix_t* P, const Matrix_t* LU, const Matrix_t* v, Matrix_t* destination);

/**
 * @brief Calculate determinant of input matrix with LU decomposition
 * 
 * @param matrix pointer to input Matrix
 * @return float returns float value of determinant
 */
float Det(const Matrix_t* matrix);

/**
 * @brief Operate diagonalization with Jacobi method (D ~= P^-1 * source * P)
 * 
 * @param source pointer to input Matrix
 * @param D pointer to output Matrix that will be assigned diagonalized matrix
 * @param P pointer to output Matrix that will be assigned eigenvector matrix (especially, P^t = P ^-1, det(P) = 1)
 * @return int error code : returns 0 if executed successfully
 */
int DiagonalizeMatrix(const Matrix_t* source, Matrix_t* D, Matrix_t* P);

/**
 * @brief Output Matrix elements to stdout
 * 
 * @param matrix pointer to input Matrix
 */
void PrintMatrix(const Matrix_t* matrix);


#endif // _LINALG_H_

