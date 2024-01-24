#include "linalg.h"

// Function to create a new Matrix_t
Matrix_t CreateMatrix(unsigned int rows, unsigned int cols)
{
    Matrix_t matrix;

    matrix.rows = rows;
    matrix.cols = cols;

    if (rows == 0 || cols == 0)
    {
        // Handle error: rows or cols is 0, cannot allocate memory
        matrix.data = NULL;
        matrix.rows = 0;
        matrix.cols = 0;

        return matrix;
    }

    matrix.data = (float*)malloc(rows * cols * sizeof(float));

    if (!matrix.data)
    {
        matrix.rows = 0;
        matrix.cols = 0;
    }

    return matrix;
}

// Function to create a same matrix
Matrix_t CopyMatrix(const Matrix_t* source)
{
    Matrix_t matrix;
    const unsigned int n = source->rows * source->cols;

    matrix.rows = source->rows;
    matrix.cols = source->cols;

    if (!n)
    {
        // Handle error: rows or cols is 0, cannot allocate memory
        matrix.data = NULL;
        matrix.rows = 0;
        matrix.cols = 0;

        return matrix;
    }

    matrix.data = (float*)malloc(n * sizeof(float));
    

    if (matrix.data)
    {
        unsigned int i;
        for (i = 0; i < n; ++i)
        {
            matrix.data[i] = source->data[i];
        }
    }
    else
    {
        matrix.rows = 0;
        matrix.cols = 0;
    }

    return matrix;
}

// Function to create a new Matrix_t and initialize with a value
Matrix_t CreateMatrixWithInitValue(unsigned int rows, unsigned int cols, float init_value)
{
    Matrix_t matrix;
    const unsigned int n = rows * cols;

    matrix.rows = rows;
    matrix.cols = cols;

    if (!n)
    {
        // Handle error: rows or cols is 0, cannot allocate memory
        matrix.data = NULL;
        matrix.rows = 0;
        matrix.cols = 0;

        return matrix;
    }

    matrix.data = (float*)malloc(n * sizeof(float));

    if (matrix.data)
    {
        unsigned int i;
        for (i = 0; i < n; ++i)
        {
            matrix.data[i] = init_value;
        }
    }
    else
    {
        matrix.rows = 0;
        matrix.cols = 0;
    }

    return matrix;
}

// Function to create a identity matrix
Matrix_t CreateIdentityMatrix(unsigned int size)
{
    Matrix_t matrix;
    const unsigned int n = size * size;

    matrix.rows = size;
    matrix.cols = size;

    if (size == 0)
    {
        // Handle error: rows or cols is 0, cannot allocate memory
        matrix.data = NULL;
        return matrix;
    }

    matrix.data = (float*)malloc(n * sizeof(float));

    if (matrix.data)
    {
        unsigned int i;
        for (i = 0; i < n; ++i)
        {
            matrix.data[i] = i%size == i/matrix.cols ? 1.0f : 0.0f;
        }
    }
    else
    {
        matrix.rows = 0;
        matrix.cols = 0;
    }

    return matrix;
}

int CreateMatrixArray(Matrix_t** pp_mat, unsigned int rows, unsigned int cols, size_t array_size)
{
    int e = 0;
    *pp_mat = (Matrix_t*)malloc(array_size * sizeof(Matrix_t));
    
    if (*pp_mat)
    {
        size_t i;
        for (i = 0; i < array_size; ++i)
        {
            (*pp_mat)[i] = CreateMatrix(rows, cols);

            e += (*pp_mat) + i ? 0 : -1;
        }
    }
    else
    {
        return -1;
    }

    if (e)
    {
        FreeMatrixArray(pp_mat, array_size);
        return -1;
    }
    else
    {
        return 0;
    }

}


// Function to assign values to a matrix.
void AssignMatrix(Matrix_t* matrix, float* _array, unsigned int size)
{
    unsigned int i;
    for (i = 0; i < size; ++i)
    {
        matrix->data[i] = _array[i];
    }
}

// Function to assign a float values to a matrix.
void AssignMatrixWithValue(Matrix_t* matrix, float value)
{
    unsigned int i;
    for (i = 0; i < matrix->rows*matrix->cols; ++i)
    {
        matrix->data[i] = value;
    }
}

// Function to assign a float values to a matrix.
void AssignIdentityMatrix(Matrix_t* matrix)
{
    unsigned int n = matrix->rows * matrix->cols;
    unsigned int i;
    for (i = 0; i < n; ++i)
    {
        matrix->data[i] = i % matrix->rows == i / matrix->cols ? 1.0f : 0.0f;
    }

}

int GetRowVector(const Matrix_t* source, Matrix_t* destination, unsigned int n_cols)
{
    if (source->cols > n_cols && source->rows == destination->rows)
    {
        unsigned int i;
        for (i = 0; i < source->rows; ++i)
        {
            Mat_at(*destination, i, 0) = Mat_at(*source, i, n_cols);
        }

        return 0;
    }
    else
    {
        return -1;
    }

}

int GetColVector(const Matrix_t* source, Matrix_t* destination, unsigned int n_rows)
{
    if (source->rows > n_rows && source->cols == destination->cols)
    {
        unsigned int i;
        for (i = 0; i < source->rows; ++i)
        {
            Mat_at(*destination, 0, i) = Mat_at(*source, n_rows, i);
        }

        return 0;
    }
    else
    {
        return -1;
    }
}

// Function to free the memory of a Matrix_t
void FreeMatrix(Matrix_t* matrix)
{
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

void FreeMatrixArray(Matrix_t** mat_array, size_t array_size)
{
    Matrix_t* p_temp = *mat_array;
    size_t i;
    for (i = 0; i < array_size; ++i)
    {
        FreeMatrix(&p_temp[i]);
    }

    free(*mat_array);
    *mat_array = NULL;
}

// Function to add two matrices (result = a + b)
int AddMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result)
{
    if (a->rows != b->rows || a->cols != b->cols || result->rows < a->rows || result->cols < a->cols || !a->data || !b->data || !result->data )
    {
        // Handle error: Matrices must have the same dimensions for addition and output matrix should not be the same as input matrices
        return -1;
    }

    unsigned int i, j;
    for (i = 0; i < a->rows; ++i)
    {
        for (j = 0; j < a->cols; ++j)
        {
            Mat_at(*result, i, j) = Mat_at(*a, i, j) + Mat_at(*b, i, j);
        }
    }

    return 0;
}


// Function to subtract two matrices (result = a - b)
int SubtractMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result)
{
    if (a->rows != b->rows || a->cols != b->cols || result->rows < a->rows || result->cols < a->cols || !a->data || !b->data || !result->data)
    {
        // Handle error: Matrices must have the same dimensions for substraction and output matrix should not be the same as input matrices
        return -1;
    }

    unsigned int i, j;
    for (i = 0; i < a->rows; ++i)
    {
        for (j = 0; j < a->cols; ++j)
        {
            Mat_at(*result, i, j) = Mat_at(*a, i, j) - Mat_at(*b, i, j);
        }
    }

    return 0;
}

// Function to subtract two matrices (target += weight * input)
int AppendMatrixWithWeight(Matrix_t* target, const Matrix_t* input, float weight)
{
    if (target->rows != input->rows || target->cols != input->cols || !target->data || !input->data)
    {
        // Handle error: Matrices must have the same dimensions for addition and output matrix should not be the same as input matrices
        return -1;
    }

    unsigned int i, j;
    for (i = 0; i < input->rows; ++i)
    {
        for (j = 0; j < input->cols; ++j)
        {
            Mat_at(*target, i, j) += weight*Mat_at(*input, i, j);
        }
    }

    return 0;
}

// Function to transpose a matrix
int TransposeMatrix(const Matrix_t* input, Matrix_t* output)
{
    if (!input->data || !output->data || output->cols * output->rows < input->cols * input->rows)
    {
        // Handle error: Invalid matrix pointers and output matrix should not be the same as input matrix, or output matrix doesn't have enough size memory.
        return -1;
    }
    
    if (input == output)
    {
        float temp;
        unsigned int i, j;

        output->rows = input->cols;
        output->cols = input->rows;

        
        for (i = 0; i < input->rows; ++i)
        {
            for (j = i+1; j < input->cols; ++j)
            {
                temp = Mat_at(*output, j, i);
                Mat_at(*output, j, i) = Mat_at(*input, i, j);
                Mat_at(*input, i, j) = temp;
            }
        }
    }
    else
    {
        unsigned int i, j;

        output->rows = input->cols;
        output->cols = input->rows;

        for (i = 0; i < input->rows; ++i)
        {
            for (j = 0; j < input->cols; ++j)
            {
                Mat_at(*output, j, i) = Mat_at(*input, i, j);
            }
        }

    }


    return 0;
}

// Function to multiply matrix a and b (result = a*b).
int MultiplyMatrices(const Matrix_t* a, const Matrix_t* b, Matrix_t* result)
{
    if (a->cols != b->rows || result->rows != a->rows || result->cols != b->cols || result->rows*result->cols < a->rows * b->cols || !a->data || !b->data || !result->data )
    {
        // Handle error: Invalid matrix dimensions for multiplication
        return -1;
    }

    if (a->data == result->data || b->data == result->data) 
    {
        // Use heap-allocated temporary matrix
        Matrix_t tempMatrix;
        unsigned int temp_size = a->rows * b->cols;
        unsigned int i, j, k;

        #ifdef USE_STATIC_ARRAY 
        float temp_array[LINALG_MATRIX_BUFFER_SIZE];


        if (LINALG_MATRIX_BUFFER_SIZE < temp_size) // when required size is larger than LINALG_MATRIX_BUFFER_SIZE
        {
            return -1;
        }

        tempMatrix.data = temp_array;

        
        #else

        tempMatrix.data = (float*)STACK_ALLOC(a->rows * b->cols * sizeof(float));

        #endif //!USE_STATIC_ARRAY

        

        if (!tempMatrix.data) 
        {
            // Handle error: Memory allocation failed
            return -1;
        }

        tempMatrix.rows = a->rows;
        tempMatrix.cols = b->cols;

        for (i = 0; i < a->rows; ++i)
        {
            for (j = 0; j < b->cols; ++j)
            {
                Mat_at(tempMatrix, i, j) = 0;
                for (k = 0; k < a->cols; ++k)
                {
                    Mat_at(tempMatrix, i, j) += Mat_at(*a, i, k) * Mat_at(*b, k, j);
                }
            }
        }

        for (i = 0; i < temp_size; ++i)
        {
            result->data[i] = tempMatrix.data[i];
        }
        
    }
    else
    {
        unsigned int i, j, k;
        // Use the result matrix directly
        result->rows = a->rows;
        result->cols = b->cols;
        

        for (i = 0; i < a->rows; ++i)
        {
            for (j = 0; j < b->cols; ++j)
            {
                Mat_at(*result, i, j) = 0;
                for (k = 0; k < a->cols; ++k)
                {
                    Mat_at(*result, i, j) += Mat_at(*a, i, k) * Mat_at(*b, k, j);
                }
            }
        }

    }

    return 0;  // Return 0 for success
}

// Function to multiply a matrix by a scalar
int MultiplyMatrixByScalar(const Matrix_t* input, float scalar, Matrix_t* output)
{
    if (!input->data || !output->data || output->rows * output->cols < input->rows * input->cols)
    {
        // Handle error: Invalid matrix pointers
        return -1;
    }

    unsigned int i, j;
    for (i = 0; i < input->rows; ++i)
    {
        for (j = 0; j < input->cols; ++j)
        {
            Mat_at(*output, i, j) = scalar * Mat_at(*input, i, j);
        }
    }

    return 0;  // Return 0 for success
}

int InvertMatrix(const Matrix_t* input, Matrix_t* output)
{
    if (input->cols != input->rows || output->cols < input->rows || output->cols < input->cols || !input->data || !output->data)
    {
        // Handle error: Invalid matrix dimensions for Invertion, or not allocated.
        return -1;
    }

    unsigned int n = input->rows;

    output->rows = input->rows;
    output->cols = input->rows;

    

    if (n == 0)
    {
        return -1;
    }
    else if (n == 1)
    {
        if (Mat_at(*input, 0, 0) == 0.0f)
        {
            return -1;
        }

        Mat_at(*output, 0, 0) = 1.0f / Mat_at(*input, 0, 0);
        return 0;
    }
    //2~2 Inverting pattern
    else if (n == 2)
    {
        float temp;
        float det = Mat_at(*input, 0, 0) * Mat_at(*input, 1, 1) - Mat_at(*input, 0, 1) * Mat_at(*input, 1, 0);

        if (det == 0.0f)
        {
            // Handle error: Matrix is singular, inverse does not exist
            return -1;
        }

        temp = Mat_at(*input, 0, 0);
        Mat_at(*output, 0, 0) = Mat_at(*input, 1, 1) * (1.0f/det);
        Mat_at(*output, 0, 1) = -Mat_at(*input, 0, 1) * (1.0f / det);
        Mat_at(*output, 1, 0) = -Mat_at(*input, 1, 0) * (1.0f / det);
        Mat_at(*output, 1, 1) = temp * (1.0f / det);

        return 0;
    }
    //3~3 Inverting pattern
    else if (n == 3)
    {

        float det = Mat_at(*input, 0, 0) * (Mat_at(*input, 1, 1) * Mat_at(*input, 2, 2) - Mat_at(*input, 1, 2) * Mat_at(*input, 2, 1)) -
            Mat_at(*input, 0, 1) * (Mat_at(*input, 1, 0) * Mat_at(*input, 2, 2) - Mat_at(*input, 1, 2) * Mat_at(*input, 2, 0)) +
            Mat_at(*input, 0, 2) * (Mat_at(*input, 1, 0) * Mat_at(*input, 2, 1) - Mat_at(*input, 1, 1) * Mat_at(*input, 2, 0));

        if (det == 0.0f)
        {
            return -1;
        }


        if (input->data == output->data)
        {
            float tmp[5];

            tmp[0] = Mat_at(*output, 0, 0);
            Mat_at(*output, 0, 0) = (Mat_at(*input, 1, 1) * Mat_at(*input, 2, 2) - Mat_at(*input, 1, 2) * Mat_at(*input, 2, 1)) * (1.0f / det);
            tmp[1] = Mat_at(*output, 0, 1);
            Mat_at(*output, 0, 1) = (Mat_at(*input, 0, 2) * Mat_at(*input, 2, 1) - Mat_at(*input, 0, 1) * Mat_at(*input, 2, 2)) * (1.0f / det);
            tmp[2] = Mat_at(*output, 0, 2);
            Mat_at(*output, 0, 2) = (tmp[1] * Mat_at(*input, 1, 2) - Mat_at(*input, 0, 2) * Mat_at(*input, 1, 1)) * (1.0f / det);

            tmp[3] = Mat_at(*output, 1, 0);
            Mat_at(*output, 1, 0) = (Mat_at(*input, 1, 2) * Mat_at(*input, 2, 0) - Mat_at(*input, 1, 0) * Mat_at(*input, 2, 2)) * (1.0f / det);
            tmp[4] = Mat_at(*output, 1, 1);
            Mat_at(*output, 1, 1) = (tmp[0] * Mat_at(*input, 2, 2) - tmp[2] * Mat_at(*input, 2, 0)) * (1.0f / det);
            Mat_at(*output, 1, 2) = (tmp[2] * tmp[3] - tmp[0] * Mat_at(*input, 1, 2)) * (1.0f / det);

            tmp[2] = Mat_at(*output, 2, 0);
            Mat_at(*output, 2, 0) = (tmp[3] * Mat_at(*input, 2, 1) - tmp[4] * Mat_at(*input, 2, 0)) * (1.0f / det);
            Mat_at(*output, 2, 1) = (tmp[1] * tmp[2] - tmp[0] * Mat_at(*input, 2, 1)) * (1.0f / det);
            Mat_at(*output, 2, 2) = (tmp[0] * tmp[4] - tmp[1] * tmp[3]) * (1.0f / det);

        }
        else
        {
            Mat_at(*output, 0, 0) = (Mat_at(*input, 1, 1) * Mat_at(*input, 2, 2) - Mat_at(*input, 1, 2) * Mat_at(*input, 2, 1)) * (1.0f / det);
            Mat_at(*output, 0, 1) = (Mat_at(*input, 0, 2) * Mat_at(*input, 2, 1) - Mat_at(*input, 0, 1) * Mat_at(*input, 2, 2)) * (1.0f / det);
            Mat_at(*output, 0, 2) = (Mat_at(*input, 0, 1) * Mat_at(*input, 1, 2) - Mat_at(*input, 0, 2) * Mat_at(*input, 1, 1)) * (1.0f / det);

            Mat_at(*output, 1, 0) = (Mat_at(*input, 1, 2) * Mat_at(*input, 2, 0) - Mat_at(*input, 1, 0) * Mat_at(*input, 2, 2)) * (1.0f / det);
            Mat_at(*output, 1, 1) = (Mat_at(*input, 0, 0) * Mat_at(*input, 2, 2) - Mat_at(*input, 0, 2) * Mat_at(*input, 2, 0)) * (1.0f / det);
            Mat_at(*output, 1, 2) = (Mat_at(*input, 0, 2) * Mat_at(*input, 1, 0) - Mat_at(*input, 0, 0) * Mat_at(*input, 1, 2)) * (1.0f / det);

            Mat_at(*output, 2, 0) = (Mat_at(*input, 1, 0) * Mat_at(*input, 2, 1) - Mat_at(*input, 1, 1) * Mat_at(*input, 2, 0)) * (1.0f / det);
            Mat_at(*output, 2, 1) = (Mat_at(*input, 0, 1) * Mat_at(*input, 2, 0) - Mat_at(*input, 0, 0) * Mat_at(*input, 2, 1)) * (1.0f / det);
            Mat_at(*output, 2, 2) = (Mat_at(*input, 0, 0) * Mat_at(*input, 1, 1) - Mat_at(*input, 0, 1) * Mat_at(*input, 1, 0)) * (1.0f / det);
        }

        return 0;
    }
    // general high order inverting pattern 
    else
    {

        Matrix_t augmented;
        unsigned int temp_size = n * 2 * n;

        #ifdef USE_STATIC_ARRAY

        float temp_array[2*LINALG_MATRIX_BUFFER_SIZE];

        if (2*LINALG_MATRIX_BUFFER_SIZE < temp_size) // when required size is larger than LINALG_MATRIX_BUFFER_SIZE
        {
            return -1;
        }

        augmented.data = temp_array;

        #else

        augmented.data = (float*)STACK_ALLOC(temp_size * sizeof(float));

        #endif // USE_STATIC_ARRAY


        if (!augmented.data)
        {
            //error: failed to allocate stack memory
            return -1;
        }

        augmented.rows = n;
        augmented.cols = 2 * n;

        unsigned int pivotRow;
        float pivotValue, factor, temp;
        
        unsigned int i, j, k;

        // init augmented matrix
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                Mat_at(augmented, i, j) = Mat_at(*input, i, j);
                Mat_at(augmented, i, j + n) = (i == j) ? 1.0f : 0.0f;
            }
        }

        // exec. gauss-jordan method
        for (i = 0; i < n; ++i)
        {

            // select a pivot
            pivotRow = i;
            for (k = i + 1; k < n; ++k)
            {
                if (fabsf(Mat_at(augmented, k, i)) > fabsf(Mat_at(augmented, pivotRow, i)))
                {
                    pivotRow = k;
                }
            }

            // swap a row of pivot
            if (pivotRow != i)
            {
                for (j = 0; j < 2 * n; ++j)
                {
                    temp = Mat_at(augmented, i, j);
                    Mat_at(augmented, i, j) = Mat_at(augmented, pivotRow, j);
                    Mat_at(augmented, pivotRow, j) = temp;
                }
            }

            // assign pivot to 1
            pivotValue = Mat_at(augmented, i, i);

            if (pivotValue == 0.0f)
            {
                //error: the matrix is not regular.
                return -1;
            }

            for (j = 0; j < 2 * n; ++j)
            {
                Mat_at(augmented, i, j) /= pivotValue;
            }

            // subs col of pivot to 0
            for (k = 0; k < n; ++k)
            {
                if (k != i) {
                    factor = Mat_at(augmented, k, i);
                    for (j = 0; j < 2 * n; ++j)
                    {
                        Mat_at(augmented, k, j) -= factor * Mat_at(augmented, i, j);
                    }
                }
            }
        }

        // copy the result
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                Mat_at(*output, i, j) = Mat_at(augmented, i, j + n);
            }
        }

        return 0;
    }
}

int CholeskyDecomposition(const Matrix_t* source, Matrix_t* destination)
{
    if (!source->data || source->rows != source->cols || destination->rows != destination->cols || source->rows != destination->rows)
    {
        // Error: not a symmetric matrix
        return -1;
    }

    unsigned int n = source->rows;
    unsigned int i, j, k;
    int e = 0;
    float sum;

    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {

            if (i >= j)
            {
                sum = Mat_at(*source, i, j);

                for (k = 0; k < j; ++k)
                {
                    sum -= Mat_at(*destination, i, k) * Mat_at(*destination, j, k);
                }

                if (sum < 0.0f)
                {
                    sum = fabsf(sum);
                    e = -1;
                }

                if (i == j)
                {
                    Mat_at(*destination, i, j) = sqrtf(sum);
                }
                else
                {
                    Mat_at(*destination, i, j) = Mat_at(*destination, j, j) != 0.0f ? sum / Mat_at(*destination, j, j) : fabsf(sum) > 1.0f ? (sum >= 0.0f ? 1.0f : -1.0f) * FLT_MAX : sum * FLT_MAX;
                }
            }
            else
            {
                Mat_at(*destination, i, j) = 0.0f;
            }

        }
    }

    return e;
}


// function to decomposition from source matrix to LU-marged matrix and Pivot matrix
int LU_Decomposition(const Matrix_t* source, Matrix_t* LU, Matrix_t* P)
{
    // size and allocation check 
    if (!source->data || !LU->data || !P->data || source->rows != source->cols || LU->rows != LU->cols || P->rows != P->cols || source->rows != LU->rows || source->rows != P->rows || LU->data == P->data)
    {
        return -1;
    }

    unsigned int n = source->rows;
    unsigned int i, j, k, max_pivot_index;
    float max_pivot;

    // initilize LU and P
    if (source->data != P->data)
    {
        for (i = 0; i < n * n; ++i)
        {
            P->data[i] = i % source->rows == i / source->cols ? 1.0f : 0.0f;
            LU->data[i] = source->data[i];
        }

        for (i = 0; i < n; ++i)
        {
            // Select pivot
            max_pivot = fabsf(Mat_at(*LU, i, i));
            max_pivot_index = i;
            for (j = i + 1; j < n; ++j)
            {
                if (fabsf(Mat_at(*LU, j, i)) > max_pivot)
                {
                    max_pivot = fabsf(Mat_at(*LU, j, i));
                    max_pivot_index = j;
                }
            }

            if (max_pivot == 0.0f)
            {
                return -1; // failed to decomposition
            }

            if (max_pivot_index != i)
            {
                float tmp;
                // swap cols
                for (k = 0; k < n; ++k)
                {
                    tmp = Mat_at(*LU, i, k);
                    Mat_at(*LU, i, k) = Mat_at(*LU, max_pivot_index, k);
                    Mat_at(*LU, max_pivot_index, k) = tmp;

                    tmp = Mat_at(*P, i, k);
                    Mat_at(*P, i, k) = Mat_at(*P, max_pivot_index, k);
                    Mat_at(*P, max_pivot_index, k) = tmp;
                }
            }

            // LU decomp
            for (j = i + 1; j < n; ++j)
            {
                Mat_at(*LU, j, i) /= Mat_at(*LU, i, i);
                for (k = i + 1; k < n; ++k)
                {
                    Mat_at(*LU, j, k) -= Mat_at(*LU, j, i) * Mat_at(*LU, i, k);
                }
            }
        }
    }
    else // if source == P
    {
        Matrix_t P_result;
#ifdef USE_STATIC_ARRAY

        float P_result_array[LINALG_MATRIX_BUFFER_SIZE];
        P_result.data = P_result_array;

        if (n*n <= LINALG_MATRIX_BUFFER_SIZE) // when the size is less than LINALG_MATRIX_BUFFER_SIZE
        {
            P_result.data = P_result_array;
        }
        else
        {
            P_result.data = NULL;
        }

#else

        P_result.data = STACK_ALLOC(source->rows * source->cols * sizeof(float));

#endif // USE_STATIC_ARRAY


        if (!P_result.data)
        {
            return -1; // failed to allocation
        }

        P_result.rows = P_result.cols = source->rows;

        for (i = 0; i < n * n; ++i)
        {
            P_result.data[i] = i % source->rows == i / source->cols ? 1.0f : 0.0f;
            LU->data[i] = source->data[i];
        }

        for (i = 0; i < n; ++i)
        {
            // Select pivot
            max_pivot = fabsf(Mat_at(*LU, i, i));
            max_pivot_index = i;
            for (j = i + 1; j < n; ++j)
            {
                if (fabsf(Mat_at(*LU, j, i)) > max_pivot)
                {
                    max_pivot = fabsf(Mat_at(*LU, j, i));
                    max_pivot_index = j;
                }
            }

            if (max_pivot == 0.0f)
            {
                return -1; // failed to decomposition
            }

            if (max_pivot_index != i)
            {
                float tmp;
                // swap cols
                for (k = 0; k < n; ++k)
                {
                    tmp = Mat_at(*LU, i, k);
                    Mat_at(*LU, i, k) = Mat_at(*LU, max_pivot_index, k);
                    Mat_at(*LU, max_pivot_index, k) = tmp;

                    tmp = Mat_at(P_result, i, k);
                    Mat_at(P_result, i, k) = Mat_at(P_result, max_pivot_index, k);
                    Mat_at(P_result, max_pivot_index, k) = tmp;
                }
            }

            // LU decomp
            for (j = i + 1; j < n; ++j)
            {
                Mat_at(*LU, j, i) /= Mat_at(*LU, i, i);
                for (k = i + 1; k < n; ++k)
                {
                    Mat_at(*LU, j, k) -= Mat_at(*LU, j, i) * Mat_at(*LU, i, k);
                }
            }
        }

        // copy from buffer to output
        for (i = 0; i < n * n; ++i)
        {
            P->data[i] = P_result.data[i];
        }

    }
    

    return 0;
}


// Function to solve linear equation with LU decomposition components
int SolveFromPLU(const Matrix_t* P, const Matrix_t* LU, const Matrix_t* v, Matrix_t* destination)
{
    if (!P->data || !LU->data || !v->data || !destination->data || LU->rows != LU->cols || P->rows != P->cols || LU->rows != P->rows || v->rows != LU->rows || destination->rows != LU->rows || v->cols != 1 || destination->cols != 1)
    {
        return -1; // check
    }

    Matrix_t y;
    unsigned int n = LU->rows;
    unsigned int i, j;
    float sum;

#ifdef USE_STATIC_ARRAY

    float y_array[LINALG_VECTOR_BUFFER_SIZE];

    y.data = y_array;

    if (LINALG_VECTOR_BUFFER_SIZE < n) // when required size is larger than LINALG_MATRIX_BUFFER_SIZE
    {
        return -1; // failed to allocate
    }

#else

    y.data = STACK_ALLOC(n * sizeof(float));


#endif // USE_STATIC_ARRAY

    if (!y.data)
    {
        return -1; // failed to allocate
    }

    y.rows = n;
    y.cols = 1;

    // y = Pv
    MultiplyMatrices(P, v, &y);

    // forward-assignmentFLy = Pv
    for (i = 0; i < n; ++i)
    {
        sum = y.data[i];
        for (j = 0; j < i; ++j)
        {
            sum -= Mat_at(*LU, i, j) * y.data[j];
        }
        y.data[i] = sum; // y[i] = Pv[i] - ƒ°Ly[0~i-1]
    }

    // backward-assignmentFUx = y
    for (i = n - 1; i != (unsigned int)(-1); --i)
    {
        sum = y.data[i];
        for (j = n - 1 ; j > i ; --j)
        {
            sum -= Mat_at(*LU, i, j) * Mat_at(*destination, j, 0);
        }
        Mat_at(*destination, i, 0) = sum / Mat_at(*LU, i, i);
    }

    return 0;
}


// function to calculate determinant of the matrix
float Det(const Matrix_t* matrix)
{
    if (!matrix->data || matrix->rows != matrix->cols)
    {
        return -1; // Error condition, non-square matrix
    }

    unsigned int n = matrix->rows, i, j, k, max_pivot_index;
    float det = 1.0f;
    float max_pivot;

#ifdef USE_STATIC_ARRAY

    if (LINALG_MATRIX_BUFFER_SIZE < n*n)
    {
        return -1;
    }


    float lu[LINALG_MATRIX_BUFFER_SIZE];

#else

    float* lu = (float*)STACK_ALLOC(n * n * sizeof(float));

#endif // USE_STATIC_ARRAY


    if (!lu)
    {
        return -1; // Error condition, memory allocation failed
    }

    // Copy the input matrix to LU matrix
    for (i = 0; i < n * n; i++)
    {
        lu[i] = matrix->data[i];
    }

    // Perform LU Decomposition with partial pivoting
    for (i = 0; i < n; i++)
    {
        // Partial pivoting
        max_pivot_index = i;
        max_pivot = fabsf(lu[i * n + i]);
        for (k = i + 1; k < n; k++)
        {
            if (fabsf(lu[k * n + i]) > max_pivot)
            {
                max_pivot = fabsf(lu[k * n + i]);
                max_pivot_index = k;
            }
        }

        // Swap rows if necessary
        if (max_pivot_index != i)
        {
            for (k = 0; k < n; k++)
            {
                float temp = lu[i * n + k];
                lu[i * n + k] = lu[max_pivot_index * n + k];
                lu[max_pivot_index * n + k] = temp;
            }
            det *= -1; // Swap changes the sign of determinant
        }

        if (lu[i * n + i] == 0.0f)
        {
            return 0.0f; // Singular matrix
        }

        // Decomposition
        for (k = i + 1; k < n; k++)
        {
            lu[k * n + i] /= lu[i * n + i];
            for (j = i + 1; j < n; j++)
            {
                lu[k * n + j] -= lu[k * n + i] * lu[i * n + j];
            }
        }
    }

    // Calculate determinant as the product of diagonal elements
    for (i = 0; i < n; i++)
    {
        det *= lu[i * n + i];
    }

    return det;
}

// Function to perform the Jacobi eigenvalue algorithm
int DiagonalizeMatrix(const Matrix_t* source, Matrix_t* D, Matrix_t* P)
{
    if (!source->data || !D->data || !P->data || source->rows != source->cols || D->rows != D->cols || P->rows != P->cols || source->rows != D->rows || source->rows != P->rows)
    {
        return -1;
    }

    unsigned int n = source->rows;
    float maxOffDiag, t, c, s, tau, temp_p, temp_q;
    unsigned int p, q, i, j, iter;

    // Initialize D and V
    for (i = 0; i < n*n; i++)
    {
            P->data[i] = i % source->rows == i / source->cols ? 1.0f : 0.0f;
            D->data[i] = source->data[i];
    }

    iter = 0;
    do {
        // Find largest off-diagonal absolute value in D
        maxOffDiag = 0.0f;
        for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                if (fabs(Mat_at(*D, i, j)) > maxOffDiag)
                {
                    maxOffDiag = fabsf(Mat_at(*D, i, j));
                    p = i;
                    q = j;
                }
            }
        }

        // Jacobi rotation
        if (maxOffDiag > LINALG_MAXOFFDIAG_EPS)
        {
            tau = (Mat_at(*D, q, q) - Mat_at(*D, p, p)) / (2.0f * Mat_at(*D, p, q));
            t = (tau > 0) ? 1.0f / (tau + sqrtf(1.0f + tau * tau)) : 1.0f / (tau - sqrtf(1.0f + tau * tau));
            c = 1.0f / sqrtf(1.0f + t * t);
            s = t * c;

            for (i = 0; i < n; i++)
            {
                if (i != p && i != q)
                {
                    temp_p = Mat_at(*D, i, p);
                    Mat_at(*D, i, p) = c * temp_p - s * Mat_at(*D, i, q);
                    Mat_at(*D, i, q) = s * temp_p + c * Mat_at(*D, i, q);
                    Mat_at(*D, p, i) = Mat_at(*D, i, p);  // Symmetric matrix
                    Mat_at(*D, q, i) = Mat_at(*D, i, q);  // Symmetric matrix
                }
            }

            temp_p = Mat_at(*D, p, p);
            Mat_at(*D, p, p) = c * c * temp_p + s * s * Mat_at(*D, q, q) - 2.0f * s * c * Mat_at(*D, p, q);
            Mat_at(*D, q, q) = s * s * temp_p + c * c * Mat_at(*D, q, q) + 2.0f * s * c * Mat_at(*D, p, q);
            Mat_at(*D, p, q) = 0.0f;
            Mat_at(*D, q, p) = 0.0f;

            for (i = 0; i < n; i++)
            {
                temp_p = Mat_at(*P, i, p);
                temp_q = Mat_at(*P, i, q);
                Mat_at(*P, i, p) = c * temp_p - s * temp_q;
                Mat_at(*P, i, q) = s * temp_p + c * temp_q;
            }
        }

        ++iter;
    } while (maxOffDiag > LINALG_MAXOFFDIAG_EPS && iter < LINALG_MAX_ITER);

    if (iter >= LINALG_MAX_ITER)
    {
        return -1;
    }

    return 0;
}


void PrintMatrix(const Matrix_t* matrix)
{
    unsigned int i, j;

    for (i = 0; i < matrix->rows; ++i)
    {
        for (j = 0; j < matrix->cols; ++j)
        {
            printf("%e ", Mat_at(*matrix, i, j));
        }
        putchar('\n');
    }
}


