#include "ssmodel.h"



int CreateSSmodel(StateSpace_t* model, Matrix_t* A, Matrix_t* B, Matrix_t* C, Matrix_t* D)
{

    if (!A || A->rows != A->cols)
    {
        //Handle error: A must be square matrix
        return -1;
    }


    if (B)
    {
        if (A->rows != B->rows )
        {
            //Handle error: rows of A and rows of B must be the same
            return -1;
        }

        if (D) //when D is given
        {
            if (B->cols != D->cols)
            {
                //Handle error: cols of B and cols of D must be the same
                return -1;
            }
        }
    }

    if (C) //When C is given
    {
        if (A->cols != C->cols)
        {
            //Handle error: cols of A and colss of C must be the same
            return -1;
        }

        if (D) //when D is given
        {
            if (C->rows != D->rows)
            {
                //Handle error: rows of C and rows of D must be the same
                return -1;
            }
        }
    }

    model->A = A ? *A : CreateMatrix(0, 0);  // Check if A is NULL
    model->B = B ? *B : CreateMatrix(0, 0);  // Check if B is NULL
    model->C = C ? *C : CreateMatrix(0, 0);  // Check if C is NULL
    model->D = D ? *D : CreateMatrix(0, 0);  // Check if D is NULL

	model->x = CreateMatrixWithInitValue(model->A.rows, 1, 0.0f);
    model->y = CreateMatrixWithInitValue(model->C.rows, 1, 0.0f);
    model->u = CreateMatrixWithInitValue(model->B.cols, 1, 0.0f);


    if (model->x.data && model->u.data && !(C && !model->y.data))
    {
        return 0;
    }
    else
    {
        FreeMatrix(&model->x);
        FreeMatrix(&model->y);
        FreeMatrix(&model->u);
        return -1;
    }
}


int SSmodelStep(StateSpace_t* model, Matrix_t* input)
{
    Matrix_t tempAx;
    unsigned int temp_size = model->x.rows * model->x.cols;
    unsigned int i;
    
    tempAx.rows = model->x.rows;
    tempAx.cols = model->x.cols;

    // Allocate memory for temp matrix
    #ifdef USE_STATIC_ARRAY
    float temp_array[LINALG_VECTOR_BUFFER_SIZE];
    if (temp_size <= LINALG_VECTOR_BUFFER_SIZE) // when the size is less than 8
    {
        tempAx.data = temp_array;
    }
    else  // This system should not be used to calculate matrices of learger size.
    {
        tempAx.data = NULL;
    }



    #else

    tempAx.data = (float*)STACK_ALLOC(temp_size * sizeof(float));

    #endif // USE_STATIC_ARRAY

    if (!tempAx.data)
    {
        //Hanlde error: Stack memory allocation failed
        return -1;
    }

    // Perform the calculations using model's matrices
    if (model->A.data)
    {
        MultiplyMatrices(&model->A, &model->x, &tempAx);
        for (i = 0; i < tempAx.rows * tempAx.cols; ++i)
        {
            model->x.data[i] = tempAx.data[i];
        }
    }
    else
    {
        //Hanlde error: failed to calculation, A must be set values
        return -1;
    }

    if (model->B.data && input)
    {
        MultiplyMatrices(&model->B, input, &model->u); // model->u is used as a temporary buffer
        for (i = 0; i < model->u.rows * model->u.cols; ++i)
        {
            model->x.data[i] += model->u.data[i]; // append (B*input) to state variables
            model->u.data[i] = input->data[i]; // set the original input values
        }
    }

    if (model->C.data)
    {
        MultiplyMatrices(&model->C, &model->x, &model->y);  // Calculate output

        if (model->D.data && input)
        {
            Matrix_t tempDu;

            #ifdef USE_STATIC_ARRAY

            float temp_array[LINALG_VECTOR_BUFFER_SIZE];
            if (model->D.rows * input->cols <= LINALG_VECTOR_BUFFER_SIZE)
            {
                tempDu.data = temp_array;
            }
            else
            {
                tempDu.data = NULL;
            }

            #else

            tempDu.data = (float*)STACK_ALLOC(model->D.rows * input->cols * sizeof(float));

            #endif // USE_STATIC_ARRAY

            
            tempDu.rows = model->D.rows;
            tempDu.cols = input->cols;

            if (!tempDu.data)
            {
                //Hanlde error: Stack memory allocation failed
                return -1;
            }

            MultiplyMatrices(&model->D, input, &tempDu);

            for (i = 0; i < model->y.rows * model->y.cols; ++i)
            {
                model->y.data[i] += tempDu.data[i];
            }

        }
    }

    return 0;
}

void FreeSSmodel(StateSpace_t* model)
{
    FreeMatrix(&model->A);
    FreeMatrix(&model->B);
    FreeMatrix(&model->C);
    FreeMatrix(&model->D);

    FreeMatrix(&model->u);
    FreeMatrix(&model->x);
    FreeMatrix(&model->y);
}

