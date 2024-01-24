#ifndef _SSMODEL_H_
#define _SSMODEL_H_

#include "linalg.h"

typedef struct
{

	Matrix_t A;
	Matrix_t B;
	Matrix_t C;
	Matrix_t D;
	Matrix_t x;
	Matrix_t y;
	Matrix_t u;


} StateSpace_t;

/**
 * @brief Initialize new digital state space object
 * 
 * @param model pointer to input StateSpace object
 * @param A pointer to state matrix (it should be n×n matrix)
 * @param B pointer to input matrix (it should be n×p matrix). If set to NULL, the input process will be ignored.(same behavior as B = O)
 * @param C pointer to output matrix (it should be q×n matrix). If set to NULL, the output process will be ignored. (same behavior as C = I)
 * @param D pointer to direct matrix (it should be q×p matrix). If set to NULL, the output process will be ignored. 
 * @return int returns 0 if executed successfully
 */
int CreateSSmodel(StateSpace_t* model, Matrix_t* A, Matrix_t* B, Matrix_t* C, Matrix_t* D);

/**
 * @brief Operate one step progression of digital state space object
 * 
 * @param model pointer to StateSpace object
 * @param input pointer to input vector (It represented by a Matrix object)
 * @return int returns 0 if executed successfully
 */
int SSmodelStep(StateSpace_t* model, Matrix_t* input);


/**
 * @brief Free allocated heap memory held StateSpace object
 * 
 * @param model pointer to StateSpace object to free
 */
void FreeSSmodel(StateSpace_t* model);


#endif // !_SSMODEL_H_
