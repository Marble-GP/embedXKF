
#ifndef _KALMAN_H_
#define _KALMAN_H_

#include <float.h>
#include <limits.h>
#include "linalg.h"
#include "ssmodel.h"

// define "USE_INCORRECT_MATRIX_SQRT_EXTENSION" when the input of CholeskyDecomposition can be non-positive definite matrix (CholeskyDecomposition is called in UKF_Operate)
#define USE_INCORRECT_MATRIX_SQRT_EXTENSION


// the initial diagonal elements of P matrix
#define P0_INIT_VAL (1e2f)

// the division rate of numerical Jacobian calculations
#define EKF_JACOBIAN_EPS ((1.0f/(float)(2U<<8)))


typedef void(*GeneralStateModelFunction_t) (Matrix_t* output, const Matrix_t* x, const Matrix_t* u, float t);

typedef struct
{
	GeneralStateModelFunction_t state_model_func;
	GeneralStateModelFunction_t state_model_jacobian_x;
	GeneralStateModelFunction_t state_model_jacobian_u;
	GeneralStateModelFunction_t output_model_func;
	GeneralStateModelFunction_t output_model_jacobian_x;
	GeneralStateModelFunction_t output_model_jacobian_u;

	StateSpace_t* linsys;

	Matrix_t P;
	Matrix_t P_adv;
	Matrix_t Q;
	Matrix_t R;
	Matrix_t G;

	Matrix_t* SP_state;
	Matrix_t* SP_output;
	float k;

	Matrix_t x_filter;

} Kalman_t;




/**
 * @brief Calculate Jacobian matrix by the input vector of x numerically (It will be used for EKF internally)
 * 
 * @param f pointer to general state model function (mathematical format of general state model function: output = f(x, u, t))
 * @param A pointer to Matrix object that will be assigned Jacobian matrix
 * @param x pointer to Matrix object that represents state vector
 * @param u pointer to Matrix object that represents input vector
 * @param t float value to represents time varying system
 * @return int error code : returns 0 if executed successfully
 * @sa _EKF_CalculateJacobian_U
 */
int _EKF_CalculateJacobian_X(GeneralStateModelFunction_t f, Matrix_t* A, const Matrix_t* x, const Matrix_t* u, float t);


/**
 * @brief Calculate Jacobian matrix by the input vector of u numerically (It will be used for EKF internally)
 * 
 * @param f pointer to general state model function (mathematical format of general state model function: output = f(x, u, t))
 * @param B pointer to Matrix object that will be assigned Jacobian matrix
 * @param x pointer to Matrix object that represents state vector
 * @param u pointer to Matrix object that represents input vector
 * @param t float value to represents time varying system
 * @return int error code : returns 0 if executed successfully
 * @sa _EKF_CalculateJacobian_X
 */
int _EKF_CalculateJacobian_U(GeneralStateModelFunction_t f, Matrix_t* B, const Matrix_t* x, const Matrix_t* u, float t);

/**
 * @brief Calculate absolute-square-root matrix using diagonalization (It will be used for UKF internally)
 * 
 * @param source pointer to input Matrix object
 * @param destination pointer to output Matrix object that will be assigned sqrt(|source|)
 * @return int error code : returns 0 if executed successfully
 */
int _UKF_AbusoluteMatrixSQRT(const Matrix_t* source, Matrix_t* destination);


/**
 * @brief Create a Kalman Filter object
 * 
 * @param kf pointer to uninitialized Kalman object
 * @param ss pointer to initialized StateSpace object
 * @param Q pointer to Matrix object that represents the system covariance matrix
 * @param R pointer to Matrix object that represents the ovservation covariance matrix
 * @return int error code : returns 0 if executed successfully
 */
int CreateKalmanFilter(Kalman_t* kf, StateSpace_t* ss, Matrix_t* Q, Matrix_t* R);


/**
 * @brief Operate one step progression of Kalman Filter
 * 
 * @param kf  pointer to Kalman object initialized as a Kalman Filter
 * @param filter_input pointer to Matrix object that reprsents the filter input (real observed value)
 * @param system_input pointer to Matrix object that reprsents the system input (real manipulated value)
 * @return int error code : returns 0 if executed successfully
 */
int KF_Operate(Kalman_t* kf, Matrix_t* filter_input, Matrix_t* system_input);

// Create an object of Extended Kalman Filter
/**
 * @brief Create a Extended Kalman Filter object
 * 
 * @param ekf  pointer to uninitialized Kalman object
 * @param Q pointer to Matrix object that represents the system covariance matrix
 * @param R pointer to Matrix object that represents the ovservation covariance matrix
 * @param f_state pointer to function that represents vector function of system model
 * @param f_output pointer to function that represents vector function of output model
 * @param f_state_jacobian_x pointer to function that represents Jacobian matrix of system model by x. If set to NULL, the Jacobian will be calculated numerically.
 * @param f_state_jacobian_u pointer to function that represents Jacobian matrix of system model by u. If set to NULL, the Jacobian will be calculated numerically.
 * @param f_output_jacobian_x pointer to function that represents Jacobian matrix of output model by x. If set to NULL, the Jacobian will be calculated numerically.
 * @param f_output_jacobian_u pointer to function that represents Jacobian matrix of output model by u. If set to NULL, the Jacobian will be calculated numerically.
 * @param state_dimension dimensions of state vector
 * @param input_dimension dimensions of input vector
 * @param output_dimension dimensions of output vector
 * @return int error code : returns 0 if executed successfully
 */
int CreateExtendedKalmanFilter(Kalman_t* ekf, Matrix_t* Q, Matrix_t* R, GeneralStateModelFunction_t f_state, GeneralStateModelFunction_t f_output, GeneralStateModelFunction_t f_state_jacobian_x, GeneralStateModelFunction_t f_state_jacobian_u, GeneralStateModelFunction_t f_output_jacobian_x, GeneralStateModelFunction_t f_output_jacobian_u, unsigned int state_dimension, unsigned int input_dimension, unsigned int output_dimension);

// Calculate a step of Extended Kalman Filter
/**
 * @brief 
 * 
 * @param ekf  pointer to Kalman object initialized as a Extended Kalman Filter
 * @param filter_input pointer to Matrix object that reprsents the filter input (real observed value)
 * @param system_input pointer to Matrix object that reprsents the system input (real manipulated value)
 * @param t float value to represents time varying system
 * @return int error code : returns 0 if executed successfully
 */
int EKF_Operate(Kalman_t* ekf, Matrix_t* filter_input, Matrix_t* system_input, float t);

// Create an object of Unscented Kalman Filter
/**
 * @brief Create a Unscented Kalman Filter object
 * 
 * @param ukf  pointer to Kalman object
 * @param Q pointer to Matrix object that represents the system covariance matrix
 * @param R pointer to Matrix object that represents the ovservation covariance matrix
 * @param f_state pointer to function that represents vector function of system model
 * @param f_output pointer to function that represents vector function of output model
 * @param state_dimension dimensions of state vector
 * @param input_dimension dimensions of input vector
 * @param output_dimension dimensions of output vector
 * @return int error code : returns 0 if executed successfully
 */
int CreateUnscentedKalmanFilter(Kalman_t* ukf, Matrix_t* Q, Matrix_t* R, GeneralStateModelFunction_t f_state, GeneralStateModelFunction_t f_output, unsigned int state_dimension, unsigned int input_dimension, unsigned int output_dimension);

// Calculate a step of Unscented Kalman Filter
/**
 * @brief 
 * 
 * @param ukf  pointer to Kalman object initialized as a Unscented Kalman Filter
 * @param filter_input pointer to Matrix object that reprsents the filter input (real observed value)
 * @param system_input pointer to Matrix object that reprsents the system input (real manipulated value)
 * @param t float value to represents time varying system
 * @return int error code : returns 0 if executed successfully
 */
int UKF_Operate(Kalman_t* ukf, Matrix_t* filter_input, Matrix_t* system_input, float t);

//Free heap memory held by an object of Kalman Filter
/**
 * @brief 
 * 
 * @param xkf  pointer to Kalman object
 */
void FreeKalman(Kalman_t* xkf);

#endif // !_KALMAN_H_
