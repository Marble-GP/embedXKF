#include "kalman.h"

// Calculate Jacobian with x numerically
int _EKF_CalculateJacobian_X(GeneralStateModelFunction_t f, Matrix_t* A, const Matrix_t* x, const Matrix_t* u, float t)
{
	if (!f || !A->data || A->cols != x->rows)
	{
		return -1;
	}

	Matrix_t eps_vect, input_buf, Ax1, Ax2;

	int error_code = 0;
	float temp_val;
	unsigned int i, j;
	unsigned long long bin_order_ref, bin_order;

	#ifdef USE_STATIC_ARRAY
	float temp_array1[LINALG_VECTOR_BUFFER_SIZE], temp_array2[LINALG_VECTOR_BUFFER_SIZE], temp_array3[LINALG_VECTOR_BUFFER_SIZE], temp_array4[LINALG_VECTOR_BUFFER_SIZE];

	if (A->rows <= LINALG_VECTOR_BUFFER_SIZE && x->rows <= LINALG_VECTOR_BUFFER_SIZE) // when the size is less than LINALG_VECTOR_BUFFER_SIZE
	{
		eps_vect.data = temp_array1;
		input_buf.data = temp_array2;
		Ax1.data = temp_array3;
		Ax2.data = temp_array4;
	}
	else
	{
		eps_vect.data = NULL;
		input_buf.data = NULL;
		Ax1.data = NULL;
		Ax2.data = NULL;
	}

	#else

	eps_vect.data = (float*)STACK_ALLOC(x->rows * sizeof(float));
	input_buf.data = (float*)STACK_ALLOC(x->rows * sizeof(float));
	Ax1.data = (float*)STACK_ALLOC(A->rows * sizeof(float));
	Ax2.data = (float*)STACK_ALLOC(A->rows * sizeof(float));

	#endif // USE_STATIC_ARRAY

	if (!eps_vect.data || !Ax1.data || !Ax2.data)
	{
		return -1;
	}

	eps_vect.rows = input_buf.rows = Ax1.rows = Ax2.rows = x->rows;
	eps_vect.cols = input_buf.cols = Ax1.cols = Ax2.cols = 1;


	for (i = 0; i < A->cols; ++i)
	{
		AssignMatrixWithValue(&eps_vect, 0.0f);

		Mat_at(eps_vect, i, 0) = 1.0f;
		temp_val = Mat_at(*x, i, 0) == 0.0f ? 1.0f : fabsf(Mat_at(*x, i, 0));

		if (temp_val > 1.0f)
		{
			while (temp_val >= 1.0f)
			{
				bin_order_ref = temp_val < (float)((ULLONG_MAX >> 1) + 1ULL) ? (unsigned long long)(temp_val) : ((ULLONG_MAX >> 1ULL) + 1ULL);
				bin_order = 2;

				while (bin_order_ref > 2)
				{
					bin_order_ref >>= 1;
					bin_order <<= 1;
				}

				temp_val /= (float)bin_order;
				Mat_at(eps_vect, i, 0) *= (float)bin_order;
			}
		}
		else
		{
			while (temp_val < 0.5f)
			{

				bin_order_ref = 1.0f/temp_val < (float)((ULLONG_MAX >> 1) + 1ULL) ? (unsigned long long)(1.0f / temp_val) : ((ULLONG_MAX >> 1) + 1ULL);
				bin_order = 2;

				while (bin_order_ref > 2)
				{
					bin_order_ref >>= 1;
					bin_order <<= 1;
				}

				temp_val *= (float)bin_order;
				Mat_at(eps_vect, i, 0) /= (float)bin_order;
			}
		}

		Mat_at(eps_vect, i, 0) *= EKF_JACOBIAN_EPS;

		error_code += AddMatrices(x, &eps_vect, &input_buf);
		f(&Ax1, &input_buf, u, t);
		error_code += SubtractMatrices(x, &eps_vect, &input_buf);
		f(&Ax2, &input_buf, u, t);

		for (j = 0; j < A->rows; ++j)
		{
			Mat_at(*A, j, i) = (Mat_at(Ax1, j, 0) - Mat_at(Ax2, j, 0)) / (2.0f * Mat_at(eps_vect, i, 0));
		}
	}

	return error_code;
}


// Calculate Jacobian with u numerically
int _EKF_CalculateJacobian_U(GeneralStateModelFunction_t f, Matrix_t* B, const Matrix_t* x, const Matrix_t* u, float t)
{
	if (!f || !B->data || B->cols != u->rows )
	{
		return -1;
	}

	Matrix_t eps_vect, input_buf, Bu1, Bu2;

	int error_code = 0;
	float temp_val;
	unsigned int i, j;
	unsigned long long bin_order_ref, bin_order;

#ifdef USE_STATIC_ARRAY
	float temp_array1[LINALG_VECTOR_BUFFER_SIZE], temp_array2[LINALG_VECTOR_BUFFER_SIZE], temp_array3[LINALG_VECTOR_BUFFER_SIZE], temp_array4[LINALG_VECTOR_BUFFER_SIZE];

	if (B->rows <= LINALG_VECTOR_BUFFER_SIZE && u->rows <= LINALG_VECTOR_BUFFER_SIZE) // when the size is less than LINALG_VECTOR_BUFFER_SIZE
	{
		eps_vect.data = temp_array1;
		input_buf.data = temp_array2;
		Bu1.data = temp_array3;
		Bu2.data = temp_array4;
	}
	else
	{
		eps_vect.data = NULL;
		input_buf.data = NULL;
		Bu1.data = NULL;
		Bu2.data = NULL;
	}

#else

	eps_vect.data = (float*)STACK_ALLOC(u->rows * sizeof(float));
	input_buf.data = (float*)STACK_ALLOC(u->rows * sizeof(float));
	Bu1.data = (float*)STACK_ALLOC(B->rows * sizeof(float));
	Bu2.data = (float*)STACK_ALLOC(B->rows * sizeof(float));

#endif // USE_STATIC_ARRAY

	if (!eps_vect.data || !Bu1.data || !Bu2.data)
	{
		return -1;
	}

	eps_vect.rows = input_buf.rows = Bu1.rows = Bu2.rows = u->rows;
	eps_vect.cols = input_buf.cols = Bu1.cols = Bu2.cols = 1;


	for (i = 0; i < B->cols; ++i)
	{
		AssignMatrixWithValue(&eps_vect, 0.0f);

		Mat_at(eps_vect, i, 0) = 1.0f;
		temp_val = Mat_at(*u, i, 0) == 0.0f ? 1.0f : fabsf(Mat_at(*u, i, 0));

		if (temp_val > 1.0f)
		{
			while (temp_val >= 1.0f)
			{
				bin_order_ref = temp_val < (float)((ULLONG_MAX >> 1) + 1ULL) ? (unsigned long long)(temp_val) : ((ULLONG_MAX >> 1) + 1ULL);
				bin_order = 2;

				while (bin_order_ref > 2)
				{
					bin_order_ref >>= 1;
					bin_order <<= 1;
				}

				temp_val /= (float)bin_order;
				Mat_at(eps_vect, i, 0) *= (float)bin_order;
			}
		}
		else
		{
			while (temp_val < 0.5f)
			{
				bin_order_ref = 1.0f / temp_val < (float)((ULLONG_MAX >> 1) + 1ULL) ? (unsigned long long)(1.0f / temp_val) : ((ULLONG_MAX >> 1) + 1ULL);
				bin_order = 2;

				while (bin_order_ref > 2)
				{
					bin_order_ref >>= 1;
					bin_order <<= 1;
				}

				temp_val *= (float)bin_order;
				Mat_at(eps_vect, i, 0) /= (float)bin_order;
			}
		}

		Mat_at(eps_vect, i, 0) *= EKF_JACOBIAN_EPS;

		error_code += AddMatrices(u, &eps_vect, &input_buf);
		f(&Bu1, x, &input_buf, t);
		error_code += SubtractMatrices(u, &eps_vect, &input_buf);
		f(&Bu2, x, &input_buf, t);
		for (j = 0; j < B->rows; ++j)
		{
			Mat_at(*B, j, i) = (Mat_at(Bu1, j, 0) - Mat_at(Bu2, j, 0)) / (2.0f * Mat_at(eps_vect, i, 0));
		}
	}

	return error_code;
}

// Calculate sqrt of matrix with diagnosing
int _UKF_AbusoluteMatrixSQRT(const Matrix_t* source, Matrix_t* destination)
{
	if (source->rows != source->cols || destination->rows != destination->cols || source->rows != destination->rows)
	{
		// Error: not a symmetric matrix
		return -1;
	}

	Matrix_t P;
	unsigned int i;
	int error_code = 0;

#ifdef USE_STATIC_ARRAY

	float P_array[LINALG_MATRIX_BUFFER_SIZE];
	P.data = LINALG_MATRIX_BUFFER_SIZE >= source->rows * source->cols ? P_array : NULL;
	
# else

	P.data = STACK_ALLOC(source->rows * source->cols * sizeof(float));

#endif

	if (!P.data)
	{
		return -1;
	}

	P.rows = P.cols = source->rows;

	error_code += DiagonalizeMatrix(source, destination, &P);

	for (i = 0; i < source->rows; ++i)
	{
		Mat_at(*destination, i, i) = sqrtf(fabsf(Mat_at(*destination, i, i)));
	}


	error_code += MultiplyMatrices(&P, destination, destination);
	error_code += TransposeMatrix(&P, &P);
	error_code += MultiplyMatrices(destination, &P, destination);


	return error_code;
}


//Create a Linear Kalman Filter
int CreateKalmanFilter(Kalman_t* kf, StateSpace_t* ss, Matrix_t* Q, Matrix_t* R)
{
	if (!kf || !ss || !Q || !R || !(Q->rows == Q->cols && Q->rows == ss->u.rows) || !(R->rows == R->cols && (R->rows == ss->y.rows || R->rows == ss->x.rows)) )
	{
		return -1;
	}
	
	// set NULL when this not used as an EKF
	kf->state_model_func = NULL;
	kf->state_model_jacobian_x = NULL;
	kf->state_model_jacobian_u = NULL;
	kf->output_model_func = NULL;
	kf->output_model_jacobian_x = NULL;
	kf->output_model_jacobian_u = NULL;

	// set NULL  this not used as an UKF
	kf->SP_state = NULL;
	kf->SP_output = NULL;

	//set model parameters
	kf->linsys = ss;
	kf->Q = *Q;
	kf->R = *R;

	//create gain matrix
	kf->G = CreateMatrixWithInitValue(ss->A.rows, ss->A.cols, 0.0f);

	//crate filtered state vector
	kf->x_filter = CopyMatrix(&ss->x);

	//create noise matrices
	kf->P_adv = CreateMatrixWithInitValue(ss->A.rows, ss->A.cols, 0.0f);
	kf->P = CreateIdentityMatrix(ss->A.rows);

	// if it is allocated successfully
	if (kf->P.data)
	{
		MultiplyMatrixByScalar(&kf->P, P0_INIT_VAL, &kf->P);
	}
	
	if (kf->G.data && kf->P.data && kf->P_adv.data)
	{
		return 0;
	}
	else
	{
		FreeKalman(kf);
		return -1;
	}
}


// Calculate Kalman Filter Step
int KF_Operate(Kalman_t* kf, Matrix_t* filter_input, Matrix_t* system_input)
{
	Matrix_t temp_nn, temp_rn, temp_pp, e_vect;
	
	int error_code = 0;

	#ifdef USE_STATIC_ARRAY 

	float temp_nn_array[LINALG_MATRIX_BUFFER_SIZE], temp_rn_array[LINALG_MATRIX_BUFFER_SIZE], temp_pp_array[LINALG_MATRIX_BUFFER_SIZE], e_vect_array[LINALG_VECTOR_BUFFER_SIZE];

	temp_nn.data = temp_nn_array;
	temp_rn.data = temp_rn_array;
	temp_pp.data = temp_pp_array;
	e_vect.data = e_vect_array;

	#else

	unsigned int P_size = kf->P.rows * kf->P.cols, R_size = kf->R.rows * kf->R.cols, rn_size = kf->P.rows * kf->Q.cols, e_size = kf->linsys->y.rows * kf->linsys->y.cols;

	temp_nn.data = (float*)STACK_ALLOC(P_size * sizeof(float));
	temp_pp.data = (float*)STACK_ALLOC(R_size * sizeof(float));
	temp_rn.data = (float*)STACK_ALLOC(rn_size * sizeof(float));
	e_vect.data = (float*)STACK_ALLOC(e_size * sizeof(float));


	#endif //!USE_STATIC_ARRAY

	temp_nn.rows = temp_nn.cols = temp_rn.cols = kf->P.rows;
	temp_rn.rows = kf->Q.rows;
	temp_pp.rows = temp_pp.cols = kf->R.rows;
	e_vect.rows = kf->linsys->C.rows;
	e_vect.cols = 1;

	// allocation check
	if (!temp_nn.data || (&kf->linsys->B && !temp_rn.data) || (kf->linsys->C.data && (!temp_pp.data || !e_vect.data)))
	{
		return -1;
	}

	// ***** Prediction step *****

	AssignMatrix(&kf->linsys->x, kf->x_filter.data, kf->linsys->x.rows);
	error_code += SSmodelStep(kf->linsys, system_input);

	error_code += TransposeMatrix(&kf->linsys->A, &temp_nn);
	error_code += MultiplyMatrices(&kf->P, &temp_nn, &temp_nn);// temp_nn = P * At
	error_code += MultiplyMatrices(&kf->linsys->A, &temp_nn, &kf->P_adv);//P_adv = A * (P * At)

	if (kf->linsys->B.data) // consider the effect of inputs on the system (B matrix exists)
	{

		error_code += TransposeMatrix(&kf->linsys->B, &temp_rn);
		error_code += MultiplyMatrices(&kf->Q, &temp_rn, &temp_rn);// temp_rn = Q*Bt
		error_code += MultiplyMatrices(&kf->linsys->B, &temp_rn, &temp_nn);//temp_nn = B*(Q*Bt)

		error_code += AddMatrices(&kf->P_adv, &temp_nn, &kf->P_adv); // complete calculating P_adv(n)
	}

	// ***** Filtering step *****

	if (kf->linsys->C.data) // consider the output system structure (C matrix exists)
	{
		error_code += TransposeMatrix(&kf->linsys->C, &kf->G); //use G as a n�~p  buffer
		error_code += MultiplyMatrices(&kf->P_adv, &kf->G, &kf->G); // G = P_adv * Ct
		error_code += MultiplyMatrices(&kf->linsys->C, &kf->G, &temp_pp);// temp_pp = C * (P_adv * Ct)

		error_code += AddMatrices(&kf->R, &temp_pp, &temp_pp);// temp_pp = C * (P_adv * Ct) + R

		error_code += InvertMatrix(&temp_pp, &temp_pp); // temp_pp = (C * (P_adv * Ct) + R)^-1  (this is like dens)

		error_code += MultiplyMatrices(&kf->G, &temp_pp, &kf->G); // G = (P_adv * Ct) * (C * (P_adv * Ct) + R)^-1

		error_code += SubtractMatrices(filter_input, &kf->linsys->y, &e_vect); // e_vect = y - y_adv
		error_code += MultiplyMatrices(&kf->G, &e_vect, &kf->x_filter); // x_filter = G * (y - y_adv)
		error_code += AddMatrices(&kf->linsys->x, &kf->x_filter, &kf->x_filter); // x_filter = x_adv + G * (y - y_adv)


		error_code += MultiplyMatrices(&kf->G, &kf->linsys->C, &kf->P); //P = G * C, use P as a n�~n  buffer

		AssignIdentityMatrix(&temp_nn);//temp_nn = I

		error_code += SubtractMatrices(&temp_nn, &kf->P, &temp_nn); //temp_nn = I - G * C
		error_code += MultiplyMatrices(&temp_nn, &kf->P_adv, &kf->P); // P = (I - G * C) * P_adv

	}
	else // if C.p_mat is NULL, it is considered y == x (i.e. C == I).
	{
		error_code += AddMatrices(&kf->R, &kf->P_adv, &temp_nn); // temp_nn = C * (P_adv * Ct) + R = P_adv + R
		error_code += InvertMatrix(&temp_nn, &temp_nn);// temp_nn = (P_adv + R)^-1
		error_code += MultiplyMatrices(&kf->P_adv, &temp_nn, &kf->G); // G = P_adv * (P_adv + R)^-1


		error_code += SubtractMatrices(filter_input, &kf->linsys->x, &kf->x_filter); // x_filter = y - y_adv, use x_filter as a buffer
		error_code += MultiplyMatrices(&kf->G, &kf->x_filter, &kf->x_filter); // x_filter = G * (y - y_adv)
		error_code += AddMatrices(&kf->linsys->x, &kf->x_filter, &kf->x_filter); // x_filter = x_adv + G * (y - y_adv)

		AssignIdentityMatrix(&temp_nn);//temp_nn = I

		error_code += SubtractMatrices(&temp_nn, &kf->G, &temp_nn); //temp_nn = I - G * C = I - G
		error_code += MultiplyMatrices(&temp_nn, &kf->P_adv, &kf->P);// P = (I - G) * P_adv
	}


	return error_code;
}

// Create a Extended Kalman Filter
int CreateExtendedKalmanFilter(Kalman_t* ekf, Matrix_t* Q, Matrix_t* R, GeneralStateModelFunction_t f_state, GeneralStateModelFunction_t f_output, GeneralStateModelFunction_t f_state_jacobian_x, GeneralStateModelFunction_t f_state_jacobian_u, GeneralStateModelFunction_t f_output_jacobian_x, GeneralStateModelFunction_t f_output_jacobian_u, unsigned int state_dimension, unsigned int input_dimension, unsigned int output_dimension)
{
	int error_code = 0;

	if (!ekf || !Q || !R || !f_state || !f_output || !(Q->rows == Q->cols && Q->rows == input_dimension) || !(R->rows == R->cols && R->rows == output_dimension))
	{
		return -1;
	}

	// set NULL  this not used as an UKF
	ekf->SP_state = NULL;
	ekf->SP_output = NULL;

	ekf->Q = *Q;
	ekf->R = *R;
	ekf->state_model_func = f_state;
	ekf->state_model_jacobian_x = f_state_jacobian_x;
	ekf->state_model_jacobian_u = f_state_jacobian_u;
	ekf->output_model_func = f_output;
	ekf->output_model_jacobian_x = f_output_jacobian_x;
	ekf->output_model_jacobian_u = f_output_jacobian_u;

	//create gain matrix
	ekf->G = CreateMatrixWithInitValue(state_dimension, output_dimension, 0.0f);

	//crate filtered state vector
	ekf->x_filter = CreateMatrixWithInitValue(state_dimension, 1, 0.0f);

	//create noise matrices
	ekf->P_adv = CreateMatrixWithInitValue(state_dimension, state_dimension, 0.0f);
	ekf->P = CreateIdentityMatrix(state_dimension);



	// if it is allocated successfully
	if (ekf->P.data)
	{
		error_code += MultiplyMatrixByScalar(&ekf->P, P0_INIT_VAL, &ekf->P);
	}

	//create linsys from heap memory
	ekf->linsys = malloc(sizeof(StateSpace_t));
	if(ekf->linsys)// create linsys automatically when this not used as a LKF
	{
		ekf->linsys->A = CreateMatrix(state_dimension, state_dimension);
		ekf->linsys->B = CreateMatrix(state_dimension, input_dimension);
		ekf->linsys->C = CreateMatrix(output_dimension, state_dimension);
		ekf->linsys->D = CreateMatrix(output_dimension, input_dimension);
		ekf->linsys->x = CreateMatrix(state_dimension, 1);
		ekf->linsys->u = CreateMatrix(input_dimension, 1);
		ekf->linsys->y = CreateMatrix(output_dimension, 1);


		if (f_state_jacobian_x) // if it is not NULL
		{
			f_state_jacobian_x(&ekf->linsys->A, &ekf->x_filter, &ekf->linsys->u, 0.0f);
		}
		else
		{
			error_code += _EKF_CalculateJacobian_X(f_state, &ekf->linsys->A, &ekf->linsys->x, &ekf->linsys->u, 0.0f);
		}

		if (f_state_jacobian_u)// if it is not NULL
		{
			f_state_jacobian_u(&ekf->linsys->B, &ekf->x_filter, &ekf->linsys->u, 0.0f);
		}
		else
		{
			error_code += _EKF_CalculateJacobian_U(f_state, &ekf->linsys->B, &ekf->linsys->x, &ekf->linsys->u, 0.0f);
		}

		if (f_output) // if it has output system
		{
			if (f_output_jacobian_x)// if it is not NULL
			{
				f_output_jacobian_x(&ekf->linsys->C, &ekf->x_filter, &ekf->linsys->u, 0.0f);

			}
			else
			{
				error_code += _EKF_CalculateJacobian_X(f_output, &ekf->linsys->C, &ekf->linsys->x, &ekf->linsys->u, 0.0f);
			}

			if (f_output_jacobian_u)// if it is not NULL
			{
				f_output_jacobian_u(&ekf->linsys->D, &ekf->x_filter, &ekf->linsys->u, 0.0f);

			}
			else
			{
				error_code += _EKF_CalculateJacobian_U(f_output, &ekf->linsys->D, &ekf->linsys->x, &ekf->linsys->u, 0.0f);
			}
		}
		

	}
	else
	{
		error_code = -1;
	}

	if (!error_code && ekf->x_filter.data && ekf->G.data && ekf->P.data && ekf->P_adv.data)
	{
		return 0;
	}
	else
	{

		FreeKalman(ekf);

		return -1;
	}
}


// Calculate a step of Extended Kalman Filter
int EKF_Operate(Kalman_t* ekf, Matrix_t* filter_input, Matrix_t* system_input, float t)
{
	Matrix_t temp_nn, temp_rn, temp_pp, e_vect;
	
	int error_code = 0;

#ifdef USE_STATIC_ARRAY 

	float temp_nn_array[LINALG_MATRIX_BUFFER_SIZE], temp_rn_array[LINALG_MATRIX_BUFFER_SIZE], temp_pp_array[LINALG_MATRIX_BUFFER_SIZE], e_vect_array[LINALG_VECTOR_BUFFER_SIZE];
	temp_nn.data = temp_nn_array;
	temp_rn.data = temp_rn_array;
	temp_pp.data = temp_pp_array;
	e_vect.data = e_vect_array;

#else
	const unsigned int P_size = ekf->P.rows * ekf->P.cols, R_size = ekf->R.rows * ekf->R.cols, rn_size = ekf->P.rows * ekf->Q.cols, e_size = ekf->linsys->y.rows * ekf->linsys->y.cols;

	temp_nn.data = (float*)STACK_ALLOC(P_size * sizeof(float));
	temp_pp.data = (float*)STACK_ALLOC(R_size * sizeof(float));
	temp_rn.data = (float*)STACK_ALLOC(rn_size * sizeof(float));
	e_vect.data = (float*)STACK_ALLOC(e_size * sizeof(float));


#endif //!USE_STATIC_ARRAY


	temp_nn.rows = temp_nn.cols = temp_rn.cols = ekf->P.rows;
	temp_rn.rows = ekf->Q.rows;
	temp_pp.rows = temp_pp.cols = ekf->R.rows;
	e_vect.rows = ekf->linsys->C.rows;
	e_vect.cols = 1;

	// allocation check
	if (!temp_nn.data || (&ekf->linsys->B && !temp_rn.data) || (ekf->linsys->C.data && (!temp_pp.data || !e_vect.data)))
	{
		return -1;
	}

	// ***** Prediction step *****
	ekf->state_model_func(&ekf->linsys->x, &ekf->x_filter, system_input, t);
	ekf->output_model_func(&ekf->linsys->y, &ekf->x_filter, system_input, t);

	if (ekf->state_model_jacobian_x)
	{
		ekf->state_model_jacobian_x(&ekf->linsys->A, &ekf->x_filter, system_input, t);
	}
	else
	{
		error_code += _EKF_CalculateJacobian_X(ekf->state_model_func, &ekf->linsys->A, &ekf->x_filter, system_input, t);
	}

	if (ekf->state_model_jacobian_u)
	{
		ekf->state_model_jacobian_u(&ekf->linsys->B, &ekf->x_filter, system_input, t);
	}
	else
	{
		error_code += _EKF_CalculateJacobian_U(ekf->state_model_func, &ekf->linsys->B, &ekf->x_filter, system_input, t);
	}

	if (ekf->output_model_func)
	{
		if (ekf->output_model_jacobian_x)
		{
			ekf->output_model_jacobian_x(&ekf->linsys->C, &ekf->x_filter, system_input, t);
		}
		else
		{
			error_code += _EKF_CalculateJacobian_X(ekf->output_model_func, &ekf->linsys->C, &ekf->x_filter, system_input, t);
		}

		if (ekf->output_model_jacobian_u)
		{
			ekf->output_model_jacobian_u(&ekf->linsys->D, &ekf->x_filter, system_input, t);
		}
		else
		{
			error_code += _EKF_CalculateJacobian_U(ekf->output_model_func, &ekf->linsys->D, &ekf->x_filter, system_input, t);
		}
	}
	

	error_code += TransposeMatrix(&ekf->linsys->A, &temp_nn);
	error_code += MultiplyMatrices(&ekf->P, &temp_nn, &temp_nn);// temp_nn = P * At
	error_code += MultiplyMatrices(&ekf->linsys->A, &temp_nn, &ekf->P_adv);//P_adv = A * (P * At)


	error_code += TransposeMatrix(&ekf->linsys->B, &temp_rn);
	error_code += MultiplyMatrices(&ekf->Q, &temp_rn, &temp_rn);// temp_rn = Q*Bt
	error_code += MultiplyMatrices(&ekf->linsys->B, &temp_rn, &temp_nn);//temp_nn = B*(Q*Bt)

	error_code += AddMatrices(&ekf->P_adv, &temp_nn, &ekf->P_adv); // complete calculating P_adv(n)

	// ***** Filtering step *****
	if (ekf->output_model_func)
	{
		error_code += TransposeMatrix(&ekf->linsys->C, &ekf->G); //use G as a n�~p  buffer
		error_code += MultiplyMatrices(&ekf->P_adv, &ekf->G, &ekf->G); // G = P_adv * Ct
		error_code += MultiplyMatrices(&ekf->linsys->C, &ekf->G, &temp_pp);// temp_pp = C * (P_adv * Ct)

		error_code += AddMatrices(&ekf->R, &temp_pp, &temp_pp);// temp_pp = C * (P_adv * Ct) + R

		error_code += InvertMatrix(&temp_pp, &temp_pp); // temp_pp = (C * (P_adv * Ct) + R)^-1 

		error_code += MultiplyMatrices(&ekf->G, &temp_pp, &ekf->G); // G = (P_adv * Ct) * (C * (P_adv * Ct) + R)^-1

		error_code += SubtractMatrices(filter_input, &ekf->linsys->y, &e_vect); // e_vect = y - y_adv
		error_code += MultiplyMatrices(&ekf->G, &e_vect, &ekf->x_filter); // x_filter = G * (y - y_adv)
		error_code += AddMatrices(&ekf->linsys->x, &ekf->x_filter, &ekf->x_filter); // x_filter = x_adv + G * (y - y_adv)


		error_code += MultiplyMatrices(&ekf->G, &ekf->linsys->C, &ekf->P); //P = G * C, use P as a n�~n  buffer

		AssignIdentityMatrix(&temp_nn);//temp_nn = I

		error_code += SubtractMatrices(&temp_nn, &ekf->P, &temp_nn); //temp_nn = I - G * C
		error_code += MultiplyMatrices(&temp_nn, &ekf->P_adv, &ekf->P); // P = (I - G * C) * P_adv
	}
	else
	{
		error_code += AddMatrices(&ekf->R, &ekf->P_adv, &temp_pp); // temp_pp = R + P_adv
		error_code += InvertMatrix(&temp_pp, &temp_pp); // temp_pp = (P_adv + R)^-1
		error_code += MultiplyMatrices(&ekf->P_adv, &temp_nn, &ekf->G); // G = P_adv * (P_adv+ R)^-1


		error_code += SubtractMatrices(filter_input, &ekf->linsys->y, &e_vect); // e_vect = y - y_adv
		error_code += MultiplyMatrices(&ekf->G, &e_vect, &ekf->x_filter); // x_filter = G * (y - y_adv)
		error_code += AddMatrices(&ekf->linsys->x, &ekf->x_filter, &ekf->x_filter); // x_filter = x_adv + G * (y - y_adv)


		AssignIdentityMatrix(&temp_nn);//temp_nn = I

		error_code += SubtractMatrices(&temp_nn, &ekf->P, &temp_nn); //temp_nn = I - G
		error_code += MultiplyMatrices(&temp_nn, &ekf->P_adv, &ekf->P); // P = (I - G) * P_adv
	}

	



	return error_code;

}

// Create an object of Unscented Kalman Filter
int CreateUnscentedKalmanFilter(Kalman_t* ukf, Matrix_t* BQBt, Matrix_t* R, GeneralStateModelFunction_t f_state, GeneralStateModelFunction_t f_output, unsigned int state_dimension, unsigned int input_dimension, unsigned int output_dimension)
{
	int error_code = 0;
	unsigned int i;
	const float SP_init_val = sqrtf(state_dimension*P0_INIT_VAL);

	// create ss automatically
	ukf->linsys = malloc(sizeof(StateSpace_t));
	if (ukf->linsys)
	{
		ukf->linsys->A = CreateMatrix(0, 0);
		ukf->linsys->B = CreateMatrix(0, 0);
		ukf->linsys->C = CreateMatrix(0, 0);
		ukf->linsys->D = CreateMatrix(0, 0);
		ukf->linsys->x = CreateMatrix(state_dimension, 1); // it is used as a buffer 
		ukf->linsys->u = CreateMatrix(input_dimension, 1); // it is used as a buffer
		ukf->linsys->y = CreateMatrix(output_dimension, 1); // it is used as a buffer
	}

	// set NULL if this not used as an KF or EKF
	ukf->state_model_jacobian_x = NULL;
	ukf->state_model_jacobian_u = NULL;
	ukf->output_model_jacobian_x = NULL;
	ukf->output_model_jacobian_u = NULL;


	ukf->Q = *BQBt;
	ukf->R = *R;
	ukf->state_model_func = f_state;
	ukf->output_model_func = f_output;
	ukf->k = 0.0f;


	//create gain matrix
	ukf->G = CreateMatrixWithInitValue(state_dimension, output_dimension, 0.0f);

	//crate filtered state vector
	ukf->x_filter = CreateMatrixWithInitValue(state_dimension, 1, 0.0f);

	//create noise matrices
	ukf->P_adv = CreateMatrixWithInitValue(state_dimension, state_dimension, 0.0f);
	ukf->P = CreateIdentityMatrix(state_dimension);

	// if it is allocated successfully
	if (ukf->P.data)
	{
		MultiplyMatrixByScalar(&ukf->P, P0_INIT_VAL, &ukf->P);
	}

	// Sigma points array
	ukf->SP_state = (Matrix_t*)malloc((2 * state_dimension + 1) * sizeof(Matrix_t));
	ukf->SP_output = (Matrix_t*)malloc((2 * state_dimension + 1) * sizeof(Matrix_t));


	// if it is allocated successfully
	if (ukf->SP_state && ukf->SP_output)
	{
		ukf->SP_state[0] = CreateMatrixWithInitValue(state_dimension, 1, 0.0f);
		ukf->SP_output[0] = CreateMatrixWithInitValue(output_dimension, 1, 0.0f);

		for (i = 1; i <= state_dimension; ++i)
		{
			ukf->SP_state[i] = CreateMatrixWithInitValue(state_dimension, 1, 0.0f);
			ukf->SP_output[i] = CreateMatrixWithInitValue(output_dimension, 1, 0.0f);
			ukf->SP_state[state_dimension + i] = CreateMatrixWithInitValue(state_dimension, 1, 0.0f);
			ukf->SP_output[state_dimension + i] = CreateMatrixWithInitValue(output_dimension, 1, 0.0f);

			error_code += (ukf->SP_state + i)->data ? 0 : -1;
			error_code += (ukf->SP_output + i)->data ? 0 : -1;
			error_code += (ukf->SP_state + state_dimension + i)->data ? 0 : -1;
			error_code += (ukf->SP_output + state_dimension + i)->data ? 0 : -1;

			if ((ukf->SP_state + i)->data && (ukf->SP_state + state_dimension + i)->data)
			{
				Mat_at(ukf->SP_state[i], i - 1, 0) = SP_init_val;
				Mat_at(ukf->SP_state[state_dimension + i], i - 1, 0) = -SP_init_val;
			}
		}
	}
	else
	{
		free(ukf->SP_state);
		free(ukf->SP_output);
		error_code = -1;
	}



	if (ukf->G.data && ukf->P.data && ukf->P_adv.data && !error_code)
	{
		return 0;
	}
	else
	{

		FreeKalman(ukf);
		
		return -1;
	}


}


// Calculate a step of Unscented Kalman Filter
int UKF_Operate(Kalman_t* ukf, Matrix_t* filter_input, Matrix_t* system_input, float t)
{
	const unsigned int state_dims = ukf->SP_state->rows, output_dims = ukf->SP_output->rows;
	const float weight_0 = ukf->k / ((float)state_dims + ukf->k), weight_i = 1.0f / (2.0f * ((float)state_dims + ukf->k)), k_sp = sqrtf((float)state_dims + ukf->k);
	Matrix_t temp_nn, temp_pn, temp_np1, temp_np2, temp_pp1, temp_pp2, temp_y_vect, temp_y_vect_T, temp_x_vect, temp_x_vect_T;
	int error_code = 0;
	unsigned int i;

#ifdef USE_STATIC_ARRAY 

	float temp_nn_array[LINALG_MATRIX_BUFFER_SIZE], temp_pn_array[LINALG_MATRIX_BUFFER_SIZE], temp_np1_array[LINALG_MATRIX_BUFFER_SIZE], temp_np2_array[LINALG_MATRIX_BUFFER_SIZE], temp_pp1_array[LINALG_MATRIX_BUFFER_SIZE], temp_pp2_array[LINALG_MATRIX_BUFFER_SIZE], temp_x_vect_array[LINALG_VECTOR_BUFFER_SIZE], temp_y_vect_array[LINALG_VECTOR_BUFFER_SIZE];

	temp_nn.data = temp_nn_array;
	temp_pn.data = temp_pn_array;
	temp_np1.data = temp_np1_array;
	temp_np2.data = temp_np2_array;
	temp_pp1.data = temp_pp1_array;
	temp_pp2.data = temp_pp2_array;
	temp_x_vect.data = temp_x_vect_T.data = temp_x_vect_array;
	temp_y_vect.data = temp_y_vect_T.data = temp_y_vect_array;
	

#else
	const unsigned int P_size = ukf->P.rows * ukf->P.cols, R_size = ukf->R.rows * ukf->R.cols, rn_size = ukf->P.rows * ukf->Q.cols;

	temp_nn.data = (float*)STACK_ALLOC(P_size * sizeof(float));
	temp_pp1.data = (float*)STACK_ALLOC(R_size * sizeof(float));
	temp_pp2.data = (float*)STACK_ALLOC(R_size * sizeof(float));
	temp_pn.data = (float*)STACK_ALLOC(rn_size * sizeof(float));
	temp_np1.data = (float*)STACK_ALLOC(rn_size * sizeof(float));
	temp_np2.data = (float*)STACK_ALLOC(rn_size * sizeof(float));
	temp_y_vect.data = temp_y_vect_T.data = (float*)STACK_ALLOC(output_dims * sizeof(float));
	temp_x_vect.data = temp_x_vect_T.data = (float*)STACK_ALLOC(state_dims * sizeof(float));


#endif //!USE_STATIC_ARRAY
	


	temp_x_vect.rows = temp_x_vect_T.cols = temp_nn.rows = temp_nn.cols = temp_pn.cols = temp_np1.rows = temp_np2.rows = state_dims;
	temp_y_vect.rows = temp_y_vect_T.cols = temp_pn.rows = temp_np1.cols = temp_np2.cols = output_dims;
	temp_pp2.rows = temp_pp2.cols = temp_pp1.rows = temp_pp1.cols = ukf->R.rows;

	temp_x_vect_T.rows = temp_x_vect.cols = temp_y_vect_T.rows = temp_y_vect.cols = 1;


	// allocation check
	if (!temp_nn.data || !temp_pn.data || !temp_np1.data || !temp_np2.data || !temp_pp1.data || !temp_pp2.data || !temp_x_vect.data || !temp_y_vect.data)
	{
		return -1;
	}

	// Calculate prediction sigma points with params that calculated in previous step. 
	ukf->state_model_func(&temp_x_vect, &(ukf->SP_state[0]), system_input, t);
	error_code += MultiplyMatrixByScalar(&temp_x_vect, weight_0, &ukf->linsys->x); //linsys.x is used as a buffer for x_prediction

	for (i = 1; i <= state_dims; ++i)
	{

		ukf->state_model_func(&temp_x_vect, &(ukf->SP_state[i]), system_input, t);
		error_code += AppendMatrixWithWeight(&ukf->linsys->x, &temp_x_vect, weight_i); // x_prediction += weight * F(SP_state[i])

		AssignMatrix(&(ukf->SP_state[i]), temp_x_vect.data, state_dims); // update SP_state[i]

		ukf->state_model_func(&temp_x_vect, &(ukf->SP_state[state_dims + i]), system_input, t);
		error_code += AppendMatrixWithWeight(&ukf->linsys->x, &temp_x_vect, weight_i); // x_prediction += weight * F(SP_state[n+i])

		AssignMatrix(&(ukf->SP_state[state_dims + i]), temp_x_vect.data, state_dims); // update SP_state[n+i]

	}

	//Calculate P_adv
	AssignMatrixWithValue(&ukf->P_adv, 0.0f); // P_adv = O
	error_code += SubtractMatrices(&(ukf->SP_state[0]), &ukf->linsys->x, &temp_x_vect); // temp_x = (SP_state[0] - x_prediction)
	error_code += MultiplyMatrices(&temp_x_vect, &temp_x_vect_T, &temp_nn); // temp_nn = temp_x * temp_x_T
	error_code += AppendMatrixWithWeight(&ukf->P_adv, &temp_nn, weight_0); // P_adv += w0 * (temp_x * temp_x_T)

	for (i = 1; i <= state_dims; ++i)
	{
		error_code += SubtractMatrices(&(ukf->SP_state[i]), &ukf->linsys->x, &temp_x_vect); // temp_x = (SP_state[i] - x_prediction)
		error_code += MultiplyMatrices(&temp_x_vect, &temp_x_vect_T, &temp_nn); // temp_nn = temp_x * temp_x_T
		error_code += AppendMatrixWithWeight(&ukf->P_adv, &temp_nn, weight_i);  // P_adv += w_i * (temp_x * temp_x_T)



		error_code += SubtractMatrices(&(ukf->SP_state[state_dims + i]), &ukf->linsys->x, &temp_x_vect); // temp_x = (SP_state[n+i] - x_prediction)
		error_code += MultiplyMatrices(&temp_x_vect, &temp_x_vect_T, &temp_nn); // temp_nn = temp_x * temp_x_T
		error_code += AppendMatrixWithWeight(&ukf->P_adv, &temp_nn, weight_i);  // P_adv += w_i * (temp_x * temp_x_T)

	}

	AddMatrices(&ukf->P_adv, &ukf->Q, &ukf->P_adv); // P_adv += BQBt


	// Re-calculate state sigma points and output sigma points with params that calculated in this step. 
	AssignMatrix(&(ukf->SP_state[0]), ukf->linsys->x.data, state_dims);  //SP_state[0] = x_prediction
	ukf->output_model_func(&(ukf->SP_output[0]), &(ukf->SP_state[0]), system_input, t); //calculate SP_output[0]
	//linsys.y is used as a buffer for y_prediction
	error_code += MultiplyMatrixByScalar(ukf->SP_output, weight_0, &ukf->linsys->y); // y_prediction = w0*SP_output[0]

	#ifdef USE_INCORRECT_MATRIX_SQRT_EXTENSION

	error_code += _UKF_AbusoluteMatrixSQRT(&ukf->P_adv, &temp_nn); //temp_nn = SQRT(P_adv)

	#else

	error_code += CholeskyDecomposition(&ukf->P_adv, &temp_nn); //temp_nn = SQRT(P_adv)

	#endif

	error_code += MultiplyMatrixByScalar(&temp_nn, k_sp, &temp_nn); // temp_nn = sqrt(n + k) * SQRT(P_adv)



	for (i = 1; i <= state_dims; ++i)
	{
		error_code += GetRowVector(&temp_nn, &temp_x_vect, i - 1); // temp_x_vect = sqrt(n + k) * SQRT(P_adv)[i]
		error_code += AddMatrices(&ukf->linsys->x, &temp_x_vect, &(ukf->SP_state[i])); // SP_state[i] = x_prediction + sqrt(n + k) * SQRT(P_adv)[i] (i = 1, 2, ... , n)
		error_code += SubtractMatrices(&ukf->linsys->x, &temp_x_vect, &(ukf->SP_state[i + state_dims]));  // SP_state[n + i] = x_prediction - sqrt(n + k) * SQRT(P_adv)[i] (i = 1, 2, ... , n)
		
		//calculate SP_output[i] and y_prediction
		ukf->output_model_func(&(ukf->SP_output[i]), &(ukf->SP_state[i]), system_input, t);
		ukf->output_model_func(&(ukf->SP_output[i + state_dims]), &(ukf->SP_state[i + state_dims]), system_input, t);
		error_code += AppendMatrixWithWeight(&ukf->linsys->y, &(ukf->SP_output[i]), weight_i); //y_prediction += w_i * SP_output[i]
		error_code += AppendMatrixWithWeight(&ukf->linsys->y, &(ukf->SP_output[i + state_dims]), weight_i); //y_prediction += w_i * SP_output[i+n]
	}

	// Calculate Pyy and Pxy
	//Pyy 
	error_code += SubtractMatrices(&(ukf->SP_output[0]), &ukf->linsys->y, &temp_y_vect); //temp_y = SP_output[0] - y_prediction
	error_code += MultiplyMatrices(&temp_y_vect, &temp_y_vect_T, &temp_pp1);
	error_code += MultiplyMatrixByScalar(&temp_pp1, weight_0, &temp_pp1);  // temp_pp1 = w0 * (SP_output[0] - y_prediction) * (SP_output[0] - y_prediction)^T
	//Pxy
	error_code += SubtractMatrices(&(ukf->SP_state[0]), &ukf->linsys->x, &temp_x_vect);
	error_code += MultiplyMatrices(&temp_x_vect, &temp_y_vect_T, &temp_np1);
	error_code += MultiplyMatrixByScalar(&temp_np1, weight_0, &temp_np1); // temp_np1 = w0 * (SP_state[0] - x_prediction) * (SP_output[0] - y_prediction)^T
	
	for (i = 1; i <= state_dims; ++i)
	{
		//Pyy (i element)
		error_code += SubtractMatrices(&(ukf->SP_output[i]), &ukf->linsys->y, &temp_y_vect); // temp_y = w_i * (SP_output[i] - y_prediction)	
		error_code += MultiplyMatrices(&temp_y_vect, &temp_y_vect_T, &temp_pp2); // temp_pp2 = w_i * (SP_output[i] - y_prediction) * (SP_output[i] - y_prediction)^T
		error_code += AppendMatrixWithWeight(&temp_pp1, &temp_pp2, weight_i); // temp_pp1 += w_i * temp_pp2
		//Pxy (i element)
		error_code += SubtractMatrices(&(ukf->SP_state[i]), &ukf->linsys->x, &temp_x_vect); // temp_x = (SP_state[i] - x_prediction)
		error_code += MultiplyMatrices(&temp_x_vect, &temp_y_vect_T, &temp_np2); // temp_np2 = (SP_state[i] - x_prediction) * (SP_output[i] - y_prediction)^T
		error_code += AppendMatrixWithWeight(&temp_np1, &temp_np2, weight_i); // temp_np1 += w_i * (SP_state[i] - x_prediction) * (SP_output[i] - y_prediction)^T


		//Pyy (n+i element)
		error_code += SubtractMatrices(&(ukf->SP_output[state_dims + i]), &ukf->linsys->y, &temp_y_vect); // temp_y = w_i * (SP_output[n+i] - y_prediction)
		error_code += MultiplyMatrices(&temp_y_vect, &temp_y_vect_T, &temp_pp2); // temp_pp2 = w_i * (SP_output[i] - y_prediction) * (SP_output[n+i] - y_prediction)^T
		error_code += AppendMatrixWithWeight(&temp_pp1, &temp_pp2, weight_i); // temp_pp1 += w_i * temp_pp2
		//Pxy (n+i element)
		error_code += SubtractMatrices(&(ukf->SP_state[state_dims + i]), &ukf->linsys->x, &temp_x_vect); // temp_x = (SP_state[n+i] - x_prediction)
		error_code += MultiplyMatrices(&temp_x_vect, &temp_y_vect_T, &temp_np2); // temp_np2 = (SP_state[n+i] - x_prediction) * (SP_output[n+i] - y_prediction)^T
		error_code += AppendMatrixWithWeight(&temp_np1, &temp_np2, weight_i); // temp_np1 += w_i * (SP_state[n+i] - x_prediction) * (SP_output[n+i] - y_prediction)^T
	

	}
	
	error_code += AddMatrices(&temp_pp1, &ukf->R, &temp_pp1); // temp_pp1 += R
	error_code += InvertMatrix(&temp_pp1, &temp_pp2); // temp_pp2 = temp_pp1 ^ -1 ( == { Pyy + R }^-1 )

	error_code += MultiplyMatrices(&temp_np1, &temp_pp2, &ukf->G); // this is the kalman gain: G = Pxy * (Pyy + R)^-1

	//Calculate the x_filter
	error_code += SubtractMatrices(filter_input, &ukf->linsys->y, &temp_y_vect); // temp_y = (y - y_prediction)
	error_code += MultiplyMatrices(&ukf->G, &temp_y_vect, &temp_x_vect); // temp_x = G * (y - y_prediction)
	error_code += AddMatrices(&ukf->linsys->x, &temp_x_vect, &ukf->x_filter); // x_filter = x_prediction + G * (y - y_prediction)


	// Calculate P
	error_code += TransposeMatrix(&temp_np1, &temp_pn); // temp_pn = Pxy^T
	error_code += MultiplyMatrices(&ukf->G, &temp_pn, &temp_nn); // temp_nn = G * Pxy^T
	error_code += SubtractMatrices(&ukf->P_adv, &temp_nn, &ukf->P);

	// Calculate sigma points for next step
	AssignMatrix(ukf->SP_state, ukf->x_filter.data, state_dims);

	#ifdef USE_INCORRECT_MATRIX_SQRT_EXTENSION

	error_code += _UKF_AbusoluteMatrixSQRT(&ukf->P, &temp_nn); //temp_nn = SQRT(P)

	#else

	error_code += CholeskyDecomposition(&ukf->P, &temp_nn); //temp_nn = SQRT(P)

	#endif


	error_code += MultiplyMatrixByScalar(&temp_nn, k_sp, &temp_nn); // temp_nn = sqrt(n + k) * SQRT(P)


	for (i = 1; i <= state_dims; ++i)
	{
		error_code += GetRowVector(&temp_nn, &temp_x_vect, i - 1); // sqrt(n + k) * SQRT(P)[i]
		error_code += AddMatrices(&ukf->x_filter, &temp_x_vect, &(ukf->SP_state[i])); // SP_state[i] = x_filter + sqrt(n + k) * SQRT(P)[i] (i = 1, 2, ... , n)
		error_code += SubtractMatrices(&ukf->x_filter, &temp_x_vect, &(ukf->SP_state[i+state_dims]));  // SP_state[n+i] = x_filter - sqrt(n + k) * SQRT(P)[i] (i = 1, 2, ... , n)
	}

	return error_code;
}

//Free heap memory held by an object of Kalman Filter
void FreeKalman(Kalman_t* xkf)
{

	// if xkf is UKF
	if (xkf->SP_state)
	{
		// FreeMatrixArray(&(xkf->SP_state), 2 * xkf->x_filter.rows + 1);
		unsigned int i;
		for (i = 0; i < 2 * xkf->x_filter.rows + 1; ++i)
			FreeMatrix(&(xkf->SP_state[i]));
		free(xkf->SP_state);
		xkf->SP_state = NULL;
	}
	// if xkf is UKF
	if (xkf->SP_output)
	{
		// FreeMatrixArray(&(xkf->SP_output), 2 * xkf->x_filter.rows + 1);
		unsigned int i;
		for (i = 0; i < 2 * xkf->x_filter.rows + 1; ++i)
			FreeMatrix(&(xkf->SP_output[i]));
		free(xkf->SP_output);
		xkf->SP_output = NULL;
	}


	// common free process
	FreeMatrix(&xkf->P);
	FreeMatrix(&xkf->P_adv);
	FreeMatrix(&xkf->Q);
	FreeMatrix(&xkf->R);
	FreeMatrix(&xkf->G);
	FreeMatrix(&xkf->x_filter);

	FreeSSmodel(xkf->linsys);


	// if xkf is EKF or UKF
	if (xkf->state_model_func || xkf->output_model_func)
	{
		free(xkf->linsys);
		xkf->linsys = NULL;
	}

	xkf->state_model_func = NULL;
	xkf->output_model_func = NULL;
	xkf->state_model_jacobian_x = NULL;
	xkf->state_model_jacobian_u = NULL;
	xkf->output_model_jacobian_x = NULL;
	xkf->output_model_jacobian_u = NULL;
	
}
