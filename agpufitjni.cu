#include <jni.h>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <limits>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <cstdlib> // for Malloc function

/* --------------------------------------------------------------------------------------------------------
NOTES: 

This gpufit code is consolidated from an open-source code. Please see https://gpufit.readthedocs.io/en/latest/ for more information. 
All codes (headers and implementations), including CUDA codes, are deliberately placed together for ease of compilation 
using nvcc compiler and for ease of distribution. We added ACF and bleach correction fitting functions. Furthermore, 
we have also placed the CUDA ACF calculations codes here. 

Some implementation details:
1. void calcacf3 kernel function calculates the block transformation values of the intensity.
2. void calcacf2a kernel function calculates the arrays according to different time bins in different parts of the correlation function.
3. void calcacf2b kernel function calculates the value of the auto or cross-correlation at every lag time. This function also performs the G1 analysis in N and B calculation.
4. void calc_data_bleach_correction kernel function is an averaging step in temporal dimension for every ave number of points, prior to performing bleach correction fitting.
5. void calc_binning kernel function performs binning of spatial data.
6. void bleachcorrection kernel function performs polynomial bleach correction given polynomial order and coefficients. It is done prior to calcacf3, calcacf2a and calcacf2b.
7. Kindly also take a look at the CPU functions for a detailed description of the variables.
8. argument float* data is the intensity array input on which the N and B or the autocorrelation or the cross-correlation has to be calculates.
9. argument double* data1 is the output array where the values of auto and cross-correlation are calculated.
--------------------------------------------------------------------------------------------------------- */

// declare as a flag during compilation, ie. -D USE_CUBLAS
// Without USE_CUBLAS, we are using solve_equation_systems_gj(), i.e. cuda_gaussjordan.cu
// #define USE_CUBLAS 

// from gpufit.h
#ifdef __linux__
    #define VISIBLE __attribute__((visibility("default")))
#else
    #define VISIBLE
#endif

/* -------------------------------------------------------------------------------------------------------
* from definitions.h START
------------------------------------------------------------------------------------------------------- */
// Precision
#ifdef GPUFIT_DOUBLE
    #define REAL double
#else
    #define REAL float
#endif // GPUFIT_DOUBLE

#ifdef USE_CUBLAS
    #include "cublas_v2.h"

    #ifdef GPUFIT_DOUBLE
        #define DECOMPOSE_LUP cublasDgetrfBatched
        #define SOLVE_LUP cublasDgetrsBatched
    #else
        #define DECOMPOSE_LUP cublasSgetrfBatched
        #define SOLVE_LUP cublasSgetrsBatched
    #endif

    #define SOLVE_EQUATION_SYSTEMS() solve_equation_systems_lup()
#else
    #define cublasHandle_t int
    #define SOLVE_EQUATION_SYSTEMS() solve_equation_systems_gj()
#endif

// Status
#define CUDA_CHECK_STATUS( cuda_function_call ) \
    if (cudaError_t const status = cuda_function_call) \
    { \
        throw std::runtime_error( cudaGetErrorString( status ) ) ; \
    }
/* -------------------------------------------------------------------------------------------------------
* from definitions.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from gpufitImFCS_GpufitImFCS.h START
------------------------------------------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    fit
 * Signature: (IILjava/nio/FloatBuffer;Ljava/nio/FloatBuffer;ILjava/nio/FloatBuffer;FIILjava/nio/IntBuffer;IILjava/nio/FloatBuffer;Ljava/nio/FloatBuffer;Ljava/nio/IntBuffer;Ljava/nio/FloatBuffer;Ljava/nio/IntBuffer;)I
 */
JNIEXPORT jint JNICALL Java_gpufitImFCS_GpufitImFCS_fit
  (JNIEnv *, jclass, jint, jint, jobject, jobject, jint, jobject, jfloat, jint, jint, jobject, jint, jint, jobject, jobject, jobject, jobject, jobject);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    getLastError
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_gpufitImFCS_GpufitImFCS_getLastError
  (JNIEnv *, jclass);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    isCudaAvailableInt
 * Signature: ()Z
 */
JNIEXPORT jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isCudaAvailableInt
  (JNIEnv *, jclass);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    getCudaVersionAsArray
 * Signature: ()[I
 */
JNIEXPORT jintArray JNICALL Java_gpufitImFCS_GpufitImFCS_getCudaVersionAsArray
  (JNIEnv *, jclass);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    resetGPU
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_gpufitImFCS_GpufitImFCS_resetGPU
  (JNIEnv *, jclass);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcDataBleachCorrection
 * Signature: ([F[FLgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
JNIEXPORT void JNICALL Java_gpufitImFCS_GpufitImFCS_calcDataBleachCorrection
  (JNIEnv *, jclass, jfloatArray, jfloatArray, jobject);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    isBinningMemorySufficient
 * Signature: (LgpufitImFCS/GpufitImFCS/ACFParameters;)Z
 */
JNIEXPORT jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isBinningMemorySufficient
  (JNIEnv *, jclass, jobject);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcBinning
 * Signature: ([F[FLgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
JNIEXPORT void JNICALL Java_gpufitImFCS_GpufitImFCS_calcBinning
  (JNIEnv *, jclass, jfloatArray, jfloatArray, jobject);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    isACFmemorySufficient
 * Signature: (LgpufitImFCS/GpufitImFCS/ACFParameters;)Z
 */
JNIEXPORT jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isACFmemorySufficient
  (JNIEnv *, jclass, jobject);

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcACF
 * Signature: ([F[D[D[D[D[D[D[D[ILgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
JNIEXPORT void JNICALL Java_gpufitImFCS_GpufitImFCS_calcACF
  (JNIEnv *, jclass, jfloatArray, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jdoubleArray, jintArray, jobject);

#ifdef __cplusplus
}
#endif
/* -------------------------------------------------------------------------------------------------------
* from gpufitImFCS_GpufitImFCS.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from constants.h START
------------------------------------------------------------------------------------------------------- */
// fitting model ID

enum ModelID {
    GAUSS_2D = 1,
    ACF_1D= 2,
    LINEAR_1D = 3
};

// estimator ID
enum EstimatorID { LSE = 0, MLE = 1 };

// fit state
enum FitState { CONVERGED = 0, MAX_ITERATION = 1, SINGULAR_HESSIAN = 2, NEG_CURVATURE_MLE = 3, GPU_NOT_READY = 4 };

// return state
enum ReturnState { OK = 0, ERROR = -1 };

enum DataLocation { HOST = 0, DEVICE = 1 };
/* -------------------------------------------------------------------------------------------------------
* from constants.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from models/gauss_2d.cuh START
------------------------------------------------------------------------------------------------------- */
/* Description of the calculate_gauss2d function
* ==============================================
*
* This function calculates the values of two-dimensional gauss model functions
* and their partial derivatives with respect to the model parameters. 
*
* No independent variables are passed to this model function.  Hence, the 
* (X, Y) coordinate of the first data value is assumed to be (0.0, 0.0).  For
* a fit size of M x N data points, the (X, Y) coordinates of the data are
* simply the corresponding array index values of the data array, starting from
* zero.
*
* Parameters:
*
* parameters: An input vector of model parameters.
*             p[0]: amplitude
*             p[1]: center coordinate x
*             p[2]: center coordinate y
*             p[3]: width (standard deviation; equal width in x and y dimensions)
*             p[4]: offset
*
* n_fits: The number of fits. (not used)
*
* n_points: The number of data points per fit.
*
* value: An output vector of model function values.
*
* derivative: An output vector of model function partial derivatives.
*
* point_index: The data point index.
*
* fit_index: The fit index. (not used)
*
* chunk_index: The chunk index. (not used)
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The size of user_info in bytes. (not used)
*
* Calling the calculate_gauss2d function
* ======================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_gauss2d(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size,
    int const num_v_coefs) // NEW
{
    // indices

    int const n_points_x = sqrt((REAL)n_points);
    int const point_index_y = point_index / n_points_x;
    int const point_index_x = point_index - point_index_y * n_points_x;

    // parameters

    REAL const * p = parameters;

    // value

    REAL const argx = (point_index_x - p[1]) * (point_index_x - p[1]) / (2 * p[3] * p[3]);
    REAL const argy = (point_index_y - p[2]) * (point_index_y - p[2]) / (2 * p[3] * p[3]);
    REAL const ex = exp(-(argx + argy));
    value[point_index] = p[0] * ex + p[4];

    // derivatives

    REAL * current_derivative = derivative + point_index;

    current_derivative[0 * n_points] = ex;
    current_derivative[1 * n_points] = p[0] * ex * (point_index_x - p[1]) / (p[3] * p[3]);
    current_derivative[2 * n_points] = p[0] * ex * (point_index_y - p[2]) / (p[3] * p[3]);
    current_derivative[3 * n_points] = ex * p[0] * ((point_index_x - p[1]) * (point_index_x - p[1]) + (point_index_y - p[2]) * (point_index_y - p[2])) / (p[3] * p[3] * p[3]);
    current_derivative[4 * n_points] = 1;
}
/* -------------------------------------------------------------------------------------------------------
* from models/gauss_2d.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from models/linear_1d.cuh START
------------------------------------------------------------------------------------------------------- */
/* Description of the calculate_linear1d function
* ===================================================
*
* This function calculates the values of one-dimensional linear model functions
* and their partial derivatives with respect to the model parameters. 
*
* This function makes use of the user information data to pass in the 
* independent variables (X values) corresponding to the data.  The X values
* must be of type REAL.
*
* Note that if no user information is provided, the (X) coordinate of the 
* first data value is assumed to be (0.0).  In this case, for a fit size of 
* M data points, the (X) coordinates of the data are simply the corresponding 
* array index values of the data array, starting from zero.
*
* There are three possibilities regarding the X values:
*
*   No X values provided: 
*
*       If no user information is provided, the (X) coordinate of the 
*       first data value is assumed to be (0.0).  In this case, for a 
*       fit size of M data points, the (X) coordinates of the data are 
*       simply the corresponding array index values of the data array, 
*       starting from zero.
*
*   X values provided for one fit:
*
*       If the user_info array contains the X values for one fit, then 
*       the same X values will be used for all fits.  In this case, the 
*       size of the user_info array (in bytes) must equal 
*       sizeof(REAL) * n_points.
*
*   Unique X values provided for all fits:
*
*       In this case, the user_info array must contain X values for each
*       fit in the dataset.  In this case, the size of the user_info array 
*       (in bytes) must equal sizeof(REAL) * n_points * nfits.
*
* Parameters:
*
* parameters: An input vector of model parameters.
*             p[0]: offset
*             p[1]: slope
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* value: An output vector of model function values.
*
* derivative: An output vector of model function partial derivatives.
*
* point_index: The data point index.
*
* fit_index: The fit index.
*
* chunk_index: The chunk index. Used for indexing of user_info.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the calculate_linear1d function
* =======================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_linear1d(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size,
    int const num_v_coefs) // NEW
{
    // indices

    REAL * user_info_float = (REAL*) user_info;
    REAL x = 0;
  //  if (!user_info_float)
  //  {
  //      x = point_index;
  //  }
  //  else if (user_info_size / sizeof(REAL) == n_points)
  //  {
        x = user_info_float[point_index];
   // }
   // else if (user_info_size / sizeof(REAL) > n_points)
  //  {
  //      int const chunk_begin = chunk_index * n_fits * n_points;
  //      int const fit_begin = fit_index * n_points;
  //      x = user_info_float[chunk_begin + fit_begin + point_index];
  //  }

    // value
    //value[point_index] = parameters[0] + parameters[1]*x + parameters[2]*x*x + parameters[3]*x*x*x + parameters[4]*x*x*x*x + parameters[5]*x*x*x*x*x + parameters[6]*x*x*x*x*x*x;
    value[point_index] = parameters[0];    

    // derivatives
    REAL * current_derivatives = derivative + point_index;
    current_derivatives[0 * n_points] = 1;

/*
    REAL tempx = 1.0;
    for (int i = 1; i < num_v_coefs; i++) {
        tempx *= x;
        current_derivatives[i * n_points] = tempx;   
	value[point_index] =value[point_index]+ parameters[i]*tempx;     
    }
*/

    double tempx = 1.0;
    double sumval = (double) parameters[0];
    for (int i = 1; i < num_v_coefs; i++) {
        tempx *= x;
        current_derivatives[i * n_points] = (float) tempx;   
        sumval += (double) parameters[i]*tempx;
    }
    value[point_index] =(float) sumval;     

}
/* -------------------------------------------------------------------------------------------------------
* from models/linear_1d.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from models/acf_1d.cuh START
------------------------------------------------------------------------------------------------------- */
__device__ void calculate_acf1d(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size,
    int const num_v_coefs) // NEW
{
    // indices

    REAL * user_info_float = (REAL*) user_info;
    REAL x = 0;
    x = user_info_float[point_index];
    double sqrpi = sqrt((double) 3.14159265359);

    double p0t = sqrt(4 * (double) parameters[1] * x + pow((double) parameters[13], 2.0));
    double p1xt = (double) parameters[11] + (double) parameters[15] - (double) parameters[2] * x;
    double p2xt = (double) parameters[11] - (double) parameters[15] + (double) parameters[2] * x;
    double p3xt = (double) parameters[15] - (double) parameters[2] * x;
    double p4xt = 2 * pow((double) parameters[11], 2.0) + 3 * pow((double) parameters[15], 2.0) - 6 * x * (double) parameters[15] * (double) parameters[2] + 3 * pow(x * (double) parameters[2], 2.0);
    double p5xt = pow(p3xt, 2.0) + pow(p1xt, 2.0);
    double p6xt = pow(p3xt, 2.0) + pow(p2xt, 2.0);
    double p7xt = 2 * (pow((double) parameters[11], 2.0) + pow((double) parameters[15], 2.0) - 2 * x * (double) parameters[15] * (double) parameters[2] + pow(x * (double) parameters[2], 2.0));
    double p1yt = (double) parameters[12] + (double) parameters[16] - (double) parameters[3] * x;
    double p2yt = (double) parameters[12] - (double) parameters[16] + (double) parameters[3] * x;
    double p3yt = (double) parameters[16] - (double) parameters[3] * x;
    double p4yt = 2 * pow((double) parameters[12], 2.0) + 3 * pow((double) parameters[16], 2.0) - 6 * x * (double) parameters[16] * (double) parameters[3] + 3 * pow(x * (double) parameters[3], 2.0);
    double p5yt = pow(p3yt, 2.0) + pow(p1yt, 2.0);
    double p6yt = pow(p3yt, 2.0) + pow(p2yt, 2.0);
    double p7yt = 2 * (pow((double) parameters[12], 2.0) + pow((double) parameters[16], 2.0) - 2 * x * (double) parameters[16] * (double) parameters[3] + pow(x * (double) parameters[3], 2.0));

    double pexpxt = exp(-pow(p1xt / p0t, 2.0)) + exp(-pow(p2xt / p0t, 2.0)) - 2 * exp(-pow(p3xt / p0t, 2.0));
    double perfxt = p1xt * erf(p1xt / p0t) + p2xt * erf(p2xt / p0t) - 2 * p3xt * erf(p3xt / p0t);
    double dDpexpxt = 2 * exp(-p4xt / pow(p0t, 2.0)) * (exp(p5xt / pow(p0t, 2.0)) + exp(p6xt / pow(p0t, 2.0)) - 2 * exp(p7xt / pow(p0t, 2.0)));
    double dvxperfxt = (erf(p2xt / p0t) + 2 * erf(p3xt / p0t) - erf(p1xt / p0t)) * x;
    double pexpyt = exp(-pow(p1yt / p0t, 2.0)) + exp(-pow(p2yt / p0t, 2.0)) - 2 * exp(-pow(p3yt / p0t, 2.0));
    double dDpexpyt = 2 * exp(-p4yt / pow(p0t, 2.0)) * (exp(p5yt / pow(p0t, 2.0)) + exp(p6yt / pow(p0t, 2.0)) - 2 * exp(p7yt / pow(p0t, 2.0)));
    double dvyperfyt = (erf(p2yt / p0t) + 2 * erf(p3yt / p0t) - erf(p1yt / p0t)) * x;
    double perfyt = p1yt * erf(p1yt / p0t) + p2yt * erf(p2yt / p0t) - 2 * p3yt * erf(p3yt / p0t);
    double pplane1 = (p0t / sqrpi * pexpxt + perfxt) * (p0t / sqrpi * pexpyt + perfyt) / (4 * (double) parameters[11]*(double) parameters[12]) * ((double) parameters[17] / ((double) parameters[11]*(double) parameters[12]));
    double pspim1 = 1 / sqrt(1 + (4 * (double) parameters[1] * x) / powf((double) parameters[14], 2));
    double acf1 = pplane1 * pspim1;

    double p0t2 = sqrt(4 * (double) parameters[6] * x + pow((double) parameters[13], 2));
    double p1xt2 = (double) parameters[11] + (double) parameters[15] - (double) parameters[2] * x;
    double p2xt2 = (double) parameters[11] - (double) parameters[15] + (double) parameters[2] * x;
    double p3xt2 = (double) parameters[15] - (double) parameters[2] * x;
    double p4xt2 = 2 * pow((double) parameters[11], 2) + 3 * pow((double) parameters[15], 2) - 6 * x * (double) parameters[15] * (double) parameters[2] + 3 * pow(x * (double) parameters[2], 2);
    double p5xt2 = pow(p3xt2, 2) + pow(p1xt2, 2);
    double p6xt2 = pow(p3xt2, 2) + pow(p2xt2, 2);
    double p7xt2 = 2 * (pow((double) parameters[11], 2) + pow((double) parameters[15], 2) - 2 * x * (double) parameters[15] * (double) parameters[2] + pow(x * (double) parameters[2], 2));
    double p1yt2 = (double) parameters[12] + (double) parameters[16] - (double) parameters[3] * x;
    double p2yt2 = (double) parameters[12] - (double) parameters[16] + (double) parameters[3] * x;
    double p3yt2 = (double) parameters[16] - (double) parameters[3] * x;
    double p4yt2 = 2 * pow((double) parameters[12], 2) + 3 * pow((double) parameters[16], 2) - 6 * x * (double) parameters[16] * (double) parameters[3] + 3 * pow(x * (double) parameters[3], 2);
    double p5yt2 = pow(p3yt2, 2) + pow(p1yt2, 2);
    double p6yt2 = pow(p3yt2, 2) + pow(p2yt2, 2);
    double p7yt2 = 2 * (pow((double) parameters[12], 2) + pow((double) parameters[16], 2) - 2 * x * (double) parameters[16] * (double) parameters[3] + pow(x * (double) parameters[3], 2));
    double pexpxt2 = exp(-pow(p1xt2 / p0t2, 2)) + exp(-pow(p2xt2 / p0t2, 2)) - 2 * exp(-pow(p3xt2 / p0t2, 2));
    double perfxt2 = p1xt2 * erf(p1xt2 / p0t2) + p2xt2 * erf(p2xt2 / p0t2) - 2 * p3xt2 * erf(p3xt2 / p0t2);
    double dDpexpxt2 = 2 * exp(-p4xt2 / pow(p0t2, 2)) * (exp(p5xt2 / pow(p0t2, 2)) + exp(p6xt2 / pow(p0t2, 2)) - 2 * exp(p7xt2 / pow(p0t2, 2)));
    double dvxperfxt2 = (erf(p2xt2 / p0t2) + 2 * erf(p3xt2 / p0t2) - erf(p1xt2 / p0t2)) * x;
    double pexpyt2 = exp(-pow(p1yt2 / p0t2, 2)) + exp(-pow(p2yt2 / p0t2, 2)) - 2 * exp(-pow(p3yt2 / p0t2, 2));
    double dDpexpyt2 = 2 * exp(-p4yt2 / pow(p0t2, 2)) * (exp(p5yt2 / pow(p0t2, 2)) + exp(p6yt2 / powf(p0t2, 2)) - 2 * exp(p7yt2 / powf(p0t2, 2)));
    double dvyperfyt2 = (erf(p2yt2 / p0t2) + 2 * erf(p3yt2 / p0t2) - erf(p1yt2 / p0t2)) * x;
    double perfyt2 = p1yt2 * erf(p1yt2 / p0t2) + p2yt2 * erf(p2yt2 / p0t2) - 2 * p3yt2 * erf(p3yt2 / p0t2);
    double pplane2 = (p0t2 / sqrpi * pexpxt2 + perfxt2) * (p0t2 / sqrpi * pexpyt2 + perfyt2) / (4 * pow((double) parameters[11]*(double) parameters[12], 2) / (double) parameters[17]);
    double pspim2 = 1 / sqrt(1 + (4 * (double) parameters[6] * x) / pow((double) parameters[14], 2));
    double acf2 = pplane2 * pspim2;

    double triplet = 1 + (double) parameters[9] / (1 - (double) parameters[9]) * exp(-x / (double) parameters[10]);

    value[point_index] = (1 / (double) parameters[0]) * ((1 - (double) parameters[5]) * acf1 + powf((double) parameters[18], 2) * (double) parameters[5] * acf2) / pow(1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5], 2) * triplet + (double) parameters[4];

    double dDplat = (1 / (sqrpi * p0t)) * (dDpexpyt * x * (p0t / sqrpi * pexpxt + perfxt) + dDpexpxt * x * (p0t / sqrpi * pexpyt + perfyt)) / (4 * powf((double) parameters[11]*(double) parameters[12], 2.0) / (double) parameters[17]);
    double dDpspim = -4 * x / (2 * pow((double) parameters[14], 2) * pow(sqrt(1 + (4 * (double) parameters[1] * x) / pow((double) parameters[14], 2)), 3));

    double dDplat2 = (1 / (sqrpi * p0t2)) * (dDpexpyt2 * x * (p0t2 / sqrpi * pexpxt2 + perfxt2) + dDpexpxt2 * x * (p0t2 / sqrpi * pexpyt2 + perfyt2)) / (4 * pow((double) parameters[11]*(double) parameters[12], 2) / (double) parameters[17]);
    double dDpspim2 = -4 * x / (2 * pow((double) parameters[14], 2) * pow(sqrt(1 + (4 * (double) parameters[6] * x) / pow((double) parameters[14], 2)), 3));

    double dtripletFtrip = exp(-x / (double) parameters[10]) * (1 / (1 - (double) parameters[9]) + (double) parameters[9] / pow(1 - (double) parameters[9], 2));
    double dtripletTtrip = exp(-x / (double) parameters[10]) * ((double) parameters[9] * x) / ((1 - (double) parameters[9]) * pow((double) parameters[10], 2));

    double pf1 = (1 - (double) parameters[5]) / (1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5]);
    double pf2 = (pow((double) parameters[18], 2) * (double) parameters[5]) / (1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5]);
    double dfnom = pow(1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5], 3);
    double df21 = 1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5] - 2 * (double) parameters[18];
    double df22 = pow((double) parameters[18], 2) * (1 + (double) parameters[5] - (double) parameters[18] * (double) parameters[5]);

    double pacf = (1 / (double) parameters[0]) * ((1 - (double) parameters[5]) * acf1 + powf((double) parameters[18], 2) * (double) parameters[5] * acf2) / pow(1 - (double) parameters[5] + (double) parameters[18] * (double) parameters[5], 2) * triplet + (double) parameters[4];

    REAL * current_derivatives = derivative + point_index;

    current_derivatives[0 * n_points] = (float) (-1 / pow((double)parameters[0], 2)) * (pf1 * acf1 + pf2 * acf2) * triplet;
    current_derivatives[1 * n_points] = (1 /  parameters[0]) * (float)(pf1 *  (pplane1 * dDpspim + pspim1 * dDplat));
    current_derivatives[2 * n_points] = (1 /   parameters[0]) * (float)((pf1 * ((p0t / sqrpi * pexpyt + perfyt) * dvxperfxt) * pspim1 / (4 * pow((double) parameters[11]*(double) parameters[12], 2) /  parameters[17]) + pf2 * ((p0t2 / sqrpi * pexpyt2 + perfyt2) * dvxperfxt2) * pspim2 / (4 * pow((double) parameters[11]*(double) parameters[12], 2) /  parameters[17])) * triplet);
    current_derivatives[3 * n_points] = (1 /   parameters[0]) * (float)((pf1 * ((p0t / sqrpi * pexpxt + perfxt) * dvyperfyt) * pspim1 / (4 * pow((double) parameters[11]*(double) parameters[12], 2) /   parameters[17]) + pf2 * ((p0t2 / sqrpi * pexpxt2 + perfxt2) * dvyperfyt2) * pspim2 / (4 * pow((double) parameters[11]*(double) parameters[12], 2) /   parameters[17])) * triplet);
    current_derivatives[4 * n_points] = 1.0;
    current_derivatives[5 * n_points] = (1 /  parameters[0]) *(float)((1 / dfnom) * (df21 * acf1 + df22 * acf2) * triplet);
    current_derivatives[6 * n_points] = (1 /  parameters[0]) * (float)(pf2 * (pplane2 * dDpspim2 + pspim2 * dDplat2) * triplet);
    current_derivatives[9 * n_points] = (float)dtripletFtrip * pacf;
    current_derivatives[10 * n_points] = (float)dtripletTtrip * pacf;
    //current_derivatives[11 * n_points] = 0.0;
    //current_derivatives[12 * n_points] = 0.0;
    //current_derivatives[13 * n_points] = 0.0;
    //current_derivatives[14 * n_points] = 0.0;
    //current_derivatives[15 * n_points] = 0.0;
    //current_derivatives[16 * n_points] = 0.0;
    //current_derivatives[17 * n_points] = 0.0;
    //current_derivatives[18 * n_points] = 0.0;
    //current_derivatives[19 * n_points] = 0.0; 
}
/* -------------------------------------------------------------------------------------------------------
* from models/acf_1d.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from models/models.cuh START
------------------------------------------------------------------------------------------------------- */
__device__ void calculate_model(
    ModelID const model_id,
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    int const user_info_size,
    int const num_v_coefs) // NEW
{
    switch (model_id)
    {
    case GAUSS_2D:
        calculate_gauss2d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size, num_v_coefs);
        break;
    case ACF_1D:
	calculate_acf1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size, num_v_coefs);
	break;
    case LINEAR_1D:
        calculate_linear1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size, num_v_coefs);
        break;
    default:
        break;
    }
}

void configure_model(ModelID const model_id, int & n_parameters, int & n_dimensions)
{
    switch (model_id)
    {
    case GAUSS_2D:              n_parameters = 5;  n_dimensions = 2;  break;
    case ACF_1D:                n_parameters = 20; n_dimensions = 1;  break;
    case LINEAR_1D:             n_parameters = 11; n_dimensions = 1;  break;
    default:                                                          break;
    }
}
/* -------------------------------------------------------------------------------------------------------
* from models/models.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from estimators/lse.cuh START
------------------------------------------------------------------------------------------------------- */
/* Description of the calculate_chi_square_lse function
* =====================================================
*
* This function calculates the chi-square values for the weighted LSE estimator.
*
* Parameters:
*
* chi_square: An output vector of chi-square values for each data point.
*
* point_index: The data point index.
*
* data: An input vector of data values.
*
* value: An input vector of fitting curve values.
*
* weight: An optional input vector of values for weighting the chi-square values.
*
* state: A pointer to a value which indicates whether the fitting
*        process was carreid out correctly or which problem occurred.
*        In this function it is not used. It can be used in functions calculating
*        other estimators than the LSE, such as MLE. It is passed into this function
*        to provide the same interface for all estimator functions.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_chi_square_lse function
* =============================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_chi_square_lse(
    volatile REAL * chi_square,
    int const point_index,
    REAL const * data,
    REAL const * value,
    REAL const * weight,
    int * state,
    char * user_info,
    std::size_t const user_info_size)
{
    REAL const deviation = value[point_index] - data[point_index];

    if (weight)
    {
        chi_square[point_index] = deviation * deviation * weight[point_index];
    }
    else
    {
        chi_square[point_index] = deviation * deviation;
    }
}

/* Description of the calculate_hessian_lse function
* ==================================================
*
* This function calculates the hessian matrix values of the weighted LSE estimator.
* The calculation is performed based on previously calculated fitting curve derivative
* values.
*
* Parameters:
*
* hessian: An output vector of values of the hessian matrix for each data point.
*
* point_index: The data point index.
*
* parameter_index_i: Index of the hessian column.
*
* parameter_index_j: Index of the hessian row.
*
* data: An input vector of data values. (not used)
*
* value: An input vector of fitting curve values. (not used)
*
* derivative: An input vector of partial derivative values of the fitting
*             curve with respect to the fitting parameters for each data point.
*
* weight: An optional input vector of values for weighting the hessian matrix values.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_hessian_lse function
* ==========================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_hessian_lse(
    double * hessian,
    int const point_index,
    int const parameter_index_i,
    int const parameter_index_j,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    if (weight)
    {
        *hessian
            += derivative[parameter_index_i] * derivative[parameter_index_j]
            * weight[point_index];
    }
    else
    {
        *hessian
            += derivative[parameter_index_i] * derivative[parameter_index_j];
    }
}

/* Description of the calculate_gradient_lse function
* ===================================================
*
* This function calculates the gradient values of the weighted LSE estimator
* based on previously calculated fitting curve derivative values.
*
* Parameters:
*
* gradient: An output vector of values of the gradient vector for each data point.
*
* point_index: The data point index.
*
* parameter_index: The parameter index.
*
* n_parameters: The number of fitting curve parameters.
*
* data: An input vector of data values.
*
* value: An input vector of fitting curve values.
*
* derivative: An input vector of partial derivative values of the fitting
*             curve with respect to the fitting parameters for each data point.
*
* weight: An optional input vector of values for weighting gradient values.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_gradient_lse function
* ===========================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_gradient_lse(
    volatile REAL * gradient,
    int const point_index,
    int const parameter_index,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    REAL const deviation = data[point_index] - value[point_index];

    if (weight)
    {
        gradient[point_index]
            = derivative[parameter_index] * deviation * weight[point_index];
    }
    else
    {
        gradient[point_index]
            = derivative[parameter_index] * deviation;
    }
}
/* -------------------------------------------------------------------------------------------------------
* from estimators/lse.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from estimators/mle.cuh START
------------------------------------------------------------------------------------------------------- */
/* Description of the calculate_chi_square_mle function
* =====================================================
*
* This function calculates the chi-square values for the MLE estimator.
*
* Parameters:
*
* chi_square: An output vector of chi-square values for each data point.
*
* point_index: The data point index.
*
* data: An input vector of data.
*
* value: An input vector of fitting curve values.
*
* weight: An input vector of values for weighting chi-square values. It is not used
*         in this function. It can be used in functions calculating other estimators
*         than the MLE, such as LSE.
*
* state: A pointer to a value which indicates whether the fitting process was carreid
*        out correctly or which problem occurred. It is set to 3 if a fitting curve
*        value is negative.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_chi_square_mle function
* =============================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_chi_square_mle(
    volatile REAL * chi_square,
    int const point_index,
    REAL const * data,
    REAL const * value,
    REAL const * weight,
    int * state,
    char * user_info,
    std::size_t const user_info_size)
{
    if (value[point_index] < 0)
    {
        *state = 3;
    }

    REAL const deviation = value[point_index] - data[point_index];

    if (data[point_index] != 0)
    {
        chi_square[point_index]
            = 2 * (deviation - data[point_index] * std::log(value[point_index] / data[point_index]));
    }
    else
    {
        chi_square[point_index] = 2 * deviation;
    }
}

/* Description of the calculate_hessian_mle function
* ==================================================
*
* This function calculates the hessian matrix values of the MLE estimator. The
* calculation is performed based on previously calculated derivative values.
* 
* Parameters:
*
* hessian: An output vector of values of the hessian matrix for each data point.
*
* point_index: The data point index.
*
* parameter_index_i: Index of the hessian column.
*
* parameter_index_j: Index of the hessian row.
*
* data: An input vector of data values.
*
* value: An input vector of fitting curve values.
*
* derivative: An input vector of partial derivative values of the fitting
*             curve with respect to the fitting parameters for each data point.
*
* weight: An input vector of values for weighting hessian matrix values. It is not
*         used in this function. It can be used in functions calculating other estimators
*         than the MLE, such as LSE.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_hessian_mle function
* ==========================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_hessian_mle(
    double * hessian,
    int const point_index,
    int const parameter_index_i,
    int const parameter_index_j,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    *hessian
        += data[point_index]
        / (value[point_index] * value[point_index])
        * derivative[parameter_index_i] * derivative[parameter_index_j];
}

/* Description of the calculate_gradient_mle function
* ===================================================
*
* This function calculates the gradient values of the MLE estimator based
* on previously calculated derivative values.
*
* Parameters:
*
* gradient: An output vector of values of the gradient vector for each data point.
*
* point_index: The data point index.
*
* parameter_index: The parameter index.
*
* data: An input vector of data values.
*
* value: An input vector of fitting curve values.
*
* derivative: An input vector of partial derivative values of the fitting
*             curve with respect to the fitting parameters for each data point.
*
* weight: An input vector of values for weighting gradient vector values. It is not
*         used in this function. It can be used in functions calculating other estimators
*         than the MLE, such as LSE.
*
* user_info: An input vector containing user information. (not used)
*
* user_info_size: The number of elements in user_info. (not used)
*
* Calling the calculate_gradient_mle function
* ===========================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
*
*/

__device__ void calculate_gradient_mle(
    volatile REAL * gradient,
    int const point_index,
    int const parameter_index,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    gradient[point_index]
        = -derivative[parameter_index]
        * (1 - data[point_index] / value[point_index]);
}
/* -------------------------------------------------------------------------------------------------------
* from estimators/mle.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from estimators/estimators.cuh START
------------------------------------------------------------------------------------------------------- */
__device__ void calculate_chi_square(
    int const estimator_id,
    volatile REAL * chi_square,
    int const point_index,
    REAL const * data,
    REAL const * value,
    REAL const * weight,
    int * state,
    char * user_info,
    std::size_t const user_info_size)
{
    switch (estimator_id)
    {
    case LSE:
        calculate_chi_square_lse(chi_square, point_index, data, value, weight, state, user_info, user_info_size);
        break;
    case MLE:
        calculate_chi_square_mle(chi_square, point_index, data, value, weight, state, user_info, user_info_size);
        break;
    default:
        break;
    }
}

__device__ void calculate_gradient(
    int const estimator_id,
    volatile REAL * gradient,
    int const point_index,
    int const parameter_index,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    switch (estimator_id)
    {
    case LSE:
        calculate_gradient_lse(gradient, point_index, parameter_index, data, value, derivative, weight, user_info, user_info_size);
        break;
    case MLE:
        calculate_gradient_mle(gradient, point_index, parameter_index, data, value, derivative, weight, user_info, user_info_size);
        break;
    default:
        break;
    }
}

__device__ void calculate_hessian(
    int const estimator_id,
    double * hessian,
    int const point_index,
    int const parameter_index_i,
    int const parameter_index_j,
    REAL const * data,
    REAL const * value,
    REAL const * derivative,
    REAL const * weight,
    char * user_info,
    std::size_t const user_info_size)
{
    switch (estimator_id)
    {
    case LSE:
        calculate_hessian_lse
        (hessian, point_index, parameter_index_i, parameter_index_j, data, value, derivative, weight, user_info,user_info_size);
        break;
    case MLE:
        calculate_hessian_mle
        (hessian, point_index, parameter_index_i, parameter_index_j, data, value, derivative, weight, user_info, user_info_size);
        break;
    default:
        break;
    }
}
/* -------------------------------------------------------------------------------------------------------
* from estimators/estimators.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from cuda_kernel.cuh START
------------------------------------------------------------------------------------------------------- */
/*
void configure_model(ModelID const model_id, int & n_parameters, int & n_dimensions);

extern __global__ void cuda_sum_chi_square_subtotals(
    REAL * chi_squares,
    REAL const * subtotals,
    int const n_blocks_per_fit,
    int const n_fits,
    int const * finished);

extern __global__ void cuda_check_fit_improvement(
    int * iteration_failed,
    REAL const * chi_squares,
    REAL const * prev_chi_squares,
    int const n_fits,
    int const * finished);

extern __global__ void cuda_calculate_chi_squares(
    REAL * chi_squares,
    int * states,
    REAL const * data,
    REAL const * values,
    REAL const * weights,
    int const n_points,
    int const n_fits,
    int const estimator_id,
    int const * finished,
    int const n_fits_per_block,
    char * user_info,
    std::size_t const user_info_size);

extern __global__ void cuda_sum_gradient_subtotals(
    REAL * gradients,
    REAL const * subtotals,
    int const n_blocks_per_fit,
    int const n_fits,
    int const n_parameters,
    int const * skip,
    int const * finished);

extern __global__ void cuda_calculate_gradients(
    REAL * gradients,
    REAL const * data,
    REAL const * values,
    REAL const * derivatives,
    REAL const * weights,
    int const n_points,
    int const n_fits,
    int const n_parameters,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
	int const estimator_id,
    int const * finished,
    int const * skip,
    int const n_fits_per_block,
    char * user_info,
    std::size_t const user_info_size);

extern __global__ void cuda_calculate_hessians(
    REAL * hessians,
    REAL const * data,
    REAL const * values,
    REAL const * derivatives,
    REAL const * weights,
    int const n_fits,
    int const n_points,
    int const n_parameters,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
    int const estimator_id,
    int const * skip,
    int const * finished,
    char * user_info,
    std::size_t const user_info_size);

extern __global__ void cuda_modify_step_widths(
    REAL * hessians,
    REAL const * lambdas,
    REAL * scaling_vectors,
    unsigned int const n_parameters,
    int const * iteration_failed,
    int const * finished,
    int const n_fits_per_block);

extern __global__ void cuda_calc_curve_values(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    int const n_parameters,
    int const * finished,
    REAL * values,
    REAL * derivatives,
    int const n_fits_per_block,
    int const n_blocks_per_fit,
    ModelID const model_id,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size,
    int const num_v_coefs); // NEW

extern __global__ void cuda_update_parameters(
    REAL * parameters,
    REAL * prev_parameters,
    REAL const * deltas,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
    int const * finished,
    int const n_fits_per_block);

extern __global__ void cuda_check_for_convergence(
    int * finished,
    REAL const tolerance,
    int * states,
    REAL const * chi_squares,
    REAL const * prev_chi_squares,
    int const iteration,
    int const max_n_iterations,
    int const n_fits);

extern __global__ void cuda_evaluate_iteration(
    int * all_finished,
    int * n_iterations,
    int * finished,
    int const iteration,
    int const * states,
    int const n_fits);

extern __global__ void cuda_prepare_next_iteration(
    REAL * lambdas,
    REAL * chi_squares,
    REAL * prev_chi_squares,
    REAL * function_parameters,
    REAL const * prev_parameters,
    int const n_fits,
    int const n_parameters);

extern __global__ void cuda_update_state_after_solving(
    int const n_fits,
    int const * singular_checks,
    int const * finished,
    int * states);
*/
/* -------------------------------------------------------------------------------------------------------
* from cuda_kernel.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from cuda_kernel.cu START
------------------------------------------------------------------------------------------------------- */
/* Description of the cuda_calc_curve_values function
* ===================================================
*
* This function calls one of the fitting curve functions depending on the input
* parameter model_id. The fitting curve function calculates the values of
* the fitting curves and its partial derivatives with respect to the fitting
* curve parameters. Multiple fits are calculated in parallel.
*
* Parameters:
*
* parameters: An input vector of concatenated sets of model parameters.
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* n_parameters: The number of curve parameters.
*
* finished: An input vector which allows the calculation to be skipped for single
*           fits.
*
* values: An output vector of concatenated sets of model function values.
*
* derivatives: An output vector of concatenated sets of model function partial
*              derivatives.
*
* n_fits_per_block: The number of fits calculated by each thread block.
*
* n_blocks_per_fit: The number of thread blocks used to calculate one fit.
*
* model_id: The fitting model ID.
*
* chunk_index: The data chunk index.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the cuda_calc_curve_values function
* ===========================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = n_points * n_fits_per_block / n_blocks_per_fit;
*   blocks.x = n_fits / n_fits_per_block * n_blocks_per_fit;
*
*   cuda_calc_curve_values<<< blocks, threads >>>(
*       parameters,
*       n_fits,
*       n_points,
*       n_parameters,
*       finished,
*       values,
*       derivatives,
*       n_fits_per_block,
*       n_blocks_per_fit,
*       model_id,
*       chunk_index,
*       user_info,
*       user_info_size,
        num_v_coefs);
*
*/

__global__ void cuda_calc_curve_values(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    int const n_parameters,
    int const * finished,
    REAL * values,
    REAL * derivatives,
    int const n_fits_per_block,
    int const n_blocks_per_fit,
    ModelID const model_id,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size,
    int const num_v_coefs) // NEW
{
    int const fit_in_block = threadIdx.x / n_points;
    int const fit_index = blockIdx.x * n_fits_per_block / n_blocks_per_fit + fit_in_block;
    int const fit_piece = blockIdx.x % n_blocks_per_fit;
    int const point_index = threadIdx.x - fit_in_block * n_points + fit_piece * blockDim.x;
    int const first_point = fit_index * n_points;

    REAL * current_values = values + first_point;
    REAL * current_derivatives = derivatives + first_point * n_parameters;
    REAL const * current_parameters = parameters + fit_index * n_parameters;

    if (finished[fit_index])
        return;
    if (point_index >= n_points)
        return;

    calculate_model(model_id, current_parameters, n_fits, n_points, current_values, current_derivatives, point_index, fit_index, chunk_index, user_info, user_info_size, num_v_coefs);
}

/* Description of the sum_up_floats function
* ==========================================
*
* This function sums up a vector of REAL values and stores the result at the
* first place of the vector.
*
* Parameters:
*
* shared_array: An input vector of REAL values. The vector must be stored
*               on the shared memory of the GPU. The size of this vector must be a
*               power of two. Use zero padding to extend it to the next highest
*               power of 2 greater than the number of elements.
*
* size: The number of elements in the input vector considering zero padding.
*
* Calling the sum_up_floats function
* ==================================
*
* This __device__ function can be only called from a __global__ function or
* an other __device__ function. When calling the function, the blocks and threads
* of the __global__ function must be set up correctly, as shown in the following
* example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = size * vectors_per_block;
*   blocks.x = n_vectors / vectors_per_block;
*
*   global_function<<< blocks, threads >>>(parameter1, ...);
*
*/

__device__ void sum_up_floats(volatile REAL* shared_array, int const size)
{
    int const fit_in_block = threadIdx.x / size;
    int const point_index = threadIdx.x - (fit_in_block*size);

    int current_n_points = size >> 1;
    __syncthreads();
    while (current_n_points)
    {
        if (point_index < current_n_points)
        {
            shared_array[point_index] += shared_array[point_index + current_n_points];
        }
        current_n_points >>= 1;
        __syncthreads();
    }
}

/* Description of the cuda_sum_chi_square_subtotals function
* ==========================================================
*
* This function sums up chi_square subtotals in place.
*
* Parameters:
*
* chi_squares: A vector of chi-square values for multiple fits.
*              in: subtotals
*              out: totals
*
* n_blocks_per_fit: The number of blocks used to calculate one fit. It is
*                   equivalent to the number of subtotals per fit.
*
* n_fits: The number of fits.
*
* finished: An input vector which allows the calculation to be skipped
*           for single fits.
*
* Calling the cuda_sum_chi_square_subtotals function
* ==================================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_sum_chi_square_subtotals<<< blocks, threads >>>(
*       chi_squares,
*       n_blocks_per_fit,
*       n_fits,
*       finished);
*
*/

__global__ void cuda_sum_chi_square_subtotals(
    REAL * chi_squares,
    REAL const * subtotals,
    int const n_blocks_per_fit,
    int const n_fits,
    int const * finished)
{
    int const index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= n_fits || finished[index])
        return;

    REAL * chi_square = chi_squares + index;
    REAL const * subtotal = subtotals + index;

    double sum = 0.0;
    for (int i = 0; i < n_blocks_per_fit; i++)
        sum += subtotal[i * n_fits];

    chi_square[0] = sum;
}

/* Description of the cuda_check_fit_improvement function
* =======================================================
*
* This function checks after each calculation of chi-square values whether the
* currently calculated chi-square values are lower than chi-square values calculated
* in the previous iteration and sets the iteration_failed flags.
*
* Parameters:
*
* iteration_failed: An output vector of flags which indicate whether the fitting
*                   process improved the fit in the last iteration. If yes it is set
*                   to 0 otherwise to 1.
*
* chi_squares: An input vector of chi-square values for multiple fits.
*
* prev_chi_squares: An input vector of chi-square values for multiple fits calculated
*                   in the previous iteration.
*
* n_fits: The number of fits.
*
* finished: An input vector which allows the calculation to be skipped
*           for single fits.
*
* Calling the cuda_check_fit_improvement function
* ===============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_check_fit_improvement <<< blocks, threads >>>(
*       iteration_failed,
*       chi_squares,
*       prev_chi_squares,
*       n_fits,
*       finished);
*
*/

__global__ void cuda_check_fit_improvement(
    int * iteration_failed,
    REAL const * chi_squares,
    REAL const * prev_chi_squares,
    int const n_fits,
    int const * finished)
{
    int const index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= n_fits || finished[index])
        return;

    bool const prev_chi_squares_initialized = prev_chi_squares[index] != 0.;
    // chi_squares[index] can be NaN which compares to false with any other number
    bool const chi_square_decreased = (chi_squares[index] < prev_chi_squares[index]);
    if (prev_chi_squares_initialized && !chi_square_decreased)
    {
        iteration_failed[index] = 1;
    }
    else
    {
        iteration_failed[index] = 0;
    }
}

/* Description of the cuda_calculate_chi_squares function
* ========================================================
*
* This function calls one of the estimator funktions depending on the input
* parameter estimator_id. The estimator function calculates the chi-square values.
* The calcluation is performed for multiple fits in parallel.
*
* Parameters:
*
* chi_squares: An output vector of concatenated chi-square values.
*
* states: An output vector of values which indicate whether the fitting process
*         was carreid out correctly or which problem occurred. In this function
*         it is only used for MLE. It is set to 3 if a fitting curve value is
*         negative. This vector includes the states for multiple fits.
*
* data: An input vector of data for multiple fits
*
* values: An input vector of concatenated sets of model function values.
*
* weights: An input vector of values for weighting chi-square, gradient and hessian,
*          while using LSE
*
* n_points: The number of data points per fit.
*
* n_fits: The number of fits.
*
* estimator_id: The estimator ID.
*
* finished: An input vector which allows the calculation to be skipped for single
*           fits.
*
* n_fits_per_block: The number of fits calculated by each thread block.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the cuda_calculate_chi_squares function
* ================================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = power_of_two_n_points * n_fits_per_block / n_blocks_per_fit;
*   blocks.x = n_fits / n_fits_per_block * n_blocks_per_fit;
*
*   int const shared_size = sizeof(REAL) * threads.x;
*
*   cuda_calculate_chi_squares<<< blocks, threads, shared_size >>>(
*       chi_squares,
*       states,
*       data,
*       values,
*       weights,
*       n_points,
*       n_fits,
*       estimator_id,
*       finished,
*       n_fits_per_block,
*       user_info,
*       user_info_size);
*
*/

__global__ void cuda_calculate_chi_squares(
    REAL * chi_squares,
    int * states,
    REAL const * data,
    REAL const * values,
    REAL const * weights,
    int const n_points,
    int const n_fits,
    int const estimator_id,
    int const * finished,
    int const n_fits_per_block,
    char * user_info,
    std::size_t const user_info_size)
{
    int const shared_size = blockDim.x / n_fits_per_block;
    int const fit_in_block = threadIdx.x / shared_size;
    int const fit_piece = blockIdx.x / n_fits;
    int const fit_index = blockIdx.x * n_fits_per_block + fit_in_block - fit_piece * n_fits;
    int const point_index = threadIdx.x - fit_in_block * shared_size + fit_piece * shared_size;
    int const first_point = fit_index * n_points;

    if (finished[fit_index])
    {
        return;
    }

    REAL const * current_data = &data[first_point];
    REAL const * current_weight = weights ? &weights[first_point] : NULL;
    REAL const * current_value = &values[first_point];
    int * current_state = &states[fit_index];

    extern __shared__ REAL extern_array[];

    volatile REAL * shared_chi_square
        = extern_array + (fit_in_block - fit_piece) * shared_size;

    if (point_index >= n_points)
    {
        shared_chi_square[point_index] = 0.;
    }

    if (point_index < n_points)
    {
        calculate_chi_square(
            estimator_id,
            shared_chi_square,
            point_index,
            current_data,
            current_value,
            current_weight,
            current_state,
            user_info,
            user_info_size);
    }
    shared_chi_square += fit_piece * shared_size;
    sum_up_floats(shared_chi_square, shared_size);
    chi_squares[fit_index + fit_piece * n_fits] = shared_chi_square[0];
}

/* Description of the cuda_sum_gradient_subtotals function
* ========================================================
*
* This function sums up the chi-square gradient subtotals in place.
*
* Parameters:
*
* gradients: A vector of gradient values for multiple fits.
*            in: subtotals
*            out: totals
*
* n_blocks_per_fit: The number of blocks used to calculate one fit
*
* n_fits: The number of fits.
*
* n_parameters_to_fit: The number of model parameters, that are not held fixed.
*
* skip: An input vector which allows the calculation to be skipped for single fits.
*
* finished: An input vector which allows the calculation to be skipped for single
*           fits.
*
* Calling the cuda_sum_gradient_subtotals function
* ================================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_sum_gradient_subtotals<<< blocks,threads >>>(
*       gradients,
*       n_blocks_per_fit,
*       n_fits,
*       n_parameters_to_fit,
*       skip,
*       finished);
*
*/

__global__ void cuda_sum_gradient_subtotals(
    REAL * gradients,
    REAL const * subtotals,
    int const n_blocks_per_fit,
    int const n_fits,
    int const n_parameters,
    int const * skip,
    int const * finished)
{
    int const index = blockIdx.x * blockDim.x + threadIdx.x;
    int const fit_index = index / n_parameters;

    if (fit_index >= n_fits || finished[fit_index] || skip[fit_index])
        return;

    REAL * gradient = gradients + index;
    REAL const * subtotal = subtotals + index;

    double sum = 0.0;
    for (int i = 0; i < n_blocks_per_fit; i++)
        sum += subtotal[i * n_fits * n_parameters];

    gradient[0] = sum;
}

/* Description of the cuda_calculate_gradients function
* =====================================================
*
* This function calls one of the gradient functions depending on the input
* parameter estimator_id. The gradient function calculates the gradient values
* of the chi-square function calling a __device__ function. The calcluation is
* performed for multiple fits in parallel.
*
* Parameters:
*
* gradients: An output vector of concatenated sets of gradient vector values.
*
* data: An input vector of data for multiple fits
*
* values: An input vector of concatenated sets of model function values.
*
* derivatives: An input vector of concatenated sets of model function partial
*              derivatives.
*
* weights: An input vector of values for weighting chi-square, gradient and hessian,
*          while using LSE
*
* n_points: The number of data points per fit.
*
* n_fits: The number of fits.
*
* n_parameters: The number of fitting curve parameters.
*
* n_parameters_to_fit: The number of fitting curve parameters, that are not held
*                      fixed.
*
* parameters_to_fit_indices: An input vector of indices of fitting curve parameters,
*                            that are not held fixed.
*
* estimator_id: The estimator ID.
*
* finished: An input vector which allows the calculation to be skipped for single
*           fits.
*
* skip: An input vector which allows the calculation to be skipped for single fits.
*
* n_fits_per_block: The number of fits calculated by each thread block.
*
* user_info: An input vector containing user information.
*
* user_info_size: The number of elements in user_info.
*
* Calling the cuda_calculate_gradients function
* =============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = power_of_two_n_points * n_fits_per_block / n_blocks_per_fit;
*   blocks.x = n_fits / n_fits_per_block * n_blocks_per_fit;
*
*   int const shared_size = sizeof(REAL) * threads.x;
*
*   cuda_calculate_gradients<<< blocks, threads, shared_size >>>(
*       gradients,
*       data,
*       values,
*       derivatives,
*       weight,
*       n_points,
*       n_fits,
*       n_parameters,
*       n_parameters_to_fit,
*       parameters_to_fit_indices,
*       estimator_id,
*       finished,
*       skip,
*       n_fits_per_block,
*       user_info,
*       user_info_size);
*
*/

__global__ void cuda_calculate_gradients(
    REAL * gradients,
    REAL const * data,
    REAL const * values,
    REAL const * derivatives,
    REAL const * weights,
    int const n_points,
    int const n_fits,
    int const n_parameters,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
    int const estimator_id,
    int const * finished,
    int const * skip,
    int const n_fits_per_block,
    char * user_info,
    std::size_t const user_info_size)
{
    int const shared_size = blockDim.x / n_fits_per_block;
    int const fit_in_block = threadIdx.x / shared_size;
    int const fit_piece = blockIdx.x / n_fits;
    int const fit_index = blockIdx.x * n_fits_per_block + fit_in_block - fit_piece * n_fits;
    int const point_index = threadIdx.x - fit_in_block * shared_size + fit_piece * shared_size;
    int const first_point = fit_index * n_points;

    if (finished[fit_index] || skip[fit_index])
    {
        return;
    }

    REAL const * current_data = &data[first_point];
    REAL const * current_weight = weights ? &weights[first_point] : NULL;
    REAL const * current_derivative = &derivatives[first_point * n_parameters];
    REAL const * current_value = &values[first_point];

    extern __shared__ REAL extern_array[];

    volatile REAL * shared_gradient = extern_array + (fit_in_block - fit_piece) * shared_size;

    if (point_index >= n_points)
    {
        shared_gradient[point_index] = 0.;
    }

    for (int parameter_index = 0; parameter_index < n_parameters_to_fit; parameter_index++)
    {
        if (point_index < n_points)
        {
            int const derivative_index = parameters_to_fit_indices[parameter_index] * n_points + point_index;

            calculate_gradient(
                estimator_id,
                shared_gradient,
                point_index,
                derivative_index,
                current_data,
                current_value,
                current_derivative,
                current_weight,
                user_info,
                user_info_size);
        }
        sum_up_floats(shared_gradient + fit_piece * shared_size, shared_size);
        gradients[(fit_index * n_parameters_to_fit + parameter_index) + fit_piece * n_fits * n_parameters_to_fit]
            = shared_gradient[fit_piece * shared_size];
    }
}

/* Description of the cuda_calculate_hessians function
* ====================================================
*
* This function calls one of the hessian function depending on the input
* parameter estimator_id. The hessian funcion calculates the hessian matrix
* values of the chi-square function calling a __device__ functions. The
* calcluation is performed for multiple fits in parallel.
*
* Parameters:
*
* hessians: An output vector of concatenated sets of hessian matrix values.
*
* data: An input vector of data for multiple fits
*
* values: An input vector of concatenated sets of model function values.
*
* derivatives: An input vector of concatenated sets of model function partial
*              derivatives.
*
* weights: An input vector of values for weighting chi-square, gradient and hessian,
*          while using LSE
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* n_parameters: The number of fitting curve parameters.
*
* n_parameters_to_fit: The number of fitting curve parameters, that are not held
*                      fixed.
*
* parameters_to_fit_indices: An input vector of indices of fitting curve parameters,
*                            that are not held fixed.
*
* estimator_id: The estimator ID.
*
* skip: An input vector which allows the calculation to be skipped for single fits.
*
* finished: An input vector which allows the calculation to be skipped for single
*           fits.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the cuda_calculate_hessians function
* ============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int n_unique_values = n_parameters_to_fit * (n_parameters_to_fit + 1) / 2;
*
*   threads.x
*       = min(n_unique_values * n_fits_per_block, max_threads_per_block);
*
*   blocks.y
*       = threads.x / max_threads_per_block
*       + int((threads.x % max_threads_per_block) > 0);
*
*   blocks.x
*       = n_fits / n_fits_per_block
*       + int((n_fits % n_fits_per_block) > 0);
*
*   cuda_calculate_hessians<<< blocks, threads >>>(
*       hessians,
*       data,
*       values,
*       derivatives,
*       weight,
*       n_fits,
*       n_points,
*       n_parameters,
*       n_parameters_to_fit,
*       parameters_to_fit_indices,
*       estimator_id,
*       skip,
*       finished,
*       user_info,
*       user_info_size);
*
*/

__global__ void cuda_calculate_hessians(
    REAL * hessians,
    REAL const * data,
    REAL const * values,
    REAL const * derivatives,
    REAL const * weights,
    int const n_fits,
    int const n_points,
    int const n_parameters,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
    int const estimator_id,
    int const * skip,
    int const * finished,
    char * user_info,
    std::size_t const user_info_size)
{
    int const n_unique_values = n_parameters_to_fit * (n_parameters_to_fit + 1) / 2;
    int const n_fits_per_block = blockDim.x * gridDim.y / n_unique_values;
    
    int const fit_in_block
        = (gridDim.y == 1)
        ? (blockIdx.y * blockDim.x + threadIdx.x) / n_unique_values
        : 0;

    int const fit_index = blockIdx.x * n_fits_per_block + fit_in_block;

    if (fit_index >= n_fits || finished[fit_index] || skip[fit_index])
    {
        return;
    }

    int const first_point = fit_index * n_points;
    int const parameter_index = (blockIdx.y * blockDim.x + threadIdx.x) - fit_in_block * n_unique_values;

    if (parameter_index >= n_unique_values)
    {
        return;
    }

    int const parameter_index_i
        = n_parameters_to_fit
        - 1.
        - std::floor(
            .5*(
                std::sqrt(
                    - 8. * (parameter_index - n_parameters_to_fit)
                    + 4. * n_parameters_to_fit * (n_parameters_to_fit - 1.)
                    - 7.
                ) - 1.
            )
        );

    int const parameter_index_j
        = parameter_index
        + parameter_index_i
        - parameter_index_i*(n_parameters_to_fit - (parameter_index_i - 1) / 2.);

    REAL * current_hessian = &hessians[fit_index * n_parameters_to_fit * n_parameters_to_fit];
    REAL const * current_data = &data[first_point];
    REAL const * current_weight = weights ? &weights[first_point] : NULL;
    REAL const * current_derivative = &derivatives[first_point*n_parameters];
    REAL const * current_value = &values[first_point];

    int const hessian_index_ij = parameter_index_i * n_parameters_to_fit + parameter_index_j;
    int const hessian_index_ji = parameter_index_j * n_parameters_to_fit + parameter_index_i;
    int const derivative_index_i = parameters_to_fit_indices[parameter_index_i] * n_points;
    int const derivative_index_j = parameters_to_fit_indices[parameter_index_j] * n_points;

    double sum = 0.0;
    for (int point_index = 0; point_index < n_points; point_index++)
    {
        calculate_hessian(
            estimator_id,
            &sum,
            point_index,
            derivative_index_i + point_index,
            derivative_index_j + point_index,
            current_data,
            current_value,
            current_derivative,
            current_weight,
            user_info,
            user_info_size);
    }
    current_hessian[hessian_index_ij] = sum;
    current_hessian[hessian_index_ji] = sum;
}

/* Description of the cuda_modify_step_widths function
* ====================================================
*
* This function midifies the diagonal elements of the hessian matrices by multiplying
* them by the factor (1+ lambda). This operation controls the step widths of the
* iteration. If the last iteration failed, befor modifying the hessian, the diagonal
* elements of the hessian are calculated back to represent unmodified values.
*
* hessians: An input and output vector of hessian matrices, which are modified by
*           the lambda values.
*
* lambdas: An input vector of values for modifying the hessians.
*
* n_parameters: The number of fitting curve parameters.
*
* iteration_failed: An input vector which indicates whether the previous iteration
*                   failed.
*
* finished: An input vector which allows the calculation to be skipped for single fits.
*
* n_fits_per_block: The number of fits calculated by each thread block.
*
* Calling the cuda_modify_step_widths function
* ============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = n_parameters_to_fit * n_fits_per_block;
*   blocks.x = n_fits / n_fits_per_block;
*
*   cuda_modify_step_width<<< blocks, threads >>>(
*       hessians,
*       lambdas,
*       n_parameters,
*       iteration_failed,
*       finished,
*       n_fits_per_block);
*
*/

__global__ void cuda_modify_step_widths(
    REAL * hessians,
    REAL const * lambdas,
    REAL * scaling_vectors,
    unsigned int const n_parameters,
    int const * iteration_failed,
    int const * finished,
    int const n_fits_per_block)
{
    int const shared_size = blockDim.x / n_fits_per_block;
    int const fit_in_block = threadIdx.x / shared_size;
    int const parameter_index = threadIdx.x - fit_in_block * shared_size;
    int const fit_index = blockIdx.x * n_fits_per_block + fit_in_block;

    if (finished[fit_index])
    {
        return;
    }

    REAL * hessian = &hessians[fit_index * n_parameters * n_parameters];
    REAL * scaling_vector = &scaling_vectors[fit_index * n_parameters];
    REAL const & lambda = lambdas[fit_index];

    int const diagonal_index = parameter_index * n_parameters + parameter_index;

    if (iteration_failed[fit_index])
    {
        hessian[diagonal_index] -= scaling_vector[parameter_index] * lambda / 10.;
    }

    // adaptive scaling
    scaling_vector[parameter_index]
        = max(scaling_vector[parameter_index], hessian[diagonal_index]);

    // continuous scaling
    //scaling_vector[parameter_index] = hessian[diagonal_index];
    
    // initial scaling
    //if (scaling_vector[parameter_index] == 0.)
    //    scaling_vector[parameter_index] = hessian[diagonal_index];

    hessian[diagonal_index] += scaling_vector[parameter_index] * lambda;
}

/* Description of the cuda_update_parameters function
* ===================================================
*
* This function stores the fitting curve parameter values in prev_parameters and
* updates them after each iteration.
*
* Parameters:
*
* parameters: An input and output vector of concatenated sets of model
*             parameters.
*
* prev_parameters: An input and output vector of concatenated sets of model
*                  parameters calculated by the previous iteration.
*
* deltas: An input vector of concatenated delta values, which are added to the
*         model parameters.
*
* n_parameters_to_fit: The number of fitted curve parameters.
*
* parameters_to_fit_indices: The indices of fitted curve parameters.
*
* finished: An input vector which allows the parameter update to be skipped for single fits.
*
* n_fits_per_block: The number of fits calculated by each threadblock.
*
* Calling the cuda_update_parameters function
* ===========================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = n_parameters * n_fits_per_block;
*   blocks.x = n_fits / n_fits_per_block;
*
*   cuda_update_parameters<<< blocks, threads >>>(
*       parameters,
*       prev_parameters,
*       deltas,
*       n_parameters_to_fit,
*       parameters_to_fit_indices,
*       finished,
*       n_fits_per_block);
*
*/

__global__ void cuda_update_parameters(
    REAL * parameters,
    REAL * prev_parameters,
    REAL const * deltas,
    int const n_parameters_to_fit,
    int const * parameters_to_fit_indices,
    int const * finished,
    int const n_fits_per_block)
{
    int const n_parameters = blockDim.x / n_fits_per_block;
    int const fit_in_block = threadIdx.x / n_parameters;
    int const parameter_index = threadIdx.x - fit_in_block * n_parameters;
    int const fit_index = blockIdx.x * n_fits_per_block + fit_in_block;

    REAL * current_parameters = &parameters[fit_index * n_parameters];
    REAL * current_prev_parameters = &prev_parameters[fit_index * n_parameters];

    current_prev_parameters[parameter_index] = current_parameters[parameter_index];

    if (finished[fit_index])
    {
        return;
    }

    if (parameter_index >= n_parameters_to_fit)
    {
        return;
    }

    REAL const * current_deltas = &deltas[fit_index * n_parameters_to_fit];

    current_parameters[parameters_to_fit_indices[parameter_index]] += current_deltas[parameter_index];
}

/* Description of the cuda_update_state_after_solving function
 * ===========================================================
 *
 * This function interprets the singular flag vector of the equation system
 * solving function according to this LM implementation.
 *
 * Parameters:
 *
 * n_fits: The number of fits.
 *
 * solution_info: An input vector used to report whether a fit is singular.
 *
 * finished: An input vector which allows the calculation to by skipped for
 *           single fits.
 *
 * gpufit_states: An output vector of values which indicate whether the fitting
 *                process was carreid out correctly or which problem occurred.
 *                If a hessian matrix of a fit is singular, it is set to 2.
 *
 * Calling the cuda_update_state_after_solving function
 * ====================================================
 *
 * When calling the function, the blocks and threads must be set up correctly,
 * as shown in the following example code.
 *
 *   dim3  threads(1, 1, 1);
 *   dim3  blocks(1, 1, 1);
 *
 *   int const example_value = 256;
 *
 *   threads.x = min(n_fits, example_value);
 *   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
 *
 *   cuda_update_state_after_solving<<< blocks, threads >>>(
 *       n_fits,
 *       solution_info,
 *       finished,
 *       gpufit_states);
 *
 */
    
__global__ void cuda_update_state_after_solving(
    int const n_fits,
    int const * cublas_info,
    int const * finished,
    int * states)
{
    int const fit_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (fit_index >= n_fits)
        return;

    if (finished[fit_index])
        return;

    if (cublas_info[fit_index] != 0)
        states[fit_index] = SINGULAR_HESSIAN;
}
    
/* Description of the cuda_check_for_convergence function
* =======================================================
*
* This function checks after each iteration whether the fits are converged or not.
* It also checks whether the set maximum number of iterations is reached.
*
* Parameters:
*
* finished: An input and output vector which allows the calculation to be skipped
*           for single fits.
*
* tolerance: The tolerance value for the convergence set by user.
*
* states: An output vector of values which indicate whether the fitting process
*         was carreid out correctly or which problem occurred. If the maximum
*         number of iterations is reached without converging, it is set to 1. If
*         the fit converged it keeps its initial value of 0.
*
* chi_squares: An input vector of chi-square values for multiple fits. Used for the
*              convergence check.
*
* prev_chi_squares: An input vector of chi-square values for multiple fits calculated
*                   in the previous iteration. Used for the convergence check.
*
* iteration: The value of the current iteration. It is compared to the value
*            of the maximum number of iteration set by user.
*
* max_n_iterations: The maximum number of iterations set by user.
*
* // NEW
* num_v_coefs: number of coefficients to be fitted. Mainly for linear_1d, ie. polynomial fit function
*
* n_fits: The number of fits.
*
* Calling the cuda_check_for_convergence function
* ===============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_check_for_convergence<<< blocks, threads >>>(
*       finished,
*       tolerance,
*       states,
*       chi_squares,
*       prev_chi_squares,
*       iteration,
*       max_n_iterations,
*       n_fits);
*
*/

__global__ void cuda_check_for_convergence(
    int * finished,
    REAL const tolerance,
    int * states,
    REAL const * chi_squares,
    REAL const * prev_chi_squares,
    int const iteration,
    int const max_n_iterations,
    int const n_fits)
{
    int const fit_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (fit_index >= n_fits)
    {
        return;
    }

    if (finished[fit_index])
    {
        return;
    }

    int const fit_found
        = abs(chi_squares[fit_index] - prev_chi_squares[fit_index])
        < tolerance * max(1., chi_squares[fit_index]);

    int const max_n_iterations_reached = iteration == max_n_iterations - 1;

    if (fit_found)
    {
        finished[fit_index] = 1;
    }
    else if (max_n_iterations_reached)
    {
        states[fit_index] = MAX_ITERATION;
    }
}

/* Description of the cuda_evaluate_iteration function
* ====================================================
*
* This function evaluates the current iteration.
*   - It marks a fit as finished if a problem occured.
*   - It saves the needed number of iterations if a fit finished.
*   - It checks if all fits finished
*
* Parameters:
*
* all_finished: An output flag, that indicates whether all fits finished.
*
* n_iterations: An output vector of needed iterations for each fit.
*
* finished: An input and output vector which allows the evaluation to be skipped
*           for single fits
*
* iteration: The values of the current iteration.
*
* states: An input vector of values which indicate whether the fitting process
*         was carreid out correctly or which problem occurred.
*
* n_fits: The number of fits.
*
* Calling the cuda_evaluate_iteration function
* ============================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_evaluate_iteration<<< blocks, threads >>>(
*       all_finished,
*       n_iterations,
*       finished,
*       iteration,
*       states,
*       n_fits);
*
*/

__global__ void cuda_evaluate_iteration(
    int * all_finished,
    int * n_iterations,
    int * finished,
    int const iteration,
    int const * states,
    int const n_fits)
{
    int const fit_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (fit_index >= n_fits)
    {
        return;
    }

    if (states[fit_index] != CONVERGED)
    {
        finished[fit_index] = 1;
    }

    if (finished[fit_index] && n_iterations[fit_index] == 0)
    {
        n_iterations[fit_index] = iteration + 1;
    }

    if (!finished[fit_index])
    {
        *all_finished = 0;
    }
}

/* Description of the cuda_prepare_next_iteration function
* ========================================================
*
* This function prepares the next iteration. It either updates previous
* chi-square values or sets currently calculated chi-square values and
* parameters to values calculated by the previous iteration. This function also
* updates lambda values.
*
* Parameters:
*
* lambdas: An output vector of values which control the step width by modifying
*          the diagonal elements of the hessian matrices.
*
* chi_squares: An input and output vector of chi-square values for multiple fits.
*
* prev_chi_squares: An input and output vector of chi-square values for multiple
*                   fits calculated in the previous iteration.
*
* parameters: An output vector of concatenated sets of model parameters.
*
* prev_parameters: An input vector of concatenated sets of model parameters
*                  calculated in the previous iteration.
*
* n_fits: The number of fits.
*
* n_parameters: The number of fitting curve parameters.
*
* Calling the cuda_prepare_next_iteration function
* ================================================
*
* When calling the function, the blocks and threads must be set up correctly,
* as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   int const example_value = 256;
*
*   threads.x = min(n_fits, example_value);
*   blocks.x = int(ceil(REAL(n_fits) / REAL(threads.x)));
*
*   cuda_prepare_next_iteration<<< blocks, threads >>>(
*       lambdas,
*       chi_squares,
*       prev_chi_squares,
*       parameters,
*       prev_parameters,
*       n_fits,
*       n_parameters);
*
*/

__global__ void cuda_prepare_next_iteration(
    REAL * lambdas,
    REAL * chi_squares,
    REAL * prev_chi_squares,
    REAL * parameters,
    REAL const * prev_parameters,
    int const n_fits,
    int const n_parameters)
{
    int const fit_index = blockIdx.x * blockDim.x + threadIdx.x;

    if (fit_index >= n_fits)
    {
        return;
    }

    if (chi_squares[fit_index] < prev_chi_squares[fit_index])
    {
        lambdas[fit_index] *= 0.1f;
        prev_chi_squares[fit_index] = chi_squares[fit_index];
    }
    else
    {
        lambdas[fit_index] *= 10.;
        chi_squares[fit_index] = prev_chi_squares[fit_index];
        for (int iparameter = 0; iparameter < n_parameters; iparameter++)
        {
            parameters[fit_index * n_parameters + iparameter] = prev_parameters[fit_index * n_parameters + iparameter];
        }
    }
}
/* -------------------------------------------------------------------------------------------------------
* from cuda_kernel.cu END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from info.h START
------------------------------------------------------------------------------------------------------- */
class Info
{
public:
    Info();
    virtual ~Info();

    void set_fits_per_block(std::size_t const n_fits);
    void set_number_of_parameters_to_fit(int const * parameters_to_fit);
    void configure();

private:
    void get_gpu_properties();
    void set_max_chunk_size();
    void set_blocks_per_fit();

public:
    int n_parameters_;
    int n_parameters_to_fit_;

    int n_points_;
    int power_of_two_n_points_;

    std::size_t n_fits_;

    std::size_t user_info_size_;

    int max_n_iterations_;
    int num_v_coefs_; // NEW
    std::size_t max_chunk_size_;

    int n_fits_per_block_;
    int n_blocks_per_fit_;
    ModelID model_id_;
    EstimatorID estimator_id_;

    bool use_weights_;

    int max_threads_;
    int warp_size_;

    DataLocation data_location_;

private:
    std::size_t max_blocks_;
    std::size_t available_gpu_memory_;
};

int getDeviceCount();
/* -------------------------------------------------------------------------------------------------------
* from info.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from info.cu START
------------------------------------------------------------------------------------------------------- */
void Info::get_gpu_properties()
{
    cudaDeviceProp devProp;
    CUDA_CHECK_STATUS(cudaGetDeviceProperties(&devProp, 0));
    max_threads_ = devProp.maxThreadsPerBlock;
    max_blocks_ = devProp.maxGridSize[0];
    warp_size_ = devProp.warpSize;

    std::size_t free_bytes;
    std::size_t total_bytes;
    CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    available_gpu_memory_ = std::size_t(double(free_bytes) * 0.1);
    
    if (double(user_info_size_) > double(free_bytes) * 0.9)
    {
        throw std::runtime_error("maximum user info size exceeded");
    }
}

int getDeviceCount()
{
	int deviceCount;
	CUDA_CHECK_STATUS(cudaGetDeviceCount(&deviceCount));
	return deviceCount;
}
/* -------------------------------------------------------------------------------------------------------
* from info.cu END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from info.cpp START
------------------------------------------------------------------------------------------------------- */
Info::Info() :
    n_parameters_(0),
    n_parameters_to_fit_(0),
    max_chunk_size_(0),
    max_n_iterations_(0),
    num_v_coefs_(0), // NEW
    n_points_(0),
    power_of_two_n_points_(0),
    n_fits_(0),
    user_info_size_(0),
    n_fits_per_block_(0),
    n_blocks_per_fit_(0),
    max_threads_(0),
    max_blocks_(0),
    warp_size_(0),
    available_gpu_memory_(0)
{
}

Info::~Info(void)
{
}

void Info::set_number_of_parameters_to_fit(int const * const parameters_to_fit)
{
    n_parameters_to_fit_ = n_parameters_;

    for (int i = 0; i < n_parameters_; i++)
    {
        if (!parameters_to_fit[i])
        {
            n_parameters_to_fit_--;
        }
    }
}

void Info::set_fits_per_block(std::size_t const current_chunk_size)
{

    n_fits_per_block_ = std::max((max_threads_ / power_of_two_n_points_), 1);

    bool is_divisible = current_chunk_size % n_fits_per_block_ == 0;

    while (!is_divisible && (n_fits_per_block_ > 1))
    {
        n_fits_per_block_ -= 1;
        is_divisible = current_chunk_size % n_fits_per_block_ == 0;
    }

}

void Info::set_blocks_per_fit()
{
    n_blocks_per_fit_ = 1;
    
    if (power_of_two_n_points_ > max_threads_)
    {
        bool enough_threads = false;
        do
        {
            n_blocks_per_fit_ *= 2;
            enough_threads = power_of_two_n_points_ / n_blocks_per_fit_ < max_threads_;
        } while (!enough_threads);
    }
}

void Info::set_max_chunk_size()
{
    int one_fit_memory
        = sizeof(REAL)
        *(1 * n_points_                                     // values
        + 1 * n_parameters_                                 // prev_parameters
        + 1 * n_parameters_to_fit_                          // gradient
        + 1 * n_parameters_to_fit_ * n_parameters_to_fit_   // hessian
        + 2 * n_parameters_to_fit_                          // delta, scaling_vector
        + 1 * n_points_*n_parameters_                       // derivatives
        + 2)                                                // prev_chi_square, lambda,
                                                            
        + sizeof(int)
        *(1 * n_parameters_to_fit_                          // indices of fitted parameters
        + 3);                                               // finished, iteration failed flag,
                                                            // solution info
    if (n_blocks_per_fit_ > 1)
    {
        one_fit_memory
            += sizeof(REAL)
             * n_parameters_to_fit_ * n_blocks_per_fit_;    // subtotals
    }

    if (data_location_ == HOST)
    {
        one_fit_memory += sizeof(REAL) * n_points_;        // data
        one_fit_memory += sizeof(REAL) * n_parameters_;    // parameters
        one_fit_memory += sizeof(REAL);                    // chi-square
        one_fit_memory += sizeof(int) * 2;                  // state, number of iterations
        if (use_weights_)
            one_fit_memory += sizeof(REAL) * n_points_;    // weights
    }

#ifdef USE_CUBLAS
    one_fit_memory
        += sizeof(REAL)
        *(2                                                 // pointer to decomposed hessian, pointer to delta
        + 1 * n_parameters_to_fit_ * n_parameters_to_fit_)  // decomposed hessian
        + sizeof(int)
        * (1 * n_parameters_to_fit_);                       // pivot vector
#endif // USE_CUBLAS
    
    std::size_t tmp_chunk_size = available_gpu_memory_ / one_fit_memory;
    
    if (tmp_chunk_size == 0)
    {
        throw std::runtime_error("not enough free GPU memory available");
    }

    tmp_chunk_size = (std::min)(tmp_chunk_size, max_blocks_ / n_blocks_per_fit_);

    std::size_t const highest_factor = n_points_ * n_parameters_;

    std::size_t const highest_size_t_value = std::numeric_limits< std::size_t >::max();

    if (tmp_chunk_size > highest_size_t_value / highest_factor)
    {
        tmp_chunk_size = highest_size_t_value / highest_factor;
    }

    max_chunk_size_ = tmp_chunk_size;

    int i = 1;
    int const divisor = 10;
    while (tmp_chunk_size > divisor)
    {
        i *= divisor;
        tmp_chunk_size /= divisor;
    }
    max_chunk_size_ = (max_chunk_size_ / i) * i;
    max_chunk_size_ = std::min(max_chunk_size_, n_fits_);
}


void Info::configure()
{
    power_of_two_n_points_ = 1;
    while (power_of_two_n_points_ < n_points_)
    {
        power_of_two_n_points_ *= 2;
    }
    
    // TODO NOTE: to address the 'too many resources requested' error due to complex calculate_acf1D function.
    // increasing the power_of_two_n_points_ by 4 times would reduce the number of threads required in Info::set_fits_per_block function
    // Info::set_blocks_per_fit functions function would also be updated accordingly.
    power_of_two_n_points_ *= 4;

    get_gpu_properties();
    set_blocks_per_fit();
    set_max_chunk_size();
}
/* -------------------------------------------------------------------------------------------------------
* from info.cpp END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from gpu_data.cuh START
------------------------------------------------------------------------------------------------------- */
template< typename Type >
struct Device_Array
{
    explicit Device_Array(std::size_t const size) : allocated_size_(size)
    {
        std::size_t const maximum_size = std::numeric_limits< std::size_t >::max();
        std::size_t const type_size = sizeof(Type);
        if (size <= maximum_size / type_size)
        {
            cudaError_t const status = cudaMalloc(&data_, size * type_size);
            if (status == cudaSuccess)
            {
                return;
            }
            else
            {
                throw std::runtime_error(cudaGetErrorString(status));
            }
        }
        else
        {
            throw std::runtime_error("maximum array size exceeded");
        }
    }

    ~Device_Array() { if (allocated_size_ > 0) cudaFree(data_); }

    operator Type * () { return static_cast<Type *>(data_); }
    operator Type const * () const { return static_cast<Type *>(data_); }

    Type const * data() const
    {
        return static_cast<Type *>(data_);
    }

    void assign(Type const * data)
    {
        data_ = const_cast<Type *>(data);
    }

    Type * copy(std::size_t const size, Type * const to) const
    {
        // TODO check size parameter
        std::size_t const type_size = sizeof(Type);
        cudaError_t const status
            = cudaMemcpy(to, data_, size * type_size, cudaMemcpyDeviceToHost);
        if (status == cudaSuccess)
        {
            return to + size;
        }
        else
        {
            throw std::runtime_error(cudaGetErrorString(status));
        }
    }

private:
    void * data_;
    std::size_t allocated_size_;
};

class GPUData
{
public:
    GPUData(Info const & info);
    ~GPUData();

    void init
    (
        int const chuk_size,
        int const chunk_index,
        REAL const * data,
        REAL const * weights,
        REAL const * initial_parameters,
        std::vector<int> const & parameters_to_fit_indices,
        int * states,
        REAL * chi_squares,
        int * n_iterations
    );
    void init_user_info(char const * user_info);

    void read(bool * dst, int const * src);
    void set(int* arr, int const value);
    void set(REAL* arr, REAL const value, int const count);
    void copy(REAL * dst, REAL const * src, std::size_t const count);

private:

    void set(int* arr, int const value, int const count);
    void write(REAL* dst, REAL const * src, int const count);
    void write(int* dst, std::vector<int> const & src);
    void write(char* dst, char const * src, std::size_t const count);
    void point_to_data_sets();

private:
    int chunk_size_;
    Info const & info_;

public:
    int chunk_index_;

    cublasHandle_t cublas_handle_;

    Device_Array< REAL > data_;
    Device_Array< REAL > weights_;
    Device_Array< REAL > parameters_;
    Device_Array< REAL > prev_parameters_;
    Device_Array< int > parameters_to_fit_indices_;
    Device_Array< char > user_info_;

    Device_Array< REAL > chi_squares_;
    Device_Array< REAL > prev_chi_squares_;
    Device_Array< REAL > gradients_;
    Device_Array< REAL > hessians_;
    Device_Array< REAL > deltas_;
    Device_Array< REAL > scaling_vectors_;
    Device_Array< REAL > subtotals_;

    Device_Array< REAL > values_;
    Device_Array< REAL > derivatives_;

    Device_Array< REAL > lambdas_;
    Device_Array< int > states_;
    Device_Array< int > finished_;
    Device_Array< int > iteration_failed_;
    Device_Array< int > all_finished_;
    Device_Array< int > n_iterations_;
    Device_Array< int > solution_info_;

#ifdef USE_CUBLAS
    Device_Array< REAL > decomposed_hessians_;
    Device_Array< REAL * > pointer_decomposed_hessians_;
    Device_Array< REAL * > pointer_deltas_;
    Device_Array< int > pivot_vectors_;
#endif // USE_CUBLAS
};
/* -------------------------------------------------------------------------------------------------------
* from gpu_data.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from gpu_data.cu START
------------------------------------------------------------------------------------------------------- */
GPUData::GPUData(Info const & info) :
    chunk_size_(0),
    info_(info),

    data_(
        (info_.data_location_ == HOST)
        ? info_.max_chunk_size_*info_.n_points_ : 0),
    weights_( 
        (info_.use_weights_ && info_.data_location_ == HOST)
        ? info_.n_points_ * info_.max_chunk_size_ : 0 ),
    parameters_(
        (info_.data_location_ == HOST)
        ? info_.max_chunk_size_*info_.n_parameters_ : 0 ),
    user_info_(
        (info_.data_location_ == HOST)
        ? info_.user_info_size_ : 0),

    prev_parameters_( info_.max_chunk_size_*info_.n_parameters_ ),
    parameters_to_fit_indices_( info_.n_parameters_to_fit_ ),

    chi_squares_(
        (info_.data_location_ == HOST)
        ? info_.max_chunk_size_ : 0),

    prev_chi_squares_( info_.max_chunk_size_ ),
    gradients_( info_.max_chunk_size_ * info_.n_parameters_to_fit_),
    hessians_( info_.max_chunk_size_ * info_.n_parameters_to_fit_ * info_.n_parameters_to_fit_ ),
    deltas_(info_.max_chunk_size_ * info_.n_parameters_to_fit_),
    scaling_vectors_(info_.max_chunk_size_ * info_.n_parameters_to_fit_),

    subtotals_(
        (info_.n_blocks_per_fit_ > 1)
        ? info_.max_chunk_size_ * info_.n_parameters_to_fit_ * info_.n_blocks_per_fit_ : 0),

    values_( info_.max_chunk_size_ * info_.n_points_ ),
    derivatives_( info_.max_chunk_size_ * info_.n_points_ * info_.n_parameters_ ),

    lambdas_( info_.max_chunk_size_ ),

    states_(
        (info_.data_location_ == HOST)
        ? info_.max_chunk_size_ : 0),
    
    finished_( info_.max_chunk_size_ ),
    iteration_failed_(info_.max_chunk_size_),
    all_finished_( 1 ),

    n_iterations_(
        (info_.data_location_ == HOST)
        ? info_.max_chunk_size_ : 0),
    
    solution_info_(info_.max_chunk_size_)

#ifdef USE_CUBLAS
    ,
    decomposed_hessians_(info_.max_chunk_size_ * info_.n_parameters_to_fit_ * info_.n_parameters_to_fit_),
    pointer_decomposed_hessians_(info_.max_chunk_size_),
    pointer_deltas_(info_.max_chunk_size_),
    pivot_vectors_(info_.max_chunk_size_ * info_.n_parameters_to_fit_)
#endif // USE_CUBLAS
{
#ifdef USE_CUBLAS
    cublasCreate(&cublas_handle_);
    point_to_data_sets();
#endif // USE_CUBLAS
}

GPUData::~GPUData()
{
#ifdef USE_CUBLAS
    cublasDestroy(cublas_handle_);
#endif // USE_CUBLAS
}

void GPUData::init
(
    int const chunk_size,
    int const chunk_index,
    REAL const * const data,
    REAL const * const weights,
    REAL const * const initial_parameters,
    std::vector<int> const & parameters_to_fit_indices,
    int * states,
    REAL * chi_squares,
    int * n_iterations)
{
    chunk_size_ = chunk_size;
    chunk_index_ = chunk_index;

    if (info_.data_location_ == HOST)
    {
        write(
            data_,
            data + chunk_index_*info_.max_chunk_size_*info_.n_points_,
            chunk_size_ * info_.n_points_);
        write(
            parameters_,
            initial_parameters + chunk_index_*info_.max_chunk_size_*info_.n_parameters_,
            chunk_size_ * info_.n_parameters_);
        if (info_.use_weights_)
            write(
                weights_,
                weights + chunk_index_*info_.max_chunk_size_*info_.n_points_,
                chunk_size_ * info_.n_points_);
    }
    else if (info_.data_location_ == DEVICE)
    {
        data_.assign(
            data + chunk_index_*info_.max_chunk_size_*info_.n_points_);
        parameters_.assign(
            initial_parameters + chunk_index_*info_.max_chunk_size_*info_.n_parameters_);
        if (info_.use_weights_)
            weights_.assign(
                weights + chunk_index_*info_.max_chunk_size_*info_.n_points_);
        states_.assign(
            states + chunk_index_ * info_.max_chunk_size_);
        chi_squares_.assign(
            chi_squares + chunk_index_ * info_.max_chunk_size_);
        n_iterations_.assign(
            n_iterations + chunk_index_ * info_.max_chunk_size_);
    }

    write(parameters_to_fit_indices_, parameters_to_fit_indices);

    set(prev_chi_squares_, 0., chunk_size_);
    set(finished_, 0, chunk_size_);
    set(scaling_vectors_, 0., chunk_size_ * info_.n_parameters_to_fit_);
    set(states_, 0, chunk_size_);
    set(lambdas_, 0.001f, chunk_size_);
}

void GPUData::init_user_info(char const * const user_info)
{
    if (info_.user_info_size_ > 0)
    {
        if (info_.data_location_ == HOST)
        {
            write(user_info_, user_info, info_.user_info_size_);
        }
        else if (info_.data_location_ == DEVICE)
        {
            user_info_.assign(user_info);
        }
    }
}

void GPUData::read(bool * dst, int const * src)
{
    int int_dst = 0;
    CUDA_CHECK_STATUS(cudaMemcpy(&int_dst, src, sizeof(int), cudaMemcpyDeviceToHost));
    * dst = (int_dst == 1) ? true : false;
}

void GPUData::write(REAL* dst, REAL const * src, int const count)
{
    CUDA_CHECK_STATUS(cudaMemcpy(dst, src, count * sizeof(REAL), cudaMemcpyHostToDevice));
}

void GPUData::write(int* dst, std::vector<int> const & src)
{
    std::size_t const size = src.size() * sizeof(int);
    CUDA_CHECK_STATUS(cudaMemcpy(dst, src.data(), size, cudaMemcpyHostToDevice));
}

void GPUData::write(char* dst, char const * src, std::size_t const count)
{
    CUDA_CHECK_STATUS(cudaMemcpy(dst, src, count * sizeof(char), cudaMemcpyHostToDevice));
}

void GPUData::copy(REAL * dst, REAL const * src, std::size_t const count)
{
    CUDA_CHECK_STATUS(cudaMemcpy(dst, src, count * sizeof(REAL), cudaMemcpyDeviceToDevice));
}

__global__ void set_kernel(int* dst, int const value, int const count)
{
    int const index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= count)
        return;

    dst[index] = value;
}

void GPUData::set(int* arr, int const value, int const count)
{
    int const tx = 256;
	int const bx = (count / tx) + 1;

    dim3  threads(tx, 1, 1);
    dim3  blocks(bx, 1, 1);

    set_kernel<<< blocks, threads >>>(arr, value, count);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void GPUData::set(int* arr, int const value)
{
    int const tx = 1;
    int const bx = 1;

    dim3  threads(tx, 1, 1);
    dim3  blocks(bx, 1, 1);

    set_kernel<<< blocks, threads >>>(arr, value, 1);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

__global__ void set_kernel(REAL* dst, REAL const value, std::size_t const count)
{
	std::size_t const index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= count)
        return;

    dst[index] = value;
}

void GPUData::set(REAL* arr, REAL const value, int const count)
{
    int const tx = 256;
	int const bx = (count / tx) + 1;

    dim3  threads(tx, 1, 1);
    dim3  blocks(bx, 1, 1);
    set_kernel<<< blocks, threads >>>(arr, value, count);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

__global__ void cuda_point_to_data_sets(
    REAL ** pointer_to_pointers,
    REAL * pointer,
    std::size_t const n_pointers,
    std::size_t const size)
{
    std::size_t const index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index >= n_pointers)
        return;

    int const begin = index * size;

    pointer_to_pointers[index] = pointer + begin;
}
#ifdef USE_CUBLAS

void GPUData::point_to_data_sets()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    std::size_t max_threads = 256;

    threads.x
        = static_cast<unsigned int>
          (std::min(info_.max_chunk_size_, max_threads));
    blocks.x
        = static_cast<unsigned int>
          (std::ceil(REAL(info_.max_chunk_size_) / REAL(threads.x)));

    cuda_point_to_data_sets <<< blocks, threads >>>(
        pointer_decomposed_hessians_,
        decomposed_hessians_,
        info_.max_chunk_size_,
        info_.n_parameters_to_fit_*info_.n_parameters_to_fit_);

    cuda_point_to_data_sets <<< blocks, threads >>> (
        pointer_deltas_,
        deltas_,
        info_.max_chunk_size_,
        info_.n_parameters_to_fit_);
}

#endif // USE_CUBLAS
/* -------------------------------------------------------------------------------------------------------
* from gpu_data.cu END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from cuda_gaussjordan.cuh START
------------------------------------------------------------------------------------------------------- */
/*
extern __global__ void cuda_gaussjordan(
    REAL * delta,
    REAL const * beta,
    REAL const * alpha,
    int const * skip_calculation,
    int * singular,
    std::size_t const n_equations,
    std::size_t const n_equations_pow2);
*/
/* -------------------------------------------------------------------------------------------------------
* from cuda_gaussjordan.cuh END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from cuda_gaussjordan.cu START
------------------------------------------------------------------------------------------------------- */
/* CUDA implementation of Gauss-Jordan elimination algorithm.
*  
* Gauss-Jordan elimination method
* ===============================
*
* This function solves a set of linear equations using the Gauss-Jordan elimination method.
* Considering a set of N equations with N unknowns, this can be written in matrix form as
* an NxN matrix of coefficients and a Nx1 column vector of right-hand side values.
*
* For example, consider the following problem with 3 equations and 3 unknowns (N=3):
* 
*   A x + B y + C z = MM
*   D x + E y + F z = NN
*   G x + H y + J z = PP
* 
* We can write this as follows in matrix form:
* 
*   [ A B C ] [ x ] = [ MM ]
*   [ D E F ] [ y ] = [ NN ] 
*   [ G H I ] [ z ] = [ PP ]
* 
* or, [A]*[X] = [B] where [A] is the matrix of coefficients and [B] is the vector of 
* right-hand side values.
*
* The Gauss Jordan elimiation method solves the system of equations in the following
* manner.  First, we form the augmented matrix (A|B):
*
*   [ A B C | MM ] 
*   [ D E F | NN ] 
*   [ G H I | PP ] 
*
* and then the augmented matrix is manipulated until its left side has the reduced
* row-echelon form.  That is to say that any individual row may be multiplied
* by a scalar factor, and any linear combination of rows may be added to another 
* row.  Finally, two rows may be swapped without affecting the solution.
* 
* When the manipulations are complete and the left side of the matrix has the desired
* form, the right side then corresponds to the solution of the system. 
*
*
* Description of the cuda_gaussjordan function
* ============================================
* 
* This algorithm is designed to perform many solutions of the Gauss Jordan elimination
* method in parallel.  One limitation of the algorithm implemented here is that for
* each solution the number of equations and unknowns (N) must be identical.  
*
* Parameters:
* 
* alpha: Coefficients matrices.  The matrix of coefficients for a single solution is 
*        a vector of NxN, where N is the number of equations.  This array stores the 
*        coefficients for the entire set of M input problems, concatenated end to end, 
*        and hence the total size of the array is MxNxN.  
*
* beta: Vector of right hand side values, concatenated together for all input problems. 
*       For a set of M inputs, the size of the vector is MxN.  Upon completion, this 
*       vector contains the results vector X for each solution.
*
* skip_calculation: An input vector which allows the calculation to be skipped for
*                   a particular solution.  For a set of M inputs, the size of this
*                   vector is M. 
*
* singular: An output vector used to report whether a given solution is singular.  For
*           a set of M inputs, this vector has size M.  Memory needs to be allocated
*           by the calling the function.
*
* n_equations: The number of equations and unknowns for a single solution.  This is
*              equal to the size N.
*
* n_equations_pow2: The next highest power of 2 greater than n_equations.
*
*
* Calling the cuda_gaussjordan function
* =====================================
*
* When calling the function, the blocks and threads must be set up correctly, as well
* as the shared memory space, as shown in the following example code.
*
*   dim3  threads(1, 1, 1);
*   dim3  blocks(1, 1, 1);
*
*   threads.x = n_equations + 1;
*   threads.y = n_equations;
*   blocks.x = n_solutions;
*   blocks.y = 1;
*
*   int const shared_size = sizeof(REAL) * 
*       ( (threads.x * threads.y) + n_parameters_pow2 + n_parameters_pow2 );
*
*   int * singular;
*   CUDA_CHECK_STATUS(cudaMalloc((void**)&singular, n_solutions * sizeof(int)));
*
*   cuda_gaussjordan<<< blocks, threads, shared_size >>>(
*       alpha,
*       beta,
*       skip_calculation,
*       singular,
*       n_equations,
*       n_equations_pow2);
*
*/

__global__ void cuda_gaussjordan(
    REAL * delta,
    REAL const * beta,
    REAL const * alpha,
    int const * skip_calculation,
    int * singular,
    std::size_t const n_equations,
    std::size_t const n_equations_pow2)
{
    extern __shared__ REAL extern_array[];     //shared memory between threads of a single block, 
    //used for storing the calculation_matrix, the 
    //abs_row vector, and the abs_row_index vector

    // In this routine we will store the augmented matrix (A|B), referred to here
    // as the calculation matrix in a shared memory space which is visible to all
    // threads within a block.  Also stored in shared memory are two vectors which 
    // are used to find the largest element in each row (the pivot).  These vectors 
    // are called abs_row and abs_row_index.
    //
    // Sizes of data stored in shared memory:
    //
    //      calculation_matrix: n_equations * (n_equations+1)
    //      abs_row:            n_equations_pow2
    //      abs_row_index:      n_equations_pow2
    //  
    // Note that each thread represents an element of the augmented matrix, with
    // the column and row indicated by the x and y index of the thread.  Each 
    // solution is calculated within one block, and the solution index is the 
    // block index x value.

    int const col_index = threadIdx.x;                  //column index in the calculation_matrix
    int const row_index = threadIdx.y;                  //row index in the calculation_matrix
    int const solution_index = blockIdx.x;

    int const n_col = blockDim.x;                       //number of columns in calculation matrix (=threads.x)
    int const n_row = blockDim.y;                       //number of rows in calculation matrix (=threads.y)
    int const alpha_size = blockDim.y * blockDim.y;     //number of entries in alpha matrix for one solution (NxN)

    if (skip_calculation[solution_index])
        return;

    REAL p;                                            //local variable used in pivot calculation

    REAL * calculation_matrix = extern_array;                          //point to the shared memory

    REAL * abs_row = extern_array + n_equations * (n_equations + 1);     //abs_row is located after the calculation_matrix
    //within the shared memory

    int * abs_row_index = (int *)(abs_row + n_equations_pow2);            //abs_row_index is located after abs_row
    //
    //note that although the shared memory is defined as
    //REAL, we are storing data of type int in this
    //part of the shared memory

    //initialize the singular vector
    if (col_index == 0 && row_index == 0)
    {
        singular[solution_index] = 0;
    }

    //initialize abs_row and abs_row_index, using only the threads on the diagonal
    if (col_index == row_index)
    {
        abs_row[col_index + (n_equations_pow2 - n_equations)] = 0.0;
        abs_row_index[col_index + (n_equations_pow2 - n_equations)] = col_index + (n_equations_pow2 - n_equations);
    }

    //initialize the calculation_matrix (alpha and beta, concatenated, for one solution)
    if (col_index != n_equations)
        calculation_matrix[row_index*n_col + col_index] = alpha[solution_index * alpha_size + row_index * n_equations + col_index];
    else
        calculation_matrix[row_index*n_col + col_index] = beta[solution_index * n_equations + row_index];

    //wait for thread synchronization

    __syncthreads();

    //start of main outer loop over the rows of the calculation matrix

    for (int current_row = 0; current_row < n_equations; current_row++)
    {

        // work in only one row, skipping the last column
        if (row_index == current_row && col_index != n_equations)
        {

            //save the absolute values of the current row
            abs_row[col_index] = abs(calculation_matrix[row_index * n_col + col_index]);

            //save the column indices
            abs_row_index[col_index] = col_index;

            __threadfence();

            //find the largest absolute value in the current row and write its index in abs_row_index[0]
            for (int n = 2; n <= n_equations_pow2; n = n * 2)
            {
                if (col_index < (n_equations_pow2 / n))
                {
                    if (abs_row[abs_row_index[col_index]] < abs_row[abs_row_index[col_index + (n_equations_pow2 / n)]])
                    {
                        abs_row_index[col_index] = abs_row_index[col_index + (n_equations_pow2 / n)];
                    }
                }
            }
        }

        __syncthreads();

        //singularity check - if all values in the row are zero, no solution exists
        if (row_index == current_row && col_index != n_equations)
        {
            if (abs_row[abs_row_index[0]] == 0.0)
            {
                singular[solution_index] = 1;
            }
        }

        //devide the row by the biggest value in the row
        if (row_index == current_row)
        {
            calculation_matrix[row_index * n_col + col_index]
                = calculation_matrix[row_index * n_col + col_index] / calculation_matrix[row_index * n_col + abs_row_index[0]];
        }

        __syncthreads();

        //The value of the largest element of the current row was found, and then current
        //row was divided by this value such that the largest value of the current row 
        //is equal to one.  
        //
        //Next, the matrix is manipulated to reduce to zero all other entries in the column 
        //in which the largest value was found.   To do this, the values in the current row
        //are scaled appropriately and substracted from the other rows of the matrix. 
        //
        //For each element of the matrix that is not in the current row, calculate the value
        //to be subtracted and let each thread store this value in the scalar variable p.

        p = calculation_matrix[current_row * n_col + col_index] * calculation_matrix[row_index * n_col + abs_row_index[0]];
        __syncthreads();

        if (row_index != current_row)
        {
            calculation_matrix[row_index * n_col + col_index] = calculation_matrix[row_index * n_col + col_index] - p;
        }
        __syncthreads();

    }

    //At this point, if the solution exists, the calculation matrix has been reduced to the 
    //identity matrix on the left side, and the solution vector on the right side.  However
    //we have not swapped rows during the procedure, so the identity matrix is out of order.
    //
    //For example, starting with the following augmented matrix as input:
    //
    //  [  3  2 -4 |  4 ]
    //  [  2  3  3 | 15 ]
    //  [  5 -3  1 | 14 ]
    //
    //we will obtain:
    //
    //  [  0  0  1 |  2 ]
    //  [  0  1  0 |  1 ]
    //  [  1  0  0 |  3 ]
    //
    //Which needs to be re-arranged to obtain the correct solution vector.  In the final
    //step, each thread checks to see if its value equals 1, and if so it assigns the value
    //in its rightmost column to the appropriate entry in the beta vector.  The solution is
    //stored in beta upon completetion.

    if (col_index != n_equations && calculation_matrix[row_index * n_col + col_index] == 1)
        delta[n_row * solution_index + col_index] = calculation_matrix[row_index * n_col + n_equations];

    __syncthreads();
}
/* -------------------------------------------------------------------------------------------------------
* from cuda_gaussjordan.cu END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from lm_fit.h START
------------------------------------------------------------------------------------------------------- */
/*
class LMFitCUDA;
*/

class LMFit
{
public:
    LMFit
    (
        REAL const * data,
        REAL const * weights,
        Info & info,
        REAL const * initial_parameters,
        int const * parameters_to_fit,
        char * user_info,
        REAL * output_parameters,
        int * output_states,
        REAL * output_chi_squares,
        int * output_n_iterations
    ) ;

    virtual ~LMFit();

    void run(REAL const tolerance);

private:
    void set_parameters_to_fit_indices();
    void get_results(GPUData const & gpu_data, int const n_fits);

    REAL const * const data_ ;
    REAL const * const weights_ ;
    REAL const * const initial_parameters_ ;
    int const * const parameters_to_fit_;
    char const * const user_info_;

    REAL * output_parameters_ ;
    int * output_states_ ;
    REAL * output_chi_squares_ ;
    int * output_n_iterations_ ;

    int ichunk_;
    int chunk_size_;
    std::size_t n_fits_left_;

    Info & info_;

    std::vector<int> parameters_to_fit_indices_;
};

class LMFitCUDA
{
public:
    LMFitCUDA(
        REAL const tolerance,
        Info const & info,
        GPUData & gpu_data,
        int const n_fits);

    virtual ~LMFitCUDA();

    void run();

private:
    void calc_curve_values();
    void calc_chi_squares();
    void calc_gradients();
    void calc_hessians();
    void evaluate_iteration(int const iteration);
    void scale_hessians();
#ifdef USE_CUBLAS
    void solve_equation_systems_lup();
#else
    void solve_equation_systems_gj();
#endif
    void update_states();
    void update_parameters();

public:

private:
    Info const & info_;
    GPUData & gpu_data_;
    int const n_fits_;

    bool all_finished_;

    REAL tolerance_;
};
/* -------------------------------------------------------------------------------------------------------
* from lm_fit.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from lm_fit.cpp START
------------------------------------------------------------------------------------------------------- */
LMFit::LMFit
(
    REAL const * const data,
    REAL const * const weights,
    Info & info,
    REAL const * const initial_parameters,
    int const * const parameters_to_fit,
    char * const user_info,
    REAL * output_parameters,
    int * output_states,
    REAL * output_chi_squares,
    int * output_n_iterations
) :
    data_( data ),
    weights_( weights ),
    initial_parameters_( initial_parameters ),
    parameters_to_fit_( parameters_to_fit ),
    user_info_( user_info ),
    output_parameters_( output_parameters ),
    output_states_( output_states ),
    output_chi_squares_( output_chi_squares ),
    output_n_iterations_( output_n_iterations ),
    info_(info),
    chunk_size_(0),
    ichunk_(0),
    n_fits_left_(info.n_fits_),
    parameters_to_fit_indices_(0)
{}

LMFit::~LMFit()
{}

void LMFit::set_parameters_to_fit_indices()
{
    int const n_parameters_to_fit = info_.n_parameters_;
    for (int i = 0; i < n_parameters_to_fit; i++)
    {
        if (parameters_to_fit_[i])
        {
            parameters_to_fit_indices_.push_back(i);
        }
    }
}

void LMFit::get_results(GPUData const & gpu_data, int const n_fits)
{
    if (info_.data_location_ == HOST)
    {
        output_parameters_
            = gpu_data.parameters_.copy(n_fits*info_.n_parameters_, output_parameters_);
        output_states_
            = gpu_data.states_.copy(n_fits, output_states_);
        output_chi_squares_
            = gpu_data.chi_squares_.copy(n_fits, output_chi_squares_);
        output_n_iterations_
            = gpu_data.n_iterations_.copy(n_fits, output_n_iterations_);
    }
}

void LMFit::run(REAL const tolerance)
{
    set_parameters_to_fit_indices();

    GPUData gpu_data(info_);
    gpu_data.init_user_info(user_info_);

    // loop over data chunks
    while (n_fits_left_ > 0)
    {
        chunk_size_ = int((std::min)(n_fits_left_, info_.max_chunk_size_));

        info_.set_fits_per_block(chunk_size_);

        gpu_data.init(
            chunk_size_,
            ichunk_,
            data_,
            weights_,
            initial_parameters_,
            parameters_to_fit_indices_,
            output_states_,
            output_chi_squares_,
            output_n_iterations_);

        LMFitCUDA lmfit_cuda(
            tolerance,
            info_,
            gpu_data,
            chunk_size_);

        lmfit_cuda.run();

        get_results(gpu_data, chunk_size_);

        n_fits_left_ -= chunk_size_;
        ichunk_++;

    }
}
/* -------------------------------------------------------------------------------------------------------
* from lm_fit.cpp END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from lm_fit_cuda.cu START
------------------------------------------------------------------------------------------------------- */
#ifdef USE_CUBLAS
    void LMFitCUDA::solve_equation_systems_lup()
    {
        dim3  threads(1, 1, 1);
        dim3  blocks(1, 1, 1);

        // initialize components of equation systems
        gpu_data_.copy(gpu_data_.decomposed_hessians_, gpu_data_.hessians_, n_fits_ * info_.n_parameters_to_fit_ * info_.n_parameters_to_fit_);

        // decompose hessians
        cublasStatus_t lu_status_decopmposition = DECOMPOSE_LUP(
            gpu_data_.cublas_handle_,
            info_.n_parameters_to_fit_,
            gpu_data_.pointer_decomposed_hessians_,
            info_.n_parameters_to_fit_,
            gpu_data_.pivot_vectors_,
            gpu_data_.solution_info_,
            n_fits_);

        // initialize deltas with values of gradients
        gpu_data_.copy(gpu_data_.deltas_, gpu_data_.gradients_, n_fits_ * info_.n_parameters_to_fit_);

        // TODO: check solution_info
        int solution_info;
    
        // solve equation systems
        cublasStatus_t lu_status_solution
            = SOLVE_LUP(
            gpu_data_.cublas_handle_,
            CUBLAS_OP_N,
            info_.n_parameters_to_fit_,
            1,
            (REAL const **)(gpu_data_.pointer_decomposed_hessians_.data()),
            info_.n_parameters_to_fit_,
            gpu_data_.pivot_vectors_,
            gpu_data_.pointer_deltas_,
            info_.n_parameters_to_fit_,
            &solution_info,
            n_fits_);
    }
#else
    void LMFitCUDA::solve_equation_systems_gj()
    {
        dim3  threads(1, 1, 1);
        dim3  blocks(1, 1, 1);

        int n_parameters_pow2 = 1;

        while (n_parameters_pow2 < info_.n_parameters_to_fit_)
        {
            n_parameters_pow2 *= 2;
        }

        //set up to run the Gauss Jordan elimination
        int const n_equations = info_.n_parameters_to_fit_;
        int const n_solutions = n_fits_;

        threads.x = n_equations + 1;
        threads.y = n_equations;
        blocks.x = n_solutions;

        //set the size of the shared memory area for each block
        int const shared_size = sizeof(REAL) * ((threads.x * threads.y) + n_parameters_pow2 + n_parameters_pow2);

        //run the Gauss Jordan elimination
        cuda_gaussjordan <<< blocks, threads, shared_size >>>(
            gpu_data_.deltas_,
            gpu_data_.gradients_,
            gpu_data_.hessians_,
            gpu_data_.finished_,
            gpu_data_.solution_info_,
            info_.n_parameters_to_fit_,
            n_parameters_pow2);
        CUDA_CHECK_STATUS(cudaGetLastError());
    }
#endif

void LMFitCUDA::update_states()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    //set up to update the lm_state_gpu_ variable with the Gauss Jordan results
    threads.x = std::min(n_fits_, 256);
    blocks.x = int(std::ceil(REAL(n_fits_) / REAL(threads.x)));

    //update the gpu_data_.states_ variable
    cuda_update_state_after_solving <<< blocks, threads >>>(
        n_fits_,
        gpu_data_.solution_info_,
        gpu_data_.finished_,
        gpu_data_.states_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::scale_hessians()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = info_.n_parameters_to_fit_*info_.n_fits_per_block_;
    blocks.x = n_fits_ / info_.n_fits_per_block_;

    cuda_modify_step_widths <<< blocks, threads >>>(
        gpu_data_.hessians_,
        gpu_data_.lambdas_,
        gpu_data_.scaling_vectors_,
        info_.n_parameters_to_fit_,
        gpu_data_.iteration_failed_,
        gpu_data_.finished_,
        info_.n_fits_per_block_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::update_parameters()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = info_.n_parameters_*info_.n_fits_per_block_;
    blocks.x = n_fits_ / info_.n_fits_per_block_;

    cuda_update_parameters <<< blocks, threads >>>(
        gpu_data_.parameters_,
        gpu_data_.prev_parameters_,
        gpu_data_.deltas_,
        info_.n_parameters_to_fit_,
        gpu_data_.parameters_to_fit_indices_,
        gpu_data_.finished_,
        info_.n_fits_per_block_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::calc_curve_values()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = info_.n_points_ * info_.n_fits_per_block_ / info_.n_blocks_per_fit_;

    if (info_.n_blocks_per_fit_ > 1)
        threads.x += info_.n_points_ % threads.x;

    threads.x = threads.x;
    blocks.x = n_fits_ / info_.n_fits_per_block_ * info_.n_blocks_per_fit_;

    cuda_calc_curve_values <<< blocks, threads >>>(
        gpu_data_.parameters_,
        n_fits_,
        info_.n_points_,
        info_.n_parameters_,
        gpu_data_.finished_,
        gpu_data_.values_,
        gpu_data_.derivatives_,
        info_.n_fits_per_block_,
        info_.n_blocks_per_fit_,
        info_.model_id_,
        gpu_data_.chunk_index_,
        gpu_data_.user_info_,
        info_.user_info_size_,
        info_.num_v_coefs_); // NEW
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::calc_chi_squares()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = info_.power_of_two_n_points_ * info_.n_fits_per_block_ / info_.n_blocks_per_fit_;
    blocks.x = n_fits_ / info_.n_fits_per_block_ * info_.n_blocks_per_fit_;

    int const shared_size = sizeof(REAL) * threads.x;

    REAL * chi_squares = 
        info_.n_blocks_per_fit_ > 1 ? gpu_data_.subtotals_ : gpu_data_.chi_squares_;

    cuda_calculate_chi_squares <<< blocks, threads, shared_size >>>(
        chi_squares,
        gpu_data_.states_,
        gpu_data_.data_,
        gpu_data_.values_,
        gpu_data_.weights_,
        info_.n_points_,
        n_fits_,
        info_.estimator_id_,
        gpu_data_.finished_,
        info_.n_fits_per_block_,
        gpu_data_.user_info_,
        info_.user_info_size_);
    CUDA_CHECK_STATUS(cudaGetLastError());

    threads.x = std::min(n_fits_, 256);
    blocks.x = int(std::ceil(REAL(n_fits_) / REAL(threads.x)));

    if (info_.n_blocks_per_fit_ > 1)
    {
        cuda_sum_chi_square_subtotals <<< blocks, threads >>> (
            gpu_data_.chi_squares_,
            gpu_data_.subtotals_,
            info_.n_blocks_per_fit_,
            n_fits_,
            gpu_data_.finished_);
        CUDA_CHECK_STATUS(cudaGetLastError());
    }

    cuda_check_fit_improvement <<< blocks, threads >>>(
        gpu_data_.iteration_failed_,
        gpu_data_.chi_squares_,
        gpu_data_.prev_chi_squares_,
        n_fits_,
        gpu_data_.finished_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::calc_gradients()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = info_.power_of_two_n_points_ * info_.n_fits_per_block_ / info_.n_blocks_per_fit_;
    blocks.x = n_fits_ / info_.n_fits_per_block_ * info_.n_blocks_per_fit_;

    int const shared_size = sizeof(REAL) * threads.x;

    REAL * gradients
        = info_.n_blocks_per_fit_ > 1 ? gpu_data_.subtotals_ : gpu_data_.gradients_;

    cuda_calculate_gradients <<< blocks, threads, shared_size >>>(
        gradients,
        gpu_data_.data_,
        gpu_data_.values_,
        gpu_data_.derivatives_,
        gpu_data_.weights_,
        info_.n_points_,
        n_fits_,
        info_.n_parameters_,
        info_.n_parameters_to_fit_,
        gpu_data_.parameters_to_fit_indices_,
        info_.estimator_id_,
        gpu_data_.finished_,
        gpu_data_.iteration_failed_,
        info_.n_fits_per_block_,
        gpu_data_.user_info_,
        info_.user_info_size_);
    CUDA_CHECK_STATUS(cudaGetLastError());

    if (info_.n_blocks_per_fit_ > 1)
    {
        int const gradients_size = n_fits_ * info_.n_parameters_to_fit_;
        threads.x = std::min(gradients_size, 256);
        blocks.x = int(std::ceil(REAL(gradients_size) / REAL(threads.x)));

        cuda_sum_gradient_subtotals <<< blocks, threads >>> (
            gpu_data_.gradients_,
            gpu_data_.subtotals_,
            info_.n_blocks_per_fit_,
            n_fits_,
            info_.n_parameters_to_fit_,
            gpu_data_.iteration_failed_,
            gpu_data_.finished_);
        CUDA_CHECK_STATUS(cudaGetLastError());
    }
}

void LMFitCUDA::calc_hessians()
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    int const n_unique_values
        = info_.n_parameters_to_fit_ * (info_.n_parameters_to_fit_ + 1) / 2;

    int n_hessians_per_block = 1;

    if (info_.n_parameters_to_fit_)
    {
        while ((n_hessians_per_block + 1) * n_unique_values < info_.warp_size_)
        {
            n_hessians_per_block++;
        }
    }

    int const temp_threads_x = n_unique_values * n_hessians_per_block;

    threads.x = std::min(temp_threads_x, info_.max_threads_);
    
    blocks.y
        = temp_threads_x / info_.max_threads_ 
        + int((temp_threads_x % info_.max_threads_) > 0);
    
    blocks.x
        = n_fits_ / n_hessians_per_block
        + int((n_fits_ % n_hessians_per_block) > 0);

    cuda_calculate_hessians <<< blocks, threads >>>(
        gpu_data_.hessians_,
        gpu_data_.data_,
        gpu_data_.values_,
        gpu_data_.derivatives_,
        gpu_data_.weights_,
        n_fits_,
        info_.n_points_,
        info_.n_parameters_,
        info_.n_parameters_to_fit_,
        gpu_data_.parameters_to_fit_indices_,
        info_.estimator_id_,
        gpu_data_.iteration_failed_,
        gpu_data_.finished_,
        gpu_data_.user_info_,
        info_.user_info_size_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}

void LMFitCUDA::evaluate_iteration(int const iteration)
{
    dim3  threads(1, 1, 1);
    dim3  blocks(1, 1, 1);

    threads.x = std::min(n_fits_, 256);
    blocks.x = int(std::ceil(REAL(n_fits_) / REAL(threads.x)));

    cuda_check_for_convergence<<< blocks, threads >>>(
        gpu_data_.finished_,
        tolerance_,
        gpu_data_.states_,
        gpu_data_.chi_squares_,
        gpu_data_.prev_chi_squares_,
        iteration,
        info_.max_n_iterations_,
        n_fits_);
    CUDA_CHECK_STATUS(cudaGetLastError());

    gpu_data_.set(gpu_data_.all_finished_, 1);

    cuda_evaluate_iteration<<< blocks, threads >>>(
        gpu_data_.all_finished_,
        gpu_data_.n_iterations_,
        gpu_data_.finished_,
        iteration,
        gpu_data_.states_,
        n_fits_);
    CUDA_CHECK_STATUS(cudaGetLastError());

    gpu_data_.read(&all_finished_, gpu_data_.all_finished_);

    cuda_prepare_next_iteration<<< blocks, threads >>>(
        gpu_data_.lambdas_,
        gpu_data_.chi_squares_,
        gpu_data_.prev_chi_squares_,
        gpu_data_.parameters_,
        gpu_data_.prev_parameters_,
        n_fits_,
        info_.n_parameters_);
    CUDA_CHECK_STATUS(cudaGetLastError());
}
/* -------------------------------------------------------------------------------------------------------
* from lm_fit_cuda.cu END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from lm_fit_cuda.cpp START
------------------------------------------------------------------------------------------------------- */
LMFitCUDA::LMFitCUDA(
    REAL const tolerance,
    Info const & info,
    GPUData & gpu_data,
    int const n_fits
    ) :
    info_(info),
    gpu_data_(gpu_data),
    n_fits_(n_fits),
    all_finished_(false),
    tolerance_(tolerance)
{
}

LMFitCUDA::~LMFitCUDA()
{
}

void LMFitCUDA::run()
{
    // initialize the chi-square values
    calc_curve_values();
    calc_chi_squares();

    if (info_.n_parameters_to_fit_ == 0)
        return;

    calc_gradients();
    calc_hessians();

    gpu_data_.copy(
        gpu_data_.prev_chi_squares_,
        gpu_data_.chi_squares_,
        n_fits_);

    // loop over the fit iterations
    for (int iteration = 0; !all_finished_; iteration++)
    {
        // modify step width
        // LUP decomposition
        // update fitting parameters
        scale_hessians();
        SOLVE_EQUATION_SYSTEMS();
        update_states();
        update_parameters();

        // calculate fitting curve values and its derivatives
        // calculate chi-squares, gradients and hessians
	calc_curve_values();
        calc_chi_squares();
        calc_gradients();
        calc_hessians();

        // check which fits have converged
        // flag finished fits
        // check whether all fits finished
        // save the number of needed iterations by each fitting process
        // check whether chi-squares are increasing or decreasing
        // update chi-squares, curve parameters and lambdas
        evaluate_iteration(iteration);
    }
}
/* -------------------------------------------------------------------------------------------------------
* from lm_fit_cuda.cpp END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from interface.h START
------------------------------------------------------------------------------------------------------- */
static_assert( sizeof( int ) == 4, "32 bit 'int' type required" ) ;

class FitInterface
{
public:
    FitInterface
    (
        REAL const * data,
        REAL const * weights,
        std::size_t n_fits,
        int n_points,
        REAL tolerance,
        int max_n_iterations,
        int num_v_coefs,
        EstimatorID estimator_id,
        REAL const * initial_parameters,
        int * parameters_to_fit,
        char * user_info,
        std::size_t user_info_size,
        REAL * output_parameters,
        int * output_states,
        REAL * output_chi_squares,
        int * output_n_iterations,
        DataLocation data_location
    ) ;
    
    virtual ~FitInterface();
    void fit(ModelID const model_id);

private:
    void check_sizes();
    void configure_info(Info & info, ModelID const model_id);

public:

private:
    //input
    REAL const * const data_ ;
    REAL const * const weights_;
    REAL const * const initial_parameters_;
    int const * const parameters_to_fit_;
    char * const user_info_;
    int n_parameters_;

    std::size_t const n_fits_;
    int const n_points_;
    REAL const  tolerance_;
    int const max_n_iterations_;
    int const num_v_coefs_; // NEW
    EstimatorID estimator_id_;
    std::size_t const user_info_size_;

    DataLocation data_location_;

    //output
    REAL * output_parameters_;
    int * output_states_;
    REAL * output_chi_squares_;
    int * output_n_iterations_;
};
/* -------------------------------------------------------------------------------------------------------
* from interface.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from interface.cpp START
------------------------------------------------------------------------------------------------------- */
FitInterface::FitInterface
(
    REAL const * data,
    REAL const * weights,
    std::size_t n_fits,
    int n_points,
    REAL tolerance,
    int max_n_iterations,
    int num_v_coefs, // NEW
    EstimatorID estimator_id,
    REAL const * initial_parameters,
    int * parameters_to_fit,
    char * user_info,
    std::size_t user_info_size,
    REAL * output_parameters,
    int * output_states,
    REAL * output_chi_squares,
    int * output_n_iterations,
    DataLocation data_location
) :
    data_( data ),
    weights_( weights ),
    initial_parameters_( initial_parameters ),
    parameters_to_fit_( parameters_to_fit ),
    user_info_( user_info ),
    n_fits_(n_fits),
    n_points_(n_points),
    tolerance_(tolerance),
    max_n_iterations_(max_n_iterations),
    num_v_coefs_(num_v_coefs), // NEW
    estimator_id_(estimator_id),
    user_info_size_(user_info_size),
    output_parameters_( output_parameters ),
    output_states_(output_states),
    output_chi_squares_(output_chi_squares),
    output_n_iterations_(output_n_iterations),
    n_parameters_(0),
    data_location_(data_location)
{}

FitInterface::~FitInterface()
{}

void FitInterface::check_sizes()
{
    std::size_t maximum_size = std::numeric_limits< std::size_t >::max();
    
    if (n_fits_ > maximum_size / n_points_ / sizeof(REAL))
    {
        throw std::runtime_error("maximum absolute number of data points exceeded");
    }
    
    if (n_fits_ > maximum_size / n_parameters_ / sizeof(REAL))
    {
        throw std::runtime_error("maximum number of fits and/or parameters exceeded");
    }
}

void FitInterface::configure_info(Info & info, ModelID const model_id)
{
    info.model_id_ = model_id;
    info.n_fits_ = n_fits_;
    info.n_points_ = n_points_;
    info.max_n_iterations_ = max_n_iterations_;
    info.num_v_coefs_ = num_v_coefs_; // NEW
    info.estimator_id_ = estimator_id_;
    info.user_info_size_ = user_info_size_;
    info.n_parameters_ = n_parameters_;
    info.use_weights_ = weights_ ? true : false;
    info.data_location_ = data_location_;

    info.set_number_of_parameters_to_fit(parameters_to_fit_);
    info.configure();
}

void FitInterface::fit(ModelID const model_id)
{
    int n_dimensions = 0;
    configure_model(model_id, n_parameters_, n_dimensions);

    check_sizes();

    Info info;
    configure_info(info, model_id);

    LMFit lmfit
    (
        data_,
        weights_,
        info,
        initial_parameters_,
        parameters_to_fit_,
        user_info_,
        output_parameters_,
        output_states_,
        output_chi_squares_,
        output_n_iterations_
    ) ;
    lmfit.run(tolerance_);
}
/* -------------------------------------------------------------------------------------------------------
* from interface.cpp END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from gpufit.h START
------------------------------------------------------------------------------------------------------- */
/*
#ifdef __cplusplus
extern "C" {
#endif

VISIBLE int gpufit
(
    size_t n_fits,
    size_t n_points,
    REAL * data,
    REAL * weights,
    int model_id,
    REAL * initial_parameters,
    REAL tolerance,
    int max_n_iterations,
    int num_v_coefs, // NEW
    int * parameters_to_fit,
    int estimator_id,
    size_t user_info_size,
    char * user_info,
    REAL * output_parameters,
    int * output_states,
    REAL * output_chi_squares,
    int * output_n_iterations
) ;

VISIBLE int gpufit_cuda_interface
(
    size_t n_fits,
    size_t n_points,
    REAL * gpu_data,
    REAL * gpu_weights,
    int model_id,
    REAL tolerance,
    int max_n_iterations,
    int num_v_coefs, // NEW
    int * parameters_to_fit,
    int estimator_id,
    size_t user_info_size,
    char * gpu_user_info,
    REAL * gpu_fit_parameters,
    int * gpu_output_states,
    REAL * gpu_output_chi_squares,
    int * gpu_output_n_iterations
);

VISIBLE char const * gpufit_get_last_error() ;

//// returns 1 if cuda is available and 0 otherwise
//VISIBLE int gpufit_cuda_available();

VISIBLE int gpufit_get_cuda_version(int * runtime_version, int * driver_version);

VISIBLE int gpufit_portable_interface(int argc, void *argv[]);

#ifdef __cplusplus
}
#endif
*/
/* -------------------------------------------------------------------------------------------------------
* from gpufit.h END
------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------
* from gpufit.cpp START
------------------------------------------------------------------------------------------------------- */
std::string last_error ;

int gpufit
(
    size_t n_fits,
    size_t n_points,
    REAL * data,
    REAL * weights,
    int model_id,
    REAL * initial_parameters,
    REAL tolerance,
    int max_n_iterations,
    int num_v_coefs, // NEW
    int * parameters_to_fit,
    int estimator_id,
    size_t user_info_size,
    char * user_info,
    REAL * output_parameters,
    int * output_states,
    REAL * output_chi_squares,
    int * output_n_iterations
)
try
{
    FitInterface fi(
        data,
        weights,
        n_fits,
        static_cast<int>(n_points),
        tolerance,
        max_n_iterations,
        num_v_coefs, // NEW
        static_cast<EstimatorID>(estimator_id),
        initial_parameters,
        parameters_to_fit,
        user_info,
        user_info_size,
        output_parameters,
        output_states,
        output_chi_squares,
        output_n_iterations,
        HOST);

    fi.fit(static_cast<ModelID>(model_id));

    return ReturnState::OK ;
}
catch( std::exception & exception )
{
    last_error = exception.what() ;

    return ReturnState::ERROR ;
}
catch( ... )
{
    last_error = "unknown error" ;

    return ReturnState::ERROR;
}

int gpufit_cuda_interface
(
    size_t n_fits,
    size_t n_points,
    REAL * gpu_data,
    REAL * gpu_weights,
    int model_id,
    REAL tolerance,
    int max_n_iterations,
    int num_v_coefs, // NEW
    int * parameters_to_fit,
    int estimator_id,
    size_t user_info_size,
    char * gpu_user_info,
    REAL * gpu_fit_parameters,
    int * gpu_output_states,
    REAL * gpu_output_chi_squares,
    int * gpu_output_n_iterations
)
try
{
    FitInterface fi(
        gpu_data,
        gpu_weights,
        n_fits,
        static_cast<int>(n_points),
        tolerance,
        max_n_iterations,
        num_v_coefs, // NEW
        static_cast<EstimatorID>(estimator_id),
        gpu_fit_parameters,
        parameters_to_fit,
        gpu_user_info,
        user_info_size,
        gpu_fit_parameters,
        gpu_output_states,
        gpu_output_chi_squares,
        gpu_output_n_iterations,
        DEVICE);

    fi.fit(static_cast<ModelID>(model_id));

    return ReturnState::OK;
}
catch (std::exception & exception)
{
    last_error = exception.what();

    return ReturnState::ERROR;
}
catch (...)
{
    last_error = "unknown error";

    return ReturnState::ERROR;
}

char const * gpufit_get_last_error()
{
    return last_error.c_str() ;
}

int gpufit_cuda_available()
{
	// Returns 1 if CUDA is available and 0 otherwise
	try
	{
		getDeviceCount();
		return 1;
	}
	catch (std::exception & exception)
	{
		last_error = exception.what();

		return 0;
	}
}

int gpufit_get_cuda_version(int * runtime_version, int * driver_version)
{
    try
    {
        cudaRuntimeGetVersion(runtime_version);
        cudaDriverGetVersion(driver_version);
        return ReturnState::OK;
    }
    catch (std::exception & exception)
    {
        last_error = exception.what();

        return ReturnState::ERROR;
    }
}

int gpufit_portable_interface(int argc, void *argv[])
{
// NOTE:
/*
0 int numberFits, 
1 int numberPoints, 
2 FloatBuffer data, 
3 FloatBuffer weights, 
4 int model_id, 
5 FloatBuffer initialParameters, 
6 float tolerance, 
7 int maxNumberIterations, 
8 int num_valid_coefs, 
9 IntBuffer parametersToFit, 
10 int estimatorID, 
11 int userInfoSize, 
12 FloatBuffer userInfo, 
13 FloatBuffer outputParameters, 
14 IntBuffer outputStates, 
15 FloatBuffer outputChiSquares, 
16 IntBuffer outputNumberIterations
*/
    return gpufit(
        *((size_t *) argv[0]),
        *((size_t *) argv[1]),
        (REAL *) argv[2],
        (REAL *) argv[3],
        *((int *) argv[4]),
        (REAL *) argv[5],
        *((REAL *) argv[6]),
        *((int *) argv[7]),
        *((int *) argv[8]),
        (int *) argv[9],
        *((int *) argv[10]),
        *((size_t *) argv[11]),
        (char *) argv[12],
        (REAL *) argv[13],
        (int *) argv[14],
        (REAL *) argv[15],
        (int *) argv[16]);

}
/* -------------------------------------------------------------------------------------------------------
* from gpufit.cpp END
------------------------------------------------------------------------------------------------------- */
/*
-------------------------------------------------------------------------------------------------------
* Calculate bleach correction input data
* from JCudaImageJExampleKernelcalcacf7.cu START
------------------------------------------------------------------------------------------------------- 
*/
__global__ void calc_data_bleach_correction(float* data, float* data1, int width, int height, int nopit, int ave)
{
    // function is an averaging step in temporal dimension for every ave number of points, prior to performing bleach correction fitting.

    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    __syncthreads();
    
    if ( (idx < width) && (idy < height) )
    {
        for (int z1 = 0; z1 < nopit; z1++) {
            double sum1 = 0;

            for (int yy = z1 * ave; yy < (z1  + 1) * ave; yy++) {
                sum1 += (float) data[yy * width * height + idy * width + idx];
            } // for yy
            data1[idy*width*nopit + idx * nopit + z1] = sum1/ave;
        } // for z1
    } //if
}
/* -------------------------------------------------------------------------------------------------------
* from JCudaImageJExampleKernelcalcacf7.cu END
------------------------------------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------------------------------------
* Calculate binning START
-------------------------------------------------------------------------------------------------------*/
__global__ void calc_binning(float* data, float* data1, int win_star, int hin_star, int w_temp, int h_temp, int framediff, int pixbinX, int pixbinY, int binningX, int binningY)
{
    // this function performs binning of spatial data.

    // NOTE: In the case overlap is OFF, we sill bin for every pixel, one pixel at a time.
    // This allows us to use cfXDistance and cfYDistance directly instead of translating these distances, which will be difficult.
 
    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    __syncthreads();

    float sum = 0.0;    

    if ( (idx < w_temp) && (idy < h_temp) )
    {
        for (int t = 0; t < framediff; t++) {
            sum = 0.0;
            for (int i = 0; i < binningX; i++) {
                for (int j = 0; j < binningY; j++) {
                    sum += data[t * win_star * hin_star + (idy + j) * win_star + (idx + i)];  
                } // for j
            } // for i
            
            data1[t * w_temp * h_temp + idy * w_temp + idx] = sum;

        } // for t
    } // if
}

/*-------------------------------------------------------------------------------------------------------
* Calculate binning END
-------------------------------------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------------------------------------
* from com_github_gpufit_Gpufit.cpp START
* NOTE: creates the gpufitJNI.dll. File is located at /Gpufit/Gpufit/java/adapter/
------------------------------------------------------------------------------------------------------- */
void * buffer_address(JNIEnv * env, jobject buffer)
{
    if (buffer == 0)
    {
        return 0;
    }
    else
    {
        return env->GetDirectBufferAddress(buffer);
    }
}

/*
* Calls gpufit(), no consistency checks on this side.
*
* Class:     gpufitImFCS_GpufitImFCS
* Method:    fit
* Signature: (IILjava/nio/FloatBuffer;Ljava/nio/FloatBuffer;ILjava/nio/FloatBuffer;FILjava/nio/IntBuffer;IILjava/nio/ByteBuffer;Ljava/nio/FloatBuffer;Ljava/nio/IntBuffer;Ljava/nio/FloatBuffer;Ljava/nio/IntBuffer;)I
*/
jint JNICALL Java_gpufitImFCS_GpufitImFCS_fit(JNIEnv * env, jclass cls, jint number_fits, jint number_points, jobject data_buffer, jobject weights_buffer, jint model_id, jobject initial_parameter_buffer, jfloat tolerance, jint max_number_iterations, jint num_valid_coefs, jobject paramters_to_fit_buffer, jint estimator_id, jint user_info_size, jobject user_info_buffer, jobject output_parameters_buffer, jobject output_states_buffer, jobject output_chi_squares_buffer, jobject output_number_iterations_buffer)
{
    // get pointer to buffers
    REAL * data = (REAL *)buffer_address(env, data_buffer);
    REAL * weights = (REAL *)buffer_address(env, weights_buffer);
    REAL * initial_parameters = (REAL *)buffer_address(env, initial_parameter_buffer);
    int * parameters_to_fit = (int *)buffer_address(env, paramters_to_fit_buffer);
    char * user_info = (char *)buffer_address(env, user_info_buffer);
    REAL * output_parameters = (REAL *)buffer_address(env, output_parameters_buffer);
    int * output_states = (int *)buffer_address(env, output_states_buffer);
    REAL * output_chi_squares = (REAL *)buffer_address(env, output_chi_squares_buffer);
    int * output_number_iterations = (int *)buffer_address(env, output_number_iterations_buffer);

    // call to gpufit
    // NOTE: Added num_valid_coefs
    int status = gpufit(number_fits, number_points, data, weights, model_id, initial_parameters, tolerance, max_number_iterations, num_valid_coefs, parameters_to_fit, estimator_id, user_info_size, user_info, output_parameters, output_states, output_chi_squares, output_number_iterations);

// lmfit_cuda

    return status;
}

/*
* Calls gpufit_get_last_error()
*
* Class:     gpufitImFCS_GpufitImFCS
* Method:    getLastError
* Signature: ()Ljava/lang/String;
*/
jstring JNICALL Java_gpufitImFCS_GpufitImFCS_getLastError(JNIEnv * env, jclass cls)
{
    char const * error = gpufit_get_last_error();
    return env->NewStringUTF(error);
}

/*
* Calls gpufit_cuda_available()
*
* Class:     gpufitImFCS_GpufitImFCS
* Method:    isCudaAvailableInt
* Signature: ()Z
*/
jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isCudaAvailableInt(JNIEnv * env, jclass cls)
{
    return gpufit_cuda_available() == 1 ? JNI_TRUE : JNI_FALSE;
}

/*
* Calls gpufit_get_cuda_version()
*
* Class:     gpufitImFCS_GpufitImFCS
* Method:    getCudaVersionAsArray
* Signature: ()[I
*/
jintArray JNICALL Java_gpufitImFCS_GpufitImFCS_getCudaVersionAsArray(JNIEnv * env, jclass cls)
{
    int runtime_version, driver_version;
    if (gpufit_get_cuda_version(&runtime_version, &driver_version) == ReturnState::OK)
    {
        // create int[2] in Java and fill with values
        jintArray array = env->NewIntArray(2);
        jint fill[2];
        fill[0] = runtime_version;
        fill[1] = driver_version;
        env->SetIntArrayRegion(array, 0, 2, fill);
        return array;
    }
    else
    {
        return 0;
    }
}

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    resetGPU
 * Signature: ()V
 */
void JNICALL Java_gpufitImFCS_GpufitImFCS_resetGPU(JNIEnv * env, jclass cls)
{
    try{
        cudaDeviceReset();
    } catch(std::runtime_error & e) {
        // see: https://www.rgagnon.com/javadetails/java-0323.html
        jclass Exception = env->FindClass("java/lang/Exception");
        env->ThrowNew(Exception, e.what());
    }
}

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcDataBleachCorrection
 * Signature: ([F[FLgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
void JNICALL Java_gpufitImFCS_GpufitImFCS_calcDataBleachCorrection(JNIEnv * env, jclass cls, jfloatArray pixels, jfloatArray outdata, jobject ACFInputParams)
{
    size_t SIZEFLOAT = sizeof(float);

    // input arrays required for calculations.
    jfloat *Cpixels;
    float *d_Cpixels;

    //output data
    float *Coutput;
    float *d_Coutput;

    try{
        Cpixels = env->GetFloatArrayElements(pixels, NULL);        

        // get parameters from the ACFInputParams object
        // we need width, height, cfXDistancegpu, cfYDistancegpu, nopit, ave
        // we also need firstframe and lastframe for setting blockSize and gridSize
        jclass ACFInputParamsCls = env->GetObjectClass(ACFInputParams);

        jfieldID w_tempId = env->GetFieldID(ACFInputParamsCls, "w_temp", "I");
        jfieldID h_tempId = env->GetFieldID(ACFInputParamsCls, "h_temp", "I");
        jfieldID firstframeId = env->GetFieldID(ACFInputParamsCls, "firstframe", "I");
        jfieldID lastframeId = env->GetFieldID(ACFInputParamsCls, "lastframe", "I");
        jfieldID nopitId = env->GetFieldID(ACFInputParamsCls, "nopit", "I");
        jfieldID aveId = env->GetFieldID(ACFInputParamsCls, "ave", "I");

        jint w_temp = env->GetIntField(ACFInputParams, w_tempId);
        jint h_temp = env->GetIntField(ACFInputParams, h_tempId);
        jint firstframe = env->GetIntField(ACFInputParams, firstframeId);
        jint lastframe = env->GetIntField(ACFInputParams, lastframeId);
        jint nopit = env->GetIntField(ACFInputParams, nopitId);
        jint ave = env->GetIntField(ACFInputParams, aveId);

        // blockSize and gridSize
        int framediff = lastframe - firstframe + 1;
        int BLKSIZEXY = 16;
        int a = ( w_temp > h_temp) ? w_temp : h_temp;
        int GRIDSIZEXY = (a + BLKSIZEXY -1) / BLKSIZEXY;

        dim3 blockSize(BLKSIZEXY, BLKSIZEXY, 1);
        dim3 gridSize(GRIDSIZEXY, GRIDSIZEXY, 1);

        // Allocate memory on GPU
        size_t size = w_temp * h_temp * framediff * SIZEFLOAT;
        cudaMalloc((void **)&d_Cpixels, size);

        // Allocate memory for Coutput and d_Coutput
        unsigned int sizeA = w_temp * h_temp * nopit;
        size_t size1 = sizeA * SIZEFLOAT;

        Coutput = (float *)malloc(size1);
        cudaMalloc((void **)&d_Coutput, size1);

        // Copy to GPU
        CUDA_CHECK_STATUS(cudaMemcpy(d_Cpixels, Cpixels, size, cudaMemcpyHostToDevice));

        cudaStream_t stream;
        CUDA_CHECK_STATUS(cudaStreamCreate( &stream ));

        calc_data_bleach_correction<<<gridSize, blockSize, 0, stream>>>(d_Cpixels, d_Coutput, w_temp, h_temp, nopit, ave);

        cudaDeviceSynchronize();
        CUDA_CHECK_STATUS(cudaGetLastError());

        // copy memory from device to host
        CUDA_CHECK_STATUS(cudaMemcpy(Coutput, d_Coutput, size1, cudaMemcpyDeviceToHost));

        CUDA_CHECK_STATUS(cudaStreamDestroy( stream ));

        //CUDA release memory
        cudaFree(d_Cpixels); cudaFree(d_Coutput);

        cudaDeviceReset();

        // copy values to Java output arrays.
        env->SetFloatArrayRegion(outdata, 0 , sizeA, Coutput);  

        //free pointers
        free(Coutput);

        // release resources
        env->ReleaseFloatArrayElements(pixels, Cpixels, 0);        
       
    } catch(std::runtime_error & e) {

        //CUDA release memory
        cudaFree(d_Cpixels); cudaFree(d_Coutput);

        //free pointers
        free(Coutput);

        // release resources
        env->ReleaseFloatArrayElements(pixels, Cpixels, 0);

        // see: https://www.rgagnon.com/javadetails/java-0323.html
        jclass Exception = env->FindClass("java/lang/Exception");
        env->ThrowNew(Exception, e.what());
    }

    return;
}

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    isBinningMemorySufficient
 * Signature: (LgpufitImFCS/GpufitImFCS/ACFParameters;)Z
 */
jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isBinningMemorySufficient(JNIEnv * env, jclass cls, jobject ACFInputParams)
{

  try {
      unsigned int SIZEFLOAT = sizeof(float);

      // get parameters from the ACFInputParams object
      jclass ACFInputParamsCls = env->GetObjectClass(ACFInputParams);
      jfieldID win_starId = env->GetFieldID(ACFInputParamsCls, "win_star", "I");
      jfieldID hin_starId = env->GetFieldID(ACFInputParamsCls, "hin_star", "I");
      jfieldID w_tempId = env->GetFieldID(ACFInputParamsCls, "w_temp", "I");
      jfieldID h_tempId = env->GetFieldID(ACFInputParamsCls, "h_temp", "I");
      jfieldID firstframeId = env->GetFieldID(ACFInputParamsCls, "firstframe", "I");
      jfieldID lastframeId = env->GetFieldID(ACFInputParamsCls, "lastframe", "I");

      jint win_star = env->GetIntField(ACFInputParams, win_starId);
      jint hin_star = env->GetIntField(ACFInputParams, hin_starId);
      jint w_temp = env->GetIntField(ACFInputParams, w_tempId);
      jint h_temp = env->GetIntField(ACFInputParams, h_tempId);
      jint firstframe = env->GetIntField(ACFInputParams, firstframeId);
      jint lastframe = env->GetIntField(ACFInputParams, lastframeId);

      // sanity check if memory on GPU is sufficient for binning
      std::size_t this_free_bytes;
      std::size_t this_total_bytes;    
      int framediff = lastframe - firstframe + 1;
      double maxmemory = (double) (win_star * hin_star + w_temp * h_temp) * framediff * SIZEFLOAT;
      CUDA_CHECK_STATUS(cudaMemGetInfo(&this_free_bytes, &this_total_bytes));

      return (maxmemory > double(this_free_bytes) * 0.9) ? JNI_FALSE : JNI_TRUE;

  } catch (std::runtime_error & e)  {
      jclass Exception = env->FindClass("java/lang/Exception");
      env->ThrowNew(Exception, e.what());
      return JNI_FALSE;
  }

}

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcBinning
 * Signature: ([F[FLgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
void JNICALL Java_gpufitImFCS_GpufitImFCS_calcBinning(JNIEnv * env, jclass cls, jfloatArray indata, jfloatArray outdata, jobject ACFInputParams)
{
    size_t SIZEFLOAT = sizeof(float);

    // input arrays required for calculations.
    jfloat *Cindata;
    float *d_Cindata;

    //output data
    float *Coutput;
    float *d_Coutput;

    try{
        Cindata = env->GetFloatArrayElements(indata, NULL);        

        // get parameters from the ACFInputParams object
        // we need w_temp, h_temp, binningX, binningY
        // we also need firstframe and lastframe for setting blockSize and gridSize
        jclass ACFInputParamsCls = env->GetObjectClass(ACFInputParams);
        jfieldID win_starId = env->GetFieldID(ACFInputParamsCls, "win_star", "I");
        jfieldID hin_starId = env->GetFieldID(ACFInputParamsCls, "hin_star", "I");
        jfieldID w_tempId = env->GetFieldID(ACFInputParamsCls, "w_temp", "I");
        jfieldID h_tempId = env->GetFieldID(ACFInputParamsCls, "h_temp", "I");
        jfieldID pixbinXId = env->GetFieldID(ACFInputParamsCls, "pixbinX", "I");
        jfieldID pixbinYId = env->GetFieldID(ACFInputParamsCls, "pixbinY", "I");
        jfieldID binningXId = env->GetFieldID(ACFInputParamsCls, "binningX", "I");
        jfieldID binningYId = env->GetFieldID(ACFInputParamsCls, "binningY", "I");
        jfieldID firstframeId = env->GetFieldID(ACFInputParamsCls, "firstframe", "I");
        jfieldID lastframeId = env->GetFieldID(ACFInputParamsCls, "lastframe", "I");

        jint win_star = env->GetIntField(ACFInputParams, win_starId);
        jint hin_star = env->GetIntField(ACFInputParams, hin_starId);
        jint w_temp = env->GetIntField(ACFInputParams, w_tempId);
        jint h_temp = env->GetIntField(ACFInputParams, h_tempId);
        jint pixbinX = env->GetIntField(ACFInputParams, pixbinXId);
        jint pixbinY = env->GetIntField(ACFInputParams, pixbinYId);
        jint binningX = env->GetIntField(ACFInputParams, binningXId);
        jint binningY = env->GetIntField(ACFInputParams, binningYId);
        jint firstframe = env->GetIntField(ACFInputParams, firstframeId);
        jint lastframe = env->GetIntField(ACFInputParams, lastframeId);

        // blockSize and gridSize
        int framediff = lastframe - firstframe + 1;
        int BLKSIZEXY = 16;
        int a = ( w_temp > h_temp) ? w_temp : h_temp;
        int GRIDSIZEXY = (a + BLKSIZEXY - 1) / BLKSIZEXY;

        dim3 blockSizeBin(BLKSIZEXY, BLKSIZEXY, 1);
        dim3 gridSizeBin(GRIDSIZEXY, GRIDSIZEXY, 1);

        // Allocate memory on GPU
        size_t size = win_star * hin_star * framediff * SIZEFLOAT;
        cudaMalloc((void **)&d_Cindata, size);

        // Allocate memory for Coutput and d_Coutput
        unsigned int sizeA = w_temp * h_temp * framediff;
        size_t size1 = sizeA * SIZEFLOAT;
        Coutput = (float *)malloc(size1);
        cudaMalloc((void **)&d_Coutput, size1);

        // Copy to GPU
        CUDA_CHECK_STATUS(cudaMemcpy(d_Cindata, Cindata, size, cudaMemcpyHostToDevice));

        cudaStream_t stream;
        CUDA_CHECK_STATUS(cudaStreamCreate( &stream ));

        calc_binning<<<gridSizeBin, blockSizeBin, 0, stream>>>(d_Cindata, d_Coutput, win_star, hin_star, w_temp, h_temp, framediff, pixbinX, pixbinY, binningX, binningY);

        cudaDeviceSynchronize();
        CUDA_CHECK_STATUS(cudaGetLastError());

        // copy memory from device to host
        CUDA_CHECK_STATUS(cudaMemcpy(Coutput, d_Coutput, size1, cudaMemcpyDeviceToHost));

        CUDA_CHECK_STATUS(cudaStreamDestroy( stream ));

        //CUDA release memory
        cudaFree(d_Cindata); cudaFree(d_Coutput);

        cudaDeviceReset();

        // copy values to Java output arrays.
        env->SetFloatArrayRegion(outdata, 0 , sizeA, Coutput);  
 
        //free pointers
        free(Coutput);

        // release resources
        env->ReleaseFloatArrayElements(indata, Cindata, 0);

    } catch(std::runtime_error & e) {
        //CUDA release memory
        cudaFree(d_Cindata); cudaFree(d_Coutput);

        //free pointers
        free(Coutput);

        // release resources
        env->ReleaseFloatArrayElements(indata, Cindata, 0);

        // see: https://www.rgagnon.com/javadetails/java-0323.html
        jclass Exception = env->FindClass("java/lang/Exception");
        env->ThrowNew(Exception, e.what());
    }

    return;

}

/*
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    isACFmemorySufficient
 * Signature: (LgpufitImFCS/GpufitImFCS/ACFParameters;)Z
 */
jboolean JNICALL Java_gpufitImFCS_GpufitImFCS_isACFmemorySufficient(JNIEnv * env, jclass cls, jobject ACFInputParams)
{
  try{
      unsigned int SIZEINT = sizeof(int);
      unsigned int SIZEFLOAT = sizeof(float);
      unsigned int SIZEDOUBLE = sizeof(double);

      double totalmemoryAll = 0.0;
      double totalmemoryCalc3 = 0.0;
      double totalmemoryCalc2 = 0.0;

      // get parameters that are required for the ACF calculations from the ACFInputParams object
      jclass ACFInputParamsCls = env->GetObjectClass(ACFInputParams);

      jfieldID widthId = env->GetFieldID(ACFInputParamsCls, "width", "I");
      jfieldID heightId = env->GetFieldID(ACFInputParamsCls, "height", "I");
      jfieldID w_tempId = env->GetFieldID(ACFInputParamsCls, "w_temp", "I");
      jfieldID h_tempId = env->GetFieldID(ACFInputParamsCls, "h_temp", "I");
      jfieldID firstframeId = env->GetFieldID(ACFInputParamsCls, "firstframe", "I");
      jfieldID lastframeId = env->GetFieldID(ACFInputParamsCls, "lastframe", "I");
      jfieldID cfXDistanceId = env->GetFieldID(ACFInputParamsCls, "cfXDistance", "I");
      jfieldID cfYDistanceId = env->GetFieldID(ACFInputParamsCls, "cfYDistance", "I");
//      jfieldID correlatorpId = env->GetFieldID(ACFInputParamsCls, "correlatorp", "D");
//      jfieldID correlatorqId = env->GetFieldID(ACFInputParamsCls, "correlatorq", "D");
      jfieldID frametimeId = env->GetFieldID(ACFInputParamsCls, "frametime", "D");
      jfieldID backgroundId = env->GetFieldID(ACFInputParamsCls, "background", "I");
      jfieldID mtab1Id = env->GetFieldID(ACFInputParamsCls, "mtab1", "D");
      jfieldID mtabchanumminus1Id = env->GetFieldID(ACFInputParamsCls, "mtabchanumminus1", "D");
      jfieldID sampchanumminus1Id = env->GetFieldID(ACFInputParamsCls, "sampchanumminus1", "D");
      jfieldID chanumId = env->GetFieldID(ACFInputParamsCls, "chanum", "I");
      jfieldID isNBcalculationId = env->GetFieldID(ACFInputParamsCls, "isNBcalculation", "Z");
      jfieldID bleachcorr_gpuId = env->GetFieldID(ACFInputParamsCls, "bleachcorr_gpu", "Z");
      jfieldID bleachcorr_orderId = env->GetFieldID(ACFInputParamsCls, "bleachcorr_order", "I");

      jint width = env->GetIntField(ACFInputParams, widthId);
      jint height = env->GetIntField(ACFInputParams, heightId);
      jint w_temp = env->GetIntField(ACFInputParams, w_tempId);
      jint h_temp = env->GetIntField(ACFInputParams, h_tempId);
      jint firstframe = env->GetIntField(ACFInputParams, firstframeId);
      jint lastframe = env->GetIntField(ACFInputParams, lastframeId);
      jint cfXDistance = env->GetIntField(ACFInputParams, cfXDistanceId);
      jint cfYDistance = env->GetIntField(ACFInputParams, cfYDistanceId);
//      jdouble correlatorpdbl = env->GetDoubleField(ACFInputParams, correlatorpId);
//      jdouble correlatorqdbl = env->GetDoubleField(ACFInputParams, correlatorqId);
      jdouble frametime = env->GetDoubleField(ACFInputParams, frametimeId);
      jint background = env->GetIntField(ACFInputParams, backgroundId);
      jdouble mtab1 = env->GetDoubleField(ACFInputParams, mtab1Id); // mtab[1], used to calculate blocknumgpu.
      jdouble mtabchanumminus1 = env->GetDoubleField(ACFInputParams, mtabchanumminus1Id); // mtab[chanum-1], used to calculate pnumgpu[counter_indexarray]
      jdouble sampchanumminus1 = env->GetDoubleField(ACFInputParams, sampchanumminus1Id); // samp[chanum-1], used to calculate pnumgpu[counter_indexarray]
      jint chanum = env->GetIntField(ACFInputParams, chanumId);
      jboolean isNBcalculation = env->GetBooleanField(ACFInputParams, isNBcalculationId);      
      jboolean bleachcorr_gpu = env->GetBooleanField(ACFInputParams, bleachcorr_gpuId);      
      jint bleachcorr_order = env->GetIntField(ACFInputParams, bleachcorr_orderId);

      // initialize parameters
//      int correlatorp = (int) correlatorpdbl;
//      int correlatorq = (int) correlatorqdbl;
      int framediff = lastframe - firstframe + 1;
      unsigned long size = w_temp * h_temp * framediff * SIZEFLOAT;
      unsigned long size1 = width * height * chanum * SIZEDOUBLE;
      unsigned long size2 = framediff * width * height * SIZEFLOAT;
      unsigned long sizeblockvararray = chanum * width * height * SIZEDOUBLE;
          
      int blocknumgpu = (int) (floor(log(mtab1)/log(2)) - 2);

      // dynamic memory allocation and/or initialization
      //------------------ common parameters ---------------------------
      totalmemoryAll = totalmemoryAll + (double) size1; //Cpixels1
      totalmemoryAll = totalmemoryAll + (double) size2; //prod
      totalmemoryAll = totalmemoryAll + (double) size; //pixels
      totalmemoryAll = totalmemoryAll + (double) (chanum * SIZEDOUBLE); //samp
      totalmemoryAll = totalmemoryAll + (double) (chanum * SIZEINT); //lag
      totalmemoryAll = totalmemoryAll + (double) width * height * bleachcorr_order * SIZEDOUBLE; 
      // totalmemoryAll = totalmemoryAll + (double) size1; // Cblocked1D copies Cpixels1 array after calcacf3 calculation and not required GPU memory.

      //------------------ calcacf3 ---------------------------
      if (!isNBcalculation) {
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu * SIZEDOUBLE); //prodnum
          totalmemoryCalc3 = totalmemoryCalc3 + (double) sizeblockvararray; //blocksd
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu*width*height*SIZEDOUBLE); //upper
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu*width*height*SIZEDOUBLE); //lower
          totalmemoryCalc3 = totalmemoryCalc3 + (double) ((blocknumgpu-1)*width*height*SIZEINT); //crt
          totalmemoryCalc3 = totalmemoryCalc3 + (double) ((blocknumgpu-2)*width*height*SIZEINT); //cr12
          totalmemoryCalc3 = totalmemoryCalc3 + (double) ((blocknumgpu-2)*width*height*SIZEINT); //cr3
          totalmemoryCalc3 = totalmemoryCalc3 + (double) ((blocknumgpu-1)*width*height*SIZEINT); //diffpos
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu*width*height*SIZEDOUBLE); //varblock0
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu*width*height*SIZEDOUBLE); //varblock1
          totalmemoryCalc3 = totalmemoryCalc3 + (double) (blocknumgpu*width*height*SIZEDOUBLE); //varblock2
      }

      //------------------ calcacf2 ---------------------------
      totalmemoryCalc2 = totalmemoryCalc2 + (double) (chanum*width*height*SIZEINT); //prodnumarray
      totalmemoryCalc2 = totalmemoryCalc2 + (double) (width*height*SIZEINT); //indexarray
      totalmemoryCalc2 = totalmemoryCalc2 + (double) sizeblockvararray; //Cblockvararray
      totalmemoryCalc2 = totalmemoryCalc2 + (double) sizeblockvararray; //blocksdarray
      totalmemoryCalc2 = totalmemoryCalc2 + (double) (width*height*SIZEINT); //pnumgpu
 
      //------------------ calculation of N & B in calcacf2 --------------------------
      if (isNBcalculation) {
          totalmemoryCalc2 = totalmemoryCalc2 + (double) ( width * height * SIZEDOUBLE * 3); //NBmeanGPU, NBmean2GPU, NBcovarianceGPU
      }

      // sanity check if memory on GPU is sufficient for calcacf3 and calcacf2
      // NOTE calcacf3 will run first. Memory of parameters related to calcacf3 only will be released after completion of calcacf3.
      std::size_t this_free_bytes;
      std::size_t this_total_bytes;
      double maxmemory = (totalmemoryCalc3 > totalmemoryCalc2)? totalmemoryCalc3 : totalmemoryCalc2;
      maxmemory = maxmemory + totalmemoryAll;
      CUDA_CHECK_STATUS(cudaMemGetInfo(&this_free_bytes, &this_total_bytes));

      return (maxmemory > double(this_free_bytes) * 0.9) ? JNI_FALSE : JNI_TRUE;

  } catch (std::runtime_error & e)  {
    jclass Exception = env->FindClass("java/lang/Exception");
    env->ThrowNew(Exception, e.what());
    return JNI_FALSE;
  }

}

/* ------------------------------------------
AUTOCORRELATION SINGLE DIMENSION ARRAY CALCULATION START
NOTE: SEE JCudaImageJExampleKernelcalcacf2.cu
------------------------------------------ */
__global__ void bleachcorrection(float * data, int w_temp, int h_temp, int d, int bleachcorr_order, double frametimegpu, double* bleachcorr_params)
{
    // function performs polynomial bleach correction given polynomial order and coefficients. It is done prior to calcacf3, calcacf2a and calcacf2b.

    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    __syncthreads();
    if ( (idx < w_temp) && (idy < h_temp) ) {
        for (int i = 0; i < d; i++) {
            float corfunc = 0;
            for (int ii = 0; ii < bleachcorr_order; ii++) {
                  corfunc += bleachcorr_params[(idy*w_temp +idx)*bleachcorr_order + ii] * powf((float)frametimegpu * (i + 0.5), (float)ii);
            } // for ii
            
            float res0 = bleachcorr_params[(idy*w_temp + idx) * bleachcorr_order];

	    data[i*w_temp*h_temp + idy*w_temp + idx] = data[i*w_temp*h_temp + idy*w_temp + idx] / sqrtf(corfunc / res0) + res0 * (1 - sqrtf(corfunc / res0));
            __syncthreads();
        } // for i
    } // if ((idx < w_temp) && (idy < h_temp))
} //bleachcorrection function

__device__ unsigned int countx = 0;
__device__ unsigned int county = 0;
__shared__ bool isLastBlockDone;
__global__ void calcacf2a(float* data, int w_temp, int h_temp, int numbin)
{
    // function calculates the arrays according to different time bins in different parts of the correlation function

    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    
    __syncthreads();
    if ( (idx < w_temp) && (idy < h_temp) ) {
        // And correct the number Of actual data points accordingly
        for (int y = 0; y < numbin; y++){ // if yes, bin the data according to the width of the current channel
            data[y*w_temp*h_temp + idy*w_temp + idx] = data[2 * y*w_temp*h_temp + idy*w_temp + idx] + data[(2 * y + 1)*w_temp*h_temp + idy*w_temp + idx];
            __syncthreads();            
        } // for int y = 0 
    }
}

__global__ void calcacf2b(float* data, int cfXDistancegpu, int cfYDistancegpu, int w, int h, int w_temp, int h_temp, int pixbinX, int pixbinY, double* data1, float* prod, int*laggpu, int* prodnumarray, int* indexarray, double* blockvararray, double* sdarray, int* pnumgpu, int x, int numbin, int currentIncrement, int ctbin, bool isNBcalculation, double* NBmeanGPU, double* NBcovarianceGPU)
{
    // function calculates the value of the auto or cross-correlation at every lag time. This function also performs the G1 analysis in N and B calculation.

    int del;				// delay Or correlation time expressed In lags
    double sumprod = 0.0;		// sum of all intensity products; divide by num to get the average <i(n)i(n+del)>

    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    double temp1=0.0, temp2=0.0;

    __syncthreads();
    if ( (idx < w) && (idy < h) ){

        del = laggpu[x] / currentIncrement;
	prodnumarray[x*w*h + idy*w + idx] = numbin - del;

	temp1 = 0.0;
	temp2 = 0.0;
	    
        for (int y = 0; y < prodnumarray[x*w*h + idy*w + idx]; y++){ // calculate the ...
	    temp1 += data[y*w_temp*h_temp + idy*pixbinY*w_temp +  idx*pixbinX];
            temp2 += data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)];
        }

	temp1 /= prodnumarray[x*w*h + idy*w + idx]; // calculate average of direct and delayed monitor, i.e. the average intensity <n(0)> and <n(tau)>
        temp2 /= prodnumarray[x*w*h + idy*w + idx];
	sumprod = 0.0;

	for (int y = 0; y < prodnumarray[x*w*h + idy*w + idx]; y++){ // calculate the correlation
            if (isNBcalculation) {
               prod[y*w*h + idy*w + idx] =  data[y*w_temp*h_temp + idy*pixbinY*w_temp + idx*pixbinX] * data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)];
           
            } else {
                prod[y*w*h + idy*w + idx] =  data[y*w_temp*h_temp + idy*pixbinY*w_temp + idx*pixbinX] * data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)] - temp2 * data[y*w_temp*h_temp + idy*pixbinY*w_temp +  idx*pixbinX] - temp1 * data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)] + temp1*temp2;
            }
            sumprod += prod[y*w*h + idy*w + idx];	
        }

        if (isNBcalculation) {
            NBmeanGPU[idy * w + idx] = temp1;
            NBcovarianceGPU[idy * w + idx] = sumprod/prodnumarray[x*w*h + idy*w + idx] -temp1*temp2;
        }

	__syncthreads();
 
        if (!isNBcalculation) {
            data1[x*w*h + idy*w + idx] = sumprod / (prodnumarray[x*w*h + idy*w + idx] * temp1 * temp2);
	    __syncthreads();

            sumprod = 0.0;
            double sumprod2 = 0.0;	// sum of all intensity products squared; divide by num to get the average <(i(n)i(n+del))^2>
            int binct = indexarray[idy*w + idx]-ctbin;
            double tempvariable =0.0;
            for (int y = 1; y <=binct; y++) {
                prodnumarray[x*w*h + idy*w + idx] = (int)floor((double)prodnumarray[x*w*h + idy*w + idx] / 2.0);
							 
                for (int z = 0; z < prodnumarray[x*w*h + idy*w + idx]; z++) {
	            prod[z*w*h + idy*w + idx]  = (prod[2 * z*w*h + idy*w + idx] + prod[(2 * z + 1)*w*h + idy*w + idx]) / 2.0;
                    __syncthreads();
	        }
            }
 			 	
	    for (int z = 0; z < pnumgpu[idy*w + idx]; z++) {
                tempvariable =  prod[z*w*h + idy*w + idx];
                sumprod += tempvariable; // calculate the sum of prod, i.e. the raw correlation value ...
                sumprod2 += powf(tempvariable, 2.0); // ... and the sum of the squares
            }

            blockvararray[x*w*h + idy*w + idx] = (sumprod2 / pnumgpu[idy*w + idx] - powf(sumprod / pnumgpu[idy*w + idx], 2.0)) / ((pnumgpu[idy*w + idx] - 1) * powf(temp1 * temp2, 2.0));
					
            sdarray[x*w*h + idy*w + idx] = sqrt(blockvararray[x*w*h + idy*w + idx]);
       } // if (!isNBcalculation)

    } // if ((idx < w) && (idy < h))
} //calcacf2b

/* ------------------------------------------
NOTE: SEE JCudaImageJExampleKernelcalcacf3.cu
------------------------------------------ */
__global__ void calcacf3(float* data, int cfXDistancegpu, int cfYDistancegpu, int blocklag, int w, int h, int w_temp, int h_temp, int pixbinX, int pixbinY, int d, int correlatorp, int correlatorq, int chanum, double frametimegpu, double* data1, float* prod, double* prodnum, double* blocksd, double* upper, double* lower, int* crt, int* cr12, int* cr3, int* diffpos, double* varblock0, double* varblock1, double* varblock2, double* sampgpu, int*laggpu)
{
    // this function calculates the block transformation values of the intensity.

    int blocknumgpu = (int)floor(log((double)d - 1.0) / log(2.0)) - 2;
    int numbin = d;		// number Of data points When they are binned
    int del;				// delay Or correlation time expressed In lags
    int currentIncrement = blocklag;
    double sumprod = 0.0;		// sum of all intensity products; divide by num to get the average <i(n)i(n+del)>
    double sumprod2 = 0.0;	// sum of all intensity products squared; divide by num to get the average <(i(n)i(n+del))^2>
    double directm = 0.0;		// direct monitor required for ACF normalization
    double delayedm = 0.0;
    int ind = 0;
    int last0 = 0;
    int idx = blockIdx.x * blockDim.x + threadIdx.x, idy = blockIdx.y * blockDim.y + threadIdx.y;
    int blockIndS = 0;

    __syncthreads();
    if ((idx < w) && (idy < h)){

        int  x = 1;
	del = laggpu[x] / currentIncrement; // calculate the delay, i.e. the correlation time
	for (int y = 0; y < numbin - del; y++) { // calculate the ...
            directm += data[y*w_temp*h_temp + idy*pixbinY*w_temp + idx*pixbinX]; // direct And ...
	    delayedm += data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)]; // delayed monitor
	}
	prodnum[0] = numbin - del; // number Of correlation products
	directm /= prodnum[0]; // calculate average Of direct And delayed monitor, 
	delayedm /= prodnum[0]; // i.e. the average intesity <n(0)> And <n(tau)>

	for (int y = 0; y < prodnum[0]; y++) { // calculate the correlation
	    prod[y*w*h + idy*w + idx] = data[y*w_temp*h_temp + idy*pixbinY*(w + cfXDistancegpu) + idx*pixbinX] * data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)] - delayedm * data[y*w_temp*h_temp + idy*pixbinY*w_temp + idx*pixbinX] - directm * data[(y + del)*w_temp*h_temp + (idy*pixbinY + cfYDistancegpu)*w_temp + (idx*pixbinX + cfXDistancegpu)] + delayedm * directm;
            __syncthreads();
	    sumprod += prod[y*w*h + idy*w + idx]; // calculate the sum Of prod, i.e. the raw correlation value ...
	    sumprod2 += powf(prod[y*w*h + idy*w + idx], 2.0); // ... And the sum Of the squares
	}

	varblock0[idy*w + idx] = currentIncrement * frametimegpu; // the time Of the block curve
	varblock1[idy*w + idx] = (sumprod2 / prodnum[0] - powf(sumprod / prodnum[0], 2.0)) / (prodnum[0] * powf(directm * delayedm, 2.0));

	for (int y = 1; y < blocknumgpu; y++) { // perform blocking operations
	    prodnum[y] = (int)floor((double)prodnum[y - 1] / 2);	// the number Of samples For the blocking curve decreases by a factor 2 With every Step
	    sumprod = 0;
	    sumprod2 = 0;
	    for (int z = 0; z < prodnum[y]; z++) { // bin the correlation data And calculate the blocked values for the SD
	        prod[z*w*h + idy*w + idx] = (prod[2 * z*w*h + idy*w + idx] + prod[(2 * z + 1)*w*h + idy*w + idx]) / 2;
                __syncthreads();
		sumprod += prod[z*w*h + idy*w + idx];
		sumprod2 += powf(prod[z*w*h + idy*w + idx], 2.0);
	    }
	
            // This is the correct one
	    varblock0[y*w*h + idy*w + idx] = (currentIncrement * powf(2, (double)y)) * frametimegpu;	// the time Of the block curve
	    varblock1[y*w*h + idy*w + idx] = (sumprod2 / prodnum[y] - powf(sumprod / prodnum[y], 2.0)) / (prodnum[y] * powf(directm * delayedm, 2.0));	// value of the block curve
	}

	for (int x = 0; x < blocknumgpu; x++) {
	    varblock1[x*w*h + idy*w + idx] = sqrt(varblock1[x*w*h + idy*w + idx]); // calculate the standard deviation
            varblock2[x*w*h + idy*w + idx] = varblock1[x*w*h + idy*w + idx] / sqrt((double)2 * (prodnum[x] - 1)); // calculate the error 
            __syncthreads();
	    upper[x*w*h + idy*w + idx] = varblock1[x*w*h + idy*w + idx] + varblock2[x*w*h + idy*w + idx]; // upper and lower quartile
	    lower[x*w*h + idy*w + idx] = varblock1[x*w*h + idy*w + idx] - varblock2[x*w*h + idy*w + idx];
	}

        // determine index where blocking criteria are fulfilled
	for (int x = 0; x < blocknumgpu - 1; x++) { // do neighboring points have overlapping error bars?
            if (upper[x*w*h + idy*w + idx] > lower[(x + 1)*w*h + idy*w + idx] && upper[(x + 1)*w*h + idy*w + idx] > lower[x*w*h + idy*w + idx]) {
                crt[x*w*h + idy*w + idx] = 1;
            }
        }

	for (int x = 0; x < blocknumgpu - 2; x++) { // do three adjacent points have overlapping error bars?
            if (crt[x*w*h + idy*w + idx] * crt[(x + 1)*w*h + idy*w + idx] == 1) {
                cr12[x*w*h + idy*w + idx] = 1;
            }
        }

	for (int x = 0; x < blocknumgpu - 1; x++) { // do neighboring points have a positive difference (increasing SD)?
            if (varblock1[(x + 1)*w*h + idy*w + idx] - varblock1[x*w*h + idy*w + idx] > 0) {
                diffpos[x*w*h + idy*w + idx] = 1;
            }
        }

        for (int x = 0; x < blocknumgpu - 2; x++) { // do three neighboring points monotonically increase?
            if (diffpos[x*w*h + idy*w + idx] * diffpos[(x + 1)*w*h + idy*w + idx] == 1) {
                cr3[x*w*h + idy*w + idx] = 1;
            }
        }

        for (int x = 0; x < blocknumgpu - 2; x++) { // find the last triple of points with monotonically increasing differences and non-overlapping error bars
            if ((cr3[x*w*h + idy*w + idx] == 1 && cr12[x*w*h + idy*w + idx] == 0)) {
                last0 = x;
            }
        }

        for (int x = 0; x <= last0; x++) { // indices of two pairs that pass criterion 1 an 2
          cr12[x*w*h + idy*w + idx] = 0;
        }

        cr12[(blocknumgpu - 3)*w*h + idy*w + idx] = 0; // criterion 3, the last two points can't be part of the blocking triple
        cr12[(blocknumgpu - 4)*w*h + idy*w + idx] = 0;

        for (int x = blocknumgpu - 5; x > 0; x--) { // index of triplet with overlapping error bars and after which no other triplet has a significant monotonic increase
            if (cr12[x*w*h + idy*w + idx] == 1) {											// or 4 increasing points
                ind = x + 1;
            }
        }

        if (ind == 0) { // if optimal blocking is not possible, use maximal blocking
            blockIndS = 0;
            if (blocknumgpu - 3 > 0) {
                ind = blocknumgpu - 3;
	    } else {
                ind = blocknumgpu - 1;
            }
        } else {
	    blockIndS = 1;
        }

        ind = (int)fmax((double)ind, (double)correlatorq - 1);
        data1[idy*w + idx] = (double)ind;
        data1[w*h + idy*w + idx] = (double)blockIndS;

    } // if ((idx < w) && (idy < h))
} // calcacf3

// initialize/set all values in the int array to value
void initializeintarr(int * array, unsigned long size, int value){
    for (unsigned long i = 0; i < size; i++){
        *(array+i) = value;
    }
}

// initialize/set all values in the double array to value
void initializedoublearr(double * array, unsigned long size, double value){
    for (unsigned long i = 0; i < size; i++){
        *(array+i) = value;
    }
}

/*
 * Calls calcacf3() and calcacf2()
 *
 * Class:     gpufitImFCS_GpufitImFCS
 * Method:    calcACF
 * Signature: ([F[D[D[D[D[D[ILgpufitImFCS/GpufitImFCS/ACFParameters;)V
 */
void JNICALL Java_gpufitImFCS_GpufitImFCS_calcACF(JNIEnv * env, jclass cls, jfloatArray pixels, jdoubleArray pixels1, jdoubleArray blockvararray, jdoubleArray NBmeanGPU, jdoubleArray NBcovarianceGPU, jdoubleArray blocked1D, jdoubleArray bleachcorr_params, jdoubleArray samp, jintArray lag, jobject ACFInputParams){
    // NOTE: outputs are stored in pixels1 and blockvararray arrays.
    // NOTE: Cpixels1 and Cblockvararray are temporary arrays to store output values before passing them to Java output arrays, by reference.

    // host arrays
    //------------------ common ---------------------------
    jfloat *Cpixels;
    double *Cpixels1;
    float *prod;
    double *Cbleachcorr_params;
    jdouble *Csamp;
    jint *Clag;
    double *Cblocked1D; // to copy Cpixels1 after calcacf3
    jboolean isNBcalculation;

    //------------------ calcacf3 ---------------------------      
    double *prodnum;
    double *blocksd;
    double *upper;
    double *lower;
    int *crt;
    int *cr12;
    int *cr3;
    int *diffpos;
    double *varblock0;
    double *varblock1;
    double *varblock2;

    //------------------ calcacf2 ---------------------------
    int *prodnumarray;
    int *indexarray;
    double *Cblockvararray;
    double *blocksdarray;
    int *pnumgpu;

    //--------------- N & B calculations in calcacf2 --------------------
    double *CNBmeanGPU;
    double *CNBcovarianceGPU;

    // CUDA arrays
    //------------------ common ---------------------------
    float *d_Cpixels;
    double *d_Cpixels1;
    double *d_Cbleachcorr_params;
    double *d_Csamp;
    int *d_Clag; 
    float *d_prod;

    //------------------ calcacf3 ---------------------------      
    double *d_prodnum;
    double *d_blocksd;
    double *d_upper;
    double *d_lower;
    int *d_crt;
    int *d_cr12;
    int *d_cr3;
    int *d_diffpos;
    double *d_varblock0;
    double *d_varblock1;
    double *d_varblock2;

    //------------------ calcacf2 ---------------------------
    int *d_prodnumarray;
    int *d_indexarray;
    double *d_Cblockvararray;
    double *d_blocksdarray;
    int *d_pnumgpu;

    //--------------- N & B calculations in calcacf2 --------------------
    double *d_NBmeanGPU;
    double *d_NBcovarianceGPU;

    try{
      int device = 0;
      CUDA_CHECK_STATUS(cudaSetDevice(device));	

      //reference: https://devblogs.nvidia.com/how-overlap-data-transfers-cuda-cc/
      // NOTE: using stream to control synchronous processing of calcacf2 and calcacf3 kernels.
      // based on the reference, we can potentially use it to speed up the transfer of data process.
      cudaStream_t stream;
      CUDA_CHECK_STATUS(cudaStreamCreate( &stream ));

      size_t SIZEINT = sizeof(int);
      size_t SIZEFLOAT = sizeof(float);
      size_t SIZEDOUBLE = sizeof(double);

      // input arrays required for calculations.
      Cpixels = env->GetFloatArrayElements(pixels, NULL);
      Cbleachcorr_params = env->GetDoubleArrayElements(bleachcorr_params, NULL);
      Csamp = env->GetDoubleArrayElements(samp, NULL);
      Clag = env->GetIntArrayElements(lag, NULL);

//      jint lensamp = env->GetArrayLength(samp);
//      jint lenlag = env->GetArrayLength(lag);
    
      // get parameters that are required for the ACF calculations from the ACFInputParams object
      jclass ACFInputParamsCls = env->GetObjectClass(ACFInputParams);

      jfieldID widthId = env->GetFieldID(ACFInputParamsCls, "width", "I");
      jfieldID heightId = env->GetFieldID(ACFInputParamsCls, "height", "I");
      jfieldID w_tempId = env->GetFieldID(ACFInputParamsCls, "w_temp", "I");
      jfieldID h_tempId = env->GetFieldID(ACFInputParamsCls, "h_temp", "I");
      jfieldID pixbinXId = env->GetFieldID(ACFInputParamsCls, "pixbinX", "I");
      jfieldID pixbinYId = env->GetFieldID(ACFInputParamsCls, "pixbinY", "I");
      jfieldID firstframeId = env->GetFieldID(ACFInputParamsCls, "firstframe", "I");
      jfieldID lastframeId = env->GetFieldID(ACFInputParamsCls, "lastframe", "I");
      jfieldID cfXDistanceId = env->GetFieldID(ACFInputParamsCls, "cfXDistance", "I");
      jfieldID cfYDistanceId = env->GetFieldID(ACFInputParamsCls, "cfYDistance", "I");
      jfieldID correlatorpId = env->GetFieldID(ACFInputParamsCls, "correlatorp", "D");
      jfieldID correlatorqId = env->GetFieldID(ACFInputParamsCls, "correlatorq", "D");
      jfieldID frametimeId = env->GetFieldID(ACFInputParamsCls, "frametime", "D");
      jfieldID backgroundId = env->GetFieldID(ACFInputParamsCls, "background", "I");
      jfieldID mtab1Id = env->GetFieldID(ACFInputParamsCls, "mtab1", "D");
      jfieldID mtabchanumminus1Id = env->GetFieldID(ACFInputParamsCls, "mtabchanumminus1", "D");
      jfieldID sampchanumminus1Id = env->GetFieldID(ACFInputParamsCls, "sampchanumminus1", "D");
      jfieldID chanumId = env->GetFieldID(ACFInputParamsCls, "chanum", "I");
      jfieldID isNBcalculationId = env->GetFieldID(ACFInputParamsCls, "isNBcalculation", "Z");
      jfieldID bleachcorr_gpuId = env->GetFieldID(ACFInputParamsCls, "bleachcorr_gpu", "Z");
      jfieldID bleachcorr_orderId = env->GetFieldID(ACFInputParamsCls, "bleachcorr_order", "I");

      jint width = env->GetIntField(ACFInputParams, widthId);
      jint height = env->GetIntField(ACFInputParams, heightId);
      jint w_temp = env->GetIntField(ACFInputParams, w_tempId);
      jint h_temp = env->GetIntField(ACFInputParams, h_tempId);
      jint pixbinX = env->GetIntField(ACFInputParams, pixbinXId);
      jint pixbinY = env->GetIntField(ACFInputParams, pixbinYId);
      jint firstframe = env->GetIntField(ACFInputParams, firstframeId);
      jint lastframe = env->GetIntField(ACFInputParams, lastframeId);
      jint cfXDistance = env->GetIntField(ACFInputParams, cfXDistanceId);
      jint cfYDistance = env->GetIntField(ACFInputParams, cfYDistanceId);
      jdouble correlatorpdbl = env->GetDoubleField(ACFInputParams, correlatorpId);
      jdouble correlatorqdbl = env->GetDoubleField(ACFInputParams, correlatorqId);
      jdouble frametime = env->GetDoubleField(ACFInputParams, frametimeId);
      jint background = env->GetIntField(ACFInputParams, backgroundId);
      jdouble mtab1 = env->GetDoubleField(ACFInputParams, mtab1Id); // mtab[1], used to calculate blocknumgpu.
      jdouble mtabchanumminus1 = env->GetDoubleField(ACFInputParams, mtabchanumminus1Id); // mtab[chanum-1], used to calculate pnumgpu[counter_indexarray]
      jdouble sampchanumminus1 = env->GetDoubleField(ACFInputParams, sampchanumminus1Id); // samp[chanum-1], used to calculate pnumgpu[counter_indexarray]
      jint chanum = env->GetIntField(ACFInputParams, chanumId);
      isNBcalculation = env->GetBooleanField(ACFInputParams, isNBcalculationId);      
      jboolean bleachcorr_gpu = env->GetBooleanField(ACFInputParams, bleachcorr_gpuId);      
      jint bleachcorr_order = env->GetIntField(ACFInputParams, bleachcorr_orderId);

      // initialize parameters
      int correlatorp = (int) correlatorpdbl;
      int correlatorq = (int) correlatorqdbl;
      int framediff = lastframe - firstframe + 1;
      size_t size = w_temp * h_temp * framediff * SIZEFLOAT;
      size_t size1 = width * height * chanum * SIZEDOUBLE;

      int blocklaggpu = 1;

      size_t size2 = framediff * width * height * SIZEFLOAT;
      size_t sizeblockvararray = chanum * width * height * SIZEDOUBLE;
          
      int blocknumgpu = (int) (floor(log(mtab1)/log(2)) - 2);

      // blockSize and gridSize
      int BLKSIZEXY = 16;
      int a = ( width > height ) ? width : height;
      int GRIDSIZEXY = (a + BLKSIZEXY -1) / BLKSIZEXY;

      dim3 blockSize(BLKSIZEXY, BLKSIZEXY, 1);
      dim3 gridSize(GRIDSIZEXY, GRIDSIZEXY, 1);

      int b = ( w_temp > h_temp ) ? w_temp : h_temp;
      int GRIDSIZEXY_Input = (b + BLKSIZEXY -1) / BLKSIZEXY;
      
      dim3 gridSize_Input(GRIDSIZEXY_Input, GRIDSIZEXY_Input, 1);

      // dynamic memory allocation and/or initialization
      //------------------ common parameters ---------------------------
      Cpixels1 = (double *)malloc(size1);
      prod = (float *)malloc(size2);
      Cblocked1D = (double *)malloc(size1); // Cblocked1D copies Cpixels1 array after calcacf3 calculation and not required GPU memory.

      //------------------ calcacf3 ---------------------------
      // Using if (!isNBcalculation) to comment out this section will cause warnings "may be used uninitialized in this function" being shown during compilation
      // instead, we will still initialize the arrays with malloc, but skip the initializedoublearr, initializeintarr functions.
      prodnum = (double *)malloc(blocknumgpu * SIZEDOUBLE);
      blocksd = (double *)malloc(sizeblockvararray);
      upper = (double *)malloc(blocknumgpu*width*height*SIZEDOUBLE);
      lower = (double *)malloc(blocknumgpu*width*height*SIZEDOUBLE);
      crt = (int *)malloc((blocknumgpu-1)*width*height*SIZEINT);
      cr12 = (int *)malloc((blocknumgpu-2)*width*height*SIZEINT);
      cr3 = (int *)malloc((blocknumgpu-2)*width*height*SIZEINT);
      diffpos = (int *)malloc((blocknumgpu-1)*width*height*SIZEINT);
      varblock0 = (double *)malloc(blocknumgpu*width*height*SIZEDOUBLE);
      varblock1 = (double *)malloc(blocknumgpu*width*height*SIZEDOUBLE);
      varblock2 = (double *)malloc(blocknumgpu*width*height*SIZEDOUBLE);

      if (!isNBcalculation) { 
          initializedoublearr(prodnum, blocknumgpu, 1.0); // initialize all values in prodnum to 1.0
          initializedoublearr(blocksd, chanum * width * height, 0.0);
          initializedoublearr(upper, blocknumgpu*width*height, 0.0);
          initializedoublearr(lower, blocknumgpu*width*height, 0.0);
          initializeintarr(crt, (blocknumgpu-1)*width*height, 0);
          initializeintarr(cr12, (blocknumgpu-2)*width*height, 0);
          initializeintarr(cr3, (blocknumgpu-2)*width*height, 0);
          initializeintarr(diffpos, (blocknumgpu-1)*width*height, 0);
          initializedoublearr(varblock0, blocknumgpu*width*height, 0.0);   
          initializedoublearr(varblock1, blocknumgpu*width*height, 0.0);   
          initializedoublearr(varblock2, blocknumgpu*width*height, 0.0);   
      }

      //------------------ calcacf2 ---------------------------
      prodnumarray = (int *)malloc(chanum*width*height*SIZEINT);
      indexarray = (int *)malloc(width*height*SIZEINT);
      Cblockvararray = (double *)malloc(sizeblockvararray);
      blocksdarray = (double *)malloc(sizeblockvararray);
      pnumgpu = (int *)malloc(width*height*SIZEINT);

      //--------------- N & B calculations in calcacf2 --------------------
      CNBmeanGPU = (double *)malloc( width * height * SIZEDOUBLE );
      CNBcovarianceGPU = (double *)malloc( width * height * SIZEDOUBLE );

      if (isNBcalculation) { 
          initializedoublearr(CNBmeanGPU, width * height, 0.0);
          initializedoublearr(CNBcovarianceGPU, width * height, 0.0);
      }

      // ------------------- perform calcacf3 calculation -------------------
      // Allocate memory on GPU for common arrays
      cudaMalloc((void **)&d_Cpixels, size);
      cudaMalloc((void **)&d_Cpixels1, size1);
      cudaMalloc((void **)&d_Cbleachcorr_params, w_temp * h_temp * bleachcorr_order * SIZEDOUBLE);
      cudaMalloc((void **)&d_Csamp, chanum * SIZEDOUBLE);
      cudaMalloc((void **)&d_Clag, chanum * SIZEINT);
      cudaMalloc((void **)&d_prod, size2);

      // Allocate memory on GPU for calcacf3
      if (!isNBcalculation) {
          cudaMalloc((void **)&d_prodnum, blocknumgpu * SIZEDOUBLE);
          cudaMalloc((void **)&d_blocksd, sizeblockvararray);
          cudaMalloc((void **)&d_upper, blocknumgpu*width*height*SIZEDOUBLE);
          cudaMalloc((void **)&d_lower, blocknumgpu*width*height*SIZEDOUBLE);
          cudaMalloc((void **)&d_crt, (blocknumgpu-1)*width*height*SIZEINT);
          cudaMalloc((void **)&d_cr12, (blocknumgpu-2)*width*height*SIZEINT);
          cudaMalloc((void **)&d_cr3, (blocknumgpu-2)*width*height*SIZEINT);
          cudaMalloc((void **)&d_diffpos, (blocknumgpu-1)*width*height*SIZEINT);
          cudaMalloc((void **)&d_varblock0, blocknumgpu*width*height*SIZEDOUBLE);
          cudaMalloc((void **)&d_varblock1, blocknumgpu*width*height*SIZEDOUBLE);
          cudaMalloc((void **)&d_varblock2, blocknumgpu*width*height*SIZEDOUBLE);
      }

      // copy to GPU for common arrays
      CUDA_CHECK_STATUS(cudaMemcpy(d_Cpixels, Cpixels, size, cudaMemcpyHostToDevice));
      CUDA_CHECK_STATUS(cudaMemcpy(d_Cpixels1, Cpixels1, size1, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_Cbleachcorr_params, Cbleachcorr_params, w_temp * h_temp * bleachcorr_order * SIZEDOUBLE, cudaMemcpyHostToDevice));
      CUDA_CHECK_STATUS(cudaMemcpy(d_Csamp, Csamp, chanum * SIZEDOUBLE, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_Clag, Clag, chanum * SIZEINT, cudaMemcpyHostToDevice)); 
      // NOTE: When pixels + prod are larger than half the memory available on CPU, this cudaMemcpy function will fail. The half memory limit is restriction in Java? If we comment this line, we will encounter a different error when the kernel calcacf3 is run.
      CUDA_CHECK_STATUS(cudaMemcpy(d_prod, prod, size2, cudaMemcpyHostToDevice)); 

      // copy to GPU for calcacf3      
      if (!isNBcalculation) {
          CUDA_CHECK_STATUS(cudaMemcpy(d_prodnum, prodnum, blocknumgpu*SIZEDOUBLE, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_blocksd, blocksd, sizeblockvararray, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_upper, upper, blocknumgpu*width*height*SIZEDOUBLE, cudaMemcpyHostToDevice));
          CUDA_CHECK_STATUS(cudaMemcpy(d_lower, lower, blocknumgpu*width*height*SIZEDOUBLE, cudaMemcpyHostToDevice));  
          CUDA_CHECK_STATUS(cudaMemcpy(d_crt, crt, (blocknumgpu-1)*width*height*SIZEINT, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_cr12, cr12, (blocknumgpu-2)*width*height*SIZEINT, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_cr3, cr3, (blocknumgpu-2)*width*height*SIZEINT, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_diffpos, diffpos, (blocknumgpu-1)*width*height*SIZEINT, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_varblock0, varblock0, blocknumgpu*width*height*SIZEDOUBLE, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_varblock1, varblock1, blocknumgpu*width*height*SIZEDOUBLE, cudaMemcpyHostToDevice)); 
          CUDA_CHECK_STATUS(cudaMemcpy(d_varblock2, varblock2, blocknumgpu*width*height*SIZEDOUBLE, cudaMemcpyHostToDevice)); 
      }

      //running kernel for calcacf3
      if (!isNBcalculation) {
          if(bleachcorr_gpu) {
              bleachcorrection<<<gridSize_Input, blockSize, 0, stream>>>(d_Cpixels, w_temp, h_temp, framediff, bleachcorr_order, frametime, d_Cbleachcorr_params);
              cudaDeviceSynchronize();
              CUDA_CHECK_STATUS(cudaGetLastError());
          }

          calcacf3<<<gridSize, blockSize, 0, stream>>>(d_Cpixels, cfXDistance, cfYDistance, blocklaggpu, width, height, w_temp, h_temp, pixbinX, pixbinY, framediff, correlatorp, correlatorq, chanum, frametime, d_Cpixels1, d_prod, d_prodnum, d_blocksd, d_upper, d_lower, d_crt, d_cr12, d_cr3, d_diffpos, d_varblock0, d_varblock1, d_varblock2, d_Csamp, d_Clag);

          cudaDeviceSynchronize();
          CUDA_CHECK_STATUS(cudaGetLastError());

          // copy memory from device to host for calcacf3
          CUDA_CHECK_STATUS(cudaMemcpy(Cpixels1, d_Cpixels1, size1, cudaMemcpyDeviceToHost));

          // CUDA release memory for calcacf3
          cudaFree(d_prodnum); cudaFree(d_blocksd); cudaFree(d_upper); cudaFree(d_lower); cudaFree(d_crt); cudaFree(d_cr12); cudaFree(d_cr3); cudaFree(d_diffpos); 
          cudaFree(d_varblock0); cudaFree(d_varblock1); cudaFree(d_varblock2); 

          // copy results in Cpixels1 to Cblocked1D. Values in Cpixels1 will be overwritten in calcacf2 calculation
          // memcpy(Cblocked1D, Cpixels1, size1);
          for (int i = 0; i < width * height * chanum; i++){
            *(Cblocked1D+i) = *(Cpixels1+i);
          }

          // initialize the values in indexarray and pnumgpu
          int counter_indexarray = 0;
          for (int y = 0; y < height; y++){
                for (int x = 0; x < width; x++){
                    *(indexarray + counter_indexarray) =  (int) *(Cpixels1 + (y*width + x));
          
                    //minimum number of products given the used correlator structure and the blockIndex ind
                    double tempval = *(indexarray + counter_indexarray) - log(sampchanumminus1) / log(2);
                    if (tempval < 0){
                        tempval = 0;
                    }
                    *(pnumgpu + counter_indexarray) = (int) floor(mtabchanumminus1/pow(2,tempval));
                    counter_indexarray = counter_indexarray + 1;
                }
          }

      } // if (!isNBcalculation)

      // ------------------- perform calcacf2 calculation -------------------

      // The pixel values in d_Cpixels has changed in calcacf3, reallocate original Cpixels to d_Cpixels
      CUDA_CHECK_STATUS(cudaMemcpy(d_Cpixels, Cpixels, size, cudaMemcpyHostToDevice));

      // Allocate memory on GPU for calcacf2
      cudaMalloc((void **)&d_prodnumarray, chanum*width*height*SIZEINT);
      cudaMalloc((void **)&d_indexarray, width*height*SIZEINT);
      cudaMalloc((void **)&d_Cblockvararray, sizeblockvararray);
      cudaMalloc((void **)&d_blocksdarray, sizeblockvararray);
      cudaMalloc((void **)&d_pnumgpu, width*height*SIZEINT);

      // Allocate memory on GPU for N & B calculations in calcacf2
      if (isNBcalculation) {
          cudaMalloc((void **)&d_NBmeanGPU, width * height * SIZEDOUBLE);
          cudaMalloc((void **)&d_NBcovarianceGPU, width * height * SIZEDOUBLE);
      }

      // copy to GPU for calcacf2
      CUDA_CHECK_STATUS(cudaMemcpy(d_indexarray, indexarray, width*height*SIZEINT, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_prodnumarray, prodnumarray, chanum*width*height*SIZEINT, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_Cblockvararray, Cblockvararray, sizeblockvararray, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_blocksdarray, blocksdarray, sizeblockvararray, cudaMemcpyHostToDevice)); 
      CUDA_CHECK_STATUS(cudaMemcpy(d_pnumgpu, pnumgpu, width*height*SIZEINT, cudaMemcpyHostToDevice)); 

      // copy to GPU for N & B calculatiosn in calcacf2
      if (isNBcalculation) {
          CUDA_CHECK_STATUS(cudaMemcpy(d_NBmeanGPU, CNBmeanGPU, width * height * SIZEDOUBLE, cudaMemcpyHostToDevice));
          CUDA_CHECK_STATUS(cudaMemcpy(d_NBcovarianceGPU, CNBcovarianceGPU, width * height * SIZEDOUBLE, cudaMemcpyHostToDevice));
      }

      //running kernel for calcacf2
      if(bleachcorr_gpu) {
          bleachcorrection<<<gridSize_Input, blockSize, 0, stream>>>(d_Cpixels, w_temp, h_temp, framediff, bleachcorr_order, frametime, d_Cbleachcorr_params);
          cudaDeviceSynchronize();
          CUDA_CHECK_STATUS(cudaGetLastError());
      }

      int numbin = framediff;		// number Of data points When they are binned
      int currentIncrement = 1;
      int ctbin = 0;
      bool runthis;

      int calcacf2_x_start = (isNBcalculation) ? 1 : 0;
      int calcacf2_x_end = (isNBcalculation) ? 2 : chanum;

      for (int x = calcacf2_x_start; x < calcacf2_x_end; x++) {
          runthis = false;
	  if (currentIncrement != Csamp[x]){ // check whether the channel width has changed
              // Set the currentIncrement accordingly
              numbin = (int)floor((double)numbin / 2.0);
              currentIncrement = (int)Csamp[x];
              ctbin++;
              runthis = true;
          }

          if (runthis){ // check whether the channel width has changed
              calcacf2a<<<gridSize_Input, blockSize, 0, stream>>>(d_Cpixels, w_temp, h_temp, numbin);
              cudaDeviceSynchronize();
              CUDA_CHECK_STATUS(cudaGetLastError());
          }

          calcacf2b<<<gridSize, blockSize, 0, stream>>>(d_Cpixels, cfXDistance, cfYDistance, width, height, w_temp, h_temp, pixbinX, pixbinY, d_Cpixels1, d_prod, d_Clag, d_prodnumarray, d_indexarray, d_Cblockvararray, d_blocksdarray, d_pnumgpu, x, numbin, currentIncrement, ctbin, isNBcalculation, d_NBmeanGPU, d_NBcovarianceGPU);
          cudaDeviceSynchronize();
          CUDA_CHECK_STATUS(cudaGetLastError());
      }

      // copy memory from device to host for calcacf2
      if (isNBcalculation) { 
          CUDA_CHECK_STATUS(cudaMemcpy(CNBmeanGPU, d_NBmeanGPU, width * height * SIZEDOUBLE, cudaMemcpyDeviceToHost));
          CUDA_CHECK_STATUS(cudaMemcpy(CNBcovarianceGPU, d_NBcovarianceGPU, width * height * SIZEDOUBLE, cudaMemcpyDeviceToHost));
      } else {
          CUDA_CHECK_STATUS(cudaMemcpy(Cpixels1, d_Cpixels1, size1, cudaMemcpyDeviceToHost));
          CUDA_CHECK_STATUS(cudaMemcpy(Cblockvararray, d_Cblockvararray, sizeblockvararray, cudaMemcpyDeviceToHost));
      }

      // CUDA release memory for calcacf2
      cudaFree(d_prodnumarray); cudaFree(d_indexarray); cudaFree(d_Cblockvararray); cudaFree(d_blocksdarray); cudaFree(d_pnumgpu);

      // CUDA release memory for N & B calculations in calcacf2
      if (isNBcalculation) { cudaFree(d_NBmeanGPU); cudaFree(d_NBcovarianceGPU); }

      // CUDA release memory for all
      cudaFree(d_Cpixels); cudaFree(d_Cpixels1); cudaFree(d_Cbleachcorr_params); cudaFree(d_Csamp); cudaFree(d_Clag); cudaFree(d_prod);

      CUDA_CHECK_STATUS(cudaStreamDestroy( stream ));

      // Reference: https://github.com/zchee/cuda-sample/blob/master/0_Simple/simpleMultiGPU/simpleMultiGPU.cu 
      // cudaDeviceReset causes the driver to clean up all state. While
      // not mandatory in normal operation, it is good practice.  It is also
      // needed to ensure correct operation when the application is being
      // profiled. Calling cudaDeviceReset causes all profile data to be
      // flushed before the application exits
      cudaDeviceReset(); 

      // copy values to Java output arrays.
      env->SetDoubleArrayRegion(pixels1, 0 , width * height * chanum, Cpixels1);  
      env->SetDoubleArrayRegion(blockvararray, 0 , chanum * width * height, Cblockvararray);
      env->SetDoubleArrayRegion(blocked1D, 0 , width * height * chanum, Cblocked1D);  

      if (isNBcalculation) {
          env->SetDoubleArrayRegion(NBmeanGPU, 0 , width * height , CNBmeanGPU); 
          env->SetDoubleArrayRegion(NBcovarianceGPU, 0 , width * height , CNBcovarianceGPU); 
      }

      // free all pointers
      free(prod); free(prodnum); free(blocksd); free(upper); free(lower); free(crt); free(cr12); free(cr3); free(diffpos); free(varblock0); free(varblock1); free(varblock2); 
      free(indexarray); free(prodnumarray); free(Cblockvararray); free(blocksdarray); free(pnumgpu);
      free(Cpixels1); free(Cblocked1D);
      free(CNBmeanGPU); free(CNBcovarianceGPU);
      
      // release resources
      env->ReleaseFloatArrayElements(pixels, Cpixels, 0);
      env->ReleaseDoubleArrayElements(bleachcorr_params, Cbleachcorr_params, 0);
      env->ReleaseDoubleArrayElements(samp, Csamp, 0);
      env->ReleaseIntArrayElements(lag, Clag, 0);

      return;  

    } catch (std::runtime_error & e)  {
      // see: https://www.rgagnon.com/javadetails/java-0323.html
      jclass Exception = env->FindClass("java/lang/Exception");
      env->ThrowNew(Exception, e.what());

      // CUDA release memory for calcacf3
      if (!isNBcalculation) {
          cudaFree(d_prodnum); cudaFree(d_blocksd); cudaFree(d_upper); cudaFree(d_lower); cudaFree(d_crt); cudaFree(d_cr12); cudaFree(d_cr3); cudaFree(d_diffpos); 
          cudaFree(d_varblock0); cudaFree(d_varblock1); cudaFree(d_varblock2); 
      }

      // CUDA release memory for calcacf2
      cudaFree(d_prodnumarray); cudaFree(d_indexarray); cudaFree(d_Cblockvararray); cudaFree(d_blocksdarray); cudaFree(d_pnumgpu);

      // CUDA release memory for N & B calculations in calcacf2
      if (isNBcalculation) { cudaFree(d_NBmeanGPU); cudaFree(d_NBcovarianceGPU); }

      // CUDA release memory for all
      cudaFree(d_Cpixels); cudaFree(d_Cpixels1); cudaFree(d_Cbleachcorr_params); cudaFree(d_Csamp); cudaFree(d_Clag); cudaFree(d_prod); 

      cudaDeviceReset(); 

      // free all pointers
      free(prod); free(prodnum); free(blocksd); free(upper); free(lower); free(crt); free(cr12); free(cr3); free(diffpos); free(varblock0); free(varblock1); free(varblock2); 
      free(indexarray); free(prodnumarray); free(Cblockvararray); free(blocksdarray); free(pnumgpu);
      free(Cpixels1); free(Cblocked1D);
      free(CNBmeanGPU); free(CNBcovarianceGPU);

      // release resources
      env->ReleaseFloatArrayElements(pixels, Cpixels, 0);
      env->ReleaseDoubleArrayElements(bleachcorr_params, Cbleachcorr_params, 0);
      env->ReleaseDoubleArrayElements(samp, Csamp, 0);
      env->ReleaseIntArrayElements(lag, Clag, 0);

      return;  
    }
}

/* ------------------------------------------
AUTOCORRELATION SINGLE DIMENSION ARRAY CALCULATION END
------------------------------------------ */

/* -------------------------------------------------------------------------------------------------------
* from com_github_gpufit_Gpufit.cpp END
------------------------------------------------------------------------------------------------------- */
