/*
 * Chain_loc v 1.0 is the source code for a mex file that will input image data and parameter estimators and output localized data
 * Calling this function in matlab will look like
 * [xf_all, yf_all, N, off_all, sigx, sigy, xf_crlb, yf_crlb, N_crlb, off_crlb, sigx_crlb, sigy_crlb, llv, framenum_all] = chain_loc( i1, numthreads, angle(in rads))
 * written by Andrew Nelson Version v 1.1 on 5/2/15
 */

/*
Version 1.2 has had substantial debugging and careful comparison between cpu version. Removed extra diagnostic code that has been commented out and has no use in regular code
	At this point we can be sured that the algorithm for calculating position, total area, and offset are working properly on the GPU
	This version is usable with the localiztion code Quhzx_01_3.m

Fixed
	Error codes should not be given off when inputs and outputs don't match expected arguments
	Fixed a problem considerably prolonging computation time by removing gpu device invokation early in main loop ( i.e. reset the gpu, synchronize threads, initialize gpu ....)
	Added redundant void __global__ loops to handle multiple sizes of input images
*/

#include <mex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#define PI 3.14159265358979323846



/*
 * Device code
 * 
 * To facilitate coding (for me) I have copied the localization algorithm to be used with multiple sized areas
 */

/*
 
 Device Functions

*/
__device__ double device_det(double Fisher[36])
{
	double det;
	det = Fisher[0] * (Fisher[7] * (Fisher[14] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) + Fisher[26] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[32] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]))) - Fisher[13] * (Fisher[8] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]))) + Fisher[19] * (Fisher[8] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[25] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) + Fisher[31] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])))) - Fisher[6] * (Fisher[1] * (Fisher[14] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) + Fisher[26] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[32] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]))) - Fisher[13] * (Fisher[2] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]))) + Fisher[19] * (Fisher[2] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])))) + Fisher[12] * (Fisher[1] * (Fisher[8] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]))) + Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) - Fisher[18] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) + Fisher[24] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) - Fisher[30] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))));
	return det;
}


/*

Global Functions

*/

// localize
__global__ void  localize(double *d_iall,
	double *d_xf_all,
	double *d_yf_all,
	double *d_N,
	double *d_off,
	double *d_sigx,
	double *d_sigy,
	double *d_xf_crlb,
	double *d_yf_crlb,
	double *d_N_crlb,
	double *d_off_crlb,
	double *d_sigx_crlb,
	double *d_sigy_crlb,
	double *d_llv,
	double ang,
	int lpcnt,
	int numi)
{

	// Declare variables
	const int pix = 23; // pixel size
	__shared__ double xgrid[pix*pix];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[pix*pix];			// this will reduce calls to global device memory
	double dudx, dudy, dudsx, dudsy, d2udx2, d2udy2, d2udsx2, d2udsy2, dudn, dudo, Ex, Ey, u;
	double d_x, d_y, d_n, d_sx, d_sy, d_o, dd_x, dd_y, dd_sx, dd_sy, dd_n, dd_o, x, y;
	// these variables will exist on the register of each thread
	double d_beta1[6] = {0, 0, 0, 0, 0, 0};
	int tx = threadIdx.x;
	int index = blockIdx.x*blockDim.x + tx;		// calculate thread index
	double d_i2[pix*pix];				// initialize data for image
	double llv;
	double fisher[36] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double det_fish = 0.0;
			// create xgrid and ygrid we want to create the grid regardless of whether the index is crunching on an image
		if (tx == 0){
			for (int i = 0; i <pix; i++){
				for(int j = 0; j <pix; j++){
					x = (double)j - ((double)pix-1)/2;
					y = (double)i - ((double)pix-1)/2;
					xgrid[j*pix + i] = x*cos(ang) - y*sin(ang);
					ygrid[j*pix + i] = x*sin(ang) + y*cos(ang);
				}
			}
		}
	
	if (index < numi){			// check to see that threads only work if an image exists
		// buffer all the variables into shared memory and registers and build guesses

		d_beta1[0] = 0.0;
		d_beta1[1] = 0.0;
		d_beta1[2] = 0.0;
		d_beta1[3] = 1.5;						 // guess on sigma widths
		d_beta1[4] = 1.5;						 // guess on sigma widths
		d_beta1[5] = 10000;
		for (int i = 0; i <pix*pix; i++) {
			d_i2[i] = d_iall[i + index*pix*pix];	         // this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
			d_beta1[0] +=xgrid[i]*d_i2[i];						 // sum of x and image weight
			d_beta1[1] +=ygrid[i]*d_i2[i];						 // sum of y and image weight
			d_beta1[2] +=d_i2[i];									 // image sum
			if (d_beta1[5] > d_i2[i]){d_beta1[5] = d_i2[i];}				 // find minimum of image
		}
		d_beta1[0] = d_beta1[0] / d_beta1[2];
		d_beta1[1] = d_beta1[1] / d_beta1[2];
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < lpcnt; counttry++){
			d_x =	0.0;
			d_y =	0.0;
			d_n =	0.0;
			d_sx =	0.0;
			d_sy =	0.0;
			d_o =	0.0;
			dd_x =	0.0; 	//wipe incremental variables each loop to give correct correction factor
			dd_y =	0.0;
			dd_n =	0.0;
			dd_sx = 0.0;
			dd_sy = 0.0;
			dd_o =	0.0;
			u =		0;
			Ey		=	0;
			Ex		=	0;
			llv		=	0.0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < pix; rowcount++){	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < pix; colcount++){	// FOR 3 loops over all columns

					// x/ygrid is col major(come from matlab) and i3 is col major 
					// these three lines help define the fitting gaussian as deined by the current iteration of parameters
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5) / sqrt(2.0 * d_beta1[4] * d_beta1[4])) - erf((ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5) / sqrt(2.0 * d_beta1[4] * d_beta1[4])));
					u = d_beta1[2] * Ex*Ey + d_beta1[5];

					// first derivatives calculations
						// these are done pixel by pixel with the sum added up in the d_x and dd_x areas
					dudx = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudy = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[4] * d_beta1[4]))*(exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))*Ex;
					dudsx = (d_beta1[2] / (sqrt(2.0*PI) * powf(d_beta1[3], 2.0)))*((xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5) * exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudsy = (d_beta1[2] / (sqrt(2.0*PI) * powf(d_beta1[4], 2.0)))*((ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5) * exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))*Ex;
					dudn = Ex*Ey;
					dudo = 1.0;

					// second derivatives
						// these are calcualted in a  similar manner to the first derivatives
					d2udx2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ey;
					d2udy2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[4], 3.0))*((ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4]))))*Ex;
					d2udsx2 = (Ey*d_beta1[2] / (sqrt(2.0 * PI)))
						*(powf(d_beta1[3], -5.0)*(powf((xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- powf((xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))
						- 2 * powf(d_beta1[3], -3.0)*((xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5)     *exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5) *exp(-powf(xgrid[rowcount + colcount*pix] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))));
					d2udsy2 = (Ex*d_beta1[2] / (sqrt(2.0 * PI)))
						*(powf(d_beta1[4], -5.0)*(powf((ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- powf((ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))
						- 2 * powf(d_beta1[4], -3.0)*((ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5)     *exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5) *exp(-powf(ygrid[rowcount + colcount*pix] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4]))));

					// summing variable to lead to correction factors
						// these variables keep track of the correction which is given by summing over the entire pixel
					d_x = d_x + dudx*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					d_sx = d_sx + dudsx*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_sx = dd_sx + d2udsx2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudsx, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_sy = d_sy + dudsy*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_sy = dd_sy + d2udsy2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudsy, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*pix] / powf(u, 2.0);

					
					if (counttry == lpcnt-1){  // on the last count, construct fisher information matrix elements
						fisher[0] += dudx*dudx / u;
						fisher[1] += dudx*dudy / u;
						fisher[2] += dudx*dudn / u;
						fisher[3] += dudx*dudo / u;
						fisher[4] += dudx*dudsx / u;
						fisher[5] += dudx*dudsy / u;

						fisher[6] += dudy*dudx / u;
						fisher[7] += dudy*dudy / u;
						fisher[8] += dudy*dudn / u;
						fisher[9] += dudy*dudo / u;
						fisher[10] += dudy*dudsx / u;;
						fisher[11] += dudy*dudsy / u;;

						fisher[12] += dudn*dudx / u;
						fisher[13] += dudn*dudy / u;
						fisher[14] += dudn*dudn / u;
						fisher[15] += dudn*dudo / u;
						fisher[16] += dudn*dudsx / u;
						fisher[17] += dudn*dudsy / u;

						fisher[18] += dudo*dudx / u;
						fisher[19] += dudo*dudy / u;
						fisher[20] += dudo*dudn / u;
						fisher[21] += dudo*dudo / u;
						fisher[22] += dudo*dudsx / u;
						fisher[23] += dudo*dudsy / u;

						fisher[24] += dudsx*dudx / u;
						fisher[25] += dudsx*dudy / u;
						fisher[26] += dudsx*dudn / u;
						fisher[27] += dudsx*dudo / u;
						fisher[28] += dudsx*dudsx / u;
						fisher[29] += dudsx*dudsy / u;

						fisher[30] += dudsy*dudx / u;
						fisher[31] += dudsy*dudy / u;
						fisher[32] += dudsy*dudn / u;
						fisher[33] += dudsy*dudo / u;
						fisher[34] += dudsy*dudsx / u;
						fisher[35] += dudsy*dudsy / u;
						llv +=  d_i2[rowcount+colcount*pix] *  log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*pix]*log(d_i2[rowcount + colcount*pix] + 0.0000000000000001) + d_i2[rowcount + colcount*pix];
					}
				} // END FOR 3
			} // END FOR2
			// correct beta1 values with tolerances
			d_beta1[0] = d_beta1[0] - d_x / dd_x;
			d_beta1[1] = d_beta1[1] - d_y / dd_y;
			d_beta1[2] = d_beta1[2] - d_n / dd_n;
			d_beta1[3] = d_beta1[3] - d_sx / dd_sx;
			d_beta1[4] = d_beta1[4] - d_sy / dd_sy;
			d_beta1[5] = d_beta1[5] - d_o / dd_o;

		}   // end FOR 1

		if (d_beta1[0] == d_beta1[0] && d_beta1[1] == d_beta1[1] && d_beta1[2] == d_beta1[2] && d_beta1[5] == d_beta1[5] && d_beta1[3] == d_beta1[3] && d_beta1[4] == d_beta1[4] && d_beta1[5] == d_beta1[5]){ // begin is numeric if statement
			if (d_beta1[2] > 0 && d_beta1[0] >= -(double)pix/2 && d_beta1[0] <= (double)pix/2 && d_beta1[1] <= (double)pix/2 && d_beta1[1] >= -(double)pix/2 ){ // was the molecule inside the grid? Was N positive? if yes then record the point
				d_xf_all[index] = d_beta1[0] ; 
				d_yf_all[index] = d_beta1[1] ;
				     d_N[index] = d_beta1[2];
				  d_sigx[index] = d_beta1[3];
				  d_sigy[index] = d_beta1[4];
				   d_off[index] = d_beta1[5];
				   d_llv[index] = llv;
				/*d_xf_all[index] = testi ; */
				
					/*d_N[index] = testi;
				  d_sigx[index] = testi;
				  d_sigy[index] = testi;
				   d_off[index] =testi;
				   d_llv[index] =testi;*/

				// calculate crlb's for estimators	
				
				// UPDATE FOR SIGMA VALUES
				det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
				d_xf_crlb[index] = (fisher[7] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[13] * (fisher[8] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) + fisher[19] * (fisher[8] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[25] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) + fisher[31] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[26] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])))) / det_fish;
				d_yf_crlb[index] = -(-(fisher[0] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[12] * (fisher[2] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[2] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) - fisher[24] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[30] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[26] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])))) / det_fish);
				d_N_crlb[index] = (fisher[0] * (fisher[7] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[31] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[1] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish;
				d_off_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish);
				d_sigx_crlb[index] = (fisher[0]*(fisher[7]*(fisher[14]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])+fisher[32]*(fisher[15]*fisher[23]-fisher[21]*fisher[17]))-fisher[13]*(fisher[8]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[23]-fisher[21]*fisher[11]))+fisher[19]*(fisher[8]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[31]*(fisher[8]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])+fisher[20]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])))-fisher[6]*(fisher[1]*(fisher[14]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])+fisher[32]*(fisher[15]*fisher[23]-fisher[21]*fisher[17]))-fisher[13]*(fisher[2]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[23]-fisher[21]*fisher[5]))+fisher[19]*(fisher[2]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])))+fisher[12]*(fisher[1]*(fisher[8]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[23]-fisher[21]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[23]-fisher[21]*fisher[5]))+fisher[19]*(fisher[2]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])-fisher[8]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])-fisher[8]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[11]-fisher[9]*fisher[5])))-fisher[18]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])-fisher[8]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])-fisher[8]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])+fisher[14]*(fisher[3]*fisher[11]-fisher[9]*fisher[5])))+fisher[30]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])+fisher[20]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])-fisher[8]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[19]*(fisher[2]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])-fisher[8]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])+fisher[14]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))))/det_fish;
				d_sigy_crlb[index] = -(-(fisher[0]*(fisher[7]*(fisher[14]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])+fisher[26]*(fisher[15]*fisher[22]-fisher[21]*fisher[16]))-fisher[13]*(fisher[8]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[22]-fisher[21]*fisher[10]))+fisher[19]*(fisher[8]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[25]*(fisher[8]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])+fisher[20]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])))-fisher[6]*(fisher[1]*(fisher[14]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])+fisher[26]*(fisher[15]*fisher[22]-fisher[21]*fisher[16]))-fisher[13]*(fisher[2]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[22]-fisher[21]*fisher[4]))+fisher[19]*(fisher[2]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])))+fisher[12]*(fisher[1]*(fisher[8]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[22]-fisher[21]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[22]-fisher[21]*fisher[4]))+fisher[19]*(fisher[2]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])-fisher[8]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])-fisher[8]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[10]-fisher[9]*fisher[4])))-fisher[18]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])-fisher[8]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])-fisher[8]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])+fisher[14]*(fisher[3]*fisher[10]-fisher[9]*fisher[4])))+fisher[24]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])+fisher[20]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])-fisher[8]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[19]*(fisher[2]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])-fisher[8]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])+fisher[14]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))))/det_fish);
				/*d_xf_crlb[index] = testi;
				d_yf_crlb[index] = testi;
				d_N_crlb[index] = testi;
				d_off_crlb[index] = testi;
				d_sigx_crlb[index] = testi;
				d_sigy_crlb[index] = testi;*/
			}
			else{						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_sigx[index] = -1;
				d_sigy[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_sigx_crlb[index] = -1;
				d_sigy_crlb[index] = -1;
				d_llv[index] = llv;
			} 
		} //end is numeric if statement
		else{
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_sigx[index] = -1;
			d_sigy[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_sigx_crlb[index] = -1;
			d_sigy_crlb[index] = -1;
			d_llv[index] = llv;
		} // end else fail statement
	} 
}   // end localize 7

/*
 * Host code
 *
 * 
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{	
	/* Declare all variables.*/
	double *iall;				// the pointer to the array of all images to be analyzed
	double *d_iall;			// Pointer to image array on gpu
	double angle;
	double *d_xf_all;
	double *d_yf_all;
	double *d_sigx_all;
	double *d_sigy_all;
	double *d_N;
	double *d_off;
	double *d_llv;
	double *d_xf_crlb;
	double *d_yf_crlb;
	double *d_N_crlb;
	double *d_off_crlb;
	double *d_sigx_crlb;
	double * d_sigy_crlb;
	double *xf, *xfc, *yf, *yfc, *n, *nc, *sigx, *sigxc, *sigy, *sigyc, *off, *offc, *llv;
	size_t threadsperblock;


	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int numi;				// number of images imported
	int lpcnt;				// number of loops to perform the MLE calculation over
	const size_t *idims;
	

	/* Throw an error if the input does not match expectations. */
	if (nrhs != 4) {
		printf("Must have 4 inputs ( i1, numthreads, angle(in rads), MLE loop count)\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
		printf("i1 must be a nxm double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
		printf("angle must be a DOUBLE\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	irow = (int)idims[0];
	numi = (int)idims[1];
	angle = (double)mxGetScalar(prhs[2]);
	lpcnt = (int)mxGetScalar(prhs[3]);
	if (numi > 1000000 || numi < 1){
		numi = 1;
	}
	int imem = irow*numi*sizeof(double);
	int vmem = numi*sizeof(double);

	// verify that the input variables are what was expected
	if (irow != 529){
		printf("Images are of incorrect size. There must be a square number of rows in the entry.\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (nlhs != 13){
		printf("You must have 13 output variables [xf_all, xf_crlb, yf_all, yf_crlb, N, N_crlb, sigx_all, sigx_crlb, sigy_all, sigy_crlb, off_all, off_crlb, llv_all]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	
	// allocate memory and copy it onto the gpu device
	// iall
	checkCudaErrors(cudaMalloc((void**)&d_iall, imem)); // allocate image memory
	checkCudaErrors(cudaMemcpy(d_iall, iall, imem, cudaMemcpyHostToDevice)); // copy images from device to host


	// allocate memory for fitted variables that will be returned from device
	checkCudaErrors(cudaMalloc((void**)&d_xf_all   , vmem)); // allocate xf_all memory
	checkCudaErrors(cudaMalloc((void**)&d_xf_crlb  , vmem)); // allocate xf_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_yf_all   , vmem)); // allocate yf_all memory
	checkCudaErrors(cudaMalloc((void**)&d_yf_crlb  , vmem)); // allocate yf_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_sigx_all , vmem)); // allocate sigx memory
	checkCudaErrors(cudaMalloc((void**)&d_sigx_crlb, vmem)); // allocate sigx_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_sigy_all , vmem)); // allocate sigy memory
	checkCudaErrors(cudaMalloc((void**)&d_sigy_crlb, vmem)); // allocate sigy_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_N		   , vmem)); // allocate N memory
	checkCudaErrors(cudaMalloc((void**)&d_N_crlb   , vmem)); // allocate N_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_off	   , vmem)); // allocate off memory
	checkCudaErrors(cudaMalloc((void**)&d_off_crlb , vmem)); // allocate N memory
	checkCudaErrors(cudaMalloc((void**)&d_llv	   , vmem)); // allocate llv memory

	/* Run GPU kernel*/
	threadsperblock = mxGetScalar(prhs[1]);  // get number of threads perblock from matlab

	localize << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_xf_all, d_yf_all, d_N, d_off, d_sigx_all, d_sigy_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, d_sigx_crlb, d_sigy_crlb, d_llv, angle, lpcnt, numi);
		


	// Allocate host side memory for output arrays at the output pointer positions
	plhs[0] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[3] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[4] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[5] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[6] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[7] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[8] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
	plhs[9] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
   plhs[10] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
   plhs[11] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);
   plhs[12] = mxCreateNumericMatrix(numi,1,mxDOUBLE_CLASS, mxREAL);


	   xf = (double *)mxGetPr(plhs[0]);
	  xfc = (double *)mxGetPr(plhs[1]);
	   yf = (double *)mxGetPr(plhs[2]);
	  yfc = (double *)mxGetPr(plhs[3]);
	    n = (double *)mxGetPr(plhs[4]);
	   nc = (double *)mxGetPr(plhs[5]);
	 sigx = (double *)mxGetPr(plhs[6]);
	sigxc = (double *)mxGetPr(plhs[7]);
	 sigy = (double *)mxGetPr(plhs[8]);
	sigyc = (double *)mxGetPr(plhs[9]);
	  off = (double *)mxGetPr(plhs[10]);
	 offc = (double *)mxGetPr(plhs[11]);
	  llv = (double *)mxGetPr(plhs[12]);
	// copy memory from device to host
	checkCudaErrors(cudaMemcpy(xf    , d_xf_all    ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(xfc   , d_xf_crlb   ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(yf    , d_yf_all    ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(yfc   , d_yf_crlb   ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(n     , d_N         ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(nc    , d_N_crlb    ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(sigx  , d_sigx_all  ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(sigxc , d_sigx_crlb ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(sigy  , d_sigy_all  ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(sigyc , d_sigy_crlb ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(off   , d_off       ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(offc  , d_off_crlb  ,vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(llv   , d_llv       ,vmem, cudaMemcpyDeviceToHost));

	// clean up
	cudaFree(d_iall);
	cudaFree(d_N);  
	cudaFree(d_xf_all); 
	cudaFree(d_yf_all); 
	cudaFree(d_off);
	cudaFree(d_sigx_all);
	cudaFree(d_sigy_all);
	cudaFree(d_xf_crlb);
	cudaFree(d_yf_crlb);
	cudaFree(d_N_crlb);
	cudaFree(d_off_crlb);
	cudaFree(d_sigx_crlb);
	cudaFree(d_sigy_crlb);
	cudaFree(d_llv);
	
}  // DONE
