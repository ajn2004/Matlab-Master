/*
 * Chain_loc v 1.0 is the source code for a mex file that will input image data and parameter estimators and output localized data
 * Calling this function in matlab will look like
 * [xf_all, yf_all, N, off_all, sigx, sigy, xf_crlb, yf_crlb, N_crlb, off_crlb, sigx_crlb, sigy_crlb, llv, framenum_all] = chain_loc(i1, beta0, framenum_temp, xpeaks, ypeaks,xpix, ypix, numthreads)
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

Still possibly broken
	CRLB generating equations have not been double / triple checked yet these values are not yet reliable
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
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

// localize 7
__global__ void  localize7(double *d_iall,
	double *d_beta0,
	double *d_x_cm,
	double *d_y_cm,
	double *d_framenum_temp,
	double *d_xf_all,
	double *d_yf_all,
	double *d_N,
	double *d_off,
	double *d_sigx,
	double *d_sigy,
	double *d_framenum_all,
	double *d_xf_crlb,
	double *d_yf_crlb,
	double *d_N_crlb,
	double *d_off_crlb,
	double *d_sigx_crlb,
	double *d_sigy_crlb,
	int irow,
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi)
{
	// Declare variables
	
	__shared__ double xgrid[7*7];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[7*7];			// this will reduce calls to global device memory

	// these variables will exist on the register of each thread
	double dudx, dudy, dudsx, dudsy, d2udx2, d2udy2, d2udsx2, d2udsy2, dudn, dudo, Ex, Ey, u, llv;
	int index = blockIdx.x*blockDim.x + threadIdx.x;		// calculate thread index
	double d_i2[7*7];				// initialize data for image
	double d_beta1[6];				// initialize data for beta1
	double d_x, d_y, d_n, d_sx, d_sy, d_o, dd_x, dd_y, dd_sx, dd_sy, dd_n, dd_o;
	double fisher[36] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double det_fish;
	
	if (index < numi){			// check to see that threads only work if an image exists
		// buffer all the variables into shared memory and registers
		for (int i = 0; i <7*7; i++) {
			xgrid[i] = d_xpix[i];
			ygrid[i] = d_ypix[i];
			d_i2[i] = d_iall[i + index * 7*7];		// this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
		}
		for (int j = 0; j<6; j++){
			d_beta1[j] = d_beta0[j + index * 6];	// similar buffering strategy employed for d_beta1
		}											// we are assuming the information comes in as (x_guess, y_guess, n_guess, sigx_guess, sigy_guess, off_guess) which gives a beta vector of 6 values

		
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < 50; counttry++){
			d_x = 0.0;
			d_y = 0.0;
			d_n = 0.0;
			d_sx = 0.0;
			d_sy = 0.0;
			d_o = 0.0;
			dd_x = 0.0; 	//wipe incremental variables each loop to give correct correction factor
			dd_y = 0.0;
			dd_n = 0.0;
			dd_sx = 0.0;
			dd_sy = 0.0;
			dd_o = 0.0;
			u = 0;
			Ey = 0;
			Ex = 0;
			llv = 0.0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < irow; rowcount++){	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < irow; colcount++){	// FOR 3 loops over all columns

					// x/ygrid is col major(come from matlab) and i3 is col major 
					// these three lines help define the fitting gaussian as deined by the current iteration of parameters
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5) / sqrt(2.0 * d_beta1[4] * d_beta1[4])) - erf((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5) / sqrt(2.0 * d_beta1[4] * d_beta1[4])));
					u = d_beta1[2] * Ex*Ey + d_beta1[5];

					// first derivatives calculations
						// these are done pixel by pixel with the sum added up in the d_x and dd_x areas
					dudx = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudy = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[4] * d_beta1[4]))*(exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))*Ex;
					dudsx = (d_beta1[2] / (sqrt(2.0*PI) * powf(d_beta1[3], 2.0)))* ((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5) * exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudsy = (d_beta1[2] / (sqrt(2.0*PI) * powf(d_beta1[4], 2.0)))* ((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5) * exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))*Ex;
					dudn = Ex*Ey;
					dudo = 1.0;

					// second derivatives
						// these are calcualted in a  similar manner to the first derivatives
					d2udx2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ey;
					d2udy2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[4], 3.0))*((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4]))))*Ex;
					d2udsx2 = (Ey*d_beta1[2] / (sqrt(2.0 * PI)))
						*(powf(d_beta1[3], -5.0)*(powf((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- powf((xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))
						- 2 * powf(d_beta1[3], -3.0)*((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5)     *exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) 
						- (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5) *exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))));
					d2udsy2 = (Ex*d_beta1[2] / (sqrt(2.0 * PI)))
						*(powf(d_beta1[4], -5.0)*(powf((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- powf((ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])))
						- 2 * powf(d_beta1[4], -3.0)*((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5)     *exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4])) 
						- (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5) *exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[4] * d_beta1[4]))));

					// summing variable to lead to correction factors
						// these variables keep track of the correction which is given by summing over the entire pixel
					d_x = d_x + dudx*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					d_sx = d_sx + dudsx*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_sx = dd_sx + d2udsx2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudsx, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_sy = d_sy + dudsy*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_sy = dd_sy + d2udsy2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudsy, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*irow] / powf(u, 2.0);

					
					if (counttry == 49){  // on the last count, construct fisher information matrix elements
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
						llv +=  d_i2[rowcount+colcount*irow] *  log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*irow]*log(d_i2[rowcount + colcount*irow] + 0.0000000000000001) + d_i2[rowcount + colcount*irow];
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
			if (d_beta1[2] > 0 && d_beta1[0] > xgrid[0] && d_beta1[0] < xgrid[irow*irow - 1] && d_beta1[1] < ygrid[irow*irow - 1] && d_beta1[1] > ygrid[0] && d_beta1[3] > 0 && d_beta1[3] < 100 && d_beta1[4] < 100 && d_beta1[4] > 0 && d_beta1[5] > 0){ // was the molecule inside the grid? Was N positive? if yes then record the point
				d_xf_all[index] = d_beta1[0] + d_x_cm[index]; // correct position for x
				d_yf_all[index] = d_beta1[1] + d_y_cm[index]; // correct position for y
				d_N[index] = d_beta1[2];
				d_sigx[index] = d_beta1[3];
				d_sigy[index] = d_beta1[4];
				d_off[index] = d_beta1[5];
				d_framenum_all[index] = d_framenum_temp[index];
				d_llv[index] = llv;

				// calculate crlb's for estimators	
				
				// UPDATE FOR SIGMA VALUES
				det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
				d_xf_crlb[index] = (fisher[7] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[13] * (fisher[8] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) + fisher[19] * (fisher[8] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[25] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) + fisher[31] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[26] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])))) / det_fish;
				d_yf_crlb[index] = -(-(fisher[0] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[12] * (fisher[2] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[2] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) - fisher[24] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[30] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[26] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])))) / det_fish);
				d_N_crlb[index] = (fisher[0] * (fisher[7] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[31] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[1] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish;
				d_off_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish);
				d_sigx_crlb[index] = (fisher[0]*(fisher[7]*(fisher[14]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])+fisher[32]*(fisher[15]*fisher[23]-fisher[21]*fisher[17]))-fisher[13]*(fisher[8]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[23]-fisher[21]*fisher[11]))+fisher[19]*(fisher[8]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[31]*(fisher[8]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])+fisher[20]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])))-fisher[6]*(fisher[1]*(fisher[14]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])+fisher[32]*(fisher[15]*fisher[23]-fisher[21]*fisher[17]))-fisher[13]*(fisher[2]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[23]-fisher[21]*fisher[5]))+fisher[19]*(fisher[2]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])))+fisher[12]*(fisher[1]*(fisher[8]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[23]-fisher[21]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[21]*fisher[35]-fisher[33]*fisher[23])-fisher[20]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[23]-fisher[21]*fisher[5]))+fisher[19]*(fisher[2]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])-fisher[8]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])-fisher[8]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[11]-fisher[9]*fisher[5])))-fisher[18]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])+fisher[32]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[35]-fisher[33]*fisher[17])-fisher[14]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[35]-fisher[33]*fisher[11])-fisher[8]*(fisher[3]*fisher[35]-fisher[33]*fisher[5])+fisher[32]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[31]*(fisher[2]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])-fisher[8]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])+fisher[14]*(fisher[3]*fisher[11]-fisher[9]*fisher[5])))+fisher[30]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])+fisher[20]*(fisher[9]*fisher[17]-fisher[15]*fisher[11]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[23]-fisher[21]*fisher[17])-fisher[14]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[17]-fisher[15]*fisher[5]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[23]-fisher[21]*fisher[11])-fisher[8]*(fisher[3]*fisher[23]-fisher[21]*fisher[5])+fisher[20]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))-fisher[19]*(fisher[2]*(fisher[9]*fisher[17]-fisher[15]*fisher[11])-fisher[8]*(fisher[3]*fisher[17]-fisher[15]*fisher[5])+fisher[14]*(fisher[3]*fisher[11]-fisher[9]*fisher[5]))))/det_fish;
				d_sigy_crlb[index] = -(-(fisher[0]*(fisher[7]*(fisher[14]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])+fisher[26]*(fisher[15]*fisher[22]-fisher[21]*fisher[16]))-fisher[13]*(fisher[8]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[22]-fisher[21]*fisher[10]))+fisher[19]*(fisher[8]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[25]*(fisher[8]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])+fisher[20]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])))-fisher[6]*(fisher[1]*(fisher[14]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])+fisher[26]*(fisher[15]*fisher[22]-fisher[21]*fisher[16]))-fisher[13]*(fisher[2]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[22]-fisher[21]*fisher[4]))+fisher[19]*(fisher[2]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])))+fisher[12]*(fisher[1]*(fisher[8]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[22]-fisher[21]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[21]*fisher[28]-fisher[27]*fisher[22])-fisher[20]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[22]-fisher[21]*fisher[4]))+fisher[19]*(fisher[2]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])-fisher[8]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])-fisher[8]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[10]-fisher[9]*fisher[4])))-fisher[18]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])+fisher[26]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[28]-fisher[27]*fisher[16])-fisher[14]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[28]-fisher[27]*fisher[10])-fisher[8]*(fisher[3]*fisher[28]-fisher[27]*fisher[4])+fisher[26]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[25]*(fisher[2]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])-fisher[8]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])+fisher[14]*(fisher[3]*fisher[10]-fisher[9]*fisher[4])))+fisher[24]*(fisher[1]*(fisher[8]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])+fisher[20]*(fisher[9]*fisher[16]-fisher[15]*fisher[10]))-fisher[7]*(fisher[2]*(fisher[15]*fisher[22]-fisher[21]*fisher[16])-fisher[14]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[16]-fisher[15]*fisher[4]))+fisher[13]*(fisher[2]*(fisher[9]*fisher[22]-fisher[21]*fisher[10])-fisher[8]*(fisher[3]*fisher[22]-fisher[21]*fisher[4])+fisher[20]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))-fisher[19]*(fisher[2]*(fisher[9]*fisher[16]-fisher[15]*fisher[10])-fisher[8]*(fisher[3]*fisher[16]-fisher[15]*fisher[4])+fisher[14]*(fisher[3]*fisher[10]-fisher[9]*fisher[4]))))/det_fish);
			}
			else{						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_sigx[index] = -1;
				d_sigy[index] = -1;
 				d_framenum_all[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_sigx_crlb[index] = -1;
				d_sigy_crlb[index] = -1;
				d_llv[index] = -1;
			} 
		} //end is numeric if statement
		else{
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_sigx[index] = -1;
			d_sigy[index] = -1;
			d_framenum_all[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_sigx_crlb[index] = -1;
			d_sigy_crlb[index] = -1;
			d_llv[index] = -1;
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
	double *beta0;				// the pointer	to all the initial estimates of parameters in order ( x coord, y coord, N, sigma, Background)
	double *framenum_temp;		// pointer to temporary array of frame numbers
	double *xpeaks;			// pointer to array of xpeaks
	double *ypeaks;			// pointer to array of y peaks
	double  *d_iall;			// Pointer to image array on gpu
	double *d_beta0;			// Pointer to initial estimates on gpu
	double *d_framenum_temp;
	double *d_framenum_all;
	double *d_x_cm;				// pointer to parameters on device
	double *d_y_cm;
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
	double *xpix;
	double *ypix;
	double *d_xpix;
	double *d_ypix;
	size_t threadsperblock;


	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const size_t *idims;
	

	/* Throw an error if the input does not match expectations. */
	if (nrhs != 8) {
		printf("Must have 8 inputs ( i1, beta0, framenum_temp, xpeaks, ypeaks,xpix, ypix, numthreads)\n");
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
		printf("i1 must be a nxnxm double array\n");
		mexErrMsgTxt("See Error above!\n");

	}
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
		printf("beta0 must be a 1x6xnumel(i1(1,1,:)) double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
		printf("framenum_temp must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])){
		printf("xpeaks must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])){
		printf("ypeaks must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	
 	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[0];
	irow = (int)idims[1];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}
	beta0 = (double *)mxGetPr(prhs[1]);
	framenum_temp = (double *)mxGetPr(prhs[2]);
	xpeaks = (double *)mxGetPr(prhs[3]);
	ypeaks = (double *)mxGetPr(prhs[4]);
	xpix = (double *)mxGetPr(prhs[5]);
	ypix = (double *)mxGetPr(prhs[6]);
	printf("irow = %d icol = %d, numi = %d\n", irow, icol, numi);
	if (mxGetM(prhs[1]) != 6){		// ensure that beta0 is a 6 rowed matrix
		printf("Your beta0 vector must be of the format beta0(:,n) if referencing the nth molecule estimators.\n");
		printf("beta0 should be of the form \nbeta0(1,:) = x guesses\nbeta0(2,:) = y guesses\nbeta0(3,:) = N guesses\nbeta0(4,:) = sigma_X guesses\nbeta0(5,:) = sigy_guesses\nbeta0(6,:) = offset guesses\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetN(prhs[1]) != numi){
		printf("Your beta0 vector must have the dimensions 6 rows x numel(i2(1,1,:)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[2]) != numi && mxGetN(prhs[2]) != numi){
		printf("Your framenum_temp vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[3]) != numi && mxGetN(prhs[3]) != numi){
		printf("Your xpeaks vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[4]) != numi && mxGetN(prhs[4]) != numi){
		printf("Your ypeaks vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[5]) != irow || mxGetN(prhs[5]) != irow){
		printf("Your  xpix matrix must have the dimensions: numel(i2(:,1,1)) rows x numel(i2(1,:,1)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[6]) != irow || mxGetN(prhs[6]) != irow){
		printf("Your  ypix matrix must have the dimensions: numel(i2(:,1,1)) rows x numel(i2(1,:,1)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// verify that the input variables are what was expected
	if (irow == 1 || icol == 1){
		printf("Images are too small must have more than 1 pixel in each dimension\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (irow != icol){			// insure inputs are square matrices
		printf("Images must have equal length and width\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (nlhs != 14){
		printf("You must have 10 output variables [xf_all, yf_all, N, off_all, sigx_all, sigy_all, xf_crlb, yf_crlb, N_crlb, off_crlb, sigx_crlb, sigy_crlb, llv_all, framenum_all]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	// allocate memory on the gpu device

	cudaError_t err0 = cudaMalloc((void**)&d_xpix, irow*irow*sizeof(double));						// allocate xpix memory
	if (err0 != cudaSuccess){
		printf("%s in \n%s \nat line %d\n", cudaGetErrorString(err0), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err1 = cudaMalloc((void**)&d_ypix, irow*irow*sizeof(double));						// allocate ypix memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_iall, irow*irow*numi*sizeof(double));				// allocate image memory
	if (err2!= cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err3 = cudaMalloc((void**)&d_beta0, 6 * numi*sizeof(double));						// allocate parameter memory
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err4 = cudaMalloc((void**)&d_x_cm, numi*sizeof(double));						// allocate x center coords
	if (err4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err5 = cudaMalloc((void**)&d_y_cm, numi*sizeof(double));						// allocate y center coords
	if (err5 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err5), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err6 = cudaMalloc((void**)&d_N, numi*sizeof(double));							// allocate N array on gpu
	if (err6 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err6), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err7 = cudaMalloc((void**)&d_off, numi*sizeof(double));						// allocate offset array on gpu
	if (err7 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err7), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err8 = cudaMalloc((void**)&d_yf_all, numi*sizeof(double));						// allocate yf_all array on gpu
	if (err8 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err8), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err9 = cudaMalloc((void**)&d_xf_all, numi*sizeof(double));						// allocate xf_all array on gpu
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err10 = cudaMalloc((void**)&d_framenum_temp, numi*sizeof(double));				// allocate framenum_temp array on gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err11 = cudaMalloc((void**)&d_framenum_all, numi*sizeof(double));				// allocate framenum_all array on gpu
	if (err11 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err11), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	cudaError_t err12 = cudaMalloc((void**)&d_sigx_all, numi*sizeof(double));				// allocate sigx_all array on gpu
	if (err12 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err12), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	cudaError_t err13 = cudaMalloc((void**)&d_sigy_all, numi*sizeof(double));				// allocate sigy_all array on gpu
	if (err13 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err13), __FILE__, __LINE__); 
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err14 = cudaMalloc((void**)&d_xf_crlb, numi*sizeof(double));					// Allocate xf_crlb array on gpu
	if (err14 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err14), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err15 = cudaMalloc((void**)&d_yf_crlb, numi*sizeof(double));					// allocate yf_crlb array on gpu
	if (err15 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err15), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err16 = cudaMalloc((void**)&d_N_crlb, numi*sizeof(double));					// allocate N_crlb array on gpu
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err17 = cudaMalloc((void**)&d_off_crlb, numi*sizeof(double));					// allocate Off_crlb array on gpu
	if (err17 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err17), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err18 = cudaMalloc((void**)&d_sigx_crlb, numi*sizeof(double));					// allocate sigx_crlb array on gpu
	if (err18 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err18), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err19 = cudaMalloc((void**)&d_sigy_crlb, numi*sizeof(double));					// allocate sigy_crlb array on gpu
	if (err19 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err19), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err20 = cudaMalloc((void**)&d_llv, numi*sizeof(double));					// allocate llv array on gpu
	if (err20 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err20), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	// copy data from host to device
	cudaError_t err21 = cudaMemcpy(d_iall, iall, icol*irow*numi*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err21 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err21), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err22 = cudaMemcpy(d_beta0, beta0, 6 * numi*sizeof(double), cudaMemcpyHostToDevice);		// copy parameter data to gpu
	if (err22 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err22), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err23 = cudaMemcpy(d_x_cm, xpeaks, numi*sizeof(double), cudaMemcpyHostToDevice);		// copy x center coords to gpu
	if (err23 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err23), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err24 = cudaMemcpy(d_y_cm, ypeaks, numi*sizeof(double), cudaMemcpyHostToDevice);		// copy y center coords to gpu
	if (err24 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err24), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err25 = cudaMemcpy(d_framenum_temp, framenum_temp, numi*sizeof(double), cudaMemcpyHostToDevice);					// copy framenum_temp array on gpu
	if (err25 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err25), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err26 = cudaMemcpy(d_xpix, xpix, irow*irow*sizeof(double), cudaMemcpyHostToDevice);
	if (err26 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err26), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err27 = cudaMemcpy(d_ypix, ypix, irow*irow*sizeof(double), cudaMemcpyHostToDevice);
	if (err27 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err27), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	/* Run GPU kernel*/
	threadsperblock = mxGetScalar(prhs[7]);  // get number of threads perblock from matlab
	printf("number of threads perblock %d\n",threadsperblock);
	localize7 << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_beta0, d_x_cm, d_y_cm, d_framenum_temp, d_xf_all, d_yf_all, d_N, d_off, d_sigx_all, d_sigy_all, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, d_sigx_crlb, d_sigy_crlb, irow, d_xpix, d_ypix, d_llv, numi);

	
	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/
	cudaError_t errk2 = cudaThreadSynchronize();
	if (errk2 != cudaSuccess){
		printf("%s in %s at line %u\n", cudaGetErrorString(errk2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t errk1 = cudaPeekAtLastError();
	
	if (errk1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errk1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	


	plhs[0] = mxDuplicateArray( prhs[3]);
	mxArray *xf_all = (mxArray *)mxGetPr(plhs[0]);
	cudaError_t err28 = cudaMemcpy( xf_all, d_xf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err28 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err28), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}
	
	plhs[1] = mxDuplicateArray(prhs[3]);
	double *yf_all = (double *)mxGetPr(plhs[1]);
	cudaError_t err29 = cudaMemcpy(yf_all, d_yf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_all data
	if (err29 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err29), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[2] = mxDuplicateArray(prhs[3]);
	double *N = (double *)mxGetPr(plhs[2]);
	cudaError_t err30 = cudaMemcpy(N, d_N, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N data
	if (err30 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err30), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[3] = mxDuplicateArray(prhs[3]);
	double *off_all = (double *)mxGetPr(plhs[3]);
	cudaError_t err31 = cudaMemcpy(off_all, d_off, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy off_all data
	if (err31 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err31), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}
	
	plhs[4] = mxDuplicateArray(prhs[3]);
	double *sig_x = (double *)mxGetPr(plhs[4]);
	cudaError_t err32 = cudaMemcpy(sig_x, d_sigx_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigx data
	if (err32 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err32), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[5] = mxDuplicateArray(prhs[3]);
	double *sig_y = (double *)mxGetPr(plhs[5]);
	cudaError_t err33 = cudaMemcpy(sig_y, d_sigy_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigy data
	if (err33 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err33), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[6] = mxDuplicateArray(prhs[3]);
	double *xf_crlb = (double *)mxGetPr(plhs[6]);
	cudaError_t err34 = cudaMemcpy(xf_crlb, d_xf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_crlb data
	if (err34 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err34), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	};

	plhs[7] = mxDuplicateArray(prhs[3]);
	double *yf_crlb = (double *)mxGetPr(plhs[7]);
	cudaError_t err35 = cudaMemcpy(yf_crlb, d_yf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_crlb data
	if (err35 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err35), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[8] = mxDuplicateArray(prhs[3]);
	double *N_crlb = (double *)mxGetPr(plhs[8]);
	cudaError_t err36 = cudaMemcpy(N_crlb, d_N_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N_crlb data
	if (err36 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err36), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[9] = mxDuplicateArray(prhs[3]);
	double *off_crlb = (double *)mxGetPr(plhs[9]);
	cudaError_t err37 = cudaMemcpy(off_crlb, d_off_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy off_crlb data
	if (err37 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err37), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[10] = mxDuplicateArray(prhs[3]);
	double *sigx_crlb = (double *)mxGetPr(plhs[10]);
	cudaError_t err38 = cudaMemcpy(sigx_crlb, d_sigx_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigx_crlb data
	if (err38 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err38), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[11] = mxDuplicateArray(prhs[3]);
	double *sigy_crlb = (double *)mxGetPr(plhs[11]);
	cudaError_t err39 = cudaMemcpy(sigy_crlb, d_sigy_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigy_crlb data
	if (err39 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err39), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[12] = mxDuplicateArray(prhs[3]);
	double *llv = (double *)mxGetPr(plhs[12]);
	cudaError_t err40 = cudaMemcpy(llv, d_llv, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy llv data
	if (err40 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err40), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[13] = mxDuplicateArray(prhs[3]);
	double *framenum_all = (double *)mxGetPr(plhs[13]);
	cudaError_t err41 = cudaMemcpy(framenum_all, d_framenum_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy framenum_all data
	if (err41 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err41), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	//cudaDeviceReset();
	
	cudaFree(d_iall);
	cudaFree(d_beta0); 
	cudaFree(d_N); 
	cudaFree(d_framenum_temp); 
	cudaFree(d_x_cm); 
	cudaFree(d_y_cm);
	cudaFree(d_framenum_all); 
	cudaFree(d_xf_all); 
	cudaFree(d_yf_all); 
	cudaFree(d_off);
	cudaFree(d_xpix);
	cudaFree(d_ypix);
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
