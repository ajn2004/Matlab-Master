/*
* full_parallel.cu is a program to take in an array of data and an array of molecular counts and a desired width for analysis
*	V 1.0
*		we expect a format of [xf_all, yf_all, N, off_all, sigx_all, sigy_all, xf_crlb, yf_crlb, N_crlb, off_crlb, sigx_crlb, sigy_crlb, llv_all, framenum_all] = full_parallel_chain_loc [i1, a1, width, xpix, ypix, mol_nums,sigma, fx_all, bkgn];
* here iall is the stack of images containing molecular emissions that need to be segmented and localized. a1 is the corresponding stack of counting images. width is the width of a segmented image typically 7, xpix and ypix are grid variables used in localization calculation, mol_nums is a scalar number of molecules to be analyzed, sigma is the initial width, bkgn is an initial offset guess, fx_all is a fake vector size of mol_num x 1 to project the output data onto. 
*	AJN 5/3/16
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define O_TILE_WIDTH 25								// variable to determine how many output tiles will be considered in a block
# define BLOCK_WIDTH (O_TILE_WIDTH + (7-1))		// block width needs to be output tiles + mask_width - 1 to ensure enough pixels are covered for calculation




/*
* Device code
*
*
*/

__device__ double device_det(double Fisher[36])
{
	double det;
	det = Fisher[0] * (Fisher[7] * (Fisher[14] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) + Fisher[26] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[32] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]))) - Fisher[13] * (Fisher[8] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]))) + Fisher[19] * (Fisher[8] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[25] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) + Fisher[31] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])))) - Fisher[6] * (Fisher[1] * (Fisher[14] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) + Fisher[26] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[32] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]))) - Fisher[13] * (Fisher[2] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]))) + Fisher[19] * (Fisher[2] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])))) + Fisher[12] * (Fisher[1] * (Fisher[8] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[21] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) + Fisher[33] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23])) - Fisher[20] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]))) + Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) - Fisher[18] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) + Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[28] * Fisher[35] - Fisher[34] * Fisher[29]) - Fisher[27] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5])) + Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) + Fisher[24] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[32] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) + Fisher[33] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[35] - Fisher[34] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) + Fisher[33] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[35] - Fisher[34] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[35] - Fisher[34] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[35] - Fisher[34] * Fisher[5]) + Fisher[33] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[32] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[31] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])))) - Fisher[30] * (Fisher[1] * (Fisher[8] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) + Fisher[20] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[26] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]))) - Fisher[7] * (Fisher[2] * (Fisher[15] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) + Fisher[27] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17])) - Fisher[14] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]))) + Fisher[13] * (Fisher[2] * (Fisher[9] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[22] * Fisher[29] - Fisher[28] * Fisher[23]) - Fisher[21] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5])) + Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) - Fisher[19] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) + Fisher[27] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[29] - Fisher[28] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[29] - Fisher[28] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[29] - Fisher[28] * Fisher[5]) + Fisher[27] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[26] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))) + Fisher[25] * (Fisher[2] * (Fisher[9] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) + Fisher[21] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11])) - Fisher[8] * (Fisher[3] * (Fisher[16] * Fisher[23] - Fisher[22] * Fisher[17]) - Fisher[15] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5])) + Fisher[14] * (Fisher[3] * (Fisher[10] * Fisher[23] - Fisher[22] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[23] - Fisher[22] * Fisher[5]) + Fisher[21] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5])) - Fisher[20] * (Fisher[3] * (Fisher[10] * Fisher[17] - Fisher[16] * Fisher[11]) - Fisher[9] * (Fisher[4] * Fisher[17] - Fisher[16] * Fisher[5]) + Fisher[15] * (Fisher[4] * Fisher[11] - Fisher[10] * Fisher[5]))));
	return det;
}


void __global__ segment and localize7(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_a1,	// makes these elements eligible for constant caching
	double sigma,
	double *xpix,
	double *ypix,
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
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi,
	int irow,
	int icol,
	double bkgn)
{
	// Declare variables

	double d_i2[7][7];		// preallocate space for shared image
	__shared__ double xgrid[7 * 7];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[7 * 7];			// this will reduce calls to global device memory
	*/
		// Coordinate building
		int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord
	int tz = threadIdx.z;
	// location of output pixel being analyzed We do not need a localization apron
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){
		if (d_a1[row_output + irow*col_output + irow*icol*imnum] >0){ // the assumption at this point is if d_a1 has a value greater than 0 the neural net said to analyze
			int row_input = row_output - 3;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
			int col_input = col_output - 3;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1
			int index = (int)d_a1[row_output + irow*col_output + irow*icol*imnum];
			// Buffer data into block

			// buffer segment into i2

			for (int row = 0; row < 7; row++){
				for (int col = 0; col < 7; col++){
					xgrid[row][col] = d_xpix[row + 7 * col]; // load x and ypix
					ygrid[row][col] = d_ypix[row + 7 * col];

					if ((row_input + row >= 0) && (row_input + row < irow) && (col_input + col >= 0) && (col_input + col < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
						d_i2[row][col] = d_iall[row_input + row + (col_input + col)*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
					}
				}
			}// end counting of rows and columns at this point the image to localize is contained in d_i2

			// at this point we should have xpix ypix and the image to localize loaded to 1 core

			// Determine the beta estimations

			// Determining X and Y guesses

			// center of mass approach
			double xsum = 0.0;
			double ysum = 0.0;
			double sumt = 0;
			
			for (int row = 0; row < 7; row++){
				for (int col = 0; col < 7; col++){
					sumt += d_i2[row][col];  // sum can also be used to determine N guess
					xsum += xgrid[row][col] * d_i2[row][col];
					ysum += ygrid[row][col] * d_i2[row][col];

				}
			} // end counting over rows

			// final estimation of xguess and yguess as xcm and ycm
			d_beta1[0] = xsum / sumt;
			d_beta1[1] = ysum / sumt;
			d_beta1[2] = sumt;
			d_beta1[3] = sigma;
			d_beta1[4] = sigma;
			d_beta1[5] = bkgn;


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
							llv += d_i2[rowcount + colcount*irow] * log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*irow] * log(d_i2[rowcount + colcount*irow] + 0.0000000000000001) + d_i2[rowcount + colcount*irow];
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
					d_xf_all[index] = d_beta1[0] + col_output; // correct position for x
					d_yf_all[index] = d_beta1[1] + row_output; // correct position for y
					d_N[index] = d_beta1[2];
					d_sigx[index] = d_beta1[3];
					d_sigy[index] = d_beta1[4];
					d_off[index] = d_beta1[5];
					d_framenum_all[index] = imnum;
					d_llv[index] = llv;

					// calculate crlb's for estimators	

					// UPDATE FOR SIGMA VALUES
					det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
					d_xf_crlb[index] = (fisher[7] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[13] * (fisher[8] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) + fisher[19] * (fisher[8] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[25] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) + fisher[31] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[26] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])))) / det_fish;
					d_yf_crlb[index] = -(-(fisher[0] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[12] * (fisher[2] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[2] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) - fisher[24] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[30] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[26] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])))) / det_fish);
					d_N_crlb[index] = (fisher[0] * (fisher[7] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[31] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[1] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish;
					d_off_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish);
					d_sigx_crlb[index] = (fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])))) / det_fish;
					d_sigy_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[25] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) + fisher[24] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])))) / det_fish);
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


	}// end if activated 
}// end if and image
} // end gpu  segment and loc 7


void __global__ segment and localize7(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_a1,	// makes these elements eligible for constant caching
	double sigma,
	double *xpix,
	double *ypix,
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
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi,
	int irow,
	int icol,
	double bkgn)
{
	// Declare variables

	double d_i2[9][9];		// preallocate space for shared image
	__shared__ double xgrid[9 * 9];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[9 * 9];			// this will reduce calls to global device memory
	*/
		// Coordinate building
		int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord
	int tz = threadIdx.z;
	// location of output pixel being analyzed We do not need a localization apron
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){
		if (d_a1[row_output + irow*col_output + irow*icol*imnum] >0){ // the assumption at this point is if d_a1 has a value greater than 0 the neural net said to analyze
			int row_input = row_output - 3;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
			int col_input = col_output - 3;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1
			int index = (int)d_a1[row_output + irow*col_output + irow*icol*imnum];
			// Buffer data into block

			// buffer segment into i2

			for (int row = 0; row < 9; row++){
				for (int col = 0; col < 9; col++){
					xgrid[row][col] = d_xpix[row + 9 * col]; // load x and ypix
					ygrid[row][col] = d_ypix[row + 9 * col];

					if ((row_input + row >= 0) && (row_input + row < irow) && (col_input + col >= 0) && (col_input + col < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
						d_i2[row][col] = d_iall[row_input + row + (col_input + col)*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
					}
				}
			}// end counting of rows and columns at this point the image to localize is contained in d_i2

			// at this point we should have xpix ypix and the image to localize loaded to 1 core

			// Determine the beta estimations

			// Determining X and Y guesses

			// center of mass approach
			double xsum = 0.0;
			double ysum = 0.0;
			double sumt = 0;
			
			for (int row = 0; row < 9; row++){
				for (int col = 0; col < 9; col++){
					sumt += d_i2[row][col];  // sum can also be used to determine N guess
					xsum += xgrid[row][col] * d_i2[row][col];
					ysum += ygrid[row][col] * d_i2[row][col];

				}
			} // end counting over rows

			// final estimation of xguess and yguess as xcm and ycm
			d_beta1[0] = xsum / sumt;
			d_beta1[1] = ysum / sumt;
			d_beta1[2] = sumt;
			d_beta1[3] = sigma;
			d_beta1[4] = sigma;
			d_beta1[5] = bkgn;


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
							llv += d_i2[rowcount + colcount*irow] * log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*irow] * log(d_i2[rowcount + colcount*irow] + 0.0000000000000001) + d_i2[rowcount + colcount*irow];
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
					d_xf_all[index] = d_beta1[0] + col_output; // correct position for x
					d_yf_all[index] = d_beta1[1] + row_output; // correct position for y
					d_N[index] = d_beta1[2];
					d_sigx[index] = d_beta1[3];
					d_sigy[index] = d_beta1[4];
					d_off[index] = d_beta1[5];
					d_framenum_all[index] = imnum;
					d_llv[index] = llv;

					// calculate crlb's for estimators	

					// UPDATE FOR SIGMA VALUES
					det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
					d_xf_crlb[index] = (fisher[7] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[13] * (fisher[8] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) + fisher[19] * (fisher[8] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[25] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) + fisher[31] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[26] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])))) / det_fish;
					d_yf_crlb[index] = -(-(fisher[0] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[12] * (fisher[2] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[2] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) - fisher[24] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[30] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[26] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])))) / det_fish);
					d_N_crlb[index] = (fisher[0] * (fisher[7] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[31] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[1] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish;
					d_off_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish);
					d_sigx_crlb[index] = (fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])))) / det_fish;
					d_sigy_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[25] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) + fisher[24] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])))) / det_fish);
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


	}// end if activated 
}// end if and image
} // end gpu  segment and loc 9

void __global__ segment and localize11(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_a1,	// makes these elements eligible for constant caching
	double sigma,
	double *xpix,
	double *ypix,
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
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi,
	int irow,
	int icol)
{
	// Declare variables

	double d_i2[11][11];		// preallocate space for shared image
	__shared__ double xgrid[11 * 11];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[11 * 11];			// this will reduce calls to global device memory
	*/
		// Coordinate building
		int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord
	int tz = threadIdx.z;
	// location of output pixel being analyzed We do not need a localization apron
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){
		if (d_a1[row_output][col_output][imnum] >0){ // the assumption at this point is if d_a1 has a value greater than 0 the neural net said to analyze
			int row_input = row_output - 3;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
			int col_input = col_output - 3;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1
			int index = (int)d_a1[row_output][col_output][imnum];
			// Buffer data into block

			// buffer segment into i2

			for (int row = 0; row < 11; row++){
				for (int col = 0; col < 11; col++){
					xgrid[row][col] = d_xpix[row + 11 * col]; // load x and ypix
					ygrid[row][col] = d_ypix[row + 11 * col];

					if ((row_input + row >= 0) && (row_input + row < irow) && (col_input + col >= 0) && (col_input + col < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
						d_i2[row][col] = d_iall[row_input + row + (col_input + col)*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
					}
				}
			}// end counting of rows and columns at this point the image to localize is contained in d_i2

			// at this point we should have xpix ypix and the image to localize loaded to 1 core

			// Determine the beta estimations

			// Determining X and Y guesses

			// center of mass approach
			double xsum = 0.0;
			double ysum = 0.0;
			double sumt = 0;
			double mina = 1000000;
			for (int row = 0; row < 11; row++){
				for (int col = 0; col < 11; col++){
					sumt += d_i2[row][col];  // sum can also be used to determine N guess
					xsum += xgrid[row][col] * d_i2[row][col];
					ysum += ygrid[row][col] * d_i2[row][col];
					if (d_i2[row][col] < mina){ // find minimum value
						mina = d_i2;
					}
				}
			} // end counting over rows

			// final estimation of xguess and yguess as xcm and ycm
			d_beta1[0] = xsum / sumt;
			d_beta1[1] = ysum / sumt;
			d_beta1[2] = sumt;
			d_beta1[3] = sigma;
			d_beta1[4] = sigma;
			d_beta1[5] = mina;


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
				for (int rowcount = 0; rowcount < 11; rowcount++){	// FOR 2  loops over all rows
					for (int colcount = 0; colcount < 11; colcount++){	// FOR 3 loops over all columns

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
							llv += d_i2[rowcount + colcount*irow] * log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*irow] * log(d_i2[rowcount + colcount*irow] + 0.0000000000000001) + d_i2[rowcount + colcount*irow];
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
					d_xf_all[index] = d_beta1[0] + col_output; // correct position for x
					d_yf_all[index] = d_beta1[1] + row_output; // correct position for y
					d_N[index] = d_beta1[2];
					d_sigx[index] = d_beta1[3];
					d_sigy[index] = d_beta1[4];
					d_off[index] = d_beta1[5];
					d_framenum_all[index] = imnum;
					d_llv[index] = llv;

					// calculate crlb's for estimators	

					// UPDATE FOR SIGMA VALUES
					det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
					d_xf_crlb[index] = (fisher[7] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[13] * (fisher[8] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) + fisher[19] * (fisher[8] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[26] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[25] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[32] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) + fisher[31] * (fisher[8] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) + fisher[20] * (fisher[9] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[26] * (fisher[9] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) + fisher[21] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])))) / det_fish;
					d_yf_crlb[index] = -(-(fisher[0] * (fisher[14] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) + fisher[26] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[32] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]))) - fisher[12] * (fisher[2] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[20] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[2] * (fisher[15] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[26] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) - fisher[24] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[33] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[15] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[32] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[30] * (fisher[2] * (fisher[15] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) + fisher[27] * (fisher[16] * fisher[23] - fisher[22] * fisher[17])) - fisher[14] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[20] * (fisher[3] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[15] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[26] * (fisher[3] * (fisher[16] * fisher[23] - fisher[22] * fisher[17]) - fisher[15] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])))) / det_fish);
					d_N_crlb[index] = (fisher[0] * (fisher[7] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[31] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[21] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) + fisher[33] * (fisher[22] * fisher[29] - fisher[28] * fisher[23])) - fisher[19] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]))) + fisher[18] * (fisher[1] * (fisher[9] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[27] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[33] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[35] - fisher[34] * fisher[23]) - fisher[21] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[9] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[33] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[9] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[27] * (fisher[10] * fisher[23] - fisher[22] * fisher[11])) - fisher[7] * (fisher[3] * (fisher[22] * fisher[29] - fisher[28] * fisher[23]) - fisher[21] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[23] - fisher[22] * fisher[5])) + fisher[19] * (fisher[3] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[9] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[27] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[3] * (fisher[10] * fisher[23] - fisher[22] * fisher[11]) - fisher[9] * (fisher[4] * fisher[23] - fisher[22] * fisher[5]) + fisher[21] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish;
					d_off_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) + fisher[25] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) + fisher[32] * (fisher[16] * fisher[29] - fisher[28] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[29] - fisher[28] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[28] * fisher[35] - fisher[34] * fisher[29]) - fisher[26] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[29] - fisher[28] * fisher[5])) + fisher[25] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) - fisher[24] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) + fisher[32] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[35] - fisher[34] * fisher[17]) - fisher[14] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[35] - fisher[34] * fisher[11]) - fisher[8] * (fisher[4] * fisher[35] - fisher[34] * fisher[5]) + fisher[32] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) + fisher[26] * (fisher[10] * fisher[17] - fisher[16] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[16] * fisher[29] - fisher[28] * fisher[17]) - fisher[14] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[17] - fisher[16] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[10] * fisher[29] - fisher[28] * fisher[11]) - fisher[8] * (fisher[4] * fisher[29] - fisher[28] * fisher[5]) + fisher[26] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])) - fisher[25] * (fisher[2] * (fisher[10] * fisher[17] - fisher[16] * fisher[11]) - fisher[8] * (fisher[4] * fisher[17] - fisher[16] * fisher[5]) + fisher[14] * (fisher[4] * fisher[11] - fisher[10] * fisher[5])))) / det_fish);
					d_sigx_crlb[index] = (fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[31] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) + fisher[32] * (fisher[15] * fisher[23] - fisher[21] * fisher[17])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[23] - fisher[21] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[35] - fisher[33] * fisher[23]) - fisher[20] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[23] - fisher[21] * fisher[5])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) + fisher[32] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[35] - fisher[33] * fisher[17]) - fisher[14] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[35] - fisher[33] * fisher[11]) - fisher[8] * (fisher[3] * fisher[35] - fisher[33] * fisher[5]) + fisher[32] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[31] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5]))) + fisher[30] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) + fisher[20] * (fisher[9] * fisher[17] - fisher[15] * fisher[11])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[23] - fisher[21] * fisher[17]) - fisher[14] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[17] - fisher[15] * fisher[5])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[23] - fisher[21] * fisher[11]) - fisher[8] * (fisher[3] * fisher[23] - fisher[21] * fisher[5]) + fisher[20] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[17] - fisher[15] * fisher[11]) - fisher[8] * (fisher[3] * fisher[17] - fisher[15] * fisher[5]) + fisher[14] * (fisher[3] * fisher[11] - fisher[9] * fisher[5])))) / det_fish;
					d_sigy_crlb[index] = -(-(fisher[0] * (fisher[7] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) + fisher[19] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[25] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]))) - fisher[6] * (fisher[1] * (fisher[14] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) + fisher[26] * (fisher[15] * fisher[22] - fisher[21] * fisher[16])) - fisher[13] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]))) + fisher[12] * (fisher[1] * (fisher[8] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[22] - fisher[21] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[21] * fisher[28] - fisher[27] * fisher[22]) - fisher[20] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[22] - fisher[21] * fisher[4])) + fisher[19] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) - fisher[18] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) + fisher[26] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[28] - fisher[27] * fisher[16]) - fisher[14] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[28] - fisher[27] * fisher[10]) - fisher[8] * (fisher[3] * fisher[28] - fisher[27] * fisher[4]) + fisher[26] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[25] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4]))) + fisher[24] * (fisher[1] * (fisher[8] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) + fisher[20] * (fisher[9] * fisher[16] - fisher[15] * fisher[10])) - fisher[7] * (fisher[2] * (fisher[15] * fisher[22] - fisher[21] * fisher[16]) - fisher[14] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[16] - fisher[15] * fisher[4])) + fisher[13] * (fisher[2] * (fisher[9] * fisher[22] - fisher[21] * fisher[10]) - fisher[8] * (fisher[3] * fisher[22] - fisher[21] * fisher[4]) + fisher[20] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])) - fisher[19] * (fisher[2] * (fisher[9] * fisher[16] - fisher[15] * fisher[10]) - fisher[8] * (fisher[3] * fisher[16] - fisher[15] * fisher[4]) + fisher[14] * (fisher[3] * fisher[10] - fisher[9] * fisher[4])))) / det_fish);
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


	}// end if activated 
}// end if and image
} // end gpu  segment and loc 11




/*

THIS IS THE SECTION FOR IDENTIFICATION

*/





/*
* Host code
*
*
*/


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	double *iall;			// the pointer to the array of all images to be analyzed
	double *a1;
	double *sigma;
	double  *d_iall;		// Pointer to image array on gpu
	double *d_a1;			// pointer to count array on gpu
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
	double *bkgn;

	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;				// n
	int numi;				// number of images imported
	int arow;
	int acol;
	int numa;
	int widths;
	const int *idims, *adims, *xdims, *ydims;


	/* Throw an error if the input does not match expectations. */
	if (nrhs != 9) {
		printf("Must have 8 inputs ( iall, a1, width, xpix, ypix, mol_num, sigs, fx_all, bkgn) line: %d\n", __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
		printf("iall must be a m x n x numel(iall(1,1,:)) double array\n");
		mexErrMsgTxt("See Error above!\n");

	}
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
		printf("a1 must be a m x n xnumel(iall(1,1,:)) double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
		printf("Width must be a l x 1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])){
		printf("xpix must be a width x witdh double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])){
		printf("ypix must be a width x witdh double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])){
		printf("mol_num must be a 1 x 1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}


	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[1];
	irow = (int)idims[0];
	numi = (int)idims[2];  // get number of images perblock from matlab
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}

	// get dimensions of activation image
	a1 = (double *)mxGetPr(prhs[2]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	adims = mxGetDimensions(prhs[2]);	// get dimensions of image array
	acol = (int)adims[1];
	arow = (int)adims[0];
	numa = (int)adims[2];  // get number of images perblock from matlab
	if (numa > 10000000 || numa < 1){
		numa = 1;
	}

	// get width pointer
	widths = (int *)mxGetPr(prhs[2]);


	// get xpix dims
	xpix = (double *)mxGetPr(prhs[3]);
	xdims = mxGetDimensions(prhs[3]);

	// get ypix dims
	ypix = (double *)mxGetPr(prhs[4]);
	ydims = mxGetDimensions(prhs[4]);

	// get number of molecules
	mol_num = (double *)mxGetPr(prhs[5]);
	sigma = (double *)mxGetPr(prhs[6]);
	bkgn = (double *)mxGetPr(prhs[8]]);
	// EVERYONE LOVES SOME GOOD VARIABLE CHECKING!!!!!!!!!!!!!!
	if (icol != acol){
		printf("a1 and iall must have same number of columns\n");
		mexErrMsgTxt("See Above Error!\n");
	}
	if (irow != arow){
		printf("a1 and iall must have same number of rows\n");
		mexErrMsgTxt("See Above Error!\n");
	}
	if (numi != numa){
		printf("a1 and iall must have same number of frames\n");
		mexErrMsgTxt("See Above Error!\n");
	}
	if (xdims[0] != ydims[0]){
		printf("xpix and ypix must have same number of columns\n");
		mexErrMsgTxt("See Above Error!\n");
	}
	if (xdims[1] != ydims[1]){
		printf("xpix and ypix must have same number of rows\n");
		mexErrMsgTxt("See Above Error!\n");
	}

	// Did the User declare an output?
	if (nlhs != 14){
		printf("You must have 14 output variables [xf_all, yf_all, N, off_all, sigx_all, sigy_all, xf_crlb, yf_crlb, N_crlb, off_crlb, sigx_crlb, sigy_crlb, llv_all, framenum_all]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	cudaDeviceReset();
	// allocate memory on the gpu device
	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*icol*(numi)*sizeof(double));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_a1, irow*icol*numa*sizeof(double));						// allocate a1 memory
	if (err2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
		
	cudaError_t err3 = cudaMalloc((void**)&d_xpix, widths*widths*sizeof(double));				// allocate xpix
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err4 = cudaMalloc((void**)&d_ypix, widths*widths*sizeof(double));				// allocate ypix
	if (err4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err5 = cudaMalloc((void**)&d_llv, mol_num*sizeof(double));					// allocate llv array on gpu
	if (err5 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err5), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err6 = cudaMalloc((void**)&d_N, mol_num*sizeof(double));							// allocate N array on gpu
	if (err6 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err6), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err7 = cudaMalloc((void**)&d_off, mol_num*sizeof(double));						// allocate offset array on gpu
	if (err7 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err7), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err8 = cudaMalloc((void**)&d_yf_all, mol_num*sizeof(double));						// allocate yf_all array on gpu
	if (err8 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err8), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err9 = cudaMalloc((void**)&d_xf_all, mol_num*sizeof(double));						// allocate xf_all array on gpu
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err10 = cudaMalloc((void**)&d_framenum_temp, mol_num*sizeof(double));				// allocate framenum_temp array on gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err11 = cudaMalloc((void**)&d_framenum_all, mol_num*sizeof(double));				// allocate framenum_all array on gpu
	if (err11 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err11), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	cudaError_t err12 = cudaMalloc((void**)&d_sigx_all, mol_num*sizeof(double));				// allocate sigx_all array on gpu
	if (err12 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err12), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	cudaError_t err13 = cudaMalloc((void**)&d_sigy_all, mol_num*sizeof(double));				// allocate sigy_all array on gpu
	if (err13 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err13), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err14 = cudaMalloc((void**)&d_xf_crlb, mol_num*sizeof(double));					// Allocate xf_crlb array on gpu
	if (err14 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err14), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err15 = cudaMalloc((void**)&d_yf_crlb, mol_num*sizeof(double));					// allocate yf_crlb array on gpu
	if (err15 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err15), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err16 = cudaMalloc((void**)&d_N_crlb, mol_num*sizeof(double));					// allocate N_crlb array on gpu
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err17 = cudaMalloc((void**)&d_off_crlb, mol_num*sizeof(double));					// allocate Off_crlb array on gpu
	if (err17 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err17), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err18 = cudaMalloc((void**)&d_sigx_crlb, mol_num*sizeof(double));					// allocate sigx_crlb array on gpu
	if (err18 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err18), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err19 = cudaMalloc((void**)&d_sigy_crlb, mol_num*sizeof(double));					// allocate sigy_crlb array on gpu
	if (err19 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err19), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}



	// copy data from host to device
	cudaError_t err20 = cudaMemcpy(d_iall, iall, irow*icol*(numi)*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err20 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err20), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err21 = cudaMemcpy(d_a1, a1, arow*acol*numa*sizeof(double), cudaMemcpyHostToDevice);		// copy a1 data to gpu
	if (err21 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err21), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err22 = cudaMemcpy(d_xpix, xpix, widths*widths*sizeof(double), cudaMemcpyHostToDevice);		// copy xpix data to gpu
	if (err22 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err22), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err23 = cudaMemcpy(d_ypix, ypix, widths*widths*sizeof(double), cudaMemcpyHostToDevice);		// copy ypix data to gpu
	if (err23 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err23), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}


	/* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH); // run 2-D gpu kernel to help with indexing
	dim3 dimGrid((icol - 1) / O_TILE_WIDTH + 1, (irow - 1) / O_TILE_WIDTH + 1, numi );
	


	switch (widths)
	{
	case 7:
		segment_and_localize7 << < dimGrid, dimBlock >> >(d_iall, d_a1, sigma, d_xpix, d_ypix, d_xf_all, d_yf_all, d_N, d_off, d_sigx_all, d_sigy_all, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, d_sigx_crlb, d_sigy_crlb, d_xpix, d_ypix, d_llv, numi, irow, icol, bkgn);
		break;
	case 9:
		segment_and_localize9 << < dimGrid, dimBlock >> >(d_iall, d_a1, sigma, d_xpix, d_ypix, d_xf_all, d_yf_all, d_N, d_off, d_sigx_all, d_sigy_all, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, d_sigx_crlb, d_sigy_crlb, d_xpix, d_ypix, d_llv, numi, irow, icol, bkgn);
		break;
	case 11:
		segment_and_localize11 << < dimGrid, dimBlock >> >(d_iall, d_a1, sigma, d_xpix, d_ypix, d_xf_all, d_yf_all, d_N, d_off, d_sigx_all, d_sigy_all, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, d_sigx_crlb, d_sigy_crlb, d_xpix, d_ypix, d_llv, numi, irow, icol, bkgn);
		break;
	default:
		printf("Image size is inappropriate please choose either 7x7, 9x9, or 11x11 size\n");
		mexErrMsgTxt("See Error Above!\n");
		break;
	}



	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/

/*
	cudaError_t errk1 = cudaPeekAtLastError();
	if (errk1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errk1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	*/
	cudaError_t err24 = cudaThreadSynchronize();
	if (err24 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err24), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}


	plhs[11] = mxDuplicateArray(prhs[7]);
	double *sigy_crlb = (double *)mxGetPr(plhs[11]);
	cudaError_t err25 = cudaMemcpy(sigy_crlb, d_sigy_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigy_crlb data
	if (err25 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err25), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[12] = mxDuplicateArray(prhs[7]);
	double *llv = (double *)mxGetPr(plhs[12]);
	cudaError_t err26 = cudaMemcpy(llv, d_llv, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy llv data
	if (err26 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err26), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[13] = mxDuplicateArray(prhs[7]);
	double *framenum_all = (double *)mxGetPr(plhs[13]);
	cudaError_t err27 = cudaMemcpy(framenum_all, d_framenum_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy framenum_all data
	if (err27 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err27), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[0] = mxDuplicateArray(prhs[7]);
	mxArray *xf_all = (mxArray *)mxGetPr(plhs[0]);
	cudaError_t err28 = cudaMemcpy(xf_all, d_xf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err28 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err28), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[1] = mxDuplicateArray(prhs[7]);
	double *yf_all = (double *)mxGetPr(plhs[1]);
	cudaError_t err29 = cudaMemcpy(yf_all, d_yf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_all data
	if (err29 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err29), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[2] = mxDuplicateArray(prhs[7]);
	double *N = (double *)mxGetPr(plhs[2]);
	cudaError_t err30 = cudaMemcpy(N, d_N, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N data
	if (err30 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err30), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[3] = mxDuplicateArray(prhs[7]);
	double *off_all = (double *)mxGetPr(plhs[3]);
	cudaError_t err31 = cudaMemcpy(off_all, d_off, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy off_all data
	if (err31 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err31), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[4] = mxDuplicateArray(prhs[7]);
	double *sig_x = (double *)mxGetPr(plhs[4]);
	cudaError_t err32 = cudaMemcpy(sig_x, d_sigx_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigx data
	if (err32 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err32), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[5] = mxDuplicateArray(prhs[7]);
	double *sig_y = (double *)mxGetPr(plhs[5]);
	cudaError_t err33 = cudaMemcpy(sig_y, d_sigy_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigy data
	if (err33 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err33), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[6] = mxDuplicateArray(prhs[7]);
	double *xf_crlb = (double *)mxGetPr(plhs[6]);
	cudaError_t err34 = cudaMemcpy(xf_crlb, d_xf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_crlb data
	if (err34 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err34), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	};

	plhs[7] = mxDuplicateArray(prhs[7]);
	double *yf_crlb = (double *)mxGetPr(plhs[7]);
	cudaError_t err35 = cudaMemcpy(yf_crlb, d_yf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_crlb data
	if (err35 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err35), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[8] = mxDuplicateArray(prhs[7]);
	double *N_crlb = (double *)mxGetPr(plhs[8]);
	cudaError_t err36 = cudaMemcpy(N_crlb, d_N_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N_crlb data
	if (err36 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err36), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[9] = mxDuplicateArray(prhs[7]);
	double *off_crlb = (double *)mxGetPr(plhs[9]);
	cudaError_t err37 = cudaMemcpy(off_crlb, d_off_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy off_crlb data
	if (err37 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err37), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[10] = mxDuplicateArray(prhs[7]);
	double *sigx_crlb = (double *)mxGetPr(plhs[10]);
	cudaError_t err38 = cudaMemcpy(sigx_crlb, d_sigx_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy sigx_crlb data
	if (err38 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err38), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	//cudaDeviceReset();

	cudaFree(d_iall);
	cudaFree(d_a1);
	cudaFree(d_N);
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

	return;
}