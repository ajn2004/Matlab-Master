/*
* image_process.cu is a program to take matlab images and process them through background subtraction, and segmentation for localization
*
*  v 0.2
Image convolution has been implemented successfully and checked with output from matlab. This can reliable convolve stacks of
*
V 0.1
we expect a format of [im_conv] = image_process [i1, i_gauss, i_ball];
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define O_TILE_WIDTH 20								// variable to determine how many output tiles will be considered in a block
# define BLOCK_WIDTH (O_TILE_WIDTH + (11-1))		// block width needs to be output tiles + mask_width - 1 to ensure enough pixels are covered for calculation


/*****************************************************************************************************************************************************************************************/
/*
*															Device code
*
*  Multiple instances of the same void __global__ are created to deal with multiple size elements of kernels
*/
/*****************************************************************************************************************************************************************************************/


/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/
/*
COVOLUTION SECTION
*/
/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/


/********CONVOLVE 3***************************************************************************************/

void __global__ convolve3(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(3)][(3)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

/*************COVOLVE 5**************************************************************************************/

void __global__ convolve5(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(5)][(5)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

//* COVOLVE 7 ***************************************************************************************************************/

void __global__ convolve7(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(7)][(7)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

/*****COVOLVE 9********************************************************************************************************/

void __global__ convolve9(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(9)][(9)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

/*****COVOLVE 11********************************************************************************************************/

void __global__ convolve11(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(11)][(11)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

/*****COVOLVE 13********************************************************************************************************/

void __global__ convolve13(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(13)][(13)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}

/*****COVOLVE 15********************************************************************************************************/

void __global__ convolve15(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_gauss,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_igauss[(15)][(15)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;

	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - gausw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - gausw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	// Buffer data into block
	for (int grow = 0; grow < gausw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
		for (int gcol = 0; gcol < gausw; gcol++){
			d_igauss[grow][gcol] = d_gauss[grow + gcol*gausw];
		}
	}

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

	// convolution calculation
	float d_res = 0.0*PI;		// initialize counting variable on thread register
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < gausw; rowcount++){
			for (int colcount = 0; colcount < gausw; colcount++){
				d_res += d_i2[rowcount + ty][colcount + tx] * d_igauss[rowcount][colcount];
			}
		}
		if (row_output < irow && col_output < icol){


			d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}

	}



}









/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/
/*
IMAGE EROSION SECTION
*/

/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/

/*****ERODE 3********************************************************************************************************/

void __global__ erode3(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(3)][(3)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion

/*****ERODE 5********************************************************************************************************/

void __global__ erode5(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(5)][(5)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion

/*****ERODE 7********************************************************************************************************/

void __global__ erode7(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(7)][(7)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion
/*****ERODE 9********************************************************************************************************/

void __global__ erode9(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(9)][(9)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion

/*****ERODE 11********************************************************************************************************/

void __global__ erode11(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(11)][(11)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion



/*****ERODE 13********************************************************************************************************/

void __global__ erode13(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(13)][(13)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion


/*****ERODE 15********************************************************************************************************/

void __global__ erode15(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(15)][(15)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 1000000.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_min = 1000.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] - d_iball[rowcount][colcount];		// subtract height of ball from pixel height at each location and check for lowest value
					if (d_res < d_res_min){
						d_res_min = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_min;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end erosion




/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/
/*
IMAGE DIALATION SECTION
*/
/*
*			Unsure how to match MATLAB here  - 5/11/15 AJN
*
*
/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/

/*****DIALATE 3********************************************************************************************************/

void __global__ dilate3(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(3)][(3)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = -1.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = -1.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){

					d_res = d_i2[rowcount + ty][colcount + tx] * d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation


/*****DIALATE 5********************************************************************************************************/

void __global__ dilate5(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(5)][(5)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = 0.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation



/*****DIALATE 7********************************************************************************************************/

void __global__ dilate7(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(7)][(7)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = 0.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation


/*****DIALATE 9********************************************************************************************************/

void __global__ dilate9(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(9)][(9)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = 0.0*PI;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation



/*****DIALATE 11********************************************************************************************************/

void __global__ dilate11(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(11)][(11)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = -1.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = -1.0;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation


/*****DIALATE 13********************************************************************************************************/

void __global__ dilate13(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(13)][(13)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = -1.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = -1.0;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation


/*****DIALATE 15********************************************************************************************************/

void __global__ dilate15(float *d_iall,   // the gaussian is a separable filter and be treated as such
	float *d_ball,	// makes these elements eligible for constant caching
	float *d_ifin,
	int irow,
	int icol,
	int ballw,				// width of guassian kernel
	int numi)
{
	// Declare variables

	__shared__ float d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	__shared__ float d_iball[(15)][(15)];					// preallocate space for image mask

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < numi){			// if there is an image to work on
		// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
		// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
		int row_input = row_output - ballw / 2;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
		int col_input = col_output - ballw / 2;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

		// Buffer data into block
		for (int grow = 0; grow < ballw; grow++){   // buffer gaussian kernel into block   **THIS HAS BEEN CHECKED AND VERIFIED TO WORK CORRECTLY **
			for (int gcol = 0; gcol < ballw; gcol++){
				d_iball[grow][gcol] = d_ball[grow + gcol*ballw];
			}
		}

		// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
		// row/col_input represents the row/col of the input pixel being considered by 
		// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
			d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else{
			d_i2[ty][tx] = -1.0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}																							// of pixels that will not contribute to the calculation

		__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here

		// convolution calculation
		float d_res_max = -1.0;		// initialize counting variable on thread register
		float d_res;
		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
			for (int rowcount = 0; rowcount < ballw; rowcount++){
				for (int colcount = 0; colcount < ballw; colcount++){
					d_res = d_i2[rowcount + ty][colcount + tx] + d_iball[rowcount][colcount];		// add height of ball from pixel height at each location and check for maximm value
					if (d_res > d_res_max){
						d_res_max = d_res;		// if found value is lower than housed value assign lowest value to lowest found
					}
				}	// end col loop
			} // end row loop
			if (row_output < irow && col_output < icol){


				d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_res_max;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

			}// end if output pixel exists inside image

		}  // end if ty and tx are withing output tile
	} // end if picture exists


}		// end dilation



/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/
/*

BACKGROUND SUBTRACTION STEP

*/
/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/



void __global__ bkgsub(float *d_iall,
	float *d_ibkg,
	int	numel)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < numel){
		d_iall[index] = d_iall[index] - d_ibkg[index];
		if (d_iall[index] < 0.0){
			d_iall[index] = 0.0;
		}
	}
}


/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/
/*
* Host code
*
*
*/
/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	float *iall;				// the pointer to the array of all images to be analyzed
	float *gauss_vec;
	float *ball_vec;
	float  *d_iall;			// Pointer to image array on gpu
	float *d_gauss;
	float *d_ball;
	float *d_ifin, *d_ifin2;
	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const size_t *idims, *gaudims, *balldims;
	cudaDeviceReset();

	/* Throw an error if the input does not match expectations. */
/*	if (nrhs != 3) {
		printf("Must have 3 inputs ( i1, i_gauss, i_ball) line: %d\n", __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	if (!mxIsfloat(prhs[0]) || mxIsComplex(prhs[0])){
		printf("i1 must be a nxnxm float array\n");
		mexErrMsgTxt("See Error above!\n");

	}
	if (!mxIsfloat(prhs[1]) || mxIsComplex(prhs[1])){
		printf("i_gauss must be a nxn float array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsfloat(prhs[2]) || mxIsComplex(prhs[2])){
		printf("number of threads per block must be an integer between 1 and 1024\n");
		mexErrMsgTxt("See Error above!\n");
	}
*/

	// get pointer to input arguments
	iall = (float *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[1];
	irow = (int)idims[0];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}

	// get gauss dims
	gauss_vec = (float *)mxGetPr(prhs[1]);
	gaudims = mxGetDimensions(prhs[1]);
	int gaurow = (int)gaudims[0];
	int gaucol = (int)gaudims[1];

	// get ball dims
	ball_vec = (float *)mxGetPr(prhs[2]);
	balldims = mxGetDimensions(prhs[1]);
	int balrow = (int)gaudims[0];
	int balcol = (int)gaudims[1];




	if (nlhs != 1){
		printf("You must have 1 output variables [i_erode]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	// allocate memory on the gpu device



	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*icol*numi*sizeof(float));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_gauss, gaurow*gaucol*sizeof(float));						// allocate gaussian memory
	if (err2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err3 = cudaMalloc((void**)&d_ifin, irow*icol*numi*sizeof(float));						// allocate completed image memory
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	
	cudaError_t err3eb = cudaMalloc((void**)&d_ifin2, irow*icol*numi*sizeof(float));						// allocate completed 2 image memory
	if (err3eb != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3eb), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err4 = cudaMalloc((void**)&d_ball, balrow*balcol*sizeof(float));						// allocate completed image memory
	if (err4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	// copy data from host to device
	cudaError_t err9 = cudaMemcpy(d_iall, iall, irow*icol*numi*sizeof(float), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err10 = cudaMemcpy(d_gauss, gauss_vec, gaurow*gaucol*sizeof(float), cudaMemcpyHostToDevice);		// copy gauss data to gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}


	cudaError_t err1a = cudaMemcpy(d_ball, ball_vec, balrow*balcol*sizeof(float), cudaMemcpyHostToDevice);		// copy gauss data to gpu
	if (err1a != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1a), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	/* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH);
	dim3 dimGrid((icol - 1) / O_TILE_WIDTH + 1, (irow - 1) / O_TILE_WIDTH + 1, numi);

	//							COVOLUTION SECTION IS COMMENTED OUT TO WORK ON EROSION
	// RUN THE CONVOLUTION WITH KERNEL 
	switch (gaurow)
	{
	case 3:
		convolve3 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 5:
		convolve5 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 7:
		convolve7 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 9:
		convolve9 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 11:
		convolve11 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 13:
		convolve13 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	case 15:
		convolve15 << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
		break;
	default:
		printf("Incorrect size of kernel please choose an odd numbered width between 3 and 15\n");
		mexErrMsgTxt("See Error Above\n");
	}



	// RUN IMAGE EROSION
	switch (balrow)
	{
	case 3:
		erode3 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate3 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 5:
		erode5 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate5 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 7:
		erode7 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate7 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 9:
		erode9 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate9 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 11:
		erode11 << <dimGrid, dimBlock >> > (d_ifin, d_ball, d_ifin2, irow, icol, balrow, numi);
		dilate11 << <dimGrid, dimBlock >> > (d_ifin2, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 13:
		erode13 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate13 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	case 15:
		erode15 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		dilate15 << <dimGrid, dimBlock >> > (d_iall, d_ball, d_ifin, irow, icol, balrow, numi);
		break;
	default:
		printf("Incorrect size of kernel please choose an odd numbered width between 3 and 15\n");
		mexErrMsgTxt("See Error Above\n");
	}

	int numel = irow*icol*numi;	// number of pixels in the entire image
	bkgsub<< <(numel-1)/1000+1, 1000 >> > (d_iall, d_ifin, numel);			// routiune vecor subtraction on GPU!!!!!













	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/
	


	plhs[0] = mxDuplicateArray(prhs[0]);
	float *ifin = (float *)mxGetPr(plhs[0]);

	//	printf("irow %d, icol %f, numi %f, line %d\n", numi, ifin[1], ifin[2], __LINE__);
	cudaError_t err16 = cudaMemcpy(ifin, d_iall, irow*icol*numi*sizeof(float), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	//	printf("irow %f, icol %f, numi %f, line %d\n", ifin[0], ifin[1], ifin[2], __LINE__);

	// cudaDeviceReset();
	 
	
	cudaFree(d_iall);
	cudaFree(d_gauss);
	cudaFree(d_ifin);
	cudaFree(d_ifin2);
	cudaFree(d_ball);
}