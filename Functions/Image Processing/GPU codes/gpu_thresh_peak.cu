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
/*
COVOLUTION SECTION
*/
/*****************************************************************************************************************************************************************************************/

/*****COVOLVE 11********************************************************************************************************/

void __global__ pick11(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_ifin,
	double *thresh,
	int irow,
	int icol,
	int numi)
{
	// Declare variables
	int pixw = 11;
	__shared__ double d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image

	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord

	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	
	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - 5;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - 5;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1

	

	// build shared image into d_i2											// THIS HAS BEEN CHECKED TO BE LOADING CORRECTLY					
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	// convolution calculation
	double d_res = 0.0*PI;		// initialize counting variable on thread register
	double d_max = 0;
	int row, col;
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		for (int rowcount = 0; rowcount < pixw; rowcount++){
			for (int colcount = 0; colcount < pixw; colcount++){
				d_res = d_i2[rowcount + ty][colcount + tx];
				if (d_res >= d_max) { 
					d_max = d_res;
					row = rowcount;
					col = colcount;
				}
			}
		}
		if (row_output < irow && col_output < icol){

			if (row == 5 & col == 5 & d_max >= *thresh) {
				d_ifin[row_output + col_output*irow + imnum*irow*icol] = 1;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin
			}
			else{ d_ifin[row_output + col_output*irow + imnum*irow*icol] = 0; }
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
	double *iall;				// the pointer to the array of all images to be analyzed
	double  *d_iall;			// Pointer to image array on gpu
	double *d_ifin;
	double *thresh;
	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const size_t *idims;
	cudaDeviceReset();

	
	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[1];
	irow = (int)idims[0];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}
	thresh = (double *)mxGetPr(prhs[1]);





	if (nlhs != 1){
		printf("You must have 1 output variables [i_peak]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	// allocate memory on the gpu device



	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*icol*numi*sizeof(double));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err3 = cudaMalloc((void**)&d_ifin, irow*icol*numi*sizeof(double));						// allocate completed image memory
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	
		// copy data from host to device
	cudaError_t err9 = cudaMemcpy(d_iall, iall, irow*icol*numi*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	
	/* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH);
	dim3 dimGrid((icol - 1) / O_TILE_WIDTH + 1, (irow - 1) / O_TILE_WIDTH + 1, numi);

	//							COVOLUTION SECTION IS COMMENTED OUT TO WORK ON EROSION
	// RUN THE CONVOLUTION WITH KERNEL 

	pick11 << < dimGrid, dimBlock >> >(d_iall, d_ifin, thresh, irow, icol, numi);




	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/
	

	plhs[0] = mxDuplicateArray(prhs[0]);
	double *ifin = (double *)mxGetPr(plhs[0]);
	printf("irow %f, icol %d, numi %d, line %d\n", thresh, icol, numi, __LINE__);
	cudaError_t err16 = cudaMemcpy(ifin, d_ifin, irow*icol*numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	//	printf("irow %f, icol %f, numi %f, line %d\n", ifin[0], ifin[1], ifin[2], __LINE__);

	// cudaDeviceReset();
	 
	
	cudaFree(d_iall);
	cudaFree(d_ifin);
	
}