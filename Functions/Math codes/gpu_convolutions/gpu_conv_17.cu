/*
* gpu_conv_3.cu is a program to take matlab images and a 3x3 kernel and return a 2D convolution
*
*  
* we expect a format of [im_conv] = image_process [i1, kernel];
* AJN 3/18/19
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define O_TILE_WIDTH 20								// variable to determine how many output tiles will be considered in a block
#define pix 17
# define BLOCK_WIDTH (O_TILE_WIDTH + (pix-1))		// block width needs to be output tiles + mask_width - 1 to ensure enough pixels are covered for calculation


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

/*****COVOLVE********************************************************************************************************/

void __global__ convolve(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_gauss,	// makes these elements eligible for constant caching
	double *d_ifin,
	int irow,
	int icol,
	int gausw,				// width of guassian kernel
	int numi)
{
	// Declare variables
	__shared__ double d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	double d_igauss[pix][pix];					// preallocate space for image mask

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
	__syncthreads();
	// convolution calculation
	double d_res = 0.0*PI;		// initialize counting variable on thread register
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
	double *gauss_vec;
	double  *d_iall;			// Pointer to image array on gpu
	double *d_gauss;
	double *d_ifin;
	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const size_t *idims, *gaudims;
	cudaDeviceReset();

    // Error Prevention Section
    if (nlhs != 1){
		printf("You must have 1 output variables [i_conv]\n");
		mexErrMsgTxt("See Error above!\n");
    }
    if (nrhs != 2){
        printf("You must have 2 input variables (i1, kernel)\n");
		mexErrMsgTxt("See Error above!\n");
    }

	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	irow = (int)idims[0];
	icol = (int)idims[1];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}

	// get gauss dims
	gauss_vec = (double *)mxGetPr(prhs[1]);
	gaudims = mxGetDimensions(prhs[1]);
	int gaurow = (int)gaudims[0];
	int gaucol = (int)gaudims[1];

	// allocate memory on the gpu device
	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*icol*numi*sizeof(double));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_gauss, gaurow*gaucol*sizeof(double));						// allocate gaussian memory
	if (err2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
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

	cudaError_t err10 = cudaMemcpy(d_gauss, gauss_vec, gaurow*gaucol*sizeof(double), cudaMemcpyHostToDevice);		// copy gauss data to gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	/* Run GPU kernel*//* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH);
	dim3 dimGrid((icol - 1) / O_TILE_WIDTH + 1, (irow - 1) / O_TILE_WIDTH + 1, numi);

	//							COVOLUTION SECTION IS COMMENTED OUT TO WORK ON EROSION
	// RUN THE CONVOLUTION WITH KERNEL 
	convolve << < dimGrid, dimBlock >> >(d_iall, d_gauss, d_ifin, irow, icol, gaurow, numi);
	plhs[0] = mxDuplicateArray(prhs[0]);
	double *ifin = (double *)mxGetPr(plhs[0]);
	
	//	printf("irow %d, icol %f, numi %f, line %d\n", numi, ifin[1], ifin[2], __LINE__);
	cudaError_t err16 = cudaMemcpy(ifin, d_ifin, irow*icol*numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}	
	cudaFree(d_iall);
	cudaFree(d_gauss);
	cudaFree(d_ifin);
}