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



void __global__ bkgsub(double *d_iall,
    double *d_ifin,
    double thrsh,
	int	numel)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if (index < numel){
		if (d_iall[index] - d_ifin[index] > thrsh){
			d_iall[index] = d_iall[index] - d_ifin[index];
        }
        else{
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
    double *iall, *i2;				// the pointer to the array of all images to be analyzed
    double *d_iall, *d_ifin, thresh;
	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const size_t *idims;
	cudaDeviceReset();

	
	// get pointer to input arguments
    iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
    i2 = (double *)mxGetPr(prhs[1]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
    thresh = mxGetScalar(prhs[2]);
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[1];
	irow = (int)idims[0];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}


	if (nlhs != 1){
		printf("You must have 1 output variables [i_diff]\n");
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

	cudaError_t err10 = cudaMemcpy(d_ifin, i2, irow*icol*numi*sizeof(double), cudaMemcpyHostToDevice);		// copy gauss data to gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	int numel = irow*icol*numi;	// number of pixels in the entire image
	bkgsub<< <(numel-1)/1000+1, 1000 >> > (d_iall, d_ifin, thresh, numel);			// routiune vecor subtraction on GPU!!!!!


	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/
	
	plhs[0] = mxDuplicateArray(prhs[0]);
	double *ifin = (double *)mxGetPr(plhs[0]);
	//	printf("irow %d, icol %f, numi %f, line %d\n", numi, ifin[1], ifin[2], __LINE__);
	cudaError_t err16 = cudaMemcpy(ifin, d_iall, irow*icol*numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}	
	cudaFree(d_iall);
	cudaFree(d_ifin);
}