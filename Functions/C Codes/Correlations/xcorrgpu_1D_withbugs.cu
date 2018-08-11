/*
 * xcorrgpu is a gpu implementation of the cross correlation function for images
 * Each thread will handle a different displacement of the images over eachother
 * This will only handle 2 images at a time so stacks of images need to be handled in matlab
 * The call for this function will be [C] = xcorrgpu(i1,i2); 
 * i2 is moved over i1
 */
#include <mex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

// includes, project
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#define	TPB 20

__global__ void xcorr(float *d_i1,
	float *d_i2,
	float *d_icorr,
	int m1,
	int n1,
	int m2,
	int n2)
{
	//grab theadID location
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// get output values based on block and thread locations
	int row_out = blockIdx.y*blockDim.y + ty;		
	int col_out = blockIdx.x*blockDim.x + tx;

	// Get starting value for the convolution as dictated by m2 and n2
	// we'll use i1 indicies as the coord syst.
	int row_st = row_out - (m2 - 1);
	int col_st = col_out - (n2 - 1);

	// correlation variable
	float corr=0; // initialize correlation variable

	if (row_out >= 0 && row_out < m1 + m2 - 1 && col_out >= 0 && col_out < n1 + n2 - 1)  // ensure output is within bounds of correlation image
	{
		// Buffering into memory would be 1 call to a global variable, From there we need 1 call for each multiplication, however we only need to make 1 call to a global
		// variable for the multiplication and move on, as such it doesn't make sense to buffer these images into local memory
		for (int i = 0; i < m2; i++) { // 
			for (int j = 0; j < n2; j++)
			{
				if (row_st + i >= 0 && row_st + i < m1 && col_st + j >= 0 && col_st + j < n1) { // if row start and col start are greater than 0 and less than the number of pixels available perform convolution
					corr += d_i1[row_st + i + (col_st + j) * m1] * d_i2[i + j * m2];
				}
				else {} // if else is invoked it's because row_st and col_st are outside of im1 bounds and the convolution should be left alone
			}
		}
		d_icorr[row_out + col_out*(m1 + m2 - 1)] = corr; // assign correlation variable to proper location in final image
	}
	else{}
	
}

// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	// Variable Declaration
	float *i1, *i2, *icorr;    // cpu pointers
	float *d_i1, *d_i2, *d_icorr; // gpu pointers

	// Get memory size of images
	const size_t *dim1, *dim2;
	dim1 = mxGetDimensions(prhs[0]);
	int m1 = (int)dim1[0];
	int n1 = (int)dim1[1];
	dim2 = mxGetDimensions(prhs[1]);
	int m2 = (int)dim2[0];
	int n2 = (int)dim2[1];

	// Determine memory requirements
	const int mem_size1 = m1*n1*sizeof(float);
	const int mem_size2 = m2*n2*sizeof(float);
	const int mem_size3 = (m1 + m2 - 1)*(n1 + n2 - 1)*sizeof(float);
	
	//error sections
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:xcorrgpu:rhs", "This function only accepts 2 images");
	}

	// Grab pointers from matlab	
	i1 = (float *)mxGetPr(prhs[0]);
	i2 = (float *)mxGetPr(prhs[1]);

	// allocate space on device for signal
	checkCudaErrors(cudaMalloc((void**)&d_i1, mem_size1));
	checkCudaErrors(cudaMalloc((void**)&d_i2, mem_size2));
	checkCudaErrors(cudaMalloc((void**)&d_icorr, mem_size3));
	// Copy data over to device
	checkCudaErrors(cudaMemcpy(d_i1, i1, mem_size1, cudaMemcpyHostToDevice)); // transfer image 1
	checkCudaErrors(cudaMemcpy(d_i2, i2, mem_size2, cudaMemcpyHostToDevice)); // transfer image 2

	/* Run GPU kernel*/
	dim3 dimBlock(TPB, TPB);
	dim3 dimGrid((m1 + m2 - 1) /TPB + 1, (n1 + n2 - 1) /TPB + 1);

	xcorr << <dimGrid, dimBlock >> > (d_i1, d_i2, d_icorr, m1, n1, m2, n2);
	
	// Copy data over from device
	plhs[0] = mxCreateNumericMatrix((m1 + m2 - 1), (n1 + n2 - 1), mxSINGLE_CLASS, mxREAL); // create numeric array for correlation image
	icorr = (float *)mxGetPr(plhs[0]); // grab pointer for this array

	checkCudaErrors(cudaMemcpy(icorr, d_icorr, mem_size3, cudaMemcpyDeviceToHost)); // copy device data to host array

	// clean GPU to prevent memory leaks
	cudaFree(d_i1);
	cudaFree(d_i2);
	cudaFree(d_icorr);
}

