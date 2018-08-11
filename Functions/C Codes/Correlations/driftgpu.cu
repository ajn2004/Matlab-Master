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


__global__ void drift(float *d_i1,
	float *d_icorr,
	int m1,
	int n1,
	int o1)
{
	//grab theadID location
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// get output values based on block and thread locations
	int row_out = blockIdx.y*blockDim.y + ty;		
	int col_out = blockIdx.x*blockDim.x + tx;
	int im_out  = blockIdx.z;

	// Get starting value for the convolution as dictated by m2 and n2
	// we'll use i1 indicies as the coord syst.
	int row_st = row_out - (m1 - 1);
	int col_st = col_out - (n1 - 1);

	// correlation variable
	float corr=0; // initialize correlation variable

	if (row_out >= 0 && row_out < 2*m1 - 1 && col_out >= 0 && col_out < 2*n1- 1 && im_out < o1-1)  // ensure output is within bounds of correlation image
	{
		// Buffering into memory would be 1 call to a global variable, From there we need 1 call for each multiplication, however we only need to make 1 call to a global
		// variable for the multiplication and move on, as such it doesn't make sense to buffer these images into local memory
		for (int i = 0; i < m1; i++) { // 
			for (int j = 0; j < n1; j++)
			{
				if (row_st + i >= 0 && row_st + i < m1 && col_st + j >= 0 && col_st + j < n1) { // if row start and col start are greater than 0 and less than the number of pixels available perform convolution
					corr += d_i1[row_st + i + (col_st + j) * m1 + im_out*m1*n1] * d_i1[i + j * m1 + (im_out+1)*m1*n1]; // shift n+1 image over n image
				}
				else {} // if else is invoked it's because row_st and col_st are outside of im1 bounds and the convolution should be left alone
			}
		}
		d_icorr[row_out + col_out*(2*m1 - 1) + im_out*(2*m1 - 1)*(2*n1 - 1)] = corr; // assign correlation variable to proper location in final image
	}
	else{}

	
}

// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	// Variable Declaration
	float *i1, *icorr;    // cpu pointers
	float *d_i1, *d_icorr; // gpu pointers
	// Get memory size of images
	const size_t *dim1;
	int TPB;
	dim1 = mxGetDimensions(prhs[0]);
	int m1 = (int)dim1[0];
	int n1 = (int)dim1[1];
	int o1 = (int)dim1[2];
	if (o1 < 1 || o1 > 10000) { o1 = 1; }
	TPB = (int)mxGetScalar(prhs[1]);
	// Determine memory requirements
	const int mem_size1 = m1*n1*o1*sizeof(float);
	const int mem_size3 = (2*m1 - 1)*(2*n1 - 1)*(o1-1)*sizeof(float);
	
	//error sections
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:driftgpu:rhs", "This function only accepts 1 images and a threadsperblock variable");
	}
	if(TPB >32){TPB = 32;}
	// Grab pointers from matlab	
	i1 = (float *)mxGetPr(prhs[0]);

	// allocate space on device for signal
	checkCudaErrors(cudaMalloc((void**)&d_i1, mem_size1));
	checkCudaErrors(cudaMalloc((void**)&d_icorr, mem_size3));

	// Copy data over to device
	checkCudaErrors(cudaMemcpy(d_i1, i1, mem_size1, cudaMemcpyHostToDevice)); // transfer image 1

	/* Run GPU kernel*/
	dim3 dimBlock(TPB, TPB,1);
	dim3 dimGrid((2*n1 - 2) /TPB + 1, (2*m1  - 2) /TPB + 1,o1);

	drift << <dimGrid, dimBlock >> > (d_i1, d_icorr, m1, n1, o1);
	
	// Copy data over from device
	size_t dims[3] = {2*m1-1, 2*n1-1, o1-1};
	plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL); // create numeric array for correlation image
	icorr = (float *)mxGetPr(plhs[0]); // grab pointer for this array

	checkCudaErrors(cudaMemcpy(icorr, d_icorr, mem_size3, cudaMemcpyDeviceToHost)); // copy device data to host array

	// clean GPU to prevent memory leaks
	cudaFree(d_i1);
	cudaFree(d_icorr);
}

