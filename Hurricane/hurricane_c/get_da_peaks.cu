/*
 * get_da_peaks is a gpu_accelerated local maxima finder
 * [iprod] = get_da_peaks(i1, r, thresh);
 * Written by Andrew Nelson 7/20/17
 * 
 * 
 *
 *
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

__global__ void da_peaks(float *d_i1,
	float thresh,
	int m,
	int n,
	int o)
{
	
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	float d_i2[25];
	// location of output pixel being analyzed
	int row_output = blockIdx.y*blockDim.y + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*blockDim.x + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z;
	if (imnum < o && row_output >=2 && row_output < m-2 && col_output >=2 && col_output <n-2)
	{
		// buffer the info into
		for(int i = 0; i <5 ; i++){
			for(int j = 0; j <5 ; j++)
			{
				d_i2[i*5 + j] = d_i1[(row_output - 2 + i) + (col_output - 2 +j)*m + imnum*m*n];
			}
		}			
		float me = d_i2[12];
		int maxi = 1;
		if(me < thresh){maxi = 0;}
		for(int k = 0; k <25; k++)
		{
			if(d_i2[k] > me){maxi = 0;}
		}
		d_i1[row_output + col_output*m + imnum*m*n] = maxi;		
	}
	else if(imnum <o){d_i1[row_output + col_output*m + imnum*m*n] = 0;}
	else{}
}

// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	
	// Variable Declaration
	float thresh;
	float *i1, *i2;
	float *d_i1;

	// Get memory size of signal
	const size_t *dims;
	dims = mxGetDimensions(prhs[0]);
	int m = (int)dims[0];
	int n = (int)dims[1];
	int o = (int)dims[2];
	if (o < 1 || o > 10000) {
		o = 1;
	}
	const int mem_size = m*n*o*sizeof(float);
	

	// allocate space on host for data	
	i1 = (float *)mxGetPr(prhs[0]);
	thresh = (float)mxGetScalar(prhs[1]);

	// allocate space on device for signal
	checkCudaErrors(cudaMalloc((void**)&d_i1, mem_size));
	// Copy data over to device
	checkCudaErrors(cudaMemcpy(d_i1, i1, mem_size, cudaMemcpyHostToDevice));

	/* Run GPU kernel*/

	dim3 dimBlock(20, 20);
	dim3 dimGrid((n - 1) / 20 + 1, (m - 1) / 20 + 1, o);

	da_peaks << <dimGrid, dimBlock >> > (d_i1, thresh, m, n, o);
	
	// Copy data over from device

	plhs[0] = mxDuplicateArray(prhs[0]);
	//plhs[0] = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
	i2 = (float *)mxGetPr(plhs[0]);

	checkCudaErrors(cudaMemcpy(i2, d_i1, mem_size, cudaMemcpyDeviceToHost));

	cudaFree(d_i1);
}

