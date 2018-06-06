/*
 * image reshaper will allow the user to divide an image up into a matrix where each row is a subimage centered around a different pixel
 * [y] = im_reshape(x,w)
 *
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
#define O_TILE_WIDTH 20								// variable to determine how many output tiles will be considered in a block
# define BLOCK_WIDTH (O_TILE_WIDTH + (7-1))

// gpu high pass filter
void __global__ imswitch(float *d_x,
	float *d_y,
	int w,
	int irow,
	int icol,
	int numi)
{
	
	__shared__ float d_xs[(BLOCK_WIDTH)][(BLOCK_WIDTH)];

	// Build GPU coordinates
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int tz = threadIdx.z;
	int v = 2 * w + 1;
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;
	int row_input = row_output - w;
	int col_input = col_output - w;
	if (tz < numi) {  // if there is an image to measure
		// Pad the region
		if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)) {		// if statement checks the row/col indices to ensure they fall onto the input image
			d_xs[ty][tx] = d_x[row_input + irow*col_input];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
		}																							// on the block
		else {
			d_xs[ty][tx] = 0;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
		}

		__syncthreads();

		if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {
			for (int i = 0; i < v; i++) {
				for (int j = 0; j < v; j++) {
					//
					if (row_output < irow && col_output < icol) {
						d_y[(i*v + j)*irow*icol*numi + tz*irow*icol + row_output + irow*col_output] = d_xs[ty + i][tx + j];
					}
				}
			}
		}
	}
}

void flt2dub(float *x,
	double *y,
	int N) 
{
	for (int i = 0; i < N; i++) {
		y[i] = (double)x[i];
	}
}


// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
	// Variable Declaration
	float *x, *dx, *dy, *y;
	double *oy;

	


		
	// Get memory size of signal
	const size_t *dims;
	dims = mxGetDimensions(prhs[0]);
	int m = (int)dims[0];
	int n = (int)dims[1];
	const int o = (int)mxGetScalar(prhs[2]);
	const int mem_size = o*m*n*sizeof(float);
	const int w = (int)mxGetScalar(prhs[1]);
	const int v = 2 * w + 1;
	// allocate space on host for data
	y = (float *)mxMalloc(mem_size*v*v);
	oy = (double *)mxMalloc(mem_size*v*v * 2);
	
	x = (float *)mxGetPr(prhs[0]);
	
	
	// allocate space on device for signal
	checkCudaErrors(cudaMalloc((void**)&dx, mem_size));	
	checkCudaErrors(cudaMalloc((void**)&dy, mem_size*v*v));
	// Copy data over to device
	
	checkCudaErrors(cudaMemcpy(dx, x, mem_size, cudaMemcpyHostToDevice));

	// at this point the memory is on the GPU ready to be manipulated
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH); // run 2-D gpu kernel to help with indexing
	dim3 dimGrid((n - 1) / O_TILE_WIDTH + 1, (m - 1) / O_TILE_WIDTH + 1, o);
	
	imswitch << <dimGrid, dimBlock >> > (dx, dy, w, m, n, o);

	plhs[0] = mxCreateDoubleMatrix(m*n*o, v*v, mxREAL);
	oy = mxGetPr(plhs[0]);

	checkCudaErrors(cudaMemcpy(y,dy,mem_size*v*v, cudaMemcpyDeviceToHost));

	flt2dub(y, oy, o*m*n*v*v);
	cudaFree(dy);	
	cudaFree(dx);
}

