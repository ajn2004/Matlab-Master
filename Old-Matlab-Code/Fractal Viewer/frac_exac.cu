/*
 * frac_exac is a gpu accelerated fractal development tool
	
	[image] frac_exac(start,mag,C,pixels)
 */
#define O_TILE_WIDTH 25
#define BLOCK_WIDTH (O_TILE_WIDTH)
#include <mex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

void __global__ g_frac(double *d_zr,
	double *d_zi,
	double *d_z0r,
	double *d_z0i,
	double *d_count,
	double thresh,
	double iter,
	int m,
	int n) {
	int indx = blockDim.x*blockIdx.x + threadIdx.x;
	int indy = blockIdx.y*blockDim.y + threadIdx.y;
	if (indx < n && indy < m) {
		double zr = d_zr[indy*n + indx];
		double zi = d_zi[indy*n + indx];
		double z0r = d_z0r[indy*n + indx];
		double z0i = d_z0i[indy*n + indx];
		double count = 0;
		double z1r, z1i;
		double thr = thresh*thresh;
//		double increm;
		double abz; 
		for (int i = 0; i < iter; i++) {
			z1r = zr;
			z1i = zi;
			zr = pow(z1r,2) - pow(z1i, 2) + z0r;
			zi = 2 * z1r*z1i + z0i;
			abz = pow(zr , 2) + pow(zi , 2) ;
			if (abz <= thr) {
				count += 1;
			}
			
		} // end for 

		d_count[indy*n + indx] = count;
	} // end check if

} // end g_frac

// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

	/*
	*
	* Variable Declaration and Setup
	*
	*
	*/ 
	// array doubles follow convention i is for input o is for output c is for crlb
	double *zr, *zi,*z0r, *z0i, thresh, iter, *count;
	// GPU Variables
	double *d_zr, *d_zi, *d_z0r, *d_z0i, *d_count;

	// Error Message Array

	// Vector Size Array
	if (nrhs != 6) {
		mexErrMsgTxt("You need 6 input variables");
	}
	if (nlhs != 1) {
		mexErrMsgTxt("You need 1 output variables");
	}
	
	// Check that variables are doubles
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}


	
	// Grab dimension size of data
	const size_t *dims;
	dims = mxGetDimensions(prhs[0]);
	int m = (int)dims[0];
	int n = (int)dims[1];
	const int mem_size = m*n*sizeof(double);

	// Get position of Data
	zr = (double *)mxGetPr(prhs[0]);
	zi = (double *)mxGetPr(prhs[1]);
	z0r = (double *)mxGetPr(prhs[2]);
	z0i = (double *)mxGetPr(prhs[3]);
	thresh = mxGetScalar(prhs[4]);
	iter = mxGetScalar(prhs[5]);
	
	
	/*
	*
	*
	*  GPU MEMORY ALLOCATION and Copying
	*  With a million molecule data set we're looking at 240 MB of data while the GeForce 1060
	*  Has ~59212 MB free or we are using ~4% of the total memory, because I've never seen a data
	*  set that big I am assuming our memory is going to be just fine
	*
	*/
	
	// cudaMalloc Array
	checkCudaErrors(cudaMalloc((void**)&d_zr, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_zi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_z0r, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_z0i, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_count, mem_size));
		
	// Data Copy Array
	checkCudaErrors(cudaMemcpy(d_zr, zr, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_zi, zi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_z0r, z0r, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_z0i, z0i, mem_size, cudaMemcpyHostToDevice));
		
	// Get Fractional error Vectors and tolerance
	/* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH); // run 2-D gpu kernel to help with indexing
	dim3 dimGrid((n - 1) / O_TILE_WIDTH + 1, (m - 1) / O_TILE_WIDTH + 1, 1);

	g_frac << <dimGrid, dimBlock >> > (d_zr, d_zi, d_z0r, d_z0i, d_count, thresh, iter, m, n);


	/*
	*
	*
	*	Copy back and free up space
	*
	*
	*
	*
	*/
	// Create Arrays at output pointers
	 plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
	

	  count = (double *)mxGetPr(plhs[0]);
	 



	// copy data array
	 checkCudaErrors(cudaMemcpy(count, d_count, mem_size, cudaMemcpyDeviceToHost));
	 
		
	// Release GPU memory
	cudaFree(d_zi);
	cudaFree(d_z0r);
	cudaFree(d_z0i);
	cudaFree(d_count);
	cudaFree(d_zr);
	
}

