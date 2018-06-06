/*
 * cgpufourier is a program that will perform fourier 
 * transforms of a given variable on the gpu and subject them to a high pass filter 
 * this will return a "background subtracted" image
 * Calling in matlab will look like 
 * [iprod] = cgpufourier(i1, sp);
 * Written by Andrew Nelson 6/1/17
 * 
 * Need much better comments than this
 *
 *
 */
#include <mex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <helper_functions.h>
#include <helper_cuda.h>


// gpu high pass filter
void __global__ highpass(cufftComplex *data,
	int sp,
	int m,
	int n)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index < n){
		/*double s = 0.0;
		double ms = 0;
		for (int i = 0; i<m; i++){
			ms = ms + powf(powf(data[index*m + i].x,2) + powf(data[index*m + i].y, 2),0.5);
		}
		ms = ms / (double)m;
		*/
		for (int i =0; i<sp*3; i++){
			data[index*m + i].x = powf(1 + exp((double)i-(double)sp),-1)*data[index*m+i].x;
			data[index*m + i].y = powf(1 + exp((double)i-(double)sp),-1)*data[index*m+i].y;
			data[index*m + (m-i-1)].y = powf(1 + exp((double)i-(double)sp),-1)*data[index*m + (m-i-1)].y;
			data[index*m + (m-i-1)].x = powf(1 + exp((double)i-(double)sp),-1)*data[index*m + (m-i-1)].x;
		}/*
		double as = 0;
		for (int i = 0; i<m; i++){
			as = as + powf(powf(data[index*m + i].x, 2) + powf(data[index*m + i].y, 2), 0.5);;
		}
		as = as/(double)m;
		s = ms/as;
		for (int i = 0; i<m; i++){
			data[index*m + i].x = s*data[index*m + i].x;
			data[index*m + i].y = s*data[index*m + i].y;
		}*/
	}
}

// gpu scale function
void __global__ scaleit(cufftComplex *data,
	int top)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index < top)
	{
		data[index].x = data[index].x/1000.0;
		data[index].y = data[index].y/1000.0;
	}
}

// convert data to handle complex
void pack_r2c(cufftComplex  *output_float,
	double *input_re,
	int Ntot)
{
	int i;
	for (i = 0; i < Ntot; i++)
	{
		output_float[i].x = (float)input_re[i];
		output_float[i].y = 0.0;
	}
}

void pack_c2c(cufftComplex  *output_float,
	double *input_re,
	double *input_im,
	int Ntot)
{
	int i;
	for (i = 0; i < Ntot; i++)
	{
		output_float[i].x = (float)input_re[i];
		output_float[i].y = (float)input_im[i];	}
}

void unpack_c2c(cufftComplex  *input_float,
	double *output_re,
	double *output_im,
	int Ntot)
{
	int i;
	for (i = 0; i < Ntot; i++)
	{
		output_re[i] = (double)input_float[i].x;
		output_im[i] = (double)input_float[i].y;
	}
}

// main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{

	// Variable Declaration
	double *yr,*yor, *yoi;
	int sp;
	sp = mxGetScalar(prhs[1]);
	cufftComplex *data, *d_data, *y_out;

	
	printf("starting Fourier...\n");

		
	// Get memory size of signal
	const size_t *dims;
	dims = mxGetDimensions(prhs[0]);
	int m = (int)dims[0];
	int n = (int)dims[1];
	const int mem_size = m*n*sizeof(cufftComplex);
	 
	
	// allocate space on host for data
	data = (cufftComplex *)mxMalloc(mem_size);
	y_out = (cufftComplex *)mxMalloc(mem_size);
	yr = (double *)mxGetPr(prhs[0]);
	// arrange the input to be complex data
	pack_r2c(data,yr,n*m);
	
	
	// allocate space on device for signal
	checkCudaErrors(cudaMalloc((void**)&d_data, mem_size));	
	// Copy data over to device
	checkCudaErrors(cudaMemcpy(d_data, data, mem_size, cudaMemcpyHostToDevice));
	
	// at this point the memory is on the GPU ready to be manipulated
	cufftHandle plan;
	checkCudaErrors(cufftPlan1d(&plan, m, CUFFT_C2C, n));
	checkCudaErrors(cufftExecC2C(plan, d_data, d_data, CUFFT_FORWARD));
	
	// at this point the fourier transform is complete and sitting on the gpu
	highpass <<<((n-1)/1024 +1) ,1024>>>(d_data,sp,m,n);
	checkCudaErrors(cufftExecC2C(plan, d_data, d_data, CUFFT_INVERSE));
	//scaleit <<<((m*n-1)/1024 +1) ,1024>>>(d_data,m*n);
	// Copy data over from device
	
	checkCudaErrors(cudaMemcpy(y_out,d_data,mem_size, cudaMemcpyDeviceToHost));
	cufftDestroy(plan);
	cudaFree(d_data);
	// create complex double in matlab
	plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);


	yor = mxGetPr(plhs[0]);
	yoi = mxGetPi(plhs[0]);
	unpack_c2c(y_out, yor, yoi, n*m);
	mxFree(data);
}

