/* Density_hist.cu
* This algorithm will populate a 3-dimensional histogram with density information
* [Histogram ] = density_hist(x,y,z,bin_size,radius)
*/

#include <mex.h>
#include <stdint.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <helper_cuda.h>
#include "make_pics.h"
#include <math.h>
#define CUDA_CHECK_ERROR() __cuda_check_errors(__FILE__, __LINE__)
#define CUDA_SAFE_CALL(err) __cuda_safe_call(err, __FILE__, __LINE__)
inline void
__cuda_check_errors (const char *filename, const int64_t line_number)
{
  cudaError err = cudaDeviceSynchronize ();
  if (err != cudaSuccess)
    {
      printf ("CUDA error %i at %s:%i: %s\n",
          err, filename, line_number, cudaGetErrorString (err));
      cudaDeviceReset();
      mexErrMsgTxt("Get Thee To A Nunnery!\n");
    }
}

inline void
__cuda_safe_call (cudaError err, const char *filename, const int64_t line_number)
{
  if (err != cudaSuccess)
    {
      printf ("CUDA error %i at %s:%i: %s\n",
          err, filename, line_number, cudaGetErrorString (err));
          cudaDeviceReset();
          mexErrMsgTxt("Get Thee To A Nunnery!\n");
    }
}



void mexFunction(int64_t nlhs, mxArray *plhs[],
    int64_t nrhs, mxArray const *prhs[])
    {
        // Declare Variables for use
        cudaEvent_t start, stop;
        float time;
        
        float *xf = (float *)mxGetPr(prhs[0]);
        float *d_i1, *d_x;
        
        int64_t *imsize;
        int64_t *bigly = (int64_t *)mxGetDimensions(prhs[0]);
        
        // imsize is an int64_t array containing the size of the histogram to be constructed
        imsize = (int64_t *)mxGetPr(prhs[1]);
        //if(imsize[2] < 0 || imsize[2] > 10000){imsize[2] = 1;} // fix glitch that occurs for 2D grids
        printf("m = %d, n = %d\n",imsize[0],imsize[1]);
        int64_t pixels = imsize[0]*imsize[1];
        printf("pixels = %d\n",pixels);
        // Allocate Space
        

        CUDA_SAFE_CALL(cudaMalloc((void**)&d_x,sizeof(float)*2*bigly[0]));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_y,sizeof(float)*bigly[0]));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_z,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_i1,sizeof(float)*pixels));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_i2,sizeof(float)*pixels));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_iy,sizeof(float)*pixels));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_iz,sizeof(float)*pixels));
        //CUDA_SAFE_CALL(cudaMalloc((void**)&d_size,sizeof(int64_t)*3));
        
        
        // Copy Data onto GPU
        
        CUDA_SAFE_CALL(cudaMemcpy(d_x,xf,sizeof(float)*2*bigly[0],cudaMemcpyHostToDevice));
        //CUDA_SAFE_CALL(cudaMemcpy(d_y,yf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        //CUDA_SAFE_CALL(cudaMemcpy(d_z,zf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        //CUDA_SAFE_CALL(cudaMemcpy(d_size,imsize,sizeof(int64_t)*3,cudaMemcpyHostToDevice));
        
        // GPU Setup and Launch
        int tpx = 32;
        int tpy = 32;  // This gives largest multiple of 32 (warp size) for maximum occupancy
        dim3 dimGrid(((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1, 1) ;
        dim3 dimBlock(tpx, tpy, 1);
        
        printf("Dimensions of the grid %d x %d x %d\n", ((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1,1);
        
        // Launch GPU with timing
        

        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start,0);
        populate_square <<<dimGrid,dimBlock>>> (d_x, d_i1, imsize[0], imsize[1], bigly[0]);
        CUDA_CHECK_ERROR();
        
        
        //CUDA_CHECK_ERROR();
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time,start ,stop);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
        printf("Kernel Execution in %f s\n",time/1000);
        // Gather Data Back
        const size_t dims[3]={imsize[0], imsize[1]};
        //cdims = (size_t)imsize;
        
        plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
        
        //printf("Requesting %d x %d x %d values\n", dims[0], dims[1], dims[2]);
        //plhs[1] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
        

        //plhs[2] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
        printf("Requesting %d x %d values\n", dims[0], dims[1]);
        float *o1 = (float *)mxGetPr(plhs[0]);
        //float *o2 = (float *)mxGetPr(plhs[1]);
        //float *o3 = (float *)mxGetPr(plhs[2]);
        
        
        CUDA_SAFE_CALL(cudaMemcpy(o1,d_i1,sizeof(float)*pixels,cudaMemcpyDeviceToHost));
        //CUDA_SAFE_CALL(cudaMemcpy(o2,d_iy,sizeof(float)*imsize[0]*imsize[1]*imsize[2],cudaMemcpyDeviceToHost));
        //CUDA_SAFE_CALL(cudaMemcpy(o3,d_iz,sizeof(float)*imsize[0]*imsize[1]*imsize[2],cudaMemcpyDeviceToHost));

        // Clean up
        
        //cudaFree(d_size);

        cudaFree(d_i1);
        cudaFree(d_x);       
    }