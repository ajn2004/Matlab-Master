/* Density_hist.cu
* This algorithm will populate a 3-dimensional histogram with density information
* [Histogram ] = density_hist(x,y,z,bin_size,radius)
*/

#include <mex.h>
#include <stdint.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <helper_cuda.h>
#include <math.h>
#define SHARE_PIX 4096
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

__global__ void populate(float *d_x,
    float *d_y,
    float *d_z,
    float *d_i1,
    int64_t *imsize,
    int64_t b,
    float radius){
        // declare variables
        // gpu position variables
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int tz = threadIdx.z;
        int bx = blockIdx.x*blockDim.x + tx;
        int by = blockIdx.y*blockDim.y + ty;
        int bz = blockIdx.z*blockDim.z + tz;
        __shared__ float xf[SHARE_PIX], yf[SHARE_PIX], zf[SHARE_PIX];  // shared logic, we get 4.1k spots of 'shared' memory
                                                        // we will use those spots to cycle over all localizations
                                                        // we can speed this up by using the parallel threads to grab data
        int crunch = blockDim.x*blockDim.y*blockDim.z;
        //if(crunch > SHARE_PIX){crunch = SHARE_PIX;}
        int dthr = tx + blockDim.x*ty + blockDim.x*blockDim.y*tz; // local thread for moving data into local
        int chunks = ((int)b-1)/crunch +1; // number of data chunks to 
        // here ind is column major so in terms of i,j,k and m,n,o being index and size respectively
        // ind = i + j*m + k*m*n
        int m = imsize[0];
        int n = imsize[1];
        int o = imsize[2];
        int ind = bx + m*by + m*n*bz;
        //int i = ind % m;
        //int j = (ind / m) % n;
        //int k = ind/(m*n);
        float cnt, dist;
        int ld, lk;
        int toplp = crunch;

        if(bx < m && by < n && bz < o){
            cnt = 0;
           //for(lk = 0; lk < chunks; lk++){
                //ld = dthr + crunch*lk; // load variable (linearized thread ID + loopvar * max_share_size)
               /* __syncthreads(); // ensure all data is loaded
                
                    if(ld < (int)b){ // if load variable has something to load, load the localization data
                        xf[ld % SHARE_PIX] = d_x[ld] - (float)by;
                        yf[ld % SHARE_PIX] = d_y[ld] - (float)bx;
                        zf[ld % SHARE_PIX] = d_z[ld] - (float)bz;
                    }
                
                __syncthreads(); // ensure all data is loaded*/
                //if(crunch*(lk+1) > (int)b-crunch*(lk)){toplp = (int)b - crunch*lk;} // if we do not utitlize full data space, ensure we only count 'new' data
                //if(toplp >0){
                    for(int ii = 0; ii < b; ii++){ // Histogram Algorithm
                        //dist = (d_x[ii] - (float)by)*(d_x[ii] - (float)by) + (d_y[ii] - (float)bx)*(d_y[ii] - (float)bx) + (d_z[ii] - (float)bz)*(d_z[ii] - (float)bz);
                        // This commented secion is 'slow' but gives correct math
                        /*dist += (d_x[ii] - (float)by)*(d_x[ii] - (float)by);    
                        dist += (d_y[ii] - (float)bx)*(d_y[ii] - (float)bx);                                    
                        dist += (d_z[ii] - (float)bz)*(d_z[ii] - (float)bz);*/
                        //if( powf(dist,0.5) <= radius){ // check if localization is within 'radius'
                        if( powf((d_x[ii] - (float)by)*(d_x[ii] - (float)by) + (d_y[ii] - (float)bx)*(d_y[ii] - (float)bx) + (d_z[ii] - (float)bz)*(d_z[ii] - (float)bz),0.5) <= radius){
                        //if( powf(xf[ii]*xf[ii] + yf[ii]*yf[ii] + zf[ii]*zf[ii],0.5) <= radius){ // check if localization is within 'radius'
                            cnt += 1; // Record found molecules
                        }
                    } 
             //   }
           // }
            d_i1[ind] = cnt; // global write to output variable
        }
    }

void mexFunction(int64_t nlhs, mxArray *plhs[],
    int64_t nrhs, mxArray const *prhs[])
    {
        // Declare Variables for use
        cudaEvent_t start, stop;
        float time;
        
        float *xf = (float *)mxGetPr(prhs[0]);
        float *yf = (float *)mxGetPr(prhs[1]);
        float *zf = (float *)mxGetPr(prhs[2]);
        float radius = (float)mxGetScalar(prhs[4]);
        float *d_i1, *d_x, *d_y, *d_z;
        
        int64_t *imsize, *d_size;
        int64_t *bigly = (int64_t *)mxGetDimensions(prhs[0]);
        
        // imsize is an int64_t array containing the size of the histogram to be constructed
        imsize = (int64_t *)mxGetPr(prhs[3]);
        if(imsize[2] < 0 || imsize[2] > 10000){imsize[2] = 1;} // fix glitch that occurs for 2D grids
        printf("m = %d, n = %d, o= %d\n",imsize[0],imsize[1],imsize[2]);
        int64_t pixels = imsize[0]*imsize[1]*imsize[2];
        printf("pixels = %d\n",pixels);
        // Allocate Space
        

        CUDA_SAFE_CALL(cudaMalloc((void**)&d_x,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_y,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_z,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_i1,sizeof(float)*pixels));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_size,sizeof(int64_t)*3));
        
        
        // Copy Data onto GPU
        
        CUDA_SAFE_CALL(cudaMemcpy(d_x,xf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_y,yf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_z,zf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_size,imsize,sizeof(int64_t)*3,cudaMemcpyHostToDevice));
        
        // GPU Setup and Launch
        int tpx = 16;
        int tpy = 8;  // This gives largest multiple of 32 (warp size) for maximum occupancy
        int tpz = 8;
        dim3 dimGrid(((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1,((int)imsize[2]-1)/tpz + 1) ;
        dim3 dimBlock(tpx, tpy, tpz);
        
        printf("Dimensions of the grid %d x %d x %d\n", ((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1,((int)imsize[2]-1)/tpz + 1);
        
        // Launch GPU with timing
        

        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start,0);
        populate <<<dimGrid,dimBlock>>> (d_x, d_y, d_z, d_i1, d_size, bigly[0], radius);
        CUDA_CHECK_ERROR();
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time,start ,stop);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
        printf("Kernel Execution in %f s\n",time/1000);
        // Gather Data Back
        size_t *dims;
        dims[0] = (size_t)imsize[0];
        dims[1] = imsize[1];
        dims[2] = imsize[2];
        plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
        float *pout = (float *)mxGetPr(plhs[0]);
        
        CUDA_SAFE_CALL(cudaMemcpy(pout,d_i1,sizeof(float)*imsize[0]*imsize[1]*imsize[2],cudaMemcpyDeviceToHost));
        
        // Clean up
        
        cudaFree(d_size);
        cudaFree(d_x);
        cudaFree(d_y);
        cudaFree(d_i1);
        cudaFree(d_z);
        
    }