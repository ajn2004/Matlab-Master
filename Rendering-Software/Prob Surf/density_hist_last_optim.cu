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

//#define O_TILE_WIDTH 20								// variable to determine how many output tiles will be considered in a block
//# define BLOCK_WIDTH (O_TILE_WIDTH + (pix-1))		// block width needs to be output tiles + mask_width - 1 to ensure enough pixels are covered for calculation



__global__ void populate(float *d_x,
    float *d_y,
    float *d_z,
    int64_t *d_i1,
    int64_t *imsize,
    int64_t b,
    float radius){
        // declare variables
        // gpu position variables
        int64_t tx = threadIdx.x;
        int64_t ty = threadIdx.y;
        int64_t tz = threadIdx.z;
        int64_t bx = blockIdx.x*blockDim.x + tx;
        int64_t by = blockIdx.y*blockDim.y + ty;
        int64_t bz = blockIdx.z*blockDim.z + tz;
        __shared__ float xf[4096], yf[4096], zf[4096];  // shared logic, we get 4.1k spots of 'shared' memory
                                                        // we will use those spots to cycle over all localizations
                                                        // we can speed this up by using the parallel threads to grab data
        int64_t crunch = blockDim.x*blockDim.y*blockDim.z;
        int64_t dthr = tx + blockDim.x*ty+blockDim.x*blockDim.y*tz; // local thread for moving data into local
        int64_t chunks = (b-1)/crunch +1; // number of data chunks to 
        // here ind is column major so in terms of i,j,k and m,n,o being index and size respectively
        // ind = i + j*m + k*m*n
        int64_t m = imsize[0];
        int64_t n = imsize[1];
        int64_t o = imsize[2];
        int64_t ind = bx + m*by + m*n*bz;
        int64_t i = ind % m;
        int64_t j = (ind / m) % n;
        int64_t k = ind/(m*n);
        float dist;
        int ld, cnt;
        int64_t toplp = crunch;

        
        cnt = 0;
        for(int lk = 0; lk < chunks; lk++){
            ld = dthr + crunch*lk; // load variable (linearized thread ID + loopvar * max_share_size)

            if(ld < b){ // if load variable has something to load, load the localization data
                xf[ld] = d_x[ld];
                yf[ld] = d_y[ld];
                zf[ld] = d_z[ld];
            }
            __syncthreads(); // ensure all data is loaded
            if(crunch*(lk+1) > b-crunch*(lk)){toplp = b - crunch*lk;} // if we do not utitlize full data space, ensure we only count 'new' data
            for(int ii = 0; ii < toplp; ii++){ // Histogram Algorithm
                dist = 0; // add distances
                dist += (xf[ii] - (float)j)*(xf[ii] - (float)j);                    
                dist += (yf[ii] - (float)i)*(yf[ii] - (float)i);                    
                dist += (zf[ii] - (float)k)*(zf[ii] - (float)k);
                /* This commented secion is 'slow' but gives correct math
                dist += (d_x[ii] - (float)j)*(d_x[ii] - (float)j);    
                dist += (d_y[ii] - (float)i)*(d_y[ii] - (float)i);                                    
                dist += (d_z[ii] - (float)k)*(d_z[ii] - (float)k);*/
                if( powf(dist,0.5) <= radius){ // check if localization is within 'radius'
                    cnt += 1; // Record found molecules
                }
            } 
        }
        if(bx < m && by < n && bz < o){
            d_i1[ind] = cnt; // global write to output variable
        }
    }

void mexFunction(int64_t nlhs, mxArray *plhs[],
    int64_t nrhs, mxArray const *prhs[])
    {
        // Declare Variables for use
        float *xf = (float *)mxGetPr(prhs[0]);
        float *yf = (float *)mxGetPr(prhs[1]);
        float *zf = (float *)mxGetPr(prhs[2]);
        float *d_x, *d_y, *d_z;
        
        int64_t *d_i1;
        float radius = (float)mxGetScalar(prhs[4]);
        int64_t *bigly = (int64_t *)mxGetDimensions(prhs[0]);
        int64_t *imsize, *d_size;
        
        // imsize is an int64_t array containing the size of the final histogram
        imsize = (int64_t *)mxGetPr(prhs[3]);
        if(imsize[2] < 0 || imsize[2] > 10000){imsize[2] = 1;}
        printf("m = %d, n = %d, o= %d\n",imsize[0],imsize[1],imsize[2]);
        int64_t pixels = imsize[0]*imsize[1]*imsize[2];
        printf("pixels = %d\n",pixels);
        // Allocate Space
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_x,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_y,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_z,sizeof(float)*bigly[0]));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_i1,sizeof(int64_t)*pixels));
        CUDA_SAFE_CALL(cudaMalloc((void**)&d_size,sizeof(int64_t)*3));

        // Copy Data onto GPU
        CUDA_SAFE_CALL(cudaMemcpy(d_x,xf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_y,yf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_z,zf,sizeof(float)*bigly[0],cudaMemcpyHostToDevice));
        CUDA_SAFE_CALL(cudaMemcpy(d_size,imsize,sizeof(int64_t)*3,cudaMemcpyHostToDevice));
      
        // GPU Setup and Launch
        int tpx = 15;
        int tpy = 8;
        int tpz = 8;
        dim3 dimGrid(((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1,((int)imsize[2]-1)/tpz + 1) ;
        dim3 dimBlock(tpx, tpy, tpz);
        printf("Dimensions of the grid %d x %d x %d\n", ((int)imsize[0]-1)/tpx + 1,((int)imsize[1]-1)/tpy + 1,((int)imsize[2]-1)/tpz + 1);
        
        // Launch GPU
        populate <<<dimGrid,dimBlock>>> (d_x, d_y, d_z, d_i1, d_size, bigly[0], radius);
        CUDA_CHECK_ERROR();
        
        // Gather Data Back
        size_t *dims;
        dims[0] = (size_t)imsize[0];
        dims[1] = imsize[1];
        dims[2] = imsize[2];
        plhs[0] = mxCreateNumericArray(3, dims, mxINT64_CLASS, mxREAL);
        int64_t *pout = (int64_t *)mxGetPr(plhs[0]);

        //printf("Maybe?");
        CUDA_SAFE_CALL(cudaMemcpy(pout,d_i1,sizeof(int64_t)*imsize[0]*imsize[1]*imsize[2],cudaMemcpyDeviceToHost));

        //printf("You don't see this");
        cudaFree(d_x);
        cudaFree(d_y);
        cudaFree(d_i1);
        cudaFree(d_z);
    }