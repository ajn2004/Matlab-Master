/* density.h 
* This header contains code necessary to generate 3D histograms of localization data

 ajn 4/17/19
*/

__global__ void populate(float *d_x,
    float *d_y,
    float *d_z,
    float *d_i1,
    int64_t m,
    int64_t n,
    int64_t o,
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
        int ind = bx + m*by + m*n*bz;
        float cnt;
        
        if(bx < m && by < n && bz < o){
            cnt = 0;
                for(int ii = 0; ii < b; ii++){ // Histogram Algorithm
                    if( powf((d_x[ii] - by)*(d_x[ii] - by) + (d_y[ii] - bx)*(d_y[ii] - bx) + (d_z[ii] - bz)*(d_z[ii] - bz),0.5) <= radius){
                        cnt += 1; // Record found molecules
                    }
                } 
            d_i1[ind] = cnt; // global write to output variable
        }
    }

__global__ void grad(float *d_i1,
    float *d_ix,
    float *d_iy,
    float *d_iz,
    int64_t m,
    int64_t n,
    int64_t o){
        // declare variables
        // gpu position variables
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int tz = threadIdx.z;
        int bx = blockIdx.x*blockDim.x + tx;
        int by = blockIdx.y*blockDim.y + ty;
        int bz = blockIdx.z*blockDim.z + tz;
        int ind = bx + m*by + m*n*bz;     
        float xs[2] = {0, 0};
        if(bx < m , by < n , bz <o){
            // load neighboring x values or pad depending on position
            if(by - 1 >-1){xs[0] = d_i1[ind - m];}else{xs[0] = 0;}
            if(by + 1 < n){xs[1] = d_i1[ind + m];}else{xs[1] = 0;}
            d_ix[ind]=0.5*(xs[1] - xs[0]); // return differential

            if(bx - 1 >-1){xs[0] = d_i1[ind - 1];}else{xs[0] = 0;}
            if(bx + 1 < m){xs[1] = d_i1[ind + 1];}else{xs[1] = 0;}
            d_iy[ind]=0.5*(xs[1] - xs[0]); // return differential

            if(bz - 1 >-1){xs[0] = d_i1[ind - m*n];}else{xs[0] = 0;}
            if(bz + 1 < o){xs[1] = d_i1[ind + m*n];}else{xs[1] = 0;}
            d_iz[ind]=0.5*(xs[1] - xs[0]); // return differential
        }
    }

__global__ void laplace(float *d_i1,
    float *d_i2,
    int64_t m,
    int64_t n,
    int64_t o){
        // declare variables
        // gpu position variables
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int tz = threadIdx.z;
        int bx = blockIdx.x*blockDim.x + tx;
        int by = blockIdx.y*blockDim.y + ty;
        int bz = blockIdx.z*blockDim.z + tz;
        int ind = bx + m*by + m*n*bz;     
        float xs[2] = {0, 0};
        float ys[2] = {0, 0};
        float zs[2] = {0, 0};
        float mp;
        if(bx < m , by < n , bz <o){
            // load neighboring x values or pad depending on position
            mp = 2*d_i1[ind];
            if(by - 2 >-1){xs[0] = d_i1[ind - 2*m];}else{xs[0] = 0;}
            if(by + 2 < n){xs[1] = d_i1[ind + 2*m];}else{xs[1] = 0;}
            if(bx - 2 >-1){ys[0] = d_i1[ind - 2];}else{ys[0] = 0;}
            if(bx + 2 < m){ys[1] = d_i1[ind + 2];}else{ys[1] = 0;}
            if(bz - 2 >-1){zs[0] = d_i1[ind - 2*m*n];}else{zs[0] = 0;}
            if(bz + 2 < o){zs[1] = d_i1[ind + 2*m*n];}else{zs[1] = 0;}

            d_i2[ind]=0.25*((xs[1] + xs[0] - mp) + (ys[1] + ys[0] - mp) + (zs[1] + zs[0] - mp)); // return differential
        }
    }