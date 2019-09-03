__global__ void populate_square(float *d_x,
    float *d_i1,
    int64_t m,
    int64_t n,
    int64_t b){
        // declare variables
        // gpu position variables
        int tx = threadIdx.x;
        int ty = threadIdx.y;
        int bx = blockIdx.x*blockDim.x + tx;
        int by = blockIdx.y*blockDim.y + ty;
        int ind = bx + m*by;
        float cnt;
        
        if(bx < m && by < n){
            cnt = 0;
                for(int ii = 0; ii < b; ii++){ // Histogram Algorithm
                    //if( powf((d_x[ii] - by)*(d_x[ii] - by) + (d_y[ii] - bx)*(d_y[ii] - bx) <= radius){
                    if((d_x[ii] - by) >= -0.5 && (d_x[ii] - by) < 0.5 && (d_x[ii+b] - bx) >= -0.5 && (d_x[ii+b] - bx) < 0.5){
                        cnt += 1; // Record found molecules
                    }
                } 
            d_i1[ind] = cnt; // global write to output variable
        }
    }