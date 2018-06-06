/*
* image_neural.cu is a program to take matlab images and process them through a neural network to determine activation values for pixels in a 7x7 region centered on each pixel
*  We are assuming a single output for the neural network, and a single hidden layer of 100 nodes this will perform the necessary calculation on each pixel in the image and return
*	an image of 0 or 1 to be used later in the computation
*	V 1.0
*		we expect a format of [im_activate] = image_neural [i1, theta1, threta2, numoframes];
*	AJN 11/2/15
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
#define O_TILE_WIDTH 26								// variable to determine how many output tiles will be considered in a block
# define BLOCK_WIDTH (O_TILE_WIDTH + (7-1))		// block width needs to be output tiles + mask_width - 1 to ensure enough pixels are covered for calculation




/*
* Device code
*
*
*/



void __global__ activate(double *d_iall,   // the gaussian is a separable filter and be treated as such
	double *d_theta1,	// makes these elements eligible for constant caching
	double *d_theta2,
	double *d_ifin,
	int irow,
	int icol,
	int numi)
{
	// Declare variables
	//double d_i3[(50)] = { 0 };
	__shared__ double d_i2[(BLOCK_WIDTH)][(BLOCK_WIDTH)];		// preallocate space for shared image
	/*__shared__ double d_th1[(50)][(100)];
	__shared__ double d_th2[(100)];
*/
	// Coordinate building
	int tx = threadIdx.x;			// local x coord
	int ty = threadIdx.y;			// local y coord
	int tz = threadIdx.z;
	// location of output pixel being analyzed
	int row_output = blockIdx.y*O_TILE_WIDTH + ty;		// gives y coordinate as a function of tile width    **these lose meaning for (ty || tx) >= O_TILE_WIDTH and the same is true for **
	int col_output = blockIdx.x*O_TILE_WIDTH + tx;		// gives x coordinate as a function of tile width
	int imnum = blockIdx.z*blockDim.z + tz;
	if (imnum < numi){
	// initialize location of apron		this forces the first pixel to take care of both the first output pixel, and loading the first input pixel
	// BLOCK_WIDTH is larger than O_TILE_WIDTH so there are more threads being used than output pixels being calculated
	int row_input = row_output - 3;	// EACH thread should load 1 input tile to the shared image as there are [BLOCK_WIDTH]x[BLOCK_WIDTH] threads in a block
	int col_input = col_output - 3;	// and BLOCK_WIDTH = O_TILE_WIDTH + MASK_WIDTH-1
/*
	// Buffer data into block
	for (int tcol = 0; tcol < 100; tcol++){   // buffer theta1 matrix into shared block space
		for (int trow = 0; trow < 50; trow++){
			d_th1[trow][tcol] = d_theta1[trow + 50 * tcol];
		}
	}

	// Buffer theta2 into shared block space
	for (int t2row = 0; t2row < 101; t2row++){
		d_th2[t2row] = d_theta2[t2row];
	} */
	// buffer shared image into d_i2															
	// row/col_input represents the row/col of the input pixel being considered by 
	// thread [blockIdx.y*BLOCK_WIDTH+ty][blockIdx.x*BLOCK_WIDTH+tx]
	if ((row_input >= 0) && (row_input < irow) && (col_input >= 0) && (col_input < icol)){		// if statement checks the row/col indices to ensure they fall onto the input image
		d_i2[ty][tx] = d_iall[row_input + col_input*irow + imnum*irow*icol];										// if true, the value of the image is written to the shared array at location d_i2[ty][tx] and stored locally
	}																							// on the block
	else{
		d_i2[ty][tx] = 0.0*PI;																	// If row/col do not satisfy boundary condtions then assign a 0 to the value to build and apron of 
	}																							// of pixels that will not contribute to the calculation

	__syncthreads();																			// each thread uploads to a shared array later accessed by all threads, it is imperative to synch threads here
	//d_i3[0] = 1.0;
	// convolution calculation
	double z1[(100)] = { 0 };
	double a1[(101)] = { 0 };
	double a = 0.0;
	double z2 = 0.0;
	a1[0] = 1;
	if (ty < O_TILE_WIDTH && tx < O_TILE_WIDTH) {										// check that the local thread should be apart of the calcualtion
		/*for (int rowcount = 0; rowcount < 7; rowcount++){
			for (int colcount = 0; colcount < 7; colcount++){
				d_i3[rowcount + 7 * colcount + 1] = d_i2[rowcount + ty][colcount + tx];					// linearize region to prep for matrix math ensure coloumn major setup for proper neural net behavior
			} // end coloumn image for loop
		}*/
		// At this point d_i3 should be linearized as matlab would do so
		// perform matrix calculation here
		for (int th1count = 0; th1count < 100; th1count++){
			z1[th1count] = d_theta1[50*th1count];
			for (int rowcount = 0; rowcount < 7; rowcount++){
				for (int colcount = 0; colcount < 7; colcount++){
					z1[th1count] += d_i2[rowcount+ty][colcount + tx] * d_theta1[rowcount + 7*colcount + 1 + 50*th1count];
				} 
			}
			a1[th1count + 1] = powf(1.0 + exp(-z1[th1count]), -1.0);
		} // this completes the first half of the calculation

		for (int th2count = 0; th2count<101; th2count++){
			z2 += a1[th2count] * d_theta2[th2count];
		}
		a = powf(1 + exp(-z2), -1.0); // activation value
		//a = z2;


		if (row_output < irow && col_output < icol && imnum < numi){

			d_ifin[row_output + col_output*irow + imnum*irow*icol] = a;			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin

		}
		/*else{
		__syncthreads();
		d_ifin[row_output + col_output*irow + imnum*irow*icol] = d_ifin[row_output + col_output*irow];			// assign to output variable  THIS SECTION WILL CORRECTLY WRITE TO d_ifin
		__syncthreads();
		}  // end if else statement tyo decide what to write to final image
		*/
	}  // end if statement to decide whether to calculate output
	}
} // end gpu  void activate 




/*

THIS IS THE SECTION FOR IDENTIFICATION

*/





/*
* Host code
*
*
*/


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	double *iall;			// the pointer to the array of all images to be analyzed
	double *theta1;			// pointer to theta1 matrix
	double *theta2;			// pointer to theta2 matrix
	double  *d_iall;		// Pointer to image array on gpu
	double *d_theta1;		// Pointer to d_theta1 on gpu
	double *d_theta2;		// Pointer to d_theta2 on gpu
	double *d_ifin;			// pointer to d_ifin on gpu
	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;				// n
	int numi;				// number of images imported
	const int *idims, *th1dims, *th2dims;


	/* Throw an error if the input does not match expectations. */
	if (nrhs != 4) {
		printf("Must have 4 inputs ( i1, theta1, theta2, numthreads) line: %d\n", __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
		printf("i1 must be a n x m x numel(i1(1,1,:)) double array\n");
		mexErrMsgTxt("See Error above!\n");

	}
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
		printf("Theta1 must be a n +1 x l double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
		printf("Theta2 must be a l x 1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])){
		printf("number of threads per block must be an integer between 1 and 1024\n");
		mexErrMsgTxt("See Error above!\n");
	}


	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[0];
	irow = (int)idims[1];
	numi = mxGetScalar(prhs[3]);  // get number of images perblock from matlab


	// get theta1 dims
	theta1 = (double *)mxGetPr(prhs[1]);
	th1dims = mxGetDimensions(prhs[1]);
	int th1row = (int)th1dims[0]; // number of rows in theta1
	int th1col = (int)th1dims[1]; // number of coloumns in theta1

	// get theta2 dims
	theta2 = (double *)mxGetPr(prhs[2]);
	th2dims = mxGetDimensions(prhs[2]);
	int th2row = (int)th2dims[0]; // number of rows in theta2
	int th2col = 1; // number of coloumns in theta2


	// EVERYONE LOVES SOME GOOD VARIABLE CHECKING!!!!!!!!!!!!!!
	if (th1row != 50){
		printf("Theta1 must have 50 rows for what you want to do\n");
		mexErrMsgTxt("See Above Error!\n");
	}

	if (th1col != 100 || th2row != 101){
		printf("Theta2 must have oone more row than Theta1 has coloumns and that number should be 100\n");
		mexErrMsgTxt("See Above Error!\n");
	}

	if (th2col != 1){
		printf("Theta2 must have 1 coloumn\n");
		mexErrMsgTxt("See Above Error!\n");
	}

	// Did the User declare an output?
	if (nlhs != 1){
		printf("Declare an output variable [im_activate] = image_neural(i1, Theta1.', Theta2.')\n"); // oh user...... TEACH THEM A LESSON!!!!
		mexErrMsgTxt("See Error above!\n");
	}
	cudaDeviceReset();
	// allocate memory on the gpu device
	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*icol*(numi)*sizeof(double));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_theta1, th1row*th1col*sizeof(double));						// allocate theta1 memory
	if (err2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err3 = cudaMalloc((void**)&d_theta2, th2row*sizeof(double));						// allocate theta2 memory
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err4 = cudaMalloc((void**)&d_ifin, irow*icol*(numi)*sizeof(double));						// allocate completed activation image memory this will be a double for convience
	if (err4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	// copy data from host to device
	cudaError_t err6 = cudaMemcpy(d_iall, iall, irow*icol*(numi)*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err6 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err6), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err7 = cudaMemcpy(d_theta1, theta1, th1row*th1col*sizeof(double), cudaMemcpyHostToDevice);		// copy theta1 data to gpu
	if (err7 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err7), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err8 = cudaMemcpy(d_theta2, theta2, th2row*sizeof(double), cudaMemcpyHostToDevice);		// copy theta2 data to gpu
	if (err8 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err8), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
/*
	cudaError_t err10 = cudaMemcpy(d_ifin, iall, irow*icol*numi*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
*/
	/* Run GPU kernel*/
	dim3 dimBlock(BLOCK_WIDTH, BLOCK_WIDTH); // run 2-D gpu kernel to help with indexing
	dim3 dimGrid((icol - 1) / O_TILE_WIDTH + 1, (irow - 1) / O_TILE_WIDTH + 1, numi );
	



	
	activate << < dimGrid, dimBlock >> >(d_iall, d_theta1, d_theta2, d_ifin, irow, icol, (numi));

	//ident << < dimGrid, dimBlock >> >(d_ifin, d_iout, irow, icol, numi);

	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/

	cudaError_t errk1 = cudaPeekAtLastError();
	if (errk1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errk1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	
	cudaError_t errk2 = cudaThreadSynchronize();
	if (errk2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errk2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}


	plhs[0] = mxDuplicateArray(prhs[0]);
	double *iout = (double *)mxGetPr(plhs[0]);
	//printf("%d\n",numi);
	cudaError_t err9 = cudaMemcpy(iout, d_ifin, irow*icol*(numi)*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	//	printf("irow %f, icol %f, numi %f, line %d\n", ifin[0], ifin[1], ifin[2], __LINE__);

	cudaDeviceReset();
/*
	cudaFree(d_iall);
	cudaFree(d_theta1);
	cudaFree(d_theta2);
	cudaFree(d_ifin);
*/
	return;
}