/*
 * gpu_tol is a program that will perform localization tolerance
 * to molecules in parallel
 * 
 * This function will be called by matlab as
 * [N, Nc, O, Oc, Sx, Sxc, Sy, Syc, X, Xc, Y, Yc, llv, Fn, lp] = gpu_tol(N, Nc, O, Oc, Sx, Sxc, Sy, Syc, X, Xc, Y, Yc, llv, Fn, lp, Nl, Nh, Ncl, Nch, Ol, Oh, Xcl, Xch, Ycl, Ych, Sl, Sh, Scl, Sch, fNh, fOh, fSh, ll, lh)
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

// gpu function to build fractional error variables
void __global__ fractionate(double *out, double *x, double *xc, int m) {
	int index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index < m) {
		if (xc[index] >= 0) {
			out[index] = sqrt(xc[index]) / x[index];
		}
		else {
			out[index] = 100.0;
		}
	}
}

void __global__ tolerate(double *d_Ni,
	double *d_No,
	double *d_Nci,
	double *d_Nco,
	double *d_Oi,
	double *d_Oo,
	double *d_Oci,
	double *d_Oco,
	double *d_Sxi,
	double *d_Sxo,
	double *d_Sxci,
	double *d_Sxco,
	double *d_Syi,
	double *d_Syo,
	double *d_Syci,
	double *d_Syco,
	double *d_Xi,
	double *d_Xo,
	double *d_Xci,
	double *d_Xco,
	double *d_Yi,
	double *d_Yo,
	double *d_Yci,
	double *d_Yco,
	double *d_llvi,
	double *d_llvo,
	double *d_fni,
	double *d_fno,
	double *d_lpi,
	double *d_lpo,
	double *d_frN,
	double *d_frO,
	double *d_frSx,
	double *d_frSy,
	double Nl,
	double Nh,
	double Ncl,
	double Nch,
	double Ol,
	double Oh,
	double Ocl,
	double Och,
	double Xcl,
	double Xch,
	double Ycl,
	double Ych,
	double Sl,
	double Sh,
	double Scl,
	double Sch,
	double fNh,
	double fOh,
	double fSh,
	double ll,
	double lh,
	int m) {
	int index = blockDim.x*blockIdx.x + threadIdx.x;
	if(index < m) 
	{ // verify you are working on a molecule
		if(d_Ni[index] >= Nl && d_Ni[index] <= Nh && d_Nci[index] >= Ncl && d_Nci[index] <= Nch && d_Oi[index] >= Ol && d_Oi[index] <= Oh && d_Oci[index] >= Ocl && d_Oci[index] <= Och && d_Xci[index] >= Xcl && d_Xci[index] <= Xch && d_Yci[index] >= Ycl && d_Yci[index] <= Ych && d_Sxi[index] >= Sl && d_Sxi[index] <= Sh && d_Syi[index] >= Sl && d_Syi[index] <= Sh && d_Sxci[index] >= Scl && d_Sxci[index] <= Sch && d_Syci[index] >= Scl && d_Syci[index] <= Sch && d_frN[index] <= fNh && d_frO[index] <= fOh && d_frSx[index] <= fSh && d_frSy[index] <= fSh && d_llvi[index] >= ll && d_llvi[index] <= lh)
		{    d_No[index] = d_Ni[index];
			d_Nco[index] = d_Nci[index];
			 d_Oo[index] = d_Oi[index];
			d_Oco[index] = d_Oci[index];
		    d_Sxo[index] = d_Sxi[index];
			d_Syo[index] = d_Syi[index];
		   d_Sxco[index] = d_Sxci[index];
		   d_Syco[index] = d_Syci[index];
			 d_Xo[index] = d_Xi[index];
		    d_Xco[index] = d_Xci[index];
		     d_Yo[index] = d_Yi[index];
		    d_Yco[index] = d_Yci[index];
		   d_llvo[index] = d_llvi[index];
			d_fno[index] = d_fni[index];
			d_lpo[index] = d_lpi[index];
		} // end tolerance if
		else {
			d_No[index] = -1;
			d_Nco[index] = -1;
			d_Oo[index] = -1;
			d_Oco[index] = -1;
			d_Sxo[index] = -1;
			d_Syo[index] = -1;
			d_Sxco[index] = -1;
			d_Syco[index] = -1;
			d_Xo[index] = -1;
			d_Xco[index] = -1;
			d_Yo[index] = -1;
			d_Yco[index] = -1;
			d_llvo[index] = -1;
			d_fno[index] = -1;
			d_lpo[index] = -1;
		}// end else
	}// end working on a molecule
}//end global


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
	double *Ni, *No, *Nci, *Nco, *Oi, *Oo, *Oci, *Oco, *Xi, *Xo, *Xci, *Xco;
	double *Yi, *Yo, *Yci, *Yco, *Sxi, *Sxo, *Sxci, *Sxco, *Syi, *Syo, *Syci, *Syco;
	double *llvi, *llvo, *lpi, *lpo, *fni, *fno;
	// GPU Variables
	double *d_Ni, *d_No, *d_Nci, *d_Nco, *d_Oi, *d_Oo, *d_Oci, *d_Oco, *d_Xi, *d_Xo, *d_Xci, *d_Xco;
	double *d_Yi, *d_Yo, *d_Yci, *d_Yco, *d_Sxi, *d_Sxo, *d_Sxci, *d_Sxco, *d_Syi, *d_Syo, *d_Syci, *d_Syco;
	double *d_llvi, *d_llvo, *d_lpi, *d_lpo, *d_fni, *d_fno, *d_frN, *d_frO, *d_frSx, *d_frSy;

	// single entry doubles
	double Nl, Nh, Ncl, Nch, Oh, Ol, Och, Ocl, Xch, Xcl, Ych, Ycl;
	double Sh, Sl, Sch, Scl, fNh, fOh, fSh, ll, lh;

	// Error Message Array

	// Vector Size Array
	if (nrhs != 36) {
		mexErrMsgTxt("You need 36 input variables");
	}
	if (nlhs != 15) {
		mexErrMsgTxt("You need 15 output variables");
	}
	if (mxGetM(prhs[0]) != mxGetM(prhs[1]) || mxGetM(prhs[0]) != mxGetM(prhs[2]) || mxGetM(prhs[0]) != mxGetM(prhs[3])) {
		mexErrMsgTxt("Your input vectors must be the same size!\n");
	}
	if (mxGetM(prhs[0]) != mxGetM(prhs[4]) || mxGetM(prhs[0]) != mxGetM(prhs[5]) || mxGetM(prhs[0]) != mxGetM(prhs[6])) {
		mexErrMsgTxt("Your input vectors must be the same size!\n");
	}
	if (mxGetM(prhs[0]) != mxGetM(prhs[7]) || mxGetM(prhs[0]) != mxGetM(prhs[8]) || mxGetM(prhs[0]) != mxGetM(prhs[9])) {
		mexErrMsgTxt("Your input vectors must be the same size!\n");
	}
	if (mxGetM(prhs[0]) != mxGetM(prhs[10]) || mxGetM(prhs[0]) != mxGetM(prhs[11]) || mxGetM(prhs[0]) != mxGetM(prhs[12])) {
		mexErrMsgTxt("Your input vectors must be the same size!\n");
	}
	if (mxGetM(prhs[0]) != mxGetM(prhs[13]) || mxGetM(prhs[0]) != mxGetM(prhs[14])) {
		mexErrMsgTxt("Your input vectors must be the same size!\n");
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
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[10]) || mxIsComplex(prhs[10])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[11]) || mxIsComplex(prhs[11])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[12]) || mxIsComplex(prhs[12])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[13]) || mxIsComplex(prhs[13])) {
		mexErrMsgTxt("Your input vectors must contain doubles\n");
	}
	if (!mxIsDouble(prhs[14]) || mxIsComplex(prhs[14])) {
		mexErrMsgTxt("Your input vectors must contain real doubles\n");
	}
	
	// Grab dimension size of data
	const size_t *dims;
	dims = mxGetDimensions(prhs[0]);
	int m = (int)dims[0];
	int n = (int)dims[1];
	const int mem_size = m*n*sizeof(double);

	// Get position of Data
	Ni = (double *)mxGetPr(prhs[0]);
	Nci = (double *)mxGetPr(prhs[1]);
	Oi = (double *)mxGetPr(prhs[2]);
	Oci = (double *)mxGetPr(prhs[3]);
	Sxi = (double *)mxGetPr(prhs[4]);
	Sxci = (double *)mxGetPr(prhs[5]);
	Syi = (double *)mxGetPr(prhs[6]);
	Syci = (double *)mxGetPr(prhs[7]);
	Xi = (double *)mxGetPr(prhs[8]);
	Xci = (double *)mxGetPr(prhs[9]);
	Yi = (double *)mxGetPr(prhs[10]);
	Yci = (double *)mxGetPr(prhs[11]);
	llvi = (double *)mxGetPr(prhs[12]);
	fni = (double *)mxGetPr(prhs[13]);
	lpi = (double *)mxGetPr(prhs[14]);
	
	// Get Tolerance Limits
	Nl = mxGetScalar(prhs[15]);
	Nh = mxGetScalar(prhs[16]);
	Ncl = mxGetScalar(prhs[17]);
	Nch = mxGetScalar(prhs[18]);

	Ol = mxGetScalar(prhs[19]);
	Oh = mxGetScalar(prhs[20]);
	Ocl = mxGetScalar(prhs[21]);
	Och = mxGetScalar(prhs[22]);

	Xcl = mxGetScalar(prhs[23]);
	Xch = mxGetScalar(prhs[24]);
	Ycl = mxGetScalar(prhs[25]);
	Ych = mxGetScalar(prhs[26]);

	Sl = mxGetScalar(prhs[27]);
	Sh = mxGetScalar(prhs[28]);
	Scl = mxGetScalar(prhs[29]);
	Sch = mxGetScalar(prhs[30]);

	fNh = mxGetScalar(prhs[31]);
	fOh = mxGetScalar(prhs[32]);
	fSh = mxGetScalar(prhs[33]);
	ll = mxGetScalar(prhs[34]);
	lh = mxGetScalar(prhs[35]);

	
	// Fairly Certain at this point all data is accessible through the program

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
	checkCudaErrors(cudaMalloc((void**)&d_Ni, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_No, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Nci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Nco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_Oi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Oo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Oci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Oco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_Sxi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Sxo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Sxci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Sxco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_Syi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Syo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Syci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Syco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_Xi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Xo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Xci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Xco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_Yi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Yo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Yci, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_Yco, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_llvi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_llvo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_fni, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_fno, mem_size));

	checkCudaErrors(cudaMalloc((void**)&d_lpi, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_lpo, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_frN, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_frO, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_frSx, mem_size));
	checkCudaErrors(cudaMalloc((void**)&d_frSy, mem_size));
	
	// Data Copy Array
	checkCudaErrors(cudaMemcpy(d_Ni, Ni, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Nci, Nci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Oi, Oi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Oci, Oci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Sxi, Sxi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Sxci, Sxci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Syi, Syi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Syci, Syci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Xi, Xi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Xci, Xci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Yi, Yi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Yci, Yci, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_llvi, llvi, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_fni, fni, mem_size, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_lpi, lpi, mem_size, cudaMemcpyHostToDevice));
	
	// Get Fractional error Vectors and tolerance
	fractionate << <((m - 1) / 1024 + 1), 1024 >> > (d_frN, d_Ni, d_Nci, m);
	fractionate << <((m - 1) / 1024 + 1), 1024 >> > (d_frO, d_Oi, d_Oci, m);
	fractionate << <((m - 1) / 1024 + 1), 1024 >> > (d_frSx, d_Sxi, d_Sxci, m);
	fractionate << <((m - 1) / 1024 + 1), 1024 >> > (d_frSy, d_Syi, d_Syci, m);

	tolerate << <((m - 1) / 1024 + 1), 1024 >> > (d_Ni, d_No, d_Nci, d_Nco, d_Oi, d_Oo, d_Oci, d_Oco, d_Sxi, d_Sxo, d_Sxci, d_Sxco, d_Syi, d_Syo, d_Syci, d_Syco, d_Xi, d_Xo, d_Xci, d_Xco, d_Yi, d_Yo, d_Yci, d_Yco, d_llvi, d_llvo, d_fni, d_fno, d_lpi, d_lpo, d_frN, d_frO, d_frSx, d_frSy, Nl, Nh, Ncl, Nch, Ol, Oh, Ocl, Och, Xcl, Xch, Ycl, Ych, Sl, Sh, Scl, Sch, fNh, fOh, fSh, ll, lh, m);


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
	 plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[3] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[4] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[5] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[6] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[7] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[8] = mxCreateDoubleMatrix(m, n, mxREAL);
	 plhs[9] = mxCreateDoubleMatrix(m, n, mxREAL);
	plhs[10] = mxCreateDoubleMatrix(m, n, mxREAL);
	plhs[11] = mxCreateDoubleMatrix(m, n, mxREAL);
	plhs[12] = mxCreateDoubleMatrix(m, n, mxREAL);
	plhs[13] = mxCreateDoubleMatrix(m, n, mxREAL);
	plhs[14] = mxCreateDoubleMatrix(m, n, mxREAL);

	  No = (double *)mxGetPr(plhs[0]);
	 Nco = (double *)mxGetPr(plhs[1]);
	  Oo = (double *)mxGetPr(plhs[2]);
	 Oco = (double *)mxGetPr(plhs[3]);
	 Sxo = (double *)mxGetPr(plhs[4]);
	Sxco = (double *)mxGetPr(plhs[5]);
	 Syo = (double *)mxGetPr(plhs[6]);
	Syco = (double *)mxGetPr(plhs[7]);
	  Xo = (double *)mxGetPr(plhs[8]);
	 Xco = (double *)mxGetPr(plhs[9]);
	  Yo = (double *)mxGetPr(plhs[10]);
	 Yco = (double *)mxGetPr(plhs[11]);
	llvo = (double *)mxGetPr(plhs[12]);
	 fno = (double *)mxGetPr(plhs[13]);
	 lpo = (double *)mxGetPr(plhs[14]);



	// copy data array
	 checkCudaErrors(cudaMemcpy(No, d_No, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Nco, d_Nco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Oo, d_Oo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Oco, d_Oco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Sxo, d_Sxo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Sxco, d_Sxco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Syo, d_Syo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Syco, d_Syco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Xo, d_Xo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Xco, d_Xco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Yo, d_Yo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(Yco, d_Yco, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(llvo, d_llvo, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(fno, d_fno, mem_size, cudaMemcpyDeviceToHost));
	 checkCudaErrors(cudaMemcpy(lpo, d_lpo, mem_size, cudaMemcpyDeviceToHost));
		
	/*


	
	
	checkCudaErrors(cudaMemcpy(y_out,d_data,mem_size, cudaMemcpyDeviceToHost));
	cufftDestroy(plan);
	cudaFree(d_data);
	// create complex double in matlab
	plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);


	yor = mxGetPr(plhs[0]);
	yoi = mxGetPi(plhs[0]);
	unpack_c2c(y_out, yor, yoi, n*m);
	mxFree(data);
	*/

	// Release GPU memory
	cudaFree(d_Ni);
	cudaFree(d_No);
	cudaFree(d_Nci);
	cudaFree(d_Nco);
	
	cudaFree(d_Oi);
	cudaFree(d_Oo);
	cudaFree(d_Oci);
	cudaFree(d_Oco);
	
	cudaFree(d_Sxi);
	cudaFree(d_Sxo);
	cudaFree(d_Sxci);
	cudaFree(d_Sxco);

	cudaFree(d_Syi);
	cudaFree(d_Syo);
	cudaFree(d_Syci);
	cudaFree(d_Syco);

	cudaFree(d_Xi);
	cudaFree(d_Xo);
	cudaFree(d_Xci);
	cudaFree(d_Xco);

	cudaFree(d_Yi);
	cudaFree(d_Yo);
	cudaFree(d_Yci);
	cudaFree(d_Yco);

	cudaFree(d_llvi);
	cudaFree(d_llvo);
	cudaFree(d_fni);
	cudaFree(d_fno);

	cudaFree(d_lpi);
	cudaFree(d_lpo);
	cudaFree(d_frN);
	cudaFree(d_frO);
	cudaFree(d_frSx);
	cudaFree(d_frSy);
}

