/*
* 3D chain_loc v 1.0
* This rendition is to build off of the chain loc code, but follow the derivation in Smith et al 2010 nat meth
*/

/*
Expected input and output
[xf_all, yf_all, zf_all, N, off_all, xf_crlb, yf_crlb, zf_crlb, N_crlb, off_crlb, llv_all] = 3d_chain_loc(ilocs, zcurve, numthreads, angle)
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

#define PI 3.14159265358979323846



/*
* Device code
*
* To facilitate coding (for me) I have copied the localization algorithm to be used with multiple sized areas
*/

/*

Device Functions

*/
__device__ float device_det(float fisher[25])  // updated determinant to 5 x 5 as of 2-6-18
{
	float det;
	det = fisher[0] * (fisher[6] * (fisher[12] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) + fisher[22] * (fisher[13] * fisher[19] - fisher[18] * fisher[14])) - fisher[11] * (fisher[7] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[19] - fisher[18] * fisher[9])) + fisher[16] * (fisher[7] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[14] - fisher[13] * fisher[9])) - fisher[21] * (fisher[7] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) + fisher[17] * (fisher[8] * fisher[14] - fisher[13] * fisher[9]))) - fisher[5] * (fisher[1] * (fisher[12] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) + fisher[22] * (fisher[13] * fisher[19] - fisher[18] * fisher[14])) - fisher[11] * (fisher[2] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[19] - fisher[18] * fisher[4])) + fisher[16] * (fisher[2] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[14] - fisher[13] * fisher[4])) - fisher[21] * (fisher[2] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[17] * (fisher[3] * fisher[14] - fisher[13] * fisher[4]))) + fisher[10] * (fisher[1] * (fisher[7] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[19] - fisher[18] * fisher[9])) - fisher[6] * (fisher[2] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[19] - fisher[18] * fisher[4])) + fisher[16] * (fisher[2] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) - fisher[7] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[9] - fisher[8] * fisher[4])) - fisher[21] * (fisher[2] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) - fisher[7] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[17] * (fisher[3] * fisher[9] - fisher[8] * fisher[4]))) - fisher[15] * (fisher[1] * (fisher[7] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[14] - fisher[13] * fisher[9])) - fisher[6] * (fisher[2] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[14] - fisher[13] * fisher[4])) + fisher[11] * (fisher[2] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) - fisher[7] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[9] - fisher[8] * fisher[4])) - fisher[21] * (fisher[2] * (fisher[8] * fisher[14] - fisher[13] * fisher[9]) - fisher[7] * (fisher[3] * fisher[14] - fisher[13] * fisher[4]) + fisher[12] * (fisher[3] * fisher[9] - fisher[8] * fisher[4]))) + fisher[20] * (fisher[1] * (fisher[7] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) + fisher[17] * (fisher[8] * fisher[14] - fisher[13] * fisher[9])) - fisher[6] * (fisher[2] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[17] * (fisher[3] * fisher[14] - fisher[13] * fisher[4])) + fisher[11] * (fisher[2] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) - fisher[7] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[17] * (fisher[3] * fisher[9] - fisher[8] * fisher[4])) - fisher[16] * (fisher[2] * (fisher[8] * fisher[14] - fisher[13] * fisher[9]) - fisher[7] * (fisher[3] * fisher[14] - fisher[13] * fisher[4]) + fisher[12] * (fisher[3] * fisher[9] - fisher[8] * fisher[4])));
	return det;
}


/*

Global Functions

*/


// localize 9 *updating for 3D 2/7/18
__global__ void  localize9(float *d_iall,  // pointer to our image variable
	float *d_xf_all,	// pointer for final x-coordinate measurement
	float *d_yf_all,	// pointer for final y-coordinate measurement
	float *d_zf_all,	// pointer for final z-coordinate measurement
	float *d_N,	// pointer for final N measurement
	float *d_off,// pointer for final offset measurement
	float *d_xf_crlb,// pointer for final x-coordinate uncertainty
	float *d_yf_crlb,// pointer for final y-coordinate uncertainty
	float *d_zf_crlb,// pointer for final z-coordinate uncertainty
	float *d_N_crlb,// pointer for final N uncertainty
	float *d_off_crlb,// pointer for final offset uncertainty
	float *d_llv,// pointer for final log likelihood value calculation
	float ang,// rotation of the fitting grid in radians
	float *d_zcurve,// pointer for defocusing constants in [sxo, syo, ax, ay, bx, by, gx, gy]
	int numi)// pointer for final x-coordinate measurement
{

	// Declare variables
	int pix = 9;   // number of pixels in the localization image
	__shared__ float xgrid[81];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ float ygrid[81];			// this will reduce calls to global device memory

											// Assign defocusing constants from zcurve
	__shared__ float sxo;
	sxo = d_zcurve[0];
	__shared__ float syo;
	syo = d_zcurve[1];
	__shared__ float axo;
	axo = d_zcurve[2];
	__shared__ float ayo;
	ayo = d_zcurve[3];
	__shared__ float bxo;
	bxo = d_zcurve[4];
	__shared__ float byo;
	byo = d_zcurve[5];
	__shared__ float dxo;
	dxo = d_zcurve[6];
	__shared__ float dyo;
	dyo = d_zcurve[7];
	__shared__ float gxo;
	gxo = d_zcurve[8];
	__shared__ float gyo;
	gyo = d_zcurve[9];

	// local register variables
	float dudx, dudy, dudz, dudsx, dudsy, d2udx2, d2udy2, d2udz2, d2udsx2, d2udsy2, dudn, dudo, Ex, Ey, u;
	float dsxdz, dsydz, d2sxdz2, d2sydz2;
	float d_x, d_y, d_z, d_n, d_o, dd_x, dd_y, dd_z, dd_n, dd_o, x, y, sx, sy;
	// fitting parameters
	float xf, yf, zf, N, b;
	int tx = threadIdx.x;
	int index = blockIdx.x*blockDim.x + tx;		// calculate thread index
	float d_i2[81];				// initialize data for image
	float llv;
	float fisher[25] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	float det_fish = 0.0;
	// create xgrid and ygrid we want to create the grid regardless of whether the index is crunching on an image
	if (tx == 0) {
		for (int i = 0; i <pix; i++) {  // here the grid is constructed by assigning x and y the pixel index, then rotating with a rotation transform by known angle ang
			for (int j = 0; j <pix; j++) {
				x = (float)j - ((float)pix - 1) / 2;
				y = (float)i - ((float)pix - 1) / 2;
				xgrid[j*pix + i] = x*cos(ang) - y*sin(ang);
				ygrid[j*pix + i] = x*sin(ang) + y*cos(ang);
			}
		}
	}

	if (index < numi) {			// check to see that threads only work if an image exists
								// buffer all the variables into shared memory and registers and build guesses

		xf = 0.0; // xf
		yf = 0.0; // yf
		N = 0.0; // N
		zf = 0;  // set z to be 0 which is close to the disk of least confusion
		b = 100000; // offset
		for (int i = 0; i <pix*pix; i++) {
			d_i2[i] = d_iall[i + index*pix*pix];	         // this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
			xf += xgrid[i] * d_i2[i];						 // sum of x and image weight
			yf += ygrid[i] * d_i2[i];						 // sum of y and image weight
			N += d_i2[i];									 // image sum
			if (b > d_i2[i]) { b = d_i2[i]; }				 // find minimum of image
		}
		xf = xf / N;
		yf = yf / N;
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < 20; counttry++) {
			d_x = 0.0;
			d_y = 0.0;
			d_z = 0.0;
			d_n = 0.0;
			d_o = 0.0;
			dd_x = 0.0; 	//wipe incremental variables each loop to give correct correction factor
			dd_y = 0.0;
			dd_z = 0.0;
			dd_n = 0.0;
			dd_o = 0.0;
			u = 0;
			Ey = 0;
			Ex = 0;
			llv = 0.0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < pix; rowcount++) {	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < pix; colcount++) {	// FOR 3 loops over all columns

					sx = sxo*powf(1 + powf((zf - gxo) / dxo, 2.0) + axo*powf((zf - gxo) / dxo, 3.0) + bxo*powf((zf - gxo) / dxo, 4.0), 0.5);
					sy = syo*powf(1 + powf((zf - gyo) / dyo, 2.0) + ayo*powf((zf - gyo) / dyo, 3.0) + byo*powf((zf - gyo) / dyo, 4.0), 0.5);
					// x/ygrid is col major(come from matlab) and i3 is col major 
					// these three lines help define the fitting gaussian as deined by the current iteration of parameters
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*pix] - xf + 0.5) / (2.0 * sx * sx)) - erf((xgrid[rowcount + colcount*pix] - xf - 0.5) / (2.0 * sx * sx)));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*pix] - yf + 0.5) / (2.0 * sy * sy)) - erf((ygrid[rowcount + colcount*pix] - yf - 0.5) / (2.0 * sy * sy)));
					u = N * Ex*Ey + b;
					// first derivatives calculations
					// these are done pixel by pixel with the sum added up in the d_x and dd_x areas
					dudx = (N / sqrt(2.0 * PI*sx * sx))*(exp(-powf(xgrid[rowcount + colcount*pix] - xf - 0.5, 2.0) / (2.0 * sx * sx))
						- exp(-powf(xgrid[rowcount + colcount*pix] - xf + 0.5, 2.0) / (2.0 * sx * sx)))*Ey;
					dudy = (N / sqrt(2.0 * PI*sy * sy))*(exp(-powf(ygrid[rowcount + colcount*pix] - yf - 0.5, 2.0) / (2.0 * sy * sy))
						- exp(-powf(ygrid[rowcount + colcount*pix] - yf + 0.5, 2.0) / (2.0 * sy * sy)))*Ex;
					dudsx = (N *Ey / (sqrt(2.0*PI) * powf(sx, 2.0)))*((xgrid[rowcount + colcount*pix] - xf - 0.5) * exp(-powf(xgrid[rowcount + colcount*pix] - xf - 0.5, 2.0) / (2.0 * sx * sx))
						- (xgrid[rowcount + colcount*pix] - xf + 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - xf + 0.5, 2.0) / (2.0 * sx * sx)));
					dudsy = (N *Ex / (sqrt(2.0*PI) * powf(sy, 2.0)))*((ygrid[rowcount + colcount*pix] - yf - 0.5) * exp(-powf(ygrid[rowcount + colcount*pix] - yf - 0.5, 2.0) / (2.0 * sy * sy))
						- (ygrid[rowcount + colcount*pix] - yf + 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - yf + 0.5, 2.0) / (2.0 * sy * sy)));
					dudn = Ex*Ey;
					dsxdz = sxo*(2 * (zf - gxo) / (dxo*dxo) + axo * 3 * powf((zf - gxo), 2) / powf(dxo, 3) + bxo * 4 * powf((zf - gxo), 3) / powf(dxo, 4)) /
						(2 * powf(1 + powf((zf - gxo) / dxo, 2.0) + axo*powf((zf - gxo) / dxo, 3.0) + bxo*powf((zf - gxo) / dxo, 4.0), 0.5));
					dsydz = syo*(2 * (zf - gyo) / (dyo*dyo) + ayo * 3 * powf((zf - gyo), 2) / powf(dyo, 3) + byo * 4 * powf((zf - gyo), 3) / powf(dyo, 4)) /
						(2 * powf(1 + powf((zf - gyo) / dyo, 2.0) + ayo*powf((zf - gyo) / dyo, 3.0) + byo*powf((zf - gyo) / dyo, 4.0), 0.5));
					dudz = dudsx*dsxdz + dudsy*dsydz;
					dudo = 1.0;

					// second derivatives
					// these are calcualted in a  similar manner to the first derivatives
					d2udx2 = (N / (sqrt(2.0 * PI)*powf(sx, 3.0))*((xgrid[rowcount + colcount*pix] - xf - 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - xf - 0.5, 2.0) / (2.0 * sx * sx))
						- (xgrid[rowcount + colcount*pix] - xf + 0.5)*exp(-powf(xgrid[rowcount + colcount*pix] - xf + 0.5, 2.0) / (2.0 * sx * sx))))*Ey;
					d2udy2 = (N / (sqrt(2.0 * PI)*powf(sy, 3.0))*((ygrid[rowcount + colcount*pix] - yf - 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - yf - 0.5, 2.0) / (2.0 * sy * sy))
						- (ygrid[rowcount + colcount*pix] - yf + 0.5)*exp(-powf(ygrid[rowcount + colcount*pix] - yf + 0.5, 2.0) / (2.0 * sy * sy))))*Ex;
					d2udsx2 = (Ey*N / (sqrt(2.0 * PI)))
						*(powf(sx, -5.0)*(powf((xgrid[rowcount + colcount*pix] - xf - 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*pix] - xf - 0.5, 2.0) / (2.0 * sx * sx))
							- powf((xgrid[rowcount + colcount*pix] - xf + 0.5), 3)*exp(-powf(xgrid[rowcount + colcount*pix] - xf + 0.5, 2.0) / (2.0 * sx * sx)))
							- 2 * powf(sx, -3.0)*((xgrid[rowcount + colcount*pix] - xf - 0.5)     *exp(-powf(xgrid[rowcount + colcount*pix] - xf - 0.5, 2.0) / (2.0 * sx * sx))
								- (xgrid[rowcount + colcount*pix] - xf + 0.5) *exp(-powf(xgrid[rowcount + colcount*pix] - xf + 0.5, 2.0) / (2.0 * sx * sx))));
					d2udsy2 = (Ex*N / (sqrt(2.0 * PI)))
						*(powf(sy, -5.0)*(powf((ygrid[rowcount + colcount*pix] - yf - 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*pix] - yf - 0.5, 2.0) / (2.0 * sy * sy))
							- powf((ygrid[rowcount + colcount*pix] - yf + 0.5), 3)*exp(-powf(ygrid[rowcount + colcount*pix] - yf + 0.5, 2.0) / (2.0 * sy * sy)))
							- 2 * powf(sy, -3.0)*((ygrid[rowcount + colcount*pix] - yf - 0.5)     *exp(-powf(ygrid[rowcount + colcount*pix] - yf - 0.5, 2.0) / (2.0 * sy * sy))
								- (ygrid[rowcount + colcount*pix] - yf + 0.5) *exp(-powf(ygrid[rowcount + colcount*pix] - yf + 0.5, 2.0) / (2.0 * sy * sy))));
					d2sxdz2 = sxo*(2 / powf(dxo, 2.0) + axo * 6 * (zf - gxo) / powf(dxo, 3) + bxo * 12 * powf(zf - gxo, 2.0) / powf(dxo, 4)) /
						(2 * powf(1 + powf((zf - gxo) / dxo, 2.0) + axo*powf((zf - gxo) / dxo, 3.0) + bxo*powf((zf - gxo) / dxo, 4.0), 0.5)) -
						sxo*powf(2 * (zf - gxo) / powf(dxo, 2) + axo * 3 * powf(zf - gxo, 2) / powf(dxo, 3) + bxo * 4 * powf(zf - gxo, 3) / powf(dxo, 4), 2) /
						(4 * powf(1 + powf((zf - gxo) / dxo, 2.0) + axo*powf((zf - gxo) / dxo, 3.0) + bxo*powf((zf - gxo) / dxo, 4.0), 1.5));
					d2sydz2 = syo*(2 / powf(dyo, 2.0) + ayo * 6 * (zf - gyo) / powf(dyo, 3) + byo * 12 * powf(zf - gyo, 2.0) / powf(dyo, 4)) /
						(2 * powf(1 + powf((zf - gyo) / dyo, 2.0) + ayo*powf((zf - gyo) / dyo, 3.0) + byo*powf((zf - gyo) / dyo, 4.0), 0.5)) -
						syo*powf(2 * (zf - gyo) / powf(dyo, 2) + ayo * 3 * powf(zf - gyo, 2) / powf(dyo, 3) + byo * 4 * powf(zf - gyo, 3) / powf(dyo, 4), 2) /
						(4 * powf(1 + powf((zf - gyo) / dyo, 2.0) + ayo*powf((zf - gyo) / dyo, 3.0) + byo*powf((zf - gyo) / dyo, 4.0), 1.5));
					d2udz2 = d2udsx2*powf(dsxdz, 2) + dudsx*d2sxdz2 + d2udsy2*powf(dsydz, 2) + dudsy*d2sydz2;

					// summing variable to lead to correction factors
					// these variables keep track of the correction which is given by summing over the entire pixel
					d_x = d_x + dudx*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_z = d_z + dudz*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_z = dd_z + d2udz2*((d_i2[rowcount + colcount*pix] / u) - 1.0) - powf(dudz, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*pix] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*pix] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*pix] / powf(u, 2.0);


					if (counttry == 19) {  // on the last count, construct fisher information matrix elements
						fisher[0] += dudx*dudx / u;
						fisher[1] += dudx*dudy / u;
						fisher[2] += dudx*dudn / u;
						fisher[3] += dudx*dudo / u;
						fisher[4] += dudx*dudz / u;

						fisher[5] += dudy*dudx / u;
						fisher[6] += dudy*dudy / u;
						fisher[7] += dudy*dudn / u;
						fisher[8] += dudy*dudo / u;
						fisher[9] += dudy*dudz / u;

						fisher[10] += dudn*dudx / u;  // the format has been updated but not the mathematics 2/7/18
						fisher[11] += dudn*dudy / u;
						fisher[12] += dudn*dudn / u;
						fisher[13] += dudn*dudo / u;
						fisher[14] += dudn*dudz / u;

						fisher[15] += dudo*dudx / u;
						fisher[16] += dudo*dudy / u;
						fisher[17] += dudo*dudn / u;
						fisher[18] += dudo*dudo / u;
						fisher[19] += dudo*dudz / u;

						fisher[20] += dudz*dudx / u;
						fisher[21] += dudz*dudy / u;
						fisher[22] += dudz*dudn / u;
						fisher[23] += dudz*dudo / u;
						fisher[24] += dudz*dudz / u;

						llv += d_i2[rowcount + colcount*pix] * log(u + 0.0000000000000001) - u - d_i2[rowcount + colcount*pix] * log(d_i2[rowcount + colcount*pix] + 0.0000000000000001) + d_i2[rowcount + colcount*pix];
					}
				} // END FOR 3
			} // END FOR2
			  // correct beta1 values with tolerances
			xf = xf - d_x / dd_x;
			yf = yf - d_y / dd_y;
			zf = zf - d_z / dd_z;
			N = N - d_n / dd_n;
			b = b - d_o / dd_o;

		}   // end FOR 1

		if (xf == xf && yf == yf && zf == zf && N == N && b == b && sx == sx && sy == sy && b == b) { // begin is numeric if statement
			if (N > 0 && xf >= -3 && xf <= 3 && yf <= 3 && yf >= -3) { // was the molecule inside the image? Was N positive? if yes then record the point
				d_xf_all[index] = xf; // correct position for x
				d_yf_all[index] = yf; // correct position for y
				d_zf_all[index] = zf;
				d_N[index] = N;
				d_off[index] = b;
				d_llv[index] = llv;

				// calculate crlb's for estimators	

				// updated for zf
				det_fish = device_det(fisher);  // these values were determined using a homemade Python code called cofacs.py and text_det.py and checking against lower rank matricies
				d_xf_crlb[index] = (fisher[6] * (fisher[12] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) + fisher[22] * (fisher[13] * fisher[19] - fisher[18] * fisher[14])) - fisher[11] * (fisher[7] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[19] - fisher[18] * fisher[9])) + fisher[16] * (fisher[7] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[22] * (fisher[8] * fisher[14] - fisher[13] * fisher[9])) - fisher[21] * (fisher[7] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) + fisher[17] * (fisher[8] * fisher[14] - fisher[13] * fisher[9]))) / det_fish;
				d_yf_crlb[index] = -(fisher[0] * (fisher[12] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) + fisher[22] * (fisher[13] * fisher[19] - fisher[18] * fisher[14])) - fisher[10] * (fisher[2] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[17] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[19] - fisher[18] * fisher[4])) + fisher[15] * (fisher[2] * (fisher[13] * fisher[24] - fisher[23] * fisher[14]) - fisher[12] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[22] * (fisher[3] * fisher[14] - fisher[13] * fisher[4])) - fisher[20] * (fisher[2] * (fisher[13] * fisher[19] - fisher[18] * fisher[14]) - fisher[12] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[17] * (fisher[3] * fisher[14] - fisher[13] * fisher[4]))) / det_fish;
				d_N_crlb[index] = +(fisher[0] * (fisher[6] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[16] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) + fisher[21] * (fisher[8] * fisher[19] - fisher[18] * fisher[9])) - fisher[5] * (fisher[1] * (fisher[18] * fisher[24] - fisher[23] * fisher[19]) - fisher[16] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[21] * (fisher[3] * fisher[19] - fisher[18] * fisher[4])) + fisher[15] * (fisher[1] * (fisher[8] * fisher[24] - fisher[23] * fisher[9]) - fisher[6] * (fisher[3] * fisher[24] - fisher[23] * fisher[4]) + fisher[21] * (fisher[3] * fisher[9] - fisher[8] * fisher[4])) - fisher[20] * (fisher[1] * (fisher[8] * fisher[19] - fisher[18] * fisher[9]) - fisher[6] * (fisher[3] * fisher[19] - fisher[18] * fisher[4]) + fisher[16] * (fisher[3] * fisher[9] - fisher[8] * fisher[4]))) / det_fish;
				d_off_crlb[index] = -(fisher[0] * (fisher[6] * (fisher[12] * fisher[24] - fisher[22] * fisher[14]) - fisher[11] * (fisher[7] * fisher[24] - fisher[22] * fisher[9]) + fisher[21] * (fisher[7] * fisher[14] - fisher[12] * fisher[9])) - fisher[5] * (fisher[1] * (fisher[12] * fisher[24] - fisher[22] * fisher[14]) - fisher[11] * (fisher[2] * fisher[24] - fisher[22] * fisher[4]) + fisher[21] * (fisher[2] * fisher[14] - fisher[12] * fisher[4])) + fisher[10] * (fisher[1] * (fisher[7] * fisher[24] - fisher[22] * fisher[9]) - fisher[6] * (fisher[2] * fisher[24] - fisher[22] * fisher[4]) + fisher[21] * (fisher[2] * fisher[9] - fisher[7] * fisher[4])) - fisher[20] * (fisher[1] * (fisher[7] * fisher[14] - fisher[12] * fisher[9]) - fisher[6] * (fisher[2] * fisher[14] - fisher[12] * fisher[4]) + fisher[11] * (fisher[2] * fisher[9] - fisher[7] * fisher[4]))) / det_fish;
				d_zf_crlb[index] = +(fisher[0] * (fisher[6] * (fisher[12] * fisher[18] - fisher[17] * fisher[13]) - fisher[11] * (fisher[7] * fisher[18] - fisher[17] * fisher[8]) + fisher[16] * (fisher[7] * fisher[13] - fisher[12] * fisher[8])) - fisher[5] * (fisher[1] * (fisher[12] * fisher[18] - fisher[17] * fisher[13]) - fisher[11] * (fisher[2] * fisher[18] - fisher[17] * fisher[3]) + fisher[16] * (fisher[2] * fisher[13] - fisher[12] * fisher[3])) + fisher[10] * (fisher[1] * (fisher[7] * fisher[18] - fisher[17] * fisher[8]) - fisher[6] * (fisher[2] * fisher[18] - fisher[17] * fisher[3]) + fisher[16] * (fisher[2] * fisher[8] - fisher[7] * fisher[3])) - fisher[15] * (fisher[1] * (fisher[7] * fisher[13] - fisher[12] * fisher[8]) - fisher[6] * (fisher[2] * fisher[13] - fisher[12] * fisher[3]) + fisher[11] * (fisher[2] * fisher[8] - fisher[7] * fisher[3]))) / det_fish;
			}
			else {						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_zf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_zf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_llv[index] = llv;
			}
		} //end is numeric if statement
		else {
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_zf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_zf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_llv[index] = llv;
		} // end else fail statement
	}
}   // end localize 9

	/*
	* Host code
	*
	*
	*/
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	float *iall;				// the pointer to the array of all images to be analyzed
	float *d_iall;			// Pointer to image array on gpu
	float *zcurve;
	float *d_zcurve;
	float angle;
	float *d_xf_all;
	float *d_yf_all;
	float *d_zf_all;
	float *d_N;
	float *d_off;
	float *d_llv;
	float *d_xf_crlb;
	float *d_yf_crlb;
	float *d_zf_crlb;
	float *d_N_crlb;
	float *d_off_crlb;
	float *xf, *xfc, *yf, *yfc, *n, *nc, *zf, *zfc, *off, *offc, *llv;
	size_t threadsperblock;


	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int numi;				// number of images imported
	const size_t *idims;


	/* Throw an error if the input does not match expectations. */
	if (nrhs != 4) {
		printf("Must have 4 inputs ( i1, numthreads, angle(in rads), defocusing constants)\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// Error statement to handle ilocs being type single
	if (!mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
		printf("i1 must be a nxm float array\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// Error statement if angle is not type single
	if (!mxIsSingle(prhs[3]) || mxIsComplex(prhs[3])) {
		printf("angle must be a single\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsSingle(prhs[1]) || mxIsComplex(prhs[1])) {
		printf("defocus constants must be of type single\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// get pointer to input arguments
	iall = (float *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	irow = (int)idims[0];
	numi = (int)idims[1];
	angle = (float)mxGetScalar(prhs[3]);
	zcurve = (float *)mxGetPr(prhs[1]);
	if (numi > 1000000 || numi < 1) {
		numi = 1;
	}

	int imem = irow*numi * sizeof(float);
	int vmem = numi * sizeof(float);

	// verify that the input variables are what was expected

	// check that iloc is a *perfect square* by num_mol array
	if (irow != 49 && irow != 81 && irow != 121) {
		printf("Images are of incorrect size. There must be a perfect square number of rows in the entry.\n");
		mexErrMsgTxt("See Error above!\n");
	}


	if (nlhs != 11) {
		printf("You must have 11 output variables [xf_all, yf_all, zf_all, N, off_all, xf_crlb, yf_crlb, zf_crlb, N_crlb, off_crlb, llv_all]\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// allocate memory and copy it onto the gpu device
	// iall
	checkCudaErrors(cudaMalloc((void**)&d_iall, imem)); // allocate image memory
	checkCudaErrors(cudaMalloc((void**)&d_zcurve, sizeof(float) * 9)); // allocate image memory
	checkCudaErrors(cudaMemcpy(d_iall, iall, imem, cudaMemcpyHostToDevice)); // copy images from device to host
	checkCudaErrors(cudaMemcpy(d_zcurve, zcurve, sizeof(float) * 9, cudaMemcpyHostToDevice));


	// allocate memory for fitted variables that will be returned from device
	checkCudaErrors(cudaMalloc((void**)&d_xf_all, vmem)); // allocate xf_all memory
	checkCudaErrors(cudaMalloc((void**)&d_xf_crlb, vmem)); // allocate xf_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_yf_all, vmem)); // allocate yf_all memory
	checkCudaErrors(cudaMalloc((void**)&d_yf_crlb, vmem)); // allocate yf_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_zf_all, vmem));   // allocate zf memory
	checkCudaErrors(cudaMalloc((void**)&d_zf_crlb, vmem));   // allocate zf_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_N, vmem)); // allocate N memory
	checkCudaErrors(cudaMalloc((void**)&d_N_crlb, vmem)); // allocate N_crlb memory
	checkCudaErrors(cudaMalloc((void**)&d_off, vmem)); // allocate off memory
	checkCudaErrors(cudaMalloc((void**)&d_off_crlb, vmem)); // allocate N memory
	checkCudaErrors(cudaMalloc((void**)&d_llv, vmem)); // allocate llv memory

													   /* Run GPU kernel*/
	threadsperblock = mxGetScalar(prhs[2]);  // get number of threads perblock from matlab

	localize9 << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_xf_all, d_yf_all, d_zf_all, d_N, d_off, d_xf_crlb, d_yf_crlb, d_zf_crlb, d_N_crlb, d_off_crlb, d_llv, angle, d_zcurve, numi);



	// Allocate host side memory for output arrays at the output pointer positions
	plhs[0] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[3] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[4] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[5] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[6] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);		 // checked so far as of 2/7/18 for updates to fit zf from full_chain_loc
	plhs[7] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[8] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[9] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);
	plhs[10] = mxCreateNumericMatrix(numi, 1, mxSINGLE_CLASS, mxREAL);

	// Copy pointers from mex array
	xf = (float *)mxGetPr(plhs[0]);
	xfc = (float *)mxGetPr(plhs[1]);
	yf = (float *)mxGetPr(plhs[2]);
	yfc = (float *)mxGetPr(plhs[3]);
	zf = (float *)mxGetPr(plhs[4]);
	zfc = (float *)mxGetPr(plhs[5]);
	n = (float *)mxGetPr(plhs[6]);
	nc = (float *)mxGetPr(plhs[7]);
	off = (float *)mxGetPr(plhs[8]);
	offc = (float *)mxGetPr(plhs[9]);
	llv = (float *)mxGetPr(plhs[10]);   // checked so far as of 2/7/18 for updates to fit zf from full_chain_loc

										// copy memory from device to host memory at mex array pointers
	checkCudaErrors(cudaMemcpy(xf, d_xf_all, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(xfc, d_xf_crlb, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(yf, d_yf_all, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(yfc, d_yf_crlb, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(zf, d_zf_all, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(zfc, d_zf_crlb, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(n, d_N, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(nc, d_N_crlb, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(off, d_off, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(offc, d_off_crlb, vmem, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(llv, d_llv, vmem, cudaMemcpyDeviceToHost)); // checked so far as of 2/7/18 for updates to fit zf from full_chain_loc

																		   // clean up
	cudaFree(d_iall);
	cudaFree(d_N);
	cudaFree(d_xf_all);
	cudaFree(d_yf_all);
	cudaFree(d_zf_all);
	cudaFree(d_off);
	cudaFree(d_xf_crlb);
	cudaFree(d_yf_crlb);
	cudaFree(d_zf_crlb);
	cudaFree(d_N_crlb);
	cudaFree(d_off_crlb);
	cudaFree(d_llv);
	cudaFree(d_zcurve);

}  // DONE
