/*
 * Chain_loc v 1.0 is the source code for a mex file that will input image data and parameter estimators and output localized data
 * Calling this function in matlab will look like
 * [xf_all, yf_all, N, off_all, xf_crlb, yf_crlb, N_crlb, off_crlb, llv] = chain_loc(i1, beta0, framenum_temp, xpeaks, ypeaks,xpix, ypix, numthreads)
 * written by Andrew Nelson Version v 1.1 on 5/2/15
 */

/*
Version 1.2 has had substantial debugging and careful comparison between cpu version. Removed extra diagnostic code that has been commented out and has no use in regular code
	At this point we can be sured that the algorithm for calculating position, total area, and offset are working properly on the GPU
	This version is usable with the localiztion code Quhzx_01_3.m

Fixed
	Error codes should not be given off when inputs and outputs don't match expected arguments
	Fixed a problem considerably prolonging computation time by removing gpu device invokation early in main loop ( i.e. reset the gpu, synchronize threads, initialize gpu ....)
	Added redundant void __global__ loops to handle multiple sizes of input images

Still possibly broken
	CRLB generating equations have not been double / triple checked yet these values are not yet reliable
*/

#include "mex.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846



/*
 * Device code
 * 
 * To facilitate coding (for me) I have copied the localization algorithm to be used with multiple sized areas
 */

void __global__ localize7(double *d_iall,
	double *d_beta0,
	double *d_x_cm,
	double *d_y_cm,
	double *d_framenum_temp,
	double *d_xf_all,
	double *d_yf_all,
	double *d_N,
	double *d_off,
	double *d_framenum_all,
	double *d_xf_crlb,
	double *d_yf_crlb,
	double *d_N_crlb,
	double *d_off_crlb,
	int irow,
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi)
{
	// Declare variables
	__shared__ double xgrid[49];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[49];			// this will reduce calls to global device memory

	// these variables will exist on the register of each thread
	double dudx, dudy, d2udx2, d2udy2, dudn, dudo, Ex, Ey, u, llv;
	int index = blockIdx.x*blockDim.x + threadIdx.x;		// calculate thread index
	double d_i2[49];				// initialize data for image
	double d_beta1[5];				// initialize data for beta1
	double d_x, d_y, d_n, d_o, dd_x, dd_y, dd_n, dd_o;
	double fisher[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double det_fish;
	if (index < numi){			// check to see that threads only work if an image exists
		// buffer all the variables into shared memory and registers
		for (int i = 0; i <49; i++) {
			xgrid[i] = d_xpix[i];
			ygrid[i] = d_ypix[i];
			d_i2[i] = d_iall[i + index * 49];		// this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
		}
		for (int j = 0; j<5; j++){
			d_beta1[j] = d_beta0[j + index * 5];	// similar buffering strategy employed for d_beta1
		}

		
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < 16; counttry++){
			d_x = 0.0;
			d_y = 0.0;
			d_n = 0.0;
			d_o = 0.0;
			dd_x = 0.0;	//wipe incremental variables
			dd_y = 0.0;
			dd_n = 0.0;
			dd_o = 0.0;
			u = 0;
			Ey = 0;
			Ex = 0;
			llv = 0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < irow; rowcount++){	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < irow; colcount++){	// FOR 3 loops over all columns
					if (d_beta1[4] < 0){
						d_beta1[4] = 0;
					}
					// x/ygrid is col major(come from matlab) and i3 is col major 
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					u = d_beta1[2] * Ex*Ey + d_beta1[4];

					// first derivatives
					dudx = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudy = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ex;
					dudn = Ex*Ey;
					dudo = 1.0;

					// second derivatives
					d2udx2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ey;
					d2udy2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ex;

					// summing variable to lead to correction factors
					d_x = d_x + dudx*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*irow] / powf(u, 2.0);

					
					if (counttry == 15){  // on the last count, construct fisher information matrix elements
						fisher[0] += dudx*dudx / u;
						fisher[4] += dudx*dudy / u;
						fisher[8] += dudx*dudn / u;
						fisher[12] += dudx*dudo / u;

						fisher[1] += dudy*dudx / u;
						fisher[5] += dudy*dudy / u;
						fisher[9] += dudy*dudn / u;
						fisher[13] += dudy*dudo / u;

						fisher[2] += dudn*dudx / u;
						fisher[6] += dudn*dudy / u;
						fisher[10] += dudn*dudn / u;
						fisher[14] += dudn*dudo / u;

						fisher[3] += dudo*dudx / u;
						fisher[7] += dudo*dudy / u;
						fisher[11] += dudo*dudn / u;
						fisher[15] += dudo*dudo / u;
						llv += llv  + d_i2[rowcount+colcount*irow] *  log(u) - u - d_i2[rowcount + colcount*irow]*log(d_i2[rowcount + colcount*irow]) + d_i2[rowcount + colcount*irow];
					}
				} // END FOR 3
			} // END FOR2
			// correct beta1 values with tolerances
			d_beta1[0] = d_beta1[0] - d_x / dd_x;
			d_beta1[1] = d_beta1[1] - d_y / dd_y;
			d_beta1[2] = d_beta1[2] - d_n / dd_n;
			d_beta1[4] = d_beta1[4] - d_o / dd_o;

		}   // end FOR 1


		if (d_beta1[0] == d_beta1[0] && d_beta1[1] == d_beta1[1] && d_beta1[2] == d_beta1[2] && d_beta1[4] == d_beta1[4]){ // begin is numeric if statement
			if (d_beta1[2] > 0 && d_beta1[0] > xgrid[0] && d_beta1[0] < xgrid[irow*irow - 1] && d_beta1[1] < ygrid[irow*irow - 1] && d_beta1[1] > ygrid[0]){ // was the molecule inside the grid? Was N positive? if yes then record the point
				//d_xf_all[index] = d_beta1[0] + d_x_cm[index];	// correct position for x
				d_xf_all[index] = d_beta1[0] + d_x_cm[index];
				d_yf_all[index] = d_beta1[1] + d_y_cm[index]; // correct position for y
				d_N[index] = d_beta1[2];
				d_off[index] = d_beta1[4];
				d_framenum_all[index] = d_framenum_temp[index];
				d_llv[index] = llv;
				
				// calculate crlb's for estimators		
								// Determinant verified
				det_fish = fisher[0] * fisher[5] * fisher[10] * fisher[15] + fisher[0] * fisher[6] * fisher[11] * fisher[13] + fisher[0] * fisher[7] * fisher[9] * fisher[14] + fisher[1] * fisher[4] * fisher[11] * fisher[14] + fisher[1] * fisher[6] * fisher[8] * fisher[15] + fisher[1] * fisher[7] * fisher[10] * fisher[12] + fisher[2] * fisher[4] * fisher[9] * fisher[15] + fisher[2] * fisher[5] * fisher[11] * fisher[12] + fisher[2] * fisher[7] * fisher[8] * fisher[13] + fisher[3] * fisher[4] * fisher[10] * fisher[13] + fisher[3] * fisher[5] * fisher[8] * fisher[14] + fisher[3] * fisher[6] * fisher[9] * fisher[12] - fisher[0] * fisher[5] * fisher[11] * fisher[14] - fisher[0] * fisher[6] * fisher[9] * fisher[15] - fisher[0] * fisher[7] * fisher[10] * fisher[13] - fisher[1] * fisher[4] * fisher[10] * fisher[15] - fisher[1] * fisher[6] * fisher[11] * fisher[12] - fisher[1] * fisher[7] * fisher[8] * fisher[14] - fisher[2] * fisher[4] * fisher[11] * fisher[13] - fisher[2] * fisher[5] * fisher[8] * fisher[15] - fisher[2] * fisher[7] * fisher[9] * fisher[12] - fisher[3] * fisher[4] * fisher[9] * fisher[14] - fisher[3] * fisher[5] * fisher[10] * fisher[12] - fisher[3] * fisher[6] * fisher[8] * fisher[13];
				// xf_crlb verified
				d_xf_crlb[index] = (fisher[5] * fisher[10] * fisher[15] + fisher[9] * fisher[14] * fisher[7] + fisher[13] * fisher[6] * fisher[11] - fisher[5] * fisher[14] * fisher[11] - fisher[9] * fisher[6] * fisher[15] - fisher[13] * fisher[10] * fisher[7]) / det_fish;
				// yf_crlb verified
				d_yf_crlb[index] = (fisher[0] * fisher[10] * fisher[15] + fisher[2] * fisher[14] * fisher[3] + fisher[12] * fisher[2] * fisher[11] - fisher[0] * fisher[14] * fisher[11] - fisher[8] * fisher[2] * fisher[15] - fisher[12] * fisher[10] * fisher[3]) / det_fish;
				// N_crlb verified
				d_N_crlb[index] = (fisher[0] * fisher[5] * fisher[15] + fisher[4] * fisher[13] * fisher[3] + fisher[12] * fisher[1] * fisher[7] - fisher[0] * fisher[13] * fisher[7] - fisher[4] * fisher[1] * fisher[15] - fisher[12] * fisher[5] * fisher[3]) / det_fish;
				// off_crlb verified
				d_off_crlb[index] = (fisher[0] * fisher[5] * fisher[10] + fisher[4] * fisher[9] * fisher[2] + fisher[8] * fisher[1] * fisher[6] - fisher[0] * fisher[9] * fisher[6] - fisher[4] * fisher[1] * fisher[10] - fisher[8] * fisher[5] * fisher[2]) / det_fish;

			}
			else{						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_framenum_all[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_llv[index] = 1000;			// separates failed llv from successfull llv values
			} 
		} //end is numeric if statement
		else{
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_framenum_all[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_llv[index] = 1000;
		}
	}
}   // end localize7


void __global__ localize9(double *d_iall,
	double *d_beta0,
	double *d_x_cm,
	double *d_y_cm,
	double *d_framenum_temp,
	double *d_xf_all,
	double *d_yf_all,
	double *d_N,
	double *d_off,
	double *d_framenum_all,
	double *d_xf_crlb,
	double *d_yf_crlb,
	double *d_N_crlb,
	double *d_off_crlb,
	int irow,
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi)
{
	// Declare variables
	__shared__ double xgrid[81];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[81];			// this will reduce calls to global device memory

	// these variables will exist on the register of each thread
	double dudx, dudy, d2udx2, d2udy2, dudn, dudo, Ex, Ey, u, llv;
	int index = blockIdx.x*blockDim.x + threadIdx.x;		// calculate thread index
	double d_i2[81];				// initialize data for image
	double d_beta1[5];				// initialize data for beta1
	double d_x, d_y, d_n, d_o, dd_x, dd_y, dd_n, dd_o;
	double fisher[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double det_fish;
	if (index < numi){			// check to see that threads only work if an image exists
		// buffer all the variables into shared memory and registers
		for (int i = 0; i <49; i++) {
			xgrid[i] = d_xpix[i];
			ygrid[i] = d_ypix[i];
			d_i2[i] = d_iall[i + index * 49];		// this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
		}
		for (int j = 0; j<5; j++){
			d_beta1[j] = d_beta0[j + index * 5];	// similar buffering strategy employed for d_beta1
		}

		
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < 16; counttry++){
			d_x = 0.0;
			d_y = 0.0;
			d_n = 0.0;
			d_o = 0.0;
			dd_x = 0.0;	//wipe incremental variables
			dd_y = 0.0;
			dd_n = 0.0;
			dd_o = 0.0;
			u = 0;
			Ey = 0;
			Ex = 0;
			llv = 0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < irow; rowcount++){	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < irow; colcount++){	// FOR 3 loops over all columns
					if (d_beta1[4] < 0){
						d_beta1[4] = 0;
					}
					// x/ygrid is col major(come from matlab) and i3 is col major 
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					u = d_beta1[2] * Ex*Ey + d_beta1[4];

					// first derivatives
					dudx = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudy = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ex;
					dudn = Ex*Ey;
					dudo = 1.0;

					// second derivatives
					d2udx2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ey;
					d2udy2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ex;

					// summing variable to lead to correction factors
					d_x = d_x + dudx*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*irow] / powf(u, 2.0);

					
					if (counttry == 15){  // on the last count, construct fisher information matrix elements
						fisher[0] += dudx*dudx / u;
						fisher[4] += dudx*dudy / u;
						fisher[8] += dudx*dudn / u;
						fisher[12] += dudx*dudo / u;

						fisher[1] += dudy*dudx / u;
						fisher[5] += dudy*dudy / u;
						fisher[9] += dudy*dudn / u;
						fisher[13] += dudy*dudo / u;

						fisher[2] += dudn*dudx / u;
						fisher[6] += dudn*dudy / u;
						fisher[10] += dudn*dudn / u;
						fisher[14] += dudn*dudo / u;

						fisher[3] += dudo*dudx / u;
						fisher[7] += dudo*dudy / u;
						fisher[11] += dudo*dudn / u;
						fisher[15] += dudo*dudo / u;
						llv += llv  + d_i2[rowcount+colcount*irow] *  log(u) - u - d_i2[rowcount + colcount*irow]*log(d_i2[rowcount + colcount*irow]) + d_i2[rowcount + colcount*irow];
					}
				} // END FOR 3
			} // END FOR2
			// correct beta1 values with tolerances
			d_beta1[0] = d_beta1[0] - d_x / dd_x;
			d_beta1[1] = d_beta1[1] - d_y / dd_y;
			d_beta1[2] = d_beta1[2] - d_n / dd_n;
			d_beta1[4] = d_beta1[4] - d_o / dd_o;

		}   // end FOR 1


		if (d_beta1[0] == d_beta1[0] && d_beta1[1] == d_beta1[1] && d_beta1[2] == d_beta1[2] && d_beta1[4] == d_beta1[4]){ // begin is numeric if statement
			if (d_beta1[2] > 0 && d_beta1[0] > xgrid[0] && d_beta1[0] < xgrid[irow*irow - 1] && d_beta1[1] < ygrid[irow*irow - 1] && d_beta1[1] > ygrid[0]){ // was the molecule inside the grid? Was N positive? if yes then record the point
				//d_xf_all[index] = d_beta1[0] + d_x_cm[index];	// correct position for x
				d_xf_all[index] = d_beta1[0] + d_x_cm[index];
				d_yf_all[index] = d_beta1[1] + d_y_cm[index]; // correct position for y
				d_N[index] = d_beta1[2];
				d_off[index] = d_beta1[4];
				d_framenum_all[index] = d_framenum_temp[index];
				d_llv[index] = llv;
				
				// calculate crlb's for estimators		
								// Determinant verified
				det_fish = fisher[0] * fisher[5] * fisher[10] * fisher[15] + fisher[0] * fisher[6] * fisher[11] * fisher[13] + fisher[0] * fisher[7] * fisher[9] * fisher[14] + fisher[1] * fisher[4] * fisher[11] * fisher[14] + fisher[1] * fisher[6] * fisher[8] * fisher[15] + fisher[1] * fisher[7] * fisher[10] * fisher[12] + fisher[2] * fisher[4] * fisher[9] * fisher[15] + fisher[2] * fisher[5] * fisher[11] * fisher[12] + fisher[2] * fisher[7] * fisher[8] * fisher[13] + fisher[3] * fisher[4] * fisher[10] * fisher[13] + fisher[3] * fisher[5] * fisher[8] * fisher[14] + fisher[3] * fisher[6] * fisher[9] * fisher[12] - fisher[0] * fisher[5] * fisher[11] * fisher[14] - fisher[0] * fisher[6] * fisher[9] * fisher[15] - fisher[0] * fisher[7] * fisher[10] * fisher[13] - fisher[1] * fisher[4] * fisher[10] * fisher[15] - fisher[1] * fisher[6] * fisher[11] * fisher[12] - fisher[1] * fisher[7] * fisher[8] * fisher[14] - fisher[2] * fisher[4] * fisher[11] * fisher[13] - fisher[2] * fisher[5] * fisher[8] * fisher[15] - fisher[2] * fisher[7] * fisher[9] * fisher[12] - fisher[3] * fisher[4] * fisher[9] * fisher[14] - fisher[3] * fisher[5] * fisher[10] * fisher[12] - fisher[3] * fisher[6] * fisher[8] * fisher[13];
				// xf_crlb verified
				d_xf_crlb[index] = (fisher[5] * fisher[10] * fisher[15] + fisher[9] * fisher[14] * fisher[7] + fisher[13] * fisher[6] * fisher[11] - fisher[5] * fisher[14] * fisher[11] - fisher[9] * fisher[6] * fisher[15] - fisher[13] * fisher[10] * fisher[7]) / det_fish;
				// yf_crlb verified
				d_yf_crlb[index] = (fisher[0] * fisher[10] * fisher[15] + fisher[2] * fisher[14] * fisher[3] + fisher[12] * fisher[2] * fisher[11] - fisher[0] * fisher[14] * fisher[11] - fisher[8] * fisher[2] * fisher[15] - fisher[12] * fisher[10] * fisher[3]) / det_fish;
				// N_crlb verified
				d_N_crlb[index] = (fisher[0] * fisher[5] * fisher[15] + fisher[4] * fisher[13] * fisher[3] + fisher[12] * fisher[1] * fisher[7] - fisher[0] * fisher[13] * fisher[7] - fisher[4] * fisher[1] * fisher[15] - fisher[12] * fisher[5] * fisher[3]) / det_fish;
				// off_crlb verified
				d_off_crlb[index] = (fisher[0] * fisher[5] * fisher[10] + fisher[4] * fisher[9] * fisher[2] + fisher[8] * fisher[1] * fisher[6] - fisher[0] * fisher[9] * fisher[6] - fisher[4] * fisher[1] * fisher[10] - fisher[8] * fisher[5] * fisher[2]) / det_fish;

			}
			else{						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_framenum_all[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_llv[index] = 1000;			// separates failed llv from successfull llv values
			} 
		} //end is numeric if statement
		else{
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_framenum_all[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_llv[index] = 1000;
		}
	}
}   // end localize9



void __global__ localize11(double *d_iall,
	double *d_beta0,
	double *d_x_cm,
	double *d_y_cm,
	double *d_framenum_temp,
	double *d_xf_all,
	double *d_yf_all,
	double *d_N,
	double *d_off,
	double *d_framenum_all,
	double *d_xf_crlb,
	double *d_yf_crlb,
	double *d_N_crlb,
	double *d_off_crlb,
	int irow,
	double *d_xpix,
	double *d_ypix,
	double * d_llv,
	int numi)
{
	// Declare variables
	__shared__ double xgrid[121];			// allocate xpix and ypix variables to the shared memory of the blocks
	__shared__ double ygrid[121];			// this will reduce calls to global device memory

	// these variables will exist on the register of each thread
	double dudx, dudy, d2udx2, d2udy2, dudn, dudo, Ex, Ey, u, llv;
	int index = blockIdx.x*blockDim.x + threadIdx.x;		// calculate thread index
	double d_i2[121];				// initialize data for image
	double d_beta1[5];				// initialize data for beta1
	double d_x, d_y, d_n, d_o, dd_x, dd_y, dd_n, dd_o;
	double fisher[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double det_fish;
	if (index < numi){			// check to see that threads only work if an image exists
		// buffer all the variables into shared memory and registers
		for (int i = 0; i <49; i++) {
			xgrid[i] = d_xpix[i];
			ygrid[i] = d_ypix[i];
			d_i2[i] = d_iall[i + index * 49];		// this buffers [49] pixels into each d_i2, the thread index determines which image is analyzed
		}
		for (int j = 0; j<5; j++){
			d_beta1[j] = d_beta0[j + index * 5];	// similar buffering strategy employed for d_beta1
		}

		
		// start the for loop iterations				FOR 1
		for (int counttry = 0; counttry < 16; counttry++){
			d_x = 0.0;
			d_y = 0.0;
			d_n = 0.0;
			d_o = 0.0;
			dd_x = 0.0;	//wipe incremental variables
			dd_y = 0.0;
			dd_n = 0.0;
			dd_o = 0.0;
			u = 0;
			Ey = 0;
			Ex = 0;
			llv = 0;
			// Calculate pixel values for derivatives, 2nd derivatives, errorfunctions and u
			for (int rowcount = 0; rowcount < irow; rowcount++){	// FOR 2  loops over all rows
				for (int colcount = 0; colcount < irow; colcount++){	// FOR 3 loops over all columns
					if (d_beta1[4] < 0){
						d_beta1[4] = 0;
					}
					// x/ygrid is col major(come from matlab) and i3 is col major 
					Ex = 0.5 * (erf((xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					Ey = 0.5 * (erf((ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])) - erf((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5) / sqrt(2.0 * d_beta1[3] * d_beta1[3])));
					u = d_beta1[2] * Ex*Ey + d_beta1[4];

					// first derivatives
					dudx = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ey;
					dudy = (d_beta1[2] / sqrt(2.0 * PI*d_beta1[3] * d_beta1[3]))*(exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])))*Ex;
					dudn = Ex*Ey;
					dudo = 1.0;

					// second derivatives
					d2udx2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5)*exp(-powf(xgrid[rowcount + colcount*irow] - d_beta1[0] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ey;
					d2udy2 = (d_beta1[2] / (sqrt(2.0 * PI)*powf(d_beta1[3], 3.0))*((ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] - 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3])) - (ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5)*exp(-powf(ygrid[rowcount + colcount*irow] - d_beta1[1] + 0.5, 2.0) / (2.0 * d_beta1[3] * d_beta1[3]))))*Ex;

					// summing variable to lead to correction factors
					d_x = d_x + dudx*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_x = dd_x + d2udx2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudx, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_y = d_y + dudy*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_y = dd_y + d2udy2*((d_i2[rowcount + colcount*irow] / u) - 1.0) - powf(dudy, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2.0);
					d_n = d_n + dudn*((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_n = dd_n - powf(dudn, 2.0) * d_i2[rowcount + colcount*irow] / powf(u, 2);
					d_o = d_o + ((d_i2[rowcount + colcount*irow] / u) - 1.0);
					dd_o = dd_o - d_i2[rowcount + colcount*irow] / powf(u, 2.0);

					
					if (counttry == 15){  // on the last count, construct fisher information matrix elements
						fisher[0] += dudx*dudx / u;
						fisher[4] += dudx*dudy / u;
						fisher[8] += dudx*dudn / u;
						fisher[12] += dudx*dudo / u;

						fisher[1] += dudy*dudx / u;
						fisher[5] += dudy*dudy / u;
						fisher[9] += dudy*dudn / u;
						fisher[13] += dudy*dudo / u;

						fisher[2] += dudn*dudx / u;
						fisher[6] += dudn*dudy / u;
						fisher[10] += dudn*dudn / u;
						fisher[14] += dudn*dudo / u;

						fisher[3] += dudo*dudx / u;
						fisher[7] += dudo*dudy / u;
						fisher[11] += dudo*dudn / u;
						fisher[15] += dudo*dudo / u;
						llv += llv  + d_i2[rowcount+colcount*irow] *  log(u) - u - d_i2[rowcount + colcount*irow]*log(d_i2[rowcount + colcount*irow]) + d_i2[rowcount + colcount*irow];
					}
				} // END FOR 3
			} // END FOR2
			// correct beta1 values with tolerances
			d_beta1[0] = d_beta1[0] - d_x / dd_x;
			d_beta1[1] = d_beta1[1] - d_y / dd_y;
			d_beta1[2] = d_beta1[2] - d_n / dd_n;
			d_beta1[4] = d_beta1[4] - d_o / dd_o;

		}   // end FOR 1


		if (d_beta1[0] == d_beta1[0] && d_beta1[1] == d_beta1[1] && d_beta1[2] == d_beta1[2] && d_beta1[4] == d_beta1[4]){ // begin is numeric if statement
			if (d_beta1[2] > 0 && d_beta1[0] > xgrid[0] && d_beta1[0] < xgrid[irow*irow - 1] && d_beta1[1] < ygrid[irow*irow - 1] && d_beta1[1] > ygrid[0]){ // was the molecule inside the grid? Was N positive? if yes then record the point
				//d_xf_all[index] = d_beta1[0] + d_x_cm[index];	// correct position for x
				d_xf_all[index] = d_beta1[0] + d_x_cm[index];
				d_yf_all[index] = d_beta1[1] + d_y_cm[index]; // correct position for y
				d_N[index] = d_beta1[2];
				d_off[index] = d_beta1[4];
				d_framenum_all[index] = d_framenum_temp[index];
				d_llv[index] = llv;
				
				// calculate crlb's for estimators		
								// Determinant verified
				det_fish = fisher[0] * fisher[5] * fisher[10] * fisher[15] + fisher[0] * fisher[6] * fisher[11] * fisher[13] + fisher[0] * fisher[7] * fisher[9] * fisher[14] + fisher[1] * fisher[4] * fisher[11] * fisher[14] + fisher[1] * fisher[6] * fisher[8] * fisher[15] + fisher[1] * fisher[7] * fisher[10] * fisher[12] + fisher[2] * fisher[4] * fisher[9] * fisher[15] + fisher[2] * fisher[5] * fisher[11] * fisher[12] + fisher[2] * fisher[7] * fisher[8] * fisher[13] + fisher[3] * fisher[4] * fisher[10] * fisher[13] + fisher[3] * fisher[5] * fisher[8] * fisher[14] + fisher[3] * fisher[6] * fisher[9] * fisher[12] - fisher[0] * fisher[5] * fisher[11] * fisher[14] - fisher[0] * fisher[6] * fisher[9] * fisher[15] - fisher[0] * fisher[7] * fisher[10] * fisher[13] - fisher[1] * fisher[4] * fisher[10] * fisher[15] - fisher[1] * fisher[6] * fisher[11] * fisher[12] - fisher[1] * fisher[7] * fisher[8] * fisher[14] - fisher[2] * fisher[4] * fisher[11] * fisher[13] - fisher[2] * fisher[5] * fisher[8] * fisher[15] - fisher[2] * fisher[7] * fisher[9] * fisher[12] - fisher[3] * fisher[4] * fisher[9] * fisher[14] - fisher[3] * fisher[5] * fisher[10] * fisher[12] - fisher[3] * fisher[6] * fisher[8] * fisher[13];
				// xf_crlb verified
				d_xf_crlb[index] = (fisher[5] * fisher[10] * fisher[15] + fisher[9] * fisher[14] * fisher[7] + fisher[13] * fisher[6] * fisher[11] - fisher[5] * fisher[14] * fisher[11] - fisher[9] * fisher[6] * fisher[15] - fisher[13] * fisher[10] * fisher[7]) / det_fish;
				// yf_crlb verified
				d_yf_crlb[index] = (fisher[0] * fisher[10] * fisher[15] + fisher[2] * fisher[14] * fisher[3] + fisher[12] * fisher[2] * fisher[11] - fisher[0] * fisher[14] * fisher[11] - fisher[8] * fisher[2] * fisher[15] - fisher[12] * fisher[10] * fisher[3]) / det_fish;
				// N_crlb verified
				d_N_crlb[index] = (fisher[0] * fisher[5] * fisher[15] + fisher[4] * fisher[13] * fisher[3] + fisher[12] * fisher[1] * fisher[7] - fisher[0] * fisher[13] * fisher[7] - fisher[4] * fisher[1] * fisher[15] - fisher[12] * fisher[5] * fisher[3]) / det_fish;
				// off_crlb verified
				d_off_crlb[index] = (fisher[0] * fisher[5] * fisher[10] + fisher[4] * fisher[9] * fisher[2] + fisher[8] * fisher[1] * fisher[6] - fisher[0] * fisher[9] * fisher[6] - fisher[4] * fisher[1] * fisher[10] - fisher[8] * fisher[5] * fisher[2]) / det_fish;

			}
			else{						// if localization failed set all parameters to -1. These can easily be identified by molecules with framenum_all -1
				d_xf_all[index] = -1;
				d_yf_all[index] = -1;
				d_N[index] = -1;
				d_off[index] = -1;
				d_framenum_all[index] = -1;
				d_xf_crlb[index] = -1;
				d_yf_crlb[index] = -1;
				d_N_crlb[index] = -1;
				d_off_crlb[index] = -1;
				d_llv[index] = 1000;			// separates failed llv from successfull llv values
			} 
		} //end is numeric if statement
		else{
			d_xf_all[index] = -1;
			d_yf_all[index] = -1;
			d_N[index] = -1;
			d_off[index] = -1;
			d_framenum_all[index] = -1;
			d_xf_crlb[index] = -1;
			d_yf_crlb[index] = -1;
			d_N_crlb[index] = -1;
			d_off_crlb[index] = -1;
			d_llv[index] = 1000;
		}
	}
}   // end localize9





/*
 * Host code
 *
 * 
 */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, mxArray const *prhs[])
{
	/* Declare all variables.*/
	double *iall;				// the pointer to the array of all images to be analyzed
	double *beta0;				// the pointer	to all the initial estimates of parameters in order ( x coord, y coord, N, sigma, Background)
	double *framenum_temp;		// pointer to temporary array of frame numbers
	double *xpeaks;			// pointer to array of xpeaks
	double *ypeaks;			// pointer to array of y peaks
	double  *d_iall;			// Pointer to image array on gpu
	double *d_beta0;			// Pointer to initial estimates on gpu
	double *d_framenum_temp;
	double *d_framenum_all;
	double *d_x_cm;				// pointer to parameters on device
	double *d_y_cm;
	double *d_xf_all;
	double *d_yf_all;
	double *d_N;
	double *d_off;
	double *d_llv;
	double *d_xf_crlb;
	double *d_yf_crlb;
	double *d_N_crlb;
	double *d_off_crlb;
	double *xpix;
	double *ypix;
	double *d_xpix;
	double *d_ypix;
	int threadsperblock;


	int irow;				// number of pixels in a row which should also be the number in a coloumn
	int icol;
	int numi;				// number of images imported
	const int *idims;
	

	/* Throw an error if the input does not match expectations. */
	if (nrhs != 8) {
		printf("Must have 8 inputs ( i1, beta0, framenum_temp, xpeaks, ypeaks,xpix, ypix, numthreads)\n");
	}

	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
		printf("i1 must be a nxnxm double array\n");
		mexErrMsgTxt("See Error above!\n");

	}
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
		printf("beta0 must be a 1x4xnumel(i1(1,1,:)) double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
		printf("framenum_temp must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])){
		printf("xpeaks must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])){
		printf("ypeaks must be a numel(i1(1,1,:))x1 double array\n");
		mexErrMsgTxt("See Error above!\n");
	}
	
	// get pointer to input arguments
	iall = (double *)mxGetPr(prhs[0]);		// matlab linearizes in a coloumn major format which affects indexing (Writing MAtlab C/MEX Code - Research Gate)
	idims = mxGetDimensions(prhs[0]);	// get dimensions of image array
	icol = (int)idims[0];
	irow = (int)idims[1];
	numi = (int)idims[2];
	if (numi > 10000000 || numi < 1){
		numi = 1;
	}
	beta0 = (double *)mxGetPr(prhs[1]);
	framenum_temp = (double *)mxGetPr(prhs[2]);
	xpeaks = (double *)mxGetPr(prhs[3]);
	ypeaks = (double *)mxGetPr(prhs[4]);
	xpix = (double *)mxGetPr(prhs[5]);
	ypix = (double *)mxGetPr(prhs[6]);
	
	if (mxGetM(prhs[1]) != 5){		// ensure that beta0 is a 5 rowed matrix
		printf("Your beta0 vector must be of the format beta0(:,n) if referencing the nth molecule estimators.\n");
		printf("beta0 should be of the form \nbeta0(1,:) = x guesses\nbeta0(2,:) = y guesses\nbeta0(3,:) = N guesses\nbeta0(4,:) = sigma guesses\nbeta0(5,:) = offset guesses\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetN(prhs[1]) != numi){
		printf("Your beta0 vector must have the dimensions 5 rows x numel(i2(1,1,:)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[2]) != numi && mxGetN(prhs[2]) != numi){
		printf("Your framenum_temp vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[3]) != numi && mxGetN(prhs[3]) != numi){
		printf("Your xpeaks vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[4]) != numi && mxGetN(prhs[4]) != numi){
		printf("Your ypeaks vector must have numel(i2(1,1,:)) elements\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[5]) != irow || mxGetN(prhs[5]) != irow){
		printf("Your  xpix matrix must have the dimensions: numel(i2(:,1,1)) rows x numel(i2(1,:,1)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}
	if (mxGetM(prhs[6]) != irow || mxGetN(prhs[6]) != irow){
		printf("Your  ypix matrix must have the dimensions: numel(i2(:,1,1)) rows x numel(i2(1,:,1)) columns\n");
		mexErrMsgTxt("See Error above!\n");
	}

	// verify that the input variables are what was expected
	if (irow == 1 || icol == 1){
		printf("Images are too small must have more than 1 pixel in each dimension\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (irow != icol){			// insure inputs are square matrices
		printf("Images must have equal length and width\n");
		mexErrMsgTxt("See Error above!\n");
	}

	if (nlhs != 10){
		printf("You must have 10 output variables [xf_all, yf_all, N, off_all, xf_crlb, yf_crlb, N_crlb, off_crlb, llv_all, framenum_all]\n");
		mexErrMsgTxt("See Error above!\n");
	}
	// allocate memory on the gpu device

	cudaError_t errxp = cudaMalloc((void**)&d_xpix, irow*irow*sizeof(double));						// allocate xpix memory
	if (errxp != cudaSuccess){
		printf("%s in \n%s \nat line %d\n", cudaGetErrorString(errxp), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erryp = cudaMalloc((void**)&d_ypix, irow*irow*sizeof(double));						// allocate ypix memory
	if (erryp != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erryp), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err1 = cudaMalloc((void**)&d_iall, irow*irow*numi*sizeof(double));				// allocate image memory
	if (err1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err2 = cudaMalloc((void**)&d_beta0, 5 * numi*sizeof(double));						// allocate parameter memory
	if (err2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err3 = cudaMalloc((void**)&d_x_cm, numi*sizeof(double));						// allocate x center coords
	if (err3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err4 = cudaMalloc((void**)&d_y_cm, numi*sizeof(double));						// allocate y center coords
	if (err4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err5 = cudaMalloc((void**)&d_N, numi*sizeof(double));							// allocate N array on gpu
	if (err5 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err5), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err6 = cudaMalloc((void**)&d_off, numi*sizeof(double));						// allocate offset array on gpu
	if (err6 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err6), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err7 = cudaMalloc((void**)&d_yf_all, numi*sizeof(double));						// allocate yf_all array on gpu
	if (err7 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err7), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err8 = cudaMalloc((void**)&d_xf_all, numi*sizeof(double));						// allocate xf_all array on gpu
	if (err8 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err8), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err14 = cudaMalloc((void**)&d_framenum_temp, numi*sizeof(double));				// allocate framenum_temp array on gpu
	if (err14 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err14), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err15 = cudaMalloc((void**)&d_framenum_all, numi*sizeof(double));				// allocate framenum_all array on gpu
	if (err15 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err15), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erra1 = cudaMalloc((void**)&d_xf_crlb, numi*sizeof(double));					// Allocate xf_crlb array on gpu
	if (erra1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erra1), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erra2 = cudaMalloc((void**)&d_yf_crlb, numi*sizeof(double));					// allocate yf_crlb array on gpu
	if (erra2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erra2), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erra3 = cudaMalloc((void**)&d_N_crlb, numi*sizeof(double));					// allocate N_crlb array on gpu
	if (erra3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erra3), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erra4 = cudaMalloc((void**)&d_off_crlb, numi*sizeof(double));					// allocate Off_crlb array on gpu
	if (erra4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erra4), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t errll = cudaMalloc((void**)&d_llv, numi*sizeof(double));					// allocate llv array on gpu
	if (errll != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errll), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	// copy data from host to device
	cudaError_t err9 = cudaMemcpy(d_iall, iall, icol*irow*numi*sizeof(double), cudaMemcpyHostToDevice);	// copy image data to gpu
	if (err9 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err9), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err10 = cudaMemcpy(d_beta0, beta0, 5 * numi*sizeof(double), cudaMemcpyHostToDevice);		// copy parameter data to gpu
	if (err10 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err10), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err11 = cudaMemcpy(d_x_cm, xpeaks, numi*sizeof(double), cudaMemcpyHostToDevice);		// copy x center coords to gpu
	if (err11 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err11), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err12 = cudaMemcpy(d_y_cm, ypeaks, numi*sizeof(double), cudaMemcpyHostToDevice);		// copy y center coords to gpu
	if (err12 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err12), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t err13 = cudaMemcpy(d_framenum_temp, framenum_temp, numi*sizeof(double), cudaMemcpyHostToDevice);					// copy framenum_temp array on gpu
	if (err13 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err13), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t errxc = cudaMemcpy(d_xpix, xpix, irow*irow*sizeof(double), cudaMemcpyHostToDevice);
	if (errxc != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errxc), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	cudaError_t erryc = cudaMemcpy(d_ypix, ypix, irow*irow*sizeof(double), cudaMemcpyHostToDevice);
	if (erryc != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(erryc), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}
	/* Run GPU kernel*/
	threadsperblock = mxGetScalar(prhs[7]);  // get number of threads perblock from matlab

	switch (irow)		// this switch statement allows us to choose the correct function for the size of the images. Sizes outside what is provided seem unreasonable for imaging
	{
	case 7: 
		localize7 << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_beta0, d_x_cm, d_y_cm, d_framenum_temp, d_xf_all, d_yf_all, d_N, d_off, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, irow, d_xpix, d_ypix, d_llv, numi);
		break;
	case 9:
		localize9 << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_beta0, d_x_cm, d_y_cm, d_framenum_temp, d_xf_all, d_yf_all, d_N, d_off, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, irow, d_xpix, d_ypix, d_llv, numi);
		break;
	case 11:
		localize11 << <((numi - 1) / threadsperblock + 1), threadsperblock >> >(d_iall, d_beta0, d_x_cm, d_y_cm, d_framenum_temp, d_xf_all, d_yf_all, d_N, d_off, d_framenum_all, d_xf_crlb, d_yf_crlb, d_N_crlb, d_off_crlb, irow, d_xpix, d_ypix, d_llv, numi);
		break;
	default:
		printf("Image size is inappropriate please choose either 7x7, 9x9, or 11x11 size\n");
		mexErrMsgTxt("See Error Above!\n");
		break;
	}
	/*		 copy data back to mxarray pointers for output
	*
	*
	*		Duplicate the input array of equal size to the output array
	*		Send the pointer to a variable
	*		copy data to place pointer points to, which is output
	*/


	plhs[0] = mxDuplicateArray( prhs[3]);
	mxArray *xf_all = (mxArray *)mxGetPr(plhs[0]);
	cudaError_t err16 = cudaMemcpy( xf_all, d_xf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_all data
	if (err16 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err16), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}
	
	plhs[1] = mxDuplicateArray(prhs[3]);
	double *yf_all = (double *)mxGetPr(plhs[1]);
	cudaError_t err17 = cudaMemcpy(yf_all, d_yf_all, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_all data
	if (err17 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err17), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[3] = mxDuplicateArray(prhs[3]);
	double *off_all = (double *)mxGetPr(plhs[3]);
	cudaError_t err19 = cudaMemcpy(off_all, d_off, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy off_all data
	if (err19 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err19), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[2] = mxDuplicateArray(prhs[3]);
	double *N = (double *)mxGetPr(plhs[2]);
	cudaError_t err18 = cudaMemcpy(N, d_N, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N data
	if (err18 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err18), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[4] = mxDuplicateArray(prhs[3]);
	double *xf_crlb = (double *)mxGetPr(plhs[4]);
	cudaError_t errb1 = cudaMemcpy(xf_crlb, d_xf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy xf_crlb data
	if (errb1 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb1), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	};

	plhs[5] = mxDuplicateArray(prhs[3]);
	double *yf_crlb = (double *)mxGetPr(plhs[5]);
	cudaError_t errb2 = cudaMemcpy(yf_crlb, d_yf_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);	// copy yf_crlb data
	if (errb2 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb2), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[6] = mxDuplicateArray(prhs[3]);
	double *N_crlb = (double *)mxGetPr(plhs[6]);
	cudaError_t errb3 = cudaMemcpy(N_crlb, d_N_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N_crlb data
	if (errb3 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb3), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[7] = mxDuplicateArray(prhs[3]);
	double *off_crlb = (double *)mxGetPr(plhs[7]);
	cudaError_t errb4 = cudaMemcpy(off_crlb, d_off_crlb, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy Off_crlb data
	if (errb4 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb4), __FILE__, __LINE__);
	mexErrMsgTxt("See Error above!\n");
	}

	plhs[8] = mxDuplicateArray(prhs[3]);
	double *llv = (double *)mxGetPr(plhs[8]);
	cudaError_t errb5 = cudaMemcpy(llv, d_llv, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N_crlb data
	if (errb5 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb5), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	plhs[9] = mxDuplicateArray(prhs[3]);
	double *framenum_all = (double *)mxGetPr(plhs[9]);
	cudaError_t errb6 = cudaMemcpy(framenum_all, d_framenum_all, numi*sizeof(double), cudaMemcpyDeviceToHost);		// copy N_crlb data
	if (errb6 != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(errb6), __FILE__, __LINE__);
		mexErrMsgTxt("See Error above!\n");
	}

	//cudaDeviceReset();
	
	cudaFree(d_iall);
	cudaFree(d_beta0); 
	cudaFree(d_N); 
	cudaFree(d_framenum_temp); 
	cudaFree(d_x_cm); 
	cudaFree(d_y_cm);
	cudaFree(d_framenum_all); 
	cudaFree(d_xf_all); 
	cudaFree(d_yf_all); 
	cudaFree(d_off);
	cudaFree(d_xpix);
	cudaFree(d_ypix);
	
}  
