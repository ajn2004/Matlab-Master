#include "mex.h"
#include <math.h>
#include <stdio.h>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *dInputLine;
  double *dOutputPixels;
  unsigned short *usOutputPixels;
  unsigned short *usInputPixels;
  unsigned short *usInputLine;
  unsigned char *ucImage;
  unsigned char *actualLine;
  
  /* to be converted to matlab inputs*/
  double *dAOIHeight;
  double *dAOIWidth;
  double *dAOIStrideLength;
  double *dImageSizeBytes;
  unsigned int aoiHeight = 0;
  unsigned int aoiWidth = 0;
  unsigned int aoiStrideLength = 0;
  unsigned int imageSizeBytes = 0;
  char *charCode;
  size_t n = 0;
  unsigned int c;
  char *input_buf;
  int   buflen, status;
  
  unsigned int iStart = 0;
  unsigned int iPixel = 0;
  unsigned int iInputPixel = 0;
  unsigned int iOutputPixel = 0;
  unsigned int iRowStart = 0;
  unsigned int iInputRowStart = 0;
  unsigned int iOutputRowStart = 0;
  unsigned int iDimension = 0;
  
  unsigned int iRow = 0;
  unsigned int iColumn = 0;
  unsigned int iHeight = 0;
  
  /* Data file details*/
  FILE *ifp;  
  mwSize adims;
  const mwSize *dims;
  mwSignedIndex outputDimensions[2];
  
  dAOIWidth = mxGetPr(prhs[1]);
  dAOIHeight = mxGetPr(prhs[2]);
  dAOIStrideLength = mxGetPr(prhs[3]);
  dImageSizeBytes = mxGetPr(prhs[4]);
  aoiHeight = (unsigned int) *dAOIHeight;
  aoiWidth = (unsigned int) *dAOIWidth;
  aoiStrideLength = (unsigned int) *dAOIStrideLength;
  imageSizeBytes = (unsigned int) *dImageSizeBytes;
  

  outputDimensions[0] = aoiWidth;
  outputDimensions[1] = aoiHeight;
  /*printf ("Got here!\n");*/
  plhs[0] = mxCreateNumericArray(2, outputDimensions, mxDOUBLE_CLASS, mxREAL);

  ucImage = (unsigned char*) mxMalloc(aoiHeight * aoiStrideLength * sizeof(unsigned char));
  usInputPixels = (unsigned short*) mxMalloc(aoiHeight * aoiStrideLength * sizeof(unsigned short));
  usOutputPixels = (unsigned short*) mxMalloc(aoiWidth * aoiHeight * sizeof(unsigned short));
  dOutputPixels = (double*) mxMalloc(aoiWidth * aoiHeight * sizeof(double));
  dOutputPixels = mxGetPr(plhs[0]);
  /* Get the length of the input string. */
  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  /* Allocate memory for input string. */
  input_buf = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[0], input_buf, buflen);
  
  /*printf ("%u %u %u %u\n", aoiWidth, aoiHeight, aoiStrideLength, imageSizeBytes);*/
  
  ifp = fopen (input_buf, "rb");
  if (ifp == NULL) {
      printf("Can't open input file %s\n", input_buf);
  }
  else {
      /* Read data file*/
      /*printf("Reading input file..."); */
      /*fread(&ucImage, aoiWidth * aoiStrideLength * sizeof(unsigned char), 1, ifp);*/
      c = fgetc(ifp); c = fgetc(ifp); c = fgetc(ifp); c = fgetc(ifp); /* Skipping 4 characters*/
      while ((c = fgetc(ifp)) != EOF) {
          ucImage[n] = (unsigned char) c;
          /*printf ("%02x ", ucImage[n]);*/
          n = n + 1;
      }
      /*printf("finished reading %d bytes\n", n); */
      
      /* Convert the byte stream data to unsigned short data */
      for (iRow = 0; iRow < aoiHeight; iRow++) {
          iRowStart = iRow * aoiStrideLength;
          for (iColumn = 0; iColumn < aoiStrideLength; iColumn++) {
              iPixel = iRowStart + iColumn;
              usInputPixels[iPixel] = (unsigned short)ucImage[iPixel];
          }
      }
      
      /* Unpack the 12MonoPacked to unsigned shorts*/
      for (iRow = 0; iRow < aoiHeight; iRow++) {
          iInputRowStart = iRow * aoiStrideLength;
          iOutputRowStart = iRow * aoiWidth;
          for (iColumn = 0; iColumn < aoiWidth; iColumn = iColumn + 2) {
              iInputPixel = iInputRowStart + (iColumn * 3)/2;
              iOutputPixel = iOutputRowStart + iColumn;
              /*printf ("iRow = %d, iColumn = %d\n", iRow, iColumn);*/
              usOutputPixels[iOutputPixel] = (usInputPixels[iInputPixel] << 4) + (usInputPixels[iInputPixel+1] & 0xF);
              /* The below comment line is as per the SDK documetation but it doesn't work the half byte is actually the MSB and full next byte is LSB*/
                usOutputPixels[iOutputPixel+1] = (usInputPixels[iInputPixel+2] << 4) + (usInputPixels[iInputPixel+1] >> 4);
              /*usOutputPixels[iOutputPixel+1] = ((usInputPixels[iInputPixel+1] & 0xF0) << 4) + (usInputPixels[iInputPixel+2]);*/
          }
      }
      
      /* Convert the unsigned short data to double*/
      for (iRow = 0; iRow < aoiHeight; iRow++) {
          iRowStart = iRow * aoiWidth;
          for (iColumn = 0; iColumn < aoiWidth; iColumn++) {
              iPixel = iRowStart + iColumn;
              dOutputPixels[iPixel] = (double)usOutputPixels[iPixel];
          }
      }
  }
  fclose(ifp); 
}
