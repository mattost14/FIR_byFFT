#include <iostream>
#include <fstream> //Stream class to both read and write from/to files.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>
#include <cmath>

#include "fft842.h"

// Filter coefficients from header file
#include "fdacoefs.h"

#define FFT_SIZE 512

using namespace std;

// Vector that holds the filter FFT 
float H_real[FFT_SIZE] = {0.0};
float H_imag[FFT_SIZE] = {0.0};


double u[BL] = {0.0};

// This buffer is designed to hold the input samples for the overlap-add implementation
double buffer[FFT_SIZE - BL + 1] = {0.0};


double filter(float x) {
    // Shift the register values.
    for(int k=BL-1; k>0; k--){
      u[k] = u[k-1];
    }
    u[0] = x;
    // The numerator
    double y = 0;
    for(int k=0; k<=BL-1; k++){
      y += B[k] * u[k];
    }
    return y;
}

int main (int argc, char *argv[]) {
  if (argc < 3)	{
		cerr << "Usage: " << argv[0] << " <input filename>  <output filename>\n";
		return -1;
	}
  
  cout << "input file = " << argv[1] << endl;
  char* inputFileName = argv[1];
  cout << "output file = " << argv[2] << endl;
  char* outputFileName = argv[2];
    

  cout << "Start FIR filtering using 512-point FFT" << endl;
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  // FILE *fx, *fy;
  ofstream fx, fy;
  fx.open (inputFileName, ios::in);
  // fx = fopen(inputFileName,"rb"); // argv[1]: Input signal file name
  if(!fx.is_open()){
    cerr << "Input file not found" << endl;
    return -1;
  }
  fy.open (outputFileName, ios::out);
  // fy = fopen(outputFileName,"wb"); // argv[2]: Output signal file name
  // fy2 = fopen("out_im.bin","wb"); // argv[2]: Output signal file name

  // Perform the FFT of the filter h[n] once
  copy(B,B+BL, H_real);
  // copy(H_real, B, BL*sizeof(real64_T));
  complx H[FFT_SIZE];
  for(int i=0; i<FFT_SIZE; i++){
    H[i].im = H_imag[i];
    H[i].re = H_real[i];
  }
  fft842(0,FFT_SIZE, H);
  for(int i=0; i<FFT_SIZE; i++){
    float mag = sqrt(pow((H[i].re),2) + pow((H[i].im),2));
    fy << mag;
  //  fwrite(*(sqrt(pow((H[i].re),2) + pow((H[i].im),2))), sizeof(float), 1, fy); // save output
  //  fwrite(&H[i].im, sizeof(float), 1, fy2); // save output
  }
  


  // float x;
  // float y;
  // while(!feof(fx)){
    
  //   fread(&x, sizeof(float), 1, fx); // read in next sample
  //   y = filter(x); 
  //   fwrite(&y, sizeof(float), 1, fy); // save output
  // }
  // chrono::steady_clock::time_point end = chrono::steady_clock::now();

  // Close all files
fx.close();
fy.close();

  // cout << "Completed in " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;

  return 0;
}