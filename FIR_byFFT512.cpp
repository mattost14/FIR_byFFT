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
int const STAGE_SIZE =  FFT_SIZE - BL + 1;

using namespace std;

// Vector that holds the filter FFT 
float H_real[FFT_SIZE] = {0.0};
float H_imag[FFT_SIZE] = {0.0};
// This buffer is designed to hold the input samples for the overlap-add implementation
float buffer[STAGE_SIZE] = {0.0};



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

  fstream fx, fy;
  fx.open(inputFileName, ios::in|ios::binary);
  if(!fx.is_open()){
    cerr << "Input file not found" << endl;
    return -1;
  }
  

  // Perform the FFT of the filter h[n] once
  copy(B,B+BL, H_real);
  complx H[FFT_SIZE];
  for(int i=0; i<FFT_SIZE; i++){
    H[i].im = H_imag[i];
    H[i].re = H_real[i];
  }
  fft842(0,FFT_SIZE, H);

  // Open output file
  fy.open (outputFileName, ios::out);
  complx X1[FFT_SIZE], y_previous[FFT_SIZE];


  int k = 0;
  while(!fx.eof()){

    for(int i=0; i<STAGE_SIZE; i++){
      buffer[i] = 0.0;
    }
    fx.read((char *)buffer, sizeof(buffer));

    for(int i=0; i<FFT_SIZE; i++){
      if(i<STAGE_SIZE){
        X1[i].im = 0.0;
        X1[i].re = buffer[i];
      }
      else {
        X1[i].im = 0.0;
        X1[i].re = 0.0;
      }
    }

    fft842(0,FFT_SIZE, X1);

    // Multiply X1 * H
    for(int i=0; i<FFT_SIZE; i++){
      float real, imag;
      real = X1[i].re*H[i].re - X1[i].im*H[i].im;
      imag = X1[i].im*H[i].re + X1[i].re*H[i].im;
      X1[i].re = real;
      X1[i].im = imag;
    }

    // Take the IFFT
    fft842(1,FFT_SIZE, X1);

    
    if(k!=0){
      // int count = 0;
      float y;
      for(int i=BL-1; i<FFT_SIZE; i++){
          if(i== FFT_SIZE-1){
            int ok=1;
          }
          if(i<STAGE_SIZE){ // No overlap region (print out)
            y = y_previous[i].re;
          }
          else { // Overlap region (add!)
            y = y_previous[i].re + X1[i-STAGE_SIZE].re;
          }
          float* var = &y;

          fy.write((char *)var, sizeof(float));
      }
    }

    // Transfer X1 to X2 and zero out X1 
    for(int i=0; i<FFT_SIZE; i++){
      y_previous[i].re = X1[i].re;
      y_previous[i].im = X1[i].im;
      X1[i].re = 0.0;
      X1[i].im = 0.0;
    }
    k++;
  }

  chrono::steady_clock::time_point end = chrono::steady_clock::now(); 

// Close all files
fx.close();
fy.close();

cout << "Completed in " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;


  return 0;
}