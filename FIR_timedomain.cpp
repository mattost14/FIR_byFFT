#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>

// Filter coefficients from header file
#include "fdacoefs.h"

using namespace std;

double u[BL] = {0.0};


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
    

  cout << "Start FIR filtering in the time domain" << endl;
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  FILE *fx, *fy;
  fx = fopen(inputFileName,"rb"); // argv[1]: Input signal file name
  if(fx == NULL){
    cerr << "Input file not found" << endl;
    return -1;
  }
  fy = fopen(outputFileName,"wb"); // argv[2]: Output signal file name


  float x;
  float y;

  while(!feof(fx)){
    
    fread(&x, sizeof(float), 1, fx); // read in next sample
    y = filter(x); 
    fwrite(&y, sizeof(float), 1, fy); // save output
  }
  chrono::steady_clock::time_point end = chrono::steady_clock::now();

  // Close all files
  fclose(fx);
  fclose(fy);

  cout << "Completed in " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;

  return 0;
}