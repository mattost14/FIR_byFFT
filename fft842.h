#pragma once

typedef struct {
	float re;
	float im;
}complx;

void fft842(int in, int n, complx *x);