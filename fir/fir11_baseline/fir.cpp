/*
	Filename: fir.cpp
		FIR lab wirtten for WES/CSE237C class at UCSD.

	INPUT:
		x: signal (chirp)

	OUTPUT:
		y: filtered output

*/

#include "fir.h"

void fir (
  data_t *y,
  data_t x
  )
{
	#pragma HLS INTERFACE mode=s_axilite port=return
	#pragma HLS INTERFACE mode=s_axilite port=x
	#pragma HLS INTERFACE mode=s_axilite port=y
	coef_t c[N] = {53, 0, -91, 0, 313, 500, 313, 0, -91, 0,53};
	// Write your code here
	static data_t shiftReg[N];
	acc_t accum = 0;
	int i;

	for( i = N -1; i >= 0; i--){
		if(i == 0){
			accum += x * c[0];
			shiftReg[0] = x;
		}
		else{
			shiftReg[i] = shiftReg[i-1];
			accum += shiftReg[i]*c[i];
		}
	}
	*y = accum;
}


