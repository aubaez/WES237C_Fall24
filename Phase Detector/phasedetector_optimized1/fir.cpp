/*
	Filename: fir.cpp
		Complex FIR or Match filter
		firI1 and firI2 share coef_t c[N]
		firQ1 and firQ2 share coef_t c[N]
		
	INPUT:
		I: signal for I sample
		I: signal for Q sample

	OUTPUT:
		X: filtered output
		Y: filtered output

*/

#include "phasedetector.h"

void firI1 (
  data_t *y,
  data_t x
  ) {

	coef_t c[N] = {1,    -1,    1,    -1,    -1,    -1,    1,    1,    -1,    -1,    -1,    1,    1,    -1,    1,    -1,    -1,    -1,    -1,    1,    1,    1,    1,    1,    -1,    -1,    1,    1,    1,    -1,    -1,    -1};

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


void firI2 (
  data_t *y,
  data_t x
  ) {

	coef_t c[N] = {1,    -1,    1,    -1,    -1,    -1,    1,    1,    -1,    -1,    -1,    1,    1,    -1,    1,    -1,    -1,    -1,    -1,    1,    1,    1,    1,    1,    -1,    -1,    1,    1,    1,    -1,    -1,    -1};

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




void firQ1 (
  data_t *y,
  data_t x
  ) {

	coef_t c[N] = {-1,    -1,    1,    -1,    1,    -1,    1,    -1,    -1,    -1,    -1,    1,    -1,    1,    -1,    1,    1,    -1,    1,    -1,    -1,    1,    -1,    1,    1,    1,    1,    -1,    1,    -1,    1,    1};


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

void firQ2 (
  data_t *y,
  data_t x
  ) {

	coef_t c[N] = {-1,    -1,    1,    -1,    1,    -1,    1,    -1,    -1,    -1,    -1,    1,    -1,    1,    -1,    1,    1,    -1,    1,    -1,    -1,    1,    -1,    1,    1,    1,    1,    -1,    1,    -1,    1,    1};


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


void fir (
  data_t I,
  data_t Q,

  data_t *X,
  data_t *Y
  ) {

	// Write your code here
	data_t Iin_Ifir, Qin_Qfir, Qin_Ifir, Iin_Qfir;

	firI1(&Iin_Ifir, I);
	firQ1(&Qin_Qfir,Q);
	firI2(&Qin_Ifir,Q);
	firQ2(&Iin_Qfir,I);
	

	//Calculate X
	*X = Iin_Ifir + Qin_Qfir;
	//Calculate Y
	*Y = Qin_Ifir - Iin_Qfir;


}


