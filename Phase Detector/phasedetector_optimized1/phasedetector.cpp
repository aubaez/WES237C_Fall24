/*
	Filename: phasedetector.cpp
		Phase detector

	INPUT:
		I: signal for I sample
		Q: signal for Q sample
		length: array size

	OUTPUT:
		R: Radius
		Theta: Angle

*/

#include "phasedetector.h"

void phasedetector (
  data_t *I,
  data_t *Q,

  data_t *R,
  data_t *Theta,

  int length
  ){

	// Write your code here
	data_t x, y;
	for(int i = 0; i < length; i++){
		fir(I[i],Q[i],&x,&y);
		cordiccart2pol(x,y,&R[i],&Theta[i]);
	}

}




