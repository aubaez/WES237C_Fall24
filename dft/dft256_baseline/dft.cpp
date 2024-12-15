#include<math.h>
#include "dft.h"
#include"coefficients256.h"

#define PI 3.141592653589

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE])
{
	//declare variables for omega, cosine value and sine value
	DTYPE w,cosine, sine;

	//declare arrays to hold temporary frequency result
	DTYPE tempReal[SIZE], tempImag[SIZE];

	for(int i = 0; i < SIZE; i++){
		//initialize values to 0
		tempReal[i] = 0.0;
		tempImag[i] = 0.0;

		// calculate value for omega
		w = (2.0*3.141592653589/SIZE) * static_cast<DTYPE>(i);

		// calculation for frequency sample at jth position
		for(int j = 0; j < SIZE; j++){
			cosine = cos(j*w);
			sine = -sin(j*w);

			// perform summation from 0 to SIZE-1 of frequency sample * input sample
			tempReal[i] += (real_sample[j]*cosine - imag_sample[j]*sine);
			tempImag[i] += (real_sample[j]*sine + imag_sample[j]*cosine);
		}
	}
	// copy over values from temp to input arrays
	for(int k = 0; k < SIZE; k++){
		real_sample[k] = tempReal[k];
		imag_sample[k] = tempImag[k];
	}
}
