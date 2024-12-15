#include<math.h>
#include "dft.h"
#include"coefficients256.h"

#define PI 3.141592653589

void dft(DTYPE real_sample[SIZE], DTYPE imag_sample[SIZE], DTYPE Out_R[SIZE], DTYPE Out_I[SIZE])
{
	//declare variables for omega, cosine value and sine value
	DTYPE cosine, sine;

	//declare arrays to hold temporary frequency result
	DTYPE tempReal[SIZE], tempImag[SIZE];

	for(int i = 0; i < SIZE; i++){
		//initialize values to 0
		tempReal[i] = 0.0;
		tempImag[i] = 0.0;


		// calculation for frequency sample as matrix multiplication
		// index of coefficient calculated as (column i * row j) % SIZE
		for(int j = 0; j < SIZE; j++){
			cosine = cos_coefficients_table[i*j%SIZE];
			sine = sin_coefficients_table[i*j%SIZE];

			// perform summation from 0 to SIZE-1 of frequency sample * input sample
			tempReal[i] += (real_sample[j]*cosine - imag_sample[j]*sine);
			tempImag[i] += (real_sample[j]*sine + imag_sample[j]*cosine);
		}
	}
	// copy over values from temp to input arrays
	for(int k = 0; k < SIZE; k++){
		Out_R[k] = tempReal[k];
		Out_I[k] = tempImag[k];
	}
}
