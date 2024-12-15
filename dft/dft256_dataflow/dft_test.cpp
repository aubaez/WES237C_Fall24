/*
This is DFT computation using matrix vector multiplication form.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal in time domain.
OUTPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal in frequency domain.

*/

#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include <math.h>
#include "dft.h"

struct Rmse
{
	int num_sq;
	float sum_sq;
	float error;

	Rmse(){ num_sq = 0; sum_sq = 0; error = 0; }

	float add_value(float d_n)
	{
		num_sq++;
		sum_sq += (d_n*d_n);
		error = sqrtf(sum_sq / num_sq);
		return error;
	}

};

Rmse rmse_R,  rmse_I;
fp_int In_R[SIZE], In_I[SIZE], OUT_R[SIZE], OUT_I[SIZE];

int main()
{
	int index;
	float gold_R, gold_I;

	FILE * fp = fopen("out.gold.dat","r");
	STREAM_T real_sample, imag_sample, Out_R,Out_I;
	AXI_T re,im;

	for(int i=0; i<SIZE; i++)
	{
		In_R[i].f = i;
		In_I[i].f = 0.0;

		re.data = In_R[i].i;
		im.data = In_I[i].i;
		if(i==SIZE-1){
			re.last = 1;
			im.last = 1;
		}
		else{
			re.last = 0;
			im.last = 0;
		}


		real_sample.write(re);
		imag_sample.write(im);
	}
	

	// DFT
	dft(real_sample, imag_sample, Out_R, Out_I);


	// comparing with golden output
	for(int i=0; i<SIZE; i++)
	{
		re = Out_R.read();
		im = Out_I.read();

		OUT_R[i].i = re.data;
		OUT_I[i].i = im.data;

		fscanf(fp, "%d %f %f", &index, &gold_R, &gold_I);
		//printf("%0.15f %0.15f %0.15f %0.15f\n",OUT_R[i].fp, OUT_I[i].fp, gold_R, gold_I);
		rmse_R.add_value((float)OUT_R[i].f - gold_R);
		rmse_I.add_value((float)OUT_I[i].f - gold_I);
	}
	fclose(fp);


	// printing error results
	printf("----------------------------------------------\n");
	printf("   RMSE(R)           RMSE(I)\n");
	printf("%0.15f %0.15f\n", rmse_R.error, rmse_I.error);
	printf("----------------------------------------------\n");

	if (rmse_R.error > 0.1 || rmse_I.error > 0.1 ) {
		fprintf(stdout, "*******************************************\n");
		fprintf(stdout, "FAIL: Output DOES NOT match the golden output\n");
		fprintf(stdout, "*******************************************\n");
	    return 1;
	}else {
		fprintf(stdout, "*******************************************\n");
		fprintf(stdout, "PASS: The output matches the golden output!\n");
		fprintf(stdout, "*******************************************\n");
	    return 0;
	}

}
