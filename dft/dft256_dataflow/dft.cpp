#include<math.h>
#include "dft.h"
#include<stdio.h>
#include <stdlib.h>
#include"coefficients256.h"
#include <ap_axi_sdata.h>
#include <hls_stream.h>

#define PI 3.141592653589

void dft(STREAM_T &real_sample, STREAM_T &imag_sample, STREAM_T &Out_R,STREAM_T &Out_I)
{
	#pragma HLS INTERFACE mode=s_axilite port=return
	#pragma HLS INTERFACE mode=axis port=real_sample
	#pragma HLS INTERFACE mode=axis port=imag_sample
	#pragma HLS INTERFACE mode=axis port=Out_R
	#pragma HLS INTERFACE mode=axis port=Out_I


	DTYPE cosine,sine;
	AXI_T real_val,imag_val;
	fp_int temp_real[SIZE], temp_imag[SIZE], real_Out[SIZE],imag_Out[SIZE];
	DTYPE cpy_re[SIZE], cpy_im[SIZE];

	//read input stream and keep a copy
	for (int i = 0; i<SIZE; i++){
		real_val = real_sample.read();
		imag_val = imag_sample.read();

		temp_real[i].i = real_val.data;
		temp_imag[i].i = imag_val.data;

		cpy_re[i] = real_val.keep;
		cpy_im[i] = imag_val.keep;

		real_Out[i].f = 0;
		imag_Out[i].f = 0;

	}

	for(int i=0; i < SIZE; i++ ){
        for(int j=0; j < SIZE; j++){
        	cosine = cos_coefficients_table[i*j%SIZE];
            sine = sin_coefficients_table[i*j%SIZE];

            real_Out[i].f += (temp_real[j].f * cosine) - (temp_imag[j].f * sine);
            imag_Out[i].f += (temp_real[j].f * sine) + (temp_imag[j].f * cosine);
        }
		real_val.data = real_Out[i].i;
		imag_val.data = imag_Out[i].i;

		if(i == SIZE - 1){
			real_val.last = 1;
			imag_val.last = 1;
		}
		else{
			real_val.last = 0;
			imag_val.last = 0;
		}


		real_val.keep = cpy_re[i];
		imag_val.keep = cpy_im[i];

		Out_R.write(real_val);
		Out_I.write(imag_val);

	}

}
