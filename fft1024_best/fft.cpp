/*
This is traditional 2-radix DIT FFT algorithm implementation.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	Out_R, Out_I[]: Real and Imag parts of Complex signal
*/

#include "fft.h"

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]);
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int STAGES, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);

// function to return bit reserved value
unsigned int reverse_bits(unsigned int n){
	int rev = 0;
	for (int i = 0; i < M; i++) {
		rev = (rev << 1) |(n & 1);
		n = n >> 1;
	}
	return rev;
}

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{
	bit_reverse(X_R, X_I);

	//Call fft
	DTYPE Stage1_R[SIZE], Stage1_I[SIZE];
	DTYPE Stage2_R[SIZE], Stage2_I[SIZE];
	DTYPE Stage3_R[SIZE], Stage3_I[SIZE];
	DTYPE Stage4_R[SIZE], Stage4_I[SIZE];
	DTYPE Stage5_R[SIZE], Stage5_I[SIZE];
	DTYPE Stage6_R[SIZE], Stage6_I[SIZE];
	DTYPE Stage7_R[SIZE], Stage7_I[SIZE];
	DTYPE Stage8_R[SIZE], Stage8_I[SIZE];
	DTYPE Stage9_R[SIZE], Stage9_I[SIZE];

	fft_stage_first(X_R, X_I, Stage1_R, Stage1_I);
	fft_stages(Stage1_R, Stage1_I, 2, Stage2_R, Stage2_I);
	fft_stages(Stage2_R, Stage2_I, 3, Stage3_R, Stage3_I);
	fft_stages(Stage3_R, Stage3_I, 4, Stage4_R, Stage4_I);
	fft_stages(Stage4_R, Stage4_I, 5, Stage5_R, Stage5_I);
	fft_stages(Stage5_R, Stage5_I, 6, Stage6_R, Stage6_I);
	fft_stages(Stage6_R, Stage6_I, 7, Stage7_R, Stage7_I);
	fft_stages(Stage7_R, Stage7_I, 8, Stage8_R, Stage8_I);
	fft_stages(Stage8_R, Stage8_I, 9, Stage9_R, Stage9_I);
	fft_stage_last(Stage9_R, Stage9_I, OUT_R, OUT_I);

}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]){
	unsigned int reversed_r;
	unsigned int reversed_i;

	DTYPE temp;
	for (unsigned int i = 0; i < SIZE; i++) {
		reversed_r = reverse_bits(X_R[i]); // Find the bit reversed index
		reversed_i = reverse_bits(X_I[i]); // Find the bit reversed index
		if (i <= reversed_r ) {
			// Swap the real values
			temp = X_R[i];
			X_R[i] = X_R[reversed_r];
			X_R[reversed_r] = temp;
		}
		if(i <= reversed_i){
			// Swap the imaginary values
			temp = X_I[i];
			X_I[i] = X_I[reversed_i];
			X_I[reversed_i] = temp;
		}
	}


}
/*=======================BEGIN: FFT=========================*/
//stage 1
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	//int step;

	int stage =1;
	int DFTpts = 1 << stage;
	int numBF = DFTpts/2;			/*Butterfly Width*/

	DTYPE c,s;

	for(int i = 0; i < SIZE; i += DFTpts){
		i_lower = i + numBF;			//index of lower point in butterfly

		c = W_real[0];					//twiddle (stage 1 uses only 1 twiddle factor)
		s = W_imag[0];					//twiddle (stage 1 uses only 1 twiddle factor)

		temp_R = X_R[i_lower]*c - X_I[i_lower]*s;
		temp_I = X_I[i_lower]*c + X_R[i_lower]*s;

		OUT_R[i_lower] = X_R[i] - temp_R;
		OUT_I[i_lower] = X_I[i] - temp_I;
		OUT_R[i] = X_R[i] + temp_R;
		OUT_I[i] = X_I[i] + temp_I;

	}


}

//stages
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i,j,k;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	int step;
	int index;			/* Twiddle Index*/

	int locStage = stage;     /* local copy of stage value */
	int DFTpts;			/* points in sub DFT */
	int numBF;			/*Butterfly Width*/

	int N2 = SIZE2;	/* N2=N>>1 */

	DFTpts = 1 << stage; /*points in sub DFT */
	numBF = DFTpts/2;       /* buttery widths in sub DFT */
	DTYPE c, s;
	index = 0; /*index of twiddle values, start at 0*/
	int twiddleFactor = SIZE/DFTpts;   /* used to calculate unique twiddle value */

	// Perform butterflies for j-th stage
		for(j = 0; j < numBF; j++){

			c = W_real[index]; 		/* twiddle factor */
			s = W_imag[index];		/* twiddle factor */

			// Compute butterflies that use same W**k
			for(i=j; i<SIZE; i += DFTpts){

				i_lower = i + numBF;			//index of lower point in butterfly

				temp_R = X_R[i_lower]*c - X_I[i_lower]*s;
				temp_I = X_I[i_lower]*c + X_R[i_lower]*s;

				OUT_R[i_lower] = X_R[i] - temp_R;
				OUT_I[i_lower] = X_I[i] - temp_I;
				OUT_R[i] = X_R[i] + temp_R;
				OUT_I[i] = X_I[i] + temp_I;

				}
			index += twiddleFactor;   /* update index to next unique twiddle val */
		}


}


//last stage
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/
	DTYPE c, s;


	int i;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	//int step;

	//int stage;
	//int DFTpts;
	//int numBF;			/*Butterfly Width*/

	for(int i = 0; i<SIZE2; i ++){ //loop will iterate Size/2 times

		i_lower = i + int(SIZE2); // lower butterfly index (not sure why doesnt work without type cast)

		c = W_real[i];		//twiddle factors index = {0, 1, 2, ...etc)
		s = W_imag[i];		//twiddler factor index = {0, 1, 2, ...etc)

		temp_R = X_R[i_lower]*c - X_I[i_lower]*s;
		temp_I = X_I[i_lower]*c + X_R[i_lower]*s;

		OUT_R[i_lower] = X_R[i] - temp_R;
		OUT_I[i_lower] = X_I[i] - temp_I;
		OUT_R[i] = X_R[i] + temp_R;
		OUT_I[i] = X_I[i] + temp_I;


	}


}
/*=======================END: FFT=========================*/
