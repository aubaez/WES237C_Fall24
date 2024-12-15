#include <hls_stream.h>
#include <ap_fixed.h>
#include <ap_axi_sdata.h>


typedef float DTYPE;
typedef ap_axis<32,2,5,6> AXI_T; //define axi stream type (model after lab def)
typedef hls::stream<AXI_T> STREAM_T; //define stream

union fp_int { //define float, int union
	int i;
	float f;
};

#define SIZE 256 		/* SIZE OF DFT */

//updated dft function for stream
void dft(STREAM_T &real_sample, STREAM_T &imag_sample, STREAM_T &Out_R,STREAM_T &Out_I);

