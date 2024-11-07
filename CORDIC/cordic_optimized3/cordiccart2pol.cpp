#include "cordiccart2pol.h"
#include "math.h"


data_t Kvalues[NO_ITER] = {1,	0.500000000000000,	0.250000000000000,	0.125000000000000,	0.0625000000000000,	0.0312500000000000,	0.0156250000000000,	0.00781250000000000,	0.00390625000000000,	0.00195312500000000,	0.000976562500000000,	0.000488281250000000,	0.000244140625000000,	0.000122070312500000,	6.10351562500000e-05,	3.05175781250000e-05};

data_t angles[NO_ITER] = {0.785398163397448,	0.463647609000806,	0.244978663126864,	0.124354994546761,	0.0624188099959574,	0.0312398334302683,	0.0156237286204768,	0.00781234106010111,	0.00390623013196697,	0.00195312251647882,	0.000976562189559320,	0.000488281211194898,	0.000244140620149362,	0.000122070311893670,	6.10351561742088e-05,	3.05175781155261e-05};

// define gain as global variable
data_t gain = 0.6072529350088812561694;

void cordiccart2pol(data_t x, data_t y, data_t * r,  data_t * theta)
{
	// Write your code here

	// define parameters for doing calculation
	my_t x2;
	data_t angle = 0.0;

	// local copy of x and y input (stored as fixed point)
	my_t temp_x = x;
	my_t temp_y = y;


	// check for signs of input to determine quadrant location and adjust angle/coordinate vals
	if( temp_x < 0){
		if(temp_y > 0){
			angle = 3.1415926535;
		}
		else if( temp_y < 0){
			angle = -3.1415926535;
		}
		temp_x = -temp_x;
		temp_y = -temp_y;
	}


	// CORDIC algorithm
	for(int i = 0; i < NO_ITER; i++){
		x2 = temp_x;
		if(temp_y < 0) {
			temp_x -= temp_y >> i;
			temp_y += x2 >> i;
			angle -=  angles[i];
		}
		else{
			temp_x += temp_y >> i;
			temp_y -= x2 >> i;
			angle = angle + angles[i];
		}

	}
	//assign calculated values of r and theta
	*r = static_cast<data_t>(temp_x)*gain;
	*theta = angle;
}
