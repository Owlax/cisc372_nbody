#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the d_hPos and d_hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(vector3* d_hPos, vector3* d_hVel, double* d_mass, vector3** accels, vector3* values){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	//first compute the pairwise accelerations.  Effect is on the first argument.
	//start parallel here: have each thread compute how two objects affect eachother and update the matrix
	//something like set i and j to the two dimensions of the resulting accel matrix and have one thread for each pair
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	int row = t/NUMENTITIES;
    int col = t%NUMENTITIES;

	if (t < NUMENTITIES * NUMENTITIES) {
        accels[row]=&values[row*NUMENTITIES];	
		//accels[row][col] = values[row * NUMENTITIES + col];

		if (row==col) {
			FILL_VECTOR(accels[row][col],0,0,0);
		}
		else{
			vector3 distance;
			for (int k=0;k<3;k++) distance[k]=d_hPos[row][k]-d_hPos[col][k];
			double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
			double magnitude=sqrt(magnitude_sq);
			double accelmag=-1*GRAV_CONSTANT*d_mass[col]/magnitude_sq;
			FILL_VECTOR(accels[row][col],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
		}
	}
}


__global__ void compute2electricboogaloo(vector3* d_hPos, vector3* d_hVel, double* d_mass, vector3** accels, vector3* values){
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	int row = t/NUMENTITIES;
    int col = t%NUMENTITIES;

	if(t<NUMENTITIES*NUMENTITIES){
		vector3 accel_sum={0,0,0};
		for (int k=0;k<3;k++){
			accel_sum[k]+=accels[row][col][k];
		}
		//compute the new velocity based on the acceleration and time interval
		//compute the new position based on the velocity and time interval
		for (int k=0;k<3;k++){
			d_hVel[row][k]+=accel_sum[k]*INTERVAL;
			d_hPos[row][k]+=d_hVel[row][k]*INTERVAL;
		}
	}
}