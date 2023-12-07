#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the d_hPos and d_hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(vector3* d_hPos, vector3* d_hVel, double* d_mass){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int k;
	//first compute the pairwise accelerations.  Effect is on the first argument.
	//start parallel here: have each thread compute how two objects affect eachother and update the matrix
	//something like set i and j to the two dimensions of the resulting accel matrix and have one thread for each pair
	int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

	//TODO: add a check that row*col<=NUMENTITIES
	if (row==col) {
		FILL_VECTOR(accels[row][col],0,0,0);
	}
	else{
		vector3 distance;
		for (k=0;k<3;k++) distance[k]=d_hPos[row][k]-d_hPos[col][k];
		double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
		double magnitude=sqrt(magnitude_sq);
		double accelmag=-1*GRAV_CONSTANT*d_mass[col]/magnitude_sq;
		FILL_VECTOR(accels[row][col],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
	}

	// Synchronize all threads in the block
    __syncthreads();

	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	//sync threads
	//have each thread add up one collumn
	vector3 accel_sum={0,0,0};
	for (k=0;k<3;k++){
		accel_sum[k]+=accels[row][col][k];
	}
	//compute the new velocity based on the acceleration and time interval
	//compute the new position based on the velocity and time interval
	for (k=0;k<3;k++){
		d_hVel[row][k]+=accel_sum[k]*INTERVAL;
		d_hPos[row][k]+=d_hVel[row][k]*INTERVAL;
	}

	//I could do another sync threads and then have each thread compute one velocity and one position

	free(accels);
	free(values);
}
