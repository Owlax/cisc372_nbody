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
	int i,j,k;
	/*vector3* values=(vector3*)malloc(sizeof(vector3)*NUMENTITIES*NUMENTITIES);	//table of num^2 vectors
	vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);	//table built in list form
	for (i=0;i<NUMENTITIES;i++)	{												//for each entity:
		accels[i]=&values[i*NUMENTITIES];		//the acceleration of that entity is the address of its row
		//this *can* be parallelized but not necessary right away
	}*/ //not needed
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
		double accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
		FILL_VECTOR(accels[row][col],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
	}

/*
	for (i=0;i<NUMENTITIES;i++){
		for (j=0;j<NUMENTITIES;j++){
			if (i==j) {
				FILL_VECTOR(accels[i][j],0,0,0);
			}
			else{
				vector3 distance;
				for (k=0;k<3;k++) distance[k]=d_hPos[i][k]-d_hPos[j][k];
				double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
				double magnitude=sqrt(magnitude_sq);
				double accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
				FILL_VECTOR(accels[i][j],accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
			}
		}
	}
	*/
	
	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	//sync threads
	//have each thread add up one collumn
	for (i=0;i<NUMENTITIES;i++){
		vector3 accel_sum={0,0,0};
		for (j=0;j<NUMENTITIES;j++){
			for (k=0;k<3;k++)
				accel_sum[k]+=accels[i][j][k];
		}
		//compute the new velocity based on the acceleration and time interval
		//compute the new position based on the velocity and time interval
		for (k=0;k<3;k++){
			d_hVel[i][k]+=accel_sum[k]*INTERVAL;
			d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
		}
	}

	//I could do another sync threads and then have each thread compute one velocity and one position

	free(accels);
	free(values);
}
