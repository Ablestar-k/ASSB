#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define DELTA_T 1

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numLi, numCl, numAl;

int readTraj(void);

FILE *fp_in;
FILE *fp_vanHove_diffusion;
FILE *fp_vanHove_Al;
FILE *fp_vanHove_rotation;

void vanHove_diffusion();
void vanHove_Al();
void vanHove_rotation();

double get_angle(double *, double *);

int main(int argc, char *argv[]){
	if (argc != 5){
		printf("USAGE : ./scattering.x ***.lammpstrj vanHove_diffusion.out vanHove_Al.out vanHove_rotation.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	fp_vanHove_diffusion = fopen(argv[2], "w");
	fp_vanHove_Al = fopen(argv[3], "w");
	fp_vanHove_rotation = fopen(argv[4], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = (box[0][0][1] - box[0][0][0]);
	box_y = (box[0][1][1] - box[0][1][0]);
	box_z = (box[0][2][1] - box[0][2][0]);

	numLi = 0; numCl = 0; numAl = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numLi += 1; break;
			case 2: numCl += 1; break;
			case 3: numAl += 1; break;
		}
	}
	printf("\tbox_x = %lf\n", box_x * 2.0);
	printf("\tbox_y = %lf\n", box_y * 2.0);
	printf("\tbox_z = %lf\n", box_z * 2.0);
	printf("\tnumLi = %d\n", numLi);
	printf("\tnumCl = %d\n", numCl);
	printf("\tnumAl = %d\n", numAl);
	printf("\n");

	vanHove_diffusion();
	vanHove_Al();
	vanHove_rotation();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void vanHove_rotation(){
	printf("\tNow Initializing Nearest List ...\n");
	double **rotation = (double **)malloc(sizeof(double *) * (numTraj / 3));

	int **nearList = (int **)malloc(sizeof(int *) * numTraj);
	double ***nearVector = (double ***)malloc(sizeof(double **) * numTraj);
	for (int t = 0; t < numTraj; t++){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		int *temp = (int *)malloc(sizeof(int) * numAl);
		double **tempVector = (double **)malloc(sizeof(double *) * numAl);
		int alIdx = 0;
		for (int i = 0; i < numAl; i++){ temp[i] = -1; }
		for (int i = 0; alIdx < numAl; i++){
			double *vector = (double *)malloc(sizeof(double) * 3);
			if (atom[t][i][1] != 3) continue;
			int nearest = 0; double nearestDist = 999.0; int clIdx = 0;
			for (int j = 0; clIdx < numCl; j++){
				if (i == j) continue;
				if (atom[t][j][1] != 2) continue;
				double dx = coord[t][j][0] - coord[t][i][0];
				double dy = coord[t][j][1] - coord[t][i][1];
				double dz = coord[t][j][2] - coord[t][i][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double distance = dx * dx + dy * dy + dz * dz;
				if (distance < nearestDist){
					nearest = j;
					nearestDist = distance; 
					vector[0] = dx; vector[1] = dy; vector[2] = dz;
				}
				clIdx++;
			}
			// printf("nearest[%d][%d] = %d\n", t, alIdx, nearest);
			tempVector[alIdx] = vector;
			temp[alIdx++] = nearest;
		}
		nearList[t] = temp;
		nearVector[t] = tempVector;
	}

	printf("\tNow Getting Rotational Correlation ...\n");
	double rotation_bin = M_PI / NUMBINS;
	for(int t = 0; t < numTraj / 3; t += DELTA_T){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double *temp = (double *)malloc(sizeof(double) * NUMBINS);
		for (int i = 0; i < NUMBINS; i++) temp[i] = 0.0;
		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int alIdx = 0;
			double normal = (double)(numTraj - t) / DELTA_T * numAl;
			for (int i = 0; alIdx < numAl; i++){
				if (atom[t][i][1] != 3) continue;
				
				int j = nearList[start][alIdx];
				double *r_0 = nearVector[start][alIdx];
				// printf("nearest[%d][%d] = %d\n", start, alIdx, nearList[start][alIdx]);
				double dx = coord[start+t][j][0] - coord[start+t][i][0];
				double dy = coord[start+t][j][1] - coord[start+t][i][1];
				double dz = coord[start+t][j][2] - coord[start+t][i][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double r_t[3] = {dx, dy, dz};

				double theta = get_angle(r_0, r_t);
				int idx = round(theta / rotation_bin);
				if (0 < idx && idx < NUMBINS){ temp[idx] += 1.0 / normal; }
				alIdx++;
			}
		}
		rotation[t / DELTA_T] = temp;
	}
	for (int i = 0; i < (numTraj / 3) / DELTA_T; i++){
		for (int j = 0; j < NUMBINS; j++){
			fprintf(fp_vanHove_rotation, " %10.6lf", rotation[i][j]);
		}
		fprintf(fp_vanHove_rotation, "\n");
	}
}

double get_angle(double *r_0, double *r_t){
	double theta;
	double x, y, z;

	x = (r_0[1] * r_t[2] - r_0[2] * r_t[1]);
	y = (r_0[2] * r_t[0] - r_0[0] * r_t[2]);
	z = (r_0[0] * r_t[1] - r_0[1] * r_t[0]);

	double outter = sqrt(x * x + y * y + z * z);
	double inner = r_0[0] * r_t[0] + r_0[1] * r_t[1] + r_0[2] * r_t[2];
	double l1 = sqrt(r_0[0] * r_0[0] + r_0[1] * r_0[1] + r_0[2] * r_0[2]);
	double l2 = sqrt(r_t[0] * r_t[0] + r_t[1] * r_t[1] + r_t[2] * r_t[2]);
	double from_outter = asin(outter / (l1 * l2));
	double from_inner = acos(inner / (l1 * l2));

	return from_inner;
}

void vanHove_diffusion(){
	printf("\tNow Gs(r,t) for Li...\n");
	double **diffusion = (double **)malloc(sizeof(double *) * (numTraj / 3));
	double binsize = 10.0 / NUMBINS;

	for(int t = 0; t < numTraj / 3; t += DELTA_T){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double *temp = (double *)malloc(sizeof(double) * NUMBINS);
		for (int i = 0; i < NUMBINS; i++) temp[i] = 0.0;
		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int liIdx = 0;
			double normal = (double)(numTraj - t) / DELTA_T * numLi * binsize;
			for (int i = 0; liIdx < numLi; i++){
				if (atom[t][i][1] != 1) continue;
				double dx = coord[start + t][i][0] - coord[start][i][0];
				double dy = coord[start + t][i][1] - coord[start][i][1];
				double dz = coord[start + t][i][2] - coord[start][i][2];

				double distance = sqrt(dx * dx + dy * dy + dz * dz);

				int idx = (int)(distance / binsize);
				if (idx < NUMBINS){ temp[idx] += 1.0 / normal; }
				liIdx++;
			}
		}
		diffusion[t / DELTA_T] = temp;
	}

	for (int i = 0; i < (numTraj / 3) / DELTA_T; i++){
		for (int j = 0; j < NUMBINS; j++){
			fprintf(fp_vanHove_diffusion, " %10.6lf", diffusion[i][j]);
		}
		fprintf(fp_vanHove_diffusion, "\n");
	}
}

void vanHove_Al(){
	printf("\tNow Gs(r,t) for Al...\n");
	double **diffusion = (double **)malloc(sizeof(double *) * (numTraj / 3));
	double binsize = 10.0 / NUMBINS;

	for(int t = 0; t < numTraj / 3; t += DELTA_T){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double *temp = (double *)malloc(sizeof(double) * NUMBINS);
		for (int i = 0; i < NUMBINS; i++) temp[i] = 0.0;
		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int alIdx = 0;
			double normal = (double)(numTraj - t) / DELTA_T * numAl * binsize;
			for (int i = 0; alIdx < numAl; i++){
				if (atom[t][i][1] != 3) continue;
				double dx = coord[start + t][i][0] - coord[start][i][0];
				double dy = coord[start + t][i][1] - coord[start][i][1];
				double dz = coord[start + t][i][2] - coord[start][i][2];

				double distance = sqrt(dx * dx + dy * dy + dz * dz);

				int idx = (int)(distance / binsize);
				if (idx < NUMBINS){ temp[idx] += 1.0 / normal; }
				alIdx++;
			}
		}
		diffusion[t / DELTA_T] = temp;
	}

	for (int i = 0; i < (numTraj / 3) / DELTA_T; i++){
		for (int j = 0; j < NUMBINS; j++){
			fprintf(fp_vanHove_Al, " %10.6lf", diffusion[i][j]);
		}
		fprintf(fp_vanHove_Al, "\n");
	}
}


int readTraj(void){
	char *iostat;
	char line[LINESIZE];

	while(1){
		iostat = fgets(line, LINESIZE, fp_in);
		// printf("s", line);

		if (!iostat) break;
		else if (strcmp(line, "ITEM: TIMESTEP\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			timestep[numTraj] = temp;
			// printf("timestep : %d\n", timestep[numTraj]);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			numAtoms[numTraj] = temp;
			// printf("numAtoms : %d\n", numAtoms[numTraj]);
		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0){
			double temp1, temp2;
			for (int i = 0; i < 3; i++){
				fscanf(fp_in, "%lf %lf", &temp1, &temp2);
				box[numTraj][i][0] = temp1;
				box[numTraj][i][1] = temp2;
			}
			/*
			printf("box :\n%f %f\n%f %f\n%f %f\n",
			        box[numTraj][0][0], box[numTraj][0][1],
			        box[numTraj][1][0], box[numTraj][1][1],
			        box[numTraj][2][0], box[numTraj][2][1]);
			*/
		}
		else if (strcmp(line, "ITEM: ATOMS id type x y z\n") == 0){
			double x, y, z;
			int type, id, ix, iy, iz;

			double boxlength[3] = {box[numTraj][0][1] - box[numTraj][0][0],
					       box[numTraj][1][1] - box[numTraj][1][0],
					       box[numTraj][2][1] - box[numTraj][2][0]};

			int **atomPerTraj;
			atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
			double **coordPerTraj;
			coordPerTraj = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

			for (int i = 0; i < numAtoms[numTraj]; i++){
				int *atomTemp;
				atomTemp = (int*)malloc(sizeof(int) * 2);
				double *coordTemp;
				coordTemp = (double*)malloc(sizeof(double) * 3);

				fscanf(fp_in, "%d %d %lf %lf %lf",
				       &id, &type, &x, &y, &z);

				atomTemp[0] = id;
				atomTemp[1] = type;
				atomPerTraj[i] = atomTemp;

				coordTemp[0] = x;
				coordTemp[1] = y;
				coordTemp[2] = z;
				coordPerTraj[i] = coordTemp;

				/*
				printf("timestep : %d id : %d, type : %d\n",
				       timestep[numTraj], atomPerTraj[i][0], atomPerTraj[i][1]);
				printf("\tx : %f y : %f z : %f\n",
				       coordPerTraj[i][0],
				       coordPerTraj[i][1],
				       coordPerTraj[i][2]);
				*/
			}
			atom[numTraj] = atomPerTraj;
			coord[numTraj] = coordPerTraj;
			numTraj += 1;
		}
	}
	fclose(fp_in);
	return 0;
}