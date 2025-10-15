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

int numNa, numTa, numO, numCl;

int readTraj_xyz(void);

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
		printf("USAGE : ./scattering.x ***.xyz vanHove_diffusion.out vanHove_Al.out vanHove_rotation.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	fp_vanHove_diffusion = fopen(argv[2], "w");
	fp_vanHove_Al = fopen(argv[3], "w");
	fp_vanHove_rotation = fopen(argv[4], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj_xyz();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = (box[0][0][1] - box[0][0][0]);
	box_y = (box[0][1][1] - box[0][1][0]);
	box_z = (box[0][2][1] - box[0][2][0]);

	numNa = 0; numTa = 0; numO = 0; numCl = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numNa++; break;
			case 2: numTa++; break;
			case 3: numO++;  break;
            case 4: numCl++; break;
		}
	}

	printf("\tbox_x = %lf\n", box_x);
	printf("\tbox_y = %lf\n", box_y);
	printf("\tbox_z = %lf\n", box_z);
	printf("\tnumNa = %d\n", numNa);
	printf("\tnumTa = %d\n", numTa);
	printf("\tnumO  = %d\n", numO);
    printf("\tnumCl = %d\n", numCl);
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
		int *temp = (int *)malloc(sizeof(int) * numNa);
		double **tempVector = (double **)malloc(sizeof(double *) * numNa);
		int alIdx = 0;
		for (int i = 0; i < numNa; i++){ temp[i] = -1; }
		for (int i = 0; alIdx < numNa; i++){
			double *vector = (double *)malloc(sizeof(double) * 3);
			if (atom[t][i][1] != 3) continue;
			int nearest = 0; double nearestDist = 999.0; int clIdx = 0;
			for (int j = 0; clIdx < numNa; j++){
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
			double normal = (double)(numTraj - t) / DELTA_T * numNa;
			for (int i = 0; alIdx < numNa; i++){
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
			double normal = (double)(numTraj - t) / DELTA_T * numNa * binsize;
			for (int i = 0; liIdx < numNa; i++){
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
			double normal = (double)(numTraj - t) / DELTA_T * numNa * binsize;
			for (int i = 0; alIdx < numNa; i++){
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

int readTraj_xyz(void){
    char line[LINESIZE];
    char symbol[16];
    double x, y, z;
    int current_natoms;

    while(fscanf(fp_in, "%d", &current_natoms) == 1) {
        if (numTraj >= MAXTIMESTEP) {
            printf("Warning: Maximum number of timesteps (%d) reached.\n", MAXTIMESTEP);
            break;
        }

        numAtoms[numTraj] = current_natoms;
        fgets(line, LINESIZE, fp_in);
        fgets(line, LINESIZE, fp_in);

        char* lattice_ptr = strstr(line, "Lattice=\"");
        if (lattice_ptr != NULL) {
            lattice_ptr += 9;
            double h[9];
            sscanf(lattice_ptr, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &h[0], &h[1], &h[2], &h[3], &h[4], &h[5], &h[6], &h[7], &h[8]);
            
            box[numTraj][0][0] = 0.0; box[numTraj][0][1] = h[0];
            box[numTraj][1][0] = 0.0; box[numTraj][1][1] = h[4];
            box[numTraj][2][0] = 0.0; box[numTraj][2][1] = h[8];
        } else {
            printf("FATAL: Lattice information not found in frame %d. Aborting.\n", numTraj);
            exit(1);
        }

        int** atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
        double** coordPerTraj = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

        for (int i = 0; i < numAtoms[numTraj]; i++) {
            atomPerTraj[i] = (int*)malloc(sizeof(int) * 2);
            coordPerTraj[i] = (double*)malloc(sizeof(double) * 3);

			char line_buffer[LINESIZE];
            if (fgets(line_buffer, LINESIZE, fp_in) == NULL) {
                printf("ERROR: Unexpected end of file or read error in frame %d\n", numTraj);
                exit(1);
            }

            if (sscanf(line_buffer, "%s %lf %lf %lf", symbol, &x, &y, &z) < 4) {
                printf("ERROR: Failed to parse atom line %d in frame %d\n", i + 1, numTraj);
                exit(1);
            }
        
            atomPerTraj[i][0] = i + 1;

           
            if (strcmp(symbol, "Na") == 0)      atomPerTraj[i][1] = 1;
            else if (strcmp(symbol, "Ta") == 0) atomPerTraj[i][1] = 2;
            else if (strcmp(symbol, "O") == 0)  atomPerTraj[i][1] = 3;
            else if (strcmp(symbol, "Cl") == 0) atomPerTraj[i][1] = 4;
            else                                atomPerTraj[i][1] = 0; 
            
            coordPerTraj[i][0] = x;
            coordPerTraj[i][1] = y;
            coordPerTraj[i][2] = z;
        }

        atom[numTraj] = atomPerTraj;
        coord[numTraj] = coordPerTraj;
        timestep[numTraj] = numTraj;

        numTraj++;
    }

    fclose(fp_in);
    return 0;
}