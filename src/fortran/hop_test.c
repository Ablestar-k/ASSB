#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 512
#define MAXTIMESTEP 100000

int numTraj = 0;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numNa, numTa, numO, numCl;

int readTraj_xyz(void);

FILE *fp_in;
FILE *fp_out;

int deltaT;
double ***avgPosA;
double ***avgPosB;
double getDistance(double *, double *);
void initialize();
void hop_function();

int main(int argc, char *argv[]){
	if (argc != 4){
		printf("USAGE : ./hop.x ***.xyz DELTA_T hop.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
    if (fp_in == NULL) {
        printf("ERROR: Cannot open input file: %s\n", argv[1]);
        exit(1);
    }
	deltaT = atoi(argv[2]);
	fp_out = fopen(argv[3], "w");

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

	initialize();
	hop_function();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}


double getDistance(double *coord1, double *coord2){
	double dx = coord1[0] - coord2[0];
	double dy = coord1[1] - coord2[1];
	double dz = coord1[2] - coord2[2];

	dx -= round(dx / box_x) * box_x;
	dy -= round(dy / box_y) * box_y;
	dz -= round(dz / box_z) * box_z;

	return dx * dx + dy * dy + dz * dz;
}


void initialize(){
	printf("\t Now Initializing Averaged Position for Na ions...\n");

	int GAP = deltaT / 2; printf("\tGAP = %d\n", GAP);

	avgPosA = (double ***)malloc(sizeof(double **) * numTraj);
	avgPosB = (double ***)malloc(sizeof(double **) * numTraj);

	for (int t = GAP; t < numTraj - GAP; t++){
		double **tempPosA = (double **)malloc(sizeof(double *) * numNa);
		double **tempPosB = (double **)malloc(sizeof(double *) * numNa);

		for (int idx = 0; idx < numNa; idx++){
			double *pA = (double *)calloc(3, sizeof(double));
			double *pB = (double *)calloc(3, sizeof(double));
            pA[0] = 0.0; pA[1] = 0.0; pA[2] = 0.0;
            pB[0] = 0.0; pB[1] = 0.0; pB[2] = 0.0;

			for (int dt = 0; dt < GAP; dt++){
				pA[0] += coord[t - dt][idx][0]; pA[1] += coord[t - dt][idx][1]; pA[2] += coord[t - dt][idx][2];
				pB[0] += coord[t + dt][idx][0]; pB[1] += coord[t + dt][idx][1]; pB[2] += coord[t + dt][idx][2];
			}

			for (int i = 0; i < 3; i++) { pA[i] /= GAP; pB[i] /= GAP; }
			
			tempPosA[idx] = pA;
			tempPosB[idx] = pB;
		}
		avgPosA[t] = tempPosA;
		avgPosB[t] = tempPosB;
	}
}


void hop_function(){
	printf("\tNow Calculating Hop Function for Na ions...\n");
	int GAP = deltaT / 2; printf("\tGAP = %d\n", GAP);
	fprintf(fp_out, "#\tt\tNa_atom_hop_values...\n");

	double **hop = (double **)malloc(sizeof(double *) * numTraj);
	for (int t = GAP + GAP; t < numTraj - GAP - GAP; t++){
		double *perTime = (double *)malloc(sizeof(double) * numNa);

		box_x = (box[t][0][1] - box[t][0][0]);
		box_y = (box[t][1][1] - box[t][1][0]);
		box_z = (box[t][2][1] - box[t][2][0]);

		for (int idx = 0; idx < numNa; idx++){
			double dA = 0.0; double dB = 0.0;
			for (int dt = 0; dt < GAP; dt++){
				dB += getDistance(coord[t + dt][idx], avgPosA[t + dt][idx]);
				dA += getDistance(coord[t - dt][idx], avgPosB[t - dt][idx]);
			}
			dA /= GAP; dB /= GAP;
			perTime[idx] = sqrt(dA * dB);
		}
		// hop[time][hop value for every atom]
		hop[t] = perTime;

		fprintf(fp_out, "%d\t", t*2);
		for (int i = 0; i < numNa; i++){
			fprintf(fp_out, "%lf\t", hop[t][i]);
			// Dump interval 2ps
		}
		fprintf(fp_out, "\n");
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