#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 1024
#define MAXTIMESTEP 100000
#define TA_O_CUTOFF  2.9
#define TA_CL_CUTOFF 3.5

int numTraj = 0;
int numAtoms[MAXTIMESTEP];
int ***atom;
double ***coord;
double ***lattice_vectors;
int numNa, numTa, numO, numCl;
int **cluster; // [t][atom_idx] -> cluster_id

int **Ta_Ta_TCF_dat; // [t][data]

int read_xyz_traj(const char *filename);
void perform_clustering();

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("USAGE: %s <input.xyz> <Ta-Ta_TCF.dat>\n", argv[0]);
        exit(1);
    }
    const char *input_file = argv[1];
    const char *tata_tcf_output_file = argv[2];

    atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
    coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);
    lattice_vectors = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

    if (!atom || !coord || !lattice_vectors) {
        fprintf(stderr, "Error: Failed to allocate memory for trajectory pointers.\n");
        exit(1);
    }

    if (read_xyz_traj(input_file) != 0) {
        fprintf(stderr, "Error reading trajectory file.\n");
        exit(1);
    }

    printf("\n\t\t< FILE SUMMARY >\n\n");
    printf("\tNumber of Timesteps Read : %d\n", numTraj);
    printf("\tNumber of Atoms per Frame: %d\n\n", numAtoms[0]);
    
    numNa = 0; numTa = 0; numO = 0; numCl = 0;
    for (int i = 0; i < numAtoms[0]; i++) {
        switch (atom[0][i][1]) {
            case 1: numNa++; break;
            case 2: numTa++; break;
            case 3: numO++; break;
            case 4: numCl++; break;
        }
    }
    printf("\tAtom Counts:\n");
    printf("\t  Na: %d\n", numNa);
    printf("\t  Ta: %d\n", numTa);
    printf("\t  O : %d\n", numO);
    printf("\t  Cl: %d\n", numCl);
    printf("\n");

    perform_clustering();

    // Free allocated memory
    for (int t = 0; t < numTraj; t++) {
        for (int i = 0; i < numAtoms[t]; i++) {
            free(atom[t][i]);
            free(coord[t][i]);
        }
        free(atom[t]);
        free(coord[t]);
        for (int i = 0; i < 3; i++) {
            free(lattice_vectors[t][i]);
        }
        free(lattice_vectors[t]);
    }
    free(atom);
    free(coord);
    free(lattice_vectors);

    printf("\tAll Tasks are Done! >:D\n");
    return 0;
}

int read_xyz_traj(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening input file");
        return -1;
    }
    printf("\tReading trajectory from %s...\n", filename);

    char line[LINESIZE];
    char species_name[10];

    while (fgets(line, LINESIZE, fp) != NULL) {
        if (numTraj >= MAXTIMESTEP) {
            fprintf(stderr, "Warning: Reached MAXTIMESTEP limit (%d). Stopping read.\n", MAXTIMESTEP);
            break;
        }

        if (sscanf(line, "%d", &numAtoms[numTraj]) != 1) continue;

        if (fgets(line, LINESIZE, fp) == NULL) break;
        
        lattice_vectors[numTraj] = (double**)malloc(sizeof(double*) * 3);
        for(int i=0; i<3; ++i) lattice_vectors[numTraj][i] = (double*)malloc(sizeof(double) * 3);

        char *lattice_str = strstr(line, "Lattice=\"");
        if (lattice_str != NULL) {
            sscanf(lattice_str + 9, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   &lattice_vectors[numTraj][0][0], &lattice_vectors[numTraj][0][1], &lattice_vectors[numTraj][0][2],
                   &lattice_vectors[numTraj][1][0], &lattice_vectors[numTraj][1][1], &lattice_vectors[numTraj][1][2],
                   &lattice_vectors[numTraj][2][0], &lattice_vectors[numTraj][2][1], &lattice_vectors[numTraj][2][2]);
        } else {
            fprintf(stderr, "Warning: 'Lattice' information not found in frame %d header. Assuming non-periodic.\n", numTraj);
            for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) lattice_vectors[numTraj][i][j] = (i==j) ? 1.0e6 : 0.0; // Large orthogonal box
        }

        atom[numTraj] = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
        coord[numTraj] = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

        for (int i = 0; i < numAtoms[numTraj]; i++) {
            if (fgets(line, LINESIZE, fp) == NULL) {
                fprintf(stderr, "Error: Unexpected end of file in frame %d.\n", numTraj);
                fclose(fp);
                return -1;
            }
            atom[numTraj][i] = (int*)malloc(sizeof(int) * 2);
            coord[numTraj][i] = (double*)malloc(sizeof(double) * 3);
            
            sscanf(line, "%s %lf %lf %lf",
                   species_name, &coord[numTraj][i][0], &coord[numTraj][i][1], &coord[numTraj][i][2]);

            atom[numTraj][i][0] = i; // Atom ID
            if (strcmp(species_name, "Na") == 0) atom[numTraj][i][1] = 1;
            else if (strcmp(species_name, "Ta") == 0) atom[numTraj][i][1] = 2;
            else if (strcmp(species_name, "O") == 0) atom[numTraj][i][1] = 3;
            else if (strcmp(species_name, "Cl") == 0) atom[numTraj][i][1] = 4;
            else atom[numTraj][i][1] = 0; 
        }
        numTraj++;
    }

    fclose(fp);
    return 0;
}

void perform_clustering() {
    printf("\tNow Performing Clustering Analysis...\n");

    const double ta_o_cutoff_sq  = TA_O_CUTOFF * TA_O_CUTOFF;
    const double ta_cl_cutoff_sq = TA_CL_CUTOFF * TA_CL_CUTOFF;

    cluster = (int **)malloc(sizeof(int *) * numTraj);

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Analyzing frame t = %d...\n", t);

        int nAtom = numAtoms[t];
        int *oldCluster = (int *)malloc(sizeof(int) * nAtom);
        int *newCluster = (int *)malloc(sizeof(int) * nAtom);

        int cluster_idx = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 1) { 
                oldCluster[i] = -1;
                newCluster[i] = -1;
            } else { // Ta, O, Cl are part of clusters
                oldCluster[i] = cluster_idx;
                newCluster[i] = cluster_idx;
                cluster_idx++;
            }
        }
        
        double inv_lattice[3][3];
        double current_lattice[3][3];
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        while (1) {
            for (int i = 0; i < nAtom; i++) {
                if (atom[t][i][1] == 1) continue; // Skip Na

                for (int j = i + 1; j < nAtom; j++) {
                    if (atom[t][j][1] == 1) continue; // Skip Na

                    int type_i = atom[t][i][1];
                    int type_j = atom[t][j][1];
                    double current_cutoff_sq;

                    if (type_i > type_j) {
                        int temp = type_i;
                        type_i = type_j;
                        type_j = temp;
                    }

                    if (type_i == 2 && type_j == 3) { // Ta-O
                        current_cutoff_sq = ta_o_cutoff_sq;
                    } else if (type_i == 2 && type_j == 4) { // Ta-Cl
                        current_cutoff_sq = ta_cl_cutoff_sq;
                    } else {
                        continue;
                    }

                    double dr[3], df[3], dr_pbc[3];
                    dr[0] = coord[t][j][0] - coord[t][i][0];
                    dr[1] = coord[t][j][1] - coord[t][i][1];
                    dr[2] = coord[t][j][2] - coord[t][i][2];

                    df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[0][1]*dr[1] + inv_lattice[0][2]*dr[2];
                    df[1] = inv_lattice[1][0]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[1][2]*dr[2];
                    df[2] = inv_lattice[2][0]*dr[0] + inv_lattice[2][1]*dr[1] + inv_lattice[2][2]*dr[2];

                    df[0] = df[0] - round(df[0]);
                    df[1] = df[1] - round(df[1]);
                    df[2] = df[2] - round(df[2]);

                    dr_pbc[0] = current_lattice[0][0]*df[0] + current_lattice[1][0]*df[1] + current_lattice[2][0]*df[2];
                    dr_pbc[1] = current_lattice[0][1]*df[0] + current_lattice[1][1]*df[1] + current_lattice[2][1]*df[2];
                    dr_pbc[2] = current_lattice[0][2]*df[0] + current_lattice[1][2]*df[1] + current_lattice[2][2]*df[2];

                    double distSq = dr_pbc[0]*dr_pbc[0] + dr_pbc[1]*dr_pbc[1] + dr_pbc[2]*dr_pbc[2];
                
                    if (distSq < current_cutoff_sq) { 
                        int id1 = newCluster[i];
                        int id2 = newCluster[j];
                        if (id1 == id2) continue;
                        
                        int min_id = (id1 < id2) ? id1 : id2;
                        int max_id = (id1 > id2) ? id1 : id2;

                        for(int k=0; k<nAtom; ++k){
                            if(newCluster[k] == max_id) newCluster[k] = min_id;
                        }
                    }
                }
            }
            if (check_for_updates(nAtom, oldCluster, newCluster) == 0) break; 

            memcpy(oldCluster, newCluster, sizeof(int) * nAtom);
        }
        
        cluster[t] = newCluster; 
        free(oldCluster);       
    }
}

void Ta_Ta_TCF(){

    Ta_Ta_TCF_dat = (int **) malloc(sizeof(int *) * numTraj);

    for (int t = 0; t < numTraj; t++) {
        Ta_Ta_TCF_dat[t] = (int *) malloc(sizeof(int) * numTa * (numTa-1) / 2);

        for (int i = 0; i < numAtoms[t]; i++) {
            for (int j = i + 1; j < numAtoms[t]; j++) {
                if (atom[t][i][1] != 2) continue; // Skip non-Ta
                if (atom[t][j][1] != 2) continue; // Skip non-Ta

                int ta_1 = atom[t][i][1];
                int ta_2 = atom[t][j][1];

                int num_Ta_pairs = 0;

                if(cluster[t][i] == cluster[t][j]){
                    Ta_Ta_TCF_dat[t][num_Ta_pairs ] = 1; // In same cluster
                } else {
                    Ta_Ta_TCF_dat[t][num_Ta_pairs] = 0; // Not in same cluster
                }
            }

        }
        
    
    }

}