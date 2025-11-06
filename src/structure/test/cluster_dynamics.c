#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 1024
#define MAXTIMESTEP 100000
#define TA_O_CUTOFF  2.9
#define TA_CL_CUTOFF 3.5

#define CLUSTER_TYPE_POLYMERIC 0
#define CLUSTER_TYPE_ISOLATED  1

#define MAX_DELTA_T 10000     
#define LIFETIME_HIST_BINS 10000 
#define DELTA_T_STEP 1           

int numTraj = 0;
int numAtoms[MAXTIMESTEP];
int ***atom;
double ***coord;
double ***lattice_vectors;
int numNa, numTa, numO, numCl;

int **cluster; 
int **na_classification; 
int **ta_classification; 

int read_xyz_traj(const char *filename);
void perform_clustering();
void classify_na_ions(); 
int check_for_updates(int n, const int *arr1, const int *arr2);
void invert_matrix_3x3(const double A[3][3], double A_inv[3][3]);

void calculate_intermittent_tcf(const char *poly_tcf_file, const char *iso_tcf_file);
void calculate_cluster_size_autocorrelation(const char *size_tcf_file);
void calculate_exchange_dynamics(const char *iso_to_poly_file, const char *poly_to_iso_file);


int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("USAGE: %s <input.xyz> <tcf_poly.dat> <tcf_iso.dat> <tcf_size_corr.dat> <life_iso_poly.dat> <life_poly_iso.dat>\n", argv[0]);
        exit(1);
    }
    const char *input_file = argv[1];
    const char *tcf_poly_output_file = argv[2];
    const char *tcf_iso_output_file = argv[3];
    const char *tcf_size_corr_output_file = argv[4];
    const char *life_iso_poly_output_file = argv[5];
    const char *life_poly_iso_output_file = argv[6];


    atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
    coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);
    lattice_vectors = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

    if (!atom || !coord || !lattice_vectors) {
        fprintf(stderr, "Error: Failed to allocate memory for trajectory pointers.\n");
        exit(1);
    }

    if (read_xyz_traj(input_file) != 0) {
        fprintf(stderr, "Error: Failed to read the trajectory file.\n");
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

    na_classification = (int **)malloc(sizeof(int *) * numTraj);
    ta_classification = (int **)malloc(sizeof(int *) * numTraj);
    if (!na_classification || !ta_classification) {
        fprintf(stderr, "Error: Failed to allocate memory for classification arrays.\n");
        exit(1);
    }

    perform_clustering();
    classify_na_ions(); 

    printf("\n\tStarting Dynamics Analysis...\n");
    
    calculate_intermittent_tcf(tcf_poly_output_file, tcf_iso_output_file);
    
    calculate_cluster_size_autocorrelation(tcf_size_corr_output_file);
    
    calculate_exchange_dynamics(life_iso_poly_output_file, life_poly_iso_output_file);

    printf("\tFreeing memory...\n");
    
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
        free(cluster[t]);
        free(na_classification[t]); 
        free(ta_classification[t]);
    }
    free(atom);
    free(coord);
    free(lattice_vectors);
    free(cluster);
    free(na_classification);
    free(ta_classification);

    printf("\tAll Dynamics Tasks are Done! >:D\n");
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
            for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) lattice_vectors[numTraj][i][j] = (i==j) ? 1.0e6 : 0.0; 
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

            atom[numTraj][i][0] = i;
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
            } else { 
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
                if (atom[t][i][1] == 1) continue; 

                for (int j = i + 1; j < nAtom; j++) {
                    if (atom[t][j][1] == 1) continue; 

                    int type_i = atom[t][i][1];
                    int type_j = atom[t][j][1];
                    double current_cutoff_sq;

                    if (type_i > type_j) { int temp = type_i; type_i = type_j; type_j = temp; }

                    if (type_i == 2 && type_j == 3) { 
                        current_cutoff_sq = ta_o_cutoff_sq;
                    } else if (type_i == 2 && type_j == 4) { 
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


void classify_na_ions() {
    printf("\tNow Classifying Na+ ions and Ta atoms...\n");

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Classifying frame t = %d...\n", t);

        int nAtom = numAtoms[t];
        na_classification[t] = (int *)malloc(sizeof(int) * numNa);
        ta_classification[t] = (int *)malloc(sizeof(int) * numTa); 
        
        if (!na_classification[t] || !ta_classification[t]) {
            fprintf(stderr, "Error: Failed to allocate na_classification[t] or ta_classification[t]\n");
            exit(1);
        }

        int max_id = 0;
        for (int i = 0; i < nAtom; i++) if (cluster[t][i] > max_id) max_id = cluster[t][i];

        int *cluster_size = (int *)calloc(max_id + 1, sizeof(int));
        int *cluster_Ta_count = (int *)calloc(max_id + 1, sizeof(int));
        int *cluster_O_count = (int *)calloc(max_id + 1, sizeof(int));
        int *cluster_Cl_count = (int *)calloc(max_id + 1, sizeof(int));

        for (int i = 0; i < nAtom; i++) {
            int cid = cluster[t][i];
            if (cid == -1) continue; 
            cluster_size[cid]++;
            switch (atom[t][i][1]) {
                case 2: cluster_Ta_count[cid]++; break;
                case 3: cluster_O_count[cid]++; break;
                case 4: cluster_Cl_count[cid]++; break;
            }
        }

        int *cluster_type = (int *)malloc(sizeof(int) * (max_id + 1));
        for (int c = 0; c <= max_id; c++) {
            if (cluster_size[c] == 0) continue;

            int is_TaCl6 = (cluster_Ta_count[c] == 1 &&
                            cluster_O_count[c] == 0 &&
                            cluster_Cl_count[c] == 6 &&
                            cluster_size[c] == 7);
            
            int is_TaOCl5 = (cluster_Ta_count[c] == 1 &&
                             cluster_O_count[c] == 1 &&
                             cluster_Cl_count[c] == 5 &&
                             cluster_size[c] == 7);
            
            if (is_TaCl6 || is_TaOCl5) {
                cluster_type[c] = CLUSTER_TYPE_ISOLATED;
            } else {
                cluster_type[c] = CLUSTER_TYPE_POLYMERIC;
            }
        }

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        int na_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 1) { 
                double min_dist_sq = 1.0e18;
                int closest_atom_idx = -1;
                for (int j = 0; j < nAtom; j++) {
                    if (atom[t][j][1] == 1) continue; 

                    double dr[3], df[3], dr_pbc[3];
                    dr[0] = coord[t][j][0] - coord[t][i][0];
                    dr[1] = coord[t][j][1] - coord[t][i][1];
                    dr[2] = coord[t][j][2] - coord[t][i][2];
                    df[0] = inv_lattice[0][0] * dr[0] + inv_lattice[0][1] * dr[1] + inv_lattice[0][2] * dr[2];
                    df[1] = inv_lattice[1][0] * dr[0] + inv_lattice[1][1] * dr[1] + inv_lattice[1][2] * dr[2];
                    df[2] = inv_lattice[2][0] * dr[0] + inv_lattice[2][1] * dr[1] + inv_lattice[2][2] * dr[2];
                    df[0] -= round(df[0]);
                    df[1] -= round(df[1]);
                    df[2] -= round(df[2]);
                    dr_pbc[0] = current_lattice[0][0] * df[0] + current_lattice[1][0] * df[1] + current_lattice[2][0] * df[2];
                    dr_pbc[1] = current_lattice[0][1] * df[0] + current_lattice[1][1] * df[1] + current_lattice[2][1] * df[2];
                    dr_pbc[2] = current_lattice[0][2] * df[0] + current_lattice[1][2] * df[1] + current_lattice[2][2] * df[2];
                    double distSq = dr_pbc[0] * dr_pbc[0] + dr_pbc[1] * dr_pbc[1] + dr_pbc[2] * dr_pbc[2];
                    if (distSq < min_dist_sq) {
                        min_dist_sq = distSq;
                        closest_atom_idx = j;
                    }
                }
                if (closest_atom_idx != -1) {
                    int closest_cluster_id = cluster[t][closest_atom_idx];
                    if(closest_cluster_id >= 0 && closest_cluster_id <= max_id) {
                        na_classification[t][na_idx_counter] = cluster_type[closest_cluster_id];
                    } else {
                        na_classification[t][na_idx_counter] = CLUSTER_TYPE_POLYMERIC; 
                    }
                } else {
                    na_classification[t][na_idx_counter] = CLUSTER_TYPE_POLYMERIC; 
                }
                na_idx_counter++;
            }
        }
        
        int ta_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 2) { 
                int cid = cluster[t][i];
                if (cid >= 0 && cid <= max_id) {
                    ta_classification[t][ta_idx_counter] = cluster_type[cid];
                } else {
                    ta_classification[t][ta_idx_counter] = CLUSTER_TYPE_POLYMERIC; 
                }
                ta_idx_counter++;
            }
        }
        free(cluster_size);
        free(cluster_Ta_count);
        free(cluster_O_count);
        free(cluster_Cl_count);
        free(cluster_type);
    }
}

int check_for_updates(int n, const int *arr1, const int *arr2) {
    for (int i = 0; i < n; i++) {
        if (arr1[i] != arr2[i]) return 1;
    }
    return 0;
}

void invert_matrix_3x3(const double A[3][3], double A_inv[3][3]) {
    double det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
                 A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                 A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    if (fabs(det) < 1.0e-12) {
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) A_inv[i][j] = (i==j) ? 1.0 : 0.0;
        return;
    }
    double inv_det = 1.0 / det;
    A_inv[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * inv_det;
    A_inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * inv_det;
    A_inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_det;
    A_inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * inv_det;
    A_inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_det;
    A_inv[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * inv_det;
    A_inv[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * inv_det;
    A_inv[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * inv_det;
    A_inv[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * inv_det;
}


void calculate_intermittent_tcf(const char *poly_tcf_file, const char *iso_tcf_file) {
    printf("\tCalculating Intermittent TCF (Poly & Iso)...\n");

    int max_dt = (numTraj / 3 < MAX_DELTA_T) ? (numTraj / 3) : MAX_DELTA_T;
    if (max_dt == 0) max_dt = 1;

    long *sum_h0_poly = (long *)calloc(max_dt, sizeof(long));
    long *sum_h0ht_poly = (long *)calloc(max_dt, sizeof(long));
    long *sum_h0_iso = (long *)calloc(max_dt, sizeof(long));
    long *sum_h0ht_iso = (long *)calloc(max_dt, sizeof(long));

    if (!sum_h0_poly || !sum_h0ht_poly || !sum_h0_iso || !sum_h0ht_iso) {
        fprintf(stderr, "Error: Failed to allocate TCF arrays.\n");
        return;
    }

    int *ta_atom_indices = (int *)malloc(sizeof(int) * numTa);
    int ta_idx_counter = 0;
    for (int i = 0; i < numAtoms[0]; i++) {
        if (atom[0][i][1] == 2) {
            ta_atom_indices[ta_idx_counter++] = i;
        }
    }

    for (int start = 0; start < numTraj; start += DELTA_T_STEP) {
        if (start % 100 == 0) printf("\t TCF t0 = %d...\n", start);

        for (int i = 0; i < numTa; i++) {
            int atom_idx = ta_atom_indices[i];
            
            if (ta_classification[start][i] == CLUSTER_TYPE_ISOLATED) {
                for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
                    if (start + dt >= numTraj) break;
                    
                    sum_h0_iso[dt]++; 
                    
                    if (ta_classification[start + dt][i] == CLUSTER_TYPE_ISOLATED) {
                        sum_h0ht_iso[dt]++; 
                    }
                }
            }
        }
        
        for (int i = 0; i < numTa; i++) {
            int atom_idx_i = ta_atom_indices[i];
            
            if (ta_classification[start][i] == CLUSTER_TYPE_POLYMERIC) {
                int cid_i = cluster[start][atom_idx_i];
                if (cid_i == -1) continue;

                for (int j = i + 1; j < numTa; j++) {
                    int atom_idx_j = ta_atom_indices[j];
                    
                    if (ta_classification[start][j] == CLUSTER_TYPE_POLYMERIC && 
                        cluster[start][atom_idx_j] == cid_i) 
                    {
                        for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
                            if (start + dt >= numTraj) break;
                            
                            sum_h0_poly[dt]++; 
                            
                            int cid_i_t = cluster[start + dt][atom_idx_i];
                            int cid_j_t = cluster[start + dt][atom_idx_j];
                            if (cid_i_t != -1 && cid_i_t == cid_j_t) {
                                sum_h0ht_poly[dt]++; 
                            }
                        }
                    }
                }
            }
        }
    } 

    FILE *fp_poly = fopen(poly_tcf_file, "w");
    fprintf(fp_poly, "# Intermittent TCF for Ta-Ta pairs in Polyanions\n");
    fprintf(fp_poly, "# dt(step)\tTCF C_poly(t) = <h_ij(0)h_ij(t)>/<h_ij(0)>\t<h(0)h(t)>\t<h(0)>\n");
    for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
        double tcf_poly = (sum_h0_poly[dt] > 0) ? (double)sum_h0ht_poly[dt] / (double)sum_h0_poly[dt] : 0.0;
        fprintf(fp_poly, "%d\t%10.6f\t%ld\t%ld\n", dt, tcf_poly, sum_h0ht_poly[dt], sum_h0_poly[dt]);
    }
    fclose(fp_poly);

    FILE *fp_iso = fopen(iso_tcf_file, "w");
    fprintf(fp_iso, "# Intermittent TCF for Isolated Ta atoms\n");
    fprintf(fp_iso, "# dt(step)\tTCF C_iso(t) = <h_i(0)h_i(t)>/<h_i(0)>\t<h(0)h(t)>\t<h(0)>\n");
    for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
        double tcf_iso = (sum_h0_iso[dt] > 0) ? (double)sum_h0ht_iso[dt] / (double)sum_h0_iso[dt] : 0.0;
        fprintf(fp_iso, "%d\t%10.6f\t%ld\t%ld\n", dt, tcf_iso, sum_h0ht_iso[dt], sum_h0_iso[dt]);
    }
    fclose(fp_iso);

    free(ta_atom_indices);
    free(sum_h0_poly);
    free(sum_h0ht_poly);
    free(sum_h0_iso);
    free(sum_h0ht_iso);
}


void calculate_cluster_size_autocorrelation(const char *size_tcf_file) {
    printf("\tCalculating Cluster Size Autocorrelation (Ta atoms)...\n");

    int max_dt = (numTraj / 3 < MAX_DELTA_T) ? (numTraj / 3) : MAX_DELTA_T;
    if (max_dt == 0) max_dt = 1;

    int **ta_cluster_size = (int **)malloc(sizeof(int *) * numTraj);
    int *ta_atom_indices = (int *)malloc(sizeof(int) * numTa);
    double avg_cluster_size = 0.0;
    long total_ta_count_valid = 0;

    int ta_idx_counter = 0;
    for (int i = 0; i < numAtoms[0]; i++) {
        if (atom[0][i][1] == 2) ta_atom_indices[ta_idx_counter++] = i;
    }

    for (int t = 0; t < numTraj; t++) {
        ta_cluster_size[t] = (int *)malloc(sizeof(int) * numTa);
        int nAtom = numAtoms[t];
        
        int max_id = 0;
        for(int i=0; i<nAtom; ++i) if(cluster[t][i] > max_id) max_id = cluster[t][i];
        int* all_cluster_sizes = (int*)calloc(max_id + 1, sizeof(int));
        for (int i = 0; i < nAtom; i++) {
            if (cluster[t][i] != -1) all_cluster_sizes[cluster[t][i]]++;
        }

        for (int i = 0; i < numTa; i++) {
            int atom_idx = ta_atom_indices[i];
            int cid = cluster[t][atom_idx];
            int size = (cid != -1) ? all_cluster_sizes[cid] : 1; 
            
            ta_cluster_size[t][i] = size;
            avg_cluster_size += (double)size;
            total_ta_count_valid++;
        }
        free(all_cluster_sizes);
    }
    
    if (total_ta_count_valid > 0) {
        avg_cluster_size /= (double)total_ta_count_valid;
    }
    printf("\tAverage Ta cluster size <N> = %f\n", avg_cluster_size);

    double *numerator = (double *)calloc(max_dt, sizeof(double));
    double denominator = 0.0;
    long *count = (long *)calloc(max_dt, sizeof(long));
    long count_denom = 0;

    for (int start = 0; start < numTraj; start += DELTA_T_STEP) {
        if (start % 100 == 0) printf("\t Size TCF t0 = %d...\n", start);
        for (int i = 0; i < numTa; i++) {
            
            double delta_n_t0 = (double)ta_cluster_size[start][i] - avg_cluster_size;
            
            denominator += delta_n_t0 * delta_n_t0;
            count_denom++;

            for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
                if (start + dt >= numTraj) break;
                
                double delta_n_t = (double)ta_cluster_size[start + dt][i] - avg_cluster_size;
                
                numerator[dt] += delta_n_t0 * delta_n_t;
                count[dt]++;
            }
        }
    }

    FILE *fp_size = fopen(size_tcf_file, "w");
    fprintf(fp_size, "# Ta Cluster Size Autocorrelation C_N(t) = <dN(0)dN(t)>/<(dN(0))^2>\n");
    fprintf(fp_size, "# <N> = %f\n", avg_cluster_size);
    fprintf(fp_size, "# dt(step)\tTCF C_N(t)\t<dN(0)dN(t)>\t<(dN(0))^2>\n");

    double denominator_avg = (count_denom > 0) ? (denominator / (double)count_denom) : 0.0;

    for (int dt = 0; dt < max_dt; dt += DELTA_T_STEP) {
        double numerator_avg = (count[dt] > 0) ? (numerator[dt] / (double)count[dt]) : 0.0;
        double tcf_size = (denominator_avg > 1e-9) ? (numerator_avg / denominator_avg) : 0.0;
        
        fprintf(fp_size, "%d\t%10.6f\t%10.6f\t%10.6f\n", dt, tcf_size, numerator_avg, denominator_avg);
    }
    fclose(fp_size);

    free(numerator);
    free(count);
    free(ta_atom_indices);
    for (int t = 0; t < numTraj; t++) free(ta_cluster_size[t]);
    free(ta_cluster_size);
}


void calculate_exchange_dynamics(const char *iso_to_poly_file, const char *poly_to_iso_file) {
    printf("\tCalculating Cluster Exchange Dynamics (Lifetimes)...\n");

    long *hist_iso_to_poly = (long *)calloc(LIFETIME_HIST_BINS, sizeof(long));
    long *hist_poly_to_iso = (long *)calloc(LIFETIME_HIST_BINS, sizeof(long));
    long total_iso_to_poly_events = 0;
    long total_poly_to_iso_events = 0;

    if (!hist_iso_to_poly || !hist_poly_to_iso) {
        fprintf(stderr, "Error: Failed to allocate lifetime histograms.\n");
        return;
    }

    for (int i = 0; i < numTa; i++) {
        
        int current_state = ta_classification[0][i];
        int last_transition_time = 0;

        for (int t = 1; t < numTraj; t++) {
            int new_state = ta_classification[t][i];

            if (new_state != current_state) {
                int lifetime = t - last_transition_time;
                
                if (lifetime < LIFETIME_HIST_BINS) {
                    if (current_state == CLUSTER_TYPE_ISOLATED) {
                        hist_iso_to_poly[lifetime]++;
                        total_iso_to_poly_events++;
                    } else {
                        hist_poly_to_iso[lifetime]++;
                        total_poly_to_iso_events++;
                    }
                }
                
                current_state = new_state;
                last_transition_time = t;
            }
        } 
    } 
    
    printf("\tTotal Iso->Poly events: %ld\n", total_iso_to_poly_events);
    printf("\tTotal Poly->Iso events: %ld\n", total_poly_to_iso_events);

    FILE *fp_iso_poly = fopen(iso_to_poly_file, "w");
    fprintf(fp_iso_poly, "# Lifetime distribution for ISOLATED state (time until transition to Poly)\n");
    fprintf(fp_iso_poly, "# Lifetime(step)\tCount\n");
    for (int i = 0; i < LIFETIME_HIST_BINS; i++) {
        if (hist_iso_to_poly[i] > 0) {
            fprintf(fp_iso_poly, "%d\t%ld\n", i, hist_iso_to_poly[i]);
        }
    }
    fclose(fp_iso_poly);

    FILE *fp_poly_iso = fopen(poly_to_iso_file, "w");
    fprintf(fp_poly_iso, "# Lifetime distribution for POLYMERIC state (time until transition to Iso)\n");
    fprintf(fp_poly_iso, "# Lifetime(step)\tCount\n");
    for (int i = 0; i < LIFETIME_HIST_BINS; i++) {
        if (hist_poly_to_iso[i] > 0) {
            fprintf(fp_poly_iso, "%d\t%ld\n", i, hist_poly_to_iso[i]);
        }
    }
    fclose(fp_poly_iso);

    free(hist_iso_to_poly);
    free(hist_poly_to_iso);
}