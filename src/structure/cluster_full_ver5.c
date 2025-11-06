#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 1024
#define MAXTIMESTEP 100000
#define TA_O_CUTOFF  2.9
#define TA_CL_CUTOFF 3.5

#define NUMBINS 1000    
#define DELTA_T 1       
#define MAX_R 10.0      

#define CLUSTER_TYPE_POLYMERIC 0
#define CLUSTER_TYPE_ISOLATED  1

#define OXYGEN_TYPE_NONBRIDGING 0 // 1-coordinate (to Ta)
#define OXYGEN_TYPE_BRIDGING    1 // 2+ coordinate (to Ta)
#define OXYGEN_TYPE_UNBOUND     2 // 0-coordinate (to Ta)

int **na_anion_classification; // 0(NBO/UO), 1(BO), 3(Cl)

int numTraj = 0;
int numAtoms[MAXTIMESTEP];
int ***atom;
double ***coord;
double ***lattice_vectors;
int numNa, numTa, numO, numCl;
int **cluster; // [t][atom_idx] -> cluster_id

int **na_classification; // [t][na_idx] -> 0 (Poly) or 1 (Iso)

int **ta_classification;      // [t][ta_idx] -> 0 (Poly) or 1 (Iso)
int **oxygen_classification;  // [t][o_idx] -> 0 (NBO), 1 (BO), 2 (UO)
int **na_oxygen_classification; // [t][na_idx] -> 0, 1, or 2 (based on closest O)


int read_xyz_traj(const char *filename);
void perform_clustering();
void classify_na_ions(); 
void classify_oxygens();
void calculate_msd_classified_ta(const char *msd_poly_filename, const char *msd_iso_filename);
void analyze_and_write_clusters(const char *cluster_size_filename, const char *na_coord_filename, const char *cluster_dist_filename);
void calculate_van_hove_classified(const char *poly_filename, const char *iso_filename); 
int check_for_updates(int n, const int *arr1, const int *arr2);
void invert_matrix_3x3(const double A[3][3], double A_inv[3][3]);
void analyze_ta_o_coordination(const char *filename);
void analyze_bridging_linkages(const char *filename);
void analyze_oxygen_classification(const char *filename);
void classify_na_by_closest_anion(const char *ratio_filename);
void calculate_van_hove_na_anion(const char *gsrt_nbo_uo_filename, const char *gsrt_bo_filename, const char *gsrt_cl_filename);

int main(int argc, char *argv[]) {
    if (argc != 16) { 
        printf("USAGE: %s <input.xyz> \n"
               "\t<cluster_size.dat> <na_coord.dat> <cluster_dist.dat> \n"
               "\t<vanHove_poly.dat> <vanHove_iso.dat> \n"
               "\t<msd_ta_poly.dat> <msd_ta_iso.dat> \n"
               "\t<na_anion_ratio.dat> \n"                
               "\t<vhf_na_nbo_uo.dat> <vhf_na_bo.dat> <vhf_na_cl.dat> \n" 
               "\t<ta_o_coord.dat> <ta_bridge.dat> <o_type_dist.dat>\n", 
               argv[0]);
        exit(1);
    }
    const char *input_file = argv[1];
    const char *cluster_size_output_file = argv[2];
    const char *na_coord_output_file = argv[3];
    const char *cluster_dist_output_file = argv[4];
    const char *vanHove_poly_output_file = argv[5];
    const char *vanHove_iso_output_file = argv[6];
    
    const char *msd_ta_poly_output_file = argv[7];
    const char *msd_ta_iso_output_file = argv[8];
    
    const char *na_anion_ratio_output_file = argv[9];
    const char *vanHove_na_nbo_uo_output_file = argv[10];
    const char *vanHove_na_bo_output_file = argv[11];
    const char *vanHove_na_cl_output_file = argv[12];
    
    const char *ta_o_coord_output_file = argv[13];
    const char *ta_bridge_linkages_output_file = argv[14];
    const char *o_type_dist_output_file = argv[15];

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

    na_anion_classification = (int **)malloc(sizeof(int *) * numTraj);
    
    ta_classification = (int **)malloc(sizeof(int *) * numTraj);
    oxygen_classification = (int **)malloc(sizeof(int *) * numTraj);
    
    if (!ta_classification || !oxygen_classification || !na_anion_classification) {
        fprintf(stderr, "Error: Failed to allocate memory for classification arrays.\n");
        exit(1);
    }

    perform_clustering();

    classify_na_ions(); 
    

    classify_oxygens();
    analyze_oxygen_classification(o_type_dist_output_file);
    
    classify_na_by_closest_anion(na_anion_ratio_output_file);
 
    analyze_ta_o_coordination(ta_o_coord_output_file);
    analyze_bridging_linkages(ta_bridge_linkages_output_file);
    analyze_and_write_clusters(cluster_size_output_file, na_coord_output_file, cluster_dist_output_file);
    calculate_van_hove_classified(vanHove_poly_output_file, vanHove_iso_output_file);
    calculate_msd_classified_ta(msd_ta_poly_output_file, msd_ta_iso_output_file);
    
    calculate_van_hove_na_anion(vanHove_na_nbo_uo_output_file, 
                                vanHove_na_bo_output_file, 
                                vanHove_na_cl_output_file);

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
        free(oxygen_classification[t]);
        free(na_anion_classification[t]); 
    }
    
    free(atom);
    free(coord);
    free(lattice_vectors);
    free(cluster);
    free(na_classification);

    free(ta_classification);
    free(oxygen_classification);
    free(na_anion_classification); 

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

                    df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                    df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                    df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];

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

    na_classification = (int **)malloc(sizeof(int *) * numTraj);

    if (!na_classification) {
        fprintf(stderr, "Error: Failed to allocate na_classification\n");
        exit(1);
    }

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
            switch (atom[t][i][1]) { // Atom type
                case 2: cluster_Ta_count[cid]++; break;
                case 3: cluster_O_count[cid]++; break;
                case 4: cluster_Cl_count[cid]++; break;
            }
        }

        // Polymeric vs Isolated
        int *cluster_type = (int *)malloc(sizeof(int) * (max_id + 1));
        for (int c = 0; c <= max_id; c++) {
            if (cluster_size[c] == 0) continue;

            // ISOLATED = cluster size == 6 or 7 && isolated Ta
            int is_iso = (cluster_Ta_count[c] == 1 &&
             (cluster_size[c] == 6 || cluster_size[c] == 7));
           
            if (is_iso) {
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
            if (atom[t][i][1] == 1) { // Na+ ion
                double min_dist_sq = 1.0e18;
                int closest_atom_idx = -1;

                // Find closest Ion in cluster
                for (int j = 0; j < nAtom; j++) {
                    if (atom[t][j][1] == 1) continue; // Na-Na skip

                    double dr[3], df[3], dr_pbc[3];
                    dr[0] = coord[t][j][0] - coord[t][i][0];
                    dr[1] = coord[t][j][1] - coord[t][i][1];
                    dr[2] = coord[t][j][2] - coord[t][i][2];

                    df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                    df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                    df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];

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
            if (atom[t][i][1] == 2) { // Ta atom
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


void analyze_and_write_clusters(const char *cluster_size_filename, const char *na_coord_filename, const char *cluster_dist_filename) {
    printf("\tAnalyzing cluster lists and writing to files...\n");

    FILE *fp_size = fopen(cluster_size_filename, "w");
    if (!fp_size) { perror("Error opening cluster size output file"); return; }
    FILE *fp_na = fopen(na_coord_filename, "w");
    if (!fp_na) { perror("Error opening Na+ coordination output file"); fclose(fp_size); return; }

    FILE *fp_dist = fopen(cluster_dist_filename, "w");
    if (!fp_dist) {
        perror("Error opening cluster distribution output file");
        fclose(fp_size);
        fclose(fp_na);
        return;
    }

    fprintf(fp_size, "# Timestep\tNum_Clusters\tMax_Cluster_Size\n");
    fprintf(fp_na, "# Timestep\tNa_Index\tClosest_Atom_Index\tClosest_Atom_Type\tCoordinated_Cluster_ID\tMin_Distance(A)\tNa_Type(0=Poly,1=Iso)\n");
    fprintf(fp_dist, "# Timestep\tCluster_Size\n"); 

    for (int t = 0; t < numTraj; t++) {
        int nAtom = numAtoms[t];

        int max_id = 0;
        for(int i=0; i<nAtom; ++i) if(cluster[t][i] > max_id) max_id = cluster[t][i];
        
        int* counter = (int*)calloc(max_id + 1, sizeof(int));
        if(!counter) {
            fprintf(stderr, "Error: Failed to allocate counter in analyze step (t=%d).\n", t);
            continue;
        }

        for (int i = 0; i < nAtom; i++) {
            if (cluster[t][i] != -1) {
                counter[cluster[t][i]]++;
            }
        }

        int numCluster = 0, maxClusterSize = 0;
        for (int i = 0; i <= max_id; i++) {
            if (counter[i] > 0) {
                numCluster++;
                if (counter[i] > maxClusterSize) {
                    maxClusterSize = counter[i];
                }
                fprintf(fp_dist, "%d\t%d\n", t, counter[i]);
            }
        }
        fprintf(fp_size, "%d\t%d\t%d\n", t, numCluster, maxClusterSize);
        free(counter);

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        int na_idx_counter = 0; // Na ion index
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 1) { // Found a Na+ ion
                double min_dist_sq = 1.0e18; 
                int closest_atom_idx = -1;

                for (int j = 0; j < nAtom; j++) {
                    if (atom[t][j][1] == 1) continue;
                    double dr[3], df[3], dr_pbc[3];
                    dr[0] = coord[t][j][0] - coord[t][i][0];
                    dr[1] = coord[t][j][1] - coord[t][i][1];
                    dr[2] = coord[t][j][2] - coord[t][i][2];

                    df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                    df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                    df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];
                    
                    df[0] -= round(df[0]);
                    df[1] -= round(df[1]);
                    df[2] -= round(df[2]);

                    dr_pbc[0] = current_lattice[0][0]*df[0] + current_lattice[1][0]*df[1] + current_lattice[2][0]*df[2];
                    dr_pbc[1] = current_lattice[0][1]*df[0] + current_lattice[1][1]*df[1] + current_lattice[2][1]*df[2];
                    dr_pbc[2] = current_lattice[0][2]*df[0] + current_lattice[1][2]*df[1] + current_lattice[2][2]*df[2];
                    
                    double distSq = dr_pbc[0]*dr_pbc[0] + dr_pbc[1]*dr_pbc[1] + dr_pbc[2]*dr_pbc[2];

                    if (distSq < min_dist_sq) {
                        min_dist_sq = distSq;
                        closest_atom_idx = j;
                    }
                }

                if (closest_atom_idx != -1) {
                    int coordinated_cluster_id = cluster[t][closest_atom_idx];
                    int closest_atom_type = atom[t][closest_atom_idx][1];
                    fprintf(fp_na, "%d\t%d\t%d\t%d\t%d\t%f\t%d\n",
                            t, i, closest_atom_idx, closest_atom_type,
                            coordinated_cluster_id, sqrt(min_dist_sq),
                            na_classification[t][na_idx_counter]); 
                }
                na_idx_counter++;
            }
        }
    }

    fclose(fp_size);
    fclose(fp_na);
    fclose(fp_dist); 
}


void calculate_van_hove_classified(const char *poly_filename, const char *iso_filename) {
    printf("\tNow Gs(r,t) for classified Na+ (Poly/Iso)...\n");

    FILE *fp_poly = fopen(poly_filename, "w");
    FILE *fp_iso = fopen(iso_filename, "w");
    if (!fp_poly || !fp_iso) {
        perror("Error opening Van Hove output files");
        return;
    }

    double binsize = MAX_R / NUMBINS;
    fprintf(fp_poly, "# Gs(r,t) for Na+ (Polymeric)\n# Row: dt (fs), Col: r (A) (bin_size=%f A)\n", binsize);
    fprintf(fp_iso, "# Gs(r,t) for Na+ (Isolated)\n# Row: dt (fs), Col: r (A) (bin_size=%f A)\n", binsize);

    int max_delta_t = numTraj / 2; 


    for (int dt = 0; dt < max_delta_t; dt += DELTA_T) {
        if (dt == 0) continue; 
        if (dt % 100 == 0) printf("\t Gs(r,t) dt = %d...\n", dt);

        double *hist_poly = (double *)calloc(NUMBINS, sizeof(double));
        double *hist_iso = (double *)calloc(NUMBINS, sizeof(double));
        
        long count_poly = 0;
        long count_iso = 0;

        for (int start = 0; start < numTraj - dt; start += DELTA_T) {
            
            int na_idx = 0; 
            for (int i = 0; i < numAtoms[start]; i++) {
                if (atom[start][i][1] != 1) continue; 

                // Unwrapped trajectory assumed
                double dx = coord[start + dt][i][0] - coord[start][i][0];
                double dy = coord[start + dt][i][1] - coord[start][i][1];
                double dz = coord[start + dt][i][2] - coord[start][i][2];

                double dist = sqrt(dx * dx + dy * dy + dz * dz);
                int bin_idx = (int)(dist / binsize);

                if (bin_idx < NUMBINS) {
                    if (na_classification[start][na_idx] == CLUSTER_TYPE_POLYMERIC) {
                        hist_poly[bin_idx] += 1.0;
                        count_poly++;
                    } else { // CLUSTER_TYPE_ISOLATED
                        hist_iso[bin_idx] += 1.0;
                        count_iso++;
                    }
                }
                na_idx++; 
            }
        } 

        fprintf(fp_poly, "%d", dt); 
        if (count_poly > 0) {
            for (int j = 0; j < NUMBINS; j++) {
                double r = (j + 0.5) * binsize;
                double dr = binsize;
                // Avoid division by zero for r=0
                if (r < 1.0e-9) {
                     fprintf(fp_poly, " 0.0");
                     continue;
                }
                double shell_volume = 4.0 * M_PI * r * r * dr;
                
                double normal_poly = shell_volume * (double)count_poly; 
                if (normal_poly > 1.0e-9) {
                    fprintf(fp_poly, " %10.6lf", hist_poly[j] / normal_poly);
                } else {
                    fprintf(fp_poly, " 0.0");
                }
            }
        } else {
            for (int j = 0; j < NUMBINS; j++) fprintf(fp_poly, " 0.0");
        }
        fprintf(fp_poly, "\n");


        fprintf(fp_iso, "%d", dt); 
        if (count_iso > 0) {
             for (int j = 0; j < NUMBINS; j++) {
                double r = (j + 0.5) * binsize;
                double dr = binsize;
                // Avoid division by zero for r=0
                if (r < 1.0e-9) {
                     fprintf(fp_iso, " 0.0");
                     continue;
                }
                double shell_volume = 4.0 * M_PI * r * r * dr;
                
                double normal_iso = shell_volume * (double)count_iso;
                if (normal_iso > 1.0e-9) {
                    fprintf(fp_iso, " %10.6lf", hist_iso[j] / normal_iso);
                } else {
                    fprintf(fp_iso, " 0.0");
                }
            }
        } else {
            for (int j = 0; j < NUMBINS; j++) fprintf(fp_iso, " 0.0");
        }
        fprintf(fp_iso, "\n");

        free(hist_poly);
        free(hist_iso);

    } 

    fclose(fp_poly);
    fclose(fp_iso);
}


void classify_oxygens() {
    printf("\tNow Classifying Oxygen atoms (BO vs NBO vs UO)...\n");
    const double ta_o_cutoff_sq  = TA_O_CUTOFF * TA_O_CUTOFF;

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Classifying Oxygens frame t = %d...\n", t);

        int nAtom = numAtoms[t];
        oxygen_classification[t] = (int *)malloc(sizeof(int) * numO);
        if (!oxygen_classification[t]) {
            fprintf(stderr, "Error: Failed to allocate oxygen_classification[t]\n");
            exit(1);
        }

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        int o_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 3) { // Oxygen atom
                int ta_neighbor_count = 0;
                
                for (int j = 0; j < nAtom; j++) {
                    if (atom[t][j][1] == 2) { // Tantalum atom
                        
                        double dr[3], df[3], dr_pbc[3];
                        dr[0] = coord[t][j][0] - coord[t][i][0];
                        dr[1] = coord[t][j][1] - coord[t][i][1];
                        dr[2] = coord[t][j][2] - coord[t][i][2];

                        df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                        df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                        df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];

                        df[0] -= round(df[0]);
                        df[1] -= round(df[1]);
                        df[2] -= round(df[2]);

                        dr_pbc[0] = current_lattice[0][0]*df[0] + current_lattice[1][0]*df[1] + current_lattice[2][0]*df[2];
                        dr_pbc[1] = current_lattice[0][1]*df[0] + current_lattice[1][1]*df[1] + current_lattice[2][1]*df[2];
                        dr_pbc[2] = current_lattice[0][2]*df[0] + current_lattice[1][2]*df[1] + current_lattice[2][2]*df[2];

                        double distSq = dr_pbc[0]*dr_pbc[0] + dr_pbc[1]*dr_pbc[1] + dr_pbc[2]*dr_pbc[2];
                        
                        if (distSq < ta_o_cutoff_sq) {
                            ta_neighbor_count++;
                        }
                    }
                } // End Ta(j) loop

                  if (ta_neighbor_count >= 2) {
                    oxygen_classification[t][o_idx_counter] = OXYGEN_TYPE_BRIDGING;    // 1
                } else if (ta_neighbor_count == 1) {
                    oxygen_classification[t][o_idx_counter] = OXYGEN_TYPE_NONBRIDGING; // 0
                } else if (ta_neighbor_count == 0) {
                    oxygen_classification[t][o_idx_counter] = OXYGEN_TYPE_UNBOUND;     // 2
                }
                o_idx_counter++;
            }
        } // End O(i) loop
    }
}


void classify_na_by_closest_anion(const char *ratio_filename) {
    printf("\tNow Classifying Na+ by closest anion (O vs Cl)...\n");
    
    FILE *fp_ratio = fopen(ratio_filename, "w");
    if (!fp_ratio) {
        perror("Error opening Na-Anion ratio output file");
        return;
    }
    fprintf(fp_ratio, "# Timestep\tRatio_NBO_UO\tRatio_BO\tRatio_Cl\n");

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Classifying Na+ by Anion type frame t = %d...\n", t);
        
        int nAtom = numAtoms[t];
        na_anion_classification[t] = (int *)malloc(sizeof(int) * numNa);
        if (!na_anion_classification[t]) {
            fprintf(stderr, "Error: Failed to allocate na_anion_classification[t]\n");
            exit(1);
        }

        // 1. Create O-index mapper (Global ID -> Local O ID)
        int *global_to_o_idx = (int *)malloc(sizeof(int) * nAtom);
        if (!global_to_o_idx) {
            fprintf(stderr, "Error: Failed to allocate global_to_o_idx map (t=%d)\n", t);
            continue;
        }
        int o_count = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 3) {
                global_to_o_idx[i] = o_count++;
            } else {
                global_to_o_idx[i] = -1;
            }
        }

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        long count_nbo_uo = 0;
        long count_bo = 0;
        long count_cl = 0;

        int na_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 1) { // Na+ atom
                double min_dist_sq = 1.0e18;
                int closest_anion_type = -1; // 0=NBO, 1=BO, 2=UO, 3=Cl

                // 2. Loop over ALL anions (O and Cl)
                for (int j = 0; j < nAtom; j++) {
                    int type_j = atom[t][j][1];
                    if (type_j != 3 && type_j != 4) continue; // Skip if not O or Cl

                    double dr[3], df[3], dr_pbc[3];
                    dr[0] = coord[t][j][0] - coord[t][i][0];
                    dr[1] = coord[t][j][1] - coord[t][i][1];
                    dr[2] = coord[t][j][2] - coord[t][i][2];

                    df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                    df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                    df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];
                    df[0] -= round(df[0]); df[1] -= round(df[1]); df[2] -= round(df[2]);
                    dr_pbc[0] = current_lattice[0][0]*df[0] + current_lattice[1][0]*df[1] + current_lattice[2][0]*df[2];
                    dr_pbc[1] = current_lattice[0][1]*df[0] + current_lattice[1][1]*df[1] + current_lattice[2][1]*df[2];
                    dr_pbc[2] = current_lattice[0][2]*df[0] + current_lattice[1][2]*df[1] + current_lattice[2][2]*df[2];

                    double distSq = dr_pbc[0]*dr_pbc[0] + dr_pbc[1]*dr_pbc[1] + dr_pbc[2]*dr_pbc[2];

                    if (distSq < min_dist_sq) {
                        min_dist_sq = distSq;
                        if (type_j == 3) { // Oxygen
                            int o_idx = global_to_o_idx[j];
                            closest_anion_type = oxygen_classification[t][o_idx]; // 0, 1, or 2
                        } else { // Chlorine
                            closest_anion_type = 3; 
                        }
                    }
                } // End Anion(j) loop

                // 3. Store classification and increment counters
                na_anion_classification[t][na_idx_counter] = closest_anion_type;
                
                if (closest_anion_type == OXYGEN_TYPE_BRIDGING) { // 1
                    count_bo++;
                } else if (closest_anion_type == 3) { // 3
                    count_cl++;
                } else { // 0 (NBO) or 2 (UO)
                    count_nbo_uo++;
                }

                na_idx_counter++;
            }
        } // End Na(i) loop

        // 4. Write ratios for this timestep
        double total_na = (double)numNa;
        if (total_na > 0) {
            fprintf(fp_ratio, "%d\t%f\t%f\t%f\n", t, 
                    (double)count_nbo_uo / total_na,
                    (double)count_bo / total_na,
                    (double)count_cl / total_na);
        } else {
            fprintf(fp_ratio, "%d\t0.0\t0.0\t0.0\n", t);
        }

        free(global_to_o_idx);
    }
    fclose(fp_ratio);
}



void calculate_msd_classified_ta(const char *msd_poly_filename, const char *msd_iso_filename) {
    printf("\tNow Calculating MSD for classified Ta+...\n");

    FILE *fp_poly = fopen(msd_poly_filename, "w");
    FILE *fp_iso = fopen(msd_iso_filename, "w");
    if (!fp_poly || !fp_iso) {
        perror("Error opening Ta MSD output files");
        return;
    }

    fprintf(fp_poly, "# MSD for Ta+ (Polymeric)\n# dt (fs)\tMSD (A^2)\tCount\n");
    fprintf(fp_iso, "# MSD for Ta+ (Isolated)\n# dt (fs)\tMSD (A^2)\tCount\n");

    int max_delta_t = numTraj / 2;

    for (int dt = 0; dt < max_delta_t; dt += DELTA_T) {
        if (dt == 0) {
             fprintf(fp_poly, "0\t0.0\t0\n");
             fprintf(fp_iso, "0\t0.0\t0\n");
             continue;
        }
        if (dt % 100 == 0) printf("\t Ta MSD dt = %d...\n", dt);

        double msd_poly = 0.0;
        double msd_iso = 0.0;
        long count_poly = 0;
        long count_iso = 0;

        for (int start = 0; start < numTraj - dt; start += DELTA_T) {
            int ta_idx = 0;
            for (int i = 0; i < numAtoms[start]; i++) {
                if (atom[start][i][1] != 2) continue; // Not a Ta atom

                // unwrapped trajectory assumed (PBC not applied)
                double dx = coord[start + dt][i][0] - coord[start][i][0];
                double dy = coord[start + dt][i][1] - coord[start][i][1];
                double dz = coord[start + dt][i][2] - coord[start][i][2];
                double distSq = dx*dx + dy*dy + dz*dz;

                if (ta_classification[start][ta_idx] == CLUSTER_TYPE_POLYMERIC) {
                    msd_poly += distSq;
                    count_poly++;
                } else { // CLUSTER_TYPE_ISOLATED
                    msd_iso += distSq;
                    count_iso++;
                }
                ta_idx++;
            }
        }

        if (count_poly > 0) {
            fprintf(fp_poly, "%d\t%10.6lf\t%ld\n", dt, msd_poly / (double)count_poly, count_poly);
        } else {
            fprintf(fp_poly, "%d\t0.0\t0\n", dt);
        }

        if (count_iso > 0) {
            fprintf(fp_iso, "%d\t%10.6lf\t%ld\n", dt, msd_iso / (double)count_iso, count_iso);
        } else {
            fprintf(fp_iso, "%d\t0.0\t0\n", dt);
        }
    }

    fclose(fp_poly);
    fclose(fp_iso);
}

void calculate_van_hove_na_anion(const char *gsrt_nbo_uo_filename, 
                                 const char *gsrt_bo_filename, 
                                 const char *gsrt_cl_filename) {
    
    printf("\tNow Gs(r,t) for classified Na+ (by NBO/BO/Cl)...\n");

    FILE *fp_nbo_uo = fopen(gsrt_nbo_uo_filename, "w");
    FILE *fp_bo = fopen(gsrt_bo_filename, "w");
    FILE *fp_cl = fopen(gsrt_cl_filename, "w");
    
    if (!fp_nbo_uo || !fp_bo || !fp_cl) {
        perror("Error opening Van Hove (Anion-classified) output files");
        return;
    }

    double binsize = MAX_R / NUMBINS;
    fprintf(fp_nbo_uo, "# Gs(r,t) for Na+ (near NBO or UO)\n# Row: dt (fs), Col: r (A) (bin_size=%f A)\n", binsize);
    fprintf(fp_bo, "# Gs(r,t) for Na+ (near BO)\n# Row: dt (fs), Col: r (A) (bin_size=%f A)\n", binsize);
    fprintf(fp_cl, "# Gs(r,t) for Na+ (near Cl)\n# Row: dt (fs), Col: r (A) (bin_size=%f A)\n", binsize);

    int max_delta_t = numTraj / 2;

    for (int dt = 0; dt < max_delta_t; dt += DELTA_T) {
        if (dt == 0) continue; 
        if (dt % 100 == 0) printf("\t Gs(r,t) Anion-type dt = %d...\n", dt);

        double *hist_nbo_uo = (double *)calloc(NUMBINS, sizeof(double));
        double *hist_bo = (double *)calloc(NUMBINS, sizeof(double));
        double *hist_cl = (double *)calloc(NUMBINS, sizeof(double));
        
        long count_nbo_uo = 0;
        long count_bo = 0;
        long count_cl = 0;

        for (int start = 0; start < numTraj - dt; start += DELTA_T) {
            
            int na_idx = 0; 
            for (int i = 0; i < numAtoms[start]; i++) {
                if (atom[start][i][1] != 1) continue; 

                // Unwrapped trajectory assumed
                double dx = coord[start + dt][i][0] - coord[start][i][0];
                double dy = coord[start + dt][i][1] - coord[start][i][1];
                double dz = coord[start + dt][i][2] - coord[start][i][2];

                double dist = sqrt(dx * dx + dy * dy + dz * dz);
                int bin_idx = (int)(dist / binsize);

                if (bin_idx < NUMBINS) {
                    // Use the new classification
                    int classification = na_anion_classification[start][na_idx];
                    
                    if (classification == OXYGEN_TYPE_BRIDGING) { // 1
                        hist_bo[bin_idx] += 1.0;
                        count_bo++;
                    } else if (classification == 3) { // 3 (Cl)
                        hist_cl[bin_idx] += 1.0;
                        count_cl++;
                    } else { // 0 (NBO) or 2 (UO)
                        hist_nbo_uo[bin_idx] += 1.0;
                        count_nbo_uo++;
                    }
                }
                na_idx++; 
            }
        } 

        // --- Write NBO/UO Data ---
        fprintf(fp_nbo_uo, "%d", dt); 
        if (count_nbo_uo > 0) {
            for (int j = 0; j < NUMBINS; j++) {
                double r = (j + 0.5) * binsize;
                double dr = binsize;
                if (r < 1.0e-9) { fprintf(fp_nbo_uo, " 0.0"); continue; }
                double shell_volume = 4.0 * M_PI * r * r * dr;
                double normal = shell_volume * (double)count_nbo_uo; 
                if (normal > 1.0e-9) {
                    fprintf(fp_nbo_uo, " %10.6lf", hist_nbo_uo[j] / normal);
                } else {
                    fprintf(fp_nbo_uo, " 0.0");
                }
            }
        } else {
            for (int j = 0; j < NUMBINS; j++) fprintf(fp_nbo_uo, " 0.0");
        }
        fprintf(fp_nbo_uo, "\n");

        // --- Write BO Data ---
        fprintf(fp_bo, "%d", dt); 
        if (count_bo > 0) {
            for (int j = 0; j < NUMBINS; j++) {
                double r = (j + 0.5) * binsize;
                double dr = binsize;
                if (r < 1.0e-9) { fprintf(fp_bo, " 0.0"); continue; }
                double shell_volume = 4.0 * M_PI * r * r * dr;
                double normal = shell_volume * (double)count_bo;
                if (normal > 1.0e-9) {
                    fprintf(fp_bo, " %10.6lf", hist_bo[j] / normal);
                } else {
                    fprintf(fp_bo, " 0.0");
                }
            }
        } else {
            for (int j = 0; j < NUMBINS; j++) fprintf(fp_bo, " 0.0");
        }
        fprintf(fp_bo, "\n");

        // --- Write Cl Data ---
        fprintf(fp_cl, "%d", dt); 
        if (count_cl > 0) {
            for (int j = 0; j < NUMBINS; j++) {
                double r = (j + 0.5) * binsize;
                double dr = binsize;
                if (r < 1.0e-9) { fprintf(fp_cl, " 0.0"); continue; }
                double shell_volume = 4.0 * M_PI * r * r * dr;
                double normal = shell_volume * (double)count_cl;
                if (normal > 1.0e-9) {
                    fprintf(fp_cl, " %10.6lf", hist_cl[j] / normal);
                } else {
                    fprintf(fp_cl, " 0.0");
                }
            }
        } else {
            for (int j = 0; j < NUMBINS; j++) fprintf(fp_cl, " 0.0");
        }
        fprintf(fp_cl, "\n");

        free(hist_nbo_uo);
        free(hist_bo);
        free(hist_cl);
    } 

    fclose(fp_nbo_uo);
    fclose(fp_bo);
    fclose(fp_cl);
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
        fprintf(stderr, "Warning: Matrix is singular, cannot invert. Using identity matrix.\n");
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

void analyze_ta_o_coordination(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening Ta-O coordination output file");
        return;
    }

    fprintf(fp, "# Timestep\tTa_Atom_Index\tTa_O_Coordination_Number\n");

    const double ta_o_cutoff_sq = TA_O_CUTOFF * TA_O_CUTOFF;

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Analyzing Ta-O CN frame t = %d...\n", t);

        int nAtom = numAtoms[t];

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] != 2) continue; // Skip if not Ta

            int ta_o_cn_count = 0;

            for (int j = 0; j < nAtom; j++) {
                if (atom[t][j][1] != 3) continue; // Skip if not O

                double dr[3], df[3], dr_pbc[3];
                dr[0] = coord[t][j][0] - coord[t][i][0];
                dr[1] = coord[t][j][1] - coord[t][i][1];
                dr[2] = coord[t][j][2] - coord[t][i][2];

                df[0] = inv_lattice[0][0]*dr[0] + inv_lattice[1][0]*dr[1] + inv_lattice[2][0]*dr[2];
                df[1] = inv_lattice[0][1]*dr[0] + inv_lattice[1][1]*dr[1] + inv_lattice[2][1]*dr[2];
                df[2] = inv_lattice[0][2]*dr[0] + inv_lattice[1][2]*dr[1] + inv_lattice[2][2]*dr[2];

                df[0] -= round(df[0]);
                df[1] -= round(df[1]);
                df[2] -= round(df[2]);

                dr_pbc[0] = current_lattice[0][0]*df[0] + current_lattice[1][0]*df[1] + current_lattice[2][0]*df[2];
                dr_pbc[1] = current_lattice[0][1]*df[0] + current_lattice[1][1]*df[1] + current_lattice[2][1]*df[2];
                dr_pbc[2] = current_lattice[0][2]*df[0] + current_lattice[1][2]*df[1] + current_lattice[2][2]*df[2];

                double distSq = dr_pbc[0]*dr_pbc[0] + dr_pbc[1]*dr_pbc[1] + dr_pbc[2]*dr_pbc[2];

                if (distSq < ta_o_cutoff_sq) {
                    ta_o_cn_count++;
                }
            } // end j (O atom) loop

            fprintf(fp, "%d\t%d\t%d\n", t, i, ta_o_cn_count);
        } // end i (Ta atom) loop
    } // end t (timestep) loop

    fclose(fp);
}

void analyze_bridging_linkages(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening Ta bridging linkages output file");
        return;
    }

    fprintf(fp, "# Timestep\tTa_Atom_Index\tTa_O_Ta_Bridges\tTa_Cl_Ta_Bridges\tTa_Classification(0=Poly,1=Iso)\n");

    const double ta_o_cutoff_sq = TA_O_CUTOFF * TA_O_CUTOFF;
    const double ta_cl_cutoff_sq = TA_CL_CUTOFF * TA_CL_CUTOFF;

    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Analyzing Ta-Bridging frame t = %d...\n", t);

        int nAtom = numAtoms[t];

        double inv_lattice[3][3];
        double current_lattice[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) current_lattice[i][j] = lattice_vectors[t][i][j];
        invert_matrix_3x3(current_lattice, inv_lattice);

        int *global_to_o_idx = (int *)malloc(sizeof(int) * nAtom);
        if (!global_to_o_idx) {
            fprintf(stderr, "Error: Failed to allocate global_to_o_idx map (t=%d)\n", t);
            continue;
        }
        int o_count = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 3) {
                global_to_o_idx[i] = o_count++;
            } else {
                global_to_o_idx[i] = -1;
            }
        }

        int ta_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] != 2) continue; // Skip if not Ta

            int ta_type = ta_classification[t][ta_idx_counter];
            int ta_o_ta_bridges = 0;
            int ta_cl_ta_bridges = 0;

            if (ta_type == CLUSTER_TYPE_POLYMERIC) {
                
                // Loop over all other atoms (j)
                for (int j = 0; j < nAtom; j++) {
                    int type_j = atom[t][j][1];

                    // --- 1. Check for Ta-O-Ta bridges ---
                    if (type_j == 3) { // j is Oxygen
                        double dr_ij[3], df_ij[3], dr_pbc_ij[3];
                        dr_ij[0] = coord[t][j][0] - coord[t][i][0];
                        dr_ij[1] = coord[t][j][1] - coord[t][i][1];
                        dr_ij[2] = coord[t][j][2] - coord[t][i][2];
                        df_ij[0] = inv_lattice[0][0]*dr_ij[0] + inv_lattice[1][0]*dr_ij[1] + inv_lattice[2][0]*dr_ij[2];
                        df_ij[1] = inv_lattice[0][1]*dr_ij[0] + inv_lattice[1][1]*dr_ij[1] + inv_lattice[2][1]*dr_ij[2];
                        df_ij[2] = inv_lattice[0][2]*dr_ij[0] + inv_lattice[1][2]*dr_ij[1] + inv_lattice[2][2]*dr_ij[2];
                        df_ij[0] -= round(df_ij[0]); df_ij[1] -= round(df_ij[1]); df_ij[2] -= round(df_ij[2]);
                        dr_pbc_ij[0] = current_lattice[0][0]*df_ij[0] + current_lattice[1][0]*df_ij[1] + current_lattice[2][0]*df_ij[2];
                        dr_pbc_ij[1] = current_lattice[0][1]*df_ij[0] + current_lattice[1][1]*df_ij[1] + current_lattice[2][1]*df_ij[2];
                        dr_pbc_ij[2] = current_lattice[0][2]*df_ij[0] + current_lattice[1][2]*df_ij[1] + current_lattice[2][2]*df_ij[2];
                        double distSq_ij = dr_pbc_ij[0]*dr_pbc_ij[0] + dr_pbc_ij[1]*dr_pbc_ij[1] + dr_pbc_ij[2]*dr_pbc_ij[2];

                        if (distSq_ij < ta_o_cutoff_sq) {
                            int o_idx = global_to_o_idx[j];
                            if (oxygen_classification[t][o_idx] == OXYGEN_TYPE_BRIDGING) {
                                ta_o_ta_bridges++;
                            }
                        }
                    } 
                    // --- 2. Check for Ta-Cl-Ta bridges ---
                    else if (type_j == 4) { // j is Chlorine
                        double dr_ij[3], df_ij[3], dr_pbc_ij[3];
                        dr_ij[0] = coord[t][j][0] - coord[t][i][0];
                        dr_ij[1] = coord[t][j][1] - coord[t][i][1];
                        dr_ij[2] = coord[t][j][2] - coord[t][i][2];
                        df_ij[0] = inv_lattice[0][0]*dr_ij[0] + inv_lattice[1][0]*dr_ij[1] + inv_lattice[2][0]*dr_ij[2];
                        df_ij[1] = inv_lattice[0][1]*dr_ij[0] + inv_lattice[1][1]*dr_ij[1] + inv_lattice[2][1]*dr_ij[2];
                        df_ij[2] = inv_lattice[0][2]*dr_ij[0] + inv_lattice[1][2]*dr_ij[1] + inv_lattice[2][2]*dr_ij[2];

                        df_ij[0] -= round(df_ij[0]); df_ij[1] -= round(df_ij[1]); df_ij[2] -= round(df_ij[2]);

                        dr_pbc_ij[0] = current_lattice[0][0]*df_ij[0] + current_lattice[1][0]*df_ij[1] + current_lattice[2][0]*df_ij[2];
                        dr_pbc_ij[1] = current_lattice[0][1]*df_ij[0] + current_lattice[1][1]*df_ij[1] + current_lattice[2][1]*df_ij[2];
                        dr_pbc_ij[2] = current_lattice[0][2]*df_ij[0] + current_lattice[1][2]*df_ij[1] + current_lattice[2][2]*df_ij[2];
                        double distSq_ij = dr_pbc_ij[0]*dr_pbc_ij[0] + dr_pbc_ij[1]*dr_pbc_ij[1] + dr_pbc_ij[2]*dr_pbc_ij[2];

                        if (distSq_ij < ta_cl_cutoff_sq) {
                            int is_cl_bridged = 0;
                            for (int k = 0; k < nAtom; k++) {
                                if (k == i) continue; 
                                if (atom[t][k][1] != 2) continue; 

                                double dr_jk[3], df_jk[3], dr_pbc_jk[3];
                                dr_jk[0] = coord[t][k][0] - coord[t][j][0];
                                dr_jk[1] = coord[t][k][1] - coord[t][j][1];
                                dr_jk[2] = coord[t][k][2] - coord[t][j][2];
                                df_jk[0] = inv_lattice[0][0]*dr_jk[0] + inv_lattice[1][0]*dr_jk[1] + inv_lattice[2][0]*dr_jk[2];
                                df_jk[1] = inv_lattice[0][1]*dr_jk[0] + inv_lattice[1][1]*dr_jk[1] + inv_lattice[2][1]*dr_jk[2];
                                df_jk[2] = inv_lattice[0][2]*dr_jk[0] + inv_lattice[1][2]*dr_jk[1] + inv_lattice[2][2]*dr_jk[2];

                                df_jk[0] -= round(df_jk[0]); df_jk[1] -= round(df_jk[1]); df_jk[2] -= round(df_jk[2]);

                                dr_pbc_jk[0] = current_lattice[0][0]*df_jk[0] + current_lattice[1][0]*df_jk[1] + current_lattice[2][0]*df_jk[2];
                                dr_pbc_jk[1] = current_lattice[0][1]*df_jk[0] + current_lattice[1][1]*df_jk[1] + current_lattice[2][1]*df_jk[2];
                                dr_pbc_jk[2] = current_lattice[0][2]*df_jk[0] + current_lattice[1][2]*df_jk[1] + current_lattice[2][2]*df_jk[2];
                                double distSq_jk = dr_pbc_jk[0]*dr_pbc_jk[0] + dr_pbc_jk[1]*dr_pbc_jk[1] + dr_pbc_jk[2]*dr_pbc_jk[2];

                                if (distSq_jk < ta_cl_cutoff_sq) {
                                    is_cl_bridged = 1;
                                    break; 
                                }
                            } // end k (other Ta) loop
                            if (is_cl_bridged) {
                                ta_cl_ta_bridges++;
                            }
                        } // if (distSq_ij < ta_cl_cutoff_sq)
                    } // else if (type_j == 4)
                } // end j (all atoms) loop

                // Polymeric case
                fprintf(fp, "%d\t%d\t%d\t%d\t%d\n", t, i, ta_o_ta_bridges, ta_cl_ta_bridges, ta_type);

            } else {
                // Isolated case
                fprintf(fp, "%d\t%d\t0\t0\t%d\n", t, i, ta_type);
            }

            ta_idx_counter++; // Increment for every Ta atom processed
        } // end i (Ta atom) loop

        free(global_to_o_idx);
    } // end t (timestep) loop

    fclose(fp);
}

void analyze_oxygen_classification(const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening oxygen classification output file");
        return;
    }
    fprintf(fp, "# Timestep\tOxygen_Global_Index\tOxygen_Local_Index\tOxygen_Type(0=NBO, 1=BO, 2=UO)\n");

    printf("\tAnalyzing and writing Oxygen Classification distribution...\n");
    for (int t = 0; t < numTraj; t++) {
        if (t % 100 == 0) printf("\t Writing O-Type frame t = %d...\n", t);
        int nAtom = numAtoms[t];
        int o_idx_counter = 0;
        for (int i = 0; i < nAtom; i++) {
            if (atom[t][i][1] == 3) { // Is Oxygen
                int o_type = oxygen_classification[t][o_idx_counter];
                // (t, Global Atom ID, Local O-only ID, O-Type)
                fprintf(fp, "%d\t%d\t%d\t%d\n", t, i, o_idx_counter, o_type);
                o_idx_counter++;
            }
        }
    }
    fclose(fp);
}