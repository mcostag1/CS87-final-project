/* 
 * SZipper - Merges point clouds by finding an appropriate
 * rotation and translation matrix to relate two overlapping
 * point clouds
 * by: Simon Bloch and Martina Costagliola
 * date: 5/12/2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>    

#include <sys/time.h>
#include <math.h>

#include "lib/eig3.h"
#include "lib/pca.h"
#include "lib/svd.h"

#define printhere printf("HERE %d\n", __LINE__)
#define NUMPOINTS 10
#define HIST_PROPORTION .85
#define CANDIDATES 10

#define DIST_THRESH 0.02

#define U 7
#define W (U - 1)/2
#define K_POINTS U*U

////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Struct Definitions //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

typedef struct pt {   
    int valid;
    float x;
    float y;
    float z;
}pt;

typedef struct candidate{
    float score;
    int i;
    int j;

}candidate;

////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Function Prototypes //////////////////////////////
////////////////////////////////////////////////////////////////////////////////


float* generate_histograms(pt* points_list, pt* normals, int n, int m);

float** get_rough_transform(candidate* scores, pt* A, pt* B, int n, int m, int n2, int m2, int pair);


pt* parse_ply_pts(char* infile_name, float*** vertices, int* n, int* m, int* v);

pt* run_pca(pt* points, int n, int m, int v);


candidate* score_hists(float* cur_hists, float* base_hists, int n, int m, int n2, int m2);

pt map(int ind, int i, int j, int n, int m, pt* mat);

pt* p_map(int ind, int i, int j, int n, int m, pt* mat);

void scale_hist(float* histograms, int i, int j, int m, int point_pair_tally);

float hist_dist(float* hA, float* hB, int i, int j);


void write_ply(float** R, pt* pts, char* outfile_name, int n, int m);

void print_mat(float** M);

void clean_normals(pt* normals, int n, int m);
     
void consider_candidate(float score, int i, int j, candidate* scores);
float Dot(pt A, pt B);
pt Cross(pt A, pt B);
pt Sub(pt A, pt B);
float Norm(pt A);
float Max(float* list, int len);

float Min(float* list, int len);
int hashFeatures(float* f, float* s);

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Main Function /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {


    // grab command line arguments
    char* infile_1_name = argv[1];
    char* infile_2_name = argv[2];
    char* outfile_1_name = argv[3];
    char* outfile_2_name = argv[4];

    // variable declarations
    int n, m, v;
    int n2, m2, v2;

    pt* cur_list;
    pt* base_list;
    pt* cur_normals = NULL;
    pt* base_normals = NULL;

    float** vertices;

    // parse ply points for both 3D point clouds
    cur_list = parse_ply_pts(infile_1_name, &vertices, &n, &m, &v);
    base_list = parse_ply_pts(infile_2_name, &vertices, &n2, &m2, &v2);

    //printf("vertices: %f %f %f\n\n", vertices[20][0], vertices[20][1], vertices[20][2]); 

    // error check to see if 3D vertices were grabbed correctly
    // from both scans

    if (cur_list==NULL) {
        printf("Cur list is NULL\n");
    }

    if (base_list==NULL) {
        printf("base list is NULL\n");
    }
    
    // grab normals for each point cloud
    cur_normals = run_pca(cur_list, n, m, v);
    base_normals = run_pca(base_list, n2, m2, v2);

    // orient the normals so that they are all facing the same direction
    // either all out of the surface, or all in towards the surface
    clean_normals(cur_normals, n, m);
    clean_normals(base_normals, n2, m2);

    float* cur_hists, *base_hists;

    // generate histograms for both point clouds
    cur_hists = generate_histograms(cur_list, cur_normals, n, m);
    base_hists = generate_histograms(base_list, base_normals, n2, m2);

    int i;
    int j;

    candidate* scores;


    // score the histograms - compare each point in scan A to a point in scan B 
    // determine the top 10 point correspondences
    scores = score_hists(cur_hists, base_hists, n, m, n2, m2);

    float** R;

    char c;
    for (i=0; i<CANDIDATES; i++) {

        if (i==0) {
            c = '0';
        }
        if (i==1) {
            c = '1';
        }
        if (i==2) {
            c = '2';
        }
        if (i==3) {
            c = '3';
        }
        if (i==4) {
            c = '4';
        }
        if (i==5) {
            c = '5';
        }
        if (i==6) {
            c = '6';
        }
        if (i==7) {
            c = '7';
        }
        if (i==8) {
            c = '8';
        }
        if (i==9) {
            c = '9';
        }

        outfile_1_name[1] = c;
        outfile_2_name[1] = c;

        // compute rotation matrix that maps scan A to scan B based on
        // the top 10 point correspondences
        R = get_rough_transform(scores, cur_list, base_list, n, m, n2, m2, i);

        // write the ply files given the rotation matrix
        write_ply(R, cur_list, outfile_1_name, n, m);

        // write the ply files given the rotation matrix
        write_ply(R, base_list, outfile_2_name, n, m);



    }

    
    // handle memory management
    // free all malloced items
    free(cur_list);
    free(base_list);
    free(cur_normals);
    free(base_normals);

    /* // free things */
    /* for (i=0; i<v; i++) { */
    /*     free(vertices[i]); */
    /* } */
    /* for (i=0; i<n; i++) {     */
    /*     free(points_list[i]); */
    /*     for (j=0; j<m; j++) { */
    /*         free(normals[i][j]); */
    /*     } */
    /*     free(normals[i]); */
    /* } */
    /* free(points_list); */
    /* free(normals);          // free(vertices); */
            
    return 0;
    
}

/* clean_normals - orient all normals at each point 
 * of a surface to face the same direction
 */

void clean_normals(pt* normals, int n, int m) {

    int cur, temp;

    cur = 0;

    while (!normals[cur].valid) {
        cur++;
    }

    printf("CUR: %d\n%d\n\n%d, %d \n\n", cur, n*m, n, m);

    int i;
    for (i=1; i < n*m; i++) {

        if (normals[i].valid) {
            if ((cur % m) > (i % m)) {
                if (i - m >= 0) {
                    cur = i - m;
                    
                    while (!normals[cur].valid) {
                        cur++;
                    }
                }
            }

            if (Dot(normals[cur], normals[i]) < 0) {

                normals[i].x *= -1;
                normals[i].y *= -1;
                normals[i].z *= -1;

                cur = i;
            }
        }
    }
}

/* write_ply - given an outfile name, points, size (nxm), and rotation
 * matrix, write out a PLY file with the new 3D point cloud 
 */

void write_ply(float** R, pt* pts, char* outfile_name, int n, int m) {
    FILE* out;
    out = fopen(outfile_name, "w");
    
    float p_in[3];
    float p_out[3];


    int i, j, k, ind=0;
    for (i=0; i < n*m; i++) {

        if (pts[i].valid) {
            p_in[0] = pts[i].x;
            p_in[1] = pts[i].y;
            p_in[2] = pts[i].z;

            memset(p_out, 0, 3*sizeof(float));
            for (j=0; j < 3; j++) {
                for (k=0; k < 3; k++) {
                    p_out[j] += R[j][k] * p_in[k];
                }
            }
            
            fprintf(out, "%f %f %f\n", p_out[0], p_out[1], p_out[2]);
            if (p_out[2] > 0.1) {
                pts[i].z -= .2;
            }
        }

    }
    for (i=0; i < n*m; i++) {
        if (pts[i].valid) {
            fprintf(out, "1 %d\n", ind);
            ind++;
        } else {
            fprintf(out, "0\n");
        }
    }
}

/* get_rough_transform - given two point clouds, and point correspondences,
 * return a rotation matrix that maps one point cloud to the orientation of
 * the other.
 * 
 */
float** get_rough_transform(candidate* scores, pt* A, pt* B, int n, int m, int n2, int m2,
                            int pair) {

    int i, j, k, ii, jj;

    float centroidA[3];
    float centroidB[3];

    memset(centroidA, 0, 3*sizeof(float));
    memset(centroidB, 0, 3*sizeof(float));

    int indA, indB;

    pt* point, *pointA, *pointB;

    for (i = 0; i < CANDIDATES; i++ ){

        for (j=0; j < K_POINTS; j++) {

            k = scores[i].i/16;

            ii= k/m;
            jj= k%m;

            point = p_map(j, ii, jj, n, m, A);
            if (point->valid) {
                if (i==pair) {
                    point->z += .2;                
                }
            }




            k = scores[i].j/16;

            ii= k/m;
            jj= k%m;

            point = p_map(j, ii, jj, n, m, B);
            if (point->valid) {
                if (i==pair) {
                    point->z += .2;
                }       
            }
        }

        /*
        for (j=0; j < K_POINTS; j++) {

            k = scores[i].i/16;

            ii= k/m;
            jj= k%m;

            pt point = map(j, ii, jj, n, m, A);

            if (i == pair) {
                int a, b;

                a = j / U;
                b = j % U;

                a = ii - W + a;
                b = jj - W + b;

                if ((a < 0) ||
                    (b < 0) ||
                    (a >= n) ||
                    (b >= m)) {
                    jj++;
                    jj--;

                } else {
                    A[a*m+b].z += .2;
                }

            }

            k = scores[i].j/16;

            ii= k/m;
            jj= k%m;

            point = map(j, ii, jj, n, m, B);

            // ^
            if (i == pair) {
                int a, b;

                a = j / U;
                b = j % U;

                a = ii - W + a;
                b = jj - W + b;

                if ((a < 0) ||
                    (b < 0) ||
                    (a >= n) ||
                    (b >= m)) {
                    jj++;
                    jj--;

                } else {
                    B[a*m+b].z += .2;
                }
            }

        }
        */
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }
    centroidA[0] /= CANDIDATES;
    centroidA[1] /= CANDIDATES;
    centroidA[2] /= CANDIDATES;
    centroidB[0] /= CANDIDATES;
    centroidB[1] /= CANDIDATES;
    centroidB[2] /= CANDIDATES;

    float** H;
    H = malloc(3*sizeof(float*));
    for (i=0; i < 3; i++) {
        H[i] = malloc(3*sizeof(float));
        memset(H[i], 0, (3*sizeof(float)));
    }

    float ptA[3];
    float ptB[3];

    for (i=0; i < CANDIDATES; i++) {
        ptA[0] = A[scores[i].i/16].x;
        ptA[1] = A[scores[i].i/16].y;
        ptA[2] = A[scores[i].i/16].z;

        ptB[0] = B[scores[i].j/16].x;
        ptB[1] = B[scores[i].j/16].y;
        ptB[2] = B[scores[i].j/16].z;

        for (j=0; j < 3; j++) {
            for (k=0; k < 3; k++) {
                H[j][k] += (ptA[j] - centroidA[j]) * (ptB[k] - centroidB[k]);
            }
        }
    }
    
    float* w;
    float** v;
    float** R;

    w = malloc(3*sizeof(float));
    v = malloc(3*sizeof(float*));
    for (i=0; i<3; i++) {
        v[i] = malloc(3*sizeof(float));
        memset(v[i], 0, 3*sizeof(float));
    }

    R = malloc(3*sizeof(float*));
    for (i=0; i<3; i++) {
        R[i] = malloc(3*sizeof(float));
        memset(R[i], 0, 3*sizeof(float));
    }


    int ret;
    printf("H Before\n");
    print_mat(H);

    ret = dsvd(H, 3, 3, w, v);

    printf("H After\n");
    print_mat(H);

    printf("v Matrix\n");
    print_mat(v);

    printf("RET: %d\n", ret);

    R[0][0] = 1;
    R[1][1] = 1;
    R[2][2] = 1;

    /*

    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            for (k=0; k < 3; k++) {
                R[i][j] += v[i][k] * H[j][k];
            }
        }
        }*/

    printf("R Matrix\n");
    print_mat(R);
    return R;
}


/* print_mat - helper function to print matrices
 *
 */
void print_mat(float** M) {
    int i;
    int j;

    printf("+ + + + + + + + + \n");
    printf("MATRIX:\n");
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            printf("%f ", M[i][j]);
        }
        printf("\n");
    }
    printf("- - - - - - - - - \n\n\n");

}

/* hist_dist - determine the similarity (or score)
 * of two given histograms.
 *
 */
float hist_dist(float* hA, float* hB, int i, int j) {
    int k;
    float sum;

    if ((hA[i] < 0) || (hB[j] < 0)) {
        return 2;
    }

    sum = 0;
    for (k=0; k<16; k++) {
        sum += (hA[i + k] - hB[j + k]) * (hA[i + k] - hB[j + k]);
    }
    
    return sqrt(sum);
}

/* score_hists - given two sets of histograms that describe two point clouds
 * find the ten best point correspondences.
 */

candidate* score_hists(float* cur_hists, float* base_hists, int n, int m, int n2, int m2) {

    int i, j;
    float score;
    long N = (long)(n*m);
    long M = (long)(n2*m2);

    long len = N*M;

    candidate* scores;

    scores = (candidate*)malloc((CANDIDATES*sizeof(candidate)));

    for(i=0; i < CANDIDATES; i++){
        scores[i].score = 2;
    }
    
    printf("%li\n", len);

    float p = .01;

    int unit = n*m*16/10;

    
    for (i=0; i < n*m*16; i += 16) {
        for (j=0; j < n*m*16; j += 16) {

            score = hist_dist(cur_hists, base_hists, i, j);

            consider_candidate(score, i, j, scores);
            
            if ((((long)i/16) * ((long)j/16)) > p*len){
                printf("%f percent:  S--> %f\n", 100*p, score);
                p+=.01;
            }
        }
    }
    
    return scores;

}

/* generate_histograms - compute the feature histogram that describes
 * each point in a given 3D point cloud.
 *
 */
float* generate_histograms(pt* points_list, pt* normals, int n, int m) {

    int i, j, k, l , idx;

    // Dependent on indices
    pt pt_k;
    pt pt_l;
    pt normal_k;
    pt normal_l;

    // Independent of indices
    pt pt_s;
    pt pt_t;
    pt normal_s;
    pt normal_t;
    pt u;
    pt v;
    pt w;

    float f[4], s[4], distances[4];

    memset(s, 0, 4*sizeof(float));

    float* histograms;

    histograms = (float*) malloc(m * n * 16 * sizeof(float));

    memset(histograms, 0, m * n * 16 * sizeof(float));

    int point_pair_tally;

    for (i=0; i< n; i++) {
        for (j=0; j<m; j++) {
            point_pair_tally = 0;

            // AT i, j
            memset(distances, 0, 4*sizeof(float));

            if ((i > W) && (j > W)) {
                distances[0] = Norm(Sub(points_list[(i*m) + j], points_list[(i-W)*m + (j-W)]));
            }
            if ((i > W) && (j + W < m)) {
                distances[1] = Norm(Sub(points_list[(i*m) + j], points_list[(i-W)*m + (j+W)]));
            }
            if ((i + W < n) && (j > W)) {
                distances[2] = Norm(Sub(points_list[(i*m) + j], points_list[(i+W)*m + (j-W)]));
            }
            if ((i + W < n) && (j + W < m)) {
                distances[3] = Norm(Sub(points_list[(i*m) + j], points_list[(i+W)*m + (j+W)])); 
            }
            
            s[1] = Max(distances, 4);

            /*
              We need neighboring points, so we can index into points list
              we need neighboring normals, so we can index into normals list
            */
            for (k = 0; k < K_POINTS - 1; k++) {
                for (l = k + 1; l < K_POINTS; l++) {

                    pt_k = map(k, i, j, n, m, points_list);
                    pt_l = map(l, i, j, n, m, points_list);
                    normal_k = map(k, i, j, n, m, normals);
                    normal_l = map(l, i, j, n, m, normals);

                    if (pt_k.valid &&
                        pt_l.valid &&
                        normal_k.valid &&
                        normal_l.valid) {

                        //Cross Dot Sub
                        if (Dot(normal_k, Sub(pt_l, pt_k)) <=
                            Dot(normal_l, Sub(pt_k, pt_l))) {
                        
                            pt_s = pt_k;
                            pt_t = pt_l;
                            normal_s = normal_k;
                            normal_t = normal_l;
                        } else {
                            pt_s = pt_l;
                            pt_t = pt_k;
                            normal_s = normal_l;
                            normal_t = normal_k;
                        }
                        u = normal_s;
                        v = Cross(Sub(pt_t, pt_s), u);
                        w = Cross(u, v);

                        f[0] = Dot(v, normal_t);
                        f[1] = Norm(Sub(pt_t, pt_s));
                        f[2] = Dot(u, Sub(pt_t, pt_s)) / f[1];
                        f[3] = atan2(Dot(w, normal_t), Dot(u, normal_t));


                        // TODO: make f array, make s array, evaluate step function, get hist!
                        /*
                          make histogram from fs and ss
                          make histogram for other point cloud
                          get point correspondences
                          get rigid transform
                          parallelize

                        */

                        idx = hashFeatures(f, s);
                        if ((idx < 0) || (idx > 15)) {
                            printf("EXITING: IMPROPER HASH\n");
                            exit(1);
                        }

                        histograms[((i*m + j) * 16) + idx]++;
                        point_pair_tally++;
                    }
                }
            } // Back out, @ (i,j)
            
            scale_hist(histograms, i, j, m, point_pair_tally);

            /*
              if (histograms[((i*m) + j)*16] >= 0) {
                
              printf("Histogram: [");
              for (k=0;k<16;k++) {
              printf("%f, ", histograms[(i*m + j) * 16 + k]);
              }
              printf("]\n");
                
              printf("(%d, %d)\n\n", i, j);
              }
            */
            
        }
    }

    return histograms;
}


/* run_pca - given a file describing a 3D point cloud,
 * run Principal Component Analysis to determine the normal
 * vector for each point in the 3D point cloud.
 *
 */
pt* run_pca(pt* points, int n, int m, int v) {
    
    int i, j, k, l;
    
    //PARALLELIZE HERE
    //
    //Generate kernel, 25x3

    pt* normals;
    normals = (pt*) malloc((n) * (m) * sizeof(pt));

    
    
    float** kernel;
    kernel = (float**) malloc(K_POINTS * sizeof(float*));

    int ind;

    float mean[3];
    float V[3][3];
    float d[3];

    for (i=0; i < n; i++) {
        for (j=0; j<m; j++) {
            // AT i, j
            ind = 0;
            
            normals[i*m+j].valid = 1;       
            

            for (k = i - 3; k < i + 4; k++) {
                for (l = j - 3; l < j + 4; l++) {
                    if ((k >= 0) && (k < n) && (l >= 0) && (l < m)) {
                        if (points[k * m + l].valid != 0) {
                            if (Norm(Sub(points[i*m+j],
                                         points[k*m + l])) < DIST_THRESH) {

                                kernel[ind] = (float*) malloc(3 * sizeof(float));

                                //printf(" points: %f %f %f\n", points[k * m + l].x,
                                //                              points[k * m + l].y,
                                //                              points[k * m + l].z );

                                kernel[ind][0] = points[k * m + l].x;
                                kernel[ind][1] = points[k * m + l].y;
                                kernel[ind][2] = points[k * m + l].z;
                                ind++;
                            }
                        }                        
                    }    
                } // Out of loop across kernel
            }
            if (ind >= 3) {
                PCA(ind, kernel, mean, V, d);

                normals[i*m+j].x = V[0][0];
                normals[i*m+j].y = V[1][0];
                normals[i*m+j].z = V[2][0];
                    
                // printf("i: %d\nj: %d\ne: %d\n", i, j, e);
                // printf("normals[i].x: %f\n", normals[i*m+j].x);                            
                // printf("normals[i].y: %f\n", normals[i*m+j].y);
                // printf("normals[i].z: %f\n", normals[i*m+j].z);
            } else {
                normals[i*m+j].valid = 0;
            }
        }
    }

    // for (i=0; i<n; i++) {
    //     for (j=0; j<m; j++) {
    //        printf("i: %d\nj: %d\n", i, j);
    //        printf("NORMAL: %f %f %f\n", normals[i*m+j].x, normals[i*m+j].y, normals[i*m+j].z);
    //     }
    // }
    

    //free(kernel);
    

    return normals;
    
    //Pass Kernel into PCA
    //get return eigenvector
    //populate output matrix with normal
}

/* parse_ply_pts - given a PLY file, parse and read the vertices
 * and indeces for those vertices, and store them in an array of point
 * structs.
 */
pt* parse_ply_pts(char* infile_name, float*** vertices, int* n, int* m, int* v) {
    
    FILE *stream;

    if ((stream = fopen(infile_name,"r")) == NULL) {
        fprintf(stderr, "p_zipper: cannot open file %s\n", infile_name);
        fprintf(stderr, "Exiting to system.");
        exit(1);
    }

    //18 then cols then 2 then rows
    char junk_val[20];
    int i;

    // Step thru irrelevant header data til #cols
    for (i = 0; i < 18; i++) {
        fscanf(stream, "%s", junk_val);
    }

    // Read #cols
    fscanf(stream, "%d", m);

    // More junk
    fscanf(stream, "%s", junk_val);
    fscanf(stream, "%s", junk_val);

    // Read #rows
    fscanf(stream, "%d", n);

    // More junk
    for (i = 0; i < 29; i++) {
        fscanf(stream, "%s", junk_val);
    }

    // Read #vertices
    fscanf(stream, "%d", v);

    // Finish reading junk
    for (i = 0; i < 18; i++) {
        fscanf(stream, "%s", junk_val);
    }

    // Now read in data.
    pt *data = NULL;

    *vertices = (float**) malloc(*v * sizeof(float *));
    float** v_array = *vertices;

    for (i = 0; i < *v; i++) {

        v_array[i] = (float*) malloc(3 * sizeof(float));

        fscanf(stream, "%f %f %f", &v_array[i][0], &v_array[i][1], &v_array[i][2]);

    }

    int is_vert, index;


    // x coordinate: 2*i
    // y coordinate: 2*i+1
    // int j;
    data = (pt*) malloc((*n) * (*m) * sizeof(pt));
    for (i = 0; i < (*n) * (*m) ; i++) {
        fscanf(stream, "%d", &is_vert);
        if (is_vert) {
            fscanf(stream, "%d", &index);
            data[i].x = v_array[index][0];
            data[i].y = v_array[index][1];
            data[i].z = v_array[index][2];
                
            free(v_array[index]);

            data[i].valid = 1;         

        } else {

            data[i].valid = 0;
        }
    }
    free(v_array);

    return data;
}
	
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Helper Functions ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Maps index to its corresponding vertex/normal
pt map(int ind, int i, int j, int n, int m, pt* mat) {

    pt invalid;
    invalid.valid = 0;

    // Have: 
    //   ind
    //   U (width of kernel)
    //   W (width of side of kernel)
    //   NB/2 (num candidate neighbors, not nec. valid)
    //   n & m
    int a, b;

    a = ind / U;
    b = ind % U;

    // a, b -- Now span from 0 to (width of kernel - 1)
    // now we need to map a and b to i,j indexes
    // 
    // if a is 0, we want i - W
    // if b is 0, we want j - W
    a = i - W + a;
    b = j - W + b;

    if ((a < 0) ||
        (b < 0) ||
        (a >= n) ||
        (b >= m)) {
        return invalid;
    }

    return mat[a*m + b];
}

// Maps index to its corresponding vertex/normal
pt* p_map(int ind, int i, int j, int n, int m, pt* mat) {

    pt invalid;
    invalid.valid = 0;

    // Have: 
    //   ind
    //   U (width of kernel)
    //   W (width of side of kernel)
    //   NB/2 (num candidate neighbors, not nec. valid)
    //   n & m
    int a, b;

    a = ind / U;
    b = ind % U;

    // a, b -- Now span from 0 to (width of kernel - 1)
    // now we need to map a and b to i,j indexes
    // 
    // if a is 0, we want i - W
    // if b is 0, we want j - W
    a = i - W + a;
    b = j - W + b;

    if ((a < 0) ||
        (b < 0) ||
        (a >= n) ||
        (b >= m)) {
        return &mat[0];
    }

    return &mat[a*m + b];
}

// Dot - performs dot product 
float Dot(pt A, pt B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

// Cross - performs cross product
pt Cross(pt A, pt B) {
    pt C;

    C.valid = 1;
    C.x = A.y * B.z - A.z * B.y;
    C.y = A.z * B.x - A.x * B.z;
    C.z = A.x * B.y - A.y * B.x;

    return C;
}

// Sub - subtracts one point from another
pt Sub(pt A, pt B) {

    B.x = A.x - B.x;
    B.y = A.y - B.y;
    B.z = A.z - B.z;

    return B;
}

// Norm - computes norm of a function
float Norm(pt A) {
    return sqrt((A.x * A.x) + (A.y * A.y) + (A.z * A.z));
}

// Max - determines the max of a list

float Max(float* list, int len) {
    
    int i, max = 0;
    for (i = 1; i < len; ++i){
        if (list[i] > list[max]){
            max = i;            
        }
    }
    return list[max];
}

float Min(float* list, int len) {

    int i, min = 0;
    for (i = 1; i < len; ++i){
        if (list[i] < list[min]){
            min = i;
        }
    }
    return min;
}

int hashFeatures(float* f, float* s) {

    int i, idx = 0;
    for (i = 0; i < 4; i++){
        if (!(f[i] < s[i])){
            idx += pow(2.0,i);
        }
    }

    return idx;
}


// scale_hist -- sorts the
void scale_hist(float* histograms, int i, int j, int m, int point_pair_tally) {
    int possibles = K_POINTS * (K_POINTS - 1)/2;
    
    int tally_threshold = possibles * HIST_PROPORTION;
        
    int p;

    if (point_pair_tally < tally_threshold) {
        histograms[(i*m + j)*16] = -1;
    } else {
        for (p=(i*m+j)*16; p < (i*m+j)*16 + 16; p++) {
            histograms[p] /= point_pair_tally;
        }
    }
}



void consider_candidate(float score, int i, int j, candidate* scores){
    int cur_ind = CANDIDATES-1;

    candidate temp;
    candidate cur;

    cur.score = score;
    cur.i = i;
    cur.j = j;

    if (cur.score < scores[CANDIDATES-1].score){
        scores[cur_ind] = cur;

        while(scores[cur_ind].score < scores[cur_ind -1].score){      
            temp = scores[cur_ind];
            scores[cur_ind] = scores[cur_ind-1];
            scores[cur_ind-1] = temp;

            if(cur_ind == 1){
                break;
            }
            cur_ind--;            
        }
    }
}

