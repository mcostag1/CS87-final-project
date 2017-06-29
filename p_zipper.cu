/* 
 * PZipper - Merges point clouds and zippers in parallel using CUDA
 *
 * by: Simon Bloch and Martina Costagliola
 * date: 4/7/2016
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>    
#include <curand_kernel.h>
#include <curand.h>
#include <sys/time.h>
#include <math.h>

#include "lib/eig3.h"
#include "lib/pca.h"
#include "lib/svd.h"
#include "lib/myopengllib.h"

#define printhere printf("HERE %d\n", __LINE__)
#define NUMPOINTS 1
#define HIST_PROPORTION .85
#define CANDIDATES 10

#define DIST_THRESH 0.02


#define CUDA_CALL(x) {cudaError_t cuda_error__ = (x); if (cuda_error__) printf("CUDA error: " #x " returned \"%s\"\n", cudaGetErrorString(cuda_error__));}
#define U 7
#define W (U - 1)/2
#define K_POINTS U*U

#define PARTITION 8
#define P (PARTITION*PARTITION)


////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Struct Definitions //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

typedef struct point {   
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


typedef struct scan_data_state {
  
    pt* base_scan;
    pt* cur_scan;
    pt* cur_normals;
    pt* base_normals;
    float* histograms_A;
    float* histograms_B;
    int* tally_points_A;
    int* tally_points_B;
    candidate* chunk_scores;
}scan_data_state;

scan_data_state scan_data;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////// Function Prototypes //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// global function prototypes

__global__ void generate_histograms(pt* points_list, pt* normals, float* histograms, int* point_pair_tally,
                                    int n, int m, int k, int l);

__global__ void scale_hist(float* histograms, int m, int* tally_points);


__global__ void score_hists(float* cur_hists, float* base_hists, candidate* full_scores,
                            int n, int m, int n2, int m2, int a_ind, int b_ind);

__global__ void rec_sort(candidate* scores, int len, int count);

// device functions prototypes
__device__ void scale_hist(float* histograms, int i, int j, int m, int point_pair_tally);
__device__ float hist_dist(float* hA, float* hB, int i, int j);
__device__ void consider_candidate(float score, int i, int j, candidate* scores);
__device__ int hashFeatures(float* f, float* s);
__device__ pt map(int ind, int i, int j, int n, int m, pt* mat);
__device__ void Drec_sort(candidate* scores, int len, int count);

 __device__ pt Cross(pt A, pt B);
 __device__ pt Sub(pt A, pt B);
 __device__ float Norm(pt A);
 __device__ float Max(float* list, int len);
 __device__ float Dot(pt A, pt B);
 __device__ float Min(float* list, int len);

candidate* generate_scores(scan_data_state scan_data, int n, int m, int n2, int m2);

void generate_hists(scan_data_state scan_data, int n, int m, int n2, int m2);

void scale_histograms(float* histograms, int* tally_points, int n, int m);

void clean_normals(pt* normals, int n, int m);

// prototypes for normal c functions

void write_ply(float** R, pt* pts, char* outfile_name, int n, int m);

float** get_rough_transform(candidate* scores, pt* A, pt* B, int n, int m, int n2, int m2);

pt* run_pca(pt* points, int n, int m, int v);

pt* parse_ply_pts(char* infile_name, float*** vertices, int* n, int* m, int* v);

void Crec_sort(candidate* scores, int len, int count);

void print_mat(float** M);

void merge(candidate* scores, candidate* input);

float CDot(pt A, pt B);
pt CSub(pt A, pt B);
float CNorm(pt A);
pt Cmap(int ind, int i, int j, int n, int m, pt* mat);

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Main Function /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {


    // grab input files
    char* infile_1_name = argv[1];
    char* infile_2_name = argv[2];
    char* outfile_1_name = argv[3];

    // declare files
    int n, m, v;
    int n2, m2, v2;

    pt* cur_list;
    pt* base_list;
    pt* cur_normals = NULL;
    pt* base_normals = NULL;

    float** vertices;
    //float* histograms;

    // Parse the ply files into an 1D, length (m*n) array
    cur_list = parse_ply_pts(infile_1_name, &vertices, &n, &m, &v);
    base_list = parse_ply_pts(infile_2_name, &vertices, &n2, &m2, &v2);

    if (cur_list==NULL) {
        printf("Cur list is NULL\n");
    }

    if (base_list==NULL) {
        printf("base list is NULL\n");
    }

    // Generate normal vectors at each point
    cur_normals = run_pca(cur_list, n, m, v);
    base_normals = run_pca(base_list, n2, m2, v2);

    // Iterate over normals (sequentially) checking that they all point the same direction
    clean_normals(cur_normals, n, m);
    clean_normals(base_normals, n2, m2);

    /* Malloc on GPU space */
    HANDLE_ERROR(cudaMalloc((void **)&scan_data.cur_scan,
                             sizeof(pt) * n * m), "malloc cur_scan");

    HANDLE_ERROR(cudaMalloc((void **)&scan_data.base_scan,
                            n2 * m2 * sizeof(pt)), "malloc base_scan");

    HANDLE_ERROR(cudaMalloc((void **)&scan_data.base_normals,
                            n2 * m2 * sizeof(pt)), "malloc space for base normal vectors");
    
    HANDLE_ERROR(cudaMalloc((void **)&scan_data.cur_normals,
                            n * m * sizeof(pt)), "malloc space for cur normal vectors");
 

    HANDLE_ERROR(cudaMalloc((void **)&scan_data.histograms_A,
                                  n * m *16* sizeof(float)), "malloc space for cur hists");
  
    HANDLE_ERROR(cudaMalloc((void **)&scan_data.histograms_B,
                                  n2 * m2 *16* sizeof(float)), "malloc space for base hists");

    HANDLE_ERROR(cudaMalloc((void **)&scan_data.tally_points_A,
                                  n * m * sizeof(int)), "malloc space for cur hist point tallies");

    HANDLE_ERROR(cudaMalloc((void **)&scan_data.tally_points_B,
                                  n2 * m2 * sizeof(int)), "malloc space for base hist point tallies");
    
    /* Copy over to the GPU*/
    HANDLE_ERROR(cudaMemcpy(scan_data.cur_scan, cur_list,
                            n * m * sizeof(pt), cudaMemcpyHostToDevice),
                 "copy cur points list to GPU");

    HANDLE_ERROR(cudaMemcpy(scan_data.base_scan, base_list,
                            n2 * m2 * sizeof(pt), cudaMemcpyHostToDevice),
                 "copy base points list to GPU");

    HANDLE_ERROR(cudaMemcpy(scan_data.cur_normals, cur_normals,
                                      n * m * sizeof(pt), cudaMemcpyHostToDevice),
                         "copy cur normals list to GPU");

    HANDLE_ERROR(cudaMemcpy(scan_data.base_normals, base_normals,
                            n2 * m2 * sizeof(pt), cudaMemcpyHostToDevice),
                 "copy base normals list to GPU");


    //////////////////////////////////////////////////////
    // PARALLELIZATION BEGINS!!! /////////////////////////
    //////////////////////////////////////////////////////

    cudaMemset(scan_data.histograms_A, 0, m * n * 16 * sizeof(float));    
    cudaMemset(scan_data.histograms_B, 0, m2 * n2 * 16 * sizeof(float));    
    cudaMemset(scan_data.tally_points_A, 0, m * n * sizeof(int));    
    cudaMemset(scan_data.tally_points_B, 0, m2 * n2 * sizeof(int));    

    // Generate the histograms for both cur and base scans
    //     (Generate hists populates scan_data with both histograms)
    generate_hists(scan_data, n, m, n2, m2);


    candidate* scores;
    scores = generate_scores(scan_data, n, m, n2, m2);

    // compute the rotation matrix that relates the two scans
    float** R;
    R = get_rough_transform(scores, cur_list, base_list, n, m, n2, m2);

    // write the completed ply files
    write_ply(R, cur_list, outfile_1_name, n, m);

    // free allocated memory on both the CPU and GPU
    free(cur_list);
    free(base_list);
    free(cur_normals);
    free(base_normals);
    free(scores);
    free(R);

    cudaFree(scan_data.cur_scan);
    cudaFree(scan_data.cur_normals);
    cudaFree(scan_data.histograms_A);
    cudaFree(scan_data.base_scan);
    cudaFree(scan_data.base_normals);
    cudaFree(scan_data.histograms_B);
    cudaFree(scan_data.chunk_scores);
            
    return 0;
    
}

// Parse the points from each ply file
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


    // Read in list of vertices
    *vertices = (float**) malloc(*v * sizeof(float *));
    float** v_array = *vertices;

    for (i = 0; i < *v; i++) {

        v_array[i] = (float*) malloc(3 * sizeof(float));

        fscanf(stream, "%f %f %f", &v_array[i][0], &v_array[i][1], &v_array[i][2]);

    }

    int is_vert, index;


    // Read vertices into data depending on whether each point is valid
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

// run PCA on sets of 3D points to aquire normals
pt* run_pca(pt* points, int n, int m, int v) {
    
    int i, j, k, l;
    
    pt* normals;
    normals = (pt*) malloc((n) * (m) * sizeof(pt));

    
    // Declare Kernel
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
            
            // At i,j loop across kernel and populate with kernel points
            for (k = i - 3; k < i + 4; k++) {
                for (l = j - 3; l < j + 4; l++) {
                    if ((k >= 0) && (k < n) && (l >= 0) && (l < m)) {
                        if (points[k * m + l].valid != 0) {
                            if (CNorm(CSub(points[i*m+j],
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

            // If we have enough points to approximate a plane, run PCA on the kernel points
            // to get normal vector V
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

    return normals;
}

// function to ensure the normals of each scan are all oriented in the same direction
void clean_normals(pt* normals, int n, int m) {

    int cur;

    cur = 0;

    while (!normals[cur].valid) {
        cur++;
    }


    // Loop across vectors
    int i;
    for (i=1; i < n*m; i++) {

        // If the vector is valid
        if (normals[i].valid) {
            if ((cur % m) > (i % m)) {


                // Look at the normal one index up in the matrix
                if (i - m >= 0) {
                    cur = i - m;
                    
                    while (!normals[cur].valid) {
                        cur++;
                    }
                }
            }

            // the dot product with the cur normal is negative, flip this one
            if (CDot(normals[cur], normals[i]) < 0) {

                normals[i].x *= -1;
                normals[i].y *= -1;
                normals[i].z *= -1;

                cur = i;
            }
        }
    }
}

// Populate scan_data with histograms
void generate_hists(scan_data_state scan_data, int n, int m, int n2, int m2){

    // declare error value for checking
    cudaError_t error;


    // specify thread and block dimensions
    dim3 blocks(n/PARTITION,m/PARTITION,1);
    dim3 threads_block(PARTITION, PARTITION, 1);

    int k, l;

    // Loop across kernel
    for(k = 0; k < K_POINTS - 1; k++){
        for(l = k + 1; l < K_POINTS; l++){
            error = cudaGetLastError();
            if (error != cudaSuccess) {
                printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
            }

            // Update histograms with the data from kernel index k,l in parallel
            generate_histograms<<<blocks, threads_block>>>(scan_data.cur_scan, scan_data.cur_normals,
                                                           scan_data.histograms_A, scan_data.tally_points_A,
                                                           n, m, k, l);

            // synchronize threads, make sure they are in sync before continuing to the next cuda kernel call
            error = cudaThreadSynchronize();
            if (error != cudaSuccess) {
                printhere;
                printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
            }
        }
    }

    scale_histograms(scan_data.histograms_A, scan_data.tally_points_A, n, m);
    for(k = 0; k < K_POINTS - 1; k++){
        for(l = k + 1; l < K_POINTS; l++){
            error = cudaGetLastError();
            if (error != cudaSuccess) {
                printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
            }

            // Call generate histograms for scan 2- second cuda kernel call
            generate_histograms<<<blocks, threads_block>>>(scan_data.base_scan, scan_data.base_normals,
                                                           scan_data.histograms_B, scan_data.tally_points_B,
                                                           n2, m2, k, l);

            // synchronize threads, make sure they are in sync before continuing to the next cuda kernel call
            error = cudaThreadSynchronize();
            if (error != cudaSuccess) {
                printhere;
                printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
            }
        }
    }

    // Use the updated tally_points array to normalize every histogram
    scale_histograms(scan_data.histograms_B, scan_data.tally_points_B, n2, m2);
}

// Normalize every histogram
void scale_histograms(float* histograms, int* tally_points, int n, int m) {
    
    cudaError_t error;

    // specify thread and block dimensions
    dim3 blocks(n/PARTITION,m/PARTITION,1);
    dim3 threads_block(PARTITION, PARTITION, 1);

    error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
    }

    // In parallel, normalize every histogram using the tally points array from scan data
    scale_hist<<<blocks, threads_block>>>(histograms, m, tally_points);

    error = cudaThreadSynchronize();
    if (error != cudaSuccess) {
        printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
    }

}

// Loops across scan A, compares hists with scan B in parallel, accumulates top 10 scores
candidate* generate_scores(scan_data_state scan_data, int n, int m, int n2, int m2) {

    cudaError_t error;

    dim3 blocks_2(PARTITION, 1);
    dim3 threads_block_2(PARTITION);

    dim3 blocks_3(1);
    dim3 threads_block_3(1);


    // Thus, idx goes from 0 to P-1

    // Chunk_scores is size 10*P
    HANDLE_ERROR(cudaMalloc((void **)&scan_data.chunk_scores,
                            CANDIDATES * P * sizeof(candidate)),
                 "malloc chunk_scores");

    int num_chunks;
    num_chunks = n*m/P;
    
    printf("n: %d\nm: %d\n", n, m);
    printf("NC: %d\n", num_chunks);
    printf("P: %d\n", P);


    candidate* scores;
    scores = (candidate*) malloc(CANDIDATES * sizeof(candidate));

    candidate* input;
    input = (candidate*) malloc(CANDIDATES * sizeof(candidate));

    candidate bad;
    bad.score = 2;

    int i, j;
    for (i=0; i < CANDIDATES; i++) {
        scores[i] = bad;
        input[i] = bad;
    }

    // Loop over chunks of A
    j = 0;
    for (i = 0; i < num_chunks; i++) {

        error = cudaGetLastError();
        if (error != cudaSuccess) {
            printhere;
            printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
        }

        score_hists<<<blocks_2, threads_block_2>>>(scan_data.histograms_A,
                                                   scan_data.histograms_B,
                                                   scan_data.chunk_scores,
                                                   n, m, n2, m2, i, j);
     
        // Now we have 10 * P scores for the first chunk, let's send these up to the CPU
   
        error = cudaThreadSynchronize();

        if (error != cudaSuccess) {
            printhere;
            printf("CudaThread Synchronize Error %s\n", cudaGetErrorString(error));
        }

            
        // Recursively sort the kernel output scores
        rec_sort<<<blocks_3, threads_block_3>>>(scan_data.chunk_scores, CANDIDATES * P, 0);

        // now do a device to host memcopy so that we can have access to the computed scores
        HANDLE_ERROR(cudaMemcpy(input,
                                scan_data.chunk_scores,
                                CANDIDATES * sizeof(candidate),
                                cudaMemcpyDeviceToHost),
                     "copy scores to CPU\n");

        // Merge the incoming array of 10 scores with the current top 10
        merge(scores, input);

    }

    free(input);

    return scores;
}


// Generate a rotation matrix from point correspondences
float** get_rough_transform(candidate* scores, pt* A, pt* B, int n, int m, int n2, int m2) {

    int i, j, k, ii, jj;

    float centroidA[3];
    float centroidB[3];

    memset(centroidA, 0, 3*sizeof(float));
    memset(centroidB, 0, 3*sizeof(float));


    pt point, pointA, pointB;
    pointA.valid = 0;
    pointB.valid = 0;

    // Get the centroid of every candidate
    for (i = 0; i < CANDIDATES; i++ ){

        for (j=0; j < K_POINTS; j++) {

            k = scores[i].i/16;

            ii= k/m;
            jj= k%m;

            point = Cmap(j, ii, jj, n, m, A);
            if (point.valid) {
                if (!pointA.valid || j <= (K_POINTS/2)) {
                    pointA = point;
                }
            }

            k = scores[i].j/16;

            ii= k/m;
            jj= k%m;

            point = Cmap(j, ii, jj, n, m, B);
            if (point.valid) {
                if (!pointB.valid || j <= (K_POINTS/2)) {
                    pointB = point;
                }
            }
        }
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

        centroidA[0] += pointA.x;
        centroidA[1] += pointA.y;
        centroidA[2] += pointA.z;

        centroidB[0] += pointB.x;
        centroidB[1] += pointB.y;
        centroidB[2] += pointB.z;

    }
    centroidA[0] /= CANDIDATES;
    centroidA[1] /= CANDIDATES;
    centroidA[2] /= CANDIDATES;
    centroidB[0] /= CANDIDATES;
    centroidB[1] /= CANDIDATES;
    centroidB[2] /= CANDIDATES;

    // H is the covariance matrix of points
    float** H;
    H = (float**)malloc(3*sizeof(float*));
    for (i=0; i < 3; i++) {
        H[i] = (float*)malloc(3*sizeof(float));
        memset(H[i], 0, (3*sizeof(float)));
    }

    float ptA[3];
    float ptB[3];


    // Populate H with covariance values from point correspondences
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


    // Malloc space for the SVD Output
    w = (float*)malloc(3*sizeof(float));
    memset(w, 0, 3*sizeof(float));

    v = (float**)malloc(3*sizeof(float*));
    for (i=0; i<3; i++) {
        v[i] = (float*)malloc(3*sizeof(float));
        memset(v[i], 0, 3*sizeof(float));
    }

    R = (float**)malloc(3*sizeof(float*));
    for (i=0; i<3; i++) {
        R[i] = (float*)malloc(3*sizeof(float));
        memset(R[i], 0, 3*sizeof(float));
    }

    int ret;
    printf("H Before\n");
    print_mat(H);


    // RUN SVD To get output matrices
    ret = dsvd(H, 3, 3, w, v);

    printf("H After\n");
    print_mat(H);

    printf("v Matrix\n");
    print_mat(v);

    printf("RET: %d\n", ret);

    //TODO: Replace this with this
    //R[0][0] = 1;
    //R[1][1] = 1;
    //R[2][2] = 1;

    // Perform matrix multiplication to get the output matrix R
    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            for (k=0; k < 3; k++) {
                R[i][j] += v[i][k] * H[j][k];
            }
        }
    }

    printf("R Matrix\n");
    print_mat(R);
    return R;
}

// write a ply file given the size of the scan data and the set of points
void write_ply(float** R, pt* pts, char* outfile_name, int n, int m) {
    FILE* out;
    out = fopen(outfile_name, "w");
    
    float p_in[3];
    float p_out[3];

    // Loop over all vertices and print them out
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

    //Print out each point's validity and index
    for (i=0; i < n*m; i++) {
        if (pts[i].valid) {
            fprintf(out, "1 %d\n", ind);
            ind++;
        } else {
            fprintf(out, "0\n");
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Kernel Functions ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// CUDA kernel function call to generate histograms for a given scan
__global__ void generate_histograms(pt* points_list, pt* normals, float* histograms, int* point_pair_tally,
                                    int n, int m, int k, int l) {

    int ind;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    //int offset = idx + idy*m;
    
    // 
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

    //for (i=0; i < some x row value; i++ )
    //  for (i=0; i < some y column dimension for the block; j++)

    // AT i, j
    memset(distances, 0, 4*sizeof(float));


    // Get the radius of the kernel's edge points
    if ((idx > W) && (idy > W)) {
        distances[0] = Norm(Sub(points_list[(idx*m) + idy], points_list[(idx-W)*m + (idy-W)]));
    }
    if ((idx > W) && (idy + W < m)) {
        distances[1] = Norm(Sub(points_list[(idx*m) + idy], points_list[(idx-W)*m + (idy+W)]));
    }
    if ((idx + W < n) && (idy > W)) {
        distances[2] = Norm(Sub(points_list[(idx*m) + idy], points_list[(idx+W)*m + (idy-W)]));
    }
    if ((idx + W < n) && (idy + W < m)) {
        distances[3] = Norm(Sub(points_list[(idx*m) + idy], points_list[(idx+W)*m + (idy+W)])); 
    }
    
    // Compute kernel radius from edge points of kernel
    s[1] = Max(distances, 4);

    /*
      We need neighboring points, so we can index into points list
      we need neighboring normals, so we can index into normals list
    */
    /* for (k = 0; k < K_POINTS - 1; k++) { */
    /*     for (l = k + 1; l < K_POINTS; l++) { */

    // Get point pair and normal pair data
    pt_k = map(k, idx, idy, n, m, points_list);
    pt_l = map(l, idx, idy, n, m, points_list);
    normal_k = map(k, idx, idy, n, m, normals);
    normal_l = map(l, idx, idy, n, m, normals);

    if (pt_k.valid && pt_l.valid) {

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
        f[2] = Dot(u, Sub(pt_t, pt_s)) / f[2];
        f[3] = atan2(Dot(w, normal_t), Dot(u, normal_t));


        ind = hashFeatures(f, s);
        if ((ind < 0) || (ind > 15)) {
            printf("EXITING: IMPROPER HASH\n");
            //exit(1);
        }

        // Update the histograms with the hashed features digit
        histograms[((idx*m + idy) * 16) + ind]++;
        point_pair_tally[(idx*m) + idy]++;
    }

}

// CUDA kernel function call which uses the tally_points array to scale each histogram
// in a given scan
__global__ void scale_hist(float* histograms, int m, int* tally_points) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    int possibles = K_POINTS * (K_POINTS - 1)/2;
    
    int tally_threshold = possibles * HIST_PROPORTION;
        
    int p;

    if (tally_points[i*m + j] < tally_threshold) {
        histograms[(i*m + j)*16] = -1;
    } else {
        for (p=(i*m+j)*16; p < (i*m+j)*16 + 16; p++) {
            //printf("%f ", histograms[p]);
            histograms[p] /= tally_points[i*m + j];
        }
    }
}



// CUDA kernel function call to score given histograms
__global__ void score_hists(float* cur_hists, float* base_hists,
                            candidate* chunk_scores,
                            int n, int m, int n2, int m2,
                            int a_ind, int b_ind) {


    
    // index appropriately into the correct block and thread within the grid
    int i, j, k;

    // Goes from 0 to P-1
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    i = 16 * (a_ind * P + idx);
    
    float score;

    // Initialize all scores to the max
    for(k = idx * CANDIDATES;
        k < idx * CANDIDATES + CANDIDATES;
        k++) {
    
        chunk_scores[k].score = 2;
    }

    // Loop across scan B, get scores
    for (j = 0;
         j < n*m*16;
         j+=16) {

        // compute the score between two histograms
        score = hist_dist(cur_hists, base_hists, i, j);

        consider_candidate(score, i, j,
                           &chunk_scores[idx * CANDIDATES]);
    }
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Helper Functions ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Merge the scores from two lists of size 10 to get the top 10
void merge(candidate* scores, candidate* input) {

    int cur_1 = 0;
    int cur_2 = 0;

    candidate temp;

    while ((cur_1 < CANDIDATES) && (cur_2 < CANDIDATES)) {
        if (input[cur_2].score > scores[cur_1].score) {
            temp = scores[cur_1];
            scores[cur_1] = input[cur_2];
            input[cur_2] = temp;
            cur_1++;

        } else {
            cur_2++;
        }
    }
}

// Recursive sort on the CPU
void Crec_sort(candidate* scores, int len, int count) {

    if (count == CANDIDATES) {
        return;
    }

    int best, i;
    best = 0;

    for (i=1; i<len; i++) {
        if (scores[best].score < scores[i].score) {
            best = i;
        }
    }

    candidate temp;

    temp = scores[0];
    scores[0] = scores[best];
    scores[best] = temp;

    Crec_sort(&(scores[1]), len - 1, count + 1);
}

// Recursive sort on the GPU
__global__ void rec_sort(candidate* scores, int len, int count) {

    if (count == CANDIDATES) {
        return;
    }

    int best, i;
    best = 0;

    // Find best score
    for (i=1; i<len; i++) {
        if (scores[best].score < scores[i].score) {
            best = i;
        }
    }
    
    

    candidate temp;

    // Swap up the best score to the top
    temp = scores[0];
    scores[0] = scores[best];
    scores[best] = temp;
    
    Drec_sort(&(scores[1]), len - 1, count + 1);
}

// Device function for recursive sort
__device__ void Drec_sort(candidate* scores, int len, int count) {

    if (count == CANDIDATES) {
        return;
    }

    int best, i;
    best = 0;

    for (i=1; i<len; i++) {
        if (scores[best].score < scores[i].score) {
            best = i;
        }
    }

    candidate temp;

    temp = scores[0];
    scores[0] = scores[best];
    scores[best] = temp;
    
    

    Drec_sort(&(scores[1]), len - 1, count + 1);
}

// Print a matrix
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

// Score two histogram distributions
__device__ float hist_dist(float* hA, float* hB, int i, int j) {
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

// Dot product
float CDot(pt A, pt B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

// Norm of a vector	
float CNorm(pt A) {
    return sqrt((A.x * A.x) + (A.y * A.y) + (A.z * A.z));
}

// Difference between vectors
pt CSub(pt A, pt B) {

    B.x = A.x - B.x;
    B.y = A.y - B.y;
    B.z = A.z - B.z;

    return B;
}

// Map a kernel index to a scan point
pt Cmap(int ind, int i, int j, int n, int m, pt* mat) {

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

    return mat[a*n + b];
}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// Device Functions ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Maps index to its corresponding vertex/normal
__device__ pt map(int ind, int i, int j, int n, int m, pt* mat) {

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

    return mat[a*n + b];
}

// Dot product
__device__ float Dot(pt A, pt B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}


// Get a vector cross product
__device__ pt Cross(pt A, pt B) {
    pt C;
    C.valid = 1;
    C.x = A.y * B.z - A.z * B.y;
    C.y = A.z * B.x - A.x * B.z;
    C.z = A.x * B.y - A.y * B.x;

    return C;
}

// Get a difference between vectors
__device__ pt Sub(pt A, pt B) {

    B.x = A.x - B.x;
    B.y = A.y - B.y;
    B.z = A.z - B.z;

    return B;
}

// Get a vector length
__device__ float Norm(pt A) {
    return sqrt((A.x * A.x) + (A.y * A.y) + (A.z * A.z));
}

// Get the max of a list
__device__ float Max(float* list, int len) {
    
    int i, max = 0;
    for (i = 1; i < len; ++i){
        if (list[i] > list[max]){
            max = i;            
        }
    }
    return list[max];
}

// Get the min of a list
__device__ float Min(float* list, int len) {

      int i, min = 0;
          for (i = 1; i < len; ++i){
              if (list[i] < list[min]){
                    min = i;
                                          }
                        }
              return min;
}

// Hash features into bits using threshold vector s
__device__ int hashFeatures(float* f, float* s) {

    int i, idx = 0;
    for (i = 0; i < 4; i++){
        if (!(f[i] < s[i])){
            idx += pow(2.0,i);
        }
    }

    return idx;
}

// Add a candidate to our scores array if its score is sufficiently high
__device__ void consider_candidate(float score, int i, int j, candidate* scores){
    int cur_ind;
    candidate temp;
    candidate cur;
    cur.score = score;
    cur.i = i;
    cur.j = j;
    if(cur.score > scores[CANDIDATES-1].score){
        cur_ind = CANDIDATES-1;
        scores[cur_ind] = cur;
        while(scores[cur_ind].score > scores[cur_ind -1].score){      
            temp = scores[cur_ind];
            scores[cur_ind] = scores[cur_ind-1];
            scores[cur_ind-1] = temp;
            if(cur_ind == 1){
                break;
            }
            
        }
    }
}

