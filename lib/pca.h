#ifndef PCA_H
#define PCA_H

#include "eig3.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void normalize(float v[3]) {
  float ss = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  float l = sqrt(ss);
  v[0] /= l;
  v[1] /= l;
  v[2] /= l;
}

void cross(float a[3], float b[3], float c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

float dot(float a[3], float b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void printvec(const float* v) {
  printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
}

float randreal() {
    return (rand() / (float)RAND_MAX) * 2 - 1;
}

void PCA(int n, float** A, float mean[3], float V[3][3], float d[3]) {

  mean[0] = mean[1] = mean[2] = 0;


  int i, j, k;
  // compute sum of points
  for (i=0; i<n; ++i) {
    for (j=0; j<3; ++j) {
      mean[j] += A[i][j];
    }
  }

  for (j=0; j<3; ++j) {
    mean[j] /= n;
  }

  // A^T A after mean shift
  float cov[3][3] = {
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 }
  };

  for (i=0; i<n; ++i) {
    for (j=0; j<3; ++j) {
      for (k=0; k<3; ++k) {
        // cov = sum of (point j - mean) * (point k - mean)^T
        cov[j][k] += (A[i][j] - mean[j])*(A[i][k] - mean[k]);
      }
    }
  }


  eigen_decomposition(cov, V, d);
  
}

#endif
