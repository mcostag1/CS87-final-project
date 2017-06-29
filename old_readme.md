# CS87Project
repo for your course project

Project Plan - Mesh Zippering Algorithm
Simon Bloch and Martina Costagliola

Week 1: 3/27 For this week, we aim to compile and run the existing sequential version of the program. There is a complete copy of this sequential version of the zipper program that is discussed in Zippered Polygon Meshes [TL94] from the Stanford - Graphics website (https://graphics.stanford.edu/software/zippack/). By the end of this week, we aim to run this program, and perform part of Experiment 1: Sequential vs. Parallel Mesh generation. We will compute average run times of the sequential mesh generator on various sized range scans and decide whether parallelizing the zippering algorithm is feasible and extensible using CUDA and MPI. 

ANNOTATION:
We spent the majority of week 1 and week 2 attempting to debug the legacy code that was used in the paper Zippered Polygon Meshes. We first believed we did not have the correct compiler to run the program, and took steps to acquire the necessary library to build it. This proved to be a dead-end, and we set on debugging the existing code for a majority of the remaining time. 
Debugging involved: 
-changing “#include” statements to accommodate more recent libraries
-fixing function calls which disagree with the function prototypes
-reconciling type mismatch for code written for a 32-bit machine
    Eventually, we switched gears, instead planning to 

Week 2: 4/3 We aim during this week to have parallelized, robust mesh generation occur simultaneously on all nodes. This is ambitious to accomplish–at the very least we aim to parallelize the creation of a mesh for one range image using CUDA. We will then complete the second half of Experiment 1, and compute the average run time of our parallel mesh generator. 

ANNOTATION:
See notes from week 2 above.

Week 3: 4/10 This week, we will finish simultaneous parallelization of mesh generation using MPI, and implement the major functionality of our reduction step (i.e. passing increasingly precise mesh scans between nodes). We will focus on the MPI portion here, and less on the zippering algorithm. 

ANNOTATION:
Since we were unsuccessful at compiling and running the sequential zipper program, we are going to redirect our focus into implementing and parallelizing the ICP algorithm, and postpone the additional components of the entire zipper program. If there is time at the end of the project, we hope to pursue these additional components, but for now our aim is to parallelize that algorithm and apply it to the existing range scan data we have.
So far we have found and incorporated a long file which runs PCA (principal components analysis) on a matrix. We will eventually parallelize this code, but for now we are integrating it as a crucial component of our ICP algorithm. It runs now on an input file, and the next step will be to parse the PLY file data for passing into the PCA program. 

Week 4: 4/17 Here, we will finish our implementation of the ICP mesh merging algorithm, and begin work on removing overlap and cleaning the mesh at the seams.

Week 5: 4/24 This final week will be dedicated to improving the completeness of our ICP algorithm, running any necessary supplemental tests, and exploring other performance improvements/necessary changes in our parallelized algorithms.



Our code is buildable and runnable. The command line should look as follows: 
    ./p_zipper filename #rows #cols option

Test values for command line arguments:

filename = pca_data.txt
#rows = 36
#columns = 8
option = R    
