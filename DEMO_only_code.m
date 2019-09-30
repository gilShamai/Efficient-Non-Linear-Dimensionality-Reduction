%{
    In this demo:
    
    Efficient approximation of non-linear dimensionality reduction by
    classical scaling (also known as multidimensional scaling (MDS) and
    Isomap).
    
    The approximation is deviated by less than ~0.001 from the full
    canonical form computation (see paper).
    
    Total computation time of a 50K vertices mesh with this approximation: 
    2-3 seconds (the full computation is impractical with this amount of
    points, see paper).
    
    If using these ideas please cite:

    [1] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. "Efficient 
    Inter-Geodesic Distance Computation and Fast Classical Scaling". 
    IEEE transactions on pattern analysis and machine intelligence (2018).
    
    [2] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. 
    "Accelerating the computation of canonical forms for 3D nonrigid 
    objects using multidimensional scaling." In Proceedings of the 
    2015 Eurographics Workshop on 3D Object Retrieval, pp. 71-78. 
    Eurographics Association, 2015.
%}

close all;
clear all;
clc;
addpath('fastmarch');

N = 20; % number of samples
dim = 3; % embedding dimension 
%% preprocessing
load 'david0.mat';
nv = length(surface.X);
[~, first_idx] = FPS(surface, 1);
[D_ext, sample2] = FPS(surface, N, first_idx);
R = D_ext';
Z = NMDS(R.^2, sample2, dim);
