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

IsGraph = false; % false for 3D triangle mesh. True for an arbitrary graph (e.g. create a graph from a 3D point by connecting near neighbors).
%% load a triangular mesh

if ~IsGraph
    fprintf('Creating shape...\n');
    load 'david0.mat';
    
    nv = length(surface.X);

    figure;
    trisurf(surface.TRIV, surface.X, surface.Y, surface.Z, zeros(nv,1)); axis equal;axis off; 
    shading interp;lighting phong;cameratoolbar;camlight headlight
end
%% instead - create a graph
if IsGraph
    fprintf('Creating graph...\n');
    load 'david0.mat';
    
    nv = length(surface.X);
    G = sparse(nv, nv);
    TRIV = surface.TRIV;
    V = [surface.X surface.Y surface.Z];
    for i=1:size(TRIV,1)        
        dd = squareform(pdist(V(TRIV(i,:),:)));    
        G(TRIV(i,:), TRIV(i,:)) = dd;
    end
end

%% farthest point sampling
fprintf('Fartherst point sampling...\n');
N = 20; % number of samples. more samples should increase accuracy, but increase complexities (quasy linearly). 
% recomended 10-200 according to application (See [1] for accuracy and time analysis).

tic
if ~IsGraph
    [~, first_idx] = FPS(surface, 1);
    [D_ext, sample2] = FPS(surface, N, first_idx);
else
    [~, first_idx] = FPS_graph(G, 1);
    [D_ext, sample2] = FPS_graph(G, N, first_idx);
end
t = toc;
fprintf('Finished in %f seconds. \n', t);
R = D_ext';

%% Canonical form
fprintf('Canonical form computation...\n');

dim = 3; % embedding dimension 
tic
Z = NMDS(R.^2, sample2, dim);
t = toc;
fprintf('Finished in %f seconds. \n', t);
% Display

figure;
trisurf(surface.TRIV, Z(:,1), Z(:,2), Z(:,3), zeros(nv,1));
axis equal;axis off;
shading interp;lighting phong;cameratoolbar;camlight headlight
title('Canonical form', 'fontsize', 20);