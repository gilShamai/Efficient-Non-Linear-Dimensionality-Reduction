%{ 
    NMDS (Nystrome MDS) - Fast multidimensional scaling solver

    Description:

     An Accelerated Multidimensional scaling solver. 

    Input: 

     sample2 - N x 1 - indices of samples.
     R - nv x N geodesic distances from samples to the rest.
     dim - The embedding dimension.

    Output: 
    
     Z       - the 3D canonical form.
    
    Reference:

    [1] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. "Efficient 
    Inter-Geodesic Distance Computation and Fast Classical Scaling". 
    IEEE transactions on pattern analysis and machine intelligence (2018).
    
    [2] Gil Shamai, Michael Zibulevsky, and Ron Kimmel. 
    "Accelerating the computation of canonical forms for 3D nonrigid 
    objects using multidimensional scaling." In Proceedings of the 
    2015 Eurographics Workshop on 3D Object Retrieval, pp. 71-78. 
    Eurographics Association, 2015.

    FOR ACADEMIC USE ONLY.
    ANY ACADEMIC USE OF THIS CODE MUST CITE THE ABOVE REFERENCE. 
    FOR ANY OTHER USE PLEASE CONTACT THE AUTHORS.
%}

function Z = NMDS(R, sample2, dim)

N = length(sample2);
N1 = round(0.5*N);

% Compute the approximation components E = R*V(:,1:N1)*D_*V(:,1:N1)'*R'
% (Subsection 3.1.2, Equation 12)

U = R(sample2, :);
[V,D] = eigs(0.5*(U + U'), N1, 'lm');
d = diag(D);
d_ = 1./d(1:N1);
D_ = diag(d_);

% compute the canonical form (Subsection 3.3)
Z = mds_get_X_nys(V(:,1:N1), D_, R, dim); 
