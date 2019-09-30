% gets V, D_, C.
% returns X such that X*X' = -0.5JCVD_V'C'J

function X = mds_get_X_nys(V, D_, C, dim)
n = size(C, 1);
B = C*V;
J_B = B-1/n*repmat(sum(B),n,1);

[Q,W] = qr(J_B, 0);

B2 = -0.5*W*D_*W';
[V2, L] = eigs(0.5*(B2+B2'),dim,'la');

X = Q*(V2*sqrt(L));
