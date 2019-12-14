%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [Bu] = mass(u,Bm,Q,Jac);

nx  = size(u,1);
ny  = size(u,2);
num = nx*ny;

uvec= reshape(u,num,1);
Bu  = Jac*Q'*Bm*Q*uvec;
Bu  = reshape(Bu,nx,ny);

end
