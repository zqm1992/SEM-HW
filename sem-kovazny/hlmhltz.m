%
%	(v,b0*u -visc*\del^2 u)
%
function [Hu] =  hlmhltz(u,visc,b0,lapu,umass)
nx   = size(u,1);
ny   = size(u,2);
num  = nx*ny;
uvec = reshape(u,num,1);

Hu   = reshape(visc*lapu*uvec,nx,ny);
Hu   = Hu +   b0*umass;
end
