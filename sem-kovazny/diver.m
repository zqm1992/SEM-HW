%
% (q,dudx+dvdy)
% JL_rs is p2v
function [DDu] = diver(ux,uy,BL_m,JL_rs,DL_r,DL_s,Qv,Qp,Jm,Jx,Jy)
nx = size(ux,1);
ny = size(ux,2);
num= nx*ny;

uxvec = reshape(ux,num,1);
uyvec = reshape(uy,num,1);

dudx  = (DL_r*(Qv*uxvec))/Jx;
dvdy  = (DL_s*(Qv*uyvec))/Jy;

DDu   = (Qp.'*(JL_rs.'*(BL_m*(dudx+dvdy))))*Jm;

nx    = round(sqrt(size(Qp,2)));
ny    = nx;
DDu   = reshape(DDu,nx,ny);
end

