%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [Cu] = advect(u,cx,cy,BL_d,JL_rs,DL_r,DL_s,Q,Jm1,Jx,Jy);

nx = size(cx,1);
ny = size(cx,2);
num= nx*ny;

uvec =reshape(u,num,1);
cxvec=reshape(cx,num,1);
cyvec=reshape(cy,num,1);

CxL_vec = JL_rs*Q*cxvec;
CyL_vec = JL_rs*Q*cyvec;

% dudx on dealising points
dudx = JL_rs*DL_r*Q*uvec/Jx;

% dudy on dealising points
dudy = JL_rs*DL_s*Q*uvec/Jy;

Cu   = Q.'*(JL_rs.'*(BL_d*(CxL_vec.*dudx + CyL_vec.*dudy)))*Jm1;
Cu   = reshape(Cu,nx,ny);
end
