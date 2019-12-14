%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [Cu] = advect(u,cx,cy,BL_d,JL_rs,DL_r,DL_s,Q);

nx = size(cx,1);
ny = size(cx,2);
num= nx*ny;

uvec =reshape(u,num,1);
cxvec=reshape(cx,num,1);
cyvec=reshape(cy,num,1);

CxL_vec = JL_rs*Q*cxvec;
CyL_vec = JL_rs*Q*cyvec;

% dudx on dealising points
dudx = JL_rs*DL_r*Q*uvec;

% dudy on dealising points
dudy = JL_rs*DL_s*Q*uvec;

Cu   = Q.'*(JL_rs.'*(BL_d*(CxL_vec.*dudx + CyL_vec.*dudy)));
Cu   = reshape(Cu,nx,ny);
end
