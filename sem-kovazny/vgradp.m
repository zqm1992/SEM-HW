% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,BL_m,JL_rs,DL_r_t,DL_s_t,Qv,Qp,Jm,Jx,Jy)

nx = size(p,1);
ny = size(p,2);
num= nx*ny;

numo=size(Qv,2);
nxo =round(sqrt(numo));
nyo =round(sqrt(numo));

pvec = reshape(p,num,1);

JBp  = BL_m*JL_rs*Qp*pvec;

px   = Qv.'*DL_r_t*JBp*Jm/Jx;
py   = Qv.'*DL_s_t*JBp*Jm/Jy;

px   = reshape(px,nxo,nyo);
py   = reshape(py,nxo,nyo);

end


