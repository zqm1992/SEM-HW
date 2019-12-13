
clear;clc;

Nex = 10;
Ney = 10;
Ex  = 3;
Ey  = 3;

nx1 = Nex;
ny1 = Ney;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

Etol = Ex*Ey;

Q=cal_Q(Ex,Ey,Nex,Ney);
R=cal_R(Ex,Ey,Nex,Ney);

f=ones(size(R,1),1);

% assume still live in [-1,1]x[-1,1]
Ele_lenx = 2/Ex;
Ele_leny = 2/Ey;

dxds = Ele_lenx/2;

[zl,wl]=zwgll(Nex);
Dh=deriv_mat(zl);
Bh=diag(wl);

[xc,yc]=cal_pos(Ex,Ey,Nex,Ney);

Ax=Dh'*Bh*Dh;
Ay=Dh'*Bh*Dh;

Bx=Bh;
By=Bh;

lap=kron(By,Ax)+kron(Ay,Bx);
lap1=lap/dxds;
lap2=kron(diag(ones(Etol,1)),lap1);
lap3=R*Q.'*lap2*Q*R.';

mass =kron(Bh,Bh);
mass1=mass*dxds;
mass2=kron(diag(ones(Etol,1)),mass1);
rhs=R*Q.'*mass2*Q*R.'*f;

numx=Ex*Nex+1;
numy=Ey*Ney+1;

x=lap3\rhs;

% dxmat_loc = kron(eye(size(Dh,1)),Dh);
% dxc = Q.'*kron(diag(ones(Etol,1)),dxmat_loc)*Q*R.'*x;
% x1=reshape(dxc,numx,numy);
% mesh(xc,yc,x1);

% x=R.'*(lap3\rhs);
% x1=reshape(x,numx,numy);
% mesh(xc,yc,x1)