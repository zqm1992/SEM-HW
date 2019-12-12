
clear;clc;

% polynomial order in each element
Nex = 6;
Ney = 6;
% Element number on each direction
Ex  = 4;
Ey  = 4;
% total number of elements
Etol = Ex*Ey;
% Number of nodes for multielement on one direction
Nnumx = Nex*Ex+1;
Nnumy = Ney*Ey+1;
Nnum  = Nnumx*Nnumy;

nx1 = Nex;
ny1 = Ney;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

% Element length
Elemx_len = 1/Ex;
Elemy_len = 1/Ey;
Elem_area = Elemx_len*Elemy_len;
Eledx_len = 1/Ex;
Eledy_len = 1/Ey;
Eled_area = Eledx_len*Eledy_len;

[zrm1,wrm1] = zwgll(nx1);
[zsm1,wsm1] = zwgll(ny1);
[zrd1,wrd1] = zwgll(nxd);
[zsd1,wsd1] = zwgll(nyd);

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drmd = dhat(zrd1);
Dsmd = dhat(zsd1);

Irm1 = eye(nx1+1);
Ism1 = eye(ny1+1);
Irmd = eye(nxd+1);
Ismd = eye(nyd+1);

[xm1,ym1]=cal_pos(Ex,Ey,nx1,ny1);
[xmd,ymd]=cal_pos(Ex,Ey,nxd,nyd);

Ndnum = size(xmd,1)*size(xmd,2);

% normal node to dealiasing node
Jr1d = interp_mat(zrd1/Ex,zrm1/Ex);
Js1d = interp_mat(zsd1/Ey,zsm1/Ey);
Jrs  = sparse(kron(diag(ones(Etol,1)),kron(Js1d,Jr1d)));

%=================================================================
% data
nu= 1e-4;
u = (xm1-0.5).^2 + (ym1-0).^2; u = exp(-u/0.03);
cx=-ymd;
cy= xmd;
f = 0*xm1;

% BC
ub = u*0;

% Q matrix for scatter gather 
Qm = cal_Q(Ex,Ey,nx1,ny1);
Qd = cal_Q(Ex,Ey,nxd,nyd);

% dir-dir boundary is used here
[R,Ry,Rx]=cal_R(Ex,Ey,Nex,Ney);

% time (steady state if T=0)
T   = 6*pi;
CFL = 0.2;

% time stepper
dx = min(min(diff(xm1)))/Ex;
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;
t  = 0;
it = 0;

% BDF3-EXT3
a = zeros(3,1);
b = zeros(4,1);

% Mass matrix
Bmh = diag(wrm1);
Bdh = diag(wrd1);
Bmk = sparse(kron(Bmh,Bmh));
Bdk = sparse(kron(Bdh,Bdh));
Bmu = sparse(Elem_area*kron(diag(ones(Etol,1)),Bmk));
Bdu = sparse(Eled_area*kron(diag(ones(Etol,1)),Bdk));
B   = Qm'*Bmu*Qm;

% Construct our laplacian operator
Ax=Drm1'*Bmh*Drm1;
Ay=Dsm1'*Bmh*Dsm1;
% kron 
lapk = sparse(kron(Bmh,Ax)+kron(Ay,Bmh));
% laplacian in physical
lapp = lapk/Elem_area;
% laplacian for multi element
A = Qm.'*kron(diag(ones(Etol,1)),lapp)*Qm;

% Advection operator
Cxd = reshape(cx, [Ndnum,1]);
Cyd = reshape(cy, [Ndnum,1]);

Cxd_block = sparse(diag(Qd*Cxd));
Cyd_block = sparse(diag(Qd*Cyd));

dxk = sparse(kron(Ism1,Drm1));
dyk = sparse(kron(Dsm1,Irm1));

dxp = kron(diag(ones(Etol,1)),dxk/Elem_area);
dyp = kron(diag(ones(Etol,1)),dyk/Elem_area);

% scatter gather need
C =  Qm'*Jrs'*Bdu*(Cxd_block*Jrs*dxp+Cyd_block*Jrs*dyp)*Qm;

t0 = 0;
t1 = 0;
t2 = 0;
u0 = u*0;
u1 = u0;
u2 = u0;
f1 = u0;
f2 = u0;

mesh(xm1,ym1,u);

for it=1:nt

	u = mask(u,Rx,Ry);

	% update histories
	t3=t2; t2=t1; t1 = t;
	u3=u2; u2=u1; u1 = u;
	f3=f2; f2=f1; 
    f1 = reshape(-C*reshape(u,Nnum,1),Nnumx,Nnumy);
    f1 = mask(f1,Rx,Ry);

	t = t + dt;

	if(it<=3)
		[a,b] = bdfext3([t t1 t2 t3]);
	end;
	
	% form BDF rhs
    rn1 = reshape(b(2)*u1+b(3)*u2+b(4)*u3,Nnumx*Nnumy,1);
	rhs = a(1)*f1 +a(2)*f2 +a(3)*f3 - reshape(Qm.'*Bmu*Qm*rn1,Nnumx,Nnumy);
	rhs = mask(rhs,Rx,Ry);

	% viscous solve
	uh = visc_slv_sem(rhs,A,B,R,b,nu,Nnumx,Nnumy);

	u  = uh;
    mesh(xm1,ym1,u);

	% vis
	%if(mod(it,floor(0.05*nt))==0)
    if(mod(it,30)==0)
		mesh(xm1,ym1,u);
		title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);
		pause(0.05)
	end

end



