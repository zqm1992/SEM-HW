
clear;clc;
ifkov = 1;
post  = 0;

fdm_use = 1;

% Node number in one element
nx1 = 8;
ny1 = 8;
% Element number on each direction
Ex  = 5;
Ey  = 5;
% total number of elements
Etol = Ex*Ey;
% Number of nodes for multielement on one direction
Nnumx = (nx1-1)*Ex+1;
Nnumy = (ny1-1)*Ey+1;
Nnum  = Nnumx*Nnumy;

slv=1; % 0: CG /todo, 1: FDM

ifvel  = 1;    % evolve velocity field
ifpres = 1;    % project velocity field onto a div-free subspace
ifps   = 0;    % evolve passive scalar per advection diffusion
%------------

nx2 = nx1 - 2;
ny2 = ny1 - 2;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);
nxp = 10*nx1;
nyp = 10*ny1;

[zrm1,wrm1] = zwgll(nx1-1); [zsm1,wsm1] = zwgll(ny1-1); % vel, scalar
[zrm2,wrm2] = zwgll(nx2-1); [zsm2,wsm2] = zwgll(ny2-1); % pres
[zrmd,wrmd] = zwgll(nxd-1); [zsmd,wsmd] = zwgll(nyd-1); % dealias
[zrmp,~   ] = zwuni(nxd-1); [zsmp,~   ] = zwgll(nyp-1); % plt

Drm1 = dhat(zrm1); Dsm1 = dhat(zsm1);
Drm2 = dhat(zrm2); Dsm2 = dhat(zsm2);
Drmd = dhat(zrmd); Dsmd = dhat(zsmd);

Irm1 = eye(nx1); Ism1 = eye(ny1);
Irm2 = eye(nx2); Ism2 = eye(ny2);
Irmd = eye(nxd); Ismd = eye(nyd);

% dx,dy operator
dxk = sparse(kron(Ism1,Drm1));
dyk = sparse(kron(Dsm1,Irm1));
dxkt= sparse(kron(Ism1,Drm1'));
dykt= sparse(kron(Dsm1',Irm1));
DxL = kron(diag(ones(Etol,1)),dxk);
DyL = kron(diag(ones(Etol,1)),dyk);
DxLt= kron(diag(ones(Etol,1)),dxkt);
DyLt= kron(diag(ones(Etol,1)),dykt);

a = -0.5; lx = 2.5; ly=2.0;
cof_lx = (lx/2.0)/Ex;
cof_ly = (ly/2.0)/Ey;

% normal node to dealiasing node
Jr_v2d = interp_mat(zrmd*cof_lx,zrm1*cof_lx);
Js_v2d = interp_mat(zsmd*cof_ly,zsm1*cof_ly);
Jrs_v2d  = sparse(kron(diag(ones(Etol,1)),kron(Js_v2d,Jr_v2d)));

% velocity node to pressure node
Jr_v2p = interp_mat(zrm2*cof_lx,zrm1*cof_lx);
Js_v2p = interp_mat(zsm2*cof_ly,zsm1*cof_ly);
Jrs_v2p  = sparse(kron(diag(ones(Etol,1)),kron(Js_v2p,Jr_v2p)));

% velocity node to pressure node
Jr_p2v = interp_mat(zrm1*cof_lx,zrm2*cof_lx);
Js_p2v = interp_mat(zsm1*cof_ly,zsm2*cof_ly);
Jrs_p2v  = sparse(kron(diag(ones(Etol,1)),kron(Js_p2v,Jr_p2v)));

%-------------------------------------------------------------------------------
% geometry
% [-1,1]x[-1,1]
[xm1,ym1]=cal_pos(Ex,Ey,nx1-1,ny1-1);
[xm2,ym2]=cal_pos(Ex,Ey,nx2-1,ny2-1);
[xmd,ymd]=cal_pos(Ex,Ey,nxd-1,nyd-1);

xm1 = a + lx/2 * (xm1+1)  ; ym1 = a + ly/2 * (ym1+1) ;
xm2 = a + lx/2 * (xm2+1)  ; ym2 = a + ly/2 * (ym2+1) ;
xmd = a + lx/2 * (xmd+1)  ; ymd = a + ly/2 * (ymd+1) ;

if(ifkov)
%-------------------------------------------------------------------------------
% kovazny

casename = 'Kovazny Flow';
cname = 'kov';

% viscosity (velocity, passive scalar)
Re = 40;
visc0 = 1/Re;
visc1 = 0e-0;

% initial condition
vx = 0*xm1;
vy = 0*xm1;
ps = 0*xm1;
pr = 0*xm2;

% exact solution
[vxe,vye] = kov_ex(xm1,ym1,Re);

% forcing
fvx = 0*xm1;
fvy = 0*xm1;
%fps = 0*xm1;

% BC
vxb = vxe;
vyb = vye;
%psb = ps;

% Restrictions
[R, Ry, Rx] =cal_R(Ex,Ey,nx1-1,ny1-1);
[Rp,Ryp,Rxp]=cal_R(Ex,Ey,nx2-1,ny2-1);
Rxvx = Rx; % vx             % dir-dir
Ryvx = Ry;                  % dir-dir
Rxvy = Rx; % vy             % dir-dir
Ryvy = Ry;                  % dir-dir

% Q matrix for scatter gather 
Qmv = cal_Q(Ex,Ey,nx1-1,ny1-1);
Qmp = cal_Q(Ex,Ey,nx2-1,ny2-1);
Qdv = cal_Q(Ex,Ey,nxd-1,nyd-1);
Qmv1d = semq(Ex,nx1-1);
Qmp1d = semq(Ex,nx2-1);

ifxperiodic = 0;
ifyperiodic = 0;

% T=0 ==> steady
T   = 20.0;
CFL = 0.2;

end

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
%mskps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

% time stepper
dx = min(min(diff(xm1)));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0)      % steady
	nt=1;
	dt=0;
else          % movie
	mov=[];
end

% jacobian (2d area ratio)
Jm1 = cof_lx*cof_ly;
Jm2 = cof_lx*cof_ly;
Jmd = cof_lx*cof_ly;
Jmx = cof_lx;
Jmy = cof_ly;

%[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
%[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Irm2,Ism2,Drm2,Dsm2);
%[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% Mass matrix
Bmh = diag(wrm1);
Bdh = diag(wrmd);
Bmk = sparse(kron(Bmh,Bmh));
Bdk = sparse(kron(Bdh,Bdh));
Bmu = sparse(kron(diag(ones(Etol,1)),Bmk));
Bdu = sparse(kron(diag(ones(Etol,1)),Bdk));

% laplacian operator preparation (on vel nodes)
dx2  = kron(Bmh,Drm1.'*Bmh*Drm1)/(Jmx*Jmx);
dy2  = kron(Dsm1.'*Bmh*Dsm1,Bmh)/(Jmy*Jmy);
lapL = kron(diag(ones(Etol,1)),(dx2+dy2)*Jm1);
lapu = Qmv.'*lapL*Qmv;

if(fdm_use == 1)
    % fast diagonalization setup
    Lx = max(max(xm1))-min(min(xm1))/Ex;
    Ly = max(max(ym1))-min(min(ym1))/Ey;
    
    % Velocity
    Bxv = (Lx/2)*diag(wrm1); Byv = (Ly/2)*diag(wsm1);
    Dxv = (2/Lx)*Drm1;       Dyv = (2/Ly)*Dsm1;
    Axv = Dxv'*Bxv*Dxv;      Ayv = Dyv'*Byv*Dyv;
    
    Bxvx = Rx*Qmv1d.'*kron(diag(ones(Ex,1)),Bxv)*Qmv1d*Rx';
    Byvx = Ry*Qmv1d.'*kron(diag(ones(Ey,1)),Byv)*Qmv1d*Ry';
    Axvx = Rx*Qmv1d.'*kron(diag(ones(Ex,1)),Axv)*Qmv1d*Rx';
    Ayvx = Ry*Qmv1d.'*kron(diag(ones(Ey,1)),Ayv)*Qmv1d*Ry';
    
    Bxvy = Rx*Qmv1d.'*kron(diag(ones(Ex,1)),Bxv)*Qmv1d*Rx';
    Byvy = Ry*Qmv1d.'*kron(diag(ones(Ey,1)),Byv)*Qmv1d*Ry';
    Axvy = Rx*Qmv1d.'*kron(diag(ones(Ex,1)),Axv)*Qmv1d*Rx';
    Ayvy = Ry*Qmv1d.'*kron(diag(ones(Ey,1)),Ayv)*Qmv1d*Ry';
    
    [Sxvx,Lxvx] = eig(Axvx,Bxvx);
    [Syvx,Lyvx] = eig(Ayvx,Byvx);
    Sxvx=Sxvx*diag(1./sqrt(diag(Sxvx'*Bxvx*Sxvx)));
    Syvx=Syvx*diag(1./sqrt(diag(Syvx'*Byvx*Syvx)));
    Lvx = visc0 * (bsxfun(@plus,diag(Lxvx),diag(Lyvx)'));
    
    [Sxvy,Lxvy] = eig(Axvy,Bxvy);
    [Syvy,Lyvy] = eig(Ayvy,Byvy);
    Sxvy=Sxvy*diag(1./sqrt(diag(Sxvy'*Bxvy*Sxvy)));
    Syvy=Syvy*diag(1./sqrt(diag(Syvy'*Byvy*Syvy)));
    Lvy = visc0 * (bsxfun(@plus,diag(Lxvy),diag(Lyvy)'));
    
    % Pressure
    Myvx = Ryvx'*Ryvx; Mxvx = Rxvx'*Rxvx;
    Myvy = Ryvy'*Ryvy; Mxvy = Rxvy'*Rxvy;
    
    Bxiv = diag(1./diag(Bxv));
    Byiv = diag(1./diag(Byv));
    
    Byp = kron(eye(Ey),Byiv);
    Axp = kron(eye(Ey),Bxiv);
    Ayp = kron(eye(Ey),Byiv);
    Bxp = kron(eye(Ey),Bxiv);
    
    % apply bc for Biv
    Byp(1,1)                    =0;
    Byp(size(Byp,1),size(Byp,2))=0;
    Axp(1,1)                    =0;
    Axp(size(Axp,1),size(Axp,2))=0;
    Ayp(1,1)                    =0;
    Ayp(size(Ayp,1),size(Ayp,2))=0;
    Bxp(1,1)                    =0;
    Bxp(size(Bxp,1),size(Bxp,2))=0;
    
    Jr21= kron(eye(Ex),Jr_p2v);
    Js21= kron(eye(Ey),Js_p2v);
    
    Byvblock = kron(eye(Ey),Byv);
    Bxvblock = kron(eye(Ex),Bxv);
    Dyvblock = kron(eye(Ey),Dyv);
    Dxvblock = kron(eye(Ex),Dxv);
    
    Byp = Js21'*Byvblock*(         (Byp)          )*Byvblock*Js21; % attack vx
    Axp = Jr21'*Bxvblock*(Dxvblock*(Axp)*Dxvblock')*Bxvblock*Jr21;
    Ayp = Js21'*Byvblock*(Dyvblock*(Ayp)*Dyvblock')*Byvblock*Js21; % attack vy
    Bxp = Jr21'*Bxvblock*(         (Bxp)          )*Bxvblock*Jr21;
    
    Byp = Qmp1d.'*Byp*Qmp1d;
    Axp = Qmp1d.'*Axp*Qmp1d;
    Ayp = Qmp1d.'*Ayp*Qmp1d;
    Bxp = Qmp1d.'*Bxp*Qmp1d;
    
    [Sxpr,Lxpr] = eig(Axp,Bxp);
    [Sypr,Lypr] = eig(Ayp,Byp);
    Sxpr=Sxpr*diag(1./sqrt(diag(Sxpr'*Bxp*Sxpr)));
    Sypr=Sypr*diag(1./sqrt(diag(Sypr'*Byp*Sypr)));
    Lpr  = bsxfun(@plus,diag(Lxpr),diag(Lypr)');
    Lipr = 1 ./ Lpr;
    Lipr(find(Lipr>1e10))  = 0;
    Lipr(find(Lipr<-1e10)) = 0;
    
    Jm_fdm = (Lx/2)*(Ly/2);
    Bm_fdm = Jm_fdm*(wrm1*wsm1');
    Bm_fdm = reshape(Qmv.'*repmat(Bm_fdm(:),Etol,1),Nnumx,Nnumy);
    Bmi_fdm= 1.0./Bm_fdm;
    
end

% time advance
time = 0; 

% initialize histories
time1 = time*0; time2 = 0; time3=0;
vx1 = vx*0; vx2 = vx1; vx3 = vx2; gvx1 = vx1; gvx2 = vx1; gvx3 = vx1;
vy1 = vy*0; vy2 = vy1; vy3 = vy2; gvy1 = vy1; gvy2 = vy1; gvy3 = vy1;
ps1 = ps*0; ps2 = ps1; ps3 = ps2; gps1 = ps1; gps2 = ps1; gps3 = ps1;
pr1 = pr*0;

%vt1 = xm1.*xm1+ym1.*ym1;

for it=1:nt

	% update time
	time3=time2; time2=time1; time1=time;
	time = time1 + dt;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
        % FDM
        Livx = 1 ./ (b(1) + Lvx);
		Livy = 1 ./ (b(1) + Lvy);
    end

    % update histories
	vx3=vx2; vx2=vx1; vx1=vx; gvx3=gvx2; gvx2=gvx1;
	vy3=vy2; vy2=vy1; vy1=vy; gvy3=gvy2; gvy2=gvy1;
	ps3=ps2; ps2=ps1; ps1=ps; gps3=gps2; gps2=gps1;
				      pr1=pr;
                      
    if(ifvel)
        gvx1 = mass(fvx,Bmu,Qmv,Jm1) - advect(vx1,vx1,vy1,Bdu,Jrs_v2d,DxL,DyL,Qmv,Jm1,Jmx,Jmy);
        gvy1 = mass(fvy,Bmu,Qmv,Jm1) - advect(vy1,vx1,vy1,Bdu,Jrs_v2d,DxL,DyL,Qmv,Jm1,Jmx,Jmy);
        
        % pressure forcing
		[px,py]=vgradp(pr1,Bmu,Jrs_p2v,DxLt,DyLt,Qmv,Qmp,Jm1,Jmx,Jmy);
        
        bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Bmu,Qmv,Jm1);
        umass = mass(vxb,Bmu,Qmv,Jm1);
		bvx = bvx - hlmhltz(vxb,visc0,b(1),lapu,umass);
		bvx = bvx + px;
        bvx = mskvx.*bvx;
		bvx = ABu(Ryvx,Rxvx,bvx);
        
        bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Bmu,Qmv,Jm1);
        umass = mass(vyb,Bmu,Qmv,Jm1);
		bvy = bvy - hlmhltz(vyb,visc0,b(1),lapu,umass);
		bvy = bvy + py;
        bvy = mskvy.*bvy;
		bvy = ABu(Ryvy,Rxvy,bvy);
        
        if(fdm_use)
            vyh = fdm(bvy,Sxvy,Syvy,Livy);
            vxh = fdm(bvx,Sxvx,Syvx,Livx);
        else
            vyh = visc_slv(bvy,lapu,Bmu,Qmv,Jm1,b(1),mskvy,visc0);
            vxh = visc_slv(bvx,lapu,Bmu,Qmv,Jm1,b(1),mskvx,visc0);
        end
        
        vx  = ABu(Ryvx',Rxvx',vxh) + vxb;
		vy  = ABu(Ryvy',Rxvy',vyh) + vyb;

		% pressure projection
		if(ifpres)
            [vx,vy,pr] = pres_proj(vx,vy,pr1,b(1),Bmi_fdm,Rxvx,Ryvx,Rxvy,Ryvy,...
                                Bmu,Jrs_p2v,DxL,DyL,DxLt,DyLt,...
                                Jm1,Jmx,Jmy,Qmv,Qmp,Sxpr,Sypr,Lipr);
        end
        if(mod(it,50)==0)
            mesh(xm1,ym1,vx);
        end

    end
    
end
