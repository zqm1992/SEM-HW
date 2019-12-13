clear;clc;

Ex=2;
Ey=3;
E = Ex*Ey;

Nx=3;
Ny=3;
Nl=(Nx+1)*(Ny+1);

Nxt = Ex*Nx+1;
Nyt = Ey*Ny+1;
Nt  = (Ex*Nx+1)*(Ey*Ny+1);

Rx  = speye(Nxt,Nxt);
Rx  = Rx(2:end-1,:);

Ry  = speye(Nyt,Nyt);
Ry  = Ry(2:end-1,:);

R   = kron(Ry,Rx);