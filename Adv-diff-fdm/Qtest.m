
% E*(N+1) x E*N+1
% E=3;
% N=3;
% Q=semq(E,N);

Ex=3;
Ey=3;
Nx=3;
Ny=3;

E = Ex*Ey;
Nl=(Nx+1)*(Ny+1);
Nt=(Ex*Nx+1)*(Ey*Ny+1);

% E*Nl x Nt
Q=spalloc(E*Nl,Nt,3*E*Nl);

for ey=1:Ey
    for ex=1:Ex
        for iy=1:Ny+1
            for ix=1:Nx+1
                Qi=((ey-1)*Ex+ex-1)*Nl + ((iy-1)*(Ny+1)+ix);
                nodex = (ex-1)*Nx + ix;
                nodey = (ey-1)*Ny + iy;
                Qj= (nodey-1)*(Ex*Nx+1) + nodex;
                Q(Qi,Qj)=1;
            end
        end
    end
end

            
        