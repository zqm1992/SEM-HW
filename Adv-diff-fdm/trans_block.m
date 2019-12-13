function [y]=trans_block(vin,Ex,Ey,Nx,Ny)
Nl = (Nx+1)*(Ny+1);
num= length(vin);
y  = speye(num);
for ii=1:Ex*Ey
    x0=(ii-1)*Nl+1;
    x1=ii*Nl;
    y(x0:x1,x0:x1)=diag(vin(x0:x1));
end

end