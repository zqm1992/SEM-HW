function [xc,yc]=cal_pos(Ex,Ey,Nx,Ny)
[zl,wl]=zwgll(Nx);

dx=2/Ex;
dy=2/Ey;
Nt=(Ex*Nx+1);

zs = (zl+1)/2*dx;
xc = zeros(Nt,Nt);
yc = zeros(Nt,Nt);

for ey=1:Ey
    for ex=1:Ex
        for iy=1:Ny+1
            for ix=1:Nx+1
                nodex = (ex-1)*Nx + ix;
                nodey = (ey-1)*Ny + iy;
                x0 = -1.0+(ex-1)*dx;
                y0 = -1.0+(ey-1)*dy;
                x1 = x0+zs(ix);
                y1 = y0+zs(iy);
                xc(nodex,nodey)=x1;
                yc(nodex,nodey)=y1;
            end
        end
    end
end

end