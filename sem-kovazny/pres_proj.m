%
% pressure project
%
function [vx,vy,pr] = pres_proj(ux,uy,pr1,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,...
                                BL_m,Jrs_p2v,DL_r,DL_s,DL_r_t,DL_s_t,...
                                Jm,Jx,Jy,Qmv,Qmp,Sxp,Syp,Lip)

	g = -diver(ux,uy,BL_m,Jrs_p2v,DL_r,DL_s,Qmv,Qmp);
    
    delp = b0 * fdm(g,Sxp,Syp,Lip);

	[px,py] = vgradp(delp,BL_m,Jrs_p2v,DL_r_t,DL_s_t,Qmv,Qmp,Jm,Jx,Jy);

	dpvx = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	dpvy = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	vx = ux + dpvx;
	vy = uy + dpvy;

	pr = pr1 + delp;

end
