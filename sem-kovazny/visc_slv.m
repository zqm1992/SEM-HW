%
% viscous solve
%
% ==============Qiming need ============== %
% need to fix the following two function
function [u] = visc_slv(b,lapu,Bm,Q,Jac,b0,mask,visc0)
	
	u = cg_visc(b,mask,0*b,lapu,Bm,Q,Jac,b0,1e-5,1e3,visc0);
	
end

function [y] = hlm_cg_use(u0,mask,b0,visc0,lapu,Bmu,Qmv,Jac)
    umass = mass(u0,Bmu,Qmv,Jac);
    y     = mask.*hlmhltz(u0,visc0,b0,lapu,umass);
end

%-------------------------------------------------------------------------------
% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
function [x,k,rsqnew] = cg_visc(b,mask,x0,lapu,Bmu,Qmv,Jac,b0,tol,maxiter,visc0);
	x = x0;
    
	r = b - hlm_cg_use(x,mask,b0,visc0,lapu,Bmu,Qmv,Jac); % r = b - Ax
	rsqold=dot(r,r);
	
	if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;
	
	p=r;
	for k=1:maxiter
		Ap = hlm_cg_use(p,mask,b0,visc0,lapu,Bmu,Qmv,Jac); % Ap = A*p
		al = rsqold / dot(p,Ap);
		x  = x + al*p;
		r  = r - al*Ap;
		rsqnew=dot(r,r); if(sqrt(rsqnew) < tol); return; end;
		be = rsqnew / rsqold;
		p  = r + be*p;
		rsqold = rsqnew;
	end
	['cg iter:',num2str(k),', residual:',num2str(sqrt(rsqnew))]
end


