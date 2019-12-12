function [uslv] = visc_slv_sem(rslv,A,B,R,b,nu,nx1,ny1)
	
		rslv = reshape(rslv,[nx1*ny1,1]);

		S    = R *(b(1)*B + nu*A)*R';
		rslv = R *rslv;
		uslv = S \rslv;
		uslv = R'*uslv;
		rslv = R'*rslv;
		
		uslv = reshape(uslv,[nx1,ny1]);
		rslv = reshape(rslv,[nx1,ny1]);

end