if broadtype == 0
    [S, S_alpha] = pl_spec_gauss(w,T,energy,n1,n2,v1,v2,dim1,dim2,tpcheck,N,v_max,v_ground,coeff,FC,broad_par,kb,w_0,threshold,Z);
    [A, A_i]     = abs_spec_gauss(w,T,energy,n1,n2,v1,v2,dim1,dim2,dim,tpcheck,N,v_max,v_ground,coeff,FC,broad_par,kb,w_0,threshold,Z);
else
    [S, S_alpha] = pl_spec_lorentz(w,T,energy,n1,n2,v1,v2,dim1,dim2,tpcheck,N,v_max,v_ground,coeff,FC,broad_par,kb,w_0,threshold,Z);
    [A, A_i]     = abs_spec_lorentz(w,T,energy,n1,n2,v1,v2,dim1,dim2,dim,tpcheck,N,v_max,v_ground,coeff,FC,broad_par,kb,w_0,threshold,Z);
end

