function [A,A_i] = abs_spec_gauss(w,T,energy,n1,n2,v1,v2,dim1,dim2,dim,tpcheck,N,v_max,v_ground,coeff,FC,gaussian_par,kb,w_0,threshold,Z)
%Function to calculate Absorption spectrum for T=0 or finite temperature, using
%Gaussian broadening

A = zeros(1,length(w));
A_i = zeros(1,length(w));
%Zero temperature
if T == 0
    for alpha=1:dim
        A = A + linestrength(0,alpha,coeff,N,v_max,FC,n1,n2,v1,v2,dim1,dim2,tpcheck) * exp(-(w-energy(alpha)).^2/gaussian_par^2);
    end
    A_i = A;
else
    %Finite temperature
    A_i = zeros(v_ground+1,length(w));
    alphin = find(energy<energy(1)+8*w_0)';
    for i = 0:v_ground
        if i*w_0 > threshold
	  break
	end
        for alpha = alphin %1:dim
            A_i(i+1,:) = A_i(i+1,:) + linestrength(i,alpha,coeff,N,v_max,FC,n1,n2,v1,v2,dim1,dim2,tpcheck) * exp(-(w-energy(alpha)+i*w_0).^2/gaussian_par^2);
        end
        A = A + A_i(i+1,:) * exp((-i*w_0)/(kb*T));
    end
    A = A/Z;
end
