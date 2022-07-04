function [S,S_alpha] = pl_spec_lorentz(w,T,energy,n1,n2,v1,v2,dim1,dim2,tpcheck,N,v_max,v_ground,coeff,FC,homwidth,kb,w_0,threshold,Z)
%Function to calculate PL spectrum for T=0 or finite temperature, using
%Lorentzian broadening
dim = dim1 + dim2;

S=zeros(1,length(w));
S_alpha=zeros(1,length(w)); %Component spectrum from alpha eigenstate
lowen= energy(1);

%Zero temperature spectrum
if T == 0
    loweigs = find(energy == lowen);
    for loweig = loweigs
        for i=0:v_ground
            S = S + linestrength(i,loweig,coeff,N,v_max,FC,n1,n2,v1,v2,dim1,dim2,tpcheck)* homwidth^2 ./ ((w-energy(loweig)+i*w_0).^2 + homwidth^2);%exp(-(w-energy(loweig)+i*w_0).^2/gaussian_par^2);
        end
    end
    S_alpha = S;
%Finite temperature spectrum
else
    S_alpha=zeros(dim,length(w)); %Component spectrum from alpha eigenstate

    for alpha=1:dim
        if energy(alpha)-lowen > threshold
            break
        end
        for i=0:v_ground
            S_alpha(alpha,:) = S_alpha(alpha,:) + linestrength(i,alpha,coeff,N,v_max,FC,n1,n2,v1,v2,dim1,dim2,tpcheck)* homwidth^2 ./ ((w-energy(alpha)+i*w_0).^2 + homwidth^2); % exp(-(w-energy(alpha)+i*w_0).^2/gaussian_par^2);
        end
        S = S + S_alpha(alpha,:) * exp((-energy(alpha)+lowen)/(kb*T));
    end
    S = S/Z;
end
end
