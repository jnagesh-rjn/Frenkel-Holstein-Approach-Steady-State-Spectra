function C_T = coherencefT(coeff,energy,n1,v1,FC,kb,T,N,v_max,w_0,dim,PBC,threshold)
%Function to calculate and return the thermal average of the Coherence
%function i.e., <C(s)>_T
%The function also plots a heatmap of the coherence function values for the
%20 lowest eigenstates where each row represents an eigenstate and the
%column represents the s value

if PBC == 0 %Open boundary conditions
    ssize = 2*N-1;
    srange = [-(N-1):(N-1)];
else        %Periodic Boundary conditions
    ssize = N;
    if mod(N,2) == 0
        srange = [-N/2+1:N/2];
    else
        srange = [-(N-1)/2:(N-1)/2];
    end
end

C_T = zeros(1,ssize);
minrange = min(srange);
Z=0;
lowen = min(energy);

if T == 0
    for s = srange
        C_T(s+1-minrange) = coherencef(s,1,coeff,n1,v1,FC,T,N,v_max,PBC);
    end
else
    for i=1:dim
        Z=Z+exp(real(lowen - energy(i))/(kb*T));
    end

    Cs_alpha = zeros(dim,ssize);

    for alpha = 1:dim
        if energy(alpha)-lowen > threshold
            break
        end
        for s = srange
            Cs_alpha(alpha,s+1-minrange) = coherencef(s,alpha,coeff,n1,v1,FC,T,N,v_max,PBC);
        end
        C_T = C_T + Cs_alpha(alpha,:)*exp(-(real(energy(alpha)-lowen))./(kb*T));
    end
    C_T = C_T/Z;
end

srange;
Cs_alpha(1,:);
%{
h = heatmap(srange,[1:10],Cs_alpha(1:10,:))
h.XLabel = 's'
h.YLabel = 'C_i(s)'
h.Title =  'Coherence Function'
%}