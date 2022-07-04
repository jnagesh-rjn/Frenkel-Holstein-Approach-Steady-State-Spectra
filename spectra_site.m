%clear
%tic
variables;

buildHam;

%Diagonalise the Hamiltonian
[coeffun, e] = eig(Ham); 
%energy = diag(e); % This turns energy into a column vector, instead of a diagonal matrix
[energy,perm] = sort(diag(e));
coeff=coeffun(:,perm);
%energy(1:3)

%calculate spectra- abs and emi
%w = w_00+w_D-4*w_0:0.0001:w_00+w_D+4*w_0;
Z=0;
if T>0
   for i=1:dim
       Z=Z+exp((energy(1)-energy(i))/(kb*T));
   end
end

cal_spec;

%plots and output files
wn=w*modefreq*au2cminverse/au2ev;
wn=wn+w00;
wDwn=w_D*modefreq/cm2ev;


