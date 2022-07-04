%Constructing the J matrix
% IMPORTANT: J(i,i) should be 0!   
J = zeros(N);

%The remainder with N is used to incorporate periodic boundary conditions
%so that the Nth momomer is coupled with the 1st monomer
if twoJ==1 && N==3
J(1,2) = J1;
J(2,1) = J1;
J(2,3) = J0;
J(3,2) = J0;
J(3,1) = J1/2;
J(1,3) = J1/2;
else 
 if PBC==1
 for i = 1:N
    J(i,rem(i,N)+1) = J0;
    J(rem(i,N)+1,i) = J0;
 end
 else
 %open boundary conditions
    for i=1:N
    J(i,i+1) = J0;
    J(i+1,i) = J0;
    end
 end
end

%Constructing the basis
%We will have two n matrices, one for the 1-particle states, and another
%for the 2-particle states
%Similarly also two v matrices
dim1 = N*(v_max+1); % Including zeros.
n1 = zeros(dim1,1); % Contains the n for each 1-particle state.
v1 = zeros(dim1,1); % Contains the ~v for each 1-particle state.
n2=double.empty;
v2=double.empty;
dim2=0;

for i = 1:dim1
    n1(i) = 1 + floor((i-1)/(v_max+1));
    % There's an (i-1) instead of i, as it increased n on the v_max 
    % phonon state, instead of the following 1 phonon state.
    
    v1(i) = rem(i,(v_max+1)); % Sets v(i) to 1,2,3,...,v_max,0.
    if v1(i) == 0
        v1(i) = v_max+1; % Sets that last 0 phonon state to v_max+1 phonons.
    end
    v1(i) = v1(i) - 1; % Brings v(i) to 0,1,2,3...v_max.
end
%The two particle part of the basis will be constructed if it is required,
%just before constructing the Hamiltonian

%Calculating the FC factors
FC = fc_difffreq(freqratio,v_max,v_ground,L);

%Adding site disorder
%siteshift = sigma*randn(N,1);
%Disorder with correlated shifts
siteshift = zeros(N,1);
%
SIGMA = zeros(N);
for i = 1:N
    for j = 1:N
        if i == j
            SIGMA(i,j) = sigma^2;
        else
        SIGMA(i,j) = (sigma^2) * exp(-abs(i-j)/l_0);
        end
    end
end
siteshift = mvnrnd(zeros(1,N),SIGMA,1);
siteshift = siteshift.';

%Constructing the Hamiltnonian matrix
%tpcheck is used to decide whether to include 2-particle states or not and
%corresponding function is called
if tpcheck == 0
    Ham = Hamiltonian1(J,FC,n1,v1,dim1,w_00,w_D,w_0,w_1,siteshift);
    dim=dim1;
else %The basis for the two particle states is also constructed here when it is needed
    vsum=v_max*(v_max+1)*0.5;
    dim2=N*(N-1)*vsum;
    n2 = zeros(dim2,2); % Contains the n,m for each 2-particle state.
    v2 = zeros(dim2,2); % Contains the ~v,v' for each 2-particle state.
    %SELF-NOTE: OPTIMISE THIS SECTION TO CONSTRUCT THE BASIS IN A MORE EFFICIENT
    %MANNER!!! This is probably not a huge problem in terms of time though
    for i=1:dim2
        n2(i,1)=ceil(i*N/dim2);
    end
    for i=1:((N-1)*vsum):dim2
        seq=setdiff(1:N,n2(i,1));
        for j=1:N-1
            for k=0:vsum-1
                n2(i+(j-1)*vsum+k,2)=seq(j);
            end
        end
    end
    for i=1:vsum:dim2
        for j=0:v_max
            for k=1:v_max-j
                v2(i-1+(j*v_max)-0.5*j*(j-1)+k,1) = j;
                v2(i-1+(j*v_max)-0.5*j*(j-1)+k,2) = k;
            end
        end
    end
    dim=dim1+dim2;
    Ham = Hamiltonian2(J,FC,n1,n2,v1,v2,dim1,dim2,w_00,w_D,w_0,w_1,siteshift);
end

