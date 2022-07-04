function periodic=per(n,N)
%This function takes in an integer n and returns its corresponding index
%between 1 & N under periodic boundary conditions

periodic = rem(n,N);

if periodic <= 0
    periodic = periodic + N;
end

end