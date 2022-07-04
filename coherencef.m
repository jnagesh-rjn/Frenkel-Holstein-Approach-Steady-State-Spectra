function C_s = coherencef(s,alpha,coeff,n1,v1,FC,T,N,v_max,PBC)
%Calculates the value of the coherence function for the Linear Aggregate at
%finite temperature for a given (alphath) eigenstate
%s is the separation. So the function returns C(s)

C_s = 0;
%Only one-particle states will contribute
if PBC == 0 %Open boundary conditions
    for n=1:N
    	if n+s <= 0 || n+s > N
            continue
        end
        for i = 0:v_max
            for j = 0:v_max
                C_s = C_s + conj(coeff((n-1)*(v_max+1)+1+j,alpha)*FC(1,j+1)) * coeff((n+s-1)*(v_max+1)+1+i,alpha) * FC(1,i+1);
            end
        end
    end
else        %Periodic Boundary Conditions
    for n = 1:N
    	for i = 0:v_max
            for j = 0:v_max
                C_s = C_s + conj(coeff((n-1)*(v_max+1)+1+j,alpha)*FC(1,1+j)) * coeff((per(n+s,N)-1)*(v_max+1)+1+i,alpha) * FC(1,i+1);
            end
        end
    end
end
