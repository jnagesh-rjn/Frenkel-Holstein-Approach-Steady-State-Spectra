function Ham=Hamiltonian1(J,FC,n1,v1,dim1,w_00,w_D,w_0,w_1,siteshift)
%This function constructs and returns the Hamiltonian in the one-particle
%basis (2-particle states are excluded)

Ham=zeros(dim1);
for i=1:dim1
    for j=1:dim1
        if n1(i) == n1(j) && v1(i) == v1(j)
            Ham(i,j) = w_00 + w_D + w_1*v1(i)+siteshift(n1(i),1);
        end
        Ham(i,j) = Ham(i,j) + J(n1(i),n1(j))*conj(FC(1,v1(i)+1))*FC(1,v1(j)+1); % The electronic coupling part.
    end
end

end
