function Ham=Hamiltonian2(J,FC,n1,n2,v1,v2,dim1,dim2,w_00,w_D,w_0,w_1,siteshift)
%This function calculates and returns the Hamiltonian matrix in the basis
%including 1-particle and 2-particle states

dim=dim1+dim2;
Ham=zeros(dim);
Ham(1:dim1,1:dim1) = Hamiltonian1(J,FC,n1,v1,dim1,w_00,w_D,w_0,w_1,siteshift);

%The two particle diagonal block
for k=dim1+1:dim
    for l=dim1+1:dim
        i=k-dim1;
        j=l-dim1;
        if n2(i,:) == n2(j,:) & v2(i,:) == v2(j,:) %Checking that it is a diagonal term
            Ham(k,l) = w_00 +  w_D + w_1*v2(i,1)+ w_0*v2(i,2) + siteshift(n2(i,1));
        end
        if n2(i,2) == n2(j,2) & v2(i,2) == v2(j,2) 
            Ham(k,l) = Ham(k,l) + J(n2(i,1),n2(j,1))*conj(FC(1,v2(i,1)+1))*FC(1,v2(j,1)+1);
        elseif n2(i,1) == n2(j,2) && n2(i,2) == n2(j,1)
            Ham(k,l) = Ham(k,l) + J(n2(i,1),n2(j,1))*conj(FC(v2(j,2)+1,v2(i,1)+1))*FC(v2(i,2)+1,v2(j,1)+1);
        end
    end
end

%Calculating the 1-p 2-p matrix elements
for i=1:dim1
    for j=dim1+1:dim
        if n1(i) == n2(j-dim1,2)
            Ham(i,j) = J(n1(i),n2(j-dim1,1))*conj(FC(v2(j-dim1,2)+1,v1(i)+1))*FC(1,v2(j-dim1,1)+1);
            Ham(j,i) = conj(Ham(i,j));
        end
    end
end

end
