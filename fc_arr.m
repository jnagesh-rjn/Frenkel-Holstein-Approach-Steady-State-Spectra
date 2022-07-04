function FC = fc_arr(v_max,v_ground,L)
%Function to calculate FC overlap matrix

FC=zeros(v_ground+1,v_max+1);
FC(1,1)= exp(-L^2/2);
for i=1:min(v_max,v_ground)
    FC(1,i+1) = exp(-L^2/2)*(L^i)/sqrt(factorial(i));
    FC(i+1,1) = (-1)^i * exp(-L^2/2)*(L^i)/sqrt(factorial(i));
end
if v_max < v_ground
    for i=v_max+1:v_ground
        FC(i+1,1) = (-1)^i * exp(-L^2/2)*(L^i)/sqrt(factorial(i));
    end
else
    for i=v_ground+1:v_max
        FC(1,i+1) = exp(-L^2/2)*(L^i)/sqrt(factorial(i));
    end
end

for u=1:v_ground
    for v=1:v_max
        FC(u+1,v+1) = (sqrt(v)*FC(u,v) - L*FC(u,v+1))/sqrt(u);
    end
end

end