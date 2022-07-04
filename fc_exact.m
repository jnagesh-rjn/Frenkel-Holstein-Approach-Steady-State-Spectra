function FC = fc_exact(u, v, L)

% u is the phonon number in the unshifted well.
% v is the phonon number in the shifted well.

if u == 0 
    FC = exp(-L^2/2)*(L^v)/sqrt(factorial(v));
elseif v == 0
    FC = fc_exact(v,u,-L);
else
    FC = (sqrt(v)*fc_exact(u-1,v-1,L) - L*fc_exact(u-1,v,L))/sqrt(u);
end    

end