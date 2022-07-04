function FC = fc_difffreq(freqratio,v_max,v_ground,L)
%This function calculates and returns an array containing the FC factors in
%the case where the S0 and S1 potentials have different curvatures. i.e.
%w_vib is different for ground and excited potential

FC = zeros(v_ground+2,v_max+1);

%Calculating <n|0~> overlaps
c = zeros(1,100);
c(1) = 1;
c(2) = -(2 * L * sqrt(freqratio) * c(1)) / ((freqratio+1) * sqrt(1));
for i = 2:100-1
    c(i+1) = -((freqratio-1) * sqrt(i-1) * c(i-1) + 2 * L * sqrt(freqratio) * c(i)) / ((freqratio+1) * sqrt(i));
end
c = c./norm(c);

FC(:,1) = c(1:v_ground+2);

%Calculating <0|n~> overlaps

Lother = -L/sqrt(freqratio);
freqreci = 1/freqratio;

c = zeros(1,100);
c(1) = 1;
c(2) = -(2 * Lother * sqrt(freqreci) * c(1)) / ((freqreci+1) * sqrt(1));
for i = 2:100-1
    c(i+1) = -((freqreci-1) * sqrt(i-1) * c(i-1) + 2 * Lother * sqrt(freqreci) * c(i)) / ((freqreci+1) * sqrt(i));
end
c = c./norm(c);

FC(1,2:v_max+1) = c(2:v_max+1);

%Calculating <m|n~> overlaps

for i = 1:v_ground
    for j = 1:v_max
        FC(i+1,j+1) = 0.5 * sqrt(i/j) * (sqrt(freqratio)+sqrt(freqreci)) * FC(i,j) + 0.5 * sqrt((i+1)/j) * (sqrt(freqratio)-sqrt(freqreci)) * FC (i+2,j) + L * FC(i+1,j)/sqrt(j);
    end
end

FC = FC(1:v_ground+1,1:v_max+1);

end