function I=linestrength(vt,alpha,coeff,N,v_max,FC,n1,n2,v1,v2,dim1,dim2,tpcheck)
%This function calculates the linestrength for the emission from the alpha
%eigenstate to the ground state with vt quanta

I=0;

if vt == 0
    %There is only contribution from the one particle states
    for i=1:dim1
        I = I + conj(coeff(i,alpha)) * conj(FC(1,v1(i)+1));
    end
    
    I = abs(I)^2;
    
else
    %The contribution from the states where there are vt quanta left on a
    %single site
    I1=0;
    for m=1:N
        I_m=0;
        for i=0:v_max
            I_m = I_m + conj(coeff((m-1)*(v_max+1)+1+i,alpha)) * conj(FC(vt+1,i+1));
        end
        if tpcheck == 1
            for i=1:dim2
                if n2(i,2) == m && v2(i,2) == vt
                    I_m = I_m + conj(coeff(i+dim1,alpha)) * conj(FC(1,v2(i,1)+1));
                end
            end
        end
        I_m = abs(I_m)^2;
        I1 = I1 + I_m;
    end
    I = I + I1;
    
    %The contribution I(2) coming from the v1 and v2 quanta spread over two
    %sites
    if tpcheck == 1 && vt ~= 1
        %vsum=v_max*(v_max+1)*0.5;
        for v_1=1:floor(vt/2)
            v_2=vt-v_1;
            I2=0;
            for m1=1:N
                for m2=1:N
                    if m1 == m2
                        continue
                    end
                    I_mm=0;
                    for i=1:dim2 %parfor
                        if n2(i,:) == [m1,m2] & v2(i,2) == v_2
                            I_mm = I_mm + conj(coeff(i+dim1,alpha) * FC(v_1+1,v2(i,1)+1));
                        %I_mm = I_mm + coeff(dim1+(m1-1)*dim2/N+(m2-1)*vsum+i*(v_max-(i-1)*0.5)+v_2,alpha) * conj(FC(v_1+1,i+1));
                        end
                        if n2(i,:) == [m2,m1] & v2(i,2) == v_1
                            I_mm = I_mm + conj(coeff(i+dim1,alpha) * FC(v_2+1,v2(i,1)+1));
                        end
                    end
                    %{
                    for i=0:v_max-v_1
                        I_mm = I_mm + coeff(dim1+(m2-1)*dim2/N+(m1-1)*vsum+i*(v_max-(i-1)*0.5)+v_1,alpha) * conj(FC(v_2+1,i+1));
                    end
                    %}
                    I_mm = abs(I_mm)^2;
                    I2 = I2 + I_mm;
                end
            end
            if v_1 == v_2
                I = I + 0.5*I2;
            else
                I = I + I2;
            end
        end
    end

end