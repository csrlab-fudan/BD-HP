function V_RF  = analog(H,Nrf,Mt,alp)
Vn = alp;
Nt = Mt;
V_RF = ones(Nt,Nrf);
[K,~,Nk] = size(H);
for nn = 1:Nk
    F(:,:,nn) = H(:,:,nn)'*H(:,:,nn);
end
F = sum(F,3)/Nk;
g = K/Nrf/Nt;
a = g/Vn;    
x = 0;
for Nloop = 1:9
    for jj = 1:Nrf
        VRF = V_RF;
        VRF(:,jj)=[];
        C = eye(Nrf-1)+a*VRF'*F*VRF;
        G = a*F-a^2*F*VRF*C^(-1)*VRF'*F;
        for i = 1:Nt
            for l = 1:Nt
                if i~=l
                    x(l)=G(i,l)*V_RF(l,jj);
                end
            end
            n = sum(x);
            if n ==0
                V_RF(i,jj)=1;
            else
                V_RF(i,jj)=n/abs(n);
            end
        end
    end
end