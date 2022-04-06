function [Rate] = WMMSEYU_OFDM(H,Nrf,snrLin)
[K,Mt,Nf] = size(H);
alp = K/snrLin;
powNoise = K/snrLin;
%%
FRF = analog(H,Nrf,Mt,alp);
% Fbb = digital(H,alp,FRF);
Fbb = digital1(H,alp,FRF);
for nn = 1:Nf
    trace(FRF*Fbb(:,:,nn)*Fbb(:,:,nn)'*FRF')
end
%%
Rate = 0;
for nn = 1:Nf
    for kk = 1:K
        Qk = powNoise;
        Hef = H(kk,:,nn)*FRF*Fbb(:,kk,nn);
        for uu = setdiff(1:K,kk)
            HB = H(kk,:,nn)*FRF*Fbb(:,uu,nn);           
            Qk = Qk + HB*HB';
        end
        Rate = Rate + log2(1+Qk^(-1)*(Hef)'*Hef);
    end
end
Rate=real(Rate)/Nf;
end
