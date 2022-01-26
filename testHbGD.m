function [Rate] = testHbGD(RR,H,Nrf,snrLin,K,Ns)
[MR,Mt,Nf] = size(H);
Mr = MR/K;
M = 10;
N = 30;
sigma2 = (K*Ns)/snrLin;
FRF = RR./abs(RR);
for mm = 1:M
    %% BD
    [URF,~]=qr(FRF,0);
    HR = zeros(MR,Nrf,Nf);
    Heff = zeros(MR,K*Ns,Nf);
    Fbb = zeros(Nrf,K*Mr,Nf);
    for nn = 1:Nf
        HR(:,:,nn) = H(:,:,nn)*URF;
        for kk = 1:K
            Hbar = [HR(1:Mr*(kk-1),:,nn);HR(kk*Mr+1:K*Mr,:,nn)];
            hk = HR(Mr*(kk-1)+1:Mr*kk,:,nn);
            HbarNew = [Hbar;hk];
            HNew = HbarNew';
            [Q,~] = qr(HNew,0);
            Fbb(:,(kk-1)*Ns+1:kk*Ns,nn) = Q(:,(K-1)*Mr+1:K*Mr);
        end
        Heff(:,:,nn) = HR(:,:,nn)*Fbb(:,:,nn);
    end
    %% Water-Filling
    V = zeros(Ns,Ns,K,Nf);
    S = zeros(Ns,Ns,K,Nf);
    QQ = zeros(K*Ns,K*Ns,Nf);
    FBb = zeros(Nrf,Ns*K,Nf);
    alp = zeros(Ns*K,Nf);
    for nn = 1:Nf
        for kk = 1:K
            [~,S(:,:,kk,nn),V(:,:,kk,nn)] = svd(Heff(Mr*(kk-1)+1:kk*Mr,Ns*(kk-1)+1:kk*Ns,nn));
            FBb(:,(kk-1)*Ns+1:kk*Ns,nn) = Fbb(:,(kk-1)*Ns+1:kk*Ns,nn)*V(:,:,kk,nn);
            s = diag(S(:,:,kk,nn));
            sv(:,kk,nn) = s(1:Ns);
        end
        QQ(:,:,nn) = FBb(:,:,nn)'*(URF)'*URF*FBb(:,:,nn);
        alp(:,nn) = diag(QQ(:,:,nn));
    end
    sv2 = (sv(:)).^2;
    q = sv2;
    b = sort(q,'descend');
    svSqur = sv.^2;
    for nn = 1:Nf
        for kk = 1:K
            for ii = 1:Ns
                ind(ii,kk,nn) = find(ismember(b,svSqur(ii,kk,nn)));
            end
        end
    end
    for ii = Nf*K*Ns:-1:1
        mu = (K*Ns*Nf)/ii+sum(sigma2./b(1:ii))/ii;
        if mu > sigma2/b(ii)
            break
        end
    end
    ldpow = zeros(Ns,K,Nf);
    b(ii+1:end)=0;
    for nn = 1:Nf
        for kk = 1:K
            ind_n = ind(:,kk,nn);
            ind_1=find(ismember(ind_n,ind_n(ind_n<ii+1)));
            ldpow(ind_1,kk,nn) = mu-sigma2./b(ind_n(ind_n<ii+1));
        end
    end
    P = zeros(Ns*K,Ns*K,Nf);
    FBB = zeros(Nrf,Ns*K,Nf);
    a = 0;
    for nn = 1:Nf
        for kk = 1:K
            Inv = 1./diag(sqrt(QQ((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn)));
            QQkTem = diag(Inv);
            P((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn) = V(:,:,kk,nn)*QQkTem*sqrt(diag(ldpow(:,kk,nn)));
        end
        FBB(:,:,nn) = Fbb(:,:,nn)*P(:,:,nn);
        a = a+real(trace((URF*FBB(:,:,nn))'*URF*FBB(:,:,nn)));
    end   
    %% Analog
    PF = FRF*inv(FRF'*FRF)*FRF';
    FBBbar = zeros(Mt,K*Ns,Nf);
    for nn = 1:Nf
        FBBbar(:,:,nn) = URF*FBB(:,:,nn);
    end
    for nnn = 1:N
        dfTheta = zeros(Mt,Nrf);
        for nn = 1:Nf
            df1 = 0;
            for kk = 1:K
                Tem = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn);
                A = eye(Mr)+ 1/sigma2*(Tem)*Tem';
                D = FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)'*PF*H((kk-1)*Mr+1:kk*Mr,:,nn)';
                df11 = (FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)*D*A^(-1)*H((kk-1)*Mr+1:kk*Mr,:,nn)*FRF*(FRF'*FRF)^(-1)+H((kk-1)*Mr+1:kk*Mr,:,nn)'*A^(-1)*D'*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)'*FRF*(FRF'*FRF)^(-1));
                df1 = df1 + df11;
            end
            dfTheta = dfTheta+df1;
        end
        %%
        vntN1 = (eye(Mt)-PF)*dfTheta;
        vntN = 2/sigma2*imag(vntN1.*conj(FRF));
        theta = angle(FRF);
        vntN = reshape(vntN,Mt,Nrf);
        thetaNew = theta + vntN;
        FRF = exp(1j*(thetaNew));
        PF = FRF*inv(FRF'*FRF)*FRF';
    end
end
%% Sum-Rate
RateA = 0;
for nn = 1:Nf
    RateS = 0;
    for kk = 1:K
        HE = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn);
        Qk = sigma2*eye(Mr);
        for gg = setdiff(1:K,kk)
            HI = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBBbar(:,(gg-1)*Ns+1:gg*Ns,nn);
            Qk = Qk + HI*HI';
        end
        RateS = RateS + real(log2(det(eye(Mr)+1/sigma2*(HE)*(HE)')));
    end
    RateA = RateA + RateS;
end
Rate = RateA/Nf;
end