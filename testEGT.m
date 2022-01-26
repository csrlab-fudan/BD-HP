function [Rate] = testEGT(H,snrLin,Nrf,K,Ns)
[MR,Mt,Nf] = size(H);
Mr = MR/K;
sigma2 = (K*Ns)/snrLin;
Hsum = zeros(MR,Mt);
for nn=1:Nf
    Hsum = Hsum + H(:,:,nn);
end
HsumBar = 1/Nf*Hsum;
FRF = exp(1j*angle(HsumBar'));
PF = FRF*(FRF'*FRF)^(-1)*FRF';
HR = zeros(MR,Mt,Nf);
Heff = zeros(MR,K*Ns,Nf);
Fbb = zeros(Mt,K*Mr,Nf);
for nn = 1:Nf
    HR(:,:,nn) = H(:,:,nn)*PF;
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
FBb = zeros(Mt,Ns*K,Nf);
for nn = 1:Nf
    for kk = 1:K
        [~,S(:,:,kk,nn),V(:,:,kk,nn)] = svd(Heff(Mr*(kk-1)+1:kk*Mr,Ns*(kk-1)+1:kk*Ns,nn));
        FBb(:,(kk-1)*Ns+1:kk*Ns,nn) = Fbb(:,(kk-1)*Ns+1:kk*Ns,nn)*V(:,:,kk,nn);
        s = diag(S(:,:,kk,nn));
        sv(:,kk,nn) = s(1:Ns);
    end
    QQ(:,:,nn) = FBb(:,:,nn)'*(PF)'*PF*FBb(:,:,nn);
end
sv2 = (sv(:)).^2;
b = sort(sv2,'descend');
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
        len = length(ind_n(ind_n<ii+1));
        ldpow(1:len,kk,nn) = mu-sigma2./b(ind(1:len,kk,nn));
    end
end
P = zeros(Ns*K,Ns*K,Nf);
FBB = zeros(Mt,Ns*K,Nf);
for nn = 1:Nf
    for kk = 1:K
        Inv = 1./diag(sqrt(QQ((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn)));
        QQkTem = diag(Inv);
        P((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn) = V(:,:,kk,nn)*QQkTem*sqrt(diag(ldpow(:,kk,nn)));
    end
    FBB(:,:,nn) = Fbb(:,:,nn)*P(:,:,nn);
end
%% Sum-Rate
PF = FRF*(FRF'*FRF)^(-1)*FRF';
RateA = 0;
for nn = 1:Nf
    RateS = 0;
    for kk = 1:K
        HE = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBB(:,(kk-1)*Ns+1:kk*Ns,nn);
        Qk = sigma2*eye(Mr);
        for gg = setdiff(1:K,kk)
            HI = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBB(:,(gg-1)*Ns+1:gg*Ns,nn);
            Qk = Qk + HI*HI';
        end
        RateS = RateS + real(log2(det(eye(Mr)+1/sigma2*(HE)*(HE)')));
    end
    RateA = RateA + RateS;
end
Rate = RateA/Nf;
end