function [Rate] = BDFullyDigitalO(H,snrLin,K,Ns)
[MR,Mt,Nf] = size(H);
Mr = MR/K;
Fbb = zeros(Mt,Ns*K,Nf);
Heff = zeros(K*Mr,K*Ns,Nf);
C2 = zeros(K*Mr,K*Ns,Nf);
sigma2 = (K*Ns)/snrLin;
%% BD
for nn = 1:Nf
    for kk = 1:K
        Hbar = [H(1:Mr*(kk-1),:,nn);H(kk*Mr+1:K*Mr,:,nn)];
        hk = H(Mr*(kk-1)+1:Mr*kk,:,nn);
        HbarNew = [Hbar;hk];
        HNew = HbarNew';
        [Q,R] = qr(HNew,0);
        C2(Mr*(kk-1)+1:kk*Mr,Ns*(kk-1)+1:kk*Ns,nn) = R((K-1)*Mr+1:K*Mr,(K-1)*Mr+1:K*Mr);
        Fbb(:,(kk-1)*Ns+1:kk*Ns,nn) = Q(:,(K-1)*Mr+1:K*Mr);
    end
    Heff(:,:,nn) = H(:,:,nn)*Fbb(:,:,nn);
end
%% Water-Filling
V = zeros(Ns,Ns,K,Nf);
S = zeros(Ns,Ns,K,Nf);
for nn = 1:Nf
    for kk = 1:K
        [~,S(:,:,kk,nn),V(:,:,kk,nn)] = svd(Heff(Mr*(kk-1)+1:kk*Mr,Ns*(kk-1)+1:kk*Ns,nn));
        s = diag(S(:,:,kk,nn));
        sv(:,kk,nn) = s(1:Ns);
    end
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
%% Test
t = 0;
P = zeros(Ns*K,Ns*K,Nf);
for nn = 1:Nf
    for kk = 1:K
        P((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn) = V(:,:,kk,nn)*sqrt(diag(ldpow(:,kk,nn)));
    end
    t = t+trace(P(:,:,nn)'*Fbb(:,:,nn)'*Fbb(:,:,nn)*P(:,:,nn));
end
%% Sum-Rate
RateA = 0;
for nn = 1:Nf
    RateS = 0;
    for kk = 1:K
        HE = H((kk-1)*Ns+1:kk*Ns,:,nn)*Fbb(:,(kk-1)*Ns+1:kk*Ns,nn)*P((kk-1)*Ns+1:kk*Ns,(kk-1)*Ns+1:kk*Ns,nn);
        Qk = 0;
        for gg = setdiff(1:K,kk)
            HI = H((kk-1)*Ns+1:kk*Ns,:,nn)*Fbb(:,(gg-1)*Ns+1:gg*Ns,nn)*P((gg-1)*Ns+1:gg*Ns,(gg-1)*Ns+1:gg*Ns,nn);
            Qk = Qk + HI*HI';
        end
        RateS = RateS + real(log2(det(eye(Ns)+1/(sigma2)* (HE)'*(HE))));
    end
    RateA = RateA + RateS;
end
Rate = RateA/Nf;
end