function [Rate,Fbb,P1,P2] = BDFullyDigitalR(H,snrLin,K,Ns,w)
[MR,Mt,Nf] = size(H);
Mr = MR/K;
Fbb = zeros(Mt,Ns*K,Nf);
Heff = zeros(K*Mr,K*Ns,Nf);
C2 = zeros(K*Mr,K*Ns,Nf);
sigma2 = (K*Ns)/snrLin;
wa = 1;
wb = 10^w;
w2=[wa,wb];
%% BD
for nn = 1:Nf
    for kk = 1:K
        Hbar = [H(1:Mr*(kk-1),:,nn);H(kk*Mr+1:K*Mr,:,nn)];
        hk = H(Mr*(kk-1)+1:Mr*kk,:,nn);
        HbarNew = [Hbar;hk];
        HNew = HbarNew';
        [Q,R] = qr(HNew,0);
        C2(Mr*(kk-1)+1:kk*Mr,Ns*(kk-1)+1:kk*Ns,nn) = R((K-1)*Mr+1:K*Mr,(K-1)*Mr+1:K*Mr);
        Fbb(:,(kk-1)*Ns+1:kk*Ns,nn) = Q(:,(K-1)*Ns+1:K*Ns);
    end
    Heff(:,:,nn) = H(:,:,nn)*Fbb(:,:,nn);
end
%% Water-Filling
V = zeros(Ns,Ns,K,Nf);
P1 = Nf*K*Ns*wa/(wa+wb);
P2 = Nf*K*Ns*wb/(wa+wb); 
sv1 = zeros(Ns,Nf);
V1 = zeros(Ns,Ns,Nf);
V2 = zeros(Ns,Ns,Nf);
for nn = 1:Nf
    [~,S1,V1(:,:,nn)] = svd(Heff(1:Mr,1:Ns,nn));
    s1 = diag(S1);
    sv1(:,nn) = s1;
end
sv12 = sv1.^2;
b = sort(sv12(:),'descend');
svSqur = sv1.^2;
for nn = 1:Nf
        for ii = 1:Ns
            ind(ii,nn) = find(ismember(b,svSqur(ii,nn)));
       end
end
for ii = Ns*Nf:-1:1
    mu = P1/ii+sum(sigma2./b(1:ii))/ii;
    if mu > sigma2/b(ii)
        break
    end 
end
ldpow1 = zeros(Ns,Nf);
b(ii+1:end)=0;
for nn = 1:Nf
    ind_n = ind(:,nn);
    len = length(ind_n(ind_n<ii+1));
    ldpow1(1:len,nn) = mu-sigma2./b(ind(1:len,nn));
end
% sum(ldpow1(:))
sv2 = zeros(Ns,Nf);
for nn = 1:Nf
    [~,S2,V2(:,:,nn)] = svd(Heff(Mr+1:2*Mr,Ns+1:2*Mr,nn));
    s2 = diag(S2);
    sv2(:,nn) = s2;
end
sv22 = sv2.^2;
b2 = sort(sv22(:),'descend');
svSqur2 = sv2.^2;
for nn = 1:Nf
        for ii = 1:Ns
            ind(ii,nn) = find(ismember(b2,svSqur2(ii,nn)));
       end
end
for ii = Ns*Nf:-1:1
    mu = P2/ii+sum(sigma2./b2(1:ii))/ii;
    if mu > sigma2/b2(ii)
        break
    end 
end
ldpow2 = zeros(Ns,Nf);
b(ii+1:end)=0;
for nn = 1:Nf
    ind_n = ind(:,nn);
    len = length(ind_n(ind_n<ii+1));
    ldpow2(1:len,nn) = mu-sigma2./b2(ind(1:len,nn));
end
% sum(ldpow2(:))
ldpow = zeros(Ns,K,Nf);
for nn = 1:Nf
    ldpow(:,1,nn) = ldpow1(:,nn);
    ldpow(:,2,nn) = ldpow2(:,nn);
end
% sum(ldpow(:))
for nn = 1:Nf
    V(:,:,1,nn) = V1(:,:,nn);
    V(:,:,2,nn) = V2(:,:,nn);
end
%% Test
t = 0;
t1=0;
P1 = zeros(Ns,Ns,Nf);
P2 = zeros(Ns,Ns,Nf);
P = zeros(Ns*K,Ns*K,Nf);
for nn = 1:Nf
    for kk = 1:K
        P1(:,:,nn) = V1(:,:,nn)*sqrt(diag(ldpow1(:,nn)));
        P2(:,:,nn) = V2(:,:,nn)*sqrt(diag(ldpow2(:,nn)));
    end
    t1 = t1+trace(P1(:,:,nn)'*P1(:,:,nn)+P2(:,:,nn)'*P2(:,:,nn));
    t = t+trace(sqrt(diag(ldpow1(:,nn)))'*sqrt(diag(ldpow1(:,nn)))+sqrt(diag(ldpow2(:,nn)))'*sqrt(diag(ldpow2(:,nn))));
    P(1:Ns,1:Ns,nn)=  P1(:,:,nn);
    P(Ns+1:2*Ns,Ns+1:2*Ns,nn)=  P2(:,:,nn);
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
        RateS = RateS + w2(kk)*real(log2(det(eye(Ns)+1/(sigma2)*(HE)'*(HE))));
    end
    RateA = RateA + RateS;
end
Rate = RateA/Nf;
end