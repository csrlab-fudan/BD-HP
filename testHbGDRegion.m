function [Rate1,Rate2,Fbb,P1,P2,FRF] = testHbGD(RR,H,Nrf,snrLin,K,Ns,w)
[MR,Mt,Nf] = size(H);
Mr = MR/K;
M = 20;
N = 30;
sigma2 = (K*Ns)/snrLin;
FRF = RR./abs(RR);
wa = 1;
wb = 10^w;
w2 = [wa,wb];
f = zeros(M,1);
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
%     sum(ldpow1(:))
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
%     sum(ldpow2(:))
    ldpow = zeros(Ns,K,Nf);
    for nn = 1:Nf
        ldpow(:,1,nn) = ldpow1(:,nn);
        ldpow(:,2,nn) = ldpow2(:,nn);
    end
%     sum(ldpow(:))
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
    FBB = zeros(Nrf,Ns*K,Nf);
    for nn = 1:Nf
        FBB(:,:,nn) = Fbb(:,:,nn)*P(:,:,nn);
    end
    %% Analog
    PF = FRF*inv(FRF'*FRF)*FRF';
    FBBbar = zeros(Mt,K*Ns,Nf);
    for nn = 1:Nf
        FBBbar(:,:,nn) = URF*FBB(:,:,nn);
    end
    g = zeros(N,1);
    for nnn = 1:N
        dfTheta = zeros(Mt,Nrf);
        for nn = 1:Nf
            df1 = 0;
            for kk = 1:K
                Tem = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn);
                A = eye(Mr)+ 1/sigma2*(Tem)*Tem';
                D = FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)'*PF*H((kk-1)*Mr+1:kk*Mr,:,nn)';
                df11 = w2(kk)*2/sigma2*(eye(Mt)-PF)*(FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)*D*A^(-1)*H((kk-1)*Mr+1:kk*Mr,:,nn)*FRF*(FRF'*FRF)^(-1)+H((kk-1)*Mr+1:kk*Mr,:,nn)'*A^(-1)*D'*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn)'*FRF*(FRF'*FRF)^(-1));
                dg11 = imag(df11.*conj(FRF));
                df1 = df1 + dg11;
            end
            dfTheta = dfTheta+df1;
        end
        %%
        vntN = dfTheta;
        theta = angle(FRF);
        vntN = reshape(vntN,Mt,Nrf);
        if w>0
            s = 10^(-w);
        else
            s = 1;
        end
        thetaNew = theta + s*vntN;
        FRF = exp(1j*(thetaNew));
        PF = FRF*inv(FRF'*FRF)*FRF';
        %% Inner test Sum-Rate
%         RateA = 0;
%         for nn = 1:Nf
%             RateS = 0;
%             for kk = 1:K
%                 HE = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn);
%                 Qk = sigma2*eye(Mr);
%                 for gg = setdiff(1:K,kk)
%                     HI = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(gg-1)*Ns+1:gg*Ns,nn);
%                     Qk = Qk + HI*HI';
%                 end
%                 RateS = RateS + real(log2(det(eye(Mr)+ 1/sigma2*(HE)*(HE)')));
% %                 RateS = RateS + real(log2(det(eye(Mr)+ 1/sigma2*(HE)*(HE)')));
%             end
%             RateA = RateA + RateS;
%         end
%         Rate = RateA/Nf;
%         g(nnn) = real(Rate);
    end
%                     figure
%                     plot([1:N],g)
    %% Outer test sum-rate
%     RateA = 0;
%             for nn = 1:Nf
%                 RateS = 0;
%                 for kk = 1:K
%                     HE = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(kk-1)*Ns+1:kk*Ns,nn);
%                     Qk = sigma2*eye(Mr);
%                     for gg = setdiff(1:K,kk)
%                         HI = H((kk-1)*Mr+1:kk*Mr,:,nn)*PF*FBBbar(:,(gg-1)*Ns+1:gg*Ns,nn);
%                         Qk = Qk + HI*HI';
%                     end
%                     RateS = RateS + real(log2(det(eye(Mr)+ 1/sigma2*(HE)*(HE)')));
%                 end
%                 RateA = RateA + RateS;
%             end
%             Rate = RateA/Nf;
%             f(mm) = real(Rate);
    %%   test power
%                 PP = zeros(K*Ns,K*Ns,Nf);
%                 pow = 0;
%                 for nn = 1:Nf
%                     PP(:,:,nn) = FBBbar(:,:,nn)'*(PF)'*PF*FBBbar(:,:,nn);
%                     pow = pow+trace(PP(:,:,nn));
%                 end
%                 f(mm) = real(pow);
end% figure
% plot([1:M],f)
PF = FRF*(FRF'*FRF)^(-1)*FRF';
%% Sum-Rate
RateA = 0;
for nn = 1:Nf
    HE = H(1:2,:,nn)*PF*FBBbar(:,1:2,nn);
    RateA = RateA + real(log2(det(eye(Mr)+1/sigma2*(HE)*(HE)')));
end
Rate1 = RateA/Nf;
RateB = 0;
for nn = 1:Nf
    HE = H(3:4,:,nn)*PF*FBBbar(:,3:4,nn);
    RateB = RateB + real(log2(det(eye(Mr)+1/sigma2*(HE)*(HE)')));
end
Rate2 = RateB/Nf;
end