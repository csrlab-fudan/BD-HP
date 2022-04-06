function V_D  = digital1(H,alp,V_RF)
[K,Mt,Nf] = size(H);
Mr = K;
Ns = 1;
powNoise = alp;
[~,Nrf] = size(V_RF);
%% effective channel
g = zeros(K,Nrf,Nf);
for nn = 1:Nf
    for kk = 1:K
        g(kk,:,nn) = H(kk,:,nn)*V_RF;
    end
end
%% Digital Initialization
V_D = zeros(Nrf,K,Nf);
V_BB = zeros(Nrf,K,Nf);
a = ones(1,Nf);
for jj = 1:Nf
    a(1,jj) = sqrt(K)/norm(g(:,:,jj)','fro');
    V_BB(:,:,jj) = g(:,:,jj)'*a(1,jj);
end
for jj = 1:Nf
    V_D(:,:,jj)= V_BB(:,:,jj).*(sqrt(K)/norm(V_RF*V_BB(:,:,jj),'fro'));
end
for nn = 1:Nf
    Pt = trace(V_RF*V_D(:,:,nn)*(V_RF*V_D(:,:,nn))');
end
b = zeros(1,Nf);
for n = 1:10
    for jj=1:Nf
        for uu1 = 1:K
            Qk = powNoise;
            Bkn = V_D(:,uu1,jj);
            Heff = g(uu1,:,jj)*Bkn;
            for uu2 = setdiff(1:K,uu1)
                HB = g(uu1,:,jj)*V_D(:,uu2,jj);
                Qk = Qk + HB*HB';
            end
            R((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns) = Heff'*Qk^(-1)*Heff;
            [VV((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns),~] = eig(R((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns));
            A(uu1,uu1,jj) = (Qk+Heff*Heff')^(-1)*Heff;
            W(uu1,uu1,jj) = (1+Qk^(-1)*(Heff)'*Heff);
        end
        v =  powNoise*trace(W(:,:,jj) * A(:,:,jj)' * A(:,:,jj))/(K*Ns);
        H1 = g(:,:,jj)'*A(:,:,jj);
        V_D(:,:,jj) = (H1 * W(:,:,jj) * H1' +  v*eye(Nrf))^(-1)*H1*W(:,:,jj);
        V_D(:,:,jj) = V_D(:,:,jj)*VV;
        b(1,jj) =   sqrt(K)/norm(V_RF*V_D(:,:,jj),'fro');
        V_D(:,:,jj) = V_D(:,:,jj)*b(1,jj);
    end
end
% for nn = 1:Nf
%     trace((V_RF*V_D(:,:,nn))*(V_RF*V_D(:,:,nn))')
% end
