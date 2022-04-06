function [Rate]= WMMSE(H,Mr,Mt,K,Nf,snrDb,Ns)
h = zeros(Mr/K,Mt,K,Nf);
A = zeros(Mr,K*Ns,Nf);
W = zeros(K*Ns,K*Ns,Nf);
B = zeros(Mt,K*Ns,Nf);
snrLin = db2pow(snrDb);
powNoise =  K/snrLin;
for jj = 1:Nf
    for ii = 1:K
        h(:,:,ii,jj) = H((ii-1)*(Mr/K)+1:ii*(Mr/K),:,jj);
    end
end
a = ones(1,Nf);
for jj = 1:Nf
    a(1,jj) = sqrt(K)/norm(H(:,:,jj)','fro');
    B(:,:,jj)= H(:,:,jj)'*a(1,jj);
end
% for jj = 1:Nf
%     [~,~,V] = svd(H(:,:,jj));
%     B(:,:,jj)= V(:,1:K*Ns);
% end
for nn = 1:Nf
    Pt = trace(B(:,:,nn)*B(:,:,nn)');
end
R = zeros(K*Ns,K*Ns);
b = zeros(1,Nf);
for n = 1:30
    for jj=1:Nf
        for uu1 = 1:K
            Qk = powNoise*eye(Mr/K);
            Bkn = B(:,(uu1-1)*Ns+1:uu1*Ns,jj);
            Heff = h(:,:,uu1,jj)*Bkn;
            for uu2 = setdiff(1:K,uu1)
                HB = h(:,:,uu1,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj);
                Qk = Qk + HB*HB';
            end
            R((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns) = Heff'*Qk^(-1)*Heff;
            [VV((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns),~] = eig(R((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns));
            A((uu1-1)*(Mr/K)+1:uu1*(Mr/K),(uu1-1)*Ns+1:uu1*Ns,jj) = (Qk+Heff*Heff')^(-1)*Heff;
            W((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns,jj) = (eye(Ns)+Heff'*Qk^(-1)*Heff);
        end
        v =  powNoise*trace(W(:,:,jj) * A(:,:,jj)' * A(:,:,jj))/(K*Ns);
        H1 = H(:,:,jj)'*A(:,:,jj);
        B(:,:,jj) = (H1 * W(:,:,jj) * H1' +  v*eye(Mt))^(-1)*H1*W(:,:,jj);
        B(:,:,jj) = B(:,:,jj)*VV;
        b(1,jj) =  sqrt(K*Ns)/norm(B(:,:,jj),'fro');
        B(:,:,jj) = B(:,:,jj)*b(1,jj);
    end
end

rate = zeros(Ns*K,1);
Rate = 0;
for jj = 1:Nf
    for uu1 = 1:K
        Qk=powNoise;
        Bk=B(:,(uu1-1)*Ns+1:uu1*Ns,jj);
        Hef=h(:,:,uu1,jj)*Bk;
        for uu2 = setdiff(1:K,uu1)
            HB = h(:,:,uu1,jj)*B(:,(uu2-1)*Ns+1:uu2*Ns,jj);           
            Qk = Qk + HB*HB';
        end
        Rate = Rate + log2(1+Qk^(-1)*Hef'*Hef);
        E((uu1-1)*Ns+1:uu1*Ns,(uu1-1)*Ns+1:uu1*Ns) = (eye(Ns)+Hef'*Qk^(-1)*Hef)^(-1);%
        mse = real(diag(E));
    end
    rate = rate + log2(1./mse);
end
Rate=real(Rate)/Nf;
for nn = 1:Nf
    Pt = trace(B(:,:,nn)*B(:,:,nn)');
end