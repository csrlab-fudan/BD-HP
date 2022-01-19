function [Rate] = fuctionCG(H,PF,FBB,Nf,K,Ns,sigma2)
RateA = 0;
Mr = Ns;
for nn = 1:Nf
    RateS = 0;
    for kk = 1:K
        HE = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBB(:,(kk-1)*Ns+1:kk*Ns,nn);
        Qk = sigma2*eye(Mr);
        for gg = setdiff(1:K,kk)
            HI = H((kk-1)*Ns+1:kk*Ns,:,nn)*PF*FBB(:,(gg-1)*Ns+1:gg*Ns,nn);
            Qk = Qk + HI*HI';
        end
         RateS = RateS + real(log2(det(eye(Ns)+(HE)'*Qk^(-1)*(HE))));
    end
    RateA = RateA + RateS;
end
Rate = RateA/Nf;
end