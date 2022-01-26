Mt = 64; Mr = 4; Nrf = 8; Ns = 2; K = 2;
Nf = 32;
numMC = 20;
snrDb = 10;
snrLin = db2pow(snrDb);
sigma2 = (K*Ns)/snrLin;
weightedSet1 =-3.5:1:3.5;
%%
chanType ='ULA';
phi = pi/3;
theta = pi/6;
zeta = 360/180*pi;
sigma = 360/180*pi;
numMp = 15;
rate1a = zeros(1,length(weightedSet1));
rate2a = zeros(1,length(weightedSet1));
rate3a = zeros(1,length(weightedSet1));
rate4a = zeros(1,length(weightedSet1));
%%
for mm = 1:numMC
    if mod(mm,10) == 1
        mm
    end
    mm
    %% channel
    chanMat = zeros(Mr,Mt,Nf);
    for kk=1:K
        chanMat(Mr/K*(kk-1)+1:Mr/K*kk,:,:) = channel_generation_ula(phi,theta,Mr/K,Mt,numMp,Nf);
    end
    %% calculate rates
    ratea =  zeros(1,length(weightedSet1));
    rateb =  zeros(1,length(weightedSet1));
    ratec =  zeros(1,length(weightedSet1));
    rated =  zeros(1,length(weightedSet1));
    for w_index = 1:length(weightedSet1)
        w = weightedSet1(w_index);
        [Rate,Fbb,P1,P2] = BDFullyDigitalR(chanMat,snrLin,K,Ns,w);
        rateA = 0;
        for nn = 1:Nf
            HE = chanMat(1:2,:,nn)*Fbb(:,1:2,nn)*P1(:,:,nn);
            rateA = rateA + log2(det(eye(Mr/K)+1/(sigma2)*HE*HE'));
        end
        rateB = 0;
        for nn = 1:Nf
            HE = chanMat(3:4,:,nn)*Fbb(:,3:4,nn)*P2(:,:,nn);
            rateB = rateB + log2(det(eye(Mr/K)+1/(sigma2)*HE*HE'));
        end
        ratea(w_index) = rateA/Nf;
        rateb(w_index) = rateB/Nf;
        A = randn(Mt,Nrf)+1j*randn(Mt,Nrf);
        [Rate1,Rate2,Fbb,P1,P2,FRF] = testHbGDRegion(A,chanMat,Nrf,snrLin,K,Ns,w);
        [URF,~]=qr(FRF,0);
        rateC = 0;
        for nn = 1:Nf
            HE = chanMat(1:2,:,nn)*URF*Fbb(:,1:2,nn)*P1(:,:,nn);
            rateC = rateC + log2(det(eye(Mr/K)+1/(sigma2)*HE*HE'));
        end
        rateD = 0;
        for nn = 1:Nf
            HE = chanMat(3:4,:,nn)*URF*Fbb(:,3:4,nn)*P2(:,:,nn);
            rateD = rateD + log2(det(eye(Mr/K)+1/(sigma2)*HE*HE'));
        end
        ratec(w_index) = rateC/Nf;
        rated(w_index) = rateD/Nf;
    end
    rate1a = rate1a + ratea;
    rate2a = rate2a + rateb;
    rate3a = rate3a + ratec;
    rate4a = rate4a + rated;
end
%%
figure
width = 2;
plot( real(rate1a)/numMC, real(rate2a)/numMC,'g+-', 'LineWidth', 2);
hold on
plot( real(rate3a)/numMC, real(rate4a)/numMC,'r+-', 'LineWidth', 2);
grid on
% axis([0 15 0 15])
xlabel('R_1[bit/channel use]');
ylabel('R_2[bit/channel use]');
title(['M_t = ' num2str(Mt) ', N_{RF} = ' num2str(Nrf) ', K = ' num2str(K) ', N = ' num2str(Nf)]);