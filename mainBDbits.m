Mt = 64; Mr = 6; Nrf = 8; Ns = 2; K = 3;
Nf = 32;
numMC = 20;
snrDbSet = -20:5:10;
%%
chanType ='ULA';
phi = pi/3;
theta = pi/6;
zeta = 360/180*pi;
sigma = 360/180*pi;
numMp = 15;
BDFD = zeros(length(snrDbSet),numMC);
BDHP = zeros(length(snrDbSet),numMC);
BDHP1 = zeros(length(snrDbSet),numMC);
BDHP2 = zeros(length(snrDbSet),numMC);
BDHP3 = zeros(length(snrDbSet),numMC);
BDHP4 = zeros(length(snrDbSet),numMC);
%%
for mm = 1:numMC
    %% channel
    chanMat = zeros(Mr,Mt,Nf);
    for kk=1:K
        chanMat(Mr/K*(kk-1)+1:Mr/K*kk,:,:) = channel_generation_ula(phi,theta,Mr/K,Mt,numMp,Nf);
    end
    %% calculate rates
    for indxSnrDb = 1:length(snrDbSet)
        snrDb = snrDbSet(indxSnrDb);
        snrLin = db2pow(snrDb);
        A = randn(Mt,Nrf)+1j*randn(Mt,Nrf);
        BDFD(indxSnrDb,mm) = BDFullyDigitalO(chanMat,snrLin,K,Ns);
        BDHP(indxSnrDb,mm) = testHbGD(A,chanMat,Nrf,snrLin,K,Ns);
        BDHP1(indxSnrDb,mm) = testHbGD1bit(A,chanMat,Nrf,snrLin,K,Ns);
        BDHP2(indxSnrDb,mm) = testHbGD2bit(A,chanMat,Nrf,snrLin,K,Ns);
        BDHP3(indxSnrDb,mm) = testHbGD3bit(A,chanMat,Nrf,snrLin,K,Ns);
    end
end
%%
figure
width = 2;
plot(snrDbSet,mean(BDFD,2),'k','LineWidth',2),hold on
plot(snrDbSet,mean(BDHP,2),'r','LineWidth',2),hold on
plot(snrDbSet,mean(BDHP1,2),'b--','LineWidth',2),hold on
plot(snrDbSet,mean(BDHP2,2),'b:','LineWidth',2),hold on
plot(snrDbSet,mean(BDHP3,2),'b-.','LineWidth',2),hold on
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
title(['M_t = ' num2str(Mt) ', N_{RF} = ' num2str(Nrf) ', K = ' num2str(K) ', N = ' num2str(Nf)]);
legend('BD-OFDM-FD','BD-OFDM-HPC','BD-OFDM-HPC (1-bit)','BD-OFDM-HPC (2-bit)','BD-OFDM-HPC (3-bit)')