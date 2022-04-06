Mt = 64; Mr = 4; Nrf = 8; Ns = 1; K = 4;
Nf = 32;
numMC = 1;
snrDbSet = -10:5:10;
%%
chanType ='ULA';
phi = pi/3;
theta = pi/6;
zeta = 360/180*pi;
sigma = 360/180*pi;
numMp = 15;
BDFD = zeros(length(snrDbSet),numMC);
BDHP = zeros(length(snrDbSet),numMC);
MMSEFD = zeros(length(snrDbSet),numMC);
MMSEHP = zeros(length(snrDbSet),numMC);
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
        [MMSEFD(indxSnrDb,mm)] = WMMSE(chanMat,Mr,Mt,K,Nf,snrDb,1);
        [MMSEHP(indxSnrDb,mm)] = WMMSEYU_OFDM(chanMat,Nrf,snrLin);
    end
end
%%
figure
width = 2;
plot(snrDbSet,mean(BDFD,2),'k','LineWidth',2),hold on
plot(snrDbSet,mean(MMSEFD,2),'b','LineWidth',2) ,hold on
plot(snrDbSet,mean(BDHP,2),'r','LineWidth',2),hold on
plot(snrDbSet,mean(MMSEHP,2),'b-.','LineWidth',2) ,hold on
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')
title(['M_t = ' num2str(Mt) ', N_{RF} = ' num2str(Nrf) ', K = ' num2str(K) ', N = ' num2str(Nf)]);
 legend('BD-OFDM-FD','WMMSE-OFDM Fully Digital [24]','Proposed BD-OFDM Hybrid （I = 35, T = 15）','WMMSE-OFDM Hybrid (A = 15, D = 20) [10]')