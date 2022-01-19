%% EGT
Mt = 64; Mr = 8; Nrf = 8; Ns =2; K = 4;
Nf = 32;
numMC = 10;
snrDbSet = -20:5:10;
%%
chanType ='ULA';
phi = pi/3;         
theta = pi/6;      
zeta = 360/180*pi;  
sigma = 360/180*pi; 
numMp = 15;       
BDHP = zeros(length(snrDbSet),numMC);
BDHPO = zeros(length(snrDbSet),numMC);
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
    for indxSnrDb = 1:length(snrDbSet)
        snrDb = snrDbSet(indxSnrDb);
        snrLin = db2pow(snrDb);
        A = randn(Mt,Nrf)+1j*randn(Mt,Nrf);
        BDHP(indxSnrDb,mm) = testEGT(chanMat,snrLin,Nrf,K,Ns);
        BDHPO(indxSnrDb,mm) = testHbGDO(A,chanMat,Nrf,snrLin,K,Ns);
    end
end
%%
figure
width = 2;
plot(snrDbSet,mean(BDHP,2),'r','LineWidth',2),hold on
plot(snrDbSet,mean(BDHPO,2),'b--','LineWidth',2),hold on
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')