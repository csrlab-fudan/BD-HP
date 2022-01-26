Mt = 64;  Ns = 2;
Nf = 32;
numMC = 20;
userSet = 1:1:8;
%%
chanType ='ULA';
phi = pi/3;         
theta = pi/6;      
zeta = 360/180*pi;  
sigma = 360/180*pi; 
numMp = 15;         
BDFD = zeros(length(userSet),numMC);
BDHP = zeros(length(userSet),numMC);
BDHP1 = zeros(length(userSet),numMC);
BDHP2 = zeros(length(userSet),numMC);
%%
for mm = 1:numMC
    if mod(mm,10) == 1
        mm
    end
    mm
    %% calculate rates
    for indxUser = 1:length(userSet)
        K = userSet(indxUser);
        Nrf = 32;
        Mr = K*2;
        %% channel
        chanMat = zeros(Mr,Mt,Nf);
        for kk=1:K
            chanMat(Mr/K*(kk-1)+1:Mr/K*kk,:,:) = channel_generation_ula(phi,theta,Mr/K,Mt,numMp,Nf);
        end
        snrDb = 0;
        snrLin = db2pow(snrDb);
        A = randn(Mt,Nrf)+1j*randn(Mt,Nrf);
        BDFD(indxUser,mm) = BDFullyDigitalO(chanMat,snrLin,K,Ns);
        BDHP(indxUser,mm) = testHbGDUsers(A,chanMat,Nrf,snrLin,K,Ns);
        Nrf1 = 24;
        B = randn(Mt,Nrf1)+1j*randn(Mt,Nrf1);
        BDHP1(indxUser,mm) = testHbGDUsers(B,chanMat,Nrf1,snrLin,K,Ns);
        Nrf2 = 16;
        C = randn(Mt,Nrf2)+1j*randn(Mt,Nrf2);
        BDHP2(indxUser,mm) = testHbGDUsers(C,chanMat,Nrf2,snrLin,K,Ns);
    end
end
%%
figure
width = 2;
plot(userSet,mean(BDFD,2),'k','LineWidth',2),hold on
plot(userSet,mean(BDHP,2),'r','LineWidth',2),hold on
plot(userSet,mean(BDHP1,2),'r-.','LineWidth',2),hold on
plot(userSet,mean(BDHP2,2),'r--','LineWidth',2),hold on
grid on
xlabel('SNR (dB)')
ylabel('Sum Rate (bps/Hz)')