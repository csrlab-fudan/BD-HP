function [H,Ar]=channel_generation_ula(phi,theta,Mr,Mt,L,Nf)
Ar = zeros(Mr,L);
deltaPhi = 2*pi*rand(1,L);
deltaTheta = 2*pi*rand(1,L);
H = zeros(Mr,Mt,Nf);
for kk = 1:Nf
    for jj = 1:L
        a = exp(1j*pi*sin(phi+deltaPhi(jj))*(0:Mr-1)');
        b = exp(1j*pi*sin(theta+deltaTheta(jj))*(0:Mt-1)');
        temp =  (randn+1i*randn)/sqrt(2*L)*a*b'*exp(-1i*2*pi/Nf*(kk-1)*(jj-1));
        Ar(:,jj) = a;
        H(:,:,kk) = H(:,:,kk)+temp;
    end
end

