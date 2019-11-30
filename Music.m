clear all
close all
derad = pi/180;        % deg -> rad
radeg = 180/pi;
twpi = 2*pi;
M = 8;               % 阵列数量
K = 3;              % number of DOA
theta = [-10 30 60];     % 角度
SNR = 10;               % input SNR (dB)
L = 1000;
A = exp(1i*pi*[0:M-1]'*sin(theta*derad));%流型矩阵
S = (randn(K,L)+1i*randn(K,L))/sqrt(2);  %randn()均值为0方差为1的正态分布
S = diag(sqrt(10^(SNR/10)./diag(1/L*(S*S'))))*S;
Noise = 1/sqrt(2)*(randn(M,L)+1i*randn(M,L));
X = A*S + Noise;

[M,N] = size(X);
Rx = (X*X')/N;
[U,S,V] = svd(Rx);

Angle = (-90:1:90);
AngScope = (-90:1:90)*derad;
ASLen = length(AngScope);

P = zeros(1,ASLen);
for ASIdx = 1:ASLen
    Atest  = exp(1i*pi*[0:M-1]'*sin(AngScope(ASIdx)));
    P(ASIdx) = 1/abs(Atest'*U(:,K+1:M)*U(:,K+1:M)'*Atest);
end

P = abs(P);
Pmax = max(P);
P = 10*log10(P/Pmax);
h = plot(Angle,P);
set(h,'Linewidth',2)
xlabel('angle (degree)')
ylabel('magnitude(dB)')
axis([-90 90 -60 0])
set(gca, 'XTick',[-90:30:90])
grid on  
