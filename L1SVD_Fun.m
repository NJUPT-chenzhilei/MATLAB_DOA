clear all
close all
%场景的设置
Omega = [1 2 3 4 5 6 7]';  %阵元数
M = length(Omega);
L = 200;    % 快拍数
DOA = [-15 5]'; % 入射角角度范围
K = length(DOA);%信源数
SNR = 10;   %信噪比
Ang1 = -90;
Ang2 = 90;
Angle = (Ang1:2:Ang2)';
N = length(Angle);
%数据集构造
A = exp(1i*pi*[0:M-1]'*sin(DOA'*pi/180));   %时延
S = (randn(K,L)+1i*randn(K,L))/sqrt(2);
S = diag(sqrt(10^(SNR/10)./diag(1/L*(S*S'))))*S;
Noise = 1/sqrt(2)*(randn(M,L)+1i*randn(M,L));
X = A * S +Noise;
%Aall的构造，稀疏恢复和非参数类解法不同的地方：在于它事先会不基于数据集构造一个流行矩阵进行X的稀疏恢复进行匹配
for m = 1:M
    for n = 1:N
        Aall(m,n) = exp(1i*pi*(m-(M+1)/2)*sin(Angle(n)*pi/180));
    end
end


%l1-svd算法
lambda_L1SVD = 0.625;
gamma_L1SVD = L1SVD_Fun(X,K,Aall,lambda_L1SVD);

%画功率谱
h = plot(Angle, 10*log10(gamma_L1SVD/max(gamma_L1SVD)));%参考MUSIC算法里功率谱的运用；plot(x,y)
hold on
plot(DOA,[0 0],'o'),
set(h,'Linewidth',2) %设置所有线宽为2
xlabel('angle (degree)')
ylabel('magnitude (dB)')
axis([-90 90 -60 0])  %axis[Xmin Xmax Ymin Ymax]
set(gca, 'XTick',[-90:30:90])%设置网格的显示格式，gca获取当前figure的句柄
grid on  

function SigEst_L1SVD = L1SVD_Fun(y,K,Aall0,lambda)
% y是接受矩阵，K是信源个数，Aall0是阵列流型矩阵，lambda是
[M,N] = size(y);
if N > K
    [U0, S0, V0] = svd(y/sqrt(N));
    d0 = diag(S0);
    K1 = min(M,K);
    Y_s = U0(:, 1:K1)*diag(d0(1:K1));
else
    Y_s = y;
end
SigEst_L1SVD = l1_soc_mljsq_joint(Y_s, Aall0, lambda);

SigEst_L1SVD = abs(SigEst_L1SVD);
