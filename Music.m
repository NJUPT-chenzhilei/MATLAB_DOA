function [P,DOAE] MUSICFun(X,K,M,Ang1,Ang2,step)   
%该函数返回功率谱P和DOA，其中X是接收数据，K是信源个数，M是天线个数，Ang1和Ang2分别是角度的范围，step是进行步长

[M,N] = size(X);%M是天线个数，N是快拍数
Rx = 1/N*(N*N');%得到协方差矩阵的估计值
[U,S,V] = svd(Rx);%svd()特征值分解，Rx = U*S*V',mean S 是一个与M*N的对角矩阵，对角元素为特征值，U为大小特征向量组成

AngScope = (Ang1:step:Ang2)*pi/180; %AngScope是一个类似于python中列表的东西，里面是要检测的角度值
ASLen = length(AngScope);%AngScope中含有多少元素

P = zeros(1,ALSen);
for AXIdx = 1:ASLen
    Atest = exp(li*pi*[0,M-1]'*sin(AngScope(ASIdx)));
    P(ASIdx) = 1/abs(Atest'*U(:,K-1:M)*U(:,K-1:M)*Atest);
end

DOAE=[];
DOAE = PeakFinding(P,K,Ang1,Ang2);%PeakFinding为matlab中自带的函数，进行谱峰搜索
