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
