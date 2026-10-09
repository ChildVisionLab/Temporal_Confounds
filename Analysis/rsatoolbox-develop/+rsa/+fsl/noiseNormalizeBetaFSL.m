function [u_hat,nu_hat] = noiseNormalizeBetaFSL(Beta,Res,Partition,nVols)

% function [u_hat,Sw_hat,resMS,beta_hat]=rsa_noiseNormalizeBeta(Y,SPM,varargin)
% Estimates beta coefficiencts beta_hat and residuals from raw time series Y
% Estimates the true activity patterns u_hat by applying noise normalization to beta_hat
% INPUT:
%    Y        raw timeseries, T by P
%    SPM:     SPM structure
% OPTIONS:
%   'normmode':   'overall': Does the multivariate noise normalisation overall (default)
%                 'runwise': Does the multivariate noise normalisation by run
%   'shrinkage':  Shrinkage coefficient. 
%                 0: No regularisation 
%                 1: Using only the diagonal - i.e. univariate noise normalisation 
%                 By default the shrinkage coeffcient is determined using
%                 the Ledoit-Wolf method. 
% OUTPUT:
%    u_hat:       estimated true activity patterns (beta_hat after multivariate noise normalization),
%    resMS:       residual mean-square - diagonal of the Var-cov matrix of the average beta weight, 1 by P  
%    Sw_hat:      overall voxel error variance-covariance matrix (PxP), before regularisation 
%    beta_hat:    estimated raw regression coefficients, K*R by P
%    shrinkage:   applied shrinkage facto 
%    trRR:        Trace of the squared residual spatial correlation matrix
%                 This number provides an estimate for the residual spatial
%                 correlation. For variance-estimation, the effective
%                 number of voxels are effVox = numVox^2/trRR 
% Alexander Walther, Joern Diedrichsen
% joern.diedrichsen@googlemail.com
% 10/2018: Updated scaling of prehwitened-betas, to take into account the
% variance of betas (bCov).

%leave nVols blank to run normal normalisation

[T,P] = size(Res);
[Q,P] = size(Beta);
Nrun= max(Partition);

if exist("nVols")

    partT = [];
    for n = 1:Nrun
        partT = [partT; (ones(1,nVols(n))*n)'];
    end
    %partT = kron(1:Nrun,ones(1,T/Nrun))'; %build own partT because run sizes are different, or don't normailse noise but check which distance measure is best to use
    partQ = Partition;
    u_hat = zeros(Q,P);
else
    partT = kron(1:Nrun,ones(1,T/Nrun))';
    partQ = Partition;
end

for i=1:Nrun
    idxT=partT==i;
    idxQ=partQ==i;
    Sw_hat(:,:,i)=rsa.stat.covdiag(Res(idxT,:));%%% regularize Sw_hat through optimal shrinkage
    u_hat(idxQ,:)=Beta(idxQ,:)*real(Sw_hat(:,:,i)^(-1/2));   %%% multivariate noise normalization %added real
end

%normalise concatenated runs
Sw_hat=rsa.stat.covdiag(Res);%%% regularize Sw_hat through optimal shrinkage
nu_hat=Beta*Sw_hat^(-1/2);   %%% multivariate noise normalization


% %normalize to avg covariance
% meanSw_hat = mean(Sw_hat,3);
% for i=1:Nrun
%     idxQ=partQ==i;
% u_hat(idxQ,:)=Beta(idxQ,:)*meanSw_hat^(-1/2);
% end

end
