function [G,dGdtheta] = pcm_calculateG_nlCorr_crossHemi_fixedCov(theta,M)
% function [G,dGdtheta] = pcm_calculateGnonlinCorr(theta,M)
% Calculate G and dGdtheta for nonlinear across-hemi correlation models. 
%
% M.Gc is the second moment of the contralater hemisphere roi during
% movement.
% M.Gc is flexibly scaled by r = (exp(2.*theta(1))-1)./(exp(2.*theta(1))+1)
% 
% Note that this scaling is equal for variances (Within-pattern
% correspondence) and covariances (structure of across-pattern
% relationships).
% This assumption is relaxed in pcm_calculateG_nlCorr_crossHemi_diffCov,
% where scaling is only applied to diagonal (variances) of M.Gc.
%
% saarbuckle 2020
error('this function needs to be completed');
z = theta(1);
r = (exp(2.*z)-1)./(exp(2.*z)+1); % Correlation is the inverse fisherz transform
idx = find(eye(5)==1);
G = M.Gc(:,:,1);
G(idx,idx) = M.Gc(:,:,2).*r; % Determine the predicted G-matrix using correlation paramter
dGdtheta(idx,idx,M.numGparams) = G./r; % Compute derivative for the correlation parameter 

