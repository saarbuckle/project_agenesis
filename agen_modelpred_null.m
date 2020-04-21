function [G,dGdtheta] = agen_modelpred_null(theta,M)
% nonlinear null model
G        = eye(M.numCond)*exp(theta(1)); 
dGdtheta = G;