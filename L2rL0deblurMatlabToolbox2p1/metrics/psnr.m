% [PSNR, error] = psnr(X,Y);
% Input Parameters:
%    X: Original
%    Y: Estimation
% Output Parameters
%    PSNR: Peak-Signal to noise ratio
%    error: Error of Estimation
%
%  Jose A. Guerrero-Colon (14/07/2004)
%  Universidad de Granada (Spain)
function [PSNR, error] = psnr(X,Y);

%% Compute the error
error = Y - X;
%%Mean Square Error
mse = mean2(error.^2);
PSNR = 10 * log10(255^2 / mse);

