% ICIP2017_Pirate:
%
% A demo for deblurring objects in images having supports with arbitrary
% shape, by using SPA.
% Here is the case (a realistic simulation) that in an image the Foreground
% is focused and Background is blurred (out of focus).
%
% This code provides Fig.2. (a), (b) (d) and (f) in 
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
%  It assumes the original image was in [0 255] range (8 bits).
%  Please, normalize the observation if you adapt this code to run with
%  other bit depth images.
%
% Requires the L2-r-L0 Deblur Matlab Toolbox 2.1:
% https://www.researchgate.net/publication/325903515_L2rL0deblurMatlabToolbox2p1
% It is included in this ZIP file. The toolbox is an implementation of the method 
% described in:
% - J. Portilla, A. Tristan-Vega, I.W. Selesnick, "Efficient and robust
%   image restoration using multiple-feature L2-relaxed sparse analysis
%   priors", IEEE Transactions on Image Processing, vol. 24, no. 12,
%   pp. 5046-5059, Dec 2015.

% Code for ICIP publication:
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
% Version 1.0, August 2018
%
% Copyright (c) 2018
% Chaoqun Dong <cdongae@connect.ust.hk>
% Javier Portilla <javier.portilla@csic.es>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.


%%
clear all;
close all;

addpath(genpath('../../data'))
addpath(genpath('../../core'))
addpath(genpath('../../functions'))
addpath(genpath('../../L2rL0deblurMatlabToolbox2p1'))


%% simulation    
background_image = double(imread('pirate512.tif'));
foreground_image = double(imread('hand.png'));
[Ny,Nx] = size(background_image);

% mask for foreground in focus object
% hand_mask.png is processed from hand.png by thresholding and 
% morphological operations.
% focused foreground pixels are zeros
mask_image = imread('hand_mask.png');   
mask_image = im2double(mask_image);    

% define blur
sigma = 0.3;   %low noise level
D = 9;         % diameter of the blurring kernel
Le = 8;        % boundary extension parameter
[Real_psf_f_o,Real_psf_f,kernel,Ly,Lx,keveny,kevenx] = disk_kernel(D,Nx,Ny,Le);   % psf

% blur mask to avoid artifacts
blur = [1 2 1]'*[1 2 1]/16; 
mask_smooth = conv2(mask_image,blur,'valid');
temp = mask_image;
temp(2:end-1,2:end-1) = mask_smooth;
mask_smooth = temp;

B_f = fft2(background_image);
randn('seed',0);
BY_f = Real_psf_f_o.*B_f + fft2(sigma*randn(Ny,Nx)); % observation model for noisy convolved image
b_simu = real(ifft2(BY_f));
y_simu = b_simu.* mask_smooth + foreground_image.*(1-mask_smooth);

y_obsv = y_simu(Ly+1:end-Ly,Lx+1:end-Lx); % cut off invalid boundary pixels

figure(1),imagescTSc(y_obsv);title('observation');   % Fig.2. (a) 


%% SPA for estimating background image
NDy = Ny + 2*Le;
NDx = Nx + 2*Le;

mask_image_back = 1-(mask_smooth<1);   % sharp background mask
% boundary handling mask
mask_boun = zeros(NDy,NDx);
mask_boun(Le+Ly-keveny+1:end-Ly-Le,Le+Lx-kevenx+1:end-Lx-Le) = 1;

% object mask
mask_back = zeros(NDy,NDx);
mask_back(Le+1:end-Le,Le+1:end-Le) = mask_image_back;

% generic mask 
mask = mask_back & mask_boun;   %knows are ones, unknowns are zeros
figure(2); imagescTSc(mask); title('unknown pixels´ mask');

y = zeros(NDy,NDx);
y(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le) = y_obsv.*mask(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le);

% PSD of the uncorrputed image
[xc, yc] = meshgrid(1:NDx,1:NDy); 
xc = xc - 1 - NDx/2;
yc = yc - 1 - NDy/2;
sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
u = xc/NDx;
v = yc/NDy;
rho = 0.65;
PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
PX_f = fftshift(PX_f);

P=4;

z = SPA(y,kernel,mask,PX_f,sigma,P);
figure(3),imagescTSc(z);title('result of SPA interpolation')   %Fig.2.(b)


%% L2-relaxed-L0 deblurring Matlab Toolbox v2.1
% please check the toolbox for details

Niter = 10;

% Load parameters for L2-r-L0 deblurring with Niter=10:
% All three representations are active.
% From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
% parameters jointly optimizing SSIM of the 24 test experiments.
load L2rL0_Niter10_parameters.mat

x_hat = deblur_L2relaxedL0(z, rep, Real_psf_f, sigma, Niter, alpha_v, sig2r0_v, beta);
figure(4),imagescTSc(x_hat);title('result of deconvolution')  %Fig.2.(d)


%% final result
final_image= x_hat(Le+1:end-Le,Le+1:end-Le).* mask_smooth + y_simu.*(1-mask_smooth);  % no sharp edge
figure(5),imagescTSc(final_image);title('final result')   %Fig.2.(f)


%% SNR only for background
x_back = double(background_image).*mask_image_back;
final_image_back = final_image.*mask_image_back;
SNR = 10*log10(mean2(x_back(Ly+1:end-Ly,Lx+1:end-Lx).^2)/mean2((final_image_back(Ly+1:end-Ly,Lx+1:end-Lx) - x_back(Ly+1:end-Ly,Lx+1:end-Lx)).^2)) 

