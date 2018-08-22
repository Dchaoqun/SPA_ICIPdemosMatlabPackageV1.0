% ICIP2017_BlockInpainting:
%
% A demo for simultaneous block inpainting and deblurring, using SPA.
%
% This is the code for getting Fig.1. (a), (b) and (d) in
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
% It assumes the original image was in [0 255] range (8 bits).
% Please, normalize the observation if you adapt the code to run with
% other bit depth images.
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
%{
% ADMM-UBC gives very good performance on simultaneous deblurring and block
% inpainting. Here we redo one experiment from these papers:
%
% [1]  Mariana S. C. Almeida and Mário A. T. Figueiredo; "Deconvolving images
%      with unknown boundaries using the alternating direction method of
%      multipliers",IEEE Transactions on Image Processing, vol. 22, No. 8,
%      pp. 3074-3086, 2013.
% [2]  Mariana S. C. Almeida and Mário A. T. Figueiredo; "Frame-based image
%      deblurring with unknown boundary conditions using the Alternating
%      Direction Method of Multipliers", IEEE Int. Conf. on Image Processing
%      Melbourne, Australia, September, 2013.
%
% The blurring kernel and noise level are the same as in Almeida and
% Figueiredo's experiment.
% The missing pixels are approximately in the same (non-reported) location
 % as in the their, and the size of the unknown block is exactly the same.
%}
im = double(imread('lena256.bmp'));
[Ny,Nx] = size(im);

x0 = im-min(im(:));
x0 = x0/max(x0(:));   % x0 in the range [0,1]

UB = 205;    % setting the upper bound
x0 = UB*x0;

BlurDim = 19 ;  % filter size (squared): 19x19 pixels
BSNR = 60;      % BSNR level (in dB)

% a uniform blurring kernel
ht = ones(BlurDim,BlurDim)/(BlurDim*BlurDim);
fsize = round( (BlurDim-1)/2 );         % fsize - size of each side of the filter
h_full = zeros(Ny,Nx);
h_full(Ny/2+1-fsize:Ny/2+1+fsize,Nx/2+1-fsize:Nx/2+1+fsize) = ht;
h_full = ifftshift(h_full);

%{
% In our ICIP 2017 paper, we compared our result with Almeida and Figueiredo's.
% So in order to make sure the degradation is identical to the one in Almeida
% and Figueiredo's TIP 2013 paper, we used their code to generate the
% blurring kernel in our experiment.
% The above implementation is an identical approach but created by us.
% One can download Almeida and Figueiredo's code from
% http://www.lx.it.pt/~mscla/ADMM_UBC.htm
% to run their original implementation as below.
%
% %Blurring filter:
n_filter = 1;   % filter type: 1 -  uniform
                %              2 -  out-of-focus
                %              3 -  linear motion at 135º
                %              4 -  Gaussian
ht2 = build_blur( BlurDim , n_filter );         % constructs the blurring filter (19x19)
h_full2 = ht2h( ht2 , size(x0,1) , size(x0,2));  % filter -> h (256x256 pixels)
%}

% In [1], they considered an observation with size 238x238, pixels affected
% by boundary conditions are discarded. In the 256x256 binary mask these
% pixels are inlcuded as zeros. Please see their papers for the details.
load lena_mask.mat

% Blurring degradation (with realistic BC):
y_blurred = ifft2(fft2(h_full).*fft2(x0));     % cyclic convolution

% Additive noise:
y_blurred_valid = y_blurred(1+fsize:end-fsize,1+fsize:end-fsize);
Py = var(y_blurred_valid(:));
sigma = sqrt(Py/10^(BSNR/10));
randn('seed',0);
noise = sigma*randn(size(y_blurred));

y_simu = (y_blurred+noise).*lena_mask;
y_obsv = y_simu(1+fsize:end-fsize,1+fsize:end-fsize);

figure(1),imagescTSc(y_obsv);title('obsevation')  %Fig.1.(a)


%% SPA for inpainting 
Le = 5;  % number of pixels for boundary extension

NDy = Ny + 2*Le;
NDx = Nx + 2*Le;

mask = zeros(NDy,NDx);
mask(Le+1:end-Le,Le+1:end-Le) = lena_mask;
figure(2); imagescTSc(mask); title('unknown pixels´ mask');

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

y = zeros(NDy,NDx);
y(Le+fsize+1:end-fsize-Le,Le+fsize+1:end-fsize-Le) = y_obsv;

P = 500;

z = SPA(y,ht,mask,PX_f,sigma,P);
figure(3),imagescTSc(z);title('result of SPA interpolation')   %Fig.1.(b)


%% restoration 
% using L2-relaxed-L0 deblurring Matlab Toolbox v2.1
% please check the toolbox for details

Niter = 10;

% Load parameters for L2-r-L0 deblurring with Niter=10:
% All three representatios are active.
% From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
% parameters jointly optimizing SSIM of the 24 test experiments.
load L2rL0_Niter10_parameters.mat

real_psf = padarray(fftshift(h_full),[Le,Le],0,'both');
Real_psf_f = real(fft2((ifftshift(real_psf))));

LB = 0;

x_hat = deblur_L2relaxedL0(z, rep, Real_psf_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
figure(4),imagescTSc(x_hat);title('result of deconvolution')  %Fig.1.(d)


%% final result (crop the external boundaries)
final_image = x_hat(Le+1:end-Le,Le+1:end-Le);
figure(5),imagescTSc(final_image);title('final result')

