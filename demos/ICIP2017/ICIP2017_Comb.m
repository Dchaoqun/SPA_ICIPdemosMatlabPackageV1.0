% ICIP2017_Comb:
%
% A demo for dealing with blurred objects having supports with arbitrary
% shapes using SPA.
% Here is the case (a real image taken by a professional camera) that in an
% image the Foreground is focused and Background is blurred.
% Camera Canon EOS 5D Mark II, 20 Mpix, CR2 to UINT16 TIFF. Focal length
% 50mm, exposure 1/2500s., ISO 100, f/3.5.
% Many thanks to Isabel Portilla for kindly providing the picture.
%
% This code recreates Fig.3. in
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
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


%% pre-processing and getting the object mask
original_image = imread('IMG_2373.tiff');
figure(1),imshow(original_image),title('observation')  %Fig.3.(a)

% get the object mask
im1 = sum(original_image,3);
im2 = im1<1.1e5;
im3 = imfill(im2,'holes');
im4 = bwmorph(im3,'dilate',[2 2]);
im5 = 1-imfill(im4,'holes');
mask_image = im5;
clear im1 im2 im3 im4 im5

% here we only process a grayscale image which is the average of RGB channel
im00 = double(mean(original_image,3));
[Ny,Nx] = size(im00);


%% define blur
D = 29.5;   % blurring kernel diameter, estimated
Le = 5;     % boundary extension
sigma = 4;  % noise level, estimated

[Real_psf_f_o,Real_psf_f,kernel,Ly,Lx,keveny,kevenx] = disk_kernel(D,Nx,Ny,Le);


%% SPA for estimating background image 
NDy = Ny + 2*Le;
NDx = Nx + 2*Le;

% blur mask to avoid artifacts
blur = [1 2 1]'*[1 2 1]/16;
mask_smooth = conv2(mask_image,blur,'valid');
temp = mask_image;
temp(2:end-1,2:end-1) = mask_smooth;
mask_smooth = temp;
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
y(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le) = im00(Ly+1:end-Ly,Lx+1:end-Lx).*mask(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le);

% PSD of the uncorrputed image
[xc, yc] = meshgrid(1:NDx,1:NDy); 
xc = xc - 1 - NDx/2;
yc = yc - 1 - NDy/2;
sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). We used it here despite the imaes are not 8bit. It works well because the deblurring part has very low contrast.
u = xc/NDx;
v = yc/NDy;
rho = 0.65;
PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
PX_f = fftshift(PX_f);

P=3;

z = SPA(y,kernel,mask,PX_f,sigma,P);
figure(3),imagescTSc(z);title('result of SPA interpolation') %Fig.3.(b)


%% restoration
% using L2-relaxed-L0 deblurring Matlab Toolbox v2.1
% please check the toolbox for details

Niter = 10;

% Load parameters for L2-r-L0 deblurring with Niter=10:
% All three representations are active.
% From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
% parameters jointly optimizing SSIM of the 24 test experiments.
load L2rL0_Niter10_parameters.mat

LB = 30000;   %lower bound for the background
UB = 54000;   %upper bound for the background

x_hat = deblur_L2relaxedL0(z, rep, Real_psf_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
figure(4),imagescTSc(x_hat);title('result of deconvolution')   %Fig.3.(c)


%% final result
final_image = x_hat(Le+1:end-Le,Le+1:end-Le).* mask_smooth + im00.*(1-mask_smooth);  % no sharp edge for the object
figure(5),imagescTSc(final_image);title('final result') %Fig.3.(d)
