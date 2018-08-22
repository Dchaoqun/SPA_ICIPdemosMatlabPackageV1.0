% ICIP2016:
%
% A demo for interpolating and restoring a blurred observation by handling
% aliasing using SPA (Spectral Pre-Adaptation), in which case the known subset
% is a uniformly sampled version of the complete set of pixels.
%
% Using this code one directly gets Fig.2 images in
% - C. Dong, J. Portilla, "Maximum likelihood interpolation for
%   aliasing-aware image restoration", 2016 IEEE International Conference
%   on Image Processing (ICIP), pp. 564-568, Sept 2016.
% And, by simply changing the parameters, one can redo all the experiments
% in the above referred paper paper and reproduce all ISNRs in Table 1.
%
% NOTES:
%   SPA.m function was previously referred to as "MLE" or "MLI", for Maximum-
%   Likelihood Extension or Interpolation. For uniformity, we enforced 
%   here to be the common central tool for all reproduced experiments in all
%   three ICIP publications - the above referred one and:
%   - J. Portilla, "Maximum likelihood extension for non-circulant
%     deconvolution", 2014 IEEE International Conference on  Image Processing
%     (ICIP), pp. 4276-4279, Oct 2014.
%   - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%     arbitrary-shape-support image restoration", 2017 IEEE International
%     Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
%   It assumes the original image was in [0 255] range (8 bits).
%   Please, normalize the observation if you adapt the code to run with
%   other bit depth images.
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
% - C. Dong, J. Portilla, "Maximum likelihood interpolation for
%   aliasing-aware image restoration", 2016 IEEE International Conference
%   on Image Processing (ICIP), pp. 564-568, Sept 2016.
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


%% set up
% In ICIP 2016 paper we tested SPA for interpolating and restoring a blurred
% and aliased observation with combinations of 2 different images, 3 different
% blurring kernels and 3 different noise levels. We also tested with two
% different restoration methods, Wiener filter and L2-r-L0.
%
% One can try different combinations of image, kernel and noise level simply
% by changing the value of 'img', 'psf', 'nsi' and'restor'.
%     img      - choose input image.
%                'pirate' for using Pirate image,
%                'barbara' for using Barbara image.
%     psf      - choose blurring kernel
%                psf = 1, Laplacian
%                psf = 2, uniform
%                psf = 3, oblique
%     nis      - choose noise level
%                nis = 1, low noise level (sigma.^2 = 0.09)
%                nis = 2, medium noise level (sigma.^2 = 0.25)
%                nis = 3, high noise level (sigma.^2 = 2.25)
%    restor    - choose restoration method
%                restor = 'Wiener' for using Wiener filter
%                restor = 'L2rL0' for using L2-relaxed-L0 restoration
%
% The results in Fig.2 in ICIP2016 are:
% 1. Barbara image with PSF3 and low noise level (sigma.^2=0.09)
%    img = 'barbara'; psf = 3; nis = 1; restor= 'L2rL0';
% 2. Pirate image with PSF2 and high level (sigma.^2=2.25):
%    img = 'pirate'; psf = 2; nis = 3; restor= 'L2rL0';

img = 'barbara';  %'pirate' for using Pirate image, 'barbara' for using Barbara image
psf = 3;   % choose kernel
nis = 1;   % choose noise level
restor= 'L2rL0'; % choose between Wiener and L2rL0

% read image
if strcmp(img,'pirate')
    x = imread('pirate512.tif');
elseif strcmp(img,'barbara')
    x = imread('barbara512.png');
end

% define kernels
if psf==1
    % PSF1, Laplacian blur
    i = meshgrid(-7:7, -7:7); j = i';
    kernel = 1./(i.^2 + j.^2 + 1);
    kernel = kernel/sum(kernel(:));
elseif psf==2
    % PSF2, uniform blur
    kernel = ones(9)/81;
    kernel = kernel/sum(kernel(:));
elseif psf==3
    % PSF3, oblique blur
    kernel = [0 0 0 0 0 0 0;
        0 0 0 0 1 1 1;
        0 0 1 2 3 2 1;
        0 1 3 4 3 1 0;
        1 2 3 2 1 0 0;
        1 1 1 0 0 0 0;
        0 0 0 0 0 0 0];
    kernel = kernel/sum(kernel(:));
end

% define noise level
if nis==1
    sigma = 0.3;  %low noise level
elseif nis==2
    sigma = 0.5;  %medium noise level
elseif nis==3
    sigma = 1.5;  % high noise level
end


%% simulation
[Ny,Nx] = size(x);

% make sure the kernel is normalized
kernel = kernel/sum(sum(kernel));
[Ky,Kx] = size(kernel);
% half kernel size
Ly = ceil((Ky-1)/2);
Lx = ceil((Kx-1)/2);
keveny = double(rem(Ky,2)==0); % odd/even kernel size flag
kevenx = double(rem(Kx,2)==0);
% Fourier transform
h = zeros(Ny,Nx);
h(Ny/2+1-Ly+keveny:Ny/2+1+Ly,Nx/2+1-Lx+kevenx:Nx/2+1+Lx) = kernel; 
H_f= fft2(fftshift(h));

X_f = fft2(single(x));
randn('seed',0);
Y_f = H_f.*X_f + fft2(sigma*randn(size(x)));
y_simu = real(ifft2(Y_f));

% ATTENTION:
% Subsampling process is done by taking one out of every two pixels
% along each dimension. 
% We also assumed here the input image has even size.
y_obsv = y_simu(1:2:end,1:2:end);

figure(1),imagescTSc(y_obsv);title('observation')

% sub-sampled kernel
% subsample by taking one out of two pixels
kernel_sub = kernel(1+mod(Ly,2):2:end,1+mod(Lx,2):2:end);
kernel_sub = kernel_sub/sum(sum(kernel_sub));
[Kys,Kxs] = size(kernel_sub);
% half kernel size
Lys = ceil((Kys-1)/2);
Lxs = ceil((Kxs-1)/2);
% Fourier transform
h = zeros(Ny/2,Nx/2); % Nx, Ny are assumed to be even!
h(Ny/4+1-Lys+keveny:Ny/4+1+Lys,Nx/4+1-Lxs+kevenx:Nx/4+1+Lxs) = kernel_sub; 
H_sub_f= fft2(fftshift(h));


%% **************************************************************************************
% ***********        Method One: Restoration + Interpolation (R+I)        ***************
% ***************************************************************************************
fprintf(1,'# Method One: Restoration + Interpolation (R+I)   ...\n');


% #############    RESTORATION   #############
if strcmp(restor,'Wiener')
    
    Nsx = Nx/2;
    Nsy = Ny/2;
    % PSD of the uncorrputed image
    [xc, yc] = meshgrid(1:Nsx,1:Nsy); 
    xc = xc - 1 - Nsx/2;
    yc = yc - 1 - Nsy/2;
    sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
    u = xc/Nsx;
    v = yc/Nsy;
    rho = 0.65;
    PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
    PX_f = fftshift(PX_f);
    
    % blurred noisy image PSD
    PY_f = abs(H_sub_f).^2.*PX_f + sigma^2;
    % Wiener filter
    Hwie_f = PX_f.*conj(H_sub_f)./PY_f;
    
    Y_OBSV_f = fft2(y_obsv);
    y_obsv_hat = real(ifft2(Y_OBSV_f.*Hwie_f));
    y_obsv_hat = max(min(y_obsv_hat, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, CoNy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    y_obsv_hat = deblur_L2relaxedL0(y_obsv, rep, H_sub_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end


% ############     INTERPOLATION   #############
[Nsy,Nsx] = size(y_obsv_hat);
[xsc, ysc] = meshgrid(1:Nsx,1:Nsy);
[xc, yc] = meshgrid(1:0.5:Nsx+1,1:0.5:Nsy+1);
xc = xc(1:end-1,1:end-1);
yc = yc(1:end-1,1:end-1);

x_hat_RI = interp2(xsc,ysc,double(y_obsv_hat),xc,yc,'spline');

% clip in range [0,255]
x_hat_RI = max(min(x_hat_RI,255),0);

figure(2),imagescTSc(x_hat_RI);title('result of R + I')


%% **************************************************************************************
% **************     Method Two: Interpolation + Restoration (I+R)         **************
% **********************************************************************%****************
fprintf(1,'# Method Two: Interpolation + Restoration (I+R)    ...\n');


% #############    INTERPOLATION   #############
[Nsy,Nsx] = size(y_obsv_hat);
[xsc, ysc] = meshgrid(1:Nsx,1:Nsy);
[xc, yc] = meshgrid(1:0.5:Nsx+1,1:0.5:Nsy+1);
xc = xc(1:end-1,1:end-1);
yc = yc(1:end-1,1:end-1);

z_IR = interp2(xsc,ysc,double(y_obsv),xc,yc,'spline');


% #############     RESTORATION    ##############
if strcmp(restor,'Wiener')
    
    % PSD of the uncorrputed image
    [xc, yc] = meshgrid(1:Nx,1:Ny); 
    xc = xc - 1 - Nx/2;
    yc = yc - 1 - Ny/2;
    sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
    u = xc/Nx;
    v = yc/Ny;
    rho = 0.65;
    PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
    PX_f = fftshift(PX_f);
    
    % blurred noisy image PSD
    PY_f = abs(H_f).^2.*PX_f + sigma^2;
    % Wiener filter
    Hwie_f = PX_f.*conj(H_f)./PY_f;
    
    Z_IR_f = fft2(z_IR);
    x_hat_IR = real(ifft2(Z_IR_f.*Hwie_f));
    x_hat_IR = max(min(x_hat_IR, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, CoNy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    x_hat_IR = deblur_L2relaxedL0(z_IR, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

figure(3),imagescTSc(x_hat_IR);title('result of I + R')



%% **************************************************************************************
% *****************            Method Three: SPA + Restoration          *****************
% **********************************************************************%****************
fprintf(1,'# Method Three: SPA + Restoration     ...\n');


% #############       SPA      ##############
% mask for SPA
mask = zeros(Ny,Nx);
mask(1:2:end,1:2:end) = 1;

% PSD of the uncorrputed image
[xc, yc] = meshgrid(1:Nx,1:Ny); 
xc = xc - 1 - Nx/2;
yc = yc - 1 - Ny/2;
sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
u = xc/Nx;
v = yc/Ny;
rho = 0.65;
PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
PX_f = fftshift(PX_f);

y = zeros(Ny,Nx);
y(1:2:end,1:2:end) = y_obsv.*mask(1:2:end,1:2:end);

P=2;

z_SPA = SPA(y,kernel,mask,PX_f,sigma,P);


% #############     RESTORATION    ##############
if strcmp(restor,'Wiener')
    
    % blurred noisy image PSD
    PY_f = abs(H_f).^2.*PX_f + sigma^2;
    % Wiener filter
    Hwie_f = PX_f.*conj(H_f)./PY_f;
    
    Z_f_SPA = fft2(z_SPA);
    x_hat_SPA = real(ifft2(Z_f_SPA.*Hwie_f));
    x_hat_SPA = max(min(x_hat_SPA, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, CoNy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    x_hat_SPA = deblur_L2relaxedL0(z_SPA, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

figure(4),imagescTSc(x_hat_SPA);title('result of SPA + Restoration')


%% Oracle
% restoration directly on simulation
fprintf(1,'Processing for Oracle  ...\n');


if strcmp(restor,'Wiener')
    % PSD of the uncorrputed image
    [xc, yc] = meshgrid(1:Nx,1:Ny); 
    xc = xc - 1 - Nx/2;
    yc = yc - 1 - Ny/2;
    sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
    u = xc/Nx;
    v = yc/Ny;
    rho = 0.65;
    PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
    PX_f = fftshift(PX_f);
    
    % blurred noisy image PSD
    PY_f = abs(H_f).^2.*PX_f + sigma^2;
    % Wiener filter
    Hwie_f = PX_f.*conj(H_f)./PY_f;
    
    Y_SIMU_f = fft2(y_simu);
    y_simu_hat = real(ifft2(Y_SIMU_f.*Hwie_f));
    y_simu_hat = max(min(y_simu_hat, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, CoNy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    y_simu_hat = deblur_L2relaxedL0(y_simu, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

figure(5),imagescTSc(y_simu_hat);title('Oracle')


%% ISNR
% observation repeated(for calculating ISNR)
y_repeat = zeros(Ny,Nx);
y_repeat(1:2:end,1:2:end) = y_obsv;
y_repeat(1:2:end,2:2:end) = y_obsv;
y_repeat(2:2:end,1:2:end) = y_obsv;
y_repeat(2:2:end,2:2:end) = y_obsv;

psnr_obsvRe = psnr(y_repeat,single(x)); % psnr of observation repeated

psnr_RI = psnr(x_hat_RI,single(x));    % psnr of R+I
isnr_RI = psnr_RI - psnr_obsvRe

psnr_IR = psnr(x_hat_IR,single(x));    % psnr of I+R
isnr_IR = psnr_IR - psnr_obsvRe

psnr_SPA = psnr(x_hat_SPA,single(x));    % psnr of SPA+Restoration
isnr_SPA = psnr_SPA - psnr_obsvRe

psnr_Ora = psnr(y_simu_hat,single(x));    % psnr of Oracle
isnr_Ora = psnr_Ora - psnr_obsvRe



%{
%% crop the result image to get the area shown in Fig.2 top row and column 3
% This part should only be used for the corresponding parameter settings
% as following:
% 1. Barbara image with PSF3 and low noise level (sigma=0.3)
%    img = 'barbara'; psf = 3; nis = 1; restor= 'L2rL0';
% 2. Pirate image with PSF2 and high level (sigma=1.5):
%    img = 'pirate'; psf = 2; nis = 3; restor= 'L2rL0';

if strcmp(img,'pirate')
    % cropped observation image area in Fig.2 top row
    crop_im_simu = y_simu(17:17+255,51:51+255);
    crop_im_obsv = crop_im_simu(1:2:end,1:2:end);
    figure(6),imagescTSc(crop_im_obsv);title(strcat('ICIP2016 Fig.2 top row',13,img))
    
    % cropped result image area in Fig.2 column 1
    crop_im_RI = x_hat_RI(17:17+255,51:51+255);
    figure(7),imagescTSc(crop_im_RI);title(strcat('ICIP2016 Fig.2 column 1',13,img))
    
    % cropped result image area in Fig.2 column 2
    crop_im_IR = x_hat_IR(17:17+255,51:51+255);
    figure(8),imagescTSc(crop_im_IR);title(strcat('ICIP2016 Fig.2 column 2',13,img))
    
    % cropped result image area in Fig.2 column 3
    crop_im_SPA = x_hat_SPA(17:17+255,51:51+255);
    figure(9),imagescTSc(crop_im_SPA);title(strcat('ICIP2016 Fig.2 column 3',13,img))
    
    % cropped result image area in Fig.2 column 4
    crop_im_Orac = y_simu_hat(17:17+255,51:51+255);
    figure(10),imagescTSc(crop_im_Orac);title(strcat('ICIP2016 Fig.2 column 4',13,img))
    
elseif strcmp(img,'barbara')
    % cropped observation image area in Fig.2 top row
    crop_im_simu = y_simu(1:256,257:end);
    crop_im_obsv = crop_im_simu(1:2:end,1:2:end);
    figure(6),imagescTSc(crop_im_obsv);title(strcat('ICIP2016 Fig.2 top row',13,img))
    
    % cropped result image area in Fig.2 column 1
    crop_im_RI = x_hat_RI(1:256,257:end);
    figure(7),imagescTSc(crop_im_RI);title(strcat('ICIP2016 Fig.2 column 1',13,img))
    
    % cropped result image area in Fig.2 column 2
    crop_im_IR = x_hat_IR(1:256,257:end);
    figure(8),imagescTSc(crop_im_IR);title(strcat('ICIP2016 Fig.2 column 2',13,img))
    
    % cropped result image area in Fig.2 column 3
    crop_im_SPA = x_hat_SPA(1:256,257:end);
    figure(9),imagescTSc(crop_im_SPA);title(strcat('ICIP2016 Fig.2 column 3',13,img))
    
    % cropped result image area in Fig.2 column 4
    crop_im_Orac = y_simu_hat(1:256,257:end);
    figure(10),imagescTSc(crop_im_Orac);title(strcat('ICIP2016 Fig.2 column 4',13,img))
end
%}
