% ICIP2014:
%
% A demo code for running the experiments in J. Portilla ICIP2014
% paper:
% - J. Portilla, "Maximum likelihood extension for non-circulant
%   deconvolution", 2014 IEEE International Conference on  Image Processing
%   (ICIP), pp. 4276-4279, Oct 2014.
%
% This is the code for getting Fig.3 images in that paper.
% Also, by changing the parameters, one can redo all the experiments in J. Portilla 
% ICIP2014 paper and reproduce all ISNRs in Table 1.
%
% NOTES:
%   In the code for ICIP2014, we were using "ConDy10" (See E. Gil-Rodrigo et al., ICIP2011)
%   for restoration, which is an antecessor of the L2-r-L0 deblurring toolbox. Here in the code
%   we are using the updated one (L2-r-L0 Deblur Matlab Toolbox 2.1), so
%   the ISNRs are different from the ones in Table 1 in J. Portilla ICIP2014
%   paper (generally a bit higher).
%
%   SPA.m function was previously referred to as "MLE" or "MLI", for Maximum-
%   Likelihood Extension or Interpolation. For uniformity, we enforced 
%   here to be the common central tool for all reproduced experiments in all three ICIP
%   publications, the above referred one and:
%   - C. Dong, J. Portilla, "Maximum likelihood interpolation for
%     aliasing-aware image restoration", 2016 IEEE International Conference
%     on Image Processing (ICIP), pp. 564-568, Sept 2016.
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
% It is included in this ZIP file. It is an implementation of the method 
% described in:
% - J. Portilla, A. Tristan-Vega, I.W. Selesnick, "Efficient and robust
%   image restoration using multiple-feature L2-relaxed sparse analysis
%   priors", IEEE Transactions on Image Processing, vol. 24, no. 12,
%   pp. 5046-5059, Dec 2015.

% Code for ICIP publication:
% - J. Portilla, "Maximum likelihood extension for non-circulant
%   deconvolution", 2014 IEEE International Conference on  Image Processing
%   (ICIP), pp. 4276-4279, Oct 2014.
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
% In the paper we tested with combinations of 3 images with 8 different
% degradations and 2 different restoration methods.
%
% One can try different combinations of image, degradation and restoration
% method simply by changing the value of 'img','deg' and 'restor'.
%  img       - choose input image
%                'house' for using House image,
%                'cameraman' for using Cameraman image
%                'straw' for using Straw image
%  deg       - choose the degradation
%                deg = 1, PSF1 + low noise (sigma.^2 = 0.25)
%                deg = 2, PSF1 + medium noise (sigma.^2 = 2)
%                deg = 3, PSF2 + low noise (sigma.^2 = 0.31)
%                deg = 4, PSF2 + medium noise (sigma.^2 = 4)
%                deg = 5, PSF3 + low noise (sigma.^2 = 1)
%                deg = 6, PSF3 + medium noise (sigma.^2 = 4)
%                deg = 7, PSF4 + low noise (sigma.^2 = 0.25)
%                deg = 8, PSF4 + medium noise (sigma.^2 = 4)
%                (PSF1, out of focus; PSF2, uniform; PSF3, motion blur; PSF4, oblique)
%  restor    - choose restoration method
%                restor = 'Wiener' for using Wiener filtering
%                restor = 'L2rL0' for using L2-relaxed-L0 restoration
%
% The examples in Fig.3 are:
%    Straw image with degradation 7 and Wiener filtering
%    Parameters are: img = 'straw'; deg = 7; restor = 'Wiener';

img = 'straw';
deg = 7;
restor = 'Wiener';

% read image
if strcmp(img,'house')
    x = imread('house.png');
elseif strcmp(img,'cameraman')
    x = imread('cameraman.png');
elseif strcmp(img,'straw')
    x = imread('straw.tif');
end

% define kernels
% PSF1, Laplacian blur
i = meshgrid(-7:7, -7:7); j = i';
kernel = 1./(i.^2 + j.^2 + 1);
kernel = kernel/sum(kernel(:));
PSF1 = kernel;
% PSF2, uniform blur
kernel = ones(9)/81;
kernel = kernel/sum(kernel(:));
PSF2 = kernel;
% PSF3, motion blur
kernel = ones(9,1);
kernel = kernel/sum(kernel(:));
PSF3 = kernel;
% PSF4, oblique blur
kernel = [0 0 0 0 0 0 0;
    0 0 0 0 1 1 1;
    0 0 1 2 3 2 1;
    0 1 3 4 3 1 0;
    1 2 3 2 1 0 0;
    1 1 1 0 0 0 0;
    0 0 0 0 0 0 0];
kernel = kernel/sum(kernel(:));
PSF4 = kernel;

% define degradations
if deg==1,
    kernel = PSF1;
    sigma = sqrt(0.25);
elseif deg==2,
    kernel = PSF1;
    sigma = sqrt(2);
elseif deg==3,
    kernel = PSF2;
    sigma = 0.56; %sqrt(0.308);
elseif deg==4,
    kernel = PSF2;
    sigma = sqrt(4);
elseif deg==5,
    kernel = PSF3;
    sigma = sqrt(1);
elseif deg==6,
    kernel = PSF3;
    sigma = sqrt(4);
elseif deg==7,
    kernel = PSF4;
    sigma = sqrt(0.25);
elseif deg==8,
    kernel = PSF4;
    sigma = sqrt(4);
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

x = single(x);
X_f = fft2(x);
randn('seed',0);
Y_f = H_f.*X_f + fft2(sigma*randn(size(x)));
y_simu = real(ifft2(Y_f));

figure(1),imagescTSc(y_simu);title('simulation')

% cut off the invalid boundary pixels
y_valid = y_simu(Ly+1:end-Ly,Lx+1:end-Lx);


%% *****************************************************************************************
% ********     Method One: boundary mirror extension and edgetaper tiling (ET)     ********
% ********                     +     Restoration                                    ********
% ******************************************************************************************
fprintf(1,'# Method One: ET + Restoration  ...\n');


% #############       ET      ###############
Le = 8;   % number of pixels for boundary extension
NDy = Ny + 2*Le;
NDx = Nx + 2*Le;

z_ET = bound_extension(y_valid,Ly+Le,Lx+Le,'mirror'); % mirror extension
PSF = fspecial('gaussian', 40, 6);
z_ET  = edgetaper(z_ET, PSF);

figure(2),imagescTSc(z_ET);title('result of ET')


% #############     RESTORATION    #############
% psf in extended size
h = zeros(NDy,NDx);
h(NDy/2+1-Ly+keveny:NDy/2+1+Ly,NDx/2+1-Lx+kevenx:NDx/2+1+Lx) = kernel;
D_f= fft2(fftshift(h));

if strcmp(restor,'Wiener')
    
    % PSD of the uncorrputed image
    [xc, yc] = meshgrid(1:NDx,1:NDy); 
    xc = xc - 1 - NDx/2;
    yc = yc - 1 - NDy/2;
    sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time), for [0 255] range images.
    u = xc/NDx;
    v = yc/NDy;
    rho = 0.65;
    PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
    PX_f = fftshift(PX_f);
    
    PY_f = abs(D_f).^2.*PX_f + sigma^2;
    Hwie_f = PX_f.*conj(D_f)./PY_f;
    
    % ET
    Z_f_ET = fft2(z_ET);
    x_hat_ET = real(ifft2(Z_f_ET.*Hwie_f));
    x_hat_ET = max(min(x_hat_ET, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    fprintf(1,'Info:  Deblurring ET result  ...\n');
    x_hat_ET = deblur_L2relaxedL0(z_ET, rep, D_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

x_hat_valid_ET = x_hat_ET(Le+1:end-Le,Le+1:end-Le); % cut off the extended pixls

figure(3),imagescTSc(x_hat_valid_ET);title('result of deconvolution on ET')


%% ******************************************************************************************
% *********************           Method Two: SPA + Restoration           *******************
% *******************************************************************************************
fprintf(1,'# Method Two: SPA + Restoration  ...\n');


% #############      SPA     #############
Le = 8;   % number of pixels for boundary extension
NDy = Ny + 2*Le;
NDx = Nx + 2*Le;

% mask for SPA
mask = zeros(NDy,NDx);
mask(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le) = 1;

y = zeros(NDy,NDx);
y(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le) = y_valid.*mask(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le);

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

z_SPA = SPA(y,kernel,mask,PX_f,sigma,P);

figure(4),imagescTSc(z_SPA);title('result of SPA interpolation')


% #############     RESTORATION    #############
% psf in extended size
h = zeros(NDy,NDx);
h(NDy/2+1-Ly+keveny:NDy/2+1+Ly,NDx/2+1-Lx+kevenx:NDx/2+1+Lx) = kernel; 
D_f= fft2(fftshift(h));

if strcmp(restor,'Wiener')
    
    PY_f = abs(D_f).^2.*PX_f + sigma^2;
    Hwie_f = PX_f.*conj(D_f)./PY_f;
    
    % SPA
    Z_f_SPA = fft2(z_SPA);
    x_hat_SPA = real(ifft2(Z_f_SPA.*Hwie_f));
    x_hat_SPA = max(min(x_hat_SPA, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    fprintf(1,'Info:  Deblurring SPA result  ...\n');
    x_hat_SPA = deblur_L2relaxedL0(z_SPA, rep, D_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

x_hat_valid_SPA = x_hat_SPA(Le+1:end-Le,Le+1:end-Le); % cut off the extended pixels

figure(5),imagescTSc(x_hat_valid_SPA);title('result of deconvolution on SPA')


%% Oracle
% directly deblurring on the simulation
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

    PY_f = abs(H_f).^2.*PX_f + sigma^2;
    Hwie_f = PX_f.*conj(H_f)./PY_f;
    
    Y_SIMU_f = fft2(y_simu);
    x_hat_Ora = real(ifft2(Y_SIMU_f.*Hwie_f));
    x_hat_Ora = max(min(x_hat_Ora, 255), 0);  % back to range [0,255]
    
elseif strcmp(restor,'L2rL0')
    
    % L2-relaxed-L0 deblurring Matlab Toolbox v2.1
    % please check the toolbox for details
    
    Niter = 10;
    
    % Load parameters for L2-r-L0 deblurring with Niter=10:
    % All three representations are active.
    % From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
    % parameters jointly optimizing SSIM of the 24 test experiments.
    load L2rL0_Niter10_parameters.mat
    
    LB = 0;
    UB = 255;
    
    fprintf(1,'Info:  Deblurring on simulation to get Oracle result  ...\n');
    x_hat_Ora = deblur_L2relaxedL0(y_simu, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
    
end

figure(6),imagescTSc(x_hat_Ora);title('Oracle')


%% ISNR
% only compute the ISNR in the valid area
psnr_simu = psnr(y_simu(Ly+1:end-Ly,Lx+1:end-Lx), x(Ly+1:end-Ly,Lx+1:end-Lx));

psnr_ET = psnr(x_hat_ET(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le), x(Ly+1:end-Ly,Lx+1:end-Lx));
isnr_ET = psnr_ET - psnr_simu

psnr_SPA = psnr(x_hat_SPA(Le+Ly+1:end-Ly-Le,Le+Lx+1:end-Lx-Le), x(Ly+1:end-Ly,Lx+1:end-Lx));
isnr_SPA = psnr_SPA - psnr_simu

psnr_Ora = psnr(x_hat_Ora(Ly+1:end-Ly,Lx+1:end-Lx),x(Ly+1:end-Ly,Lx+1:end-Lx));
isnr_Ora = psnr_Ora - psnr_simu

%{
%% crop the result image to get the area shown in Fig.3
% This part should only be used for the corresponding parameter settings
% as following:
%    Straw image with degradation 7 and Wiener filtering
%    Paramters are: img = 'straw'; deg = 7; rest = 'Wiener';

% cropped simulation image area in Fig.3 left image
crop_im_simu = y_simu(4:123,4:123);
figure(7),imagescTSc(crop_im_simu);title(strcat('ICIP2014 Fig.3 left image'))

% cropped result image area in Fig.3 right image
crop_im_ET = x_hat_valid_ET(4:123,4:123);
figure(8),imagescTSc(crop_im_ET);title(strcat('ICIP2014 Fig.3 middle image'))

% cropped result image area in Fig.3 right image
crop_im_SPA = x_hat_valid_SPA(4:123,4:123);
figure(9),imagescTSc(crop_im_SPA);title(strcat('ICIP2014 Fig.3 right image'))
%}
