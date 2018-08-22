% ICIP2017_Pen:
%
% A demo for deblurring objects in images having supports with arbitrary
% shapes, by using SPA.
% Here is the case (a real image taken by a mobile phone) that in an
% image the Foreground is focused and Background is blurred.
%
% This code recreates Fig.4. in
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
im = imread('abstract_pen.png');
figure(1),imshow(im),title('observation')  %Fig.4. left image

% get the object mask
im1 =  im(:,:,3)-im(:,:,1);
im2 = histeq(im1);
im3 = im2<65;
im4 = imfill(im3,'holes');
im5 = bwmorph(im4,'dilate',[3 3]);
im6 = 1-imfill(im5,'holes');
mask_image = im6;
clear im1 im2 im3 im4 im5 im6

im = double(im);


%% process the average of RGB channel
%  here we only process a grayscale image which is the average of RGB
%  channel, later we do color correction to get the final RGB image

im00 = mean(im,3); % The average
[Ny,Nx] = size(im00);


%% define blur
D = 14.5;   % blurring kernel diameter, estimated
Le = 5;     % boundary extension
sigma = 1.25;  % noise level

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
sig2x = 30^2; % sig2x and rho work well for Wiener (both boundary processing and Wiener "optimized" at the same time). Adapted to the [0 255] range.
u = xc/NDx;
v = yc/NDy;
rho = 0.65;
PX_f = sig2x*4*log(rho)^2./((log(rho)^2 + 4*pi^2*u.^2).*(log(rho)^2 + 4*pi^2*v.^2));
PX_f = fftshift(PX_f);

P=3;

z = SPA(y,kernel,mask,PX_f,sigma,P);
figure(3),imagescTSc(z);title('result of SPA interpolation')


%% restoration
% using L2-relaxed-L0 deblurring Matlab Toolbox v2.0
% please check the toolbox for details

Niter = 10;

% Load parameters for L2-r-L0 deblurring with Niter=10:
% All three representations are active.
% From Table I, ConDy 10 (number of iterations) in Portilla TIP2015,
% parameters jointly optimizing SSIM of the 24 test experiments.
load L2rL0_Niter10_parameters.mat

LB = 120; % Lower and upper bounds for the background paper
UB = 210;

x_hat = deblur_L2relaxedL0(z, rep, Real_psf_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB);
figure(4),imagescTSc(x_hat);title('result of deconvolution')


%% final result
final_image = x_hat(Le+1:end-Le,Le+1:end-Le).* mask_smooth + im00.*(1-mask_smooth);  % no sharp edge
figure(5),imagescTSc(final_image);title('final grayscale result')


%% color correction for getting color image result
for color = 1:3
    
    im00 = im(:,:,color);
    
    % compute background color mean for color correction
    im0 = im00 - mean2(im00);
    [Ny,Nx] = size(im0);
    
    % First, it extracts a big chunk (as big as possible) of background area,
    % from which to train the predictors
    L = 3;
    im0_ext = bound_extension(im0,L,L,'mirror');
    
    % High-pass filter the image to mark non-background area(in focus area)
    
    im1 = conv2([-1 2 -1]/4,[ 0 1 0], im0_ext, 'valid');
    im2 = conv2([0 1 0],[-1 2 -1].'/4,im0_ext, 'valid');
    
    res = im1.^2 + im2.^2;
    
    res = res(1+L:end-L,1+L:end-L);
    
    mu = mean2(res);
    
    aux = zeros(Ny,Nx);
    aux(2:end-1,2:end-1) = res;
    res = aux;
    
    T = 6;
    
    Bin = double(res>T*mu);
    
    % Remove isolated outliers
    Bin = conv2([1 1 1]/3,[1 1 1]/3,Bin,'same');
    Bin = double(Bin>0.5);
    
    % Now try to find the largest area background rectangle
    % First define a MxM grid of posible shapes/sizes
    % binning the image
    
    N = 2; % 2^N scale zoom out
    Bin_small = blurDn(Bin,N,ones(2,2)/4);  %blur and downsample the image
    
    [Nsy,Nsx] = size(Bin_small);
    
    smin = 3;
    step = 2; % Watch out: it must be an even number, so keeping the rectangle kernel dimensions odd
    sxsy = 0;
    n = 0;
    for sx = Nsx-1:-step:smin,
        for sy = Nsy-1:-step:smin,
            if sy*sx>=sxsy,
                aux = conv2(ones(sy,1),ones(sx,1),Bin_small,'valid');
                valid = double(abs(aux)<1e-5);
                if sum(valid(:))>0,
                    n = n+1;
                    syv(n) = sy;
                    sxv(n) = sx;
                    [fooy foox] = find(valid);
                    xv(n) = foox(1) + (sx-1)/2;
                    yv(n) = fooy(1) + (sy-1)/2;
                    ind = 0;
                    ind = ind + double(yv(n) - (sy-1)/2 ==1);
                    ind = ind + double(xv(n) - (sx-1)/2 ==1);
                    ind = ind + double(yv(n) + (sy-1)/2 == Nsy);
                    ind = ind + double(xv(n) + (sx-1)/2 == Nsx);
                    if (sx*sy > sxsy) && (ind>=2),
                        sxsy = sx*sy;
                        %showIm(Bin); title('background area');hold on
                        xvmax = 2^N*(xv(n)-1)+1;
                        yvmax = 2^N*(yv(n)-1)+1;
                        sxmax = 2^N*sxv(n);
                        symax = 2^N*syv(n);
                        %plot([xvmax-(sxmax-2^N)/2 xvmax+(sxmax-2^N)/2 xvmax+(sxmax-2^N)/2 xvmax-(sxmax-2^N)/2 xvmax-(sxmax-2^N)/2],[yvmax-(symax-2^N)/2 yvmax-(symax-2^N)/2 yvmax+(symax-2^N)/2 yvmax+(symax-2^N)/2 yvmax-(symax-2^N)/2],'y');
                        %hold off
                        %drawnow; %pause
                    end
                    
                end
            end
        end
    end
    
    bg_chunk = im00(max(yvmax-symax/2,1):min(yvmax+symax/2,Ny),max(xvmax-sxmax/2,1):min(xvmax+sxmax/2,Nx));
    
    % Now obtain the mean
    mu_bg = mean2(bg_chunk);
    
    mu_rgb_correction = mean2(x_hat(max(yvmax-symax/2,1):min(yvmax+symax/2,Ny),max(xvmax-sxmax/2,1):min(xvmax+sxmax/2,Nx)));
    
    final_image_RGB(:,:,color) = mu_bg/mu_rgb_correction.*final_image;
    mask_RGB(:,:,color) = mask_smooth;
    
end

final_image_RGB = final_image_RGB .* mask_RGB + double(im).*(1-mask_RGB);
figure(6),imshow(uint8(final_image_RGB));title('final RGB result') %Fig.4 right image

