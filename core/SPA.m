 function [z, z_init] = SPA(y, kernel, mask, PX_f, sigma, P)
% SPA Complete (extend and/or interpolate) an incomplete blurry observation.  
%   Z = SPA(Y, KERNEL, MASK, PX_F, SIGMA, P) returns a maximum-likelihood
%   interpolation Z according to a Gaussian spectral model.
%
%   "SPA" stands for Spectral Pre-Adaptation. 
%   It does a Maximum-Likelihood interpolation as a likely guess of the unknown 
%   pixels of an incomplete observation of a blurry image ("y"), based on 
%   consistency with the blurring kernel ("kernel"), the noise standard deviation
%   ("sigma", assumed Gaussian), assumed both known, and a simple AR-1 image model.
%   The known and unknown values of "y" correspond to "TRUE"s (or 1s) and "FALSE"s
%   (or 0s) in the "mask" array, respectively.
%   P is the number of iterations (e.g., 4 for typical situations, although it may be
%   necessary to increase this up to 50 or more for big kernels and little noise).
%
%   [Z, Z_INIT] = SPA(Y, KERNEL, MASK, PX_F, SIGMA, P) returns the maximum-likelihood
%   interpolation Z and a first guess Z_INIT. Z_INIT is an optional output 
%   giving the first guess obtained with the "extrapol.m" method. It is 
%   used as the starting point for maximum-likelihood interpolation to reduce
%   the number of iterations.   
%
%   INPUT:
%     y           - an incomplete observation of a blurry image
%                   (preferrably in the [0 255] range - otherwise please,
%                   normalize it to that range)
%     kernel      - blurring kernel
%     mask        - a mask indicating the unknown pixels 
%     PX_f        - power spectral density model of the uncorrupted image
%     sigma       - standard deviation of noise
%     P           - number of iterations
%
%   OUTPUT:
%     z           - the interpolated image
%     z_init      - a first guess for 'z' obtained with 'extrapol.m'
%
%   It is useful for:
%   (A) Boundary handling (by extending the boundaries and imposing circular
%   boundaries condition) (see ICIP 2014)
%   (B) Deblurring aliased observations, in which case the known
%   subset is a uniformly sampled version of the complete set of pixels
%   (see ICIP 2016).
%   (C) Deblurring blurred objects having supports with arbitrary shapes
%   (e.g., out-of-focus background in presence of in-focus foreground)
%   (see ICIP 2017).
%
% RELATED PUBLICATIONS:
% - J. Portilla, "Maximum likelihood extension for non-circulant
%   deconvolution", Image Processing (ICIP) 2014 IEEE International
%   Conference, pp. 4276-4279, Oct 2014.
% - C. Dong, J. Portilla, "Maximum likelihood interpolation for
%   aliasing-aware image restoration", 2016 IEEE International Conference
%   on Image Processing (ICIP), pp. 564-568, Sept 2016.
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.

% Version 1.0, August 2018
%
% Copyright (c) 2018
% Javier Portilla <javier.portilla@csic.es>
% Chaoqun Dong <cdongae@connect.ust.hk>
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
%fprintf(1,'Info:  Processing by SPA...\n');
    
mask = logical(mask); % in case it was given as numbers, in binary, etc.

% blurring kernel h
[Ky,Kx] = size(kernel);

Ly = ceil((Ky-1)/2);  % half size of the kernel
Lx = ceil((Kx-1)/2);
keveny = double(rem(Ky,2)==0);
kevenx = double(rem(Kx,2)==0);

%observation
y00_obsv = y;
mu = mean2(y00_obsv.*mask)/mean2(mask);
y00 = (y00_obsv - mu).*mask;

% the blurring kernel, in the Fourier domain.
[NDy, NDx] = size(y);
h = zeros(NDy,NDx);
h(NDy/2+1-Ly+keveny:NDy/2+1+Ly,NDx/2+1-Lx+kevenx:NDx/2+1+Lx) = kernel; % Warning: it assumes Nx, Ny are even numbers!
H_f = fft2(fftshift(h));  %% FT of the kernel

% PSD model for the (incomplete) observation
PZ_f = abs(H_f).^2.*PX_f + sigma^2;

Q = 25;  % Length of the linear prediction. This value works well for all tried examples

% Initiallize surrounded by zeros

% Watch out! Here it is assuming "bytes"
z = mu*ones(NDy,NDx);
z(mask) = y00_obsv(mask);
z = extrapol(z, mask);
z_init = z;

% Pre-processing

% First compute q'
aux = zeros(NDy,NDx);
% part of b(y): Fo*y
aux(mask) = y00(mask);
Y_f = fft2(aux);
% part of b(y): F.D_pz^-1(Fo*.y), here we have Fi* = F*.E^T, so Fi = E.F
% So q = -b(y)
q = real(ifft2(Y_f./PZ_f));
q(mask) = zeros;

z0 = z - mu; % We NEED to remove the mean
ze0 = -z0;
ze0(mask) = zeros;
ze = ze0;

lambda = 2*sigma^2;   % fastest convergence

S = NDy*NDx - sum(sum(double(mask))); % Number of degrees of freedom to be optimized

U = zeros(S,Q);
Z = zeros(NDy*NDx,Q);
mask2 = ones(NDy,NDx);
mask2(mask) = zeros;
mask2 = mask2(:);

for m = 1:P,
       
    % Applies the original method (see ICIP 2014) but saving the vector
    % values, so then can do a linear prediction of the convergence limit
    for n = 1:Q,
        % part of Az,cause A= Fi.D_pz^-1.Fi*, Ze_f = Fi*z, because it is -Az, so
        % ze0 = -z0
        Ze_f = fft2(ze);
        % part of Az, F.D_pz^-1.Ze_f
        aux = real(ifft2(Ze_f./PZ_f));
        % part of Az, E(...)
        aux(mask) = zeros;
        ze_old = ze;
        ze = ze + lambda*(q-aux);
        u = ze - ze_old;
        u = u(:);
        U(:,n) = u(mask2==1);
        Z(:,n) = ze(:);
    end
    
    % Estimate u(N) from the previous Q-1 ones
    % (see ICIP 2016)
    
    Ur = U(:,1:Q-1);
    c = - pinv(Ur)*U(:,Q);
    c = [c; 1];
    c = c/sum(c);
    
    Zr = Z(:,1:Q);
    ze = Zr*c;
    ze = reshape(ze,NDy,NDx);
  
    fprintf(1,'Info:  Processing by SPA, current p = %d out of %d ...\n',m,P);
    
end 

z = -ze;
z(mask) = y00(mask);   % get back known

z = z + mu;



