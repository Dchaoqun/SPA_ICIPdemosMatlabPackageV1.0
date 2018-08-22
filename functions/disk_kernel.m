function [psf_f_o,psf_f,kernel,Ly,Lx,keveny,kevenx,psf_sub_f] = disk_kernel(D,Ny,Nx,Le)
% DISK_KERNEL Return a disk blurring kernel (a rough approximation of a defocus PSF)
% with/without boundary handling.
%
%   SYNTAX:
%     [psf_f_o,psf_f,kernel,Ly,Lx,keveny,kevenx,psf_sub_f] = disk_kernel(D,Ny,Nx,Le)
%
%   INPUT:
%    D            - diameter of the disk kernel
%    Nx, Ny       - size of the test image
%    Le           - number of pixels for boundary extension
%
%   OUTPUT:
%    psf_f_o           - PSF image with original size, in Fourier domain
%    psf_f             - PSF image with boundary extension, in Fourier domain
%    kernel            - blurring kernel in spatial domian
%    Ly,Lx             - half size of the kernel, for boundary extension
%    keveny,kevenx     - odd/even flag for the size of the kernel
%    psf_sub_f         - subsampled PSF, in Fourier domain

% Used in the code for ICIP publications:
% - C. Dong, J. Portilla, "Maximum likelihood interpolation for
%   aliasing-aware image restoration", 2016 IEEE International Conference
%   on Image Processing (ICIP), pp. 564-568, Sept 2016.
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
% Version 1.0, August 2018
%
% Copyright (c) 2018
% Chaoqun Dong <chaodong999@gmail.com>
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

R = D/2;
cy = Ny/2+1; % Ny, Nx assumed to be even
cx = Nx/2+1;
[u,v] = meshgrid(-Nx/2:Nx/2-1,-Ny/2:Ny/2-1);
u = u/Nx;
v = v/Ny;

f = sqrt(u.^2 + v.^2);

Ky = ceil(D/2)*2 + 1;  % size of the kernel
Kx = ceil(D/2)*2 + 1;  % assume kernel is square because the psf is a disk

Ly = ceil((Ky-1)/2);  % half size of the kernel
Lx = ceil((Kx-1)/2);

keveny = double(rem(Ky,2)==0);  % odd/even size handling
kevenx = double(rem(Kx,2)==0);


%% psf

z = 2*pi*f*R;
Disk_f = 2*besselj(1,z)./(z +double(z==0));
Disk_f(cy,cx) = 1;              % set the center of the disk to be 1

disk = fftshift(real(ifft2(fftshift(Disk_f))));

% a mask in the spatial domain to smooth the uniform disk
mask = [1 2 1]'*[1 2 1]/16;
disk_smooth = conv2(disk,mask,'same');

kernel = disk_smooth(Ny/2+1-Ly+keveny:Ny/2+1+Ly,Nx/2+1-Lx+kevenx:Nx/2+1+Lx);
kernel = kernel/sum(sum(kernel));
psf = zeros(Ny,Nx);
psf(Ny/2+1-Ly+keveny:Ny/2+1+Ly,Nx/2+1-Lx+kevenx:Nx/2+1+Lx) = kernel;

% FT
psf_f_o = fft2(fftshift(psf));

%boundary handling
NDy = Ny + 2*Le;
NDx = Nx + 2*Le;
D_f = zeros(NDy,NDx);
D_f(NDy/2+1-Ly+keveny:NDy/2+1+Ly,NDx/2+1-Lx+kevenx:NDx/2+1+Lx) = kernel; % Warning: Nx, Ny are assummed to be even numbers!
D_f = fft2(fftshift(D_f));
psf_f = D_f;

%% sub psf

kernel_sub = kernel(1+mod(Ly,2):2:end,1+mod(Lx,2):2:end);
kernel_sub = kernel_sub/sum(sum(kernel_sub));

Nsy = size(Ly-keveny+1:2:Ny-Ly,2);
Nsx = size(Lx-kevenx+1:2:Nx-Lx,2);

psf_sub = zeros(Nsy,Nsx);
psf_sub(floor(Nsy/2)+1-floor(Ly/2)+keveny:floor(Nsy/2)+1+floor(Ly/2),floor(Nsx/2)+1-floor(Lx/2)+kevenx:floor(Nsx/2)+1+floor(Lx/2)) = kernel_sub;

% FT
psf_sub_f = fft2(fftshift(psf_sub));

