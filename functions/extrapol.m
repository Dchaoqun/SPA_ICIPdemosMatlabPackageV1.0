function imx = extrapol(im, mask)
% EXTRAPOL Return an inter/extrapolation of an incomplete image.
%   It is based on using a very peaky local integration kernel (r.^(-7)), such that,
%   for close-by known values the inter/extrapolated pixels are an average of those
%   neighbor pixels, whereas for far-away pixels the inter/extrapolated pixels are
%   a much more extended area average.
%   The trick consists of dividing always by the area of the kernel on the
%   known-pixels' area. This ensures that it is a true average.
%
%   SYNTAX:
%     imx = extrapol(im,mask);
%
%   INPUT: 
%     im     - image (incomplete)
%     mask   - binary mask indicating the pixels to be kept (1's) and the ones
%              to be filled in (0's)
%
%   OUTPUT:
%     imx    - an image keeping the original pixels (with mask ==1) and with new
%              values in the rest of the pixels such that the result is smooth

% Used in the code for ICIP publication:
% - C. Dong, J. Portilla, "Spectral pre-adaptation for two-step
%   arbitrary-shape-support image restoration", 2017 IEEE International
%   Conference on Image Processing (ICIP), pp. 3515-3519, Sept 2017.
%
% Version 1.0, August 2018
%
% Copyright (c) 2018
% Javier Portilla <javier.portilla@csic.es>
% Chaoqun Dong <chaodong999@gmail.com>
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

Np = size(im,1);
N = size(mask,1);
img = zeros(N);
img((N-Np)/2+1:(N+Np)/2,(N-Np)/2+1:(N+Np)/2) = im;
img = img.*mask;
int = double(mask);

% Define the kernel
fil = 1./(meshgrid(-N/2:N/2-1,-N/2:N/2-1).^2+(meshgrid(-N/2:N/2-1,-N/2:N/2-1).').^2).^3.5;
fil = fftshift(fil);
fil(1,1) = 1000;	
F = real(fft2(fil));

% Filter the image, doing a local average
imx = ifft2(F.*fft2(img))./ifft2(F.*fft2(int));  %divide the filtered image by the filtered mask
imx = imx.*(1-mask)+img.*mask;

