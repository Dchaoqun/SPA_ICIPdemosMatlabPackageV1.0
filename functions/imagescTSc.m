function imagescTSc(im, scale, sc, cmap)
% IMAGESCTSC Display an image in True Size or scaled True Size.
%	Displaying is performed by imagesc (if sc==1), and colormap (cmap) is
%	set by default to 'gray'. It also works for full color images.
%
%   SYNTAX:
%     imagescTSc(im, scale, sc, cmap);
%
%   INPUT:
%     im      - image to display
%     scale   - image resize scale (an integer >=1)
%     sc      - when 1 the display values range from the minimum to the maximum
%     cmap    - colormap, default value is 'gray'
%
%   Example:
%       % Display a grayscale image twice the size
%       im = imread('pout.tif');
%       imagescTSc(im,2);   

% Used in the code for ICIP publications:
% - J. Portilla, "Maximum likelihood extension for non-circulant
%   deconvolution", Image Processing (ICIP) 2014 IEEE International
%   Conference, pp. 4276-4279, Oct 2014.
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

if nargin<2
    scale=1;
end
if nargin<3
    sc=1;
end
if nargin<4
    cmap = 'gray'; %default value
end

s = size(im);
s(1) = s(1)*scale;
s(2) = s(2)*scale;

if length(s)==2, % gray-scale image
    if scale>1,
        im_o = zeros(s(1),s(2));
        for nsx = 1:scale,
            for nsy = 1:scale,
                im_o(nsy:scale:end,nsx:scale:end) = im;
            end
        end
        im = im_o;
    end
else % true color image
    if scale>1,
        im_o = zeros(s(1),s(2),3);
        for nsx = 1:scale,
            for nsy = 1:scale,
                im_o(nsy:scale:end,nsx:scale:end,:) = im;
            end
        end
        im = im_o;
    end
end

% Display image
if length(s)==3,
    if sc
        im = im - min(min(min(im)));
        im = im/max(max(max(im)));
        imshow(im);
    else
        imshow(im);
    end
else
    if sc
        imagesc(im);
    else
        image(im);
    end
end

% set axes (and colormap, if applies)

axis('image'); axis('off');
if length(s)==2,
    colormap(cmap);
end
truesize;

