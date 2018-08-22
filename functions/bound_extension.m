function im_ext = bound_extension(im,By,Bx,type)
% BOUND_EXTENSION Extend an image for avoiding boundary artifacts,
%
%   SYNTAX:
%     im_ext = bound_extension(im,By,Bx,type)
%
%   INPUT:
%     By, Bx     - widths of the added stripes.
%     type:   
%           'mirror'            - Mirror extension
%           'mirror_nr'         - Mirror without repeating the last pixel
%           'circular'          - fft2-like
%           'double_mirror'     - Mirror spatially and in values
%           'constant'          - repeat the last row/column
%           'zeros'             - all zeros
%
%   OUTPUT:
%      im_ext    - extended image

% Used in the code for ICIP publications:
% - J. Portilla, "Maximum likelihood extension for non-circulant
%   deconvolution", Image Processing (ICIP) 2014 IEEE International
%   Conference, pp. 4276-4279, Oct 2014.
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

[Ny,Nx,Nc] = size(im);

im_ext = zeros(Ny+2*By,Nx+2*Bx,Nc);
im_ext(By+1:Ny+By,Bx+1:Nx+Bx,:) = im;

if strcmp(type,'mirror'),

    im_ext(1:By,:,:) = im_ext(2*By:-1:By+1,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,2*Bx:-1:Bx+1,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(Ny+By:-1:Ny+1,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Nx+Bx:-1:Nx+1,:);
    im_ext(1:By,1:Bx,:) = im_ext(2*By:-1:By+1,2*Bx:-1:Bx+1,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+By:-1:Ny+1,Nx+Bx:-1:Nx+1,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(2*By:-1:By+1,Nx+Bx:-1:Nx+1,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(Ny+By:-1:Ny+1,2*Bx:-1:Bx+1,:);

elseif strcmp(type,'mirror_nr'),    
        
    im_ext(1:By,:,:) = im_ext(2*By+1:-1:By+2,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,2*Bx+1:-1:Bx+2,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(Ny+By-1:-1:Ny,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Nx+Bx-1:-1:Nx,:);
    im_ext(1:By,1:Bx,:) = im_ext(2*By+1:-1:By+2,2*Bx+1:-1:Bx+2,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+By-1:-1:Ny,Nx+Bx-1:-1:Nx,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(2*By+1:-1:By+2,Nx+Bx-1:-1:Nx,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(Ny+By-1:-1:Ny,2*Bx+1:-1:Bx+2,:);
        
elseif strcmp(type,'circular'),        
        
    im_ext(1:By,:,:) =  im_ext(Ny+1:Ny+By,:,:);
    im_ext(:,1:Bx,:) = im_ext(:,Nx+1:Nx+Bx,:);
    im_ext(Ny+1+By:Ny+2*By,:,:) = im_ext(By+1:2*By,:,:);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = im_ext(:,Bx+1:2*Bx,:);
    im_ext(1:By,1:Bx,:) = im_ext(Ny+1:Ny+By,Nx+1:Nx+Bx,:);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(By+1:2*By,Bx+1:2*Bx,:);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = im_ext(Ny+1:Ny+By,Bx+1:2*Bx,:);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = im_ext(By+1:2*By,Nx+1:Nx+Bx,:);

elseif strcmp(type,'constant'),        

    im_ext(1:By,:,:) = repmat(im_ext(By+1,:,:),[By 1]);
    im_ext(:,1:Bx,:) = repmat(im_ext(:,Bx+1,:),[1 Bx]);
    im_ext(Ny+1+By:Ny+2*By,:,:) = repmat(im_ext(By+Ny,:,:),[By 1]);
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = repmat(im_ext(:,Bx+Nx,:),[1 Bx]);
    im_ext(1:By,1:Bx,:) = repmat(im_ext(By+1,Bx+1,:),[By Bx]);
    im_ext(Ny+1+By:Ny+2*By,Nx+1+Bx:Nx+2*Bx,:) = repmat(im_ext(By+Ny,Bx+Nx,:),[By Bx]);
    im_ext(1:By,Nx+1+Bx:Nx+2*Bx,:) = repmat(im_ext(By+1,Bx+Nx,:),[By Bx]);
    im_ext(Ny+1+By:Ny+2*By,1:Bx,:) = repmat(im_ext(By+Ny,Bx+1,:),[By Bx]);
    
elseif strcmp(type,'zeros'),        
        
    im_ext(1:By,:,:) =  zeros;
    im_ext(:,1:Bx,:) = zeros;
    im_ext(Ny+1+By:Ny+2*By,:,:) = zeros;
    im_ext(:,Nx+1+Bx:Nx+2*Bx,:) = zeros;
   
end    
        