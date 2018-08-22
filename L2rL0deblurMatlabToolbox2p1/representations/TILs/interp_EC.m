function im = interp_EC(im0, typ);

% function im = interp_EC(im0, typ);
%
% It enlarges an image in a factor of 2 along each dimension.
% It uses an efficient convolution (EC) technique of mine.
% It provides a performance close to a standard spline, but
% it is significantly faster.
%
%       typ:
%           1: Linear
%           2: Cubic
%           3: The next one ;-) (default)
%           4: The other :-)
%
% Javier Portilla
% Madrid, 22 October 2008.

if ~exist('typ'),
    typ = 3;
end    

% First, it extends the original image so the interpolated pixels have the
% expected dimensions
im0 = [im0 im0(:,end); im0(end,:) im0(end,end)];

im = zeros(2*size(im0,1)-1,2*size(im0,2)-1);
% x -
% - -
im(1:2:end,1:2:end) = im0;

aux = im0;
if typ==1,
    fil = [1 1]/2;   % Linear filter
elseif typ==2,
    fil = [-1 9 9 -1]/16;   % Cubic filter
elseif typ==3,
    fil = [3 -25 150 150 -25 3]/256;   % Next one
elseif typ==4,
    fil = [-5 49 -245 1225 1225 -245 49 -5]/2048;
else    
    %fil = design_interp_filter(typ)';
    error('Interpolation filter parameter ranges from 1 to 4 in this implementation');
    return;
end

L = typ; % = length(fil)/2;   It must be even!

if L>1,
% Boundary extension for filtering in the spatial domain
aux = [aux(L-1:-1:1,:); aux; aux(end:-1:end-L+2,:)];
aux = [aux(:,L-1:-1:1) aux aux(:,end:-1:end-L+2)];
end

% x O
% - -
im_int = conv2(aux,fil,'valid');
im(1:2:end,2:2:end) = im_int(L:end-L+1,:); 

% x -
% O -
im_int = conv2(aux,fil','valid');
im(2:2:end,1:2:end) = im_int(:,L:end-L+1); 

% x -
% - O
im_int = conv2(im_int,fil,'valid');
im(2:2:end,2:2:end) = im_int; 


% Remove the file/row added at the beginning
im = im(1:end-1,1:end-1); 
