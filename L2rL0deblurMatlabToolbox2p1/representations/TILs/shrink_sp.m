function y = shrink_sp(x, N);

% Efficient spatial shrinkage of an image by low-pass filtering it in the spatial
% domain and subsampling with a factor 2 along each dimension.
%
%   y = shrink_sp(x, N);
%
% It applies a low-pass filter with N vanishing moments and cut-off
% frequency of f = 0.25 cyc/pixel, symmetric in frequency domain,
% which has the structure [...0 d 0 c 0 b 1 b 0 c 0 d 0...]/2.
% Then it subsamples the image in a factor of 2x2.
% Given this structure the filtering + subsampling can be done very
% fast, typically twice as fast as using the direct method (separable
% convolution with the given filter and then subsampling).
% Boundary extension is done by mirror reflection.
%
% It requires the image dimensions to be even.
%
% Javier Portilla, Madrid, October 2008

[Ny,Nx] = size(x);

if N==1, % linear
    fil = [1 2 1]/4;
elseif N==2, % cubic    
    fil = [-1 0 9 16 9 0 -1]/32;
elseif N==3,    
    fil = [3 0 -25 0 150 256 150 0 -25 0 3]/512;
elseif N==4,
    fil = [-5 0 49 0 -245 0 1225 2048 1225 0 -245 0 49 0 -5]/4096;
else    
    fil = design_interp_filter(N)';
    aux = zeros(4*N-1,1);
    aux(1:2:end) = fil;
    aux(2*N) = 0.5;
    fil = aux;
end    
 
%y = blurDn(x,1,fil);   % + 0.45 dB, I do not know why (?)

x1 = x(1:2:Ny,:);
x2 = x(2:2:Ny,:);
% Boundary extension for filtering in the spatial domain
x2 = [x2(N:-1:1,:);x2;x2(end:-1:end-N+2,:)];
aux = conv2(x2,fil(1:2:end)','valid');
yx2 = x1/2 + aux;

yy1 = yx2(:,1:2:Nx);
yy2 = yx2(:,2:2:Nx);
% Boundary extension for filtering in the spatial domain
yy2 = [yy2(:,N:-1:1) yy2 yy2(:,end:-1:end-N+2)];
aux = conv2(yy2,fil(1:2:end),'valid');
y = yy1/2 + aux;
