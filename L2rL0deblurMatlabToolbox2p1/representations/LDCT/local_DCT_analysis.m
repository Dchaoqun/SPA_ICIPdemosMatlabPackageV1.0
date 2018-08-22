function [local_dct,pind,ovl] = local_DCT_analysis(im,N,rf,ovl);

% Implementation of Local DCT analysis
% [local_dct,pind] = local_DCT_analysis(im,N,rf,ovl);
% INPUT:
%       im : image to be processed
%       N  : block side size (square)    (default = 8)
%       rf : redundancy factor along each dimension (default = 2)
%       ovl: overlapping mask (OPTIONAL). It contains the square root of
%            the number of blocks contributing to each pixel in the image. 
% OUTPUT:
%       local_dct: 2-D DCT vector of the possibly overlapping blocks,
%                   sequentialized row-wise.
%       pind: vector containing, in this order:
%           - Ny, Nx: image size
%           - N: block side size
%           - rf: redundancy factor
%
% NOTE: blocks are properly weighted in the image domain so the
% representation is a Parseval frame (it preserves the Euclidean norm)
% NOTE2: the DC component is removed from the image and included as the
% last sample in the representation.
%
% See also: local_DCT_synthesis, my_DCT, my_IDCT
%
% Javier Portilla, CSIC, Madrid
% April 2009

%im = single(im);

if~exist('N'),
    N = 8;
end

if~exist('rf'),
    rf = 2;
end

[Ny,Nx] = size(im);

DC = sum(im(:))/sqrt(Ny*Nx);
im = im - DC/sqrt(Ny*Nx); % Subtract the mean

Nby = rf*Ny/N;
Nbx = rf*Nx/N;

% Use compute_mask_LDCT.m, instead!
%
% if~exist('ovl'),
% 
% % First it computes the overlapping mask, telling us the redundancy of each
% % pixel in the image
% ovl = zeros(size(im));
% for nby = 1:Nby-(rf-1),
%     y = (nby-1)*N/rf + 1;
%     for nbx = 1:Nbx-(rf-1),
%         x = (nbx-1)*N/rf + 1;
%         ovl(y:y+N-1,x:x+N-1) = ovl(y:y+N-1,x:x+N-1) + ones(N,N);
%     end
% end
% ovl = sqrt(ovl);
%
% end % if ovl does not exist

%ovl = single(ovl);

% Now it computes the transform itself
local_dct = zeros((Nby-(rf-1))*(Nbx-(rf-1))*N^2+1,1);

% % Compute weights to multiply DFT coefficients
ww = (exp(-i*(0:N-1)*pi/(2*N))/sqrt(2*N)).';
ww(1) = ww(1) / sqrt(2);
ww = 2*ww;  % Double the weights for even-length case  
W = ww(:,ones(1,N));

n = 1;
for nby = 1:Nby-(rf-1),
    y = (nby-1)*N/rf + 1;
    for nbx = 1:Nbx-(rf-1),
        x = (nbx-1)*N/rf + 1;
        block = im(y:y+N-1,x:x+N-1);
        block = block./ovl(y:y+N-1,x:x+N-1);
        %dct_b = dct2(block); % SLOW
        %dct_b = dct(dct(block).').'; % a bit faster, equivalent
        dct_b = my_DCT(my_DCT(block, W).',W); % much faster! equivalent
        local_dct(n:n+N^2-1) = dct_b(:);
        n = n + N^2;
    end
end    
local_dct(end) = DC;
pind = [Ny; Nx; N; rf];