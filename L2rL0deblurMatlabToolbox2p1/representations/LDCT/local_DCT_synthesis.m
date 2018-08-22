function im = local_DCT_synthesis(local_dct,pind,ovl);

% Implementation of Local DCT synthesis
%
% im = local_DCT_synthesis(local_dct,pind,ovl);
% INPUT:
%       local_dct: 2-D DCT vector of the possibly overlapping blocks,
%                   sequentialized row-wise.
%       pind: vector containing, in this order:
%           - Ny, Nx: image size
%           - N: block side size
%           - rf: redundancy factor
%       ovl: overlapping mask (OPTIONAL). It contains the square root of
%            the number of blocks contributing to each pixel in the image. 
% OUTPUT:
%       im : reconstructed image
%
% NOTE: Exact Parseval frame (perfect reconstruction, preserves the Euclidean
%       norm).
%
% See also: local_DCT_analysis, my_DCT, my_IDCT
%
% Javier Portilla, CSIC, Madrid
% April 2009


Ny = pind(1); Nx = pind(2); % Image dimensions
N = pind(3);    % Block dimensions (N x N)
rf = pind(4);   % Redundancy factor

Nby = rf*Ny/N;
Nbx = rf*Nx/N;

% Use compute_mask_LDCT.m, instead!
% if~exist('ovl'),
% 
% % It computes the overlapping mask, which tells us the redundancy of each
% % pixel in the image
% ovl = zeros(Ny,Nx);
% for nby = 1:Nby-(rf-1),
%     y = (nby-1)*N/rf + 1;
%     for nbx = 1:Nbx-(rf-1),
%         x = (nbx-1)*N/rf + 1;
%         ovl(y:y+N-1,x:x+N-1) = ovl(y:y+N-1,x:x+N-1) + ones(N,N);
%     end
% end
% ovl = sqrt(ovl);
% 
% end % ovl does not exist

% Compute weights
ww = sqrt(2*N) * exp(j*(0:N-1)*pi/(2*N)).';
% Compute precorrection factor
ww(1) = ww(1)/sqrt(2);
W = ww(:,ones(1,N));

% Now it reverts the transform itself
im = zeros(Ny,Nx);
n = 1; N2 = N^2;
for nby = 1:Nby-(rf-1),
    y = (nby-1)*N/rf + 1;
    for nbx = 1:Nbx-(rf-1),
        x = (nbx-1)*N/rf + 1;
        dct_b = local_dct(n:n+N2-1);
        dct_b = reshape(dct_b,N,N);
        %block = idct((idct(dct_b)).').';
        block = my_IDCT((my_IDCT(dct_b,W)).',W); % much faster (40% less time), equivalent
        im(y:y+N-1,x:x+N-1)= im(y:y+N-1,x:x+N-1) + block./ovl(y:y+N-1,x:x+N-1);
        n = n + N2;
    end
end    

DC = local_dct(end);
im = im + DC/sqrt(Ny*Nx); % Add back the mean
