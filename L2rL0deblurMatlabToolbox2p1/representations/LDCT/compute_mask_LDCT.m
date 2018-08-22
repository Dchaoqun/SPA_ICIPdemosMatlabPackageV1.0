function ovl = compute_mask_LDCT(im, N, rf);

% It computes the overlapping mask, telling us the redundancy of each
% pixel in the image
%
% ovl = compute_mask_LDCT(im, N, rf);
%
% im: image of the same dimensions
% N: block size
% rf = redundancy factor along each dimension
%
% Used by local_DCT_analysis, local_DCT_synthesis
%
% Javier Portilla
% Madrid 2009

[Ny,Nx] = size(im);

Nby = rf*Ny/N;
Nbx = rf*Nx/N;

ovl = zeros(size(im));
for nby = 1:Nby-(rf-1),
    y = (nby-1)*N/rf + 1;
    for nbx = 1:Nbx-(rf-1),
        x = (nbx-1)*N/rf + 1;
        ovl(y:y+N-1,x:x+N-1) = ovl(y:y+N-1,x:x+N-1) + ones(N,N);
    end
end
ovl = sqrt(ovl);