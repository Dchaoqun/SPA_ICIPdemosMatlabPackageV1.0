function [pyr,pind] = buildTDP1(x,Nsc, LPR);

% Triadic-dyadic pyramid, or Translation Invariant Laplacian Separable
% (TILs) representation.
%
%   [pyr,pind] = buildTDP1(x,Nsc, LPR);
%           x:         Image to be decomposed (double, grayscale)
%           Nsc:     Number of scales
%           LPR:    Including (1) or not (0) Low Pass Residual
%           pyr:      pyramid
%           pind:    subbands' sizes
%
% Javier Portilla, Madrid, 20 de Octubre de 2011

[Ny,Nx] = size(x);

% The orthogonal basis:

kh = [-1 2 -1]'/(sqrt(3*6));  % Bar detector, high frequency
km = [-1 0 1]'/(sqrt(3*2));    % Edge detector, medium frequency
kl = [1 1 1]'/(sqrt(3*3));     % Low pass

pyr = [];
pind = [];  

xf_ll = x;

for nsc = 1:Nsc,

    xf_ll = [xf_ll(end,:); xf_ll; xf_ll(1,:)];
    xf_ll = [xf_ll(:,end) xf_ll xf_ll(:,1)];

    
% 9 bands, one of them is the low-pass
xf_h = conv2(xf_ll,kh,'valid');
xf_hh = conv2(xf_h,kh','valid');
xf_hm = conv2(xf_h,km','valid');
xf_hl = conv2(xf_h,kl','valid');
xf_m = conv2(xf_ll,km,'valid');
xf_mh = conv2(xf_m,kh','valid');
xf_mm = conv2(xf_m,km','valid');
xf_ml = conv2(xf_m,kl','valid');
xf_l = conv2(xf_ll,kl,'valid');
xf_lh = conv2(xf_l,kh','valid');
xf_lm = conv2(xf_l,km','valid');
xf_ll = conv2(xf_l,kl','valid');

pyr = [pyr; vector(xf_hh)];
pyr = [pyr; vector(xf_hm)];
pyr = [pyr; vector(xf_mh)];
pyr = [pyr; vector(xf_hl)];
pyr = [pyr; vector(xf_lh)];
pyr = [pyr; vector(xf_mm)];
pyr = [pyr; vector(xf_ml)];
pyr = [pyr; vector(xf_lm)];

pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];
pind = [pind; Ny Nx];

if Nsc>1,
xf_ll = 2*shrink_sp(xf_ll,4);

Ny = Ny/2;
Nx = Nx/2;
end

end    

if LPR,
    pyr = [pyr; vector(xf_ll)];
    pind = [pind; Ny Nx];
end
