function x_r = reconTDP1(pyr,pind)

% Reconstruction of an image from the Truadic-dyadic pyramid (TDP) or Translation Invariant Laplacian Separable (TILs)
% representation
%
% x_r = reconTDP1(pyr,pind)
%
% See buildTDP1.m
%
% Javier Portilla, Madrid, 20 de Octubre de 2011


% The orthogonal basis (conjugated, for synthesis):

kh = [-1 2 -1]'/(sqrt(3*6));  % Bar detector, high frequency
km = [1 0 -1]'/(sqrt(3*2));    % Edge detector, medium frequency
kl = [1 1 1]'/(sqrt(3*3));     % Low pass


Ny = pind(1,1); Nx = pind(1,2);
Nb = 8;
Nband = size(pind,1);
LPR = (Nband>=8)&&(Nband/2 ~=floor(Nband/2));
Nsc = (Nband-LPR)/Nb;
  
Ny = Ny/2^(Nsc-1);
Nx = Nx/2^(Nsc-1);
if LPR,
    x_r = pyrBand(pyr,pind,Nband);
else
    x_r = zeros(Ny,Nx);  % No LPR!
end


for nsc = Nsc:-1:1,    
   
    if Nsc>1.
    x_r = interp_EC(x_r,4)/2;
    end
    xf_ll = x_r;
    
    xf_lm = pyrBand(pyr,pind,nsc*Nb);
    xf_ml = pyrBand(pyr,pind,nsc*Nb - 1);
    xf_mm = pyrBand(pyr,pind,nsc*Nb - 2);   
    xf_lh = pyrBand(pyr,pind,nsc*Nb - 3);
    xf_hl = pyrBand(pyr,pind,nsc*Nb - 4);
    xf_mh = pyrBand(pyr,pind,nsc*Nb - 5);   
    xf_hm = pyrBand(pyr,pind,nsc*Nb - 6);
    xf_hh = pyrBand(pyr,pind,nsc*Nb - 7);
   
    ba = xf_lm;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = conv2(kl,km,ba,'valid');
    
    ba = xf_ml;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(km,kl,ba,'valid');

    ba = xf_mm;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(km,km,ba,'valid');

    ba = xf_lh;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(kl,kh,ba,'valid');
    
    ba = xf_hl;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(kh,kl,ba,'valid');

    ba = xf_mh;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(km,kh,ba,'valid');

    ba = xf_hm;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(kh,km,ba,'valid');
    
    ba = xf_hh;
    ba = [ba(end,:); ba; ba(1,:)];
    ba = [ba(:,end) ba ba(:,1)];
    x_r = x_r + conv2(kh,kh,ba,'valid');

    if Nsc>1 || LPR,
        ba = xf_ll;
        ba = [ba(end,:); ba; ba(1,:)];
        ba = [ba(:,end) ba ba(:,1)];
        x_r = x_r + conv2(kl,kl,ba,'valid');
    end
     
    Ny = 2*Ny; Nx = 2*Nx;
    
end    


