function x_hat = deblur_L2relaxedL0(y, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB, sig2r_v_floor)

%  x_hat = deblur_L2relaxedL0(y, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta, LB, UB, sig2r_v_floor);
%
% Given a blurred and noisy image, its degradation parameters (the blurring filter, and the standard
% deviation of the noise) and the deblurring parameters of the L2-relaxed
% L0 method (J. Portilla, A. Tristán-Vega and I.W.Selesnick, Efficient and Robust Image Restoration
% Using Multiple-Feature L2-Relaxed Sparse Analysis Priors, IEEE TIP, vol
% 24(12), pp 5046-5059, 2015), performs the image estimation (deblurring).
%
% OUTPUT:
%       x_hat:                  Deblurred image estimation.
%
% INPUT
%       y:                      A blurred and noisy image.
%       rep:                    Binary representation vector, indicating which representations are active;
%                               rep = [TILs DTCWT LDCT]; Usually all active
%                               (TILs = 1; DTCWT = 1; LDCT = 1;)
%       H_f:                    Fourier transform of the blurring filter.
%       sigma:                  Standard deviation of the noise.
%       Niter:                  Number of iterations of the estimation loop. E.g. 10, for ConDy10
%       alpha_v:                Sparsity prior parameters, [alpha1, ..,alphaN], row vector.
%       sig2r0_v:               Initial prior low-contrast texture variance, [sig2r01, ..,sig2r0N], row vector.
%       beta:                   Exponential decay factor for dynamic shrinkage.
%       LB,UB:                  Lower and upper bound for the non-degraded image estimation.
%                               Default values are LB = -Inf, UB = Inf.
%       sig2r_v_floor:          Final prior low-contrast texture variance floor. Default value is [0; 0; 0].

% Javier Portilla, Instituto de Óptica, CSIC, Madrid, 2015.
% Last update: June 2018.

if ~exist('sig2r_v_floor'),
    sig2r_v_floor = zeros(3,1);
end   

if ~exist('LB'),
    LB = -Inf;
    warning('No lower bound is set. You probably would get better results using LB.')
end

if ~exist('UB'),
    UB = Inf;
    warning('No upper bound is set. You probably would get better results using UB.')
end

[Ny,Nx] = size(y);

z = y;

Y_f = fft2(y);
A_f = conj(H_f).*Y_f;
H2_f = abs(H_f).^2;

% Combination of representations
% rep = [TILs LDCT DTCWT]; % each is a binary value active/not active
TILs = rep(1);
DTCWT = rep(2);
LDCT = rep(3);

alpha1 = alpha_v(1);
alpha2 = alpha_v(2);
alpha3 = alpha_v(3);

sig2r1 = sig2r0_v(1);
sig2r2 = sig2r0_v(2);
sig2r3 = sig2r0_v(3);

sig2r1_f = sig2r_v_floor(1);
sig2r2_f = sig2r_v_floor(2);
sig2r3_f = sig2r_v_floor(3);


Nr = 256; % Reference size for saved spectral profiles (from which any size can be computed)

% Load the PhiPhiT Fourier responses if necessary, and adapt them to the
% image support size

if TILs,
    load PhiPhiTf_TILs % TILs (no LPR, no subsampling)
    PhiPhiT = fftshift(real(ifft2(PhiPhiTf)));
    L = 2; % 3 - 1
    PhiPhiT = PhiPhiT(Nr/2+1-L:Nr/2+1+L,Nr/2+1-L:Nr/2+1+L);
    aux = zeros(Ny,Nx);
    aux(Ny/2+1-L:Ny/2+1+L,Nx/2+1-L:Nx/2+1+L) = PhiPhiT;
    PhiPhiTf_TILs = real(fft2(fftshift(aux)));
end    
    
if DTCWT,
    Nsc = 3; par2 = 'near_sym_a'; par3 = 'qshift_b';
end

if LDCT,
    Nb = 16; % Block size
    rf = 2; % Redundancy factor along each dimension
    ovl = compute_mask_LDCT(y, Nb, rf);
end


for n = 1: Niter,   

    % Sparsify z in the chosen representations

    x_hat = zeros(size(y));

    % TILs
    
    if TILs,

        theta1 = sqrt(sig2r1*2/alpha1);
        [pyr,pind] = buildTDP1(z, 1,0);
        
        pyr = pyr.*(abs(pyr)>=theta1); 
        aux = reconTDP1(pyr, pind)/(2*sig2r1);       
        x_hat = x_hat + aux;
        
    end  
   
    
    % DTCWT
       
    if DTCWT,

        theta2 = sqrt(sig2r2*2/alpha2);
        [Yl,Yh] = dtwavexfm2(z,Nsc,par2,par3);
               
        for s=1:Nsc
            Yh{s} = single(Yh{s}.*(abs(Yh{s})>=theta2));
        end
        Yl = single(Yl.*(abs(Yl)>=theta2));
        aux = dtwaveifm2(Yl,Yh,par2,par3)/(2*sig2r2);
        x_hat = x_hat + aux;
                                      
    end   
    
        
    % LDCT - Local DCT
    
    if LDCT,
        
        theta3 = sqrt(sig2r3*2/alpha3);
        [local_dct,pind] = local_DCT_analysis(z,Nb,rf,ovl);
           
        local_dct = local_dct.*(abs(local_dct)>= theta3);
        aux = local_DCT_synthesis(local_dct,pind,ovl)/(2*sig2r3);
        x_hat = x_hat + aux;
    
    end    
        
    
    % The solution is close to the sparsified z, but also consistent to the
    % observation
    
    X_hat_f = fft2(single(x_hat));
       
    Z_f = (X_hat_f + A_f/(2*sigma^2)) ./ ( TILs*PhiPhiTf_TILs/(2*sig2r1) + DTCWT/(2*sig2r2) + LDCT/(2*sig2r3) + H2_f/(2*sigma^2) );
    
    z = real(ifft2(single(Z_f)));
    z = max(min(z,UB),LB); % Enforces valid pixels to be in [LB UB] range
    
    
    % Prepared for dynamic convergent version (having a "floor" value for
    % sig2r), if not (sig2rj_f = 0) it does not have an impact.
    sig2r1 = beta*(sig2r1 - sig2r1_f) + sig2r1_f;
    sig2r2 = beta*(sig2r2 - sig2r2_f) + sig2r2_f;
    sig2r3 = beta*(sig2r3 - sig2r3_f) + sig2r3_f;

    fprintf(1,'Info:  Deblurring with L2-relaxed-L0, current iter = %d out of %d ...\n',n,Niter);
    
end    
x_hat = z;

