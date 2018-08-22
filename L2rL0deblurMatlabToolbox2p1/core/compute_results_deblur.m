function [err2 results] = compute_results_deblur(params, Niter, rep, provide_results, mode, merit, nim0, deg0);

%  [err2 results] = compute_results_deblur(params, Niter, rep, provide_results, mode, merit, nim0, deg0);
%
% Given a set of degradation experiments (a set of blurring kernels
% combined with added Gaussian noise of given variance, on a set of test
% images), a set of deblurring parameters, and a performance criterion,
% compute the global (average) and individual (for each experiment)
% performance of the L2-r-L0 deblurring method.
%
% OUTPUT:
%       err2:                     Global quality measurement.
%                                   For merit = 'MSEIR', err2 is the average Mean
%                                       Square Error Improvement Rate <MSE_out/MSE_in>
%                                   For merit = 'SSIM', err2 is 1 - <SSIM>.
%       results:                Array with the performance measurement for each
%                                   experiment.
%                                   For merit = 'MSEIR',  MSE_out/MSE_in.
%                                   For merit = 'SSIM', 1 - <SSIM>.
% INPUT
%       params:               Model-method parameters (params = [sig2r01, .., sig2r0N,
%                                                       alpha1, .., alphaN, beta], N<=3, row vector)
%       Niter:                  Number of iterations of the estimation loop. E.g. 10, for ConDy10
%       rep:                      Binary representation vector, indicating which representations are active;
%                                   rep = [TILs DTCWT LDCT]; e.g., TILs = 1; DTCWT = 1; LDCT = 1;
%       provide_results:  1 for displaying individual results, 0 for silent mode.  
%       mode:                  'training' for testing the performance on the training images set.
%                                    'test' for measuring the performance  on the test images set.
%       merit:                   Merit/cost function. 2 options: 'MSEIR' or 'SSIM'
%                                   (other performance criteria can be easily implemented)
%       nim0, deg0:         Use this additional parameters for the case of
%                                   wanting to test the deblurring of a specific image (index nim0, from 1 to 3) 
%                                   affected by a specific degradation (index deg0, from 1 to 8).
%
%
% Example: Compute the results on the test images for the SSIM-optimized
% parameters obtained for Niter = 10 (ConDy10).
% Uncomment, copy and paste on the command line the code lines below: 
%       TILs = 1; DTCWT = 1; LDCT = 1; 
%       rep = [TILs DTCWT LDCT];
%       Niter = 10;
%       %Values from Table II (ConDy10, SSIM)
%       alpha1 = 12.03; sig2r01 = 2.26e3; % TILs
%       alpha2 = 18.02; sig2r02 = 6.77e3;  % DTCWT
%       alpha3 = 11.39;  sig2r03 = 4.11e3; % LDCT
%       beta = 0.6200;
%       alpha_v = [alpha1 alpha2 alpha3];
%       sig2r0_v = [sig2r01 sig2r02 sig2r03];
%       params = [sig2r0_v alpha_v beta];
%       mode = 'test';
%       merit = 'SSIM';
%       provide_results = 1;
%       [err2 results] = compute_results_deblur(params, Niter, rep, provide_results, mode, merit);
%       ssim_r = 1 - results
%       SSIM_global = 1 - err2
%
% Javier Portilla
% Instituto de Optica, CSIC
% Madrid, September 2015


if~exist('merit'),
    merit = 'MSEIR'; % Mean Square Error Improvement Ratio
end

% %
% Parameter assignation according to the active representations
% (information coded in "rep")
% %


% rep = [TILs LDCT DTCWT]; % each is a binary value active/not active
TILs = rep(1);
DTCWT = rep(2);
LDCT = rep(3);

% Count nr: number of active representations
nr = DTCWT + TILs + LDCT;

sig2r01 = 1;  % They must be initialized to a non-zero value (regardless of they are used or not).
sig2r02 = 1;
sig2r03 = 1;
alpha1 = 1;
alpha2 = 1;
alpha3 = 1;

nc = 1;
if TILs,
    sig2r01 = abs(params(nc));    % Only positive values
    alpha1 = abs(params(nc+nr));
    nc = nc + 1;
end
if DTCWT,
    sig2r02 = abs(params(nc));
    alpha2 = abs(params(nc+nr));
    nc = nc + 1;
end
if LDCT,
    sig2r03 = abs(params(nc));
    alpha3 = abs(params(nc+nr));
end
if length(params) == 2*nr+1,
    beta = exp(-abs(log(abs(params(2*nr+1))))); % Only values between 0 and 1
elseif length(params) == 2*nr,
    beta = 1;  % static
end    

err2 = 0;
Ny = 256; Nx = 256; % For now

Nim = 3;
Ndeg = 8; % 3x8 = 24 experiments, in total

if exist('nim0'),  % Specific image
    a_im = nim0;
    b_im = nim0;
else
    a_im = 1;
    b_im = Nim;
end    

if exist('deg0'),  % Specific de
    a_deg = deg0;
    b_deg = deg0;
else
    a_deg = 1;
    b_deg = Ndeg;
end    

Nnim = b_im/a_im;   % Effective number of images to process (Nim or 1)
Nndeg = b_deg/a_deg; % Effective number of degradations to process (Ndeg or 1)

results = zeros(Nnim, Nndeg);


% SIMULATE DEGRADATION (BLUR + NOISE)

% Blurring kernel definition

i = meshgrid(-7:7, -7:7); j = i';
kernel = 1./(i.^2 + j.^2 + 1);    % PSF 1
PSF1 = kernel/sum(kernel(:));

kernel = ones(9)/81;                   % PSF 2    
PSF2 = kernel;

kernel = [1 4 6 4 1]'*[1 4 6 4 1]/256; % PSF 3
PSF3 = kernel;

N     = 10;         % PSF4
sigma = 3;
i     = -N:N;
G     = exp( -i.*i./(2*sigma*sigma) );
G     = G(:)./sum(G(:));
PSF4  = G*G';

N     =  4; % PSF5; warning, cropped Gaussian (small 9x9 support considering large sigma)
sigma = 4;
i     = -N:N;
G     = exp( -i.*i./(2*sigma*sigma) );
G     = G(:)./sum(G(:));
PSF5  = G*G';

% It is assumed that current directory is where the Readme.rtf file is

impath_training = 'images/training/';
impath_test = 'images/test/';

for nim = a_im:b_im,

    if strcmp(mode,'training'),
        if nim==1,
            x = imread([impath_training,'house240.png']);
        elseif nim==2,    
            x = imread([impath_training,'pout240.png']);
        elseif nim==3,    
            x = imread([impath_training,'wcop240.png']);
        end
    elseif strcmp(mode,'test'),
        if nim==1,
            x = imread([impath_test,'lena256.png']);
        elseif nim==2,    
            x = imread([impath_test,'barbara256.png']);
        elseif nim==3,    
            x = imread([impath_test,'cameraman.png']);
        end
    else error('Invalid parameter ''mode'': it must be "training" or "test".');
    end    

        x = double(x);
        [Ny,Nx] = size(x);
        
% Degradation experiments' loop
for deg = a_deg:b_deg,            

if deg==1,
    kernel = PSF1;
    sigma = sqrt(2);
elseif deg==2, 
    kernel = PSF1;
    sigma = sqrt(8);
elseif deg==3,
    kernel = PSF2;    
    sigma = 0.56; %sqrt(0.308);
elseif deg==4,    
    kernel = PSF2;
    sigma = sqrt(2);
elseif deg==5,    
    kernel = PSF2;
    sigma = sqrt(4);
elseif deg==6,
    kernel = PSF3;
    sigma = sqrt(8);
elseif deg==7,
    kernel = PSF4;
    sigma = sqrt(4);
elseif deg==8,
    kernel = PSF5;
    sigma = 0.25; 
end    
    
H_f = zeros(Ny,Nx);
kernel = kernel/sum(sum(kernel));
L = (size(kernel,1)-1)/2;
H_f(Ny/2+1-L:Ny/2+1+L,Nx/2+1-L:Nx/2+1+L) = kernel; % Warning: Nx and Ny are assumed to be even numbers
H_f = fft2(fftshift(H_f));

X_f = fft2(x);
randn('seed',0);
Y_f = H_f.*X_f + fft2(sigma*randn(size(x)));
y = real(ifft2(Y_f));


% PERFORM DEBLURRING

PSNR_0 = psnr(y,x);
SSIM_0 = ssim(double(y), double(x));
    
alpha_v = [alpha1 alpha2 alpha3];
sig2r0_v = [sig2r01 sig2r02 sig2r03];

if provide_results,
    tic;
end    
x_hat = deblur_L2relaxedL0(y, rep, H_f, sigma, Niter, alpha_v, sig2r0_v, beta);
if provide_results,
    toc
end    

aux = psnr(x_hat, x);
ISNR = aux - PSNR_0;
if merit(1:3)=='MSE',
    results(min(nim,Nnim),min(deg,Nndeg)) = 10^(-ISNR/10);
    if provide_results,
        ISNR
    end    
elseif merit=='SSIM',
    SSIM = ssim(double(x_hat), double(x));
    results(min(nim,Nnim),min(deg,Nndeg)) = 1 - SSIM;   
    if provide_results,
        SSIM
    end    
end    
    
if provide_results,
    imagesc([255*ones(Ny+2,1) [255*ones(1,2*Nx); double([y x_hat]); zeros(1,2*Nx)]  zeros(Ny+2,1)]); axis off; axis('image'); colormap gray
    drawnow; %pause
end

end  % for deg

end % for nim


err2 = mean2(results)
