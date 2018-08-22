function params = optimize_params(params0, dynamic, Niter, MaxIter, rep, merit, test_results);

% params = optimize_params(params0, dynamic, Niter, MaxIter, rep, merit, test_results)
%
% Optimize the set of 6 + 1(if dynamic) parameters ({alpha_j,sig2r0j, j = 1..3}, beta) of
% the L2-r-L0 method (TIP 2015).
%
% OUTPUT:
%       params:        Optimized parameters
%
% INPUT
%       params0:      Initial parameters (params = [sig2r01, .., sig2r0N,
%                                                       alpha1, .., alphaN, beta], N<=3, row vector)
%       dynamic:     if dynamic = 1, it assumes Constrained Dynamic optimization (optimizes beta).
%                           if dynamic = 0, it fixes beta = 1 (Static L2-r-L0)
%       Niter:           Number of iterations of the estimation loop. E.g. 10, for ConDy10
%       MaxIter:       Maximum number of iterations allowed in the optimization
%                           loop (e.g., 100).
%       rep:              Binary representation vector, indicating which are
%                           the active representations; rep = [TILs DTCWT LDCT];
%                           e.g., TILs = 1; DTCWT = 1; LDCT = 1;
%       merit:          Merit/cost function. 2 options: 'MSEIR' or 'SSIM'
%                           (other performance criteria can be easily implemented)
%       test_results: If set to 1 it runs the experiments with the image tests using the optimized
%                           parameters, after having optimized the parameters using the 
%                           training set. 1 = Yes, 0 = No.
%
% WARNING: Due to each evaluation of the cost function requires performing 24 deblurring tests,
% optimizing parameters may range from few minutes (example above, TILS, Niter = 5) to several hours
% (e.g., static mode, 3 representations, Niter = 50).
%
%
% Example 1: Optimize parameters for only TILs, Niter = 5, merit = 'MSEIR'.
% The deblurring using this configuration is very fast (e.g. < 0.3 s for 256x256).
% The parameters' optimization itself may take around 10 minutes.
%
%      TILs = 1; DTCWT = 0; LDCT = 0;  % Only TILs active
%       rep = [TILs DTCWT LDCT];
%       Niter = 5;
%       sig2r01 = 5e3; alpha1 = 5; beta = 0.3;
%       params0 = [sig2r01, alpha1, beta];
%       dynamic = 1;
%       MaxIter = 50;
%       merit = 'MSEIR';
%       test_results = 1;
%
% Example 2: Optimize parameters for the three representations, Niter = 10, merit = 'SSIM'.
% The deblurring using this configuration is still pretty fast (around 3 s for
% 256x256), but the parameters' optimization may take several hours
% (typically from 2 to 5 h. depending on your platform).
%
%      TILs = 1; DTCWT = 1; LDCT = 1;  
%      rep = [TILs DTCWT LDCT];
%      Niter = 10;
%      sig2r01 = 5e3; sig2r02 = 5e3; sig2r03 = 5e3; alpha1 = 5; alpha2 = 5; alpha3 = 5; beta = 0.5;
%      params0 = [sig2r01, sig2r02, sig2r03, alpha1, alpha2, alpha3, beta];
%      dynamic = 1;
%      MaxIter = 150;
%      merit = 'SSIM';
%      test_results = 1;
%
%
% Javier Portilla
% Instituto de Optica, CSIC
% Madrid, September 2015



options = optimset('MaxIter',MaxIter);

mode = 'training';
provide_results = 0;

params = params0;

tic;
if dynamic,
    params = fminsearch(@(params) compute_results_deblur(params, Niter, rep, provide_results, mode, merit), params, options)
else
    params = fminsearch(@(params) compute_results_deblur([params 1], Niter, rep, provide_results, mode, merit), params, options)
end    
toc

if test_results,

mode = 'test';
provide_results = 1;
if dynamic,
    [err2, results] = compute_results_deblur(params, Niter, rep, provide_results, mode, merit);
else
    [err2, results] = compute_results_deblur([params 1], Niter, rep, provide_results, mode, merit);
end
if strcmp(merit,'MSEIR'),
    ISNR_table = -10*log10(results.')
elseif strcmp(merit,'SSIM'),
    SSIM_table = 1 - results.'    
end

end