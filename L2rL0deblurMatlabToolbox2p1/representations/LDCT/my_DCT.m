function dct = my_DCT(x, W);

% Fast 1-D DCT computation
%
% Javier Portilla
% Madrid, 2009

% % Compute weights to multiply DFT coefficients
% ww = (exp(-i*(0:N-1)*pi/(2*N))/sqrt(2*N)).';
% ww(1) = ww(1) / sqrt(2);
% ww = 2*ww;  % Double the weights for even-length case  
% W = ww(:,ones(1,N));

[n,m] = size(x);

x = [ x(1:2:n,:); x(n:-2:2,:) ];
dct = fft(x);  

% Multiply FFT by weights:
dct = real(W .* dct);