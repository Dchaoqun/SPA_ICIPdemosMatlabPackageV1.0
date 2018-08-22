function x = my_IDCT(dct,W);

% Fast 1-D inverse DCT computation
%
% Javier Portilla
% Madrid, 2009

[n,m] = size(dct);

% % Compute weights
% ww = sqrt(2*n) * exp(j*(0:n-1)*pi/(2*n)).';
% 
% % Compute precorrection factor
% ww(1) = ww(1)/sqrt(2);
% W = ww(:,ones(1,m));



dct = W.*dct;

% Compute x tilde using equation (5.93) in Jain
y = real(ifft(dct));

% Re-order elements of each column according to equations (5.93) and
  % (5.94) in Jain
x = zeros(n,m);
x(1:2:n,:) = y(1:n/2,:);
x(2:2:n,:) = y(n:-1:n/2+1,:);

