function viewHarmonic(wv, xscale, yscale)
% viewHarmonic(wv, xscale, yscale)
%
% Helper function that visualizes Fourier harmonic on 2 plane.
%
% wv - 2 x 1 wavevector
% xscale - max value of x
% yscale - max value of y

% set up the domain
[X,Y] = meshgrid(-1:1e-2:1);
X = xscale*X;
Y = yscale*Y;

% compute the values
wvx = wv(1);
wvy = wv(2);
val = exp(1j*pi*(wvx*X/xscale + wvy*Y/yscale)); % fourier harmoni

%visuslize
subplot(1,2,1);
pcolor( X,Y, real(val)); shading flat;
title('Real');

subplot(1,2,2);
pcolor( X,Y, imag(val)); shading flat;
title('Imaginary')
