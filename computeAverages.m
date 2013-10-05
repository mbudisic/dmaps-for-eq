function [avgs_real, avgs_imag] = computeAverages( t, xy, wv, xscale, yscale )
%  [avgs_real, avgs_imag] = computeAveragesSingle( t, xy, wv, xscale, yscale )
%
% Evaluate averages of Fourier modes along a single 2D trajectory
%
% t - time vector of Nsteps length
% xy - Nsteps x 2 trajectory
% wv - wavevectors - 2 x K matrix
% xscale - width of the state space domain (assumed symmetric around 0)
% yscale - height of the state space domain (assumed symmetric around 0)
%
% Real and imaginary components are computed separately because 
% Matlab Coder (MEX-generator) does not like to return complex variables.
% I don't know why.

K = size(wv, 2);
Nsteps = size(xy,1);
avgs_real = zeros(K, 1); % output is complex valued
avgs_imag = zeros(K, 1); % output is complex valued

x = xy(:,1);
y = xy(:,2);

for k = 1:K
    
    wvx = wv(1,k);
    wvy = wv(2,k);

    val_real = cos(pi*(wvx*x/xscale + wvy*y/yscale)); % fourier harmonic
    val_imag = sin(pi*(wvx*x/xscale + wvy*y/yscale)); % fourier harmonic
    
    % 1st order integral divided by timespan
    avgs_real(k) = sum( (val_real(1:end-1) + val_real(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );
    avgs_imag(k) = sum( (val_imag(1:end-1) + val_imag(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );
    
end
