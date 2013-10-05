function [avgs_real, avgs_imag] = computeAverages( t, xy, wv, scales )
%  [avgs_real, avgs_imag] = computeAverages( t, xy, wv, scales )
%
% Evaluate averages of Fourier modes along a single 2D trajectory
%
% t - time vector of Nsteps length
% xy - Nsteps x 2 trajectory
% wv - wavevectors - 2 x K matrix
% scales - size of each coordinate in the state space, assumed symmetric around 0
%        - e.g. if state space is [-2,2] x [-5,5]
%        -      scale should be [4, 10]
%
% Real and imaginary components are computed separately because 
% Matlab Coder (MEX-generator) does not like to return complex variables.
% I don't know why.

K = size(wv, 2);
D = size(wv, 1);

Nsteps = size(xy,1);

assert( D == size(xy,2), 'Dimensions of wavevectors and the trajectory do not match. Wavevectors should be a D x K matrix and the trajectory a Nsteps x D matrix.')

avgs_real = zeros(K, 1); % output is complex valued
avgs_imag = zeros(K, 1); % output is complex valued

x = xy(:,1);
y = xy(:,2);

for k = 1:K

	% compute the argument of the Fourier harmonic
	argument = zeros( Nsteps,1);
	for d = 1:D
		argument = argument + wv(d,k) * xy(:,d)/scales(d);
	end
    
    val_real = cos(2*pi*argument); % fourier harmonic
    val_imag = sin(2*pi*argument); % fourier harmonic
    
    % 1st order integral divided by timespan
    avgs_real(k) = sum( (val_real(1:end-1) + val_real(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );
    avgs_imag(k) = sum( (val_imag(1:end-1) + val_imag(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );
    
end
