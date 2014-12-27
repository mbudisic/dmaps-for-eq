function avgs = computeAverages_mat( t, xy, wv, scales )
%  [avgs_real, avgs_imag] = computeAverages( t, xy, wv, scales )
%
% Evaluate averages of Fourier modes along a single 2D trajectory
%
% t - time vector of Nsteps length
% xy - Nsteps x 2 trajectory
% wv - wavevectors - K x 2 matrix
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

avgs = zeros(K, 1, 'like',1+1j); % output is complex valued

for k = 1:K

	% compute the argument of the Fourier harmonic
	argument = zeros( Nsteps,1);
	for d = 1:D
		argument = argument + wv(d,k) * xy(:,d)/scales(d);
	end
    
    val = exp(2j*pi*argument); % fourier harmonic
    
    % 1st order integral divided by timespan
    avgs(k) = sum( (val(1:end-1) + val(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );
    
end

