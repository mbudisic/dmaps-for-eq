function avgs = computeAverages( t, xy, wv, domain )
%  [avgs] = computeAverages( t, xy, wv, domain )
%
% Evaluate averages of Fourier modes along a single trajectory in a
% D-dimensional state space.
%
% Inputs:
% t - time vector of Nsteps length
% xy - Nsteps x D trajectory - each row is a point on trajectory
% wv - wavevectors - D x K matrix  - each column is a D-dimensional wavevector
% domain - width of the domain in each state variable ( 1 x D )
%
% Returns:
% avgs - K x 1 vector of averages
%
% Function invokes either MEX version of code if it is compiled.

D = size(wv, 1);
K = size(wv, 2);

Nsteps = size(xy,1);

assert( D == size(xy,2), ['Dimensions of wavevectors and the trajectory ' ...
                    'do not match. Wavevectors should be a D x K ' ...
                    'matrix and the trajectory a Nsteps x D matrix.'])

avgs = zeros(K, 1, 'like',1+1j); % output is complex valued

for k = 1:K

	% compute the argument of the Fourier harmonic
	argument = zeros( Nsteps,1);
	for d = 1:D
		argument = argument + wv(d,k) * xy(:,d)/domain(d);
	end

    val = exp(2j*pi*argument); % fourier harmonic

    % 1st order integral divided by timespan
    avgs(k) = sum( (val(1:end-1) + val(2:end) )/2 .* diff(t(:)) ) / ( t(end) - t(1) );

end
