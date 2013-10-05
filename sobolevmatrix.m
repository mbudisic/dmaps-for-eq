function D = sobolevMatrix( V, wv, s )
% SOBOLEVMATRIX( V, wv, s )
%
% Compute squared-Sobolev-distance NxN matrix between 'vectors' of observables of length K.
%
% V - K x N complex matrix - each column is a vector of observables
%           at a different point
% wv - d x K - wavevectors - each column is a wavevector
% s - order of Sobolev distance 

K = size(wv, 2);
N = size(V, 2);

assert( size(V, 1) == K, 'Number of rows in V should match length of wv' );


D = zeros(N);

% compute wavenumbers
wavenumbers = sum(abs(wv),1).';

% compute sobolev weights
weights = (1 + (2*pi*wavenumbers).^2).^s;

% populate distance matrix by Sobolev distance
for i = 1:N
    for j = 1:(i-1)
        dv = sum( abs(V(:,i) - V(:,j)).^2 .* weights );
        D(i,j) = dv;
        D(j,i) = dv;
    end
end


