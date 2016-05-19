function D2 = sobolevMatrix( V, wv, order )
% SOBOLEVMATRIX( V, wv, order )
%
% Compute squared-Sobolev-distance NxN matrix between 'vectors' of observables of length K.
%
% V - K x N complex matrix - each column corresponds to a different point
% wv - D x K - wavevectors - each column is a D-dimensional wavevector
% order - order of Sobolev distance 

import DiffusionMaps.*


D = size(wv, 1); % number of states
K = size(wv, 2); % number of observables

validateattributes( V, {'numeric'}, {'nrows', K} );
N = size(V, 2); % number of points

% populate distance matrix by Sobolev distance - result is a vector
% storing just the lower triangle of the matrix
D2vec = pdist( V.', @(Xi, Xj)sobdist2( Xi, Xj, wv, order ) );
% this should be simplified using 'mahalanobis' in pdist

assert( ~any( isnan(D2vec(:))), 'There is a NaN in the matrix');

% create square matrix from the lower triangle
D2 = squareform(D2vec);

% clear persistent variables in sobdist2
sobdist2;

end

function d2 = sobdist2( onerow, multirow, wv, order )
% SOBDIST2
%
% Computes row-wise Sobolev distance squared between
% 'onerow' and each of the rows in 'multirow'. Sobolev distance is
% really a weighted Euclidean distance with weights computed based
% on wavenumbers.
%
% wv and order are used to compute weights of the sobolev distance.
%
% SOBDIST2 stores weights in an internal persistent variable
% Call it without arguments to clear the internal state.

persistent weights;

if nargin < 1
  clear weights
  return
end

% compute Sobolev matrix weights on the first invocation of the function
if isempty(weights)
  wavenumbers = sum(abs(wv),1);
  weights = (1 + (2*pi*wavenumbers).^2).^order;
  weights = weights(:);
end

% vectorized weighted euclidean distance squared
Nrows = size(multirow, 1);
d2 = abs( repmat( onerow, [Nrows,1] ) - multirow ).^2 * weights;
end
