function [evectors, evalues] = dist2diff(D, Nvec, h)
% DIST2DIFF(D, Nvec, h)
%
% Compute diffusion eigenfunctions from a distance matrix.
%
% D - squared-distance matrix
% Nvec - number of vectors computed
% h - diffusion bandwidth (if omitted or negative, it will be estimated)
%
% returns:
% evectors - eigenvectors of the diffusion (without the trivial
%            eigenvector)
% evalues - eigenvalues of the diffusion (without the trivial one)

% number of points
N = size(D,1);

%% estimate diffusion bandwidth
if ~exist('h', 'var') || h < 0
    Nsize = min( fix(N*5e-3), N-1); % neighborhood size - this is heuristic, sets to .5% of dataset
    h = nss(D, Nsize); % estimate using neighborhood size stability
    fprintf(1,'Estimated bandwidth: %.2e\n',h);
end


%% unbiased heat kernel
A = exp( -D/(4*h)); % heat kernel evaluation

p = sum(A,1);    % estimated sampling density
Ahat = A ./ (p.' * p); % remove sampling density bias from heat kernel


%% convert unbiased heat kernel to a symmetrized Markov chain matrix
symmass_col = sqrt(sum(Ahat,1));
symmass_row = sqrt(sum(Ahat,2));

scaling = (symmass_row * symmass_col);

assert( size(scaling,1) == size(Ahat,1) );
assert( size(scaling,2) == size(Ahat,2) );

S = Ahat ./ scaling;
assert( norm(S - S.') < 1e-16, 'Symmetrized matrix not symmetric!')

%% compute eigenvectors of the heat Markov chain
opts.issym = true;
opts.isreal = true;
opts.disp = 0;
opts.v0 = [1; zeros(N-1,1)]; % remove randomness by fixing initial vector
[V, Sigma] = eigs(S, Nvec+1, 'LA',opts); % add +1 to account for 1st trivial evector


%% make sure vectors are sorted by descending eigenvalue order
evalues = diag(Sigma);
[evalues, eind] = sort(evalues,'descend');
evectors = V(:, eind);

%% rescale with dominant eigenvector to account for symmetrization
evectors = evectors ./ repmat(evectors(:,1), [1, size(evectors,2)]);

