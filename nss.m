function h = nss(D, Nsize)
% h = nss(D, Nneigh)
%
% Estimate bandwidth using Neighborhood Size Stability.
%
% D - square of distances
% Nsize - size of the neighborhood


N = size(D,1);
ds = zeros(1,N);

Nsize = max(Nsize,1); % Nsize cannot be less than 1

% find Nsize-smallest distance
for k = 1:N
    y = sort(D(:,k));
    ds(k) = y(Nsize+1); % +1 accounts for diagonal which is always 0
end

% bandwidth h such that each point has Nsize neighbors within sqrt(2h)
h = max(ds)/2;
