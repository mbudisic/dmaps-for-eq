function setname = exampleDynamics
% exampleDynamics
%
% Demonstration of use of Diffusion Maps for analysis of dynamical systems.

%% preamble
Ngrid = 30; % dimension of grid of initial conditions per axis
Tmax = 10;   % trajectory time length
Wmax = 5;   % max wavevector used (results in (2 Wmax + 1)^2 observables used
hband = 0;  % diffusion bandwidth - <= 0 to autodetect (see nss.m)
fwdbwd = 0; % averaging direction -- forward when > 0, backward
            % when < 0, time-symmetric when == 0

% if you want to force the use of a certain vector for clustering, 
% enter its index in k1, k2, or k3
% otherwise leave as NaN
k1 = 1;
k2 = 2;
k3 = NaN;
kvec = [k1,k2,k3];
clustersel = [true, true, false];

% cluster into 2 (one true element), 4 (two true elements), or 8
% clusters (three true elements)
assert(any(clustersel),'Set at least one clustering selection');
fprintf('Clustering into %d clusters.\n', 2^sum(double(clustersel)));


%% Compute or load trajectories from a file
demofile = sprintf('exampleDynamicsTrajectories_T%.1f.mat',Tmax);
setname = demofile(1:end-4);
if exist(demofile,'file')
    disp(['Loading trajectories. Erase ' demofile ' to recompute.']);
    load(demofile);
else
    disp('Computing trajectories.');
    [xy, t, icgridX, icgridY,ic] = computeTrajectories(Ngrid, Tmax, 0);
    % parameters set at the beginning of the file
    disp(['Saving trajectories to ' demofile]);    
    save(demofile, 'xy', 't', 'icgridX', 'icgridY', 'ic');
end
Npoints = size(xy,3);
% xy contains trajectories in format
% Nsteps x 2 x Npoints
% so, e.g., xy(:,:,1) is a matrix containing the trajectory of the first
% initial condition

%% generate all relevant wavevector pairs up to Wmax harmonic in each dimension
disp('Generating wavevectors')
[Wx,Wy] = meshgrid(-Wmax:Wmax); % Wmax is set at the beginning
wv = [Wx(:), Wy(:)].';
% wv is a 2 x K matrix of wavevectors -
% K = (2*Wmax+1)^2
K = size(wv,2);

%% compute averages of Fourier functions along trajectories
avgs = zeros( K, Npoints, 'like', 1+1j );

xscale = 2*max(max(abs(xy(:,1,:))));
yscale = 2*max(max(abs(xy(:,2,:))));

% select Matlab Coder MEX if it exists

disp('Computing averages')
tic;
for n = 1:Npoints
    avgs(:,n) = computeAverages( t, xy(:,:,n), wv, ...
                               [xscale, yscale] );
            
end

% avgs is a K x Npoints complex matrix in which each column
% is a vector of averages computed along a single trajectory

%% compute sobolev distances between trajectories
disp('Computing distance matrix');
spaceDim = 2;

if exist('sobolevMatrix_mex') == 3
    disp('Using MEX distance function')
    distance = @sobolevMatrix_mex;
else
    disp('Using Matlab distance function. Run "deploytool -build sobolevMatrix.prj" to speed up computation.')
    distance = @sobolevmatrix;
end

D = distance( avgs, wv, -(spaceDim + 1)/2 );
% D is a Npoints x Npoints real matrix with positive entries

%% compute diffusion coordinates for trajectories
disp('Computing diffusion coordinates');
Nvec = 10; % we need only a few eigenvectors
[evectors, evalues] = dist2diff(D, Nvec, hband); % % h is set at the beginning of the file
% evalues is not really important for visualization
% each column in evectors is Npoints long - elements give diffusion
% coordinates for the corresponding trajectory.

% Heuristic: find indices of three "most independent" coordinates.
[k1,k2,k3] = threeIndependent(evectors, [k1, k2, k3]);

%% THE END OF COMPUTATION
disp('Visualizing')
%% visualization of results (just for purposes of demonstration)
[X,Y] = meshgrid( icgridX, icgridY );

% evaluate the Hamiltonian of the f
icx = ic(1,:);
icy = ic(2,:);
k = 1;
b = 2;
H = icy.^2/2 - k*(icx.^2/2 - b*icx.^4/4);

% embed into three dominant eigenvectors and color using original
% parameters
figure('name','Embedding dynamics in three independent diffusion coordinates')
scatter3(evectors(:,k1), evectors(:,k2), evectors(:,k3), 5, H.', 'fill');
xlabel(sprintf('Coordinate %d',k1)); 
ylabel(sprintf('Coordinate %d',k2)); 
zlabel(sprintf('Coordinate %d',k3));
set(gcf,'Color','white');

axis equal
axis square
title({'Embedding dynamics in three independent';'diffusion coordinates';
 'Color is the';'value of stream function'})

colormap(hot)
set(gca,'Color',[0.7,0.7,0.7])
%caxis( [-1,1]*max(abs(H)) );
a1 = gca;
% color the state space using diffusion coordinates
figure
for n = 2:10
    subplot(3,3,n-1);
    sel = n-1;
    colorfield = reshape( evectors(:,sel), size(X) );
    pcolor(X, Y, colorfield); shading flat;
    caxis( [-1,1]*max(abs(colorfield(:))) );
    axis square;
    xlabel('x'); ylabel('y');
    title(sprintf('Diffusion coordinate %d', sel));
    overlay(ic, X, Y);    
    colorbar
end
set(gcf,'color','white');
subtitle('Coloring of the state space by diffusion coordinates');

figure
pl = 1;
for n = ([k1,k2,k3]+1)
    subplot(1,3,pl);

    sel = n-1;
    colorfield = reshape( evectors(:,sel), size(X) );
    pcolor(X, Y, sign(colorfield)); shading flat;
    axis square;
    xlabel('x'); ylabel('y');
    overlay(ic, X, Y);        
    title(sprintf('Diff. coordinate (%d)', sel));
    pl = pl+1;    
    %colorbar
end
set(gcf,'color','white');
subtitle('State space colored by signs of independent coordinates');



%% Clustering

% The clustering is performed very crudely: using the signs of the
% three "most independent" coordinate functions (see function
% threeIndependent at the end).
%
% For this reason, there are always 8 clusters, corresponding to
% combinations of signs of coordinates k1 k2 k3, e.g., +++, ++-,
% +-+, etc
% 
% This is performed for demonstration purposes only. In a more
% serious implementation, this clustering step might be omitted to
% retain the fine-grained classification, or replaced with a more
% sophisticated algorithm.
figure;

disp('Coarse clustering')
zeromean = evectors - repmat( mean(evectors,1), [Npoints,1] );
clusterweight = 2.^[2,1,0].';
clusterweight(~clustersel) = 0;

clusters = round( [sign(zeromean(:,[k1,k2,k3])) + 1] * clusterweight );

subplot(1,2,1)
colorfield = reshape( clusters, size(X) );
pcolor(X, Y, colorfield); shading flat;
axis square;
xlabel('x'); ylabel('y');
title('Initial conditions labeled by clusters');
overlay(ic, X, Y);
colormap(jet)
colorbar

subplot(1,2,2)
scatter3(evectors(:,k1), evectors(:,k2), evectors(:,k3), 5, clusters, 'fill');
xlabel(sprintf('Evector %d',k1)); ylabel(sprintf('Evector %d',k2)); zlabel(sprintf('Evector %d',k3));
axis equal
axis square
title({'Embedding into three independent';['coordinates colored by ' ...
                    'clusters']})
set(gca,'Color',repmat(0.7,[1,3]))
colormap(jet)

set(gcf,'color','white');
subtitle('Sign-based clusters');

end
%% Auxiliaries


function [xy, t, icgridX, icgridY, ic] = computeTrajectories(Ngrid, Tmax, ...
                                                    fwdbwd)
% Compute an ensemble of trajectories started at a uniform grid of 
% points on the state space.

% generate a uniform grid of initial conditions
icgridX = linspace(1/2/Ngrid - 1, 1 - 1/2/Ngrid, Ngrid);
icgridY = icgridX;
[X,Y] = meshgrid(icgridX, icgridY);
ic = [X(:), Y(:)].';
Npoints = size(ic, 2);

% time vector
tpos = [];
tneg = [];
if fwdbwd == 0
  disp('Symmetric time average')
  tpos = 0:1e-2:(Tmax/2);
  tneg = -( 0:1e-2:(Tmax/2) );
  t = [tneg(end:-1:2), tpos];
elseif fwdbwd > 0
  disp('Forward average')  
  tpos = 0:1e-2:Tmax;
  tneg = [];
  t = tpos;
else
  disp('Backward average')    
  tpos = [];
  tneg = -( 0:1e-2:(Tmax/2) );  
  t = tneg(end:-1:2);
end

Nsteps = numel(t);
xy = zeros( Nsteps, 2, Npoints );

for n = 1:Npoints
    p = ic(:, n);
    if ~isempty(tpos) && isempty(tneg)
      [~,ypos] = ode23t( @vf, tpos, p );
      XYt = ypos;      
    elseif isempty(tpos) && ~isempty(tneg)    
      [~,yneg] = ode23t( @vf, tneg, p );
      XYt = yneg;
    else
      [~,ypos] = ode23t( @vf, tpos, p );  
      [~,yneg] = ode23t( @vf, tneg, p );      
      XYt = [yneg(end:-1:2,:); ypos];      
    end

    xy(:,:,n) = XYt;
end



end

function u = vf(t,p)
% vector field -- double well potential
k = 1;
b = 2;

u = zeros(2, numel(t) );

% H(x,y) = y^2/2 - k(x^2/2 - b*x^4/4)
x = p(1,:);
y = p(2,:);

u(1,:) = -y;
u(2,:) = -k*(x - b*x.^3);
end

function [k1,k2,k3] = threeIndependent( coord, kvec )
%
% Retrieve indices of three "most independent" coordinate functions
% (each column is an evaluation of a coordinate function on the
% points)
%
% This is a heuristic, but a sensible one in the context of diffusion
% maps.  The idea behind the heuristic is to extract coordinate
% functions which vary along independent spatial directions, e.g., for
% planar Fourier functions, k1 could be the first harmonic in x
% direction k2 first harmonic in y direction, and k3 the first mixed
% xy harmonic.
%
% k1 is always one, as this is assumed to be the
% "large scale" mode
%
% Second coordinate is the first coordinate that shows a large
% variation in gradient against k1(larger than mean of all the
% coordinates)
% 
% Third coordinate is the coordinate (different from k1 and k2)
% that has the largest variation in gradient against both k1 and
% k2.
%
  
  Npoints = size(coord,1);
  Nvec = size(coord,2);
  
  if isnan( kvec(1) )
  k1 = 1;
  else
    k1 = kvec(1);
  end
    
  sorted1 = sortrows(coord,k1);
  D1 = sum(abs(diff(sorted1,1,1))/Npoints,1);
  sorted1 = sortrows(coord,k1);
  D1 = sum(abs(diff(sorted1,1,1))/Npoints,1);
  
  if isnan( kvec(2) )
    k2 = find( D1 > mean(D1) & (1:Nvec > k1), 1, 'first' );
  else
    k2 = kvec(2);
  end
  

  sorted2 = sortrows(coord,k2);
  D2 = sum(abs(diff(sorted2,1,1))/Npoints,1);

  D2(k1) = -inf;
  D2(k2) = -inf;
  
  if isnan(kvec(3))
    [~,k3] = max( D1 + D2 );
  else
    k3 = kvec(3);
  end

end

function [ax,h]=subtitle(text)
%
%Centers a title over a group of subplots.
%Returns a handle to the title and the handle to an axis.
% [ax,h]=subtitle(text)
%           returns handles to both the axis and the title.
% ax=subtitle(text)
%           returns a handle to the axis only.
ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
set(get(ax,'Title'),'Visible','on');
uistack(ax,'bottom') % push title axis to the bottom
title(text);
if (nargout < 2)
    return
end
h=get(ax,'Title');
end

% overlay function to plot on top of color plots
% e.g., simulated trajectories or hamiltonian level curves
function overlay(ic, X, Y)
  hold on
  % evaluate the Hamiltonian of the f
  icx = ic(1,:);
  icy = ic(2,:);
  k = 1;
  b = 2;
  H = icy.^2/2 - k*(icx.^2/2 - b*icx.^4/4);
  hold on; 
  contour(X,Y,reshape(H,size(X)),...
          'Color', 0.7*ones(1,3)); 
  hold off;
end
