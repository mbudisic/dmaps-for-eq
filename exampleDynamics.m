%% Diffusion Maps as coordinates of Ergodic Partition/Quotient
% In this file we will demonstrate how Diffusion Maps can be used
% to give a set of time-invariant coordinates to ergodic sets in
% the state space of a time-independent dynamical system.
%%
function setname = exampleDynamics

%% Compute or load trajectories from a file
Ngrid = 30; % dimension of grid of initial conditions per axis
Tmax = 10;  % trajectory time length
fwdbwd = 0; % averaging direction -- forward when > 0, backward
            % when < 0, time-symmetric when == 0

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

Npoints = size(xy,3); % number of trajectories
% xy contains trajectories in format
% Nsteps x 2 x Npoints
% so, e.g., xy(:,:,1) is a matrix containing the trajectory of the first
% initial condition

%% Average observables (2D Fourier functions) along trajectories
% Maximum wavevector used in either spatial dimension.
% Number of observables used is K = (2 Wmax + 1)^2
Wmax = 5;   
disp('Generating wavevectors')
[Wx,Wy] = meshgrid(-Wmax:Wmax); 
wv = [Wx(:), Wy(:)].'; % wv is a 2 x K matrix of wavevectors
K = size(wv,2);
D = size(wv,1);

% Average observables along trajectories and store to avgs
avgs = zeros( K, Npoints, 'like', 1+1j );

% Rescaling factors which ensure that the initial harmonic of
% Fourier functions varies sufficiently over the box including all orbits.
xscale = 2*max(max(abs(xy(:,1,:))));
yscale = 2*max(max(abs(xy(:,2,:))));

disp('Computing averages')
for n = 1:Npoints    
    avgs(:,n) = computeAverages( t, xy(:,:,n), wv, ...
                               [xscale, yscale] );
                
end
% avgs is a K x Npoints complex matrix in which each column
% is a vector of averages computed along a single trajectory

%% Compute pairwise distance matrix between trajectories
% Pairwise distance is computed using vectors of observables
% computed in the previous step. It corresponds to a Sobolev
% distance on the space of Fourier coefficients. Matrix of pairwise
% distances is the input to Diffusion Maps algorithm.
%
% Order -( D + 1 )/2 where D is the state dimension was used in
% (Budisic 2012). 
% Order -1/2 was used as a mix norm in (Mathew, 2005)
% Order -1 was used as a mix norm in (Lin, 2011)
disp('Computing distance matrix');
sobolevOrder = -(D + 1)/2;
D2 = sobolevMatrix( avgs, wv, sobolevOrder );

%% Compute Diffusion Coordinates
% Diffusion maps treats the Pairwise Distance matrix as a set of
% discrete samples on a Riemann manifold. Its goal is to extract
% the intrinsic distance on the manifold from the samples by using
% numerical diffusion. The algorithm returns the coordinates of
% samples using the Diffusion Coordinates, in which Euclidean
% distance corresponds to, informally, mean diffusion distance on
% the sampled manifold.
%
% The main parameter of Diffusion Maps algorithm is the diffusion
% bandwidth which governs how strongly will the neighboring points in the
% ergodic quotient be bridged by diffusion. Too small bandwidth will
% result in many disconnected components, too large bandwidth may
% artificially create one large connected component in the ergodic
% quotient.
%
% Heuristically, the calculation does not seem to be too sensitive
% to the choice of bandwidth - changes of even an order of
% magnitude may not influence the final result. The value of the
% bandwidth can be inferred from data to which algorithm is applied
% (Lee, 2009). One of such algorithms are invoked by setting hband
% at 0.
hband = 0;  % diffusion bandwidth - <= 0 to autodetect (see nss.m)

disp('Computing diffusion coordinates');
Ncoord = 10; % Coordinates are sorted, so retain just Ncoord 
[evectors, evalues] = dist2diff(D2, Ncoord, hband); 

% evalues is not really important for visualization
% each column in evectors is Npoints long - elements give diffusion
% coordinates for the corresponding trajectory.

% At this step, we have representation of trajectories in the
% diffusion coordinate space. When averaging length is long enough
% (ergodic averages), and initial conditions cover the state space densely
% enough, these coordinates are the coordinate system for the
% Ergodic Quotient. Everything that follows is post-processing.

%% Parameterize clustering algorithm
% The clustering is performed very crudely: using the signs of the
% three "most independent" coordinate functions (see function
% threeIndependent at the end).  For this reason, using 3 diffusion
% coordinates to cluster will result in at most 8 clusters,
% corresponding to combinations of signs of coordinates at a point,
% e.g., +++, ++-, +-+, etc
% 
% This is performed for demonstration purposes only. In a more
% serious implementation, this clustering step might be omitted to
% retain the fine-grained classification, or replaced with a more
% sophisticated algorithm.
% 
% Choice of diffusion coordinates used for clustering. Ideally,
% coordinates 1, 2, 3 would give the best layout. However, if a
% user wants to leave it open to the algorithm to choose, set NaN
% instead one of the indices.
kvec = threeIndependent(evectors, [1,2,NaN]);

% Number of clusters to split data in:
% [true, false, false] - 2
% [true, true, false] - 4
% [true, true, true] - 8
clustersel = [true, true, false];
assert(any(clustersel),'Set at least one clustering selection');
fprintf('Clustering into %d clusters.\n', ...
        2^sum(double(clustersel)));

disp('Coarse clustering')
zeromean = evectors - repmat( mean(evectors,1), [Npoints,1] );
clusterweight = 2.^[2,1,0].';
clusterweight(~clustersel) = 0;

clusters = round( [sign(zeromean(:,kvec)) + 1] * clusterweight );

%% Visualizing
% There are several ways in which Diffusion Coordinates can be used
% to visualize information about the dynamical system.
%

%%
% *Embedding trajectories in diffusion coordinate space*
% First, we can visualize each trajectory as a point in the space
% of Diffusion Coordinates. In the ergodic limit, this is a
% graphical representation of the quotient of the state space by
% ergodic partition, where each ergodic set is represented by a
% single point.
%
% Points are colored by the value of the stream function at
% the initial condition.
hf = figure(1);
hf.Name = 'Embedding dynamics in three independent diffusion coordinates';
scatter3(evectors(:,kvec(1)), ...
         evectors(:,kvec(2)), ...
         evectors(:,kvec(3)), 5, ...
         streamfunction(ic), 'fill');
xlabel(sprintf('Coordinate %d',kvec(1))); 
ylabel(sprintf('Coordinate %d',kvec(2))); 
zlabel(sprintf('Coordinate %d',kvec(3)));
axis equal
axis square
title({'Embedding dynamics in three independent';'diffusion coordinates';
 'Color is the';'value of stream function'})
colormap(hot)
set(gca,'Color',[0.7,0.7,0.7])

%%
% *Coloring state space using diffusion coordinates*
[X,Y] = meshgrid( icgridX, icgridY );
hf = figure(2);
hf.Name = 'What is this?';
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

%%
% *Coloring state space using signs of diffusion coordinates used for clustering*
hf = figure(3);
hf.Name = 'What is this, a biggun?';
pl = 1;
for n = (kvec+1)
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

%%
% *Coloring state space using clusters*
hf = figure(4);
hf.Name = 'Colored by clusters?'

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
scatter3(evectors(:,kvec(1)), evectors(:,kvec(2)), evectors(:,kvec(3)), 5, clusters, 'fill');
xlabel(sprintf('Evector %d',kvec(1))); ylabel(sprintf('Evector %d',kvec(2))); zlabel(sprintf('Evector %d',kvec(3)));
axis equal
axis square
title({'Embedding into three independent';['coordinates colored by ' ...
                    'clusters']})
set(gca,'Color',repmat(0.7,[1,3]))
colormap(jet)

set(gcf,'color','white');
subtitle('Sign-based clusters');

%% Auxiliaries

function H = streamfunction( ic )
% evaluate the Streamfunction of the f
  icx = ic(1,:);
  icy = ic(2,:);
  k = 1;
  b = 2;
  H = icy(:).^2/2 - k*(icx(:).^2/2 - b*icx(:).^4/4);
  
function [xy, t, icgridX, icgridY, ic] = computeTrajectories(Ngrid, Tmax, fwdbwd)
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


function kvec = threeIndependent( coord, kvec )
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
  Ncoord = size(coord,2);
  
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
    k2 = find( D1 > mean(D1) & (1:Ncoord > k1), 1, 'first' );
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
  
  kvec = [k1, k2, k3];

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

% overlay function to plot on top of color plots
% e.g., simulated trajectories or streamfunction level curves
function overlay(ic, X, Y)
  % evaluate the Streamfunction of the f
  H = streamfunction(ic);

  hold on
  contour(X,Y,reshape(H,size(X)),...
          'Color', 0.7*ones(1,3)); 
  hold off;

