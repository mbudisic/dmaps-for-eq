function exampleDynamics
% exampleDynamics
%
% Demonstration of use of Diffusion Maps for analysis of dynamical systems.

%% preamble
Ngrid = 10; % dimension of grid of initial conditions per axis
Tmax = 5;   % trajectory time length
Wmax = 3;   % max wavevector used (results in (2 Wmax + 1)^2 observables used
h = -1;     % diffusion bandwidth - <= 0 to autodetect (see nss.m)

%% Compute or load trajectories from a file
demofile = 'exampleDynamicsTrajectories.mat';
if exist(demofile,'file')
    disp(['Loading trajectories. Erase ' demofile ' to recompute.']);
    load(demofile);
else
    disp('Computing trajectories.');
    [xy, t, icgridX, icgridY] = computeTrajectories(Ngrid, Tmax); % parameters set at the beginning of the file
    disp(['Saving trajectories to ' demofile]);    
    save(demofile, 'xy', 't', 'icgridX', 'icgridY');
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
xscale = 2;
yscale = 2;

% select Matlab Coder MEX if it exists
if exist('computeAverages_mex') == 3
    disp('Using MEX averaging function')
    average = @computeAverages_mex;
else
    disp('Using Matlab averaging function. Run "deploytool -build computeAverages.prj" to speed up computation.')
    average = @computeAverages;
end

disp('Computing averages')
if matlabpool('size') < 2
    warning('Open parallel threads by running "matlabpool open" to speed up computation of averages, if desired.')
end
parfor n = 1:Npoints
    [myavg_real, myavg_imag] = average( t, xy(:,:,n), wv, [xscale, yscale] );
    avgs(:,n) = complex(myavg_real, myavg_imag);
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
    distance = @sobolevMatrix;
end

D = distance( avgs, wv, -(spaceDim + 1)/2 );
% D is a Npoints x Npoints real matrix with positive entries

%% compute diffusion coordinates for trajectories
disp('Computing diffusion coordinates');
Nvec = 10; % we need only a few eigenvectors
[evectors, evalues] = dist2diff(D, Nvec, h); % % h is set at the beginning of the file
% evalues is not really important for visualization
% each column in evectors is Npoints long - elements give diffusion coordinates for
% corresponding trajectory

%% THE END OF COMPUTATION
disp('Visualizing')
%% visualization of results (just for purposes of demonstration)
[X,Y] = meshgrid( icgridX, icgridY );
colorind = 1;
color = evectors(:,colorind);

% embed into three dominant eigenvectors and color using original
% parameters
figure('name','Embedding dynamics in three diffusion eigenvectors')
subplot(2,3,1);
scatter3(evectors(:,1), evectors(:,2), evectors(:,3), 5, color, 'fill');
xlabel('Evector 1'); ylabel('Evector 2'); zlabel('Evector 3');
axis equal
axis square
title(sprintf('Color: Diffusion coordinate %d',colorind))
colormap(jet)
set(gca,'Color',[0.5,0.5,0.5])
caxis( [-1,1]*max(abs(color)) );

% color the state space using diffusion coordinates
for n = 2:6
    subplot(2,3,n);
    sel = n-1;
    colorfield = reshape( evectors(:,sel), size(X) );
    pcolor(X, Y, colorfield); shading flat;
    caxis( [-1,1]*max(abs(colorfield(:))) );
    axis square;
    xlabel('x'); ylabel('y');
    title(sprintf('Diffusion evector %d', sel));
    colorbar
end

end

%% Auxiliaries


function [xy, t, icgridX, icgridY] = computeTrajectories(Ngrid, Tmax)
% Compute an ensemble of trajectories started at a uniform grid of 
% points on the state space.

% generate a uniform grid of initial conditions
icgridX = linspace(1/2/Ngrid - 1, 1 - 1/2/Ngrid, Ngrid);
icgridY = icgridX;
[X,Y] = meshgrid(icgridX, icgridY);
ic = [X(:), Y(:)].';
Npoints = size(ic, 2);

% time vector
t = 0:1e-2:Tmax;
Nsteps = numel(t);

xy = zeros( Nsteps, 2, Npoints );

if matlabpool('size') < 2
    warning('Open parallel threads by running "matlabpool open" to speed up computation of trajectories, if desired.')
end
parfor n = 1:Npoints
    p = ic(:, n);
    [~,yout] = ode23t( @vf, t, p );
    xy(:,:,n) = yout;
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

