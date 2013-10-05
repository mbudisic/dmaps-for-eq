function [evectors, evalues] = exampleTorus

%% generate points on a 2-torus
disp('Generating points on a torus... ')

% vals will contain values of Fourier harmonics, coord is just for
% visualization
[vals, wv, coord] = generatePoints(50,5);

%% parametrize algorithms
d = 3; % data is parametrized by two parameters
s = -(d+1)/2.; % set to Sobolev norm used for dynamical systems 
h = -1; % negative value for bandwidth turns on NSS estimation
Nvec = 6; % number of diffusion eigenvectors to return eigenvectors

% compute distance matrix
disp('Computing Sobolev distance matrix...')

try
    D = sobolevmatrix_mex( vals.', wv.', s );
catch
    disp('Something is wrong with _mex file -- does not exist? Run "deploytool -build sobolevmatrix.prj" in Matlab prompt.')
    try
        D = sobolevmatrix( vals.', wv.', s );
    catch ME
        disp('Error with .m file as well')
        rethrow(ME)
    end
end

% compute diffusion eigenvectors
disp('Computing diffusion eigenvectors...')
[evectors, evalues] = dist2diff(D, Nvec, h); % h = -1 indicates that the heat bandwidth is to be estimated

%% THE END OF COMPUTATION

%% visualization of results (just for purposes of demonstration)

% embed into three dominant eigenvectors and color using original
% parameters
figure('name','Embedding into three eigenvectors')
viewlist(1) = subplot(2,2,1);
scatter3(evectors(:,1), evectors(:,2), evectors(:,3), 5, coord.phi, 'fill');
xlabel('Evector 1'); ylabel('Evector 2'); zlabel('Evector 3');
axis equal
axis square
title('Color: phi')
colormap(hsv)
set(gca,'Color',[0.5,0.5,0.5])

viewlist(2) = subplot(2,2,2);
scatter3(evectors(:,1), evectors(:,2), evectors(:,3), 5, coord.theta, 'fill');
xlabel('Evector 1'); ylabel('Evector 2'); zlabel('Evector 3');
axis equal
title('Color: theta')
axis square
colormap(hsv)
set(gca,'Color',[0.5,0.5,0.5])

global linkview1;
linkview1 = linkprop(viewlist, {'View'});

viewlist(1) = subplot(2,2,3);
scatter3(evectors(:,1), evectors(:,2), evectors(:,4), 5, coord.phi, 'fill');
xlabel('Evector 1'); ylabel('Evector 2'); zlabel('Evector 4');
axis equal
axis square
title('Color: phi')
colormap(hsv)
set(gca,'Color',[0.5,0.5,0.5])

viewlist(2) = subplot(2,2,4);
scatter3(evectors(:,1), evectors(:,2), evectors(:,4), 5, coord.theta, 'fill');
xlabel('Evector 1'); ylabel('Evector 2'); zlabel('Evector 4');
axis equal
title('Color: theta')
axis square
colormap(hsv)
set(gca,'Color',[0.5,0.5,0.5])

global linkview2;
linkview2 = linkprop(viewlist, {'View'});

% embed into  (theta,phi) original coordinates and color using eigenvectors
figure('name', 'evector scatter color')
views = zeros(1,Nvec);
for k = 1:Nvec
    views(k) = subplot(3,ceil(Nvec/3),k);
    scatter(coord.theta, coord.phi, 5, evectors(:,k),'fill');
    title(sprintf('Evector %d',k))
    xlabel('Theta'); ylabel('Phi');
    axis equal
    axis square
    axis tight
end

global linkview3;
linkview3 = linkprop(views, {'View'});

% sort points according to one of the first four eigenvectors
% and visualize other eigenvectors in that order
figure('name', 'eigenvector re-sorted plots');
[~,id1] = sort(evectors(:,1),1,'ascend');
[~,id2] = sort(evectors(:,2),1,'ascend');
[~,id3] = sort(evectors(:,3),1,'ascend');
[~,id4] = sort(evectors(:,4),1,'ascend');


for k = 1:Nvec
    subplot(4,Nvec,k)
    plot(evectors(id1,1), evectors(id1,k), '.', 'MarkerSize',1);
    title(sprintf('Evector %d',k))
    xlabel('Eigenvector 1')
    
    subplot(4,Nvec,k+Nvec)
    plot(evectors(id2,2), evectors(id2,k), '.', 'MarkerSize',1);
    title(sprintf('Evector %d',k))
    xlabel('Eigenvector 2')

    subplot(4,Nvec,k+2*Nvec)
    plot(evectors(id3,3), evectors(id3,k), '.', 'MarkerSize',1);
    title(sprintf('Evector %d',k))
    xlabel('Eigenvector 3')
    
    
    subplot(4,Nvec,k+3*Nvec)
    plot(evectors(id4,4), evectors(id4,k), '.', 'MarkerSize',1);
    title(sprintf('Evector %d',k))
    xlabel('Eigenvector 4')
    
end

function [retval, wvs, coord] = generatePoints(N, K)
% generate Np = NxN points on a 2-torus
% and evaluate Nk = (K+1) x (2K+1) 2d Fourier harmonics on them
%
% returns:
% retval (Np, Nk) - each row is a vector of evaluations (corresponding to a
%                   different point
% wvs (Nk, 2)     - each row is a different wavevector
% coord           - coordinates of points

ax1 = linspace(0,1,N+1);
ax2 = linspace(0,1,N+1);

[th, phi] = meshgrid( ax1(2:end), ax2(2:end) );
points = [th(:), phi(:)];
Np = size(points, 1)

% generate K+1 x 2K+1 harmonics (this is enough due to data being realc)
[kx, ky] = meshgrid( -K:K, -K:K );
wvs = [kx(:), ky(:)];
Nk = size(wvs, 1)

% evaluate 2d Fourier harmonics on given points
harmonic = @(x,y, kx, ky) exp( 2j*pi*(kx.*x + ky.*y) );

% compute evaluations of harmonics
retval = zeros(Np, Nk);
for i = 1:Nk
    retval(:,i) = harmonic( points(:,1), points(:,2), wvs(i,1), wvs(i,2) );
end

[x,y,z] = tor2cart( points(:,1),points(:,2), 5, 1);
coord.x = x; coord.y = y; coord.z = z; coord.theta = points(:,1); coord.phi = points(:,2);

function [x,y,z] = tor2cart( theta,phi, Rmajor, aminor)
% convert from torus to cartesian coordinates

x = (Rmajor + aminor.*cos(phi)) .* cos(theta);
y = (Rmajor + aminor.*cos(phi)) .* sin(theta);
z = aminor.*sin(phi);



