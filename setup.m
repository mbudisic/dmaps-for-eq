disp('Setting up Diffusion Maps.');

disp('Compiling MEX file computeAverages.cpp'); 
try
    mex computeAverages.cpp
    disp('MEX compiled successfully');
catch me
    warning(me.message)
    disp('It appears that something is wrong with the MEX setup.')
    disp('Proceeding to use Matlab code');
end

disp('Displaying exampleDynamics.m');
disp('This may take a few minutes if this is your first run.')
disp('Subsequent runs load trajectories from a file.')
web(publish('exampleDynamics'));
close all;