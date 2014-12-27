function avgs = computeAverages( t, xy, wv, scales )
%  [avgs_real, avgs_imag] = computeAverages( t, xy, wv, scales )
%
% Evaluate averages of Fourier modes along a single 2D trajectory
%
% t - time vector of Nsteps length
% xy - Nsteps x 2 trajectory
% wv - wavevectors - K x 2 matrix
% scales - size of each coordinate in the state space, assumed symmetric around 0
%        - e.g. if state space is [-2,2] x [-5,5]
%        -      scale should be [4, 10]
%
% Real and imaginary components are computed separately because 
% Matlab Coder (MEX-generator) does not like to return complex variables.
% I don't know why.

global DIFFMAPS_nomex
global DIFFMAPS_test

if isempty(DIFFMAPS_test)
    DIFFMAPS_test = false;
end
if isempty(DIFFMAPS_nomex)
    DIFFMAPS_nomex = false;
end

if exist('computeAverages_mex') == 3 && DIFFMAPS_test
    tic
    avgs_c = computeAverages_mex( t, xy, wv, scales );
    fprintf(1, 'Avgs using MEX: %.3f msec\n', toc*1000);
    
    tic
    avgs_m = computeAverages_mat( t, xy, wv, scales );
    fprintf(1, 'Avgs using MAT: %.3f msec\n', toc*1000);
    fprintf(1, 'Error: %.3f \n', log10(max(abs(avgs_m(:) - avgs_c(:)))) );
end


if exist('computeAverages_mex') == 3 && ~DIFFMAPS_nomex
    avgs = computeAverages_mex( t, xy, wv, scales );
else
    avgs = computeAverages_mat( t, xy, wv, scales );
end
    
