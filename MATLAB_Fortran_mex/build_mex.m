% build_mex.m - Compile the Riemann zeta MEX functions
%
% Place Riemann_zeta_module.f90, Riemann_zeta_mex.F90 and
% Riemann_zeta_prime_mex.F90 in the current directory, then run:
%   >> build_mex

% Use -R2017b for the separated complex API (mxGetPr/mxGetPi)
% This flag is needed for MATLAB R2018a and later
flags = {'-R2017b', 'FOPTIMFLAGS=-O3'};

fprintf('Compiling Riemann_zeta_mex...\n');
mex(flags{:}, '-output', 'Riemann_zeta_mex', ...
    'Riemann_zeta_module.f90', 'Riemann_zeta_mex.F90');

fprintf('Compiling Riemann_zeta_prime_mex...\n');
mex(flags{:}, '-output', 'Riemann_zeta_prime_mex', ...
    'Riemann_zeta_module.f90', 'Riemann_zeta_prime_mex.F90');


