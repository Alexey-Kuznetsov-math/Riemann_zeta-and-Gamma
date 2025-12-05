clear all;
%**************************************************************
N=1000;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test (a): computing the relative error Riemann_zeta_mex(s)/zeta(s)-1')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:')
s=15*rand(1,N)-5+100i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:')
tic
f1=Riemann_zeta_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test (b): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:')
s=15*rand(1,N)-5+1000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:')
tic
f1=Riemann_zeta_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test (c): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:')
s=15*rand(1,N)-5+10000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:')
tic
f1=Riemann_zeta_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test (d): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100 000, -5<Re(s)<10:')
s=15*rand(1,N)-5+100000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:')
tic
f1=Riemann_zeta_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
