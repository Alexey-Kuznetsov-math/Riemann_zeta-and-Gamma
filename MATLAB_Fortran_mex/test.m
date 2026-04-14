function test
%**************************************************************
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_prime_mex function: computing the first non-trivial zero with Newtons method')
% testing the function Riemann_zeta_prime_mex
s1=1/2+14.13472514173469379i;% the first non-trivial zero of zeta(s), downloaded from https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros2
s=1/2+1i*14.1;% initial approximation to this zero
for i=1:5
    [f,f1]=Riemann_zeta_prime_mex(s);
    s=s-f/f1 % use Newton's method to find better approximation
end
error=s-s1
disp('Testing the Riemann_zeta_prime_mex function: computing the non-trivial zero close to 1/2+401.8i with Newtons method')
s1=1/2+401.8392286005332165399113i; 
s=1/2+1i*401.8;% initial approximation to this zero
for i=1:5
    [f,f1]=Riemann_zeta_prime_mex(s);
    s=s-f/f1 % use Newton's method to find better approximation
end
error=s-s1
%**************************************************************
N=1000;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test (a): comparing the values of Riemann_zeta_mex(s) and zeta(s)')
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
disp('the maximum mixed absolute/relative error is:')
disp(mixed_abs_rel_error(f1,f2))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test (b): comparing the values of Riemann_zeta_mex(s) and zeta(s)')
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
disp('the maximum mixed absolute/relative error is:')
disp(mixed_abs_rel_error(f1,f2))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test (c): comparing the values of Riemann_zeta_mex(s) and zeta(s)')
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
disp('the maximum mixed absolute/relative error is:')
disp(mixed_abs_rel_error(f1,f2))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test (d): comparing the values of Riemann_zeta_mex(s) and zeta(s)')
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
disp('the maximum mixed absolute/relative error is:')
disp(mixed_abs_rel_error(f1,f2))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)')
disp(' ')
disp('Test (e): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:')
s=15*rand(1,N)-5+100i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:')
tic
[f,f1]=Riemann_zeta_prime_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:')
tic
g=zeta(s);
g1=zeta(1,s);
T2=toc;
disp(T2);
disp('the maximum mixed absolute/relative error for zeta(s) is:')
disp(mixed_abs_rel_error(f,g))
disp('the maximum mixed absolute/relative error for derivative of zeta(s) is:')
disp(mixed_abs_rel_error(f1,g1))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)')
disp(' ')
disp('Test (f): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:')
s=15*rand(1,N)-5+1000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:')
tic
[f,f1]=Riemann_zeta_prime_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:')
tic
g=zeta(s);
g1=zeta(1,s);
T2=toc;
disp(T2);
disp('the maximum mixed absolute/relative error for zeta(s) is:')
disp(mixed_abs_rel_error(f,g))
disp('the maximum mixed absolute/relative error for derivative of zeta(s) is:')
disp(mixed_abs_rel_error(f1,g1))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)')
disp(' ')
disp('Test (g): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)')
disp('for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:')
s=15*rand(1,N)-5+10000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:')
tic
[f,f1]=Riemann_zeta_prime_mex(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:')
tic
g=zeta(s);
g1=zeta(1,s);
T2=toc;
disp(T2);
disp('the maximum mixed absolute/relative error for zeta(s) is:')
disp(mixed_abs_rel_error(f,g))
disp('the maximum mixed absolute/relative error for derivative of zeta(s) is:')
disp(mixed_abs_rel_error(f1,g1))
%**************************************************************
function e=mixed_abs_rel_error(f,g)
% compute the maximum componentwise mixed absolute/relative error between 
% f=[f_1,...,f_n] and g=[g_1,...,g_n] 
% we define e=max(h_1,...,h_n), where h_i=|f_i-g_i| if |g_i|<=1 
% and h_i=|f_i/g_i-1| if |g_i|>1 
i1=find(abs(g)<=1);
i2=find(abs(g)>1);
e=max([abs(f(i1)-g(i1)), abs(f(i2)/g(i2)-1)]);
