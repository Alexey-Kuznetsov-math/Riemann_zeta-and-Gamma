clear all;
%**************************************************************
% Test 1: testing ln_gamma(z) againt MATLAB built-in function gamma(z)
z=0.5:0.01:10;
f_exact=log(gamma(z));
f=ln_gamma(z);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 1: the maximum absolute error of ln_gamma(z)-log(gamma(z)) in the range 1/2<z<10 is')
disp(max(abs(f-f_exact)))

%**************************************************************
% Test 2: testing the functional equation ln(Gamma(z+1))=ln(Gamma(z))+ln(z)
N=10000;
z=20*(2*rand(1,N)-1)+20i*(2*rand(1,N)-1);
error=ln_gamma(z+1)-ln_gamma(z)-log(z);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 2: the maximum absolute value of ln_gamma(z+1)-ln_gamma(z)-ln(z)')
disp('for 10 000 random points, uniformly distributed in the square |Im(z)|<20, |Re(z)|<20:')
disp(max(abs(error)))

%**************************************************************
% Test 3: testing the Riemann zeta function against MATLAB built-in function zeta(s)
N=100;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 3: testing the Riemann zeta function against the MATLAB built-in function zeta(s)')
disp(' ')
disp('Test 3(a): computing the relative error Riemann_zeta(s)/zeta(s)-1')
disp('for 100 random points, uniformly distributed in the rectangle |Im(s)|<100, 0<Re(s)<10:')
s=10*rand(1,N)+100i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta(s) for these 100 random points:')
tic
f1=Riemann_zeta(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 100 random points:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 3(b): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1')
disp('for 100 random points, uniformly distributed in the rectangle |Im(s)|<1000, 0<Re(s)<10:')
s=10*rand(1,N)+1000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta(s) for these 100 random points:')
tic
f1=Riemann_zeta(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 100 random points:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 3(c): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1')
disp('for 100 random points, uniformly distributed in the rectangle |Im(s)|<10000, 0<Re(s)<10:')
s=10*rand(1,N)+10000i*(2*rand(1,N)-1);
disp(' ')
disp('Computation time (in seconds) of Riemann_zeta(s) for these 100 random points:')
tic
f1=Riemann_zeta(s);
T1=toc;
disp(T1);
disp(' ')
disp('Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 100 random points:')
tic
f2=zeta(s);
T2=toc;
disp(T2);
disp('the maximum relative error is:')
disp(max(abs(f1./f2-1)))
%**************************************************************
% Test 4: plot the Riemann-Siegel function Z(t) for 0<t<100, 1000<t<1100 and 10000<t<10100
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Test 4: plot the Riemann-Siegel function Z(t) for 0<t<100, 1000<t<1100 and 10000<t<10100')
t=0:0.01:100;
f0=real(Riemann_zeta(0.5+1i*t).*exp(1i*imag(ln_gamma(0.25+0.5i*t))-0.5i*t*log(pi)));
t1=t+1000;
f1=real(Riemann_zeta(0.5+1i*t1).*exp(1i*imag(ln_gamma(0.25+0.5i*t1))-0.5i*t1*log(pi)));
t2=t+10000;
f2=real(Riemann_zeta(0.5+1i*t2).*exp(1i*imag(ln_gamma(0.25+0.5i*t2))-0.5i*t2*log(pi)));
subplot(3,1,1);
plot(t,f0);grid;
title('Z(t) for 0<t<100');
subplot(3,1,2);
plot(t1,f1);grid;
title('Z(t) for 1000<t<1100');
subplot(3,1,3);
plot(t2,f2);grid;
title('Z(t) for 10000<t<10100');



