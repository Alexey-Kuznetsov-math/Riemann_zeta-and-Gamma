function f=Riemann_zeta(s)
% Riemann_zeta approximates the Riemann zeta function in the entire complex plane
%
%   f = Riemann_zeta(s) returns an approximation to zeta(s) for complex input s
%   The input s can be a scalar, a vector or an array
%
% This function requires ln_gamma function (which computes the logarithm of the Gamma function 
% in the entire complex plane)
% -------------------------------------------------------------------------
% Author: Alexey Kuznetsov
% York University, Toronto, Canada
% Website: https://kuznetsovmath.ca/
% Email: akuznets@yorku.ca
%
% Created: 18-Nov-2025
% Last updated: 5-Dec-2025
%
% License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
%--------------------------------------------------------------------------
    i1=find(real(s)>=0.5);
    i2=find((real(s)<0.5)&(imag(s)>=0));
    i3=find((real(s)<0.5)&(imag(s)<0));
    f=zeros(size(s));
    if (~isempty(i1)) 
        f(i1)=Riemann_zeta_half_plane(s(i1));
    end
    if (~isempty(i2)) % use reflection formula if Re(s)<0.5
        w2=s(i2);
        f(i2)=1i*(2*pi).^(w2-1).*(1-exp(1i*pi*w2)).*exp(-0.5i*pi*w2+ln_gamma(1-w2)).*Riemann_zeta_half_plane(1-w2);
    end
    if (~isempty(i3)) % use reflection formula if Re(s)<0.5
        w3=s(i3);
        f(i3)=-1i*(2*pi).^(w3-1).*(1-exp(-1i*pi*w3)).*exp(0.5i*pi*w3+ln_gamma(1-w3)).*Riemann_zeta_half_plane(1-w3);
    end
end
%##########################################################################
function f=Riemann_zeta_half_plane(s)
% approximate zeta(s) in the half-plane Re(s)>0.5
%--------------------------------------------------------------------------
    i1=find(real(s)>=5);
    i2=find((real(s)<5)&(abs(imag(s))<=200));
    i3=find((real(s)<5)&(imag(s)>200));
    i4=find((real(s)<5)&(imag(s)<-200)); 
    f=zeros(size(s));
    if (~isempty(i1)) % in the half-plane Re(s)>=5 we compute zeta(s) by summing \sum_{k=1}^{\infty} k^(-s) 
        f(i1)=zeta_summation(s(i1));
    end
    if (~isempty(i2)) % if Re(s)<5 and |Im(s)|<=200, we use Euler-Maclaurin method to compute zeta(s)
        f(i2)=zeta_Euler_Maclaurin(s(i2));
    end
    if (~isempty(i3)) % if Re(s)<5 and Im(s)>200, we use zeta_8(s) approximation
        f(i3)=zeta_8(s(i3));
    end
    if (~isempty(i4)) % if Re(s)<5 and Im(s)<-200, we use zeta_8(s) approximation
        f(i4)=conj(zeta_8(conj(s(i4))));
    end 
end   
%##########################################################################
function f=zeta_8(s)
% zeta_8 computes the approximation to the Riemann zeta function described in the paper
% A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function", 2025, https://arxiv.org/abs/2503.09519
%--------------------------------------------------------------------------
% the cofficients lambda_j for j=1,2,...,8
    lambda=[0.152845417613666702426-0.119440685603870510384i
            0.302346225128945757427-0.243989695504400621268i          
	        0.451119584531782942888-0.378479770209444563858i 
	        0.604563710297226464637-0.523486888629095259770i
	        0.765965706759629396959-0.678405572413543444272i
	        0.938371150977889047740-0.845332361280975174880i
	        1.128148837845288402558-1.030737947568157685685i
	        1.353030558654668162533-1.252503278108132307164i];               
    % the cofficients omega_j for j=0,1,2,...,8
    omega0=1.926019633029103199063e-1+2.472986965795651842299e-2i;
    omega=[1.582954327321094104502e-1+4.149113569204600502105e-2i   
	       7.826728293587305110862e-2+5.215518667623989653254e-2i   
	       1.940595049247490540621e-2+2.977286598777633378610e-2i  
	       1.691184771902755036966e-3+8.938933548999206800196e-3i     
	      -2.994777986686168319731e-4+1.567541981830224487301e-3i   
	      -9.837202592542590210980e-5+1.502108057352792742070e-4i  
          -9.346989286415688998740e-6+5.793852209955845432028e-6i  
	      -2.451577304299235983015e-7+6.134784898751456953524e-9i]; 
%--------------------------------------------------------------------------
    if (length(s)==1) 
        % compute chi(s)=(2*pi)^s/(2*cos(pi*s/2)*gamma(s))
        chi=exp(s*log(2*pi)+0.5i*pi*s-log(1+exp(pi*1i*s))-ln_gamma(s));
        N=floor(sqrt(imag(s)/(2*pi)));
        M=N+0.5;
        lnn=log(2:N);
        f=1+sum(exp(-s*lnn))+chi*(1+sum(exp((s-1)*lnn))); %compute the main sum 
        I1=exp(-s*log(M))*(omega0+sum(omega.*(exp(-2*pi*M*lambda-s*log(1+1i*lambda/M))+exp(2*pi*M*lambda-s*log(1-1i*lambda/M)))));
        s1=1-conj(s);
        I2=conj(exp(-s1*log(M))*(omega0+sum(omega.*(exp(-2*pi*M*lambda-s1*log(1+1i*lambda/M))+exp(2*pi*M*lambda-s1*log(1-1i*lambda/M))))));
        f=f-0.5*(-1)^N*(I1+chi*I2); 
    else
        chi=exp(s*log(2*pi)+0.5i*pi*s-log(1+exp(pi*1i*s))-ln_gamma(s));
        f=1+chi;
        N=floor(sqrt(imag(s)/(2*pi)));
        for n=2:max(N)  %compute the main sum 
            u=n.^(-s);
            f=f+(n<=N).*(u+chi./(n*u)); 
        end
        M=N+0.5;
        s1=1-conj(s);
        I1=omega0*ones(size(s));
        I2=I1;
        for j=1:8
            I1=I1+omega(j)*(exp(-2*pi*M*lambda(j)-s.*log(1+1i*lambda(j)./M))+exp(2*pi*M*lambda(j)-s.*log(1-1i*lambda(j)./M)));
            I2=I2+omega(j)*(exp(-2*pi*M*lambda(j)-s1.*log(1+1i*lambda(j)./M))+exp(2*pi*M*lambda(j)-s1.*log(1-1i*lambda(j)./M)));
        end
        I1=exp(-s.*log(M)).*I1;
        I2=conj(exp(-s1.*log(M)).*I2);
        f=f-0.5*(-1).^N.*(I1+chi.*I2); 
    end
end
%##########################################################################
function f=zeta_summation(s)
% computes zeta(s)=\sum_{n=1}^{\infty} n^{-s}  
% We remove from this sum all numbers divisible by 2,3,5,7,11,13,17,19 and truncate the resulting sum at n=500
% This results in a more efficient algorithm (fewer terms in the main sum)
% This function gives a good approximation to zeta(s) in the half-plane Re(s)>5
    n=[23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 ...
       179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263 269 271 277 281 283 293 307 311 313 317 331 337 ...
       347 349 353 359 367 373 379 383 389 397 401 409 419 421 431 433 439 443 449 457 461 463 467 479 487 491 499];
    f=ones(size(s));
    for k=1:87
        f=f+n(k).^(-s);
    end
    f=f./((1-2.^(-s)).*(1-3.^(-s)).*(1-5.^(-s)).*(1-7.^(-s)).*(1-11.^(-s)).*(1-13.^(-s)).*(1-17.^(-s)).*(1-19.^(-s)));
end
%##########################################################################
function f=zeta_Euler_Maclaurin(s)
% computes zeta(s) via Euler-Maclaurin summation 
% we use formula (25.2.10) at https://dlmf.nist.gov/25.2 with N=100 and n=15
%--------------------------------------------------------------------------
% precompute b(k)=B_{2k}/(2k)!, where B_{k} are Bernoulli numbers
    b=[8.33333333333333e-2,-1.38888888888889e-3,3.306878306878306e-05,-8.267195767195768e-07,2.087675698786810e-08, ...
      -5.284190138687493e-10,1.338253653068468e-11,-3.389680296322583e-13,8.586062056277845e-15,-2.174868698558061e-16, ...
       5.509002828360229e-18, -1.395446468581252e-19, 3.534707039629467e-21,-8.953517427037546e-23,2.267952452337683e-24];
    N=100;
    m=s/N;
    f=zeros(size(s));
    for k=1:15
        f=f+b(k)*m;
        m=m.*(s+2*k-1).*(s+2*k)/N^2;
    end
    f=1+N.^(-s).*(f+0.5+N./(s-1));
    for k=2:(N-1)
        f=f+k.^(-s);
    end
end




