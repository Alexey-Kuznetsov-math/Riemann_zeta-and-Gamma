function f=ln_gamma(z)
% ln_gamma  Computes the logarithm of the Gamma function in the entire complex plane.
%
%   f = ln_gamma(z) returns log(Gamma(z)) for any complex input z.
%   The input z can be a scalar, vector, or array.
%
% -------------------------------------------------------------------------
% Author: Alexey Kuznetsov
% York University, Toronto, Canada
% Website: https://kuznetsovmath.ca/
% Email: akuznets@yorku.ca
%
% Created: 11-Nov-2025
% Last updated: 29-Nov-2025
%
% License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
% -------------------------------------------------------------------------
    i1=find(real(z)>=1.5);
    i2=find((real(z)>=0.5)&(real(z)<1.5));
    i3=find((real(z)<0.5)&(imag(z)>=0));
    i4=find((real(z)<0.5)&(imag(z)<0));
    f=zeros(size(z));
    if (~isempty(i1))
        f(i1)=ln_gamma_half_plane(z(i1));
    end
    if (~isempty(i2)) % use functional equation log(Gamma(z))=log(Gamma(z+1))-log(z)
        f(i2)=ln_gamma_half_plane(z(i2)+1)-log(z(i2));
    end
    if (~isempty(i3))  % use reflection formula
        w3=z(i3);
        f(i3)=log(1-w3)-ln_gamma_half_plane(2-w3)+1.83787706640935+1i*pi*(w3-0.5)-log(1-exp(2i*pi*w3));
    end
    if (~isempty(i4))  % use reflection formula
        w4=z(i4);
        f(i4)=log(1-w4)-ln_gamma_half_plane(2-w4)+1.83787706640935-1i*pi*(w4-0.5)-log(1-exp(-2i*pi*w4));
    end
end
% -------------------------------------------------------------------------
function f=ln_gamma_half_plane(z)
% approximates the logarithm of the Gamma function in the half-plane Re(z)>=1.5
% -------------------------------------------------------------------------
    c=[-5.3035486658210424e-9+8.47546314707099368e-10i
        1.0093321787093660e-6+2.72761804689954947e-7i
       -2.9946104278440048e-5+3.68867702224184848e-5i
        5.1324321142829496e-4+1.31127415401525033e-3i];
    l=[3.036399609619316394+2.183538146270238962i
       2.619619673091566707+1.285513733113457465i 
       2.236468161644125811-6.451943678012780725e-1i
       1.822691064348029853+1.260996493971840621e-1i];
    cr=[ 1.91267259422448627e-12
        -1.34922683153592188e-3
        -8.28182554958193984e-4
        -1.80081775867242537e-4];
    lr=[7.862008471885223854
        1.278757134063109466   
        1.129114403699685475   
        1.037917494906115128];
    f=zeros(size(z));
    for j=1:4
        f=f+c(j)./(z-1+l(j)).^3+conj(c(j))./(z-1+conj(l(j))).^3+cr(j)./(z-1+lr(j)).^3;
    end
    f=2*f+0.918938533204673+(z-0.5).*log(z)-z+0.083333333333333./z;
end
% -------------------------------------------------------------------------
