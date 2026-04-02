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
% Last updated: 2-April-2026
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
        f(i3)=log(1-w3)-ln_gamma_half_plane(2-w3)+1.83787706640934548+1i*pi*(w3-0.5)-log(1-exp(2i*pi*w3));
    end
    if (~isempty(i4))  % use reflection formula
        w4=z(i4);
        f(i4)=log(1-w4)-ln_gamma_half_plane(2-w4)+1.83787706640934548-1i*pi*(w4-0.5)-log(1-exp(-2i*pi*w4));
    end
end
% -------------------------------------------------------------------------
function f=ln_gamma_half_plane(z)
% approximates the logarithm of the Gamma function in the half-plane Re(z)>=1.5
% -------------------------------------------------------------------------
    c=[ 7.8489806116985789e-8+1.69933072157308014e-8i
       -2.0915898156004773e-5-1.57315657837253062e-5i
        8.0382728578667684e-4+1.46123728240593489e-3i
       -1.7781861010073780e-2-3.91857202978460444e-2i];
    l=[3.089224411647167808+2.346905763041727244i
       2.676547234010318606+1.370743226374250769i
       2.266206347366603493+6.752681709263032345e-1i
       1.829497420906493754+1.564512477809790312e-1i];
    cr=[-3.32926108478200017e-6
         4.79008236266177865e-2
         4.62832185218743242e-2
         2.31503627111999874e-2];
    lr=[3.754580109023743872
        1.253693011087479497   
        1.099995202351846400 
        1.018671901578566989];
    f=zeros(size(z));
    for j=1:4
        f=f+c(j)./(z-1+l(j))+conj(c(j))./(z-1+conj(l(j)))+cr(j)./(z-1+lr(j));
    end
    f=f+0.91893853320467274+(z-0.5).*log(z)-z;
end
% -------------------------------------------------------------------------
