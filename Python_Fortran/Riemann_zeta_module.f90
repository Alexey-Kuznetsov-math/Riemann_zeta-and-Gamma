module Riemann_zeta_module
! Riemann_zeta_module provides functions for computing the Riemann zeta function and its derivative
! -------------------------------------------------------------------------
! Author: Alexey Kuznetsov
! York University, Toronto, Canada
! Website: https://kuznetsovmath.ca/
! Email: akuznets@yorku.ca
! 
! Created: 28-Nov-2025
! Last updated: 14-April-2026
! 
! License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
!--------------------------------------------------------------------------
    implicit none
    private
    public :: Riemann_zeta, Riemann_zeta_prime

    integer, parameter            :: qp=selected_real_kind(33, 4931)
    real(kind=qp), parameter      :: pi_qp=3.141592653589793238462643383279502884_qp
    complex(kind=qp), parameter   :: i_qp=(0.0_qp,1.0_qp)
    real(kind=qp), parameter      :: log_2pi=1.8378770664093454835606594728112352797_qp   

    contains
!#######################################################################
    function Riemann_zeta(s) result(f)
!    f=Riemann_zeta(s) computes zeta(s) for complex input s
!    The input s must be a scalar (this code is not vectorized)
!-------------------------------------------------------------------------- 
    implicit none
    complex (kind=qp), intent(in) :: s    
    complex (kind=qp)             :: f    
    real (kind=qp)                :: d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (s%re>=0.5_qp) then
        f=Riemann_zeta_half_plane(s)
    elseif (abs(s)<1.0e-3_qp) then ! use Taylor series if |s| is small
        f=1/(s-1)+0.5_qp+s*(8.106146679532725821967026359438e-2_qp+s*(-3.178227954292425605050013e-3_qp&
        +s*(-7.851944770424079601768e-4_qp+s*(1.207004994288350422e-4_qp+s*(-1.940896320456e-6_qp&
        +s*(-1.301146014e-6_qp+s*(1.686158264e-7_qp-s*5.764676e-9_qp))))))) 
    else ! use reflection formula for the Riemann zeta function for Re(s)<0.5
        d=sign(1.0_qp,s%im)
        f=d*i_qp*exp(log_2pi*(s-1)-d*0.5_qp*i_qp*pi_qp*s+ln_gamma(1-s))*Riemann_zeta_half_plane(1-s)
        if (abs(s%im)<26.0_qp) f=f*(1-exp(d*i_qp*pi_qp*s))
    end if  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end function Riemann_zeta
!#######################################################################
    function Riemann_zeta_half_plane(s) result(f)
! computes zeta(s) in the half-plane Re(s)>=0.5
    implicit none
    complex (kind=qp), intent(in)  :: s    
    complex (kind=qp)              :: f
    real(kind=qp)                  :: N_sum, N_EM, N_12
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    N_sum=exp(40/max(s%re-1,0.5_qp))/4   ! computational complexity of zeta_summation function -- counting the number of evaluations of k^{-s}
    N_EM=abs(s)/8+5                      ! computational complexity of zeta_Euler_Maclaurin function         
    N_12=sqrt(abs(s%im)/(2*pi_qp))/4+50  ! computational complexity of zeta_12 function
    if (abs(s%im)<400) then 
        if (N_sum<N_EM) then
            f=zeta_summation(s)
        else
            f=zeta_Euler_Maclaurin(s)
        end if
    elseif (N_sum<N_12) then 
        f=zeta_summation(s) 
    else
        if (s%im>0) then
            f=zeta_12(s)
        else
            f=conjg(zeta_12(conjg(s)))
        end if              
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end function Riemann_zeta_half_plane
!#######################################################################    
    function zeta_summation(s) result(f)
! approximates the Riemann zeta function f=zeta(s)=\sum_{n=1}^{\infty} n^{-s}   
! to have more efficient computation, we instead compute g(s)=((1-2^{-s})(1-3^{-s})(1-5^{-s})) zeta(s) 
! which is the Dirichlet series g(z)=\sum_{n>=1, gcd(n,30)=1} n^{-s} 
! this trick removes 11/15 (or 73.3%) of the terms from the original sum \sum_{n=1}^{\infty} n^{-s}
! we truncate the Dirichlet series for g(z) at n=N*30
! for Re(s)>2 the error is smaller than (30*N)^(1-Re(s))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)  :: s
    complex (kind=qp)              :: f
    integer, parameter      :: res30(1:8)=(/1, 7, 11, 13, 17, 19, 23, 29/) ! residue classes mod 30 that are co-prime to 30
    integer                 :: N, j, l
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    N=floor(exp(40/(s%re-1))/30)  ! choose N so that the error is smaller than 10^{-17}
    f=0
    do j=1,8
        do l=0,N
            f=f+(res30(j)+30*l)**(-s)
        end do
    end do
    f=f/((1-2**(-s))*(1-3**(-s))*(1-5**(-s)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function zeta_summation
!#######################################################################    
    function RS_main_sum(s,x) result(f)
! This function computes f(1)=\sum_{1\le n \le x} n^{-s} and f(2)=\sum_{1\le n \le x} n^{s-1}.
! We reduce the number of evaluations of the exponential function n^{-s} 
! by using factorization n=m 2^a 3^b 5^c, where m is coprime to 30.
! The algorithm evaluates n^{-s} only for n=2,3,5 and for integers m coprime to 30,
! and generates all multiples of m by repeated multiplication by 2, 3, and 5.
! The number of exponential function evaluations is decreased by 73.3%,
! compared with the simple summation approach to compute f(1) and f(2). 
! This optimization allows to compute zeta(s) faster when Im(s) is very large.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    integer, intent(in)           :: x
    complex(kind=qp), intent(in)  :: s
    integer, parameter            :: res30(1:8) = (/1, 7, 11, 13, 17, 19, 23, 29/)
    integer                       :: b, r, m, n2, n3, n5, x2, x3, x5
    complex(kind=qp)              :: u2, u3, u5, v2, v3, v5, t0, t2, t3, t5
    complex(kind=qp)              :: g0, g2, g3, g5, sum1, sum2, f(1:2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    u2=2**(-s); u3=3**(-s); u5=5**(-s)
    v2=1/(2*u2); v3=1/(3*u3); v5=1/(5*u5)
    x2=x/2; x3=x/3; x5=x/5
    sum1=0.0_qp
    sum2=0.0_qp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do b=0,((x-1)/30)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do r=1,8
        m=30*b+res30(r)
        if (m>x) exit
        if (m==1) then
            t0=1.0_qp
            g0=1.0_qp
        else
            t0=m**(-s)
            g0=1/(m*t0)
        end if
        n2=m; t2=t0; g2=g0
        do
            n3=n2; t3=t2; g3=g2
            do
                n5=n3; t5=t3; g5=g3
                    do
                    sum1=sum1+t5
                    sum2=sum2+g5
                    if (n5>x5) exit
                    n5=5*n5; t5=t5*u5; g5=g5*v5
                end do
                if (n3>x3) exit
                n3=3*n3; t3=t3*u3; g3=g3*v3 
            end do
            if (n2>x2) exit
            n2=2*n2; t2=t2*u2; g2=g2*v2
        end do
    end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f(1)=sum1
    f(2)=sum2
    end function RS_main_sum            
!#######################################################################
    function zeta_Euler_Maclaurin(s) result(f)
    implicit none
    complex (kind=qp), intent(in)   :: s
    complex (kind=qp)               :: f, m, g(1:2)
    integer                         :: k, N
    real (kind=qp)                  :: N2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! computes zeta(s) via Euler-Maclaurin summation 
! we use formula (25.2.10) at https://dlmf.nist.gov/25.2 with n=25
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! precompute b(k)=B_{2k}/(2k)!, where B_{k} are Bernoulli numbers
   real(kind=qp), parameter :: b(1:25)=&    
   (/8.33333333333333333333333333333333293e-2_qp,&
    -1.38888888888888888888888888888888893e-3_qp,&
    3.30687830687830687830687830687830844e-5_qp,&
    -8.26719576719576719576719576719577148e-7_qp,&
    2.08767569878680989792100903212014440e-8_qp,&
    -5.28419013868749318484768220217955963e-10_qp,&
    1.33825365306846788328269809751291339e-11_qp,&
    -3.38968029632258286683019539124944525e-13_qp,&
    8.58606205627784456413590545042563573e-15_qp,&
    -2.17486869855806187304151642386592068e-16_qp,&
    5.50900282836022951520265260890226040e-18_qp,&
    -1.39544646858125233407076862640635689e-19_qp,&
    3.53470703962946747169322997780380424e-21_qp,&
    -8.95351742703754685040261131811275565e-23_qp,&
    2.26795245233768306031095073886817000e-24_qp,&
    -5.74479066887220244526388198760702843e-26_qp,&
    1.45517247561486490186626486727133221e-27_qp,&
    -3.68599494066531017818178247990866798e-29_qp,&
    9.33673425709504467203255515278564467e-31_qp,&
    -2.36502241570062993455963519636984353e-32_qp,&
    5.99067176248213430465991239681967181e-34_qp,&
    -1.51745488446829026171081313586472285e-35_qp,&
    3.84375812545418823222944529099024291e-37_qp,&
    -9.73635307264669103526762127925048019e-39_qp,&
    2.46624704420068095710640028028884984e-40_qp/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    N=floor(abs(s)/2)+20
    N2=1.0_qp/N**2  
    m=s/N
    f=0
    do k=1,25
        f=f+b(k)*m
        m=m*(s+2*k-1)*(s+2*k)*N2
    end do
    f=N**(-s)*(f+0.5_qp+N/(s-1))
    g=RS_main_sum(s,N-1)
    f=f+g(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    end function zeta_Euler_Maclaurin
!#######################################################################
    function zeta_12(s) result(f)
! computes the approximation zeta_{12}(s) to the Riemann zeta function described in 
! [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function", 2025,  https://arxiv.org/abs/2503.09519  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)   :: s
    complex (kind=qp)               :: f, s1, w, chi, I1, I2, g(1:2), l1(1:12), l2(1:12)
    real (kind=qp)                  :: M
    integer                         :: N    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    complex(kind=qp), parameter	:: lambda(1:12)=&	
    (/(0.12666507030478299171328588879538_qp,-0.10087636358615151558524566443202_qp),&
    (0.25122099390287200219569569142124_qp,-0.20405074409349581911004051031180_qp),&
    (0.37405911972784972054636868358418_qp,-0.31221281154087242842137513443055_qp),&
    (0.49776929817534881045635955701099_qp,-0.42635960060301833679056058340141_qp),&
    (0.62435778140713997700002673039044_qp,-0.54572495512799403188584815398318_qp),&
    (0.75466390867013388720749045161642_qp,-0.66972935838659877528943507209541_qp),&
    (0.88933701315173701350361581796676_qp,-0.79869551590918500550023908781514_qp),&         
    (1.02947832738015104394795537702100_qp,-0.93374897103422687909199888577669_qp),&              
    (1.17700666949730216301023441132388_qp,-1.07684193858043548080679607396027_qp),&
    (1.33525480422098724064095160659601_qp,-1.23127339545247342505983081362199_qp),&  
    (1.51071045858258192077553791057342_qp,-1.40345091272949325437536665778861_qp),&
    (1.72033314454678436509164855529537_qp,-1.61022353225906768661496722514114_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    complex(kind=qp), parameter	:: omega0=(1.6098168692729432925346560422022541e-1_qp,1.87823072277176154607044212350895e-2_qp)
    complex(kind=qp), parameter	:: omega(1:12)=&
    (/(1.4106875185404046965105645712788459e-1_qp,3.02465842133145100446327863373405e-2_qp),&   
    (8.8469767244257187866040969336291841e-2_qp,4.46960624042180008692385694241609e-2_qp),&   
    (3.5150310372194186343936296492476377e-2_qp,3.66197739954210998843410055490889e-2_qp),& 
    (8.1172398241584947575282044650705794e-3_qp,1.78939741414190712284813334727467e-2_qp),& 
    (6.7293124348065999744683947100109004e-4_qp,5.92406929893300257148499876708815e-3_qp),&  
    (-2.6231285561276628016088792543058718e-4_qp,1.40952189832532584552807898212223e-3_qp),&     
    (-1.2407125021539101474191351741328731e-4_qp,2.37524498834962686368608234375620e-4_qp),& 
    (-2.6548327944815849601557456533448347e-5_qp,2.63187311599115470745428256223651e-5_qp),& 
    (-3.2872693233000221851679206259221078e-6_qp,1.60583114552664837555774293393550e-6_qp),&                
    (-2.264607013792395805037591106325525e-7_qp,2.30109385236452501550096833543338e-8_qp),&    
    (-7.184384829748545261872001768132292e-9_qp,-2.17226972533644575270298136976476e-9_qp),&                
    (-6.032583864104558568614884445020897e-11_qp,-5.85674035801747810623084784839811e-11_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    ! compute chi(s)=(2*pi)^s/(2*cos(pi*s/2)*gamma(s))
    chi=exp(s*log_2pi+0.5_qp*pi_qp*i_qp*s-ln_gamma(s))
    if (s%im<26.0_qp) chi=chi/(1+exp(pi_qp*i_qp*s))
    !compute the main sum
    N=floor(sqrt(s%im/(2*pi_qp)))
    g=RS_main_sum(s,N)  ! compute the main sums in the Riemann-Siegel formula
    !compute I1=I_{M,12}(s) and  I2=conjg(I_{M,12}(1-conjg(s)))
    M=N+0.5_qp
    l1=log(1+i_qp*lambda/M)
    l2=log(1-i_qp*lambda/M)  
    w=M**(-s)
    I1=(omega0+sum(omega*(exp(-2*pi_qp*M*lambda-s*l1)+exp(2*pi_qp*M*lambda-s*l2))))*w
    s1=1-conjg(s)   
    I2=conjg(omega0+sum(omega*(exp(-2*pi_qp*M*lambda-s1*l1)+exp(2*pi_qp*M*lambda-s1*l2))))/(M*w)
    f=g(1)+chi*g(2)-0.5_qp*(-1)**N*(I1+chi*I2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    end function zeta_12
!#######################################################################
    function Riemann_zeta_prime(s) result(f)
!    f=Riemann_zeta_prime(s) computes f(1)=zeta(s) and f(2)=zeta'(s) for complex input s
!    The input s must be a scalar (this code is not vectorized)
!-------------------------------------------------------------------------- 
    implicit none
    complex (kind=qp), intent(in)   :: s    
    complex (kind=qp)               :: f(1:2), g(1:2), chi(1:2) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (s%re>=0.5_qp) then
        f=Riemann_zeta_prime_half_plane(s)
    elseif (abs(s)<1.0e-2_qp) then ! use Taylor series if |s| is small
        f=Taylor_series(s)          
    else ! use reflection formula for the Riemann zeta function for Re(s)<0.5
        chi=chi_prime(s)
        g=Riemann_zeta_prime_half_plane(1-s)
        f(1)=chi(1)*g(1)
        f(2)=chi(2)*g(1)-chi(1)*g(2)
    end if  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end function Riemann_zeta_prime 
!#######################################################################
    function Riemann_zeta_prime_half_plane(s) result(f)
! computes f(1)=zeta(s) and f(2)=zeta'(s) in the half-plane Re(s)>=0.5
    implicit none
    complex (kind=qp), intent(in)   :: s    
    complex (kind=qp)               :: f(1:2), g(1:4)
    real(kind=qp)                   :: N_sum, N_EM, N_12    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    N_sum=exp(40/max(s%re-1,0.5_qp))/4    ! computational complexity of RS_main_sums_with_derivatives function -- counting the number of evaluations of k^{-s}
    N_EM=abs(s)/8+5                       ! computational complexity of zeta_prime_Euler_Maclaurin function         
    N_12=sqrt(abs(s%im)/(2*pi_qp))/4+50   ! computational complexity of zeta_prime_12 function
    if (abs(s%im)<400) then 
        if (N_sum<N_EM) then
            g=RS_main_sums_with_derivatives(s,floor(4*N_sum))
            f=g(1:2)
        else
            f=zeta_prime_Euler_Maclaurin(s)
        end if
    elseif (N_sum<N_12) then 
        g=RS_main_sums_with_derivatives(s,floor(4*N_sum))
        f=g(1:2)    
    else
        if (s%im>0) then
            f=zeta_prime_12(s)
        else
            f=conjg(zeta_prime_12(conjg(s)))
        end if              
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end function Riemann_zeta_prime_half_plane
!#######################################################################
    function Taylor_series(s) result(f)
! computes f(1)=zeta(s) and f(2)=zeta'(s) via Taylor series expansion of g(s)=zeta(s)-1/(s-1) centered at s=0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)   :: s
    complex (kind=qp)               :: f(1:2)
    integer                         :: k
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! precompute Taylor coefficients of g(s)=zeta(s)-1/(s-1)
   real(kind=qp), parameter :: a(0:12)=&    
   (/ 0.5_qp, &
   8.10614667953272582196702635943823843e-2_qp, &
  -3.17822795429242560505001336498022119e-3_qp, &
  -7.85194477042407960176802227729143456e-4_qp, &
   1.20700499428835042199186344124721842e-4_qp, &
  -1.94089632045603779988198163183913663e-6_qp, &
  -1.30114601395962431150487297976478201e-6_qp, &
   1.68615826389220069782941984547082385e-7_qp, &
  -5.76467597994939441606374162706121810e-9_qp, &
  -9.11016489231416570921867354622218771e-10_qp, &
   1.49700759419011373520720083150847656e-10_qp, &
  -9.40689566566617690964790248200256190e-12_qp, &
  -4.09258263041583154762866605291246775e-14_qp /)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
! use the extended Horner's method
         f(1)=a(12)
         f(2)=0.0_qp
         do k=11,0,(-1) 
            f(2)=f(1)+s*f(2)
            f(1)=a(k)+s*f(1)
         end do
         f(1)=f(1)+1/(s-1)
         f(2)=f(2)-1/(s-1)**2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
        end function Taylor_series  
!#######################################################################    
    function RS_main_sums_with_derivatives(s,x) result(f)
! This function computes f(1)=\sum_{1\le n \le x} n^{-s},  f(2)=-\sum_{1\le n \le x}  \ln(n) n^{-s}, 
! f(3)=\sum_{1\le n \le x} n^{s-1}, and f(4)=\sum_{1\le n \le x}  \ln(n) n^{s-1}.
! We reduce the number of evaluations of the exponential function n^{-s} 
! by using factorization n=m 2^a 3^b 5^c, where m is coprime to 30.
! The algorithm evaluates n^{-s} only for n=2,3,5 and for integers m coprime to 30,
! and generates all multiples of m by repeated multiplication by 2, 3, and 5.
! The number of exponential function evaluations is decreased by 73.3%,
! compared with the simple summation approach to compute f(1) and f(2). 
! This optimization allows to compute zeta(s) faster when Im(s) is large.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    integer, intent(in)             :: x
    complex(kind=qp), intent(in)    :: s
    integer, parameter              :: res30(1:8) = (/1, 7, 11, 13, 17, 19, 23, 29/)
    integer                         :: b, r, m, n2, n3, n5, x2, x3, x5
    complex(kind=qp)                :: u2, u3, u5, v2, v3, v5, t0, t2, t3, t5
    complex(kind=qp)                :: g0, g2, g3, g5, sum1, sum2, sum3, sum4, f(1:4)
    real(kind=qp)                   :: l0, l2, l3, l5, ln2, ln3, ln5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    u2=2**(-s); u3=3**(-s); u5=5**(-s)
    v2=1/(2*u2); v3=1/(3*u3); v5=1/(5*u5)
    ln2=log(2.0_qp); ln3=log(3.0_qp); ln5=log(5.0_qp)
    x2=x/2; x3=x/3; x5=x/5
    sum1=0.0_qp; sum2=0.0_qp; sum3=0.0_qp; sum4=0.0_qp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do b=0,((x-1)/30)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do r=1,8
        m=30*b+res30(r)
        if (m>x) exit
        if (m==1) then
            t0=1.0_qp
            g0=1.0_qp
            l0=0.0_qp
        else
            l0=log(m+0.0_qp)
            t0=exp(-s*l0)
            g0=1.0_qp/(m*t0)
        end if
        n2=m; t2=t0; g2=g0; l2=l0
        do
            n3=n2; t3=t2; g3=g2; l3=l2
            do
                n5=n3; t5=t3; g5=g3; l5=l3
                    do
                    sum1=sum1+t5
                    sum2=sum2+l5*t5
                    sum3=sum3+g5
                    sum4=sum4+l5*g5
                    if (n5>x5) exit
                    n5=5*n5; t5=t5*u5; g5=g5*v5; l5=l5+ln5
                end do
                if (n3>x3) exit
                n3=3*n3; t3=t3*u3; g3=g3*v3; l3=l3+ln3 
            end do
            if (n2>x2) exit
            n2=2*n2; t2=t2*u2; g2=g2*v2; l2=l2+ln2
        end do
    end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    f(1)=sum1; f(2)=-sum2; f(3)=sum3; f(4)=sum4
    end function RS_main_sums_with_derivatives  
!#######################################################################
    function zeta_prime_Euler_Maclaurin(s) result(f)
! computes f(1)=zeta(s) and f(2)=zeta'(s) via Euler-Maclaurin summation 
! we use formula (25.2.10) at https://dlmf.nist.gov/25.2 with n=25
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)   :: s
    complex (kind=qp)               :: f(1:2), m, m1, w, g(1:4)
    integer                         :: k, N
    real (kind=qp)                  :: N2, v
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! precompute b(k)=B_{2k}/(2k)!, where B_{k} are Bernoulli numbers
   real(kind=qp), parameter :: b(1:25)=&    
   (/8.33333333333333333333333333333333293e-2_qp,&
    -1.38888888888888888888888888888888893e-3_qp,&
    3.30687830687830687830687830687830844e-5_qp,&
    -8.26719576719576719576719576719577148e-7_qp,&
    2.08767569878680989792100903212014440e-8_qp,&
    -5.28419013868749318484768220217955963e-10_qp,&
    1.33825365306846788328269809751291339e-11_qp,&
    -3.38968029632258286683019539124944525e-13_qp,&
    8.58606205627784456413590545042563573e-15_qp,&
    -2.17486869855806187304151642386592068e-16_qp,&
    5.50900282836022951520265260890226040e-18_qp,&
    -1.39544646858125233407076862640635689e-19_qp,&
    3.53470703962946747169322997780380424e-21_qp,&
    -8.95351742703754685040261131811275565e-23_qp,&
    2.26795245233768306031095073886817000e-24_qp,&
    -5.74479066887220244526388198760702843e-26_qp,&
    1.45517247561486490186626486727133221e-27_qp,&
    -3.68599494066531017818178247990866798e-29_qp,&
    9.33673425709504467203255515278564467e-31_qp,&
    -2.36502241570062993455963519636984353e-32_qp,&
    5.99067176248213430465991239681967181e-34_qp,&
    -1.51745488446829026171081313586472285e-35_qp,&
    3.84375812545418823222944529099024291e-37_qp,&
    -9.73635307264669103526762127925048019e-39_qp,&
    2.46624704420068095710640028028884984e-40_qp/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    N=floor(abs(s)/2)+20
    N2=1.0_qp/N**2  
    m=s/N
    m1=1.0_qp/N 
    f=0.0_qp
    do k=1,25
        f(1)=f(1)+b(k)*m
        f(2)=f(2)+b(k)*m1
        w=(s+2*k-1)*(s+2*k)
        m1=(m1*w+m*(2*s+4*k-1))*N2
        m=m*w*N2
    end do
    v=log(N+0.0_qp)
    w=exp(-s*v)
    f(1)=w*(f(1)+0.5_qp+N/(s-1))
    f(2)=-v*f(1)+w*(f(2)-N/(s-1)**2)
    g=RS_main_sums_with_derivatives(s,N-1)
    f=f+g(1:2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    end function zeta_prime_Euler_Maclaurin 
!#######################################################################    
    function zeta_prime_12(s) result(f)
! computes the approximations f(1)=zeta_{12}(s) and f(2)=zeta_{12}'(s) to the Riemann zeta function and its derivative described in 
! [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function", 2025,  https://arxiv.org/abs/2503.09519  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)   :: s
    complex (kind=qp)               :: f(1:2), s1, chi(1:2), I1, I2, I1_prime, I2_prime
    complex (kind=qp)               :: l1(1:12), l2(1:12), e1(1:12), e2(1:12), w, g(1:4)
    real (kind=qp)                  :: M, v
    integer                         :: N    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    complex(kind=qp), parameter	:: lambda(1:12)=&	
    (/(0.12666507030478299171328588879538_qp,-0.10087636358615151558524566443202_qp),&
    (0.25122099390287200219569569142124_qp,-0.20405074409349581911004051031180_qp),&
    (0.37405911972784972054636868358418_qp,-0.31221281154087242842137513443055_qp),&
    (0.49776929817534881045635955701099_qp,-0.42635960060301833679056058340141_qp),&
    (0.62435778140713997700002673039044_qp,-0.54572495512799403188584815398318_qp),&
    (0.75466390867013388720749045161642_qp,-0.66972935838659877528943507209541_qp),&
    (0.88933701315173701350361581796676_qp,-0.79869551590918500550023908781514_qp),&         
    (1.02947832738015104394795537702100_qp,-0.93374897103422687909199888577669_qp),&              
    (1.17700666949730216301023441132388_qp,-1.07684193858043548080679607396027_qp),&
    (1.33525480422098724064095160659601_qp,-1.23127339545247342505983081362199_qp),&  
    (1.51071045858258192077553791057342_qp,-1.40345091272949325437536665778861_qp),&
    (1.72033314454678436509164855529537_qp,-1.61022353225906768661496722514114_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    complex(kind=qp), parameter	:: omega0=(1.6098168692729432925346560422022541e-1_qp,1.87823072277176154607044212350895e-2_qp)
    complex(kind=qp), parameter	:: omega(1:12)=&
    (/(1.4106875185404046965105645712788459e-1_qp,3.02465842133145100446327863373405e-2_qp),&   
    (8.8469767244257187866040969336291841e-2_qp,4.46960624042180008692385694241609e-2_qp),&   
    (3.5150310372194186343936296492476377e-2_qp,3.66197739954210998843410055490889e-2_qp),& 
    (8.1172398241584947575282044650705794e-3_qp,1.78939741414190712284813334727467e-2_qp),& 
    (6.7293124348065999744683947100109004e-4_qp,5.92406929893300257148499876708815e-3_qp),&  
    (-2.6231285561276628016088792543058718e-4_qp,1.40952189832532584552807898212223e-3_qp),&     
    (-1.2407125021539101474191351741328731e-4_qp,2.37524498834962686368608234375620e-4_qp),& 
    (-2.6548327944815849601557456533448347e-5_qp,2.63187311599115470745428256223651e-5_qp),& 
    (-3.2872693233000221851679206259221078e-6_qp,1.60583114552664837555774293393550e-6_qp),&                
    (-2.264607013792395805037591106325525e-7_qp,2.30109385236452501550096833543338e-8_qp),&    
    (-7.184384829748545261872001768132292e-9_qp,-2.17226972533644575270298136976476e-9_qp),&                
    (-6.032583864104558568614884445020897e-11_qp,-5.85674035801747810623084784839811e-11_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    chi=chi_prime(s)
    N=floor(sqrt(s%im/(2*pi_qp)))
    M=N+0.5_qp
    g=RS_main_sums_with_derivatives(s,N)
    l1=log(1+i_qp*lambda/M)
    l2=log(1-i_qp*lambda/M)
    v=log(M)
    w=exp(-v*s)
    e1=exp(-2*pi_qp*M*lambda-l1*s)
    e2=exp(2*pi_qp*M*lambda-l2*s)   
    I1=w*(omega0+sum(omega*(e1+e2)))
    I1_prime=-v*I1-w*sum(omega*(l1*e1+l2*e2))
    s1=1-conjg(s)   
    w=exp(-v*s1)
    e1=exp(-2*pi_qp*M*lambda-l1*s1)
    e2=exp(2*pi_qp*M*lambda-l2*s1)
    I2=conjg(w*(omega0+sum(omega*(e1+e2))))
    I2_prime=v*I2+conjg(w*sum(omega*(l1*e1+l2*e2)))
    f(1)=g(1)+chi(1)*g(3)-0.5_qp*(-1)**N*(I1+chi(1)*I2)
    f(2)=g(2)+chi(2)*g(3)+chi(1)*g(4)-0.5_qp*(-1)**N*(I1_prime+chi(2)*I2+chi(1)*I2_prime)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    end function zeta_prime_12  
!#######################################################################    
    function chi_prime(s) result(f)
! computes f(1)=chi(s) and f(2)=chi'(s), where chi(s)=2^s pi^{s-1} sin(pi s/2) Gamma(1-s)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex (kind=qp), intent(in)   :: s    
    complex (kind=qp)               :: w, f(1:2)        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    if (aimag(s)>=0) then
        w=0.0_qp
        if (s%im<26) w=exp(pi_qp*i_qp*s)
        f(1)=i_qp*(1-w)*exp(log_2pi*(s-1)-0.5_qp*pi_qp*i_qp*s+ln_gamma(1-s))
        f(2)=f(1)*(log_2pi-0.5_qp*pi_qp*i_qp*(1+w)/(1-w)-psi(1-s))
    else
        w=0.0_qp
        if (s%im>-26) w=exp(-pi_qp*i_qp*s)
        f(1)=-i_qp*(1-w)*exp(log_2pi*(s-1)+0.5_qp*pi_qp*i_qp*s+ln_gamma(1-s))
        f(2)=f(1)*(log_2pi+0.5_qp*pi_qp*i_qp*(1+w)/(1-w)-psi(1-s))  
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    end function chi_prime      
!#######################################################################
    function ln_gamma(z) result(f)
! f=ln_gamma(z)  Computes the logarithm of the Gamma function in the entire complex plane to relative accuracy 10^{-33} 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    implicit none
    complex(kind=qp), intent(in)    :: z
    complex(kind=qp)                :: f
    real(kind=qp)                   :: d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (z%re>=0.5_qp) then
        f=ln_gamma_half_plane(z)
    else   ! use reflection formula
        d=sign(1.0_qp,z%im)
        f=log_2pi-ln_gamma_half_plane(1-z)+d*pi_qp*i_qp*(z-0.5_qp)
        if (abs(z%im)<13.0_qp) f=f-log(1-exp(2*d*pi_qp*i_qp*z))
    end if      
    end function ln_gamma   
!#######################################################################    
    function ln_gamma_half_plane(z) result(f)
! Computes the logarithm of the Gamma function in the half-plane Re(z)>=0.5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex(kind=qp), intent(in)    :: z
    complex(kind=qp)                :: f, w
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  complex(kind=qp), parameter   :: c(1:10)=&  
  (/(5.24616021279591452486366931886507231e-18_qp,-7.28858517089204357378285631098968017e-19_qp),&
    (-1.03645727946895950202973370696517980e-14_qp,-4.91235660742904730225402198655414545e-16_qp),&
    (2.83116216383233958589589637252505638e-12_qp,1.77570368606396873629631397951997996e-12_qp),&
    (-1.40944973572617784675246296142057149e-10_qp,-3.53502539895230844271539953919813180e-10_qp),&
    (-5.20787854697287794169176796369336582e-9_qp,1.94617502015990090067084029863278257e-8_qp),&
    (3.47966411442394778440342480541775312e-7_qp,-3.33535383403921029829086583875030695e-7_qp),&
    (-3.18843642498486831790083135535258358e-6_qp,-1.82958014796132699877266562780625476e-6_qp),&
    (-5.57562016424223636703328227099497974e-5_qp,4.92679042937267561935689700370687005e-5_qp),&
    (1.38643273928272474218148558674608392e-3_qp,4.83806934765748596342211686330260287e-4_qp),&
    (-1.59803102233662871107318165146306651e-2_qp,-1.52516955883481628555080732879825731e-2_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   real(kind=qp), parameter     :: c_r(1:10)=&  
   (/8.38343225530011439082309752725882187e-2_qp,&
   -6.74430082338650078957962313610703220e-2_qp,&
   -2.70213900992268839762892525434748157e-2_qp,&
    1.57803808804528694981326494005746286e-2_qp,&
    2.33954804824275664453236597530851904e-2_qp,&
    2.52284396609808043688885398749700275e-2_qp,&
    2.30809809521001645590635555499840879e-2_qp,&
    1.83501055870866307420764135274564410e-2_qp,&
    1.21423902396302907922628772466779253e-2_qp,&
    5.29059031423024471955612225192891902e-3_qp/)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   complex(kind=qp), parameter  :: lambda(1:10)=&  
   (/(6.02971195313860842539370931670450114_qp,6.11437996129723446402665525970436514_qp),& 
     (5.56804054759171092241042892117635281_qp,4.85457153733643031933210712482122880_qp),& 
     (5.15844569491674055604349517690042181_qp,3.90327112121204740598345590376013661_qp),& 
     (4.76871586259457801768396785533919680_qp,3.12556145061703120318850853810674408_qp),& 
     (4.38341245174507563739820419799248207_qp,2.47337891615845342867937611036767555_qp),& 
     (3.99133784291107400227962001230393159_qp,1.94137541852531086840355444808637554_qp),& 
     (3.66958804329562827823962808418181359_qp,1.56236796729048577916994524331369589_qp),& 
     (3.46159050180126136182926728376695136_qp,1.12410611018322832754089391206660747_qp),& 
     (3.15765276731885888010820360687494277_qp,6.55786691120793104689870709013343090e-1_qp),& 
     (2.76850563819759711075932558400262024_qp,2.34600416035934699487584444243166950e-1_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
   real(kind=qp), parameter :: lambda_r(1:10)=&  
   (/2.18124521405284135318527507076287705_qp,& 
    1.84252251506495844651969378580123083_qp,& 
    1.65768924207144731540686634527932492_qp,& 
    1.36817642421453973550980633810027734_qp,& 
    1.25980243023253236269819828718601771_qp,& 
    1.17271186374413119803174350058265146_qp,& 
    1.10507741943276034796403254535216293_qp,& 
    1.05533682765507104656603577401315009_qp,& 
    1.02210293012682038807462088453101857_qp,& 
    1.00414527528635943285444465013050254_qp/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (abs(z)>200.0_qp) then  ! use Stirling's asymptotic formula
        w=z**(-2)
        f=(z-0.5_qp)*log(z)-z+log_2pi/2+(1.0_qp/12+w*(-1.0_qp/360+w*(1.0_qp/1260+&
            w*(-1.0_qp/1680+w*(1.0_qp/1188+w*(-691.0_qp/360360+w*(1.0_qp/156)))))))/z           
    elseif (z%re>=1.5_qp) then ! use Binet's formula
        w=z-1
        f=sum(c/(w+lambda)+conjg(c)/(w+conjg(lambda))+c_r/(w+lambda_r))+(z-0.5_qp)*log(z)-z+log_2pi/2       
    else ! in the strip 0.5<=Re(z)<1.5 we use the functional equation log(Gamma(z))=log(Gamma(z+1))-log(z)  
        f=sum(c/(z+lambda)+conjg(c)/(z+conjg(lambda))+c_r/(z+lambda_r))+(z+0.5_qp)*log(z+1)-z-1+log_2pi/2-log(z)                    
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function ln_gamma_half_plane
!#######################################################################
    function psi(z) result(f)
! f=psi(z)  Computes the digamma function psi(z)=Gamma'(z)/Gamma(z) in the entire complex plane 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    implicit none
    complex(kind=qp), intent(in)    :: z
    complex(kind=qp)                :: f
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (z%re>=0.5_qp) then
        f=psi_half_plane(z)
    else   ! use reflection formula
        if (abs(z%im)<13.0_qp) then
            f=psi_half_plane(1-z)-pi_qp/tan(pi_qp*z)
        else
            f=psi_half_plane(1-z)+sign(1.0_qp,z%im)*pi_qp*i_qp
        end if
    end if      
    end function psi    
!#######################################################################    
    function psi_half_plane(z) result(f)
! Computes the digamma function psi(z)=Gamma'(z)/Gamma(z) in the half-plane Re(z)>=0.5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    implicit none
    complex(kind=qp), intent(in)    :: z
    complex(kind=qp)                :: f, w
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  complex(kind=qp), parameter   :: c(1:10)=&  
  (/(5.24616021279591452486366931886507231e-18_qp,-7.28858517089204357378285631098968017e-19_qp),&
    (-1.03645727946895950202973370696517980e-14_qp,-4.91235660742904730225402198655414545e-16_qp),&
    (2.83116216383233958589589637252505638e-12_qp,1.77570368606396873629631397951997996e-12_qp),&
    (-1.40944973572617784675246296142057149e-10_qp,-3.53502539895230844271539953919813180e-10_qp),&
    (-5.20787854697287794169176796369336582e-9_qp,1.94617502015990090067084029863278257e-8_qp),&
    (3.47966411442394778440342480541775312e-7_qp,-3.33535383403921029829086583875030695e-7_qp),&
    (-3.18843642498486831790083135535258358e-6_qp,-1.82958014796132699877266562780625476e-6_qp),&
    (-5.57562016424223636703328227099497974e-5_qp,4.92679042937267561935689700370687005e-5_qp),&
    (1.38643273928272474218148558674608392e-3_qp,4.83806934765748596342211686330260287e-4_qp),&
    (-1.59803102233662871107318165146306651e-2_qp,-1.52516955883481628555080732879825731e-2_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   real(kind=qp), parameter     :: c_r(1:10)=&  
   (/8.38343225530011439082309752725882187e-2_qp,&
   -6.74430082338650078957962313610703220e-2_qp,&
   -2.70213900992268839762892525434748157e-2_qp,&
    1.57803808804528694981326494005746286e-2_qp,&
    2.33954804824275664453236597530851904e-2_qp,&
    2.52284396609808043688885398749700275e-2_qp,&
    2.30809809521001645590635555499840879e-2_qp,&
    1.83501055870866307420764135274564410e-2_qp,&
    1.21423902396302907922628772466779253e-2_qp,&
    5.29059031423024471955612225192891902e-3_qp/)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   complex(kind=qp), parameter  :: lambda(1:10)=&  
   (/(6.02971195313860842539370931670450114_qp,6.11437996129723446402665525970436514_qp),& 
     (5.56804054759171092241042892117635281_qp,4.85457153733643031933210712482122880_qp),& 
     (5.15844569491674055604349517690042181_qp,3.90327112121204740598345590376013661_qp),& 
     (4.76871586259457801768396785533919680_qp,3.12556145061703120318850853810674408_qp),& 
     (4.38341245174507563739820419799248207_qp,2.47337891615845342867937611036767555_qp),& 
     (3.99133784291107400227962001230393159_qp,1.94137541852531086840355444808637554_qp),& 
     (3.66958804329562827823962808418181359_qp,1.56236796729048577916994524331369589_qp),& 
     (3.46159050180126136182926728376695136_qp,1.12410611018322832754089391206660747_qp),& 
     (3.15765276731885888010820360687494277_qp,6.55786691120793104689870709013343090e-1_qp),& 
     (2.76850563819759711075932558400262024_qp,2.34600416035934699487584444243166950e-1_qp)/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
   real(kind=qp), parameter :: lambda_r(1:10)=&  
   (/2.18124521405284135318527507076287705_qp,& 
    1.84252251506495844651969378580123083_qp,& 
    1.65768924207144731540686634527932492_qp,& 
    1.36817642421453973550980633810027734_qp,& 
    1.25980243023253236269819828718601771_qp,& 
    1.17271186374413119803174350058265146_qp,& 
    1.10507741943276034796403254535216293_qp,& 
    1.05533682765507104656603577401315009_qp,& 
    1.02210293012682038807462088453101857_qp,& 
    1.00414527528635943285444465013050254_qp/)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (abs(z)>200.0_qp) then  ! use asymptotic formula 5.11.2 in https://dlmf.nist.gov/5.11
        w=z**(-2)
        f=log(z)-0.5_qp/z-w*(1.0_qp/12+w*(-1.0_qp/120+w*(1.0_qp/252+&
            w*(-1.0_qp/240+w*(1.0_qp/132+w*(-691.0_qp/32760+w*(1.0_qp/12-w*(3617.0_qp/8160))))))))          
    elseif (z%re>=1.5_qp) then ! use derivative of Binet's formula (see the approximation for ln_gamma)
        w=z-1
        f=-sum(c/(w+lambda)**2+conjg(c)/(w+conjg(lambda))**2+c_r/(w+lambda_r)**2)+log(z)-0.5_qp/z       
    else ! in the strip 0.5<=Re(z)<1.5 we use the functional equation psi(z)=psi(z+1)-1/z   
        f=-sum(c/(z+lambda)**2+conjg(c)/(z+conjg(lambda))**2+c_r/(z+lambda_r)**2)+log(z+1)-0.5_qp/(z+1)-1.0_qp/z                    
    end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function psi_half_plane 
!#######################################################################
end module Riemann_zeta_module
