! to compile for python on linux with ifort do
! dof2py2 -x '-openmp -D__OPENMP -fast ' --open_mp interp.F90
module interp
  IMPLICIT NONE
  REAL*8 :: grav  = 9.80665
  REAL*8 :: rn_alpha = 2.e-4
  REAL*8 ::rau0 = 1035.d0
  REAL*8 :: rn_beta = 7.7e-4
  INTEGER*4 :: neos = 0
contains
  subroutine interpolate8(kmt,T,S,d3,d0_in,km,jm,im, &
       & k_below_s,r_above_s,T_s,S_s,outcropmask,groundmask)
    INTEGER , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*8 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),d3(km,im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    INTEGER , INTENT(OUT) :: k_below_s(im,jm)
    REAL*8 , INTENT(OUT) :: r_above_s(im,jm),T_s(im,jm),S_s(im,jm)
    LOGICAL(KIND=1) , INTENT(OUT) :: outcropmask(im,jm),groundmask(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,d3,d0_in
    !f2py intent (out) k_below,r_above_s,T_s,S_s,outcropmask,groundmask

    INTEGER, allocatable :: kin(:)
    INTEGER :: i,j,k, kmc
    REAL*8 :: rdd,r_below,r_above
    REAL*8 :: d0

    outcropmask(:,:) = .false.
    groundmask(:,:) = .false.

    allocate ( kin(im) )


    kin(:) = 10000
    d0 = d0_in
    do j=1,jm
!!!     !$omp parallel do private(kmc,k,r_above,rdd,r_below) shared(im,j,kmt,kin,d3,outcropmask,groundmask)
       do i=1,im
          kmc = kmt(i,j)
          if (kmc==0) cycle
          call huntc18(kmc,d3(1:kmc,i,j),kin(i),k,d0)
          k_below_s(i,j) = k

          if (k==0) then
             outcropmask(i,j) = .true.

          else if (k>=kmc) then
             groundmask(i,j) = .true.
             k=255

          else
             rdd = 1./(d3(k+1,i,j)-d3(k,i,j))
             r_above = (d3(k+1,i,j) - d0)*rdd
             r_below = 1. - r_above
             r_above_s(i,j) = r_above

             S_s(i,j) = r_above*S(k,i,j)+r_below*S(k+1,i,j)
             T_s(i,j) = r_above*T(k,i,j)+r_below*T(k+1,i,j)
          end if

          kin(i) = k
       end do
    end do
  end subroutine interpolate8

  subroutine interpolate4(kmt,T,S,rho,drho0,d0_in,km,jm,im, &
       & k_below_s,r_above_s,T_s,S_s,outcropmask,groundmask)
    INTEGER , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*4 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),rho(km,im,jm),drho0(km,im,jm)
    INTEGER , INTENT(OUT) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    REAL*8 , INTENT(OUT) :: r_above_s(im,jm),T_s(im,jm),S_s(im,jm)
    LOGICAL(KIND=1), INTENT(OUT) :: outcropmask(im,jm),groundmask(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,rho,drho0,d0_in
    !f2py intent (out) k_below,r_above_s,T_s,S_s,outcropmask,groundmask

    INTEGER, allocatable :: kin(:)
    ! REAL*4, allocatable :: d(:)
    REAL*4 :: d(100)
    INTEGER :: i,j,k,kk,kmc
    REAL*8 :: dd,r_below,r_above
    REAL*4 :: d0

    outcropmask(:,:) = .false.
    groundmask(:,:) = .false.

    allocate ( kin(im) )
    ! allocate ( d(km) )


    kin(:) = 10000
    d0 = d0_in
    do j=1,jm
       !$omp parallel do private(i,kmc,k,d,r_above,dd,r_below) shared(im,j,kmt,kin,rho,drho0,d0,outcropmask,groundmask,T_s,S_s)
       do i=1,im
          kmc = kmt(i,j)
          if (kmc==0) cycle
          d(1:kmc) = rho(1:kmc,i,j) - drho0(1:kmc,i,j)
          call huntc14(kmc,d(1:kmc),kin(i),k,d0)
          if (k<0) then
             do kk=2,kmc
                if(d(kk)>=d0) then
                   k=kk-1
                   exit
                endif
             end do
          end if

          if (k==0) then
             outcropmask(i,j) = .true.

          else if (k>=kmc) then
             groundmask(i,j) = .true.
             k = 127

          else
             dd = (d(k+1)-d(k))
             if (abs(dd)>1.d-20) then
                r_above = (d(k+1) - d0)/dd
                r_below = 1. - r_above
             else
                r_above = 0.5d0
                r_below = 0.5d0
             endif
             r_above_s(i,j) = r_above

             S_s(i,j) = r_above*S(k,i,j)+r_below*S(k+1,i,j)
             T_s(i,j) = r_above*T(k,i,j)+r_below*T(k+1,i,j)
          end if

          kin(i) = k
          k_below_s(i,j) = k

       end do
       !$omp  end parallel do
    end do
  end subroutine interpolate4

  subroutine mginterpolate4(kmt,T,S,rho,drho0,ssh,dzw,depth, &
       & k_below_s,r_above_s,active,depth_km,d0_in,km,jm,im, &
       & z_s,Mg)
    INTEGER*4 , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*4 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),rho(km,im,jm),drho0(km,im,jm),ssh(im,jm)
    REAL*8 , INTENT(IN) :: dzw(km,im,jm)
    REAL*4 , INTENT(IN) :: depth(km,im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    INTEGER*4 , INTENT(IN) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: r_above_s(im,jm), depth_km
    LOGICAL(KIND=1) , INTENT(IN) :: active(im,jm)
    REAL*8 , INTENT(OUT) :: z_s(im,jm),Mg(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,rho,drho0,ssh,dzw,depth,k_below,r_above_s,depth_km,d0_in,active
    !f2py intent (out) z_s,Mg

    INTEGER*4 :: i,j,k, kmc
    REAL*8 :: r_below,r_above,p_above,p,z
    REAL*8 :: rho00 = 1035.d0
    REAL*8 :: rrho0 = 1.d0/1035.d0
    REAL*8 :: grav = 9.80665d0
    REAL*8 :: buoy(100)
    ! REAL*8, ALLOCATABLE :: rhbar(:)
    ! REAL*4, ALLOCATABLE :: rho(:)
    REAL*4 :: b0

    b0 = -d0_in*grav*rrho0
    do j=1,jm
       !$omp parallel do private(kmc,k,i,r_above,r_below,z,buoy,p_above,p) shared(active,im,j,kmt,k_below_s,r_above_s,depth_km,rho,drho0,b0,z_s,Mg)
       do i=1,im
          if (active(i,j)) then
             kmc = kmt(i,j)
             k = k_below_s(i,j)
             r_above = r_above_s(i,j)
             r_below = 1.d0 - r_above
             z = r_above*depth(k,i,j) + r_below*depth(k+1,i,j)
             z_s(i,j) = z
             buoy(1:k) = -grav*rrho0*(rho(1:k,i,j) - drho0(1:k,i,j))
             p_above = ssh(i,j)*grav - dzw(1,i,j)*buoy(1)
             if(k>1) then
                p_above = p_above - dot_product(.5d0*(buoy(1:k-1)+buoy(2:k)),dzw(2:k,i,j))
             end if
             p = p_above - .5d0*(buoy(k) + b0)*r_below*dzw(k+1,i,j)
             Mg(i,j) = p + b0*z
          end if
       end do
       !$omp  end parallel do
    end do
  end subroutine mginterpolate4

  subroutine siginterpolate4(T,S,k_below_s,r_above_s,active,depth_km,median_depth_km,km,jm,im, &
       & sig_s,median_sig_s)
    INTEGER*4 , INTENT(IN) :: km,jm,im
    REAL*4 , INTENT(IN) :: T(km,im,jm),S(km,im,jm)
    INTEGER*4 , INTENT(IN) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: r_above_s(im,jm), depth_km, median_depth_km
    LOGICAL(KIND=1) , INTENT(IN) :: active(im,jm)
    REAL*8 , INTENT(OUT) :: sig_s(im,jm),median_sig_s(im,jm)
    !f2py intent (in) km,jm,im,T,S,k_below,r_above_s,active,depth_km,median_depth_km
    !f2py intent (out) sig_s,median_sig_s

    INTEGER*4 :: i,j,k
    REAL*8 :: r_below,r_above !,sig8(2),median_sig8(2),T8(2),S8(2)
    REAL*4 :: sig(2),median_sig(2)

    do j=1,jm
       !$omp parallel do private(k,i,r_above,r_below,sig,median_sig) shared(active,im,j,k_below_s,r_above_s,depth_km,median_depth_km,sig_s,median_sig_s)
       do i=1,im
          if (active(i,j)) then
             k = k_below_s(i,j)
             r_above = r_above_s(i,j)
             r_below = 1.d0 - r_above
             ! T8 = T(k:k+1,i,j)
             ! S8 = S(k:k+1,i,j)
             ! call sigma_n(T8,S8,2,depth_km,sig8)
             ! call sigma_n(T8,S8,2,median_depth_km,median_sig8)
             call eos_sigman4(T(k:k+1,i,j),S(k:k+1,i,j),depth_km,sig)
             call eos_sigman4(T(k:k+1,i,j),S(k:k+1,i,j),median_depth_km,median_sig)
             sig_s(i,j) = r_above*sig(1) + r_below*sig(2)
             median_sig_s(i,j) = r_above*median_sig(1) + r_below*median_sig(2)
          end if
       end do
       !$omp  end parallel do
    end do
  end subroutine siginterpolate4

  subroutine mginterpolate8(kmt,T,S,d3,ssh,dzw,depth, &
       & k_below_s,r_above_s,active,depth_km,d0_in,km,jm,im, &
       & z_s,d_s,Mg)
    INTEGER*4 , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*8 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),d3(km,im,jm),ssh(im,jm)
    REAL*8 , INTENT(IN) :: dzw(km,im,jm)
    REAL*4 , INTENT(IN) :: depth(km,im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    INTEGER*4 , INTENT(IN) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: r_above_s(im,jm), depth_km
    LOGICAL(KIND=1) , INTENT(IN) :: active(im,jm)
    REAL*8 , INTENT(OUT) :: z_s(im,jm),d_s(im,jm),Mg(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,d3,ssh,dzw,depth,k_below,r_above_s,depth_km,d0_in,active
    !f2py intent (out) z_s,d_s,Mg

    INTEGER*4 :: i,j,k, kmc
    REAL*8 :: r_below,r_above,p_above,p,z,sig8(2),T8(2),S8(2)
    REAL*8 :: rho0 = 1035.d0
    REAL*8 :: rrho0 = 1.d0/1035.d0
    REAL*8 :: grav = 9.80665d0
    REAL*8 :: buoy(100)
    ! REAL*8, ALLOCATABLE :: rhbar(:)
    ! REAL*4, ALLOCATABLE :: rho(:)
    REAL*4 :: b0

    ! ALLOCATE (rhbar(km))
    ! ALLOCATE (rho(km),rhbar(km))

    b0 = -d0_in*grav*rrho0
    do j=1,jm
!!!     !$omp parallel do private(kmc,k,i,r_above,r_below,z,T8,S8,sig8,buoy,p_above,p) shared(active,im,j,kmt,k_below_s,r_above_s,depth_km,d3,b0,z_s,d_s,Mg)
       do i=1,im
          if (active(i,j)) then
             kmc = kmt(i,j)
             k = k_below_s(i,j)
             r_above = r_above_s(i,j)
             r_below = 1.d0 - r_above
             z = r_above*depth(k,i,j) + r_below*depth(k+1,i,j)
             z_s(i,j) = z
             T8 = T(k:k+1,i,j)
             S8 = S(k:k+1,i,j)
             call sigma_n(T8,S8,2,depth_km,sig8)
             d_s(i,j) = r_above*sig8(1) + r_below*sig8(2)
             buoy(1:k) = -grav*rrho0*d3(1:k,i,j)
             p_above = ssh(i,j)*grav - dzw(1,i,j)*buoy(1)
             if(k>1) then
                p_above = p_above - dot_product(.5d0*(buoy(1:k-1)+buoy(2:k)),dzw(2:k,i,j))
             end if
             p = p_above - .5d0*(buoy(k) + b0)*r_below*dzw(k+1,i,j)
             Mg(i,j) = p + b0*z
          end if
       end do
    end do
    end subroutine mginterpolate8


    SUBROUTINE eos_sigman4(T, S, depth_km, rho, n)
    !!----------------------------------------------------------------------
    !!                   ***  ROUTINE eos_insitu  ***
    !!
    !! ** Purpose :   Compute the in unmasked situ density deviation from 1000kg/m^3 ( rho(t,s,z)) from
    !!       potential temperature salinity and depth using an equation of state
    !!       selected in the nameos namelist
    !!
    !! ** Method  :   rho(t,s,z)
    !!                t      TEOS10: CT or EOS80: PT      Celsius
    !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
    !!                depth  reference depth                    km
    !!                rho    potential density              kg/m^3
    !!
    !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
    !!               Note global mean r0(z) is added to perturbation r(t,s,z) = rho(t,s,z) - r0(z)
    !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
    !!
    !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
    !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
    !!
    !!     ln_seos : simplified equation of state
    !!              rho(t,s,z) = -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0)
    !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
    !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
    !!              Vallis like equation: use default values of coefficients
    !!
    !! ** Action  :   compute rho , the in situ density (kg/m^3)
    !!
    !! References :   Roquet et al, Ocean Modelling (2015)
    !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
    !!                TEOS-10 Manual, 2010
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN)                        ::  n
    REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
    REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
    REAL*4, INTENT(in)                           ::  depth_km! reference depth                    [km]
    REAL*4, DIMENSION(n), INTENT(out)            ::  rho     ! deviation of rho from 1000          [kg/m^3]
    !
    INTEGER  ::  i
    REAL*8  :: zt, zh ! local scalars
    REAL*8  :: zs! local scalars
    REAL*8  :: zn1, zn2!   -      -
    REAL*8  :: zn, zn0, zn3!   -      -
    REAL*8  :: zr0
    ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
    REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
    !!----------------------------------------------------------------------
    !
    SELECT CASE( neos )
    !
    CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
       !
       zh = depth_km * 1.d3 * r1_Z0
       ! Define reference profile zr0 to be added to anomaly
       zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
       !$omp parallel do private(zt,zs,zn,zn0,zn1,zn2,zn3)
       DO i=1,n
          !
          zt  = T (i) * r1_T0                           ! temperature
          zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
          !
          zn3 = EOS013*zt   &
             &   + EOS103*zs+EOS003
             !
          zn2 = (EOS022*zt   &
             &   + EOS112*zs+EOS012)*zt   &
             &   + (EOS202*zs+EOS102)*zs+EOS002
             !
          zn1 = (((EOS041*zt   &
             &   + EOS131*zs+EOS031)*zt   &
             &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
             &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
             &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
             !
          zn0 = (((((EOS060*zt   &
             &   + EOS150*zs+EOS050)*zt   &
             &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
             &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
             &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
             &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
             &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
             !
          zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
          !!  Subtract 1000. to give density anomaly
          !!  Add reference profile zr0 to anomaly
             rho(i) = REAL(zn - 1000.d0 + zr0, KIND=4)
    !!----------------------------------------------------------------------
          !
          ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
          !
       END DO
       !$omp end parallel do
       !
    CASE( np_seos )                !==  simplified EOS  ==!
       !
       zh = depth_km * 1.d3
       DO i=1,n
          zt  = T(i) - 10.d0
          zs  = S(i) - 35.d0
          !
          zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
             &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
             &  - rn_nu * zt * zs
             !
          ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
          rho(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
       END DO
       !
    CASE(np_old_eos80)
       zh = depth_km*1000.d0                  ! depth
       !$omp parallel do private(zt,zs,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
       DO i=1,n
           zt = T(i)
           zs = S(i)
           zsr= SQRT( ABS( zs) )        ! square root salinity
           !
           ! compute volumic mass pure water at atm pressure
           zr1= ( ( ( ( 6.536332d-9*zt-1.120083d-6 )*zt+1.001685d-4 )*zt   &
                &                          -9.095290d-3 )*zt+6.793952d-2 )*zt+999.842594d0
           ! seawater volumic mass atm pressure
           zr2= ( ( ( 5.3875d-9*zt-8.2467d-7 ) *zt+7.6438d-5 ) *zt   &
                &                                         -4.0899d-3 ) *zt+0.824493d0
           zr3= ( -1.6546d-6*zt+1.0227d-4 )    *zt-5.72466d-3
           zr4= 4.8314d-4
           !
           ! potential volumic mass (reference to the surface)
           zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
           !
           ! add the compression terms
           ze = ( -3.508914d-8*zt-1.248266d-8 ) *zt-2.595994d-6
           zbw= (  1.296821d-6*zt-5.782165d-9 ) *zt+1.045941d-4
           zb = zbw + ze * zs
           !
           zd = -2.042967d-2
           zc =   (-7.267926d-5*zt+2.598241d-3 ) *zt+0.1571896
           zaw= ( ( 5.939910d-6*zt+2.512549d-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
           za = ( zd*zsr + zc ) *zs + zaw
           !
           zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
           za1= ( (   2.326469d-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
           zkw= ( ( (-1.361629d-4*zt-1.852732d-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
           zk0= ( zb1*zsr + za1 )*zs + zkw
           !
           ! unmasked in situ density anomaly
           rho(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0, KIND=4)
        END DO
       !$omp end parallel do
    CASE DEFAULT
       PRINT *,'eos_sigman called without neos set by previous call of eos_init'
       STOP
    END SELECT
    !
 END SUBROUTINE eos_sigman4


  subroutine hunt4(kmt,nmax,xx,nlo,y)
    INTEGER , INTENT(IN) :: kmt,nmax
    REAL*4 , INTENT(IN) :: xx(nmax),y(kmt)
    INTEGER, INTENT(INOUT) :: nlo(kmt)
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.4.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    do k=1,kmt
       n=nlo(k)
       inc = 1
       ! #ifdef debug
       !        print *,'k ',k,'n ',n,'y on entry',y(k)
       !        print *,'nlo on entry',nlo(k)
       ! #endif
       if (y(k).ge.xx(n)) then
1         nhi = n+inc
          if(nhi.gt.nmax)then
             nhi=nmax
          else if (y(k).ge.xx(nhi)) then
             n=nhi
             inc=inc+inc
             goto 1
          end if
       else
          nhi=n
2         n=nhi-inc
          if(n.lt.1)then
             n=1
          else if (y(k).lt.xx(n)) then
             nhi=n
             inc=inc+inc
             goto 2
          endif
       end if
3      if(nhi.ne.n+1)then
          nmean=(nhi+n)/2
          if(y(k).gt.xx(nmean))then
             n=nmean
          else
             nhi=nmean
          endif
          goto 3
       endif
       nlo(k)=n
       ! #ifdef debug
       !        print *,'nlo on exit',nlo(k)
       ! #endif
    end do
  end subroutine hunt4

  subroutine huntfrac4(kmt,nmax,xx,nlo,y,dfy)
    INTEGER , INTENT(IN) :: kmt,nmax
    REAL*4 , INTENT(IN) :: xx(nmax),y(kmt)
    INTEGER, INTENT(INOUT) :: nlo(kmt)
    REAL*4 , INTENT(OUT) :: dfy(kmt)
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !f2py intent (out) dfy
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) and fractions dfy(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !     and y(k) = xx(nlo(k))*(1.-dfy(k)) + xx(nlo(k)+1)*dfy(k)
    !
    !     Based on  HUNT routine. Numerical recipes,3.4.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    do k=1,kmt
       n=nlo(k)
       inc = 1
       ! #ifdef debug
       !        print *,'k ',k,'n ',n,'y on entry',y(k)
       !        print *,'nlo on entry',nlo(k)
       ! #endif
       if (y(k).ge.xx(n)) then
1         nhi = n+inc
          if(nhi.gt.nmax)then
             nhi=nmax
          else if (y(k).ge.xx(nhi)) then
             n=nhi
             inc=inc+inc
             goto 1
          end if
       else
          nhi=n
2         n=nhi-inc
          if(n.lt.1)then
             n=1
          else if (y(k).lt.xx(n)) then
             nhi=n
             inc=inc+inc
             goto 2
          endif
       end if
3      if(nhi.ne.n+1)then
          nmean=(nhi+n)/2
          if(y(k).gt.xx(nmean))then
             n=nmean
          else
             nhi=nmean
          endif
          goto 3
       endif
       nlo(k)=n
       dfy(k) = (y(k) - xx(n))/(xx(n+1) - xx(n))
       ! #ifdef debug
       !        print *,'nlo on exit',nlo(k)
       ! #endif
    end do
  end subroutine huntfrac4


  subroutine huntfrac48(kmt,nmax,xx,nlo,y,dfy)
    INTEGER , INTENT(IN) :: kmt,nmax
    REAL*8 , INTENT(IN) :: xx(nmax)
    REAL*4 , INTENT(IN) :: y(kmt)
    INTEGER, INTENT(INOUT) :: nlo(kmt)
    REAL*4 , INTENT(OUT) :: dfy(kmt)
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !f2py intent (out) dfy
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) and fractions dfy(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !     and y(k) = xx(nlo(k))*(1.-dfy(k)) + xx(nlo(k)+1)*dfy(k)
    !
    !     Based on  HUNT routine. Numerical recipes,3.4.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    do k=1,kmt
       n=nlo(k)
       inc = 1
       ! #ifdef debug
       !        print *,'k ',k,'n ',n,'y on entry',y(k)
       !        print *,'nlo on entry',nlo(k)
       ! #endif
       if (y(k).ge.xx(n)) then
1         nhi = n+inc
          if(nhi.gt.nmax)then
             nhi=nmax
          else if (y(k).ge.xx(nhi)) then
             n=nhi
             inc=inc+inc
             goto 1
          end if
       else
          nhi=n
2         n=nhi-inc
          if(n.lt.1)then
             n=1
          else if (y(k).lt.xx(n)) then
             nhi=n
             inc=inc+inc
             goto 2
          endif
       end if
3      if(nhi.ne.n+1)then
          nmean=(nhi+n)/2
          if(y(k).gt.xx(nmean))then
             n=nmean
          else
             nhi=nmean
          endif
          goto 3
       endif
       nlo(k)=n
       dfy(k) = (y(k) - xx(n))/(xx(n+1) - xx(n))
       ! #ifdef debug
       !        print *,'nlo on exit',nlo(k)
       ! #endif
    end do
  end subroutine huntfrac48

  subroutine hunt14(nmax,xx,nlo,y)
    INTEGER , INTENT(IN) :: nmax
    REAL*4 , INTENT(IN) :: xx(nmax),y
    INTEGER, INTENT(INOUT) :: nlo
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.4.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    n=nlo
    inc = 1
    if (y.ge.xx(n)) then
1      nhi = n+inc
       if(nhi.gt.nmax)then
          nhi=nmax
       else if (y.ge.xx(nhi)) then
          n=nhi
          inc=inc+inc
          goto 1
       end if
    else
       nhi=n
2      n=nhi-inc
       if(n.lt.1)then
          n=1
       else if (y.lt.xx(n)) then
          nhi=n
          inc=inc+inc
          goto 2
       endif
    end if
3   if(nhi.ne.n+1)then
       nmean=(nhi+n)/2
       if(y.gt.xx(nmean))then
          n=nmean
       else
          nhi=nmean
       endif
       goto 3
    endif
    nlo=n
  end subroutine hunt14

  subroutine hunt18(nmax,xx,nlo,y)
    INTEGER , INTENT(IN) :: nmax
    REAL*8 , INTENT(IN) :: xx(nmax),y
    INTEGER, INTENT(INOUT) :: nlo
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.8.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    n=nlo
    inc = 1
    if (y.ge.xx(n)) then
1      nhi = n+inc
       if(nhi.gt.nmax)then
          nhi=nmax
       else if (y.ge.xx(nhi)) then
          n=nhi
          inc=inc+inc
          goto 1
       end if
    else
       nhi=n
2      n=nhi-inc
       if(n.lt.1)then
          n=1
       else if (y.lt.xx(n)) then
          nhi=n
          inc=inc+inc
          goto 2
       endif
    end if
3   if(nhi.ne.n+1)then
       nmean=(nhi+n)/2
       if(y.gt.xx(nmean))then
          n=nmean
       else
          nhi=nmean
       endif
       goto 3
    endif
    nlo=n
  end subroutine hunt18

  subroutine huntc18(nmax,xx,nlo_in,n,y)
    INTEGER , INTENT(IN) :: nmax
    REAL*8 , INTENT(IN) :: xx(nmax),y
    INTEGER, INTENT(IN) :: nlo_in
    INTEGER, INTENT(OUT) :: n
    !f2py intent (in) kmt,nmax,xx,y,nlo_in
    !f2py intent (out) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.8.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: inc,nhi,nmean,k,iterations,itmax=100
    if (y.gt.xx(nmax)) then
       n = nmax
       return
    else if (y.lt.xx(1)) then
       n = 0
       return
    else
       n=nlo_in+1
       if ((n.lt.1).or.(n.ge.nmax)) then
          n = nmax/2
       end if
    end if
    iterations = 0
    inc = 1
    if (y.ge.xx(n)) then
1      nhi = n+inc
       iterations = iterations + 1
       if(nhi.gt.nmax)then
          nhi=nmax
       else if (y.ge.xx(nhi)) then
          iterations = iterations + 1
          if (iterations>itmax) then
             n=-1
             return
          endif
          n=nhi
          inc=inc+inc
          goto 1
       end if
    else
       nhi=n
2      n=nhi-inc
       if(n.lt.1)then
          n=1
       else if (y.lt.xx(n)) then
          iterations = iterations + 1
          if (iterations>itmax) then
             n=-1
             return
          endif
          nhi=n
          inc=inc+inc
          goto 2
       endif
    end if
3   if(nhi.ne.n+1)then
       nmean=(nhi+n)/2
       if(y.gt.xx(nmean))then
          iterations = iterations + 1
          if (iterations>itmax) then
             n=-1
             return
          endif
          n=nmean
       else
          nhi=nmean
       endif
       goto 3
    endif
    return
  end subroutine huntc18

  subroutine huntc14(nmax,xx,nlo_in,n,y)
    INTEGER , INTENT(IN) :: nmax
    REAL*4 , INTENT(IN) :: xx(nmax),y
    INTEGER, INTENT(IN) :: nlo_in
    INTEGER, INTENT(OUT) :: n
    !f2py intent (in) kmt,nmax,xx,y,nlo_in
    !f2py intent (out) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.8.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: inc,nhi,nmean,k,iterations,itmax=100
    if (y.gt.xx(nmax)) then
       n = nmax
       return
    else if (y.lt.xx(1)) then
       n = 0
       return
    else
       n=nlo_in
       if ((n.lt.1).or.(n.ge.nmax)) then
          n = nmax/2
       end if
    end if
    iterations = 0
    inc = 1
    if (y.ge.xx(n)) then
1      nhi = n+inc
       iterations = iterations + 1
       if (iterations>itmax) then
          n=-1
          return
       endif
       if(nhi.gt.nmax)then
          nhi=nmax
       else if (y.ge.xx(nhi)) then
          n=nhi
          inc=inc+inc
          goto 1
       end if
    else
       nhi=n
2      n=nhi-inc
       iterations = iterations + 1
       if (iterations>itmax) then
          n=-1
          return
       endif
       if(n.lt.1)then
          n=1
       else if (y.lt.xx(n)) then
          nhi=n
          inc=inc+inc
          goto 2
       endif
    end if
3   if(nhi.ne.n+1)then
       iterations = iterations + 1
       if (iterations>itmax) then
          n=-1
          return
       endif
       nmean=(nhi+n)/2
       if(y.gt.xx(nmean))then
          n=nmean
       else
          nhi=nmean
       endif
       goto 3
    endif
    return
  end subroutine huntc14

  subroutine hunt8(kmt,nmax,xx,nlo,y)
    INTEGER , INTENT(IN) :: kmt,nmax
    REAL*8 , INTENT(IN) :: xx(nmax),y(kmt)
    INTEGER, INTENT(INOUT) :: nlo(kmt)
    !f2py intent (in) kmt,nmax,xx,y
    !f2py intent (inout) nlo
    !------------------------------------------------------------------------------------------------------
    !     finds array of indices nlo(1...kmt) such that
    !     increasing table xx(nlo(k)),(k+1)  brackets y(k) in input array  y(1...kmt)
    !
    !     Based on  HUNT routine. Numerical recipes,3.4.
    !     but assumes ascending xx, and works on vector array y(1...kmt)
    !
    !     Uses input values of nlo(..)  as first guess for new values.
    !     Note no error checking, so nlo(..) needs to be initialised
    !     to nonzero values
    !------------------------------------------------------------------------------------------------------

    INTEGER :: n,inc,nhi,nmean,k
    do k=1,kmt
       n=nlo(k)
       inc = 1
       ! #ifdef debug
       !        print *,'k ',k,'n ',n,'y on entry',y(k)
       !        print *,'nlo on entry',nlo(k)
       ! #endif
       if (y(k).ge.xx(n)) then
1         nhi = n+inc
          if(nhi.gt.nmax)then
             nhi=nmax
          else if (y(k).ge.xx(nhi)) then
             n=nhi
             inc=inc+inc
             goto 1
          end if
       else
          nhi=n
2         n=nhi-inc
          if(n.lt.1)then
             n=1
          else if (y(k).lt.xx(n)) then
             nhi=n
             inc=inc+inc
             goto 2
          endif
       end if
3      if(nhi.ne.n+1)then
          nmean=(nhi+n)/2
          if(y(k).gt.xx(nmean))then
             n=nmean
          else
             nhi=nmean
          endif
          goto 3
       endif
       nlo(k)=n
       ! #ifdef debug
       !        print *,'nlo on exit',nlo(k)
       ! #endif
    end do
  end subroutine hunt8


  subroutine convect8(dens,dh,pass,kmc)
    integer , INTENT(IN) :: kmc
    real*8 , intent(inout) :: dens(kmc), pass(kmc)
    real*8 , intent(in) :: dh(kmc)

    integer kt,kb,k,lb,la
    logical chk_la,chk_lb
    real*8 ru,rl,dz1,dz2,tsm,tmx,zsm,psm,pmx
    !
    !     search for unstable regions starting from the top
    !
    if (minval(dens(2:kmc)-dens(1:kmc-1)).ge.0.d0) return
    kt = 1
    kb = 2
    do while (kt .lt. kmc)
       ru = dens(kt)
       rl = dens(kb)
       !
       !     sum the first pair found in an unstable region
       !
       if (ru .gt. rl) then
          chk_la = .true.
          chk_lb = .true.
          dz1 = dh(kt)
          dz2 = dh(kb)
          zsm = dz1 + dz2
          tsm = ru*dz1 + rl*dz2
          tmx = tsm/zsm
          psm = pass(kt)*dz1 + pass(kb)*dz2
          pmx = psm/zsm
          !
          do while (chk_lb .or. chk_la)
             !
             !     check for an unstable level (lb) below kb
             !
             if (kb .ge. kmc) chk_lb = .false.
             do while (chk_lb)
                chk_lb = .false.
                lb = kb + 1
                ru = tmx
                rl = dens(lb)
                if (ru .gt. rl) then
                   !
                   !     add new level to sums
                   !
                   kb = lb
                   dz2 = dh(kb)
                   zsm = zsm + dz2
                   tsm = tsm + dens(kb)*dz2
                   tmx = tsm/zsm
                   psm = psm + pass(kb)*dz2
                   pmx = psm/zsm
                   chk_la = .true.
                   if (kb .lt. kmc) chk_lb = .true.
                endif
             enddo
             !
             !     check for an unstable level (la) above kt
             !     to get equivalent of original Rahmstorf code, uncomment the next line
             !     chk_la = .true.
             if (kt .le. 1) chk_la = .false.
             do while (chk_la)
                chk_la = .false.
                la = kt - 1
                ru = dens(la)
                rl = tmx
                if (ru .gt. rl) then
                   !
                   !     add new level to sums
                   !
                   kt = la
                   dz1 = dh(kt)
                   zsm = zsm + dz1
                   tsm = tsm + dens(kt)*dz1
                   tmx = tsm/zsm
                   psm = psm + pass(kt)*dz1
                   pmx = psm/zsm
                   chk_lb = .true.
                   !
                   !     to get equivalent of original Rahmstorf code, comment out the next line
                   !
                   if (kt .gt. 1) chk_la = .true.
                endif
             enddo
          enddo
          !
          !     mix all tracers from kt to kb
          !
          do k=kt,kb
             dens(k) = tmx
             pass(k) = pmx
          enddo

          !
          !     some possible diagnostics
          !     ktot = ktot + kb - kt + 1
          !     if (kt .eq. 1) kven = kb
          !
          kt = kb + 1
       else
          kt = kb
       endif
       !
       !     continue the search for other unstable regions
       !
       kb = kt + 1
    enddo
    !
    return
  end subroutine convect8

  subroutine convect4(dens,dh,pass,kmc)
    integer , INTENT(IN) :: kmc
    real*4 , intent(inout) :: dens(kmc), pass(kmc)
    real*8 , intent(in) :: dh(kmc)
    !f2py intent(in) dh,kmc
    !f2py intent(inout) dens,pass

    integer kt,kb,k,lb,la
    logical chk_la,chk_lb
    real*8 ru,rl,dz1,dz2,tsm,tmx,zsm,psm,pmx
    !
    !     search for unstable regions starting from the top
    !
    if (minval(dens(2:kmc)-dens(1:kmc-1)).ge.0.d0) return
    kt = 1
    kb = 2
    do while (kt .lt. kmc)
       ru = dens(kt)
       rl = dens(kb)
       !
       !     sum the first pair found in an unstable region
       !
       if (ru .gt. rl) then
          chk_la = .true.
          chk_lb = .true.
          dz1 = dh(kt)
          dz2 = dh(kb)
          zsm = dz1 + dz2
          tsm = ru*dz1 + rl*dz2
          tmx = tsm/zsm
          psm = pass(kt)*dz1 + pass(kb)*dz2
          pmx = psm/zsm
          !
          do while (chk_lb .or. chk_la)
             !
             !     check for an unstable level (lb) below kb
             !
             if (kb .ge. kmc) chk_lb = .false.
             do while (chk_lb)
                chk_lb = .false.
                lb = kb + 1
                ru = tmx
                rl = dens(lb)
                if (ru .gt. rl) then
                   !
                   !     add new level to sums
                   !
                   kb = lb
                   dz2 = dh(kb)
                   zsm = zsm + dz2
                   tsm = tsm + dens(kb)*dz2
                   tmx = tsm/zsm
                   psm = psm + pass(kb)*dz2
                   pmx = psm/zsm
                   chk_la = .true.
                   if (kb .lt. kmc) chk_lb = .true.
                endif
             enddo
             !
             !     check for an unstable level (la) above kt
             !     to get equivalent of original Rahmstorf code, uncomment the next line
             !     chk_la = .true.
             if (kt .le. 1) chk_la = .false.
             do while (chk_la)
                chk_la = .false.
                la = kt - 1
                ru = dens(la)
                rl = tmx
                if (ru .gt. rl) then
                   !
                   !     add new level to sums
                   !
                   kt = la
                   dz1 = dh(kt)
                   zsm = zsm + dz1
                   tsm = tsm + dens(kt)*dz1
                   tmx = tsm/zsm
                   psm = psm + pass(kt)*dz1
                   pmx = psm/zsm
                   chk_lb = .true.
                   !
                   !     to get equivalent of original Rahmstorf code, comment out the next line
                   !
                   if (kt .gt. 1) chk_la = .true.
                endif
             enddo
          enddo
          !
          !     mix all tracers from kt to kb
          !
          do k=kt,kb
             dens(k) = tmx
             pass(k) = pmx
          enddo

          !
          !     some possible diagnostics
          !     ktot = ktot + kb - kt + 1
          !     if (kt .eq. 1) kven = kb
          !
          kt = kb + 1
       else
          kt = kb
       endif
       !
       !     continue the search for other unstable regions
       !
       kb = kt + 1
    enddo
    !
    return
  end subroutine convect4

   SUBROUTINE eos_init(neos_in)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_init  ***
      !!
      !! ** Purpose :   initializations for the equation of state
      !!
      !! ** Method  :   Read the namelist nameos and control the parameters
      !!----------------------------------------------------------------------
      ! NAMELIST/nameos/ rn_a0, rn_b0, rn_lambda1, rn_mu1,   &
      !    &                                             rn_lambda2, rn_mu2, rn_nu
      !!----------------------------------------------------------------------
      !
      !
      INTEGER*4 neos_in


      neos = neos_in
      SELECT CASE( neos )         ! check option
      !
      CASE( np_teos10 )                       !==  polynomial TEOS-10  ==!
         rdeltaS = 32.d0
         r1_S0  = 0.875d0/35.16504d0
         r1_T0  = 1.d0/40.d0
         r1_Z0  = 1.d-4
         !
         EOS000 = 8.0189615746d+02
         EOS100 = 8.6672408165d+02
         EOS200 = -1.7864682637d+03
         EOS300 = 2.0375295546d+03
         EOS400 = -1.2849161071d+03
         EOS500 = 4.3227585684d+02
         EOS600 = -6.0579916612d+01
         EOS010 = 2.6010145068d+01
         EOS110 = -6.5281885265d+01
         EOS210 = 8.1770425108d+01
         EOS310 = -5.6888046321d+01
         EOS410 = 1.7681814114d+01
         EOS510 = -1.9193502195d0
         EOS020 = -3.7074170417d+01
         EOS120 = 6.1548258127d+01
         EOS220 = -6.0362551501d+01
         EOS320 = 2.9130021253d+01
         EOS420 = -5.4723692739d0
         EOS030 = 2.1661789529d+01
         EOS130 = -3.3449108469d+01
         EOS230 = 1.9717078466d+01
         EOS330 = -3.1742946532d0
         EOS040 = -8.3627885467d0
         EOS140 = 1.1311538584d+01
         EOS240 = -5.3563304045d0
         EOS050 = 5.4048723791d-01
         EOS150 = 4.8169980163d-01
         EOS060 = -1.9083568888d-01
         EOS001 = 1.9681925209d+01
         EOS101 = -4.2549998214d+01
         EOS201 = 5.0774768218d+01
         EOS301 = -3.0938076334d+01
         EOS401 = 6.6051753097d0
         EOS011 = -1.3336301113d+01
         EOS111 = -4.4870114575d0
         EOS211 = 5.0042598061d0
         EOS311 = -6.5399043664d-01
         EOS021 = 6.7080479603d0
         EOS121 = 3.5063081279d0
         EOS221 = -1.8795372996d0
         EOS031 = -2.4649669534d0
         EOS131 = -5.5077101279d-01
         EOS041 = 5.5927935970d-01
         EOS002 = 2.0660924175d0
         EOS102 = -4.9527603989d0
         EOS202 = 2.5019633244d0
         EOS012 = 2.0564311499d0
         EOS112 = -2.1311365518d-01
         EOS022 = -1.2419983026d0
         EOS003 = -2.3342758797d-02
         EOS103 = -1.8507636718d-02
         EOS013 = 3.7969820455d-01
         !
         ALP000 = -6.5025362670d-01
         ALP100 = 1.6320471316d0
         ALP200 = -2.0442606277d0
         ALP300 = 1.4222011580d0
         ALP400 = -4.4204535284d-01
         ALP500 = 4.7983755487d-02
         ALP010 = 1.8537085209d0
         ALP110 = -3.0774129064d0
         ALP210 = 3.0181275751d0
         ALP310 = -1.4565010626d0
         ALP410 = 2.7361846370d-01
         ALP020 = -1.6246342147d0
         ALP120 = 2.5086831352d0
         ALP220 = -1.4787808849d0
         ALP320 = 2.3807209899d-01
         ALP030 = 8.3627885467d-01
         ALP130 = -1.1311538584d0
         ALP230 = 5.3563304045d-01
         ALP040 = -6.7560904739d-02
         ALP140 = -6.0212475204d-02
         ALP050 = 2.8625353333d-02
         ALP001 = 3.3340752782d-01
         ALP101 = 1.1217528644d-01
         ALP201 = -1.2510649515d-01
         ALP301 = 1.6349760916d-02
         ALP011 = -3.3540239802d-01
         ALP111 = -1.7531540640d-01
         ALP211 = 9.3976864981d-02
         ALP021 = 1.8487252150d-01
         ALP121 = 4.1307825959d-02
         ALP031 = -5.5927935970d-02
         ALP002 = -5.1410778748d-02
         ALP102 = 5.3278413794d-03
         ALP012 = 6.2099915132d-02
         ALP003 = -9.4924551138d-03
         !
         BET000 = 1.0783203594d+01
         BET100 = -4.4452095908d+01
         BET200 = 7.6048755820d+01
         BET300 = -6.3944280668d+01
         BET400 = 2.6890441098d+01
         BET500 = -4.5221697773d0
         BET010 = -8.1219372432d-01
         BET110 = 2.0346663041d0
         BET210 = -2.1232895170d0
         BET310 = 8.7994140485d-01
         BET410 = -1.1939638360d-01
         BET020 = 7.6574242289d-01
         BET120 = -1.5019813020d0
         BET220 = 1.0872489522d0
         BET320 = -2.7233429080d-01
         BET030 = -4.1615152308d-01
         BET130 = 4.9061350869d-01
         BET230 = -1.1847737788d-01
         BET040 = 1.4073062708d-01
         BET140 = -1.3327978879d-01
         BET050 = 5.9929880134d-03
         BET001 = -5.2937873009d-01
         BET101 = 1.2634116779d0
         BET201 = -1.1547328025d0
         BET301 = 3.2870876279d-01
         BET011 = -5.5824407214d-02
         BET111 = 1.2451933313d-01
         BET211 = -2.4409539932d-02
         BET021 = 4.3623149752d-02
         BET121 = -4.6767901790d-02
         BET031 = -6.8523260060d-03
         BET002 = -6.1618945251d-02
         BET102 = 6.2255521644d-02
         BET012 = -2.6514181169d-03
         BET003 = -2.3025968587d-04
         !
         PEN000 = -9.8409626043d0
         PEN100 = 2.1274999107d+01
         PEN200 = -2.5387384109d+01
         PEN300 = 1.5469038167d+01
         PEN400 = -3.3025876549d0
         PEN010 = 6.6681505563d0
         PEN110 = 2.2435057288d0
         PEN210 = -2.5021299030d0
         PEN310 = 3.2699521832d-01
         PEN020 = -3.3540239802d0
         PEN120 = -1.7531540640d0
         PEN220 = 9.3976864981d-01
         PEN030 = 1.2324834767d0
         PEN130 = 2.7538550639d-01
         PEN040 = -2.7963967985d-01
         PEN001 = -1.3773949450d0
         PEN101 = 3.3018402659d0
         PEN201 = -1.6679755496d0
         PEN011 = -1.3709540999d0
         PEN111 = 1.4207577012d-01
         PEN021 = 8.2799886843d-01
         PEN002 = 1.7507069098d-02
         PEN102 = 1.3880727538d-02
         PEN012 = -2.8477365341d-01
         !
         APE000 = -1.6670376391d-01
         APE100 = -5.6087643219d-02
         APE200 = 6.2553247576d-02
         APE300 = -8.1748804580d-03
         APE010 = 1.6770119901d-01
         APE110 = 8.7657703198d-02
         APE210 = -4.6988432490d-02
         APE020 = -9.2436260751d-02
         APE120 = -2.0653912979d-02
         APE030 = 2.7963967985d-02
         APE001 = 3.4273852498d-02
         APE101 = -3.5518942529d-03
         APE011 = -4.1399943421d-02
         APE002 = 7.1193413354d-03
         !
         BPE000 = 2.6468936504d-01
         BPE100 = -6.3170583896d-01
         BPE200 = 5.7736640125d-01
         BPE300 = -1.6435438140d-01
         BPE010 = 2.7912203607d-02
         BPE110 = -6.2259666565d-02
         BPE210 = 1.2204769966d-02
         BPE020 = -2.1811574876d-02
         BPE120 = 2.3383950895d-02
         BPE030 = 3.4261630030d-03
         BPE001 = 4.1079296834d-02
         BPE101 = -4.1503681096d-02
         BPE011 = 1.7676120780d-03
         BPE002 = 1.7269476440d-04
         !
      CASE( np_eos80 )                        !==  polynomial EOS-80 formulation  ==!
         !
         rdeltaS = 20.d0
         r1_S0  = 1.d0/40.d0
         r1_T0  = 1.d0/40.d0
         r1_Z0  = 1.d-4
         !
         EOS000 = 9.5356891948d+02
         EOS100 = 1.7136499189d+02
         EOS200 = -3.7501039454d+02
         EOS300 = 5.1856810420d+02
         EOS400 = -3.7264470465d+02
         EOS500 = 1.4302533998d+02
         EOS600 = -2.2856621162d+01
         EOS010 = 1.0087518651d+01
         EOS110 = -1.3647741861d+01
         EOS210 = 8.8478359933d0
         EOS310 = -7.2329388377d0
         EOS410 = 1.4774410611d0
         EOS510 = 2.0036720553d-01
         EOS020 = -2.5579830599d+01
         EOS120 = 2.4043512327d+01
         EOS220 = -1.6807503990d+01
         EOS320 = 8.3811577084d0
         EOS420 = -1.9771060192d0
         EOS030 = 1.6846451198d+01
         EOS130 = -2.1482926901d+01
         EOS230 = 1.0108954054d+01
         EOS330 = -6.2675951440d-01
         EOS040 = -8.0812310102d0
         EOS140 = 1.0102374985d+01
         EOS240 = -4.8340368631d0
         EOS050 = 1.2079167803d0
         EOS150 = 1.1515380987d-01
         EOS060 = -2.4520288837d-01
         EOS001 = 1.0748601068d+01
         EOS101 = -1.7817043500d+01
         EOS201 = 2.2181366768d+01
         EOS301 = -1.6750916338d+01
         EOS401 = 4.1202230403d0
         EOS011 = -1.5852644587d+01
         EOS111 = -7.6639383522d-01
         EOS211 = 4.1144627302d0
         EOS311 = -6.6955877448d-01
         EOS021 = 9.9994861860d0
         EOS121 = -1.9467067787d-01
         EOS221 = -1.2177554330d0
         EOS031 = -3.4866102017d0
         EOS131 = 2.2229155620d-01
         EOS041 = 5.9503008642d-01
         EOS002 = 1.0375676547d0
         EOS102 = -3.4249470629d0
         EOS202 = 2.0542026429d0
         EOS012 = 2.1836324814d0
         EOS112 = -3.4453674320d-01
         EOS022 = -1.2548163097d0
         EOS003 = 1.8729078427d-02
         EOS103 = -5.7238495240d-02
         EOS013 = 3.8306136687d-01
         !
         ALP000 = -2.5218796628d-01
         ALP100 = 3.4119354654d-01
         ALP200 = -2.2119589983d-01
         ALP300 = 1.8082347094d-01
         ALP400 = -3.6936026529d-02
         ALP500 = -5.0091801383d-03
         ALP010 = 1.2789915300d0
         ALP110 = -1.2021756164d0
         ALP210 = 8.4037519952d-01
         ALP310 = -4.1905788542d-01
         ALP410 = 9.8855300959d-02
         ALP020 = -1.2634838399d0
         ALP120 = 1.6112195176d0
         ALP220 = -7.5817155402d-01
         ALP320 = 4.7006963580d-02
         ALP030 = 8.0812310102d-01
         ALP130 = -1.0102374985d0
         ALP230 = 4.8340368631d-01
         ALP040 = -1.5098959754d-01
         ALP140 = -1.4394226233d-02
         ALP050 = 3.6780433255d-02
         ALP001 = 3.9631611467d-01
         ALP101 = 1.9159845880d-02
         ALP201 = -1.0286156825d-01
         ALP301 = 1.6738969362d-02
         ALP011 = -4.9997430930d-01
         ALP111 = 9.7335338937d-03
         ALP211 = 6.0887771651d-02
         ALP021 = 2.6149576513d-01
         ALP121 = -1.6671866715d-02
         ALP031 = -5.9503008642d-02
         ALP002 = -5.4590812035d-02
         ALP102 = 8.6134185799d-03
         ALP012 = 6.2740815484d-02
         ALP003 = -9.5765341718d-03
         !
         BET000 = 2.1420623987d0
         BET100 = -9.3752598635d0
         BET200 = 1.9446303907d+01
         BET300 = -1.8632235232d+01
         BET400 = 8.9390837485d0
         BET500 = -1.7142465871d0
         BET010 = -1.7059677327d-01
         BET110 = 2.2119589983d-01
         BET210 = -2.7123520642d-01
         BET310 = 7.3872053057d-02
         BET410 = 1.2522950346d-02
         BET020 = 3.0054390409d-01
         BET120 = -4.2018759976d-01
         BET220 = 3.1429341406d-01
         BET320 = -9.8855300959d-02
         BET030 = -2.6853658626d-01
         BET130 = 2.5272385134d-01
         BET230 = -2.3503481790d-02
         BET040 = 1.2627968731d-01
         BET140 = -1.2085092158d-01
         BET050 = 1.4394226233d-03
         BET001 = -2.2271304375d-01
         BET101 = 5.5453416919d-01
         BET201 = -6.2815936268d-01
         BET301 = 2.0601115202d-01
         BET011 = -9.5799229402d-03
         BET111 = 1.0286156825d-01
         BET211 = -2.5108454043d-02
         BET021 = -2.4333834734d-03
         BET121 = -3.0443885826d-02
         BET031 = 2.7786444526d-03
         BET002 = -4.2811838287d-02
         BET102 = 5.1355066072d-02
         BET012 = -4.3067092900d-03
         BET003 = -7.1548119050d-04
         !
         PEN000 = -5.3743005340d0
         PEN100 = 8.9085217499d0
         PEN200 = -1.1090683384d+01
         PEN300 = 8.3754581690d0
         PEN400 = -2.0601115202d0
         PEN010 = 7.9263222935d0
         PEN110 = 3.8319691761d-01
         PEN210 = -2.0572313651d0
         PEN310 = 3.3477938724d-01
         PEN020 = -4.9997430930d0
         PEN120 = 9.7335338937d-02
         PEN220 = 6.0887771651d-01
         PEN030 = 1.7433051009d0
         PEN130 = -1.1114577810d-01
         PEN040 = -2.9751504321d-01
         PEN001 = -6.9171176978d-01
         PEN101 = 2.2832980419d0
         PEN201 = -1.3694684286d0
         PEN011 = -1.4557549876d0
         PEN111 = 2.2969116213d-01
         PEN021 = 8.3654420645d-01
         PEN002 = -1.4046808820d-02
         PEN102 = 4.2928871430d-02
         PEN012 = -2.8729602515d-01
         !
         APE000 = -1.9815805734d-01
         APE100 = -9.5799229402d-03
         APE200 = 5.1430784127d-02
         APE300 = -8.3694846809d-03
         APE010 = 2.4998715465d-01
         APE110 = -4.8667669469d-03
         APE210 = -3.0443885826d-02
         APE020 = -1.3074788257d-01
         APE120 = 8.3359333577d-03
         APE030 = 2.9751504321d-02
         APE001 = 3.6393874690d-02
         APE101 = -5.7422790533d-03
         APE011 = -4.1827210323d-02
         APE002 = 7.1824006288d-03
         !
         BPE000 = 1.1135652187d-01
         BPE100 = -2.7726708459d-01
         BPE200 = 3.1407968134d-01
         BPE300 = -1.0300557601d-01
         BPE010 = 4.7899614701d-03
         BPE110 = -5.1430784127d-02
         BPE210 = 1.2554227021d-02
         BPE020 = 1.2166917367d-03
         BPE120 = 1.5221942913d-02
         BPE030 = -1.3893222263d-03
         BPE001 = 2.8541225524d-02
         BPE101 = -3.4236710714d-02
         BPE011 = 2.8711395266d-03
         BPE002 = 5.3661089288d-04
         !
         IF(neos == np_teos10) r1_S0  = 0.875d0/35.16504d0   ! Used to convert CT to potential temperature when using bulk formulae
                                         !   (eos_pot_from_CT)
      CASE( np_seos )                        !==  Simplified EOS     ==!

      END SELECT
      !
      r1_rho0     = 1.d0 / rho0
      ! rcp         = 3991.86795711963d0      !: heat capacity     [J/K]
      ! rho0_rcp    = rho0 * rcp
      ! r1_rcp      = 1.d0 / rcp
      ! r1_rho0_rcp = 1.d0 / rho0_rcp
      ! !
      !
   END SUBROUTINE eos_init


end module interp
