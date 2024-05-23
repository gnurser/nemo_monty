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

  subroutine interpolate4(kmt,T,S,rho,rho0,d0_in,km,jm,im, &
       & k_below_s,r_above_s,T_s,S_s,outcropmask,groundmask)
    INTEGER , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*4 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),rho(km,im,jm),rho0(km,im,jm)
    INTEGER , INTENT(OUT) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    REAL*8 , INTENT(OUT) :: r_above_s(im,jm),T_s(im,jm),S_s(im,jm)
    LOGICAL(KIND=1), INTENT(OUT) :: outcropmask(im,jm),groundmask(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,rho,rho0,d0_in
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
       !$omp parallel do private(i,kmc,k,d,r_above,dd,r_below) shared(im,j,kmt,kin,rho,rho0,d0,outcropmask,groundmask,T_s,S_s)
       do i=1,im
          kmc = kmt(i,j)
          if (kmc==0) cycle
          d(1:kmc) = rho(1:kmc,i,j) - rho0(1:kmc,i,j)
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

  subroutine mginterpolate4(kmt,T,S,rho,rho0,ssh,dzw,depth, &
       & k_below_s,r_above_s,active,depth_km,d0_in,km,jm,im, &
       & z_s,Mg)
    INTEGER*4 , INTENT(IN) :: kmt(im,jm),km,jm,im
    REAL*4 , INTENT(IN) :: T(km,im,jm),S(km,im,jm),rho(km,im,jm),rho0(km,im,jm),ssh(im,jm)
    REAL*8 , INTENT(IN) :: dzw(km,im,jm)
    REAL*4 , INTENT(IN) :: depth(km,im,jm)
    REAL*8 , INTENT(IN) :: d0_in
    INTEGER*4 , INTENT(IN) :: k_below_s(im,jm)
    REAL*8 , INTENT(IN) :: r_above_s(im,jm), depth_km
    LOGICAL(KIND=1) , INTENT(IN) :: active(im,jm)
    REAL*8 , INTENT(OUT) :: z_s(im,jm),Mg(im,jm)
    !f2py intent (in) kmt,km,jm,im,T,S,rho,rho0,ssh,dzw,depth,k_below,r_above_s,depth_km,d0_in,active
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
       !$omp parallel do private(kmc,k,i,r_above,r_below,z,buoy,p_above,p) shared(active,im,j,kmt,k_below_s,r_above_s,depth_km,rho,rho0,b0,z_s,Mg)
       do i=1,im
          if (active(i,j)) then
             kmc = kmt(i,j)
             k = k_below_s(i,j)
             r_above = r_above_s(i,j)
             r_below = 1.d0 - r_above
             z = r_above*depth(k,i,j) + r_below*depth(k+1,i,j)
             z_s(i,j) = z
             buoy(1:k) = -grav*rrho0*(rho(1:k,i,j) - rho0(1:k,i,j))
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
    REAL*8 :: r_below,r_above,sig8(2),median_sig8(2),T8(2),S8(2)

    do j=1,jm
       !$omp parallel do private(k,i,r_above,r_below,T8,S8,sig8,median_sig8) shared(active,im,j,k_below_s,r_above_s,depth_km,median_depth_km,sig_s,median_sig_s)
       do i=1,im
          if (active(i,j)) then
             k = k_below_s(i,j)
             r_above = r_above_s(i,j)
             r_below = 1.d0 - r_above
             T8 = T(k:k+1,i,j)
             S8 = S(k:k+1,i,j)
             call sigma_n(T8,S8,2,depth_km,sig8)
             call sigma_n(T8,S8,2,median_depth_km,median_sig8)
             sig_s(i,j) = r_above*sig8(1) + r_below*sig8(2)
             median_sig_s(i,j) = r_above*median_sig8(1) + r_below*median_sig8(2)
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


  SUBROUTINE sigma_n( theta,S,n,depth_km,rho)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN) ::   theta(n),S(n)
    REAL*8, INTENT(IN) ::   depth_km
    REAL*8, INTENT(OUT) ::   rho(n)
    !f2py intent (in) n
    !f2py intent (in) theta,S
    !f2py intent (in) depth_km
    !f2py intent (out) rho
    INTEGER ::   i  ! loop counter
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    if (neos==0) then   !==  Jackett and McDougall (1994) formulation  ==!
       do i=1,n
          zt = theta(i)
          zs = S(i)
          zh = depth_km*1000.d0                  ! depth
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
          ! masked in situ density anomaly
          rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
       end do
       ! nemo uses prd = (rho-rau0)/rau0
    else if (neos==1) then  !==  Linear formulation = F( temperature )  ==!
       zs = 1.d0 + 0.0285d0
       ! prd =  0.0285 - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(zs - rn_alpha *theta(i) ) - 1000.d0
       end do
    else if (neos==2) then  !==  Linear formulation = F( temperature , salinity )  ==!
       ! prd =  rn_beta * S - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(1.d0 + rn_beta*S(i) - rn_alpha *theta(i) ) - 1000.d0
       end do
    else
       print *,'unknown eos'
       stop
    end if
  END SUBROUTINE sigma_n

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
end module interp
