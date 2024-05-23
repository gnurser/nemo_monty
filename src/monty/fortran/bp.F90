module bp
  IMPLICIT NONE
  REAL*8 :: grav  = 9.80665d0
  REAL*8 :: sigma0  = 35.d0
  REAL*8 :: rho0  = 1035.d0
  REAL*8 :: fillvalue = 0.d0
  real*8, allocatable :: e3w(:,:,:), dp(:)
  INTEGER*4, allocatable :: kmt(:,:)

contains
  SUBROUTINE setup_bp_mn(fillvalue_in,kmt_in,e3w_in,nz,ny,nx)

    INTEGER*4, INTENT(IN) ::   nz,ny,nx
    INTEGER*4, INTENT(IN) :: kmt_in(nx,ny)
    REAL*4, INTENT(IN) ::   fillvalue_in
    REAL*8, INTENT(IN) ::   e3w_in(nz,ny,nx)
    !f2py intent (in) nz,ny,nx,fillvalue_in,kmt_in,e3w_in

    allocate (kmt(nx,ny),e3w(nz,ny,nx),dp(nz))
    e3w(:,:,:) = e3w_in(:,:,:)
    kmt(:,:) = kmt_in(:,:)

    fillvalue = fillvalue_in
  end subroutine setup_bp_mn


  SUBROUTINE bp_mn4(rho,ssh,bp,nz,ny,nx)
    INTEGER*4, INTENT(IN) ::   nz,ny,nx
    REAL*4, INTENT(IN) ::   rho(nz,ny,nx),ssh(nx,ny)
    REAL*4, INTENT(OUT) ::  bp(nx,ny)
    !f2py intent (in) rho,ssh,nz,ny,nx
    !f2py intent (out) bp

    INTEGER*4 :: i,j,k, km

    bp(:,:) = fillvalue

    do j=1,ny
       do i=1,nx
          km = kmt(i,j)
          if (km < 2) cycle
          dp(1) = grav*(ssh(i,j)*rho0 + .5d0*e3w(1,j,i)*(rho(1,j,i)) - sigma0)
          do k=2,km
             dp(k) = grav*e3w(k,j,i)*(.5d0*(rho(k,j,i)+rho(k-1,j,i)) - sigma0)
             end do
             ! if (i==150) then
             !    print *,'j=',j
             !    print *
             !    print *,'k=',km
             !    print *,dp(1:km)
             !    print *
             !    print *, shape(rho),shape(e3w)
             !    print *, rho(1:km,j,i)
             !    print *, e3w(1:km,j,i)
             !    print *
             !    print *
             ! end if
             bp(i,j) = 1.d-3*sum(dp(1:km))
          end do
       end do

  end subroutine bp_mn4

  SUBROUTINE bp_mn8(rho,ssh,bp,nz,ny,nx)
    INTEGER*4, INTENT(IN) ::   nz,ny,nx
    REAL*8, INTENT(IN) ::   rho(nz,ny,nx),ssh(nx,ny)
    REAL*8, INTENT(OUT) ::  bp(nx,ny)
    !f2py intent (in) rho,ssh,nz,ny,nx
    !f2py intent (out) bp

    INTEGER*4 :: i,j,k, km

    bp(:,:) = fillvalue

    do j=1,ny
       do i=1,nx
          km = kmt(i,j)
          if (km < 2) cycle
          dp(1) = grav*(ssh(i,j)*rho0 + .5d0*e3w(1,j,i)*(rho(1,j,i)) - sigma0)
          do k=2,km
             dp(k) = grav*e3w(k,j,i)*(.5d0*(rho(k,j,i)+rho(k-1,j,i)) - sigma0)
             end do
             ! if (i==150) then
             !    print *,'j=',j
             !    print *
             !    print *,'k=',km
             !    print *,dp(1:km)
             !    print *
             !    print *, shape(rho),shape(e3w)
             !    print *, rho(1:km,j,i)
             !    print *, e3w(1:km,j,i)
             !    print *
             !    print *
             ! end if
             bp(i,j) = 1.d-3*sum(dp(1:km))
          end do
       end do

  end subroutine bp_mn8

  SUBROUTINE rho_mn( fillvalue,mask,theta,S,depth,rho,n)
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
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta(n),S(n),depth(n)
    REAL*4, INTENT(OUT) ::  rho(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) rho

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    do i=1,n
       if (mask(i)) then
          rho(i) = fillvalue
          cycle
       end if

       zt = theta(i)
       zs = S(i)
       zh = depth(i)                 ! depth
       zsr= SQRT( ABS( zs) )        ! square root salinity
       !
       ! compute volumic mass pure water at atm pressure
       zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
            &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594d0
       ! seawater volumic mass atm pressure
       zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
            &                                         -4.0899e-3 ) *zt+0.824493d0
       zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
       zr4= 4.8314e-4
       !
       ! potential volumic mass (reference to the surface)
       zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
       !
       ! add the compression terms
       ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
       zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
       zb = zbw + ze * zs
       !
       zd = -2.042967e-2
       zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
       zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
       za = ( zd*zsr + zc ) *zs + zaw
       !
       zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
       za1= ( (   2.326469e-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
       zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
       zk0= ( zb1*zsr + za1 )*zs + zkw
       !
       ! masked in situ density anomaly
       rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
    end do
    !
  END SUBROUTINE rho_mn
end module bp
