MODULE module_gad_nrel5

REAL, PARAMETER :: pi=3.1415926535897932384626433

CONTAINS

!===============================================================================
!
!  NREL 5 MW Turbine
!
!===============================================================================

!=========================================================================

FUNCTION calc_rpm_nrel5( x ) result( rpm )

!_________________________________________________________________________
!
! Data from NREL 5 MW
! x =  wind speed
!_________________________________________________________________________


  REAL, INTENT( in ) :: x

  INTEGER, PARAMETER :: nws = 50

  REAL, DIMENSION( 0:nws-1 ) :: ws
  REAL, DIMENSION( 0:nws-1 ) :: rp

  REAL :: rpm
  INTEGER :: i

  rpm = 0.0

  rp(0) = 0.05
  rp(1) = 0.31
  rp(2) = 3.29
  rp(3) = 5.23
  rp(4) = 6.57
  rp(5) = 6.97
  rp(6) = 7.06
  rp(7) = 7.18
  rp(8) = 7.32
  rp(9) = 7.50
  rp(10) = 7.70
  rp(11) = 7.93
  rp(12) = 8.18
  rp(13) = 8.46
  rp(14) = 8.76
  rp(15) = 9.15
  rp(16) = 9.72
  rp(17) = 10.29
  rp(18) = 10.86
  rp(19) = 11.43
  rp(20) = 11.77
  rp(21) = 11.89
  rp(22) = 12.10
  rp(23) = 12.10
  rp(24) = 12.10
  rp(25) = 12.10
  rp(26) = 12.10
  rp(27) = 12.10
  rp(28) = 12.10
  rp(29) = 12.10
  rp(30) = 12.10
  rp(31) = 12.10
  rp(32) = 12.10
  rp(33) = 12.10
  rp(34) = 12.10
  rp(35) = 12.10
  rp(36) = 12.10
  rp(37) = 12.10
  rp(38) = 12.10
  rp(39) = 12.10
  rp(40) = 12.10
  rp(41) = 12.10
  rp(42) = 12.10
  rp(43) = 12.10
  rp(44) = 12.10
  rp(45) = 12.10
  rp(46) = 12.10
  rp(47) = 12.10
  rp(48) = 12.10
  rp(49) = 12.10

  ws(0) = 0.5
  ws(1) = 1.0
  ws(2) = 1.5
  ws(3) = 2.0
  ws(4) = 2.5
  ws(5) = 3.0
  ws(6) = 3.5
  ws(7) = 4.0
  ws(8) = 4.5
  ws(9) = 5.0
  ws(10) = 5.5
  ws(11) = 6.0
  ws(12) = 6.5
  ws(13) = 7.0
  ws(14) = 7.5
  ws(15) = 8.0
  ws(16) = 8.5
  ws(17) = 9.0
  ws(18) = 9.5
  ws(19) = 10.0
  ws(20) = 10.5
  ws(21) = 11.0
  ws(22) = 11.5
  ws(23) = 12.0
  ws(24) = 12.5
  ws(25) = 13.0
  ws(26) = 13.5
  ws(27) = 14.0
  ws(28) = 14.5
  ws(29) = 15.0
  ws(30) = 15.5
  ws(31) = 16.0
  ws(32) = 16.5
  ws(33) = 17.0
  ws(34) = 17.5
  ws(35) = 18.0
  ws(36) = 18.5
  ws(37) = 19.0
  ws(38) = 19.5
  ws(39) = 20.0
  ws(40) = 20.5
  ws(41) = 21.0
  ws(42) = 21.5
  ws(43) = 22.0
  ws(44) = 22.5
  ws(45) = 23.0
  ws(46) = 23.5
  ws(47) = 24.0
  ws(48) = 24.5
  ws(49) = 25.0

  DO i = 1,nws-1
     
     IF ( ( x .GE. ws(i-1) ) .AND. ( x .LT. ws(i) ) ) THEN
        rpm = rp(i-1) + (x - ws(i-1))*( ( rp(i) - rp(i-1) )/( ws(i) - ws(i-1) ) )
        GOTO 1
     ENDIF

  END DO
  
  IF ( x .LT. ws(0) ) rpm = rp(0)
  IF ( x .EQ. ws(nws-1) ) rpm = rp(nws-1)

1 CONTINUE
 
!PRINT*,x,rpm
!_________________________________________________________________________

END FUNCTION calc_rpm_nrel5
!_________________________________________________________________________
 
!=========================================================================

FUNCTION calc_pitch_nrel5( x ) RESULT( pitch )

!_________________________________________________________________________
!
! Data from NREL 5 MW
! x = wind speeds
!_________________________________________________________________________
 
  REAL, INTENT( IN ) :: x

  INTEGER, PARAMETER :: nws = 50

  REAL, DIMENSION( 0:nws-1 ) :: ws 
  REAL, DIMENSION( 0:nws-1 ) :: p

  REAL :: pitch
  INTEGER :: i

  pitch = 0.0

  p(0) = 0.00
  p(1) = 0.00
  p(2) = 0.00
  p(3) = 0.00
  p(4) = 0.00
  p(5) = 0.00
  p(6) = 0.00
  p(7) = 0.00
  p(8) = 0.00
  p(9) = 0.00
  p(10) = 0.00
  p(11) = 0.00
  p(12) = 0.00
  p(13) = 0.00
  p(14) = 0.00
  p(15) = 0.00
  p(16) = 0.00
  p(17) = 0.00
  p(18) = 0.00
  p(19) = 0.00
  p(20) = 0.00
  p(21) = 0.00
  p(22) = 1.49
  p(23) = 3.71
  p(24) = 5.21
  p(25) = 6.45
  p(26) = 7.54
  p(27) = 8.53
  p(28) = 9.46
  p(29) = 10.32
  p(30) = 11.15
  p(31) = 11.94
  p(32) = 12.7
  p(33) = 13.44
  p(34) = 14.13
  p(35) = 14.81
  p(36) = 15.46
  p(37) = 16.1
  p(38) = 16.73
  p(39) = 17.34
  p(40) = 17.95
  p(41) = 18.54
  p(42) = 19.13
  p(43) = 19.71
  p(44) = 20.28
  p(45) = 20.84
  p(46) = 21.39
  p(47) = 21.93
  p(48) = 22.47
  p(49) = 22.99
  
  ws(0) = 0.5
  ws(1) = 1.0
  ws(2) = 1.5
  ws(3) = 2.0
  ws(4) = 2.5
  ws(5) = 3.0
  ws(6) = 3.5
  ws(7) = 4.0
  ws(8) = 4.5
  ws(9) = 5.0
  ws(10) = 5.5
  ws(11) = 6.0
  ws(12) = 6.5
  ws(13) = 7.0
  ws(14) = 7.5
  ws(15) = 8.0
  ws(16) = 8.5
  ws(17) = 9.0
  ws(18) = 9.5
  ws(19) = 10.0
  ws(20) = 10.5
  ws(21) = 11.0
  ws(22) = 11.5
  ws(23) = 12.0
  ws(24) = 12.5
  ws(25) = 13.0
  ws(26) = 13.5
  ws(27) = 14.0
  ws(28) = 14.5
  ws(29) = 15.0
  ws(30) = 15.5
  ws(31) = 16.0
  ws(32) = 16.5
  ws(33) = 17.0
  ws(34) = 17.5
  ws(35) = 18.0
  ws(36) = 18.5
  ws(37) = 19.0
  ws(38) = 19.5
  ws(39) = 20.0
  ws(40) = 20.5
  ws(41) = 21.0
  ws(42) = 21.5
  ws(43) = 22.0
  ws(44) = 22.5
  ws(45) = 23.0
  ws(46) = 23.5
  ws(47) = 24.0
  ws(48) = 24.5
  ws(49) = 25.0

  DO i = 1,nws-1

     IF ( ( x .GT. ws(i-1) ) .AND. ( x .LE. ws(i) ) ) THEN 
        pitch = p(i-1) + (x - ws(i-1))*( ( p(i) - p(i-1) )/( ws(i) - ws(i-1) ) )
        GOTO 2
     ENDIF
     
  ENDDO

  IF ( x .LE. ws(0) ) pitch = p(0)
  IF ( x .GT. ws(nws-1) ) pitch = p(nws-1)

2 CONTINUE
  
!PRINT*,x,pitch
!_________________________________________________________________________

END FUNCTION calc_pitch_nrel5

!=========================================================================

FUNCTION calc_fnft_nrel5( v0, rho, blades, bladelen, diameter,           &
                         radius, pitchdeg, rpm, dscale, cscale)         &
                         result( fnft )
!_________________________________________________________________________
!
!  --v0       - "free stream" velocity normal to the turbine rotor plane, [m/s]
!  --rho      - air density, [kg/m^3]
!  --blades   - number of turbine blades, [non-dimensional]
!  --diameter - turbine rotor diameter, [m]
!  --xcoord   - streamwise coordinate of the grid point with respect to
!               local coordinate system with the origin at the turbine hub, [m]
!  --ycoord   - cross-stream coordinate of the grid point with respect to
!               local coordinate system with the origin at the turbine hub, [m]
!  --zcoord   - vertical coordinate of the grid point with respect to
!               local coordinate system with the origin at the turbine hub, [m]
!
!  --fnft     - one-dimensional array (two members) with FUNCTION output
!               specific normal and tangential forces per unit area, [N/m^2]
!
     real, intent( in ) :: v0, rho, blades, bladelen, diameter
     real, intent( in ) :: dscale, cscale 
     real, intent( in ) :: radius
     real, intent( in ) :: rpm
     real, intent( in ) :: pitchdeg
     real, dimension( 2 )  :: fnft
     real               :: rps, omega 
     real               :: radmax, radmin
     real               :: an, at
     real               :: anold, atold
     real               :: cord
     real               :: twist
     real               :: xi, xideg, beta, betadeg, psi, psideg
     real               :: vn, vt
     real               :: cl, cd
     real               :: cn, ct
     real               :: solid
     real               :: fac1, fac2, tiploss, hubloss, lossfac
     real               :: vr
     real               :: lift, drag
     real               :: fn, ft
     real               :: epsi, epscrit
     real               :: cd_nac
     integer            :: iter, itmax
!
     rps   = rpm/60.0
     omega = 2.0*pi*rps
!
     radmax = 0.5*diameter
     radmin = radmax - bladelen
!

  fnft(1) = 0.0
  fnft(2) = 0.0

  IF (radius .LT. radmax) THEN !JDM
     
     IF (radius .GT. radmin) THEN !JDM ! MSG CHANGED FROM .GE. TO .GT.

        an = 0.0
        at = 0.0
!
        iter = 0
        epsi = 1.0        
        epscrit = 1.e-5
        itmax   = 20
!
!    --Incidence angle

        xideg = calc_twist_nrel5( radius/dscale ) + pitchdeg
        xi    = xideg*pi/180.0
!
!    --Cord length 

        cord  = cscale * calc_chord_nrel5( radius/dscale )
!
!    --Compute induction factors an and at iteratively

        do while ( ( epsi .gt. epscrit ) .and. ( iter .lt. itmax ) ) 
!
           iter = iter + 1
!
!       --Normal and tangential velocity components

           an = 0.02 !JDM This is just the tip loss component of an. This was computed using the
                     !    IDL code as the difference between an computed with and without tip loss
                     !    across the power curve. The value of 0.02 was nearly constant.

           vn = v0*( 1.0 - an )
           vt = omega*radius*( 1.0 + at ) 
 
!       --Advance angle 
  
           psi    = atan2(vn,vt)
           psideg = psi*180.0/pi
! print*,'iter',iter
! print*,'v0,  omega, radius',v0,omega,radius
! print*,'vn, at, vt',vn,at,vt
! print*,'psi',psi
! print*,'------'

!       --Incidence angle  
 
           beta    = psi - xi
           betadeg = beta*180.0/pi

!       --Lift and drag coefficients  
     
           cl = calc_cl_nrel5( radius/dscale, betadeg )
           cd = calc_cd_nrel5( radius/dscale, betadeg )
  
!       --Normal and tangential force coefficients     
  
           cn = cl*cos( psi ) + cd*sin( psi )
           ct = cl*sin( psi ) - cd*cos( psi )

!       --Solidity factor

           solid = cord*blades / ( 2.0*pi*radius )

!       --Tip and hub loss

!          When psi < 0, fac1 and fac 2 can give bizarre answers, and resulting in tip and hub loss terms - NaN.
!          Some logic should be added to keep things from going crazy
           
!          RSA update: the argument of acos should have an abs. value <= 1. When psi
!          is negative, fac1 and/or fac2 are negative, causing exp(-fac) to be a
!          large positive number, causing acos to return NaN. When psi is positive, fac1 and fac2 are
!          positive, such that exp(-fac) is small and positive, which is ok.
!          Thus, we have added logic around the calls to calc_fnft above such
!          that when v0 < 0 (ie psi < 0), this function is not called and the
!          force is set to 0. The turbine should still yaw into place eventually.

           if ( abs( psi ) .gt. 0.001 ) then
             if ( psi .lt. (pi-0.01)) then  !! if STATEMENT ADDED BY MSG TO AVOID GETTING NAN IN tiploss AND hubloss
              fac1 = blades / 2.0*( radmax - radius ) / ( radius*sin( psi ) )
              fac2 = blades / 2.0*( radius - radmin ) / ( radmin*sin( psi ) )

              tiploss = (2.0/pi)*acos( exp( -fac1 ) )
              hubloss = (2.0/pi)*acos( exp( -fac2 ) )

              lossfac = tiploss*hubloss

!             --Induction factors

              anold = an
              atold = at

!JDM             an = ( 1. + 4.0*lossfac*sin( psi )*sin( psi )/ ( solid*cn ) )**(-1.0)

              at = ( 4.0*lossfac*sin( psi )*cos( psi )/ ( solid*ct ) - 1.0 )**(-1.0)

! print*,'inside of iteration'
! print*,'iter',iter
! print*,'omega, radius, radmin, radmax',omega,radius,radmin,radmax
! print*,'v0, vn, vt, an, at',v0,vn,vt,an,at
! print*,'psi, fac1, fac2',psi,fac1,fac2
! print*,'tiploss, hubloss, lossfac',tiploss,hubloss,lossfac
! print*,'------'

            else

             anold = an
             atold = at

             an = 0.0
             at = 0.0


            endif
          else

             anold = an
             atold = at

             an = 0.0
             at = 0.0


          endif
!       --Convergence criteria

           epsi = sqrt( ( anold - an )**2.0 + ( atold - at )**2.0 )
  
        end do
        if (isnan(at)) then
                print*,'at is NaN lossfac=',lossfac,'psi=',psi,'solid=',solid,'ct=',ct, & 
                        'tiploss=',tiploss,'hubloss=',hubloss,'fac1=',fac1,'fac2=',fac2,'radius=',radius
        endif
!        if (at .GT. 0) then
!                print*,'at is > 0 at=',at,'psi=',psi,'fac1=',fac1,'fac2=',fac2,'radius=',radius
!        endif
        if (isnan(at)) then
                print*,'ADD LOGIC SO THAT PSI > 0 deg AND PSI < 180 deg (add error margins to each limit) '
        endif
!print*,'radius,radmax,radmin',radius,radmax,radmin
!print*,'radius - radmin',radius - radmin
!print*,'fac1,fac2',fac1,fac2
!print*,'psi,sin(psi)',psi,sin(psi)
!print*,'tiploss,hubloss',tiploss,hubloss
!print*,'cl,cd',cl,cd
!print*,'cn,ct',cn,ct

!    --Normal and tangential velocity components

       vn = v0*( 1.0 - an )

       vt = omega*radius*( 1.0 + at ) 

!    --Relative velocity

       vr = sqrt( vn*vn + vt*vt ) 
! 
!    --Advance angle   

       psi    = atan2(vn,vt)
       psideg = psi*180.0/pi
        if (isnan(psi)) then
                print*,'vn=',vn,'vt=',vt,'v0=',v0,'an=',an,'omega=',omega,'radius=',radius,'at=',at
        endif
!
!    --Incidence angle   

       beta    = psi - xi
       betadeg = beta*180.0/pi

! print*,'outside of iteration'
! print*,'v0, vn, at, vt',v0,vn,at,vt
! print*,'psideg, betadeg',psideg,betadeg
! print*,'------'

!    --Lift and drag coefficients     
  
       cl = calc_cl_nrel5( radius/dscale, betadeg )
       cd = calc_cd_nrel5( radius/dscale, betadeg )

!    --Lift and drag forces

       lift = 0.5*rho*vr*vr*cord*cl 
       drag = 0.5*rho*vr*vr*cord*cd 

!    --Specific normal and tangential forces

       fn = lift*cos( psi ) + drag*sin( psi )
       ft = lift*sin( psi ) - drag*cos( psi )

!       print*,'psi',psideg

       fn = blades*fn/(2.0*pi*radius)
       ft = blades*ft/(2.0*pi*radius)

       if( fn .lt. 0.0 ) fn = 0.0
       if( ft .lt. 0.0 ) ft = 0.0
        if (isnan(fn)) then
                print*,'fn is NaN, psi=',psi,'lift=',lift,'drag=',drag
        endif
        if (isnan(ft)) then
                print*,'ft is NaN, psi=',psi,'lift=',lift,'drag=',drag
        endif
       fnft(1) = fn
       fnft(2) = ft              

!print*,'radius,radmax,radmin',radius,radmax,radmin
!print*,'radius - radmin',radius - radmin
!print*,'fac1,fac2',fac1,fac2
!print*,'psi,sin(psi)',psi,sin(psi)
!print*,'tiploss,hubloss',tiploss,hubloss
!print*,'cl,cd',cl,cd
!print*,'cn,ct',cn,ct
!print*,'fn,ft',fn,ft


     ENDIF !radius .GE. radmin

  ENDIF   !radius .LT. radmax

END FUNCTION calc_fnft_nrel5

!=========================================================================

FUNCTION calc_twist_nrel5( x ) RESULT( twist )
!_________________________________________________________________________
!
! Computes blade twist as a function of radial distance from hub center
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: x

  INTEGER, PARAMETER :: ns = 19

  REAL, DIMENSION( 0:ns-1 ) :: r
  REAL, DIMENSION( 0:ns-1 ) :: t

  REAL :: twist 
  INTEGER :: i

  twist = 0.0
  
  r(0) = 1.5
  r(1) = 2.8667
  r(2) = 5.6
  r(3) = 8.3333
  r(4) = 11.75
  r(5) = 15.85
  r(6) = 19.95
  r(7) = 24.05
  r(8) = 28.15
  r(9) = 32.25
  r(10) = 36.35
  r(11) = 40.45
  r(12) = 44.55
  r(13) = 48.65
  r(14) = 52.75
  r(15) = 56.1667
  r(16) = 58.9
  r(17) = 61.6333
  r(18) = 63.000
 
  t(0) = 13.308
  t(1) = 13.308
  t(2) = 13.308
  t(3) = 13.308
  t(4) = 13.308
  t(5) = 11.48
  t(6) = 10.162
  t(7) = 9.011
  t(8) = 7.795
  t(9) = 6.544
  t(10) = 5.361
  t(11) = 4.188
  t(12) = 3.125
  t(13) = 2.319
  t(14) = 1.526
  t(15) = 0.863
  t(16) = 0.37
  t(17) = 0.106
  t(18) = 0.106
  
  DO i = 1,ns-1
 
     IF ( ( x .GT. r(i-1) ) .AND. ( x .LE. r(i) ) ) THEN
        twist = t(i-1) + (x - r(i-1))*( ( t(i) - t(i-1) )/( r(i) - r(i-1) ) )
        goto 1
     ENDIF

1 CONTINUE

     IF ( x .LT. r(0) ) THEN
        PRINT*,'—————————————'
        PRINT*,'radius < blade root in calc_twist_nrel5.pro'
	PRINT*, 'radius: ', x
        PRINT*,'—————————————'
        STOP
     ENDIF
     
     IF ( x .GT. r(ns-1) ) THEN
        PRINT*,'—————————————'
        PRINT*,'radius >  blade length in calc_twist_nrel5.pro'
	PRINT*,'radius: ', x
        PRINT*,'—————————————'
        STOP
     ENDIF 
     
  ENDDO
  
!PRINT*,x,twist
  
!_________________________________________________________________________
  
END FUNCTION calc_twist_nrel5
!_________________________________________________________________________


!=========================================================================

FUNCTION calc_chord_nrel5( x ) RESULT( chord )
!_________________________________________________________________________
!
! Data from NREL 5 MW
! x = radial distance from hub center
!_________________________________________________________________________
 

  REAL, INTENT( IN ) :: x

  INTEGER, PARAMETER :: ns = 19

  REAL, DIMENSION( 0:ns-1 ) :: r
  REAL, DIMENSION( 0:ns-1 ) :: c

  REAL :: chord 
  INTEGER :: i

  chord = 0.0
  
  r(0) = 1.5
  r(1) = 2.8667
  r(2) = 5.6
  r(3) = 8.3333
  r(4) = 11.75
  r(5) = 15.85
  r(6) = 19.95
  r(7) = 24.05
  r(8) = 28.15
  r(9) = 32.25
  r(10) = 36.35
  r(11) = 40.45
  r(12) = 44.55
  r(13) = 48.65
  r(14) = 52.75
  r(15) = 56.1667
  r(16) = 58.9
  r(17) = 61.6333
  r(18) = 63.000

  c(0) = 3.542
  c(1) = 3.542
  c(2) = 3.854
  c(3) = 4.167
  c(4) = 4.557
  c(5) = 4.652
  c(6) = 4.458
  c(7) = 4.249
  c(8) = 4.007
  c(9) = 3.748
  c(10) = 3.502
  c(11) = 3.256
  c(12) = 3.01
  c(13) = 2.764
  c(14) = 2.518
  c(15) = 2.313
  c(16) = 2.086
  c(17) = 1.419
  c(18) = 1.419

  DO i = 1,ns-1
     IF ( ( x .GE. r(i-1) ) .AND. ( x .LT. r(i) ) ) THEN
        chord = c(i-1) + (x - r(i-1))*( ( c(i) - c(i-1) )/( r(i) - r(i-1) ) )
        GOTO 1
     ENDIF
  ENDDO

  IF ( x .LT. r(0) ) chord = 0.0 !No airfoil here, just a cylinder
  IF ( x .EQ. r(ns-1) ) chord = c(ns-1)
  
  IF ( x .GT. r(ns-1) ) THEN 
     PRINT*,'—————————————'
     PRINT*,'radius > blade length in calc_chord_nrel5.pro'
     PRINT*,'radius: ', x
     PRINT*,'—————————————'
     STOP
  ENDIF
  
1 CONTINUE
  
!PRINT*,x,chord
   
!_________________________________________________________________________
  
END FUNCTION calc_chord_nrel5
!_________________________________________________________________________


!=========================================================================

FUNCTION calc_cl_nrel5( x,y ) RESULT( cl )
!_________________________________________________________________________
!
! Calls appropriate airfoil section table depending upon radial distance
! from blade root. Data from NREL 5 MW
!
! calc_cl1_nrel5 uses data from Cylinder1.dat
! calc_cl2_nrel5 uses data from Cylinder2.dat
! calc_cl3_nrel5 uses data from DU40_A17.dat
! calc_cl4_nrel5 uses data from DU35_A17.dat
! calc_cl5_nrel5 uses data from DU30_A17.dat
! calc_cl6_nrel5 uses data from DU25_A17.dat
! calc_cl7_nrel5 uses data from DU21_A17.dat
! calc_cl8_nrel5 uses data from NACA64_A17.dat
!
! x = radial distance from hub center
! y = angle of attack
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: x,y

  INTEGER, PARAMETER :: ns = 9

  REAL, DIMENSION( 0:ns-1 ) :: r

  REAL :: cl

  INTEGER :: i

  cl = 0.0
  
  r(0) = 1.5
  r(1) = 8.3333
  r(2) = 11.75
  r(3) = 15.85
  r(4) = 24.05
  r(5) = 28.15
  r(6) = 36.35
  r(7) = 44.55
  r(8) = 63.000
  
  IF ( ( x .GT. r( 0) ) .AND. ( x .LE. r( 1) ) ) THEN
       cl = calc_cl1_nrel5(y)
     ELSEIF ( ( x .GT. r( 1) ) .AND. ( x .LE. r( 2) ) ) THEN
       cl = calc_cl2_nrel5(y)
     ELSEIF ( ( x .GT. r( 2) ) .AND. ( x .LE. r( 3) ) ) THEN
       cl = calc_cl3_nrel5(y)
     ELSEIF ( ( x .GT. r( 3) ) .AND. ( x .LE. r( 4) ) ) THEN
       cl = calc_cl4_nrel5(y)
     ELSEIF ( ( x .GT. r( 4) ) .AND. ( x .LE. r( 5) ) ) THEN
       cl = calc_cl5_nrel5(y)
     ELSEIF ( ( x .GT. r( 5) ) .AND. ( x .LE. r( 6) ) ) THEN
       cl = calc_cl6_nrel5(y)
     ELSEIF ( ( x .GT. r( 6) ) .AND. ( x .LE. r( 7) ) ) THEN
       cl = calc_cl7_nrel5(y)
     ELSEIF ( ( x .GT. r( 7) ) .AND. ( x .LE. r( 8) ) ) THEN
       cl = calc_cl8_nrel5(y)  
  ENDIF
  
  IF ( x .LT. r(0) ) THEN 
     PRINT*,'—————————————'
     PRINT*,'radius less than blade root in calc_cl_nrel5.pro'
     PRINT*,'radius: ', x
     PRINT*,'—————————————'
     STOP
  ENDIF
  IF ( x .GT. r(ns-1) ) THEN 
     PRINT*,'—————————————'
     PRINT*,'radius greater than blade length in calc_cl_nrel5.pro'
     PRINT*,'radius: ', x
     PRINT*,'—————————————'
     STOP
  ENDIF
  
!_________________________________________________________________________
  
END FUNCTION calc_cl_nrel5
!_________________________________________________________________________


FUNCTION calc_cl1_nrel5( y ) RESULT( cl1 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: Cylinder1
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  REAL :: cl1

  cl1 = 0.0001
!_________________________________________________________________________
  
END FUNCTION calc_cl1_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl2_nrel5( y ) RESULT( cl2 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: Cylinder2
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  REAL :: cl2

  cl2 = 0.0001
!_________________________________________________________________________
  
END FUNCTION calc_cl2_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl3_nrel5( y ) RESULT( cl3 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU40_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 136

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl3
  INTEGER :: i

  cl3 = 0.0

!                alpha    CL   
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.2180]
  cl(0:1, 2  ) = [-170.0, 0.3970]
  cl(0:1, 3  ) = [-160.0, 0.6420]
  cl(0:1, 4  ) = [-155.0, 0.7150]
  cl(0:1, 5  ) = [-150.0, 0.7570]
  cl(0:1, 6  ) = [-145.0, 0.7720]
  cl(0:1, 7  ) = [-140.0, 0.7620]
  cl(0:1, 8  ) = [-135.0, 0.7310]
  cl(0:1, 9  ) = [-130.0, 0.6800]
  cl(0:1, 10 ) = [-125.0, 0.6130]
  cl(0:1, 11 ) = [-120.0, 0.5320]
  cl(0:1, 12 ) = [-115.0, 0.4390]
  cl(0:1, 13 ) = [-110.0, 0.3370]
  cl(0:1, 14 ) = [-105.0, 0.2280]
  cl(0:1, 15 ) = [-100.0, 0.1140]
  cl(0:1, 16 ) = [-95.0, -0.0020]
  cl(0:1, 17 ) = [-90.0, -0.1200]
  cl(0:1, 18 ) = [-85.0, -0.2360]
  cl(0:1, 19 ) = [-80.0, -0.3490]
  cl(0:1, 20 ) = [-75.0, -0.4560]
  cl(0:1, 21 ) = [-70.0, -0.5570]
  cl(0:1, 22 ) = [-65.0, -0.6470]
  cl(0:1, 23 ) = [-60.0, -0.7270]
  cl(0:1, 24 ) = [-55.0, -0.7920]
  cl(0:1, 25 ) = [-50.0, -0.8420]
  cl(0:1, 26 ) = [-45.0, -0.8740]
  cl(0:1, 27 ) = [-40.0, -0.8860]
  cl(0:1, 28 ) = [-35.0, -0.8750]
  cl(0:1, 29 ) = [-30.0, -0.8390]
  cl(0:1, 30 ) = [-25.0, -0.7770]
  cl(0:1, 31 ) = [-24.0, -0.7610]
  cl(0:1, 32 ) = [-23.0, -0.7440]
  cl(0:1, 33 ) = [-22.0, -0.7250]
  cl(0:1, 34 ) = [-21.0, -0.7060]
  cl(0:1, 35 ) = [-20.0, -0.6850]
  cl(0:1, 36 ) = [-19.0, -0.6620]
  cl(0:1, 37 ) = [-18.0, -0.6350]
  cl(0:1, 38 ) = [-17.0, -0.6050]
  cl(0:1, 39 ) = [-16.0, -0.5710]
  cl(0:1, 40 ) = [-15.0, -0.5340]
  cl(0:1, 41 ) = [-14.0, -0.4940]
  cl(0:1, 42 ) = [-13.0, -0.4520]
  cl(0:1, 43 ) = [-12.0, -0.4070]
  cl(0:1, 44 ) = [-11.0, -0.3600]
  cl(0:1, 45 ) = [-10.0, -0.3110]
  cl(0:1, 46 ) = [-8.0, -0.2080]
  cl(0:1, 47 ) = [-6.0, -0.1110]
  cl(0:1, 48 ) = [-5.5, -0.0900]
  cl(0:1, 49 ) = [-5.0, -0.0720]
  cl(0:1, 50 ) = [-4.5, -0.0650]
  cl(0:1, 51 ) = [-4.0, -0.0540]
  cl(0:1, 52 ) = [-3.5, -0.0170]
  cl(0:1, 53 ) = [-3.0, 0.0030]
  cl(0:1, 54 ) = [-2.5, 0.0140]
  cl(0:1, 55 ) = [-2.0, 0.0090]
  cl(0:1, 56 ) = [-1.5, 0.0040]
  cl(0:1, 57 ) = [-1.0, 0.0360]
  cl(0:1, 58 ) = [-0.5, 0.0730]
  cl(0:1, 59 ) = [0.0, 0.1370]
  cl(0:1, 60 ) = [0.5, 0.2130]
  cl(0:1, 61 ) = [1.0, 0.2920]
  cl(0:1, 62 ) = [1.5, 0.3690]
  cl(0:1, 63 ) = [2.0, 0.4440]
  cl(0:1, 64 ) = [2.5, 0.5140]
  cl(0:1, 65 ) = [3.0, 0.5800]
  cl(0:1, 66 ) = [3.5, 0.6450]
  cl(0:1, 67 ) = [4.0, 0.7100]
  cl(0:1, 68 ) = [4.5, 0.7760]
  cl(0:1, 69 ) = [5.0, 0.8410]
  cl(0:1, 70 ) = [5.5, 0.9040]
  cl(0:1, 71 ) = [6.0, 0.9670]
  cl(0:1, 72 ) = [6.5, 1.0270]
  cl(0:1, 73 ) = [7.0, 1.0840]
  cl(0:1, 74 ) = [7.5, 1.1400]
  cl(0:1, 75 ) = [8.0, 1.1930]
  cl(0:1, 76 ) = [8.5, 1.2420]
  cl(0:1, 77 ) = [9.0, 1.2870]
  cl(0:1, 78 ) = [9.5, 1.3330]
  cl(0:1, 79 ) = [10.0, 1.3680]
  cl(0:1, 80 ) = [10.5, 1.4000]
  cl(0:1, 81 ) = [11.0, 1.4250]
  cl(0:1, 82 ) = [11.5, 1.4490]
  cl(0:1, 83 ) = [12.0, 1.4730]
  cl(0:1, 84 ) = [12.5, 1.4940]
  cl(0:1, 85 ) = [13.0, 1.5130]
  cl(0:1, 86 ) = [13.5, 1.5380]
  cl(0:1, 87 ) = [14.5, 1.5870]
  cl(0:1, 88 ) = [15.0, 1.6140]
  cl(0:1, 89 ) = [15.5, 1.6310]
  cl(0:1, 90 ) = [16.0, 1.6490]
  cl(0:1, 91 ) = [16.5, 1.6660]
  cl(0:1, 92 ) = [17.0, 1.6810]
  cl(0:1, 93 ) = [17.5, 1.6990]
  cl(0:1, 94 ) = [18.0, 1.7190]
  cl(0:1, 95 ) = [19.0, 1.7510]
  cl(0:1, 96 ) = [19.5, 1.7670]
  cl(0:1, 97 ) = [20.5, 1.7980]
  cl(0:1, 98 ) = [21.0, 1.8100]
  cl(0:1, 99 ) = [22.0, 1.8300]
  cl(0:1, 100) = [23.0, 1.8470]
  cl(0:1, 101) = [24.0, 1.8610]
  cl(0:1, 102) = [25.0, 1.8720]
  cl(0:1, 103) = [26.0, 1.8810]
  cl(0:1, 104) = [28.0, 1.8940]
  cl(0:1, 105) = [30.0, 1.9040]
  cl(0:1, 106) = [32.0, 1.9150]
  cl(0:1, 107) = [35.0, 1.9290]
  cl(0:1, 108) = [40.0, 1.9030]
  cl(0:1, 109) = [45.0, 1.8200]
  cl(0:1, 110) = [50.0, 1.6900]
  cl(0:1, 111) = [55.0, 1.5220]
  cl(0:1, 112) = [60.0, 1.3230]
  cl(0:1, 113) = [65.0, 1.1060]
  cl(0:1, 114) = [70.0, 0.8800]
  cl(0:1, 115) = [75.0, 0.6580]
  cl(0:1, 116) = [80.0, 0.4490]
  cl(0:1, 117) = [85.0, 0.2670]
  cl(0:1, 118) = [90.0, 0.1240]
  cl(0:1, 119) = [95.0, 0.0020]
  cl(0:1, 120) = [100.0, -0.1180]
  cl(0:1, 121) = [105.0, -0.2350]
  cl(0:1, 122) = [110.0, -0.3480]
  cl(0:1, 123) = [115.0, -0.4530]
  cl(0:1, 124) = [120.0, -0.5490]
  cl(0:1, 125) = [125.0, -0.6330]
  cl(0:1, 126) = [130.0, -0.7020]
  cl(0:1, 127) = [135.0, -0.7540]
  cl(0:1, 128) = [140.0, -0.7870]
  cl(0:1, 129) = [145.0, -0.7970]
  cl(0:1, 130) = [150.0, -0.7820]
  cl(0:1, 131) = [155.0, -0.7390]
  cl(0:1, 132) = [160.0, -0.6640]
  cl(0:1, 133) = [170.0, -0.4100]
  cl(0:1, 134) = [175.0, -0.2260]
  cl(0:1, 135) = [180.0, 0.0000]
  
  IF ( y .EQ. cl(0,0) ) cl3 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl3 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl3 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl3_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl3 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl3_nrel5.pro'
  ENDIF
  
!print,y,cl3
  
!_________________________________________________________________________
  
END FUNCTION calc_cl3_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl4_nrel5( y ) RESULT( cl4 )
!_________________________________________________________________________
!
! Data from NREL 5 DU35_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 135

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl4
  INTEGER :: i

  cl4 = 0.0

!                alpha    CL   
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.2230]
  cl(0:1, 2  ) = [-170.0, 0.4050]
  cl(0:1, 3  ) = [-160.0, 0.6580]
  cl(0:1, 4  ) = [-155.0, 0.7330]
  cl(0:1, 5  ) = [-150.0, 0.7780]
  cl(0:1, 6  ) = [-145.0, 0.7950]
  cl(0:1, 7  ) = [-140.0, 0.7870]
  cl(0:1, 8  ) = [-135.0, 0.7570]
  cl(0:1, 9  ) = [-130.0, 0.7080]
  cl(0:1, 10 ) = [-125.0, 0.6410]
  cl(0:1, 11 ) = [-120.0, 0.5600]
  cl(0:1, 12 ) = [-115.0, 0.4670]
  cl(0:1, 13 ) = [-110.0, 0.3650]
  cl(0:1, 14 ) = [-105.0, 0.2550]
  cl(0:1, 15 ) = [-100.0, 0.1390]
  cl(0:1, 16 ) = [-95.0, 0.0210]
  cl(0:1, 17 ) = [-90.0, -0.0980]
  cl(0:1, 18 ) = [-85.0, -0.2160]
  cl(0:1, 19 ) = [-80.0, -0.3310]
  cl(0:1, 20 ) = [-75.0, -0.4410]
  cl(0:1, 21 ) = [-70.0, -0.5440]
  cl(0:1, 22 ) = [-65.0, -0.6380]
  cl(0:1, 23 ) = [-60.0, -0.7200]
  cl(0:1, 24 ) = [-55.0, -0.7880]
  cl(0:1, 25 ) = [-50.0, -0.8400]
  cl(0:1, 26 ) = [-45.0, -0.8750]
  cl(0:1, 27 ) = [-40.0, -0.8890]
  cl(0:1, 28 ) = [-35.0, -0.8800]
  cl(0:1, 29 ) = [-30.0, -0.8460]
  cl(0:1, 30 ) = [-25.0, -0.7840]
  cl(0:1, 31 ) = [-24.0, -0.7680]
  cl(0:1, 32 ) = [-23.0, -0.7510]
  cl(0:1, 33 ) = [-22.0, -0.7330]
  cl(0:1, 34 ) = [-21.0, -0.7140]
  cl(0:1, 35 ) = [-20.0, -0.6930]
  cl(0:1, 36 ) = [-19.0, -0.6710]
  cl(0:1, 37 ) = [-18.0, -0.6480]
  cl(0:1, 38 ) = [-17.0, -0.6240]
  cl(0:1, 39 ) = [-16.0, -0.6010]
  cl(0:1, 40 ) = [-15.0, -0.5790]
  cl(0:1, 41 ) = [-14.0, -0.5590]
  cl(0:1, 42 ) = [-13.0, -0.5390]
  cl(0:1, 43 ) = [-12.0, -0.5190]
  cl(0:1, 44 ) = [-11.0, -0.4990]
  cl(0:1, 45 ) = [-10.0, -0.4800]
  cl(0:1, 46 ) = [-5.5, -0.3850]
  cl(0:1, 47 ) = [-5.0, -0.3590]
  cl(0:1, 48 ) = [-4.5, -0.3600]
  cl(0:1, 49 ) = [-4.0, -0.3550]
  cl(0:1, 50 ) = [-3.5, -0.3070]
  cl(0:1, 51 ) = [-3.0, -0.2460]
  cl(0:1, 52 ) = [-3.0, -0.2400]
  cl(0:1, 53 ) = [-2.5, -0.1630]
  cl(0:1, 54 ) = [-2.0, -0.0910]
  cl(0:1, 55 ) = [-1.5, -0.0190]
  cl(0:1, 56 ) = [-1.0, 0.0520]
  cl(0:1, 57 ) = [-0.5, 0.1210]
  cl(0:1, 58 ) = [0.0, 0.1960]
  cl(0:1, 59 ) = [0.5, 0.2650]
  cl(0:1, 60 ) = [1.0, 0.3350]
  cl(0:1, 61 ) = [1.5, 0.4040]
  cl(0:1, 62 ) = [2.0, 0.4720]
  cl(0:1, 63 ) = [2.5, 0.5400]
  cl(0:1, 64 ) = [3.0, 0.6080]
  cl(0:1, 65 ) = [3.5, 0.6740]
  cl(0:1, 66 ) = [4.0, 0.7420]
  cl(0:1, 67 ) = [4.5, 0.8090]
  cl(0:1, 68 ) = [5.0, 0.8750]
  cl(0:1, 69 ) = [5.5, 0.9410]
  cl(0:1, 70 ) = [6.0, 1.0070]
  cl(0:1, 71 ) = [6.5, 1.0710]
  cl(0:1, 72 ) = [7.0, 1.1340]
  cl(0:1, 73 ) = [7.5, 1.1980]
  cl(0:1, 74 ) = [8.0, 1.2600]
  cl(0:1, 75 ) = [8.5, 1.3180]
  cl(0:1, 76 ) = [9.0, 1.3680]
  cl(0:1, 77 ) = [9.5, 1.4220]
  cl(0:1, 78 ) = [10.0, 1.4750]
  cl(0:1, 79 ) = [10.5, 1.5230]
  cl(0:1, 80 ) = [11.0, 1.5700]
  cl(0:1, 81 ) = [11.5, 1.6090]
  cl(0:1, 82 ) = [12.0, 1.6420]
  cl(0:1, 83 ) = [12.5, 1.6750]
  cl(0:1, 84 ) = [13.0, 1.7000]
  cl(0:1, 85 ) = [13.5, 1.7170]
  cl(0:1, 86 ) = [14.0, 1.7120]
  cl(0:1, 87 ) = [14.5, 1.7030]
  cl(0:1, 88 ) = [15.5, 1.6710]
  cl(0:1, 89 ) = [16.0, 1.6490]
  cl(0:1, 90 ) = [16.5, 1.6210]
  cl(0:1, 91 ) = [17.0, 1.5980]
  cl(0:1, 92 ) = [17.5, 1.5710]
  cl(0:1, 93 ) = [18.0, 1.5490]
  cl(0:1, 94 ) = [19.0, 1.5440]
  cl(0:1, 95 ) = [19.5, 1.5490]
  cl(0:1, 96 ) = [20.0, 1.5650]
  cl(0:1, 97 ) = [21.0, 1.5650]
  cl(0:1, 98 ) = [22.0, 1.5630]
  cl(0:1, 99 ) = [23.0, 1.5580]
  cl(0:1, 100) = [24.0, 1.5520]
  cl(0:1, 101) = [25.0, 1.5460]
  cl(0:1, 102) = [26.0, 1.5390]
  cl(0:1, 103) = [28.0, 1.5270]
  cl(0:1, 104) = [30.0, 1.5220]
  cl(0:1, 105) = [32.0, 1.5290]
  cl(0:1, 106) = [35.0, 1.5440]
  cl(0:1, 107) = [40.0, 1.5290]
  cl(0:1, 108) = [45.0, 1.4710]
  cl(0:1, 109) = [50.0, 1.3760]
  cl(0:1, 110) = [55.0, 1.2490]
  cl(0:1, 111) = [60.0, 1.0970]
  cl(0:1, 112) = [65.0, 0.9280]
  cl(0:1, 113) = [70.0, 0.7500]
  cl(0:1, 114) = [75.0, 0.5700]
  cl(0:1, 115) = [80.0, 0.3960]
  cl(0:1, 116) = [85.0, 0.2370]
  cl(0:1, 117) = [90.0, 0.1010]
  cl(0:1, 118) = [95.0, -0.0220]
  cl(0:1, 119) = [100.0, -0.1430]
  cl(0:1, 120) = [105.0, -0.2610]
  cl(0:1, 121) = [110.0, -0.3740]
  cl(0:1, 122) = [115.0, -0.4800]
  cl(0:1, 123) = [120.0, -0.5750]
  cl(0:1, 124) = [125.0, -0.6590]
  cl(0:1, 125) = [130.0, -0.7270]
  cl(0:1, 126) = [135.0, -0.7780]
  cl(0:1, 127) = [140.0, -0.8090]
  cl(0:1, 128) = [145.0, -0.8180]
  cl(0:1, 129) = [150.0, -0.8000]
  cl(0:1, 130) = [155.0, -0.7540]
  cl(0:1, 131) = [160.0, -0.6770]
  cl(0:1, 132) = [170.0, -0.4170]
  cl(0:1, 133) = [175.0, -0.2290]
  cl(0:1, 134) = [180.0, 0.0000]
 
  IF ( y .EQ. cl(0,0) ) cl4 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl4 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl4 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl4_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl4 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl4_nrel5.pro'
  ENDIF
  
!print,y,cl4
  
!_________________________________________________________________________
  
END FUNCTION calc_cl4_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl5_nrel5( y ) RESULT( cl5 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU30_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 143

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl5
  INTEGER :: i

  cl5 = 0.0

!                alpha    CL  
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.2740]
  cl(0:1, 2  ) = [-170.0, 0.5470]
  cl(0:1, 3  ) = [-160.0, 0.6850]
  cl(0:1, 4  ) = [-155.0, 0.7660]
  cl(0:1, 5  ) = [-150.0, 0.8160]
  cl(0:1, 6  ) = [-145.0, 0.8360]
  cl(0:1, 7  ) = [-140.0, 0.8320]
  cl(0:1, 8  ) = [-135.0, 0.8040]
  cl(0:1, 9  ) = [-130.0, 0.7560]
  cl(0:1, 10 ) = [-125.0, 0.6900]
  cl(0:1, 11 ) = [-120.0, 0.6090]
  cl(0:1, 12 ) = [-115.0, 0.5150]
  cl(0:1, 13 ) = [-110.0, 0.4110]
  cl(0:1, 14 ) = [-105.0, 0.3000]
  cl(0:1, 15 ) = [-100.0, 0.1820]
  cl(0:1, 16 ) = [-95.0, 0.0610]
  cl(0:1, 17 ) = [-90.0, -0.0610]
  cl(0:1, 18 ) = [-85.0, -0.1830]
  cl(0:1, 19 ) = [-80.0, -0.3020]
  cl(0:1, 20 ) = [-75.0, -0.4160]
  cl(0:1, 21 ) = [-70.0, -0.5230]
  cl(0:1, 22 ) = [-65.0, -0.6220]
  cl(0:1, 23 ) = [-60.0, -0.7080]
  cl(0:1, 24 ) = [-55.0, -0.7810]
  cl(0:1, 25 ) = [-50.0, -0.8380]
  cl(0:1, 26 ) = [-45.0, -0.8770]
  cl(0:1, 27 ) = [-40.0, -0.8950]
  cl(0:1, 28 ) = [-35.0, -0.8890]
  cl(0:1, 29 ) = [-30.0, -0.8580]
  cl(0:1, 30 ) = [-25.0, -0.8320]
  cl(0:1, 31 ) = [-24.0, -0.8520]
  cl(0:1, 32 ) = [-23.0, -0.8820]
  cl(0:1, 33 ) = [-22.0, -0.9190]
  cl(0:1, 34 ) = [-21.0, -0.9630]
  cl(0:1, 35 ) = [-20.0, -1.0130]
  cl(0:1, 36 ) = [-19.0, -1.0670]
  cl(0:1, 37 ) = [-18.0, -1.1250]
  cl(0:1, 38 ) = [-17.0, -1.1850]
  cl(0:1, 39 ) = [-16.0, -1.2450]
  cl(0:1, 40 ) = [-15.3, -1.2900]
  cl(0:1, 41 ) = [-14.2, -1.2290]
  cl(0:1, 42 ) = [-13.2, -1.1480]
  cl(0:1, 43 ) = [-12.2, -1.0520]
  cl(0:1, 44 ) = [-11.2, -0.9650]
  cl(0:1, 45 ) = [-10.2, -0.8670]
  cl(0:1, 46 ) = [-9.7, -0.8220]
  cl(0:1, 47 ) = [-9.2, -0.7690]
  cl(0:1, 48 ) = [-8.2, -0.7560]
  cl(0:1, 49 ) = [-7.2, -0.6900]
  cl(0:1, 50 ) = [-6.7, -0.6160]
  cl(0:1, 51 ) = [-6.1, -0.5420]
  cl(0:1, 52 ) = [-6.0, -0.5250]
  cl(0:1, 53 ) = [-5.5, -0.4510]
  cl(0:1, 54 ) = [-5.0, -0.3820]
  cl(0:1, 55 ) = [-4.5, -0.3140]
  cl(0:1, 56 ) = [-4.0, -0.2510]
  cl(0:1, 57 ) = [-3.5, -0.1890]
  cl(0:1, 58 ) = [-3.0, -0.1200]
  cl(0:1, 59 ) = [-2.5, -0.0510]
  cl(0:1, 60 ) = [-2.0, 0.0170]
  cl(0:1, 61 ) = [-1.5, 0.0850]
  cl(0:1, 62 ) = [-1.0, 0.1520]
  cl(0:1, 63 ) = [-0.5, 0.2190]
  cl(0:1, 64 ) = [0.0, 0.2880]
  cl(0:1, 65 ) = [0.5, 0.3540]
  cl(0:1, 66 ) = [1.0, 0.4210]
  cl(0:1, 67 ) = [1.5, 0.4870]
  cl(0:1, 68 ) = [2.0, 0.5540]
  cl(0:1, 69 ) = [2.5, 0.6190]
  cl(0:1, 70 ) = [3.0, 0.6850]
  cl(0:1, 71 ) = [3.5, 0.7490]
  cl(0:1, 72 ) = [4.0, 0.8150]
  cl(0:1, 73 ) = [4.5, 0.8790]
  cl(0:1, 74 ) = [5.0, 0.9440]
  cl(0:1, 75 ) = [5.5, 1.0080]
  cl(0:1, 76 ) = [6.0, 1.0720]
  cl(0:1, 77 ) = [6.5, 1.1350]
  cl(0:1, 78 ) = [7.0, 1.1970]
  cl(0:1, 79 ) = [7.5, 1.2560]
  cl(0:1, 80 ) = [8.0, 1.3050]
  cl(0:1, 81 ) = [9.0, 1.3900]
  cl(0:1, 82 ) = [9.5, 1.4240]
  cl(0:1, 83 ) = [10.0, 1.4580]
  cl(0:1, 84 ) = [10.5, 1.4880]
  cl(0:1, 85 ) = [11.0, 1.5120]
  cl(0:1, 86 ) = [11.5, 1.5330]
  cl(0:1, 87 ) = [12.0, 1.5490]
  cl(0:1, 88 ) = [12.5, 1.5580]
  cl(0:1, 89 ) = [13.0, 1.4700]
  cl(0:1, 90 ) = [13.5, 1.3980]
  cl(0:1, 91 ) = [14.0, 1.3540]
  cl(0:1, 92 ) = [14.5, 1.3360]
  cl(0:1, 93 ) = [15.0, 1.3330]
  cl(0:1, 94 ) = [15.5, 1.3260]
  cl(0:1, 95 ) = [16.0, 1.3290]
  cl(0:1, 96 ) = [16.5, 1.3260]
  cl(0:1, 97 ) = [17.0, 1.3210]
  cl(0:1, 98 ) = [17.5, 1.3310]
  cl(0:1, 99 ) = [18.0, 1.3330]
  cl(0:1, 100) = [18.5, 1.3400]
  cl(0:1, 101) = [19.0, 1.3620]
  cl(0:1, 102) = [19.5, 1.3820]
  cl(0:1, 103) = [20.0, 1.3980]
  cl(0:1, 104) = [20.5, 1.4260]
  cl(0:1, 105) = [21.0, 1.4370]
  cl(0:1, 106) = [22.0, 1.4180]
  cl(0:1, 107) = [23.0, 1.3970]
  cl(0:1, 108) = [24.0, 1.3760]
  cl(0:1, 109) = [25.0, 1.3540]
  cl(0:1, 110) = [26.0, 1.3320]
  cl(0:1, 111) = [28.0, 1.2930]
  cl(0:1, 112) = [30.0, 1.2650]
  cl(0:1, 113) = [32.0, 1.2530]
  cl(0:1, 114) = [35.0, 1.2640]
  cl(0:1, 115) = [40.0, 1.2580]
  cl(0:1, 116) = [45.0, 1.2170]
  cl(0:1, 117) = [50.0, 1.1460]
  cl(0:1, 118) = [55.0, 1.0490]
  cl(0:1, 119) = [60.0, 0.9320]
  cl(0:1, 120) = [65.0, 0.7990]
  cl(0:1, 121) = [70.0, 0.6570]
  cl(0:1, 122) = [75.0, 0.5090]
  cl(0:1, 123) = [80.0, 0.3620]
  cl(0:1, 124) = [85.0, 0.2210]
  cl(0:1, 125) = [90.0, 0.0920]
  cl(0:1, 126) = [95.0, -0.0300]
  cl(0:1, 127) = [100.0, -0.1500]
  cl(0:1, 128) = [105.0, -0.2670]
  cl(0:1, 129) = [110.0, -0.3790]
  cl(0:1, 130) = [115.0, -0.4830]
  cl(0:1, 131) = [120.0, -0.5780]
  cl(0:1, 132) = [125.0, -0.6600]
  cl(0:1, 133) = [130.0, -0.7270]
  cl(0:1, 134) = [135.0, -0.7770]
  cl(0:1, 135) = [140.0, -0.8070]
  cl(0:1, 136) = [145.0, -0.8150]
  cl(0:1, 137) = [150.0, -0.7970]
  cl(0:1, 138) = [155.0, -0.7500]
  cl(0:1, 139) = [160.0, -0.6730]
  cl(0:1, 140) = [170.0, -0.5470]
  cl(0:1, 141) = [175.0, -0.2740]
  cl(0:1, 142) = [180.0, 0.0000]

  IF ( y .EQ. cl(0,0) ) cl5 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl5 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl5 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl5_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl5 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl5_nrel5.pro'
  ENDIF
  
!print,y,cl5
  
!_________________________________________________________________________
  
END FUNCTION calc_cl5_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl6_nrel5( y ) RESULT( cl6 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU25_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 140

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl6
  INTEGER :: i

  cl6 = 0.0

!                alpha    CL  
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.3680]
  cl(0:1, 2  ) = [-170.0, 0.7350]
  cl(0:1, 3  ) = [-160.0, 0.6950]
  cl(0:1, 4  ) = [-155.0, 0.7770]
  cl(0:1, 5  ) = [-150.0, 0.8280]
  cl(0:1, 6  ) = [-145.0, 0.8500]
  cl(0:1, 7  ) = [-140.0, 0.8460]
  cl(0:1, 8  ) = [-135.0, 0.8180]
  cl(0:1, 9  ) = [-130.0, 0.7710]
  cl(0:1, 10 ) = [-125.0, 0.7050]
  cl(0:1, 11 ) = [-120.0, 0.6240]
  cl(0:1, 12 ) = [-115.0, 0.5300]
  cl(0:1, 13 ) = [-110.0, 0.4260]
  cl(0:1, 14 ) = [-105.0, 0.3140]
  cl(0:1, 15 ) = [-100.0, 0.1950]
  cl(0:1, 16 ) = [-95.0, 0.0730]
  cl(0:1, 17 ) = [-90.0, -0.0500]
  cl(0:1, 18 ) = [-85.0, -0.1730]
  cl(0:1, 19 ) = [-80.0, -0.2940]
  cl(0:1, 20 ) = [-75.0, -0.4090]
  cl(0:1, 21 ) = [-70.0, -0.5180]
  cl(0:1, 22 ) = [-65.0, -0.6170]
  cl(0:1, 23 ) = [-60.0, -0.7060]
  cl(0:1, 24 ) = [-55.0, -0.7800]
  cl(0:1, 25 ) = [-50.0, -0.8390]
  cl(0:1, 26 ) = [-45.0, -0.8790]
  cl(0:1, 27 ) = [-40.0, -0.8980]
  cl(0:1, 28 ) = [-35.0, -0.8930]
  cl(0:1, 29 ) = [-30.0, -0.8620]
  cl(0:1, 30 ) = [-25.0, -0.8030]
  cl(0:1, 31 ) = [-24.0, -0.7920]
  cl(0:1, 32 ) = [-23.0, -0.7890]
  cl(0:1, 33 ) = [-22.0, -0.7920]
  cl(0:1, 34 ) = [-21.0, -0.8010]
  cl(0:1, 35 ) = [-20.0, -0.8150]
  cl(0:1, 36 ) = [-19.0, -0.8330]
  cl(0:1, 37 ) = [-18.0, -0.8540]
  cl(0:1, 38 ) = [-17.0, -0.8790]
  cl(0:1, 39 ) = [-16.0, -0.9050]
  cl(0:1, 40 ) = [-15.0, -0.9320]
  cl(0:1, 41 ) = [-14.0, -0.9590]
  cl(0:1, 42 ) = [-13.0, -0.9850]
  cl(0:1, 43 ) = [-12.0, -0.9530]
  cl(0:1, 44 ) = [-11.0, -0.9000]
  cl(0:1, 45 ) = [-10.0, -0.8270]
  cl(0:1, 46 ) = [-9.0, -0.7530]
  cl(0:1, 47 ) = [-8.5, -0.6910]
  cl(0:1, 48 ) = [-7.5, -0.5550]
  cl(0:1, 49 ) = [-6.4, -0.4130]
  cl(0:1, 50 ) = [-5.4, -0.2710]
  cl(0:1, 51 ) = [-5.0, -0.2200]
  cl(0:1, 52 ) = [-4.5, -0.1520]
  cl(0:1, 53 ) = [-4.0, -0.0840]
  cl(0:1, 54 ) = [-3.5, -0.0180]
  cl(0:1, 55 ) = [-3.0, 0.0490]
  cl(0:1, 56 ) = [-2.5, 0.1150]
  cl(0:1, 57 ) = [-2.0, 0.1810]
  cl(0:1, 58 ) = [-1.5, 0.2470]
  cl(0:1, 59 ) = [-1.0, 0.3120]
  cl(0:1, 60 ) = [-0.5, 0.3770]
  cl(0:1, 61 ) = [0.0, 0.4440]
  cl(0:1, 62 ) = [0.5, 0.5080]
  cl(0:1, 63 ) = [1.0, 0.5730]
  cl(0:1, 64 ) = [1.5, 0.6360]
  cl(0:1, 65 ) = [2.0, 0.7010]
  cl(0:1, 66 ) = [2.5, 0.7650]
  cl(0:1, 67 ) = [3.0, 0.8270]
  cl(0:1, 68 ) = [3.5, 0.8900]
  cl(0:1, 69 ) = [4.0, 0.9520]
  cl(0:1, 70 ) = [4.5, 1.0130]
  cl(0:1, 71 ) = [5.0, 1.0620]
  cl(0:1, 72 ) = [6.0, 1.1610]
  cl(0:1, 73 ) = [6.5, 1.2080]
  cl(0:1, 74 ) = [7.0, 1.2540]
  cl(0:1, 75 ) = [7.5, 1.3010]
  cl(0:1, 76 ) = [8.0, 1.3360]
  cl(0:1, 77 ) = [8.5, 1.3690]
  cl(0:1, 78 ) = [9.0, 1.4000]
  cl(0:1, 79 ) = [9.5, 1.4280]
  cl(0:1, 80 ) = [10.0, 1.4420]
  cl(0:1, 81 ) = [10.5, 1.4270]
  cl(0:1, 82 ) = [11.0, 1.3740]
  cl(0:1, 83 ) = [11.5, 1.3160]
  cl(0:1, 84 ) = [12.0, 1.2770]
  cl(0:1, 85 ) = [12.5, 1.2500]
  cl(0:1, 86 ) = [13.0, 1.2460]
  cl(0:1, 87 ) = [13.5, 1.2470]
  cl(0:1, 88 ) = [14.0, 1.2560]
  cl(0:1, 89 ) = [14.5, 1.2600]
  cl(0:1, 90 ) = [15.0, 1.2710]
  cl(0:1, 91 ) = [15.5, 1.2810]
  cl(0:1, 92 ) = [16.0, 1.2890]
  cl(0:1, 93 ) = [16.5, 1.2940]
  cl(0:1, 94 ) = [17.0, 1.3040]
  cl(0:1, 95 ) = [17.5, 1.3090]
  cl(0:1, 96 ) = [18.0, 1.3150]
  cl(0:1, 97 ) = [18.5, 1.3200]
  cl(0:1, 98 ) = [19.0, 1.3300]
  cl(0:1, 99 ) = [19.5, 1.3430]
  cl(0:1, 100) = [20.0, 1.3540]
  cl(0:1, 101) = [20.5, 1.3590]
  cl(0:1, 102) = [21.0, 1.3600]
  cl(0:1, 103) = [22.0, 1.3250]
  cl(0:1, 104) = [23.0, 1.2880]
  cl(0:1, 105) = [24.0, 1.2510]
  cl(0:1, 106) = [25.0, 1.2150]
  cl(0:1, 107) = [26.0, 1.1810]
  cl(0:1, 108) = [28.0, 1.1200]
  cl(0:1, 109) = [30.0, 1.0760]
  cl(0:1, 110) = [32.0, 1.0560]
  cl(0:1, 111) = [35.0, 1.0660]
  cl(0:1, 112) = [40.0, 1.0640]
  cl(0:1, 113) = [45.0, 1.0350]
  cl(0:1, 114) = [50.0, 0.9800]
  cl(0:1, 115) = [55.0, 0.9040]
  cl(0:1, 116) = [60.0, 0.8100]
  cl(0:1, 117) = [65.0, 0.7020]
  cl(0:1, 118) = [70.0, 0.5820]
  cl(0:1, 119) = [75.0, 0.4560]
  cl(0:1, 120) = [80.0, 0.3260]
  cl(0:1, 121) = [85.0, 0.1970]
  cl(0:1, 122) = [90.0, 0.0720]
  cl(0:1, 123) = [95.0, -0.0500]
  cl(0:1, 124) = [100.0, -0.1700]
  cl(0:1, 125) = [105.0, -0.2870]
  cl(0:1, 126) = [110.0, -0.3990]
  cl(0:1, 127) = [115.0, -0.5020]
  cl(0:1, 128) = [120.0, -0.5960]
  cl(0:1, 129) = [125.0, -0.6770]
  cl(0:1, 130) = [130.0, -0.7430]
  cl(0:1, 131) = [135.0, -0.7920]
  cl(0:1, 132) = [140.0, -0.8210]
  cl(0:1, 133) = [145.0, -0.8260]
  cl(0:1, 134) = [150.0, -0.8060]
  cl(0:1, 135) = [155.0, -0.7580]
  cl(0:1, 136) = [160.0, -0.6790]
  cl(0:1, 137) = [170.0, -0.7350]
  cl(0:1, 138) = [175.0, -0.3680]
  cl(0:1, 139) = [180.0, 0.0000]

  IF ( y .EQ. cl(0,0) ) cl6 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl6 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl6 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl6_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl6 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl6_nrel5.pro'
  ENDIF
  
!print,y,cl6
  
!_________________________________________________________________________
  
END FUNCTION calc_cl6_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cl7_nrel5( y ) RESULT( cl7 )
!_________________________________________________________________________
!
! Data from  NREL 5 MW: DU21_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 142

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl7
  INTEGER :: i

  cl7 = 0.0

!                alpha    CL  
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.3940]
  cl(0:1, 2  ) = [-170.0, 0.7880]
  cl(0:1, 3  ) = [-160.0, 0.6700]
  cl(0:1, 4  ) = [-155.0, 0.7490]
  cl(0:1, 5  ) = [-150.0, 0.7970]
  cl(0:1, 6  ) = [-145.0, 0.8180]
  cl(0:1, 7  ) = [-140.0, 0.8130]
  cl(0:1, 8  ) = [-135.0, 0.7860]
  cl(0:1, 9  ) = [-130.0, 0.7390]
  cl(0:1, 10 ) = [-125.0, 0.6750]
  cl(0:1, 11 ) = [-120.0, 0.5960]
  cl(0:1, 12 ) = [-115.0, 0.5050]
  cl(0:1, 13 ) = [-110.0, 0.4030]
  cl(0:1, 14 ) = [-105.0, 0.2940]
  cl(0:1, 15 ) = [-100.0, 0.1790]
  cl(0:1, 16 ) = [-95.0, 0.0600]
  cl(0:1, 17 ) = [-90.0, -0.0600]
  cl(0:1, 18 ) = [-85.0, -0.1790]
  cl(0:1, 19 ) = [-80.0, -0.2950]
  cl(0:1, 20 ) = [-75.0, -0.4070]
  cl(0:1, 21 ) = [-70.0, -0.5120]
  cl(0:1, 22 ) = [-65.0, -0.6080]
  cl(0:1, 23 ) = [-60.0, -0.6930]
  cl(0:1, 24 ) = [-55.0, -0.7640]
  cl(0:1, 25 ) = [-50.0, -0.8200]
  cl(0:1, 26 ) = [-45.0, -0.8570]
  cl(0:1, 27 ) = [-40.0, -0.8750]
  cl(0:1, 28 ) = [-35.0, -0.8690]
  cl(0:1, 29 ) = [-30.0, -0.8380]
  cl(0:1, 30 ) = [-25.0, -0.7910]
  cl(0:1, 31 ) = [-24.0, -0.7940]
  cl(0:1, 32 ) = [-23.0, -0.8050]
  cl(0:1, 33 ) = [-22.0, -0.8210]
  cl(0:1, 34 ) = [-21.0, -0.8430]
  cl(0:1, 35 ) = [-20.0, -0.8690]
  cl(0:1, 36 ) = [-19.0, -0.8990]
  cl(0:1, 37 ) = [-18.0, -0.9310]
  cl(0:1, 38 ) = [-17.0, -0.9640]
  cl(0:1, 39 ) = [-16.0, -0.9990]
  cl(0:1, 40 ) = [-15.0, -1.0330]
  cl(0:1, 41 ) = [-14.5, -1.0500]
  cl(0:1, 42 ) = [-12.0, -0.9530]
  cl(0:1, 43 ) = [-11.0, -0.9000]
  cl(0:1, 44 ) = [-10.0, -0.8270]
  cl(0:1, 45 ) = [-8.1, -0.5360]
  cl(0:1, 46 ) = [-7.6, -0.4670]
  cl(0:1, 47 ) = [-7.1, -0.3930]
  cl(0:1, 48 ) = [-6.6, -0.3230]
  cl(0:1, 49 ) = [-6.5, -0.3110]
  cl(0:1, 50 ) = [-6.0, -0.2450]
  cl(0:1, 51 ) = [-5.5, -0.1780]
  cl(0:1, 52 ) = [-5.0, -0.1130]
  cl(0:1, 53 ) = [-4.5, -0.0480]
  cl(0:1, 54 ) = [-4.0, 0.0160]
  cl(0:1, 55 ) = [-3.5, 0.0800]
  cl(0:1, 56 ) = [-3.0, 0.1450]
  cl(0:1, 57 ) = [-2.5, 0.2080]
  cl(0:1, 58 ) = [-2.0, 0.2700]
  cl(0:1, 59 ) = [-1.5, 0.3330]
  cl(0:1, 60 ) = [-1.0, 0.3960]
  cl(0:1, 61 ) = [-0.5, 0.4580]
  cl(0:1, 62 ) = [0.0, 0.5210]
  cl(0:1, 63 ) = [0.5, 0.5830]
  cl(0:1, 64 ) = [1.0, 0.6450]
  cl(0:1, 65 ) = [1.5, 0.7060]
  cl(0:1, 66 ) = [2.0, 0.7680]
  cl(0:1, 67 ) = [2.5, 0.8280]
  cl(0:1, 68 ) = [3.0, 0.8880]
  cl(0:1, 69 ) = [3.5, 0.9480]
  cl(0:1, 70 ) = [4.0, 0.9960]
  cl(0:1, 71 ) = [4.5, 1.0460]
  cl(0:1, 72 ) = [5.0, 1.0950]
  cl(0:1, 73 ) = [5.5, 1.1450]
  cl(0:1, 74 ) = [6.0, 1.1920]
  cl(0:1, 75 ) = [6.5, 1.2390]
  cl(0:1, 76 ) = [7.0, 1.2830]
  cl(0:1, 77 ) = [7.5, 1.3240]
  cl(0:1, 78 ) = [8.0, 1.3580]
  cl(0:1, 79 ) = [8.5, 1.3850]
  cl(0:1, 80 ) = [9.0, 1.4030]
  cl(0:1, 81 ) = [9.5, 1.4010]
  cl(0:1, 82 ) = [10.0, 1.3580]
  cl(0:1, 83 ) = [10.5, 1.3130]
  cl(0:1, 84 ) = [11.0, 1.2870]
  cl(0:1, 85 ) = [11.5, 1.2740]
  cl(0:1, 86 ) = [12.0, 1.2720]
  cl(0:1, 87 ) = [12.5, 1.2730]
  cl(0:1, 88 ) = [13.0, 1.2730]
  cl(0:1, 89 ) = [13.5, 1.2730]
  cl(0:1, 90 ) = [14.0, 1.2720]
  cl(0:1, 91 ) = [14.5, 1.2730]
  cl(0:1, 92 ) = [15.0, 1.2750]
  cl(0:1, 93 ) = [15.5, 1.2810]
  cl(0:1, 94 ) = [16.0, 1.2840]
  cl(0:1, 95 ) = [16.5, 1.2960]
  cl(0:1, 96 ) = [17.0, 1.3060]
  cl(0:1, 97 ) = [17.5, 1.3080]
  cl(0:1, 98 ) = [18.0, 1.3080]
  cl(0:1, 99 ) = [18.5, 1.3080]
  cl(0:1, 100) = [19.0, 1.3080]
  cl(0:1, 101) = [19.5, 1.3070]
  cl(0:1, 102) = [20.0, 1.3110]
  cl(0:1, 103) = [20.5, 1.3250]
  cl(0:1, 104) = [21.0, 1.3240]
  cl(0:1, 105) = [22.0, 1.2770]
  cl(0:1, 106) = [23.0, 1.2290]
  cl(0:1, 107) = [24.0, 1.1820]
  cl(0:1, 108) = [25.0, 1.1360]
  cl(0:1, 109) = [26.0, 1.0930]
  cl(0:1, 110) = [28.0, 1.0170]
  cl(0:1, 111) = [30.0, 0.9620]
  cl(0:1, 112) = [32.0, 0.9370]
  cl(0:1, 113) = [35.0, 0.9470]
  cl(0:1, 114) = [40.0, 0.9500]
  cl(0:1, 115) = [45.0, 0.9280]
  cl(0:1, 116) = [50.0, 0.8840]
  cl(0:1, 117) = [55.0, 0.8210]
  cl(0:1, 118) = [60.0, 0.7400]
  cl(0:1, 119) = [65.0, 0.6460]
  cl(0:1, 120) = [70.0, 0.5400]
  cl(0:1, 121) = [75.0, 0.4250]
  cl(0:1, 122) = [80.0, 0.3040]
  cl(0:1, 123) = [85.0, 0.1790]
  cl(0:1, 124) = [90.0, 0.0530]
  cl(0:1, 125) = [95.0, -0.0730]
  cl(0:1, 126) = [100.0, -0.1980]
  cl(0:1, 127) = [105.0, -0.3190]
  cl(0:1, 128) = [110.0, -0.4340]
  cl(0:1, 129) = [115.0, -0.5410]
  cl(0:1, 130) = [120.0, -0.6370]
  cl(0:1, 131) = [125.0, -0.7200]
  cl(0:1, 132) = [130.0, -0.7870]
  cl(0:1, 133) = [135.0, -0.8360]
  cl(0:1, 134) = [140.0, -0.8640]
  cl(0:1, 135) = [145.0, -0.8690]
  cl(0:1, 136) = [150.0, -0.8470]
  cl(0:1, 137) = [155.0, -0.7950]
  cl(0:1, 138) = [160.0, -0.7110]
  cl(0:1, 139) = [170.0, -0.7880]
  cl(0:1, 140) = [175.0, -0.3940]
  cl(0:1, 141) = [180.0, 0.0000]

  IF ( y .EQ. cl(0,0) ) cl7 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl7 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl7 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl7_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl7 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl7_nrel5.pro'
  ENDIF
  
!print,y,cl7
  
!_________________________________________________________________________
  
END FUNCTION calc_cl7_nrel5
!_________________________________________________________________________

!=========================================================================

Function calc_cl8_nrel5( y ) RESULT( cl8 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: NACA64_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 127

  REAL, DIMENSION( 0:1,0:na-1 ) :: cl

  REAL :: cl8
  INTEGER :: i

  cl8 = 0.0

!                alpha    CL  
  cl(0:1, 0  ) = [-180.0, 0.0000]
  cl(0:1, 1  ) = [-175.0, 0.3740]
  cl(0:1, 2  ) = [-170.0, 0.7490]
  cl(0:1, 3  ) = [-160.0, 0.6590]
  cl(0:1, 4  ) = [-155.0, 0.7360]
  cl(0:1, 5  ) = [-150.0, 0.7830]
  cl(0:1, 6  ) = [-145.0, 0.8030]
  cl(0:1, 7  ) = [-140.0, 0.7980]
  cl(0:1, 8  ) = [-135.0, 0.7710]
  cl(0:1, 9  ) = [-130.0, 0.7240]
  cl(0:1, 10 ) = [-125.0, 0.6600]
  cl(0:1, 11 ) = [-120.0, 0.5810]
  cl(0:1, 12 ) = [-115.0, 0.4910]
  cl(0:1, 13 ) = [-110.0, 0.3900]
  cl(0:1, 14 ) = [-105.0, 0.2820]
  cl(0:1, 15 ) = [-100.0, 0.1690]
  cl(0:1, 16 ) = [-95.0, 0.0520]
  cl(0:1, 17 ) = [-90.0, -0.0670]
  cl(0:1, 18 ) = [-85.0, -0.1840]
  cl(0:1, 19 ) = [-80.0, -0.2990]
  cl(0:1, 20 ) = [-75.0, -0.4090]
  cl(0:1, 21 ) = [-70.0, -0.5120]
  cl(0:1, 22 ) = [-65.0, -0.6060]
  cl(0:1, 23 ) = [-60.0, -0.6890]
  cl(0:1, 24 ) = [-55.0, -0.7590]
  cl(0:1, 25 ) = [-50.0, -0.8140]
  cl(0:1, 26 ) = [-45.0, -0.8500]
  cl(0:1, 27 ) = [-40.0, -0.8660]
  cl(0:1, 28 ) = [-35.0, -0.8600]
  cl(0:1, 29 ) = [-30.0, -0.8290]
  cl(0:1, 30 ) = [-25.0, -0.8530]
  cl(0:1, 31 ) = [-24.0, -0.8700]
  cl(0:1, 32 ) = [-23.0, -0.8900]
  cl(0:1, 33 ) = [-22.0, -0.9110]
  cl(0:1, 34 ) = [-21.0, -0.9340]
  cl(0:1, 35 ) = [-20.0, -0.9580]
  cl(0:1, 36 ) = [-19.0, -0.9820]
  cl(0:1, 37 ) = [-18.0, -1.0050]
  cl(0:1, 38 ) = [-17.0, -1.0820]
  cl(0:1, 39 ) = [-16.0, -1.1130]
  cl(0:1, 40 ) = [-15.0, -1.1050]
  cl(0:1, 41 ) = [-14.0, -1.0780]
  cl(0:1, 42 ) = [-13.5, -1.0530]
  cl(0:1, 43 ) = [-13.0, -1.0150]
  cl(0:1, 44 ) = [-12.0, -0.9040]
  cl(0:1, 45 ) = [-11.0, -0.8070]
  cl(0:1, 46 ) = [-10.0, -0.7110]
  cl(0:1, 47 ) = [-9.0, -0.5950]
  cl(0:1, 48 ) = [-8.0, -0.4780]
  cl(0:1, 49 ) = [-7.0, -0.3750]
  cl(0:1, 50 ) = [-6.0, -0.2640]
  cl(0:1, 51 ) = [-5.0, -0.1510]
  cl(0:1, 52 ) = [-4.0, -0.0170]
  cl(0:1, 53 ) = [-3.0, 0.0880]
  cl(0:1, 54 ) = [-2.0, 0.2130]
  cl(0:1, 55 ) = [-1.0, 0.3280]
  cl(0:1, 56 ) = [0.0, 0.4420]
  cl(0:1, 57 ) = [1.0, 0.5560]
  cl(0:1, 58 ) = [2.0, 0.6700]
  cl(0:1, 59 ) = [3.0, 0.7840]
  cl(0:1, 60 ) = [4.0, 0.8980]
  cl(0:1, 61 ) = [5.0, 1.0110]
  cl(0:1, 62 ) = [6.0, 1.1030]
  cl(0:1, 63 ) = [7.0, 1.1810]
  cl(0:1, 64 ) = [8.0, 1.2570]
  cl(0:1, 65 ) = [8.5, 1.2930]
  cl(0:1, 66 ) = [9.0, 1.3260]
  cl(0:1, 67 ) = [9.5, 1.3560]
  cl(0:1, 68 ) = [10.0, 1.3820]
  cl(0:1, 69 ) = [10.5, 1.4000]
  cl(0:1, 70 ) = [11.0, 1.4150]
  cl(0:1, 71 ) = [11.5, 1.4250]
  cl(0:1, 72 ) = [12.0, 1.4340]
  cl(0:1, 73 ) = [12.5, 1.4430]
  cl(0:1, 74 ) = [13.0, 1.4510]
  cl(0:1, 75 ) = [13.5, 1.4530]
  cl(0:1, 76 ) = [14.0, 1.4480]
  cl(0:1, 77 ) = [14.5, 1.4440]
  cl(0:1, 78 ) = [15.0, 1.4450]
  cl(0:1, 79 ) = [15.5, 1.4470]
  cl(0:1, 80 ) = [16.0, 1.4480]
  cl(0:1, 81 ) = [16.5, 1.4440]
  cl(0:1, 82 ) = [17.0, 1.4380]
  cl(0:1, 83 ) = [17.5, 1.4390]
  cl(0:1, 84 ) = [18.0, 1.4480]
  cl(0:1, 85 ) = [18.5, 1.4520]
  cl(0:1, 86 ) = [19.0, 1.4480]
  cl(0:1, 87 ) = [19.5, 1.4380]
  cl(0:1, 88 ) = [20.0, 1.4280]
  cl(0:1, 89 ) = [21.0, 1.4010]
  cl(0:1, 90 ) = [22.0, 1.3590]
  cl(0:1, 91 ) = [23.0, 1.3000]
  cl(0:1, 92 ) = [24.0, 1.2200]
  cl(0:1, 93 ) = [25.0, 1.1680]
  cl(0:1, 94 ) = [26.0, 1.1160]
  cl(0:1, 95 ) = [28.0, 1.0150]
  cl(0:1, 96 ) = [30.0, 0.9260]
  cl(0:1, 97 ) = [32.0, 0.8550]
  cl(0:1, 98 ) = [35.0, 0.8000]
  cl(0:1, 99 ) = [40.0, 0.8040]
  cl(0:1, 100) = [45.0, 0.7930]
  cl(0:1, 101) = [50.0, 0.7630]
  cl(0:1, 102) = [55.0, 0.7170]
  cl(0:1, 103) = [60.0, 0.6560]
  cl(0:1, 104) = [65.0, 0.5820]
  cl(0:1, 105) = [70.0, 0.4950]
  cl(0:1, 106) = [75.0, 0.3980]
  cl(0:1, 107) = [80.0, 0.2910]
  cl(0:1, 108) = [85.0, 0.1760]
  cl(0:1, 109) = [90.0, 0.0530]
  cl(0:1, 110) = [95.0, -0.0740]
  cl(0:1, 111) = [100.0, -0.1990]
  cl(0:1, 112) = [105.0, -0.3210]
  cl(0:1, 113) = [110.0, -0.4360]
  cl(0:1, 114) = [115.0, -0.5430]
  cl(0:1, 115) = [120.0, -0.6400]
  cl(0:1, 116) = [125.0, -0.7230]
  cl(0:1, 117) = [130.0, -0.7900]
  cl(0:1, 118) = [135.0, -0.8400]
  cl(0:1, 119) = [140.0, -0.8680]
  cl(0:1, 120) = [145.0, -0.8720]
  cl(0:1, 121) = [150.0, -0.8500]
  cl(0:1, 122) = [155.0, -0.7980]
  cl(0:1, 123) = [160.0, -0.7140]
  cl(0:1, 124) = [170.0, -0.7490]
  cl(0:1, 125) = [175.0, -0.3740]
  cl(0:1, 126) = [180.0, 0.0000]
  
  IF ( y .EQ. cl(0,0) ) cl8 = cl(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cl(0,i-1) ) .AND. (y .LE. cl(0,i) ) ) THEN
        cl8 = cl(1,i-1) + (y - cl(0,i-1))*( ( cl(1,i) - cl(1,i-1) )/( cl(0,i) - cl(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cl(0,0) ) THEN 
     cl8 = cl(1,0)
!     PRINT*,'alpha < tabulated range in calc_cl8_nrel5.pro'
  ENDIF
  IF ( y .GT. cl(0,na-1) ) THEN
     cl8 = cl(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cl8_nrel5.pro'
  ENDIF
  
!print,y,cl8
  
!_________________________________________________________________________
  
END FUNCTION calc_cl8_nrel5
!_________________________________________________________________________


!=========================================================================

FUNCTION calc_cd_nrel5( x,y ) RESULT( cd )
!_________________________________________________________________________
!
! Calls appropriate airfoil section table depending upon radial distance
! from blade root. Data from NREL 5 MW
!
! calc_cd1_nrel5 uses data from Cylinder1.dat
! calc_cd2_nrel5 uses data from Cylinder2.dat
! calc_cd3_nrel5 uses data from DU40_A17.dat
! calc_cd4_nrel5 uses data from DU35_A17.dat
! calc_cd5_nrel5 uses data from DU30_A17.dat
! calc_cd6_nrel5 uses data from DU25_A17.dat
! calc_cd7_nrel5 uses data from DU21_A17.dat
! calc_cd8_nrel5 uses data from NACA64_A17.dat
!
! x = radial distance from hub center
! y = angle of attack
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: x,y

  INTEGER, PARAMETER :: ns = 9

  REAL, DIMENSION( 0:ns-1 ) :: r

  REAL :: cd
  INTEGER :: i

  cd = 0.0

  r(0)  = 1.5
  r(1)  = 8.3333
  r(2)  = 11.75
  r(3)  = 15.85
  r(4)  = 24.05
  r(5)  = 28.15
  r(6)  = 36.35
  r(7)  = 44.55
  r(8)  = 63.000

  IF ( ( x .GT. r( 0) ) .AND. (x .LE. r( 1) ) ) THEN
       cd = calc_cd1_nrel5(y)
     ELSEIF ( ( x .GT. r( 1) ) .AND. (x .LE. r( 2) ) ) THEN
       cd = calc_cd2_nrel5(y)
     ELSEIF ( ( x .GT. r( 2) ) .AND. (x .LE. r( 3) ) ) THEN
       cd = calc_cd3_nrel5(y)
     ELSEIF ( ( x .GT. r( 3) ) .AND. (x .LE. r( 4) ) ) THEN
       cd = calc_cd4_nrel5(y)
     ELSEIF ( ( x .GT. r( 4) ) .AND. (x .LE. r( 5) ) ) THEN
       cd = calc_cd5_nrel5(y)
     ELSEIF ( ( x .GT. r( 5) ) .AND. (x .LE. r( 6) ) ) THEN
       cd = calc_cd6_nrel5(y)
     ELSEIF ( ( x .GT. r( 6) ) .AND. (x .LE. r( 7) ) ) THEN
       cd = calc_cd7_nrel5(y)
     ELSEIF ( ( x .GT. r( 7) ) .AND. (x .LE. r( 8) ) ) THEN
       cd = calc_cd8_nrel5(y)  
  ENDIF
  
  IF ( x .LT. r(0) ) THEN
     PRINT*,'—————————————'
     PRINT*,'radius less than blade root in calc_cd_nrel5.pro'
     PRINT*,'radius: ', x
     PRINT*,'—————————————'
     STOP
  ENDIF
  
  IF ( x .GT. r(ns-1) ) THEN
     PRINT*,'—————————————'
     PRINT*,'radius greater than blade length in calc_cd_nrel5.pro'
     PRINT*,'radius: ', x
     PRINT*,'—————————————'
     STOP
  ENDIF
!_________________________________________________________________________
  
END FUNCTION calc_cd_nrel5
!_________________________________________________________________________

FUNCTION calc_cd1_nrel5( y ) RESULT( cd1 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: Cylinder1
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  REAL :: cd1

  cd1 = 0.5
!_________________________________________________________________________
  
END FUNCTION calc_cd1_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd2_nrel5( y ) RESULT( cd2 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: Cylinder2
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  REAL :: cd2

  cd2 = 0.35
!_________________________________________________________________________
  
END FUNCTION calc_cd2_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd3_nrel5( y ) RESULT( cd3 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU40_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 136

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd3
  INTEGER :: i

  cd3 = 0.0

!               alpha    CD   
  cd(0:1, 0  ) = [-180.0, 0.0602]
  cd(0:1, 1  ) = [-175.0, 0.0699]
  cd(0:1, 2  ) = [-170.0, 0.1107]
  cd(0:1, 3  ) = [-160.0, 0.3045]
  cd(0:1, 4  ) = [-155.0, 0.4179]
  cd(0:1, 5  ) = [-150.0, 0.5355]
  cd(0:1, 6  ) = [-145.0, 0.6535]
  cd(0:1, 7  ) = [-140.0, 0.7685]
  cd(0:1, 8  ) = [-135.0, 0.8777]
  cd(0:1, 9  ) = [-130.0, 0.9788]
  cd(0:1, 10 ) = [-125.0, 1.0700]
  cd(0:1, 11 ) = [-120.0, 1.1499]
  cd(0:1, 12 ) = [-115.0, 1.2174]
  cd(0:1, 13 ) = [-110.0, 1.2716]
  cd(0:1, 14 ) = [-105.0, 1.3118]
  cd(0:1, 15 ) = [-100.0, 1.3378]
  cd(0:1, 16 ) = [-95.0, 1.3492]
  cd(0:1, 17 ) = [-90.0, 1.3460]
  cd(0:1, 18 ) = [-85.0, 1.3283]
  cd(0:1, 19 ) = [-80.0, 1.2964]
  cd(0:1, 20 ) = [-75.0, 1.2507]
  cd(0:1, 21 ) = [-70.0, 1.1918]
  cd(0:1, 22 ) = [-65.0, 1.1204]
  cd(0:1, 23 ) = [-60.0, 1.0376]
  cd(0:1, 24 ) = [-55.0, 0.9446]
  cd(0:1, 25 ) = [-50.0, 0.8429]
  cd(0:1, 26 ) = [-45.0, 0.7345]
  cd(0:1, 27 ) = [-40.0, 0.6215]
  cd(0:1, 28 ) = [-35.0, 0.5067]
  cd(0:1, 29 ) = [-30.0, 0.3932]
  cd(0:1, 30 ) = [-25.0, 0.2849]
  cd(0:1, 31 ) = [-24.0, 0.2642]
  cd(0:1, 32 ) = [-23.0, 0.2440]
  cd(0:1, 33 ) = [-22.0, 0.2242]
  cd(0:1, 34 ) = [-21.0, 0.2049]
  cd(0:1, 35 ) = [-20.0, 0.1861]
  cd(0:1, 36 ) = [-19.0, 0.1687]
  cd(0:1, 37 ) = [-18.0, 0.1533]
  cd(0:1, 38 ) = [-17.0, 0.1398]
  cd(0:1, 39 ) = [-16.0, 0.1281]
  cd(0:1, 40 ) = [-15.0, 0.1183]
  cd(0:1, 41 ) = [-14.0, 0.1101]
  cd(0:1, 42 ) = [-13.0, 0.1036]
  cd(0:1, 43 ) = [-12.0, 0.0986]
  cd(0:1, 44 ) = [-11.0, 0.0951]
  cd(0:1, 45 ) = [-10.0, 0.0931]
  cd(0:1, 46 ) = [-8.0,  0.0930]
  cd(0:1, 47 ) = [-6.0,  0.0689]
  cd(0:1, 48 ) = [-5.5,  0.0614]
  cd(0:1, 49 ) = [-5.0,  0.0547]
  cd(0:1, 50 ) = [-4.5,  0.0480]
  cd(0:1, 51 ) = [-4.0,  0.0411]
  cd(0:1, 52 ) = [-3.5,  0.0349]
  cd(0:1, 53 ) = [-3.0, 0.0299]
  cd(0:1, 54 ) = [-2.5, 0.0255]
  cd(0:1, 55 ) = [-2.0, 0.0198]
  cd(0:1, 56 ) = [-1.5, 0.0164]
  cd(0:1, 57 ) = [-1.0, 0.0147]
  cd(0:1, 58 ) = [-0.5, 0.0137]
  cd(0:1, 59 ) = [0.0, 0.0113]
  cd(0:1, 60 ) = [0.5, 0.0114]
  cd(0:1, 61 ) = [1.0, 0.0118]
  cd(0:1, 62 ) = [1.5, 0.0122]
  cd(0:1, 63 ) = [2.0, 0.0124]
  cd(0:1, 64 ) = [2.5, 0.0124]
  cd(0:1, 65 ) = [3.0, 0.0123]
  cd(0:1, 66 ) = [3.5, 0.0120]
  cd(0:1, 67 ) = [4.0, 0.0119]
  cd(0:1, 68 ) = [4.5, 0.0122]
  cd(0:1, 69 ) = [5.0, 0.0125]
  cd(0:1, 70 ) = [5.5, 0.0129]
  cd(0:1, 71 ) = [6.0, 0.0135]
  cd(0:1, 72 ) = [6.5, 0.0144]
  cd(0:1, 73 ) = [7.0, 0.0158]
  cd(0:1, 74 ) = [7.5, 0.0174]
  cd(0:1, 75 ) = [8.0, 0.0198]
  cd(0:1, 76 ) = [8.5, 0.0231]
  cd(0:1, 77 ) = [9.0, 0.0275]
  cd(0:1, 78 ) = [9.5, 0.0323]
  cd(0:1, 79 ) = [10.0, 0.0393]
  cd(0:1, 80 ) = [10.5, 0.0475]
  cd(0:1, 81 ) = [11.0, 0.0580]
  cd(0:1, 82 ) = [11.5, 0.0691]
  cd(0:1, 83 ) = [12.0, 0.0816]
  cd(0:1, 84 ) = [12.5, 0.0973]
  cd(0:1, 85 ) = [13.0, 0.1129]
  cd(0:1, 86 ) = [13.5, 0.1288]
  cd(0:1, 87 ) = [14.5, 0.1650]
  cd(0:1, 88 ) = [15.0, 0.1845]
  cd(0:1, 89 ) = [15.5, 0.2052]
  cd(0:1, 90 ) = [16.0, 0.2250]
  cd(0:1, 91 ) = [16.5, 0.2467]
  cd(0:1, 92 ) = [17.0, 0.2684]
  cd(0:1, 93 ) = [17.5, 0.2900]
  cd(0:1, 94 ) = [18.0, 0.3121]
  cd(0:1, 95 ) = [19.0, 0.3554]
  cd(0:1, 96 ) = [19.5, 0.3783]
  cd(0:1, 97 ) = [20.5, 0.4212]
  cd(0:1, 98 ) = [21.0, 0.4415]
  cd(0:1, 99 ) = [22.0, 0.4830]
  cd(0:1, 100) = [23.0, 0.5257]
  cd(0:1, 101) = [24.0, 0.5694]
  cd(0:1, 102) = [25.0, 0.6141]
  cd(0:1, 103) = [26.0, 0.6593]
  cd(0:1, 104) = [28.0, 0.7513]
  cd(0:1, 105) = [30.0, 0.8441]
  cd(0:1, 106) = [32.0, 0.9364]
  cd(0:1, 107) = [35.0, 1.0722]
  cd(0:1, 108) = [40.0, 1.2873]
  cd(0:1, 109) = [45.0, 1.4796]
  cd(0:1, 110) = [50.0, 1.6401]
  cd(0:1, 111) = [55.0, 1.7609]
  cd(0:1, 112) = [60.0, 1.8360]
  cd(0:1, 113) = [65.0, 1.8614]
  cd(0:1, 114) = [70.0, 1.8347]
  cd(0:1, 115) = [75.0, 1.7567]
  cd(0:1, 116) = [80.0, 1.6334]
  cd(0:1, 117) = [85.0, 1.4847]
  cd(0:1, 118) = [90.0, 1.3879]
  cd(0:1, 119) = [95.0, 1.3912]
  cd(0:1, 120) = [100.0, 1.3795]
  cd(0:1, 121) = [105.0, 1.3528]
  cd(0:1, 122) = [110.0, 1.3114]
  cd(0:1, 123) = [115.0, 1.2557]
  cd(0:1, 124) = [120.0, 1.1864]
  cd(0:1, 125) = [125.0, 1.1041]
  cd(0:1, 126) = [130.0, 1.0102]
  cd(0:1, 127) = [135.0, 0.9060]
  cd(0:1, 128) = [140.0, 0.7935]
  cd(0:1, 129) = [145.0, 0.6750]
  cd(0:1, 130) = [150.0, 0.5532]
  cd(0:1, 131) = [155.0, 0.4318]
  cd(0:1, 132) = [160.0, 0.3147]
  cd(0:1, 133) = [170.0, 0.1144]
  cd(0:1, 134) = [175.0, 0.0702]
  cd(0:1, 135) = [180.0, 0.0602]
  
  IF ( y .EQ. cd(0,0) ) cd3 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd3 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd3 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd3_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd3 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd3_nrel5.pro'
  ENDIF
  
!print,y,cd3
  
!_________________________________________________________________________
  
END FUNCTION calc_cd3_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd4_nrel5( y ) RESULT( cd4 )
!_________________________________________________________________________
!
! Data from NREL 5 MW DU35_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 135

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd4
  INTEGER :: i

  cd4 = 0.0

!                 alpha    CD   
  cd(0:1, 0  ) = [-180.0, 0.0407]
  cd(0:1, 1  ) = [-175.0, 0.0507]
  cd(0:1, 2  ) = [-170.0, 0.1055]
  cd(0:1, 3  ) = [-160.0, 0.2982]
  cd(0:1, 4  ) = [-155.0, 0.4121]
  cd(0:1, 5  ) = [-150.0, 0.5308]
  cd(0:1, 6  ) = [-145.0, 0.6503]
  cd(0:1, 7  ) = [-140.0, 0.7672]
  cd(0:1, 8  ) = [-135.0, 0.8785]
  cd(0:1, 9  ) = [-130.0, 0.9819]
  cd(0:1, 10 ) = [-125.0, 1.0756]
  cd(0:1, 11 ) = [-120.0, 1.1580]
  cd(0:1, 12 ) = [-115.0, 1.2280]
  cd(0:1, 13 ) = [-110.0, 1.2847]
  cd(0:1, 14 ) = [-105.0, 1.3274]
  cd(0:1, 15 ) = [-100.0, 1.3557]
  cd(0:1, 16 ) = [-95.0, 1.3692]
  cd(0:1, 17 ) = [-90.0, 1.3680]
  cd(0:1, 18 ) = [-85.0, 1.3521]
  cd(0:1, 19 ) = [-80.0, 1.3218]
  cd(0:1, 20 ) = [-75.0, 1.2773]
  cd(0:1, 21 ) = [-70.0, 1.2193]
  cd(0:1, 22 ) = [-65.0, 1.1486]
  cd(0:1, 23 ) = [-60.0, 1.0660]
  cd(0:1, 24 ) = [-55.0, 0.9728]
  cd(0:1, 25 ) = [-50.0, 0.8705]
  cd(0:1, 26 ) = [-45.0, 0.7611]
  cd(0:1, 27 ) = [-40.0, 0.6466]
  cd(0:1, 28 ) = [-35.0, 0.5299]
  cd(0:1, 29 ) = [-30.0, 0.4141]
  cd(0:1, 30 ) = [-25.0, 0.3030]
  cd(0:1, 31 ) = [-24.0, 0.2817]
  cd(0:1, 32 ) = [-23.0, 0.2608]
  cd(0:1, 33 ) = [-22.0, 0.2404]
  cd(0:1, 34 ) = [-21.0, 0.2205]
  cd(0:1, 35 ) = [-20.0, 0.2011]
  cd(0:1, 36 ) = [-19.0, 0.1822]
  cd(0:1, 37 ) = [-18.0, 0.1640]
  cd(0:1, 38 ) = [-17.0, 0.1465]
  cd(0:1, 39 ) = [-16.0, 0.1300]
  cd(0:1, 40 ) = [-15.0, 0.1145]
  cd(0:1, 41 ) = [-14.0, 0.1000]
  cd(0:1, 42 ) = [-13.0, 0.0867]
  cd(0:1, 43 ) = [-12.0, 0.0744]
  cd(0:1, 44 ) = [-11.0, 0.0633]
  cd(0:1, 45 ) = [-10.0, 0.0534]
  cd(0:1, 46 ) = [-5.5, 0.0245]
  cd(0:1, 47 ) = [-5.0, 0.0225]
  cd(0:1, 48 ) = [-4.5, 0.0196]
  cd(0:1, 49 ) = [-4.0, 0.0174]
  cd(0:1, 50 ) = [-3.5, 0.0162]
  cd(0:1, 51 ) = [-3.0, 0.0144]
  cd(0:1, 52 ) = [-3.0, 0.0240]
  cd(0:1, 53 ) = [-2.5, 0.0188]
  cd(0:1, 54 ) = [-2.0, 0.0160]
  cd(0:1, 55 ) = [-1.5, 0.0137]
  cd(0:1, 56 ) = [-1.0, 0.0118]
  cd(0:1, 57 ) = [-0.5, 0.0104]
  cd(0:1, 58 ) = [0.0, 0.0094]
  cd(0:1, 59 ) = [0.5, 0.0096]
  cd(0:1, 60 ) = [1.0, 0.0098]
  cd(0:1, 61 ) = [1.5, 0.0099]
  cd(0:1, 62 ) = [2.0, 0.0100]
  cd(0:1, 63 ) = [2.5, 0.0102]
  cd(0:1, 64 ) = [3.0, 0.0103]
  cd(0:1, 65 ) = [3.5, 0.0104]
  cd(0:1, 66 ) = [4.0, 0.0105]
  cd(0:1, 67 ) = [4.5, 0.0107]
  cd(0:1, 68 ) = [5.0, 0.0108]
  cd(0:1, 69 ) = [5.5, 0.0109]
  cd(0:1, 70 ) = [6.0, 0.0110]
  cd(0:1, 71 ) = [6.5, 0.0113]
  cd(0:1, 72 ) = [7.0, 0.0115]
  cd(0:1, 73 ) = [7.5, 0.0117]
  cd(0:1, 74 ) = [8.0, 0.0120]
  cd(0:1, 75 ) = [8.5, 0.0126]
  cd(0:1, 76 ) = [9.0, 0.0133]
  cd(0:1, 77 ) = [9.5, 0.0143]
  cd(0:1, 78 ) = [10.0, 0.0156]
  cd(0:1, 79 ) = [10.5, 0.0174]
  cd(0:1, 80 ) = [11.0, 0.0194]
  cd(0:1, 81 ) = [11.5, 0.0227]
  cd(0:1, 82 ) = [12.0, 0.0269]
  cd(0:1, 83 ) = [12.5, 0.0319]
  cd(0:1, 84 ) = [13.0, 0.0398]
  cd(0:1, 85 ) = [13.5, 0.0488]
  cd(0:1, 86 ) = [14.0, 0.0614]
  cd(0:1, 87 ) = [14.5, 0.0786]
  cd(0:1, 88 ) = [15.5, 0.1173]
  cd(0:1, 89 ) = [16.0, 0.1377]
  cd(0:1, 90 ) = [16.5, 0.1600]
  cd(0:1, 91 ) = [17.0, 0.1814]
  cd(0:1, 92 ) = [17.5, 0.2042]
  cd(0:1, 93 ) = [18.0, 0.2316]
  cd(0:1, 94 ) = [19.0, 0.2719]
  cd(0:1, 95 ) = [19.5, 0.2906]
  cd(0:1, 96 ) = [20.0, 0.3085]
  cd(0:1, 97 ) = [21.0, 0.3447]
  cd(0:1, 98 ) = [22.0, 0.3820]
  cd(0:1, 99 ) = [23.0, 0.4203]
  cd(0:1, 100) = [24.0, 0.4593]
  cd(0:1, 101) = [25.0, 0.4988]
  cd(0:1, 102) = [26.0, 0.5387]
  cd(0:1, 103) = [28.0, 0.6187]
  cd(0:1, 104) = [30.0, 0.6978]
  cd(0:1, 105) = [32.0, 0.7747]
  cd(0:1, 106) = [35.0, 0.8869]
  cd(0:1, 107) = [40.0, 1.0671]
  cd(0:1, 108) = [45.0, 1.2319]
  cd(0:1, 109) = [50.0, 1.3747]
  cd(0:1, 110) = [55.0, 1.4899]
  cd(0:1, 111) = [60.0, 1.5728]
  cd(0:1, 112) = [65.0, 1.6202]
  cd(0:1, 113) = [70.0, 1.6302]
  cd(0:1, 114) = [75.0, 1.6031]
  cd(0:1, 115) = [80.0, 1.5423]
  cd(0:1, 116) = [85.0, 1.4598]
  cd(0:1, 117) = [90.0, 1.4041]
  cd(0:1, 118) = [95.0, 1.4053]
  cd(0:1, 119) = [100.0, 1.3914]
  cd(0:1, 120) = [105.0, 1.3625]
  cd(0:1, 121) = [110.0, 1.3188]
  cd(0:1, 122) = [115.0, 1.2608]
  cd(0:1, 123) = [120.0, 1.1891]
  cd(0:1, 124) = [125.0, 1.1046]
  cd(0:1, 125) = [130.0, 1.0086]
  cd(0:1, 126) = [135.0, 0.9025]
  cd(0:1, 127) = [140.0, 0.7883]
  cd(0:1, 128) = [145.0, 0.6684]
  cd(0:1, 129) = [150.0, 0.5457]
  cd(0:1, 130) = [155.0, 0.4236]
  cd(0:1, 131) = [160.0, 0.3066]
  cd(0:1, 132) = [170.0, 0.1085]
  cd(0:1, 133) = [175.0, 0.0510]
  cd(0:1, 134) = [180.0, 0.0407]
  
  IF ( y .EQ. cd(0,0) ) cd4 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd4 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd4 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd4_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd4 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd4_nrel5.pro'
  ENDIF
  
!print,y,cd4
  
!_________________________________________________________________________
  
END FUNCTION calc_cd4_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd5_nrel5( y ) RESULT( cd5 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU30_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 143

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd5
  INTEGER :: i

  cd5 = 0.0

!                alpha    CD  
  cd(0:1, 0  ) = [-180.0, 0.0267]
  cd(0:1, 1  ) = [-175.0, 0.0370]
  cd(0:1, 2  ) = [-170.0, 0.0968]
  cd(0:1, 3  ) = [-160.0, 0.2876]
  cd(0:1, 4  ) = [-155.0, 0.4025]
  cd(0:1, 5  ) = [-150.0, 0.5232]
  cd(0:1, 6  ) = [-145.0, 0.6454]
  cd(0:1, 7  ) = [-140.0, 0.7656]
  cd(0:1, 8  ) = [-135.0, 0.8807]
  cd(0:1, 9  ) = [-130.0, 0.9882]
  cd(0:1, 10 ) = [-125.0, 1.0861]
  cd(0:1, 11 ) = [-120.0, 1.1730]
  cd(0:1, 12 ) = [-115.0, 1.2474]
  cd(0:1, 13 ) = [-110.0, 1.3084]
  cd(0:1, 14 ) = [-105.0, 1.3552]
  cd(0:1, 15 ) = [-100.0, 1.3875]
  cd(0:1, 16 ) = [-95.0, 1.4048]
  cd(0:1, 17 ) = [-90.0, 1.4070]
  cd(0:1, 18 ) = [-85.0, 1.3941]
  cd(0:1, 19 ) = [-80.0, 1.3664]
  cd(0:1, 20 ) = [-75.0, 1.3240]
  cd(0:1, 21 ) = [-70.0, 1.2676]
  cd(0:1, 22 ) = [-65.0, 1.1978]
  cd(0:1, 23 ) = [-60.0, 1.1156]
  cd(0:1, 24 ) = [-55.0, 1.0220]
  cd(0:1, 25 ) = [-50.0, 0.9187]
  cd(0:1, 26 ) = [-45.0, 0.8074]
  cd(0:1, 27 ) = [-40.0, 0.6904]
  cd(0:1, 28 ) = [-35.0, 0.5703]
  cd(0:1, 29 ) = [-30.0, 0.4503]
  cd(0:1, 30 ) = [-25.0, 0.3357]
  cd(0:1, 31 ) = [-24.0, 0.3147]
  cd(0:1, 32 ) = [-23.0, 0.2946]
  cd(0:1, 33 ) = [-22.0, 0.2752]
  cd(0:1, 34 ) = [-21.0, 0.2566]
  cd(0:1, 35 ) = [-20.0, 0.2388]
  cd(0:1, 36 ) = [-19.0, 0.2218]
  cd(0:1, 37 ) = [-18.0, 0.2056]
  cd(0:1, 38 ) = [-17.0, 0.1901]
  cd(0:1, 39 ) = [-16.0, 0.1754]
  cd(0:1, 40 ) = [-15.3, 0.1649]
  cd(0:1, 41 ) = [-14.2, 0.1461]
  cd(0:1, 42 ) = [-13.2, 0.1263]
  cd(0:1, 43 ) = [-12.2, 0.1051]
  cd(0:1, 44 ) = [-11.2, 0.0886]
  cd(0:1, 45 ) = [-10.2, 0.0740]
  cd(0:1, 46 ) = [-9.7, 0.0684]
  cd(0:1, 47 ) = [-9.2, 0.0605]
  cd(0:1, 48 ) = [-8.2, 0.0270]
  cd(0:1, 49 ) = [-7.2, 0.0180]
  cd(0:1, 50 ) = [-6.7, 0.0166]
  cd(0:1, 51 ) = [-6.1, 0.0152]
  cd(0:1, 52 ) = [-6.0, 0.0117]
  cd(0:1, 53 ) = [-5.5, 0.0105]
  cd(0:1, 54 ) = [-5.0, 0.0097]
  cd(0:1, 55 ) = [-4.5, 0.0092]
  cd(0:1, 56 ) = [-4.0, 0.0091]
  cd(0:1, 57 ) = [-3.5, 0.0089]
  cd(0:1, 58 ) = [-3.0, 0.0089]
  cd(0:1, 59 ) = [-2.5, 0.0088]
  cd(0:1, 60 ) = [-2.0, 0.0088]
  cd(0:1, 61 ) = [-1.5, 0.0088]
  cd(0:1, 62 ) = [-1.0, 0.0088]
  cd(0:1, 63 ) = [-0.5, 0.0088]
  cd(0:1, 64 ) = [0.0, 0.0087]
  cd(0:1, 65 ) = [0.5, 0.0087]
  cd(0:1, 66 ) = [1.0, 0.0088]
  cd(0:1, 67 ) = [1.5, 0.0089]
  cd(0:1, 68 ) = [2.0, 0.0090]
  cd(0:1, 69 ) = [2.5, 0.0091]
  cd(0:1, 70 ) = [3.0, 0.0092]
  cd(0:1, 71 ) = [3.5, 0.0093]
  cd(0:1, 72 ) = [4.0, 0.0095]
  cd(0:1, 73 ) = [4.5, 0.0096]
  cd(0:1, 74 ) = [5.0, 0.0097]
  cd(0:1, 75 ) = [5.5, 0.0099]
  cd(0:1, 76 ) = [6.0, 0.0101]
  cd(0:1, 77 ) = [6.5, 0.0103]
  cd(0:1, 78 ) = [7.0, 0.0107]
  cd(0:1, 79 ) = [7.5, 0.0112]
  cd(0:1, 80 ) = [8.0, 0.0125]
  cd(0:1, 81 ) = [9.0, 0.0155]
  cd(0:1, 82 ) = [9.5, 0.0171]
  cd(0:1, 83 ) = [10.0, 0.0192]
  cd(0:1, 84 ) = [10.5, 0.0219]
  cd(0:1, 85 ) = [11.0, 0.0255]
  cd(0:1, 86 ) = [11.5, 0.0307]
  cd(0:1, 87 ) = [12.0, 0.0370]
  cd(0:1, 88 ) = [12.5, 0.0452]
  cd(0:1, 89 ) = [13.0, 0.0630]
  cd(0:1, 90 ) = [13.5, 0.0784]
  cd(0:1, 91 ) = [14.0, 0.0931]
  cd(0:1, 92 ) = [14.5, 0.1081]
  cd(0:1, 93 ) = [15.0, 0.1239]
  cd(0:1, 94 ) = [15.5, 0.1415]
  cd(0:1, 95 ) = [16.0, 0.1592]
  cd(0:1, 96 ) = [16.5, 0.1743]
  cd(0:1, 97 ) = [17.0, 0.1903]
  cd(0:1, 98 ) = [17.5, 0.2044]
  cd(0:1, 99 ) = [18.0, 0.2186]
  cd(0:1, 100) = [18.5, 0.2324]
  cd(0:1, 101) = [19.0, 0.2455]
  cd(0:1, 102) = [19.5, 0.2584]
  cd(0:1, 103) = [20.0, 0.2689]
  cd(0:1, 104) = [20.5, 0.2814]
  cd(0:1, 105) = [21.0, 0.2943]
  cd(0:1, 106) = [22.0, 0.3246]
  cd(0:1, 107) = [23.0, 0.3557]
  cd(0:1, 108) = [24.0, 0.3875]
  cd(0:1, 109) = [25.0, 0.4198]
  cd(0:1, 110) = [26.0, 0.4524]
  cd(0:1, 111) = [28.0, 0.5183]
  cd(0:1, 112) = [30.0, 0.5843]
  cd(0:1, 113) = [32.0, 0.6492]
  cd(0:1, 114) = [35.0, 0.7438]
  cd(0:1, 115) = [40.0, 0.8970]
  cd(0:1, 116) = [45.0, 1.0402]
  cd(0:1, 117) = [50.0, 1.1686]
  cd(0:1, 118) = [55.0, 1.2779]
  cd(0:1, 119) = [60.0, 1.3647]
  cd(0:1, 120) = [65.0, 1.4267]
  cd(0:1, 121) = [70.0, 1.4621]
  cd(0:1, 122) = [75.0, 1.4708]
  cd(0:1, 123) = [80.0, 1.4544]
  cd(0:1, 124) = [85.0, 1.4196]
  cd(0:1, 125) = [90.0, 1.3938]
  cd(0:1, 126) = [95.0, 1.3943]
  cd(0:1, 127) = [100.0, 1.3798]
  cd(0:1, 128) = [105.0, 1.3504]
  cd(0:1, 129) = [110.0, 1.3063]
  cd(0:1, 130) = [115.0, 1.2481]
  cd(0:1, 131) = [120.0, 1.1763]
  cd(0:1, 132) = [125.0, 1.0919]
  cd(0:1, 133) = [130.0, 0.9962]
  cd(0:1, 134) = [135.0, 0.8906]
  cd(0:1, 135) = [140.0, 0.7771]
  cd(0:1, 136) = [145.0, 0.6581]
  cd(0:1, 137) = [150.0, 0.5364]
  cd(0:1, 138) = [155.0, 0.4157]
  cd(0:1, 139) = [160.0, 0.3000]
  cd(0:1, 140) = [170.0, 0.1051]
  cd(0:1, 141) = [175.0, 0.0388]
  cd(0:1, 142) = [180.0, 0.0267]

  IF ( y .EQ. cd(0,0) ) cd5 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd5 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd5 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd5_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd5 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd5_nrel5.pro'
  ENDIF
  
!print,y,cd5
  
!_________________________________________________________________________
  
END FUNCTION calc_cd5_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd6_nrel5( y ) RESULT( cd6 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: DU25_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 140

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd6
  INTEGER :: i

  cd6 = 0.0

!                alpha    CD  
  cd(0:1, 0  ) = [-180.0, 0.0202]
  cd(0:1, 1  ) = [-175.0, 0.0324]
  cd(0:1, 2  ) = [-170.0, 0.0943]
  cd(0:1, 3  ) = [-160.0, 0.2848]
  cd(0:1, 4  ) = [-155.0, 0.4001]
  cd(0:1, 5  ) = [-150.0, 0.5215]
  cd(0:1, 6  ) = [-145.0, 0.6447]
  cd(0:1, 7  ) = [-140.0, 0.7660]
  cd(0:1, 8  ) = [-135.0, 0.8823]
  cd(0:1, 9  ) = [-130.0, 0.9911]
  cd(0:1, 10 ) = [-125.0, 1.0905]
  cd(0:1, 11 ) = [-120.0, 1.1787]
  cd(0:1, 12 ) = [-115.0, 1.2545]
  cd(0:1, 13 ) = [-110.0, 1.3168]
  cd(0:1, 14 ) = [-105.0, 1.3650]
  cd(0:1, 15 ) = [-100.0, 1.3984]
  cd(0:1, 16 ) = [-95.0, 1.4169]
  cd(0:1, 17 ) = [-90.0, 1.4201]
  cd(0:1, 18 ) = [-85.0, 1.4081]
  cd(0:1, 19 ) = [-80.0, 1.3811]
  cd(0:1, 20 ) = [-75.0, 1.3394]
  cd(0:1, 21 ) = [-70.0, 1.2833]
  cd(0:1, 22 ) = [-65.0, 1.2138]
  cd(0:1, 23 ) = [-60.0, 1.1315]
  cd(0:1, 24 ) = [-55.0, 1.0378]
  cd(0:1, 25 ) = [-50.0, 0.9341]
  cd(0:1, 26 ) = [-45.0, 0.8221]
  cd(0:1, 27 ) = [-40.0, 0.7042]
  cd(0:1, 28 ) = [-35.0, 0.5829]
  cd(0:1, 29 ) = [-30.0, 0.4616]
  cd(0:1, 30 ) = [-25.0, 0.3441]
  cd(0:1, 31 ) = [-24.0, 0.3209]
  cd(0:1, 32 ) = [-23.0, 0.2972]
  cd(0:1, 33 ) = [-22.0, 0.2730]
  cd(0:1, 34 ) = [-21.0, 0.2485]
  cd(0:1, 35 ) = [-20.0, 0.2237]
  cd(0:1, 36 ) = [-19.0, 0.1990]
  cd(0:1, 37 ) = [-18.0, 0.1743]
  cd(0:1, 38 ) = [-17.0, 0.1498]
  cd(0:1, 39 ) = [-16.0, 0.1256]
  cd(0:1, 40 ) = [-15.0, 0.1020]
  cd(0:1, 41 ) = [-14.0, 0.0789]
  cd(0:1, 42 ) = [-13.0, 0.0567]
  cd(0:1, 43 ) = [-12.0, 0.0271]
  cd(0:1, 44 ) = [-11.0, 0.0303]
  cd(0:1, 45 ) = [-10.0, 0.0287]
  cd(0:1, 46 ) = [-9.0, 0.0271]
  cd(0:1, 47 ) = [-8.5, 0.0264]
  cd(0:1, 48 ) = [-7.5, 0.0114]
  cd(0:1, 49 ) = [-6.4, 0.0094]
  cd(0:1, 50 ) = [-5.4, 0.0086]
  cd(0:1, 51 ) = [-5.0, 0.0073]
  cd(0:1, 52 ) = [-4.5, 0.0071]
  cd(0:1, 53 ) = [-4.0, 0.0070]
  cd(0:1, 54 ) = [-3.5, 0.0069]
  cd(0:1, 55 ) = [-3.0, 0.0068]
  cd(0:1, 56 ) = [-2.5, 0.0068]
  cd(0:1, 57 ) = [-2.0, 0.0068]
  cd(0:1, 58 ) = [-1.5, 0.0067]
  cd(0:1, 59 ) = [-1.0, 0.0067]
  cd(0:1, 60 ) = [-0.5, 0.0067]
  cd(0:1, 61 ) = [0.0, 0.0065]
  cd(0:1, 62 ) = [0.5, 0.0065]
  cd(0:1, 63 ) = [1.0, 0.0066]
  cd(0:1, 64 ) = [1.5, 0.0067]
  cd(0:1, 65 ) = [2.0, 0.0068]
  cd(0:1, 66 ) = [2.5, 0.0069]
  cd(0:1, 67 ) = [3.0, 0.0070]
  cd(0:1, 68 ) = [3.5, 0.0071]
  cd(0:1, 69 ) = [4.0, 0.0073]
  cd(0:1, 70 ) = [4.5, 0.0076]
  cd(0:1, 71 ) = [5.0, 0.0079]
  cd(0:1, 72 ) = [6.0, 0.0099]
  cd(0:1, 73 ) = [6.5, 0.0117]
  cd(0:1, 74 ) = [7.0, 0.0132]
  cd(0:1, 75 ) = [7.5, 0.0143]
  cd(0:1, 76 ) = [8.0, 0.0153]
  cd(0:1, 77 ) = [8.5, 0.0165]
  cd(0:1, 78 ) = [9.0, 0.0181]
  cd(0:1, 79 ) = [9.5, 0.0211]
  cd(0:1, 80 ) = [10.0, 0.0262]
  cd(0:1, 81 ) = [10.5, 0.0336]
  cd(0:1, 82 ) = [11.0, 0.0420]
  cd(0:1, 83 ) = [11.5, 0.0515]
  cd(0:1, 84 ) = [12.0, 0.0601]
  cd(0:1, 85 ) = [12.5, 0.0693]
  cd(0:1, 86 ) = [13.0, 0.0785]
  cd(0:1, 87 ) = [13.5, 0.0888]
  cd(0:1, 88 ) = [14.0, 0.1000]
  cd(0:1, 89 ) = [14.5, 0.1108]
  cd(0:1, 90 ) = [15.0, 0.1219]
  cd(0:1, 91 ) = [15.5, 0.1325]
  cd(0:1, 92 ) = [16.0, 0.1433]
  cd(0:1, 93 ) = [16.5, 0.1541]
  cd(0:1, 94 ) = [17.0, 0.1649]
  cd(0:1, 95 ) = [17.5, 0.1754]
  cd(0:1, 96 ) = [18.0, 0.1845]
  cd(0:1, 97 ) = [18.5, 0.1953]
  cd(0:1, 98 ) = [19.0, 0.2061]
  cd(0:1, 99 ) = [19.5, 0.2170]
  cd(0:1, 100) = [20.0, 0.2280]
  cd(0:1, 101) = [20.5, 0.2390]
  cd(0:1, 102) = [21.0, 0.2536]
  cd(0:1, 103) = [22.0, 0.2814]
  cd(0:1, 104) = [23.0, 0.3098]
  cd(0:1, 105) = [24.0, 0.3386]
  cd(0:1, 106) = [25.0, 0.3678]
  cd(0:1, 107) = [26.0, 0.3972]
  cd(0:1, 108) = [28.0, 0.4563]
  cd(0:1, 109) = [30.0, 0.5149]
  cd(0:1, 110) = [32.0, 0.5720]
  cd(0:1, 111) = [35.0, 0.6548]
  cd(0:1, 112) = [40.0, 0.7901]
  cd(0:1, 113) = [45.0, 0.9190]
  cd(0:1, 114) = [50.0, 1.0378]
  cd(0:1, 115) = [55.0, 1.1434]
  cd(0:1, 116) = [60.0, 1.2333]
  cd(0:1, 117) = [65.0, 1.3055]
  cd(0:1, 118) = [70.0, 1.3587]
  cd(0:1, 119) = [75.0, 1.3922]
  cd(0:1, 120) = [80.0, 1.4063]
  cd(0:1, 121) = [85.0, 1.4042]
  cd(0:1, 122) = [90.0, 1.3985]
  cd(0:1, 123) = [95.0, 1.3973]
  cd(0:1, 124) = [100.0, 1.3810]
  cd(0:1, 125) = [105.0, 1.3498]
  cd(0:1, 126) = [110.0, 1.3041]
  cd(0:1, 127) = [115.0, 1.2442]
  cd(0:1, 128) = [120.0, 1.1709]
  cd(0:1, 129) = [125.0, 1.0852]
  cd(0:1, 130) = [130.0, 0.9883]
  cd(0:1, 131) = [135.0, 0.8818]
  cd(0:1, 132) = [140.0, 0.7676]
  cd(0:1, 133) = [145.0, 0.6481]
  cd(0:1, 134) = [150.0, 0.5264]
  cd(0:1, 135) = [155.0, 0.4060]
  cd(0:1, 136) = [160.0, 0.2912]
  cd(0:1, 137) = [170.0, 0.0995]
  cd(0:1, 138) = [175.0, 0.0356]
  cd(0:1, 139) = [180.0, 0.0202]
  
  IF ( y .EQ. cd(0,0) ) cd6 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd6 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd6 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd6_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd6 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd6_nrel5.pro'
  ENDIF
  
!print,y,cd6
  
!_________________________________________________________________________
  
END FUNCTION calc_cd6_nrel5
!_________________________________________________________________________

!=========================================================================

FUNCTION calc_cd7_nrel5( y ) RESULT( cd7 )
!_________________________________________________________________________
!
! Data from  NREL 5 MW: DU21_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 142

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd7
  INTEGER :: i

  cd7 = 0.0

!                alpha    CD  
  cd(0:1, 0  ) = [-180.0, 0.0185]
  cd(0:1, 1  ) = [-175.0, 0.0332]
  cd(0:1, 2  ) = [-170.0, 0.0945]
  cd(0:1, 3  ) = [-160.0, 0.2809]
  cd(0:1, 4  ) = [-155.0, 0.3932]
  cd(0:1, 5  ) = [-150.0, 0.5112]
  cd(0:1, 6  ) = [-145.0, 0.6309]
  cd(0:1, 7  ) = [-140.0, 0.7485]
  cd(0:1, 8  ) = [-135.0, 0.8612]
  cd(0:1, 9  ) = [-130.0, 0.9665]
  cd(0:1, 10 ) = [-125.0, 1.0625]
  cd(0:1, 11 ) = [-120.0, 1.1476]
  cd(0:1, 12 ) = [-115.0, 1.2206]
  cd(0:1, 13 ) = [-110.0, 1.2805]
  cd(0:1, 14 ) = [-105.0, 1.3265]
  cd(0:1, 15 ) = [-100.0, 1.3582]
  cd(0:1, 16 ) = [-95.0, 1.3752]
  cd(0:1, 17 ) = [-90.0, 1.3774]
  cd(0:1, 18 ) = [-85.0, 1.3648]
  cd(0:1, 19 ) = [-80.0, 1.3376]
  cd(0:1, 20 ) = [-75.0, 1.2962]
  cd(0:1, 21 ) = [-70.0, 1.2409]
  cd(0:1, 22 ) = [-65.0, 1.1725]
  cd(0:1, 23 ) = [-60.0, 1.0919]
  cd(0:1, 24 ) = [-55.0, 1.0002]
  cd(0:1, 25 ) = [-50.0, 0.8990]
  cd(0:1, 26 ) = [-45.0, 0.7900]
  cd(0:1, 27 ) = [-40.0, 0.6754]
  cd(0:1, 28 ) = [-35.0, 0.5579]
  cd(0:1, 29 ) = [-30.0, 0.4405]
  cd(0:1, 30 ) = [-25.0, 0.3256]
  cd(0:1, 31 ) = [-24.0, 0.3013]
  cd(0:1, 32 ) = [-23.0, 0.2762]
  cd(0:1, 33 ) = [-22.0, 0.2506]
  cd(0:1, 34 ) = [-21.0, 0.2246]
  cd(0:1, 35 ) = [-20.0, 0.1983]
  cd(0:1, 36 ) = [-19.0, 0.1720]
  cd(0:1, 37 ) = [-18.0, 0.1457]
  cd(0:1, 38 ) = [-17.0, 0.1197]
  cd(0:1, 39 ) = [-16.0, 0.0940]
  cd(0:1, 40 ) = [-15.0, 0.0689]
  cd(0:1, 41 ) = [-14.5, 0.0567]
  cd(0:1, 42 ) = [-12.0, 0.0271]
  cd(0:1, 43 ) = [-11.0, 0.0303]
  cd(0:1, 44 ) = [-10.0, 0.0287]
  cd(0:1, 45 ) = [-8.1,  0.0124]
  cd(0:1, 46 ) = [-7.6,  0.0109]
  cd(0:1, 47 ) = [-7.1,  0.0092]
  cd(0:1, 48 ) = [-6.6,  0.0083]
  cd(0:1, 49 ) = [-6.5,  0.0089]
  cd(0:1, 50 ) = [-6.0,  0.0082]
  cd(0:1, 51 ) = [-5.5,  0.0074]
  cd(0:1, 52 ) = [-5.0,  0.0069]
  cd(0:1, 53 ) = [-4.5,  0.0065]
  cd(0:1, 54 ) = [-4.0, 0.0063]
  cd(0:1, 55 ) = [-3.5, 0.0061]
  cd(0:1, 56 ) = [-3.0, 0.0058]
  cd(0:1, 57 ) = [-2.5, 0.0057]
  cd(0:1, 58 ) = [-2.0, 0.0057]
  cd(0:1, 59 ) = [-1.5, 0.0057]
  cd(0:1, 60 ) = [-1.0, 0.0057]
  cd(0:1, 61 ) = [-0.5, 0.0057]
  cd(0:1, 62 ) = [0.0, 0.0057]
  cd(0:1, 63 ) = [0.5, 0.0057]
  cd(0:1, 64 ) = [1.0, 0.0058]
  cd(0:1, 65 ) = [1.5, 0.0058]
  cd(0:1, 66 ) = [2.0, 0.0059]
  cd(0:1, 67 ) = [2.5, 0.0061]
  cd(0:1, 68 ) = [3.0, 0.0063]
  cd(0:1, 69 ) = [3.5, 0.0066]
  cd(0:1, 70 ) = [4.0, 0.0071]
  cd(0:1, 71 ) = [4.5, 0.0079]
  cd(0:1, 72 ) = [5.0, 0.0090]
  cd(0:1, 73 ) = [5.5, 0.0103]
  cd(0:1, 74 ) = [6.0, 0.0113]
  cd(0:1, 75 ) = [6.5, 0.0122]
  cd(0:1, 76 ) = [7.0, 0.0131]
  cd(0:1, 77 ) = [7.5, 0.0139]
  cd(0:1, 78 ) = [8.0, 0.0147]
  cd(0:1, 79 ) = [8.5, 0.0158]
  cd(0:1, 80 ) = [9.0, 0.0181]
  cd(0:1, 81 ) = [9.5, 0.0211]
  cd(0:1, 82 ) = [10.0, 0.0255]
  cd(0:1, 83 ) = [10.5, 0.0301]
  cd(0:1, 84 ) = [11.0, 0.0347]
  cd(0:1, 85 ) = [11.5, 0.0401]
  cd(0:1, 86 ) = [12.0, 0.0468]
  cd(0:1, 87 ) = [12.5, 0.0545]
  cd(0:1, 88 ) = [13.0, 0.0633]
  cd(0:1, 89 ) = [13.5, 0.0722]
  cd(0:1, 90 ) = [14.0, 0.0806]
  cd(0:1, 91 ) = [14.5, 0.0900]
  cd(0:1, 92 ) = [15.0, 0.0987]
  cd(0:1, 93 ) = [15.5, 0.1075]
  cd(0:1, 94 ) = [16.0, 0.1170]
  cd(0:1, 95 ) = [16.5, 0.1270]
  cd(0:1, 96 ) = [17.0, 0.1368]
  cd(0:1, 97 ) = [17.5, 0.1464]
  cd(0:1, 98 ) = [18.0, 0.1562]
  cd(0:1, 99 ) = [18.5, 0.1664]
  cd(0:1, 100) = [19.0, 0.1770]
  cd(0:1, 101) = [19.5, 0.1878]
  cd(0:1, 102) = [20.0, 0.1987]
  cd(0:1, 103) = [20.5, 0.2100]
  cd(0:1, 104) = [21.0, 0.2214]
  cd(0:1, 105) = [22.0, 0.2499]
  cd(0:1, 106) = [23.0, 0.2786]
  cd(0:1, 107) = [24.0, 0.3077]
  cd(0:1, 108) = [25.0, 0.3371]
  cd(0:1, 109) = [26.0, 0.3664]
  cd(0:1, 110) = [28.0, 0.4246]
  cd(0:1, 111) = [30.0, 0.4813]
  cd(0:1, 112) = [32.0, 0.5356]
  cd(0:1, 113) = [35.0, 0.6127]
  cd(0:1, 114) = [40.0, 0.7396]
  cd(0:1, 115) = [45.0, 0.8623]
  cd(0:1, 116) = [50.0, 0.9781]
  cd(0:1, 117) = [55.0, 1.0846]
  cd(0:1, 118) = [60.0, 1.1796]
  cd(0:1, 119) = [65.0, 1.2617]
  cd(0:1, 120) = [70.0, 1.3297]
  cd(0:1, 121) = [75.0, 1.3827]
  cd(0:1, 122) = [80.0, 1.4202]
  cd(0:1, 123) = [85.0, 1.4423]
  cd(0:1, 124) = [90.0, 1.4512]
  cd(0:1, 125) = [95.0, 1.4480]
  cd(0:1, 126) = [100.0, 1.4294]
  cd(0:1, 127) = [105.0, 1.3954]
  cd(0:1, 128) = [110.0, 1.3464]
  cd(0:1, 129) = [115.0, 1.2829]
  cd(0:1, 130) = [120.0, 1.2057]
  cd(0:1, 131) = [125.0, 1.1157]
  cd(0:1, 132) = [130.0, 1.0144]
  cd(0:1, 133) = [135.0, 0.9033]
  cd(0:1, 134) = [140.0, 0.7845]
  cd(0:1, 135) = [145.0, 0.6605]
  cd(0:1, 136) = [150.0, 0.5346]
  cd(0:1, 137) = [155.0, 0.4103]
  cd(0:1, 138) = [160.0, 0.2922]
  cd(0:1, 139) = [170.0, 0.0969]
  cd(0:1, 140) = [175.0, 0.0334]
  cd(0:1, 141) = [180.0, 0.0185]
  
  IF ( y .EQ. cd(0,0) ) cd7 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd7 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd7 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd7_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd7 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd7_nrel5.pro'
  ENDIF
  
!print,y,cd7
  
!_________________________________________________________________________
  
END FUNCTION calc_cd7_nrel5
!_________________________________________________________________________

!=========================================================================

Function calc_cd8_nrel5( y ) RESULT( cd8 )
!_________________________________________________________________________
!
! Data from NREL 5 MW: NACA64_A17
!_________________________________________________________________________
  
  REAL, INTENT( IN ) :: y

  INTEGER, PARAMETER :: na = 127

  REAL, DIMENSION( 0:1,0:na-1 ) :: cd

  REAL :: cd8
  INTEGER :: i

  cd8 = 0.0

!              alpha    CD  
  cd(0:1, 0  ) = [-180.0, 0.0198]
  cd(0:1, 1  ) = [-175.0, 0.0341]
  cd(0:1, 2  ) = [-170.0, 0.0955]
  cd(0:1, 3  ) = [-160.0, 0.2807]
  cd(0:1, 4  ) = [-155.0, 0.3919]
  cd(0:1, 5  ) = [-150.0, 0.5086]
  cd(0:1, 6  ) = [-145.0, 0.6267]
  cd(0:1, 7  ) = [-140.0, 0.7427]
  cd(0:1, 8  ) = [-135.0, 0.8537]
  cd(0:1, 9  ) = [-130.0, 0.9574]
  cd(0:1, 10 ) = [-125.0, 1.0519]
  cd(0:1, 11 ) = [-120.0, 1.1355]
  cd(0:1, 12 ) = [-115.0, 1.2070]
  cd(0:1, 13 ) = [-110.0, 1.2656]
  cd(0:1, 14 ) = [-105.0, 1.3104]
  cd(0:1, 15 ) = [-100.0, 1.3410]
  cd(0:1, 16 ) = [-95.0, 1.3572]
  cd(0:1, 17 ) = [-90.0, 1.3587]
  cd(0:1, 18 ) = [-85.0, 1.3456]
  cd(0:1, 19 ) = [-80.0, 1.3181]
  cd(0:1, 20 ) = [-75.0, 1.2765]
  cd(0:1, 21 ) = [-70.0, 1.2212]
  cd(0:1, 22 ) = [-65.0, 1.1532]
  cd(0:1, 23 ) = [-60.0, 1.0731]
  cd(0:1, 24 ) = [-55.0, 0.9822]
  cd(0:1, 25 ) = [-50.0, 0.8820]
  cd(0:1, 26 ) = [-45.0, 0.7742]
  cd(0:1, 27 ) = [-40.0, 0.6610]
  cd(0:1, 28 ) = [-35.0, 0.5451]
  cd(0:1, 29 ) = [-30.0, 0.4295]
  cd(0:1, 30 ) = [-25.0, 0.3071]
  cd(0:1, 31 ) = [-24.0, 0.2814]
  cd(0:1, 32 ) = [-23.0, 0.2556]
  cd(0:1, 33 ) = [-22.0, 0.2297]
  cd(0:1, 34 ) = [-21.0, 0.2040]
  cd(0:1, 35 ) = [-20.0, 0.1785]
  cd(0:1, 36 ) = [-19.0, 0.1534]
  cd(0:1, 37 ) = [-18.0, 0.1288]
  cd(0:1, 38 ) = [-17.0, 0.1037]
  cd(0:1, 39 ) = [-16.0, 0.0786]
  cd(0:1, 40 ) = [-15.0, 0.0535]
  cd(0:1, 41 ) = [-14.0, 0.0283]
  cd(0:1, 42 ) = [-13.5, 0.0158]
  cd(0:1, 43 ) = [-13.0, 0.0151]
  cd(0:1, 44 ) = [-12.0, 0.0134]
  cd(0:1, 45 ) = [-11.0, 0.0121]
  cd(0:1, 46 ) = [-10.0, 0.0111]
  cd(0:1, 47 ) = [-9.0, 0.0099]
  cd(0:1, 48 ) = [-8.0, 0.0091]
  cd(0:1, 49 ) = [-7.0, 0.0086]
  cd(0:1, 50 ) = [-6.0, 0.0082]
  cd(0:1, 51 ) = [-5.0, 0.0079]
  cd(0:1, 52 ) = [-4.0, 0.0072]
  cd(0:1, 53 ) = [-3.0, 0.0064]
  cd(0:1, 54 ) = [-2.0, 0.0054]
  cd(0:1, 55 ) = [-1.0, 0.0052]
  cd(0:1, 56 ) = [0.0, 0.0052]
  cd(0:1, 57 ) = [1.0, 0.0052]
  cd(0:1, 58 ) = [2.0, 0.0053]
  cd(0:1, 59 ) = [3.0, 0.0053]
  cd(0:1, 60 ) = [4.0, 0.0054]
  cd(0:1, 61 ) = [5.0, 0.0058]
  cd(0:1, 62 ) = [6.0, 0.0091]
  cd(0:1, 63 ) = [7.0, 0.0113]
  cd(0:1, 64 ) = [8.0, 0.0124]
  cd(0:1, 65 ) = [8.5, 0.0130]
  cd(0:1, 66 ) = [9.0, 0.0136]
  cd(0:1, 67 ) = [9.5, 0.0143]
  cd(0:1, 68 ) = [10.0, 0.0150]
  cd(0:1, 69 ) = [10.5, 0.0267]
  cd(0:1, 70 ) = [11.0, 0.0383]
  cd(0:1, 71 ) = [11.5, 0.0498]
  cd(0:1, 72 ) = [12.0, 0.0613]
  cd(0:1, 73 ) = [12.5, 0.0727]
  cd(0:1, 74 ) = [13.0, 0.0841]
  cd(0:1, 75 ) = [13.5, 0.0954]
  cd(0:1, 76 ) = [14.0, 0.1065]
  cd(0:1, 77 ) = [14.5, 0.1176]
  cd(0:1, 78 ) = [15.0, 0.1287]
  cd(0:1, 79 ) = [15.5, 0.1398]
  cd(0:1, 80 ) = [16.0, 0.1509]
  cd(0:1, 81 ) = [16.5, 0.1619]
  cd(0:1, 82 ) = [17.0, 0.1728]
  cd(0:1, 83 ) = [17.5, 0.1837]
  cd(0:1, 84 ) = [18.0, 0.1947]
  cd(0:1, 85 ) = [18.5, 0.2057]
  cd(0:1, 86 ) = [19.0, 0.2165]
  cd(0:1, 87 ) = [19.5, 0.2272]
  cd(0:1, 88 ) = [20.0, 0.2379]
  cd(0:1, 89 ) = [21.0, 0.2590]
  cd(0:1, 90 ) = [22.0, 0.2799]
  cd(0:1, 91 ) = [23.0, 0.3004]
  cd(0:1, 92 ) = [24.0, 0.3204]
  cd(0:1, 93 ) = [25.0, 0.3377]
  cd(0:1, 94 ) = [26.0, 0.3554]
  cd(0:1, 95 ) = [28.0, 0.3916]
  cd(0:1, 96 ) = [30.0, 0.4294]
  cd(0:1, 97 ) = [32.0, 0.4690]
  cd(0:1, 98 ) = [35.0, 0.5324]
  cd(0:1, 99 ) = [40.0, 0.6452]
  cd(0:1, 100) = [45.0, 0.7573]
  cd(0:1, 101) = [50.0, 0.8664]
  cd(0:1, 102) = [55.0, 0.9708]
  cd(0:1, 103) = [60.0, 1.0693]
  cd(0:1, 104) = [65.0, 1.1606]
  cd(0:1, 105) = [70.0, 1.2438]
  cd(0:1, 106) = [75.0, 1.3178]
  cd(0:1, 107) = [80.0, 1.3809]
  cd(0:1, 108) = [85.0, 1.4304]
  cd(0:1, 109) = [90.0, 1.4565]
  cd(0:1, 110) = [95.0, 1.4533]
  cd(0:1, 111) = [100.0, 1.4345]
  cd(0:1, 112) = [105.0, 1.4004]
  cd(0:1, 113) = [110.0, 1.3512]
  cd(0:1, 114) = [115.0, 1.2874]
  cd(0:1, 115) = [120.0, 1.2099]
  cd(0:1, 116) = [125.0, 1.1196]
  cd(0:1, 117) = [130.0, 1.0179]
  cd(0:1, 118) = [135.0, 0.9064]
  cd(0:1, 119) = [140.0, 0.7871]
  cd(0:1, 120) = [145.0, 0.6627]
  cd(0:1, 121) = [150.0, 0.5363]
  cd(0:1, 122) = [155.0, 0.4116]
  cd(0:1, 123) = [160.0, 0.2931]
  cd(0:1, 124) = [170.0, 0.0971]
  cd(0:1, 125) = [175.0, 0.0334]
  cd(0:1, 126) = [180.0, 0.0198]

  IF ( y .EQ. cd(0,0) ) cd8 = cd(1,0)
  DO i = 1,na-1
     IF ( ( y .GT. cd(0,i-1) ) .AND. (y .LE. cd(0,i) ) ) THEN
        cd8 = cd(1,i-1) + (y - cd(0,i-1))*( ( cd(1,i) - cd(1,i-1) )/( cd(0,i) - cd(0,i-1) ) )
        GOTO 1
     ENDIF
  ENDDO 
1 CONTINUE
  
  IF ( y .LT. cd(0,0) ) THEN 
     cd8 = cd(1,0)
!     PRINT*,'alpha < tabulated range in calc_cd8_nrel5.pro'
  ENDIF
  IF ( y .GT. cd(0,na-1) ) THEN
     cd8 = cd(1,na-1)
!     PRINT*,'alpha greater than tabulated range in calc_cd8_nrel5.pro'
  ENDIF
  
!print,y,cd8
  
!_________________________________________________________________________
  
END FUNCTION calc_cd8_nrel5
!_________________________________________________________________________

END MODULE module_gad_nrel5

