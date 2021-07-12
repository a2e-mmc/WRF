!WRF
!
MODULE module_horizontal_slices

USE module_model_constants

CONTAINS

SUBROUTINE compute_horizontal_slices( slices_u, slices_v,                 &
                                      slices_w, slices_t,                 &
                                      slices_m13, slices_m23,             &
                                      slices_mh3, slices_tke,             & 
                                      slices_z,                           &
                                      u, v, w, t, m13, m23, mh3, tke,     &
                                      z, z_at_w, ht,                      &
                                      num_slices,                         &
                                      slice_heights,                      &
                                      ids, ide, jds, jde, kds, kde,       &
                                      ims, ime, jms, jme, kms, kme,       &
                                      its, ite, jts, jte, kts, kte         )

!=======================================================================
!
!  This subroutine calculates num_slices of horizontal planes of 
!  u, v, w, and t, interpolated to  slice_heights, for output.
!
!  J. Mirocha (LLNL) and P. Hawbecker (NCAR), January, 2021
!
!======================================================================= 

   IMPLICIT NONE

   REAL, DIMENSION(ims:ime,1:num_slices,jms:jme), INTENT( OUT ) :: slices_u, slices_v,     &
                                                                   slices_w, slices_t,     &
                                                                   slices_m13, slices_m23, &
                                                                   slices_mh3, slices_tke

   
   REAL, DIMENSION(1:num_slices), INTENT( OUT ) :: slices_z
   
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(  IN ) :: u, v, w, t, m13, m23, mh3, tke, z, z_at_w

   REAL, DIMENSION(ims:ime,jms:jme), INTENT(  IN ) :: ht
   
   INTEGER, INTENT(  IN ) :: num_slices
   
   REAL , DIMENSION (1:num_slices), INTENT( IN )  :: slice_heights
   
   INTEGER , INTENT( IN  ) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte

   ! local variables
   
   INTEGER :: i, j, k, k_slice

   REAL :: d_z, d_1, d_2, fac_1, fac_2
   REAL :: slice_height_abs
   REAL :: u_m_1, u_m_2, v_m_1, v_m_2
   
!-----------------------------------------------------------------------
! Executable code starts here
!

   DO k = 1, num_slices
 
     slices_z(k) = slice_heights(k)
     
   END DO

   
   DO k_slice = 1, num_slices

!      print*,k_slice,slice_heights(k_slice)

      DO i = its, MIN(ite,ide-1)
         DO j = jts, MIN(jte,jde-1)

            slice_height_abs = slice_heights(k_slice) + ht(i,j)
            
            DO k = kts+1, MIN(kte,kde-1)
               
               IF ( z_at_w(i,k,j) .GE. slice_height_abs ) THEN

                   d_z = z_at_w(i,k,j) - z_at_w(i,k-1,j)
                   d_2 = z_at_w(i,k,j) - slice_height_abs
                   d_1 = slice_height_abs - z_at_w(i,k-1,j)
                   fac_2 = d_2/d_z
                   fac_1 = d_1/d_z
             
                   slices_w(i,k_slice,j) = fac_2*w(i,k-1,j) + fac_1*w(i,k,j)
                   slices_m13(i,k_slice,j) = fac_2*( 0.5*( m13(i,k-1,j) + m13(i+1,k-1,j) ) ) + fac_1*( 0.5*( m13(i,k,j) + m13(i+1,k,j) ) )
                   slices_m23(i,k_slice,j) = fac_2*( 0.5*( m23(i,k-1,j) + m23(i,k-1,j+1) ) ) + fac_1*( 0.5*( m23(i,k,j) + m23(i,k,j+1) ) )
                   slices_mh3(i,k_slice,j) = fac_2*mh3(i,k-1,j) + fac_1*mh3(i,k,j)
                   
!                   IF ( (i .EQ. its+2) .AND. (j .EQ. jts+2) )print*,k,fac_1*z_at_w(i,k,j)+fac_2*z_at_w(i,k-1,j)
                   
                   GOTO 1
               
               ENDIF

            END DO !k

1 CONTINUE

            DO k = kts+1, MIN(kte,kde-1)
            
               IF ( z(i,k,j) .GE. slice_height_abs ) THEN

                  d_z = z(i,k,j) - z(i,k-1,j)
                  d_2 = z(i,k,j) - slice_height_abs
                  d_1 = slice_height_abs - z(i,k-1,j)

                  fac_2 = d_2/d_z
                  fac_1 = d_1/d_z
               
                  u_m_2 = 0.5*(u(i+1,k,j) + u(i,k,j) ) ! u at cell center level k
                  u_m_1 = 0.5*(u(i+1,k-1,j) + u(i,k-1,j) ) ! u at cell center level k-1
                  v_m_2 = 0.5*(v(i,k,j+1) + v(i,k,j) ) ! v at cell center level k
                  v_m_1 = 0.5*(v(i,k-1,j+1) + v(i,k-1,j) ) ! v at cell center level k-1

                  slices_u(i,k_slice,j) = fac_2*u_m_1 + fac_1*u_m_2
                  slices_v(i,k_slice,j) = fac_2*v_m_1 + fac_1*v_m_2
                  slices_t(i,k_slice,j) = fac_2*t(i,k-1,j) + fac_1*t(i,k,j)
                  slices_tke(i,k_slice,j) = fac_2*tke(i,k-1,j) + fac_1*tke(i,k,j)
                  
!                   IF ( (i .EQ. its+2) .AND. (j .EQ. jts+2) )print*,k,fac_1*z(i,k,j)+fac_2*z(i,k-1,j)
                  
                  GOTO 2
               
               ENDIF

            END DO !k
         
2       CONTINUE
            
         END DO
      END DO
   
   END DO !k_slice
   
   RETURN

END SUBROUTINE compute_horizontal_slices

  
END MODULE module_horizontal_slices
