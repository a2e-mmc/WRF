
! ======================================================================================
! This file was generated by the version 4.4.0 of ADG on 12/23/2010. The Adjoint Code
! Generator (ADG) was developed and sponsored by LASG of IAP (1999-2010)
! The Copyright of the ADG system was declared by Walls at LASG, 1999-2010
! ======================================================================================
! corrrected by zzma on 01/10/2011

MODULE a_module_microphysics_zero_out

CONTAINS

   SUBROUTINE a_microphysics_zero_outa(moist_new,a_moist_new,n_moist,config_flags, &
   ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

!PART I: DECLARATION OF VARIABLES

   USE module_state_description
   USE module_configure
   USE module_wrf_error

   IMPLICIT NONE

   INTEGER :: K0_ADJ,K1_ADJ,K2_ADJ,K3_ADJ
   TYPE(grid_config_rec_type) :: config_flags
   INTEGER :: ids,ide,jds,jde,kds,kde
   INTEGER :: ims,ime,jms,jme,kms,kme
   INTEGER :: its,ite,jts,jte,kts,kte
   INTEGER :: n_moist
   REAL,DIMENSION(ims:ime,kms:kme,jms:jme,n_moist) :: moist_new,a_moist_new
   INTEGER i,j,k,n

   IF( config_flags%mp_zero_out .EQ. 0 ) THEN
   ELSE IF( config_flags%mp_zero_out .EQ. 1 ) THEN
         CALL wrf_debug ( 100 , 'zero out small condensates, vapor not included')
     DO n =PARAM_FIRST_SCALAR, n_moist
       IF( n .NE. P_QV ) THEN
         DO j =jts, jte
         DO k =kts, kte
         DO i =its, ite
           IF( moist_new(i,k,j,n) .LT. config_flags%mp_zero_out_thresh ) THEN
             a_moist_new(i,k,j,n) =0.
           END IF
         ENDDO
         ENDDO
         ENDDO
       END IF
     ENDDO

   ELSE IF( config_flags%mp_zero_out .EQ. 2 ) then
         CALL wrf_debug ( 100 , 'zero out small condensates, zero out negative vapor')
   DO n =PARAM_FIRST_SCALAR, n_moist
     IF( n .NE. P_QV ) THEN

       DO j =jte, jts, -1
       DO k =kte, kts, -1
       DO i =ite, its, -1
       IF( moist_new(i,k,j,n) .LT. config_flags%mp_zero_out_thresh ) THEN
       a_moist_new(i,k,j,n) =0.0
       END IF
       ENDDO
       ENDDO
       ENDDO

     ELSE IF( n .EQ. P_QV ) THEN

       DO j =jte, jts, -1
       DO k =kte, kts, -1
       DO i =ite, its, -1
       a_moist_new(i,k,j,n) =(1.0 +(1.0)*sign(1.0,moist_new(i,k,j,n)-0.))*0.5*a_moist_new(i,k,j,n)
       ENDDO
       ENDDO
       ENDDO

     END IF
   ENDDO

   END IF

   END SUBROUTINE a_microphysics_zero_outa

   SUBROUTINE a_microphysics_zero_outb(moist_new,a_moist_new,n_moist,config_flags, &
   ids,ide,jds,jde,kds,kde,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

!PART I: DECLARATION OF VARIABLES

   USE module_state_description
   USE module_configure
   USE module_wrf_error

   IMPLICIT NONE

   INTEGER :: K0_ADJ,K1_ADJ,K2_ADJ,K3_ADJ
   TYPE(grid_config_rec_type) :: config_flags
   INTEGER :: ids,ide,jds,jde,kds,kde
   INTEGER :: ims,ime,jms,jme,kms,kme
   INTEGER :: its,ite,jts,jte,kts,kte
   INTEGER :: n_moist
   REAL,DIMENSION(ims:ime,kms:kme,jms:jme,n_moist) :: moist_new,a_moist_new
   INTEGER i,j,k,n

 IF( config_flags%mp_zero_out .NE. 0 ) THEN
   DO n =n_moist, PARAM_FIRST_SCALAR, -1

   i = ide-1
   IF( ( i .GE. its ) .AND. ( i .LE. MIN( ite , ide-1) ) ) THEN
   DO j =min(jte, jde-1), jts, -1
   DO k =kte, kts, -1
   a_moist_new(i,k,j,n) =(1.0 +(1.0)*sign(1.0,moist_new(i,k,j,n)-0.))*0.5*a_moist_new(i,k,j,n)
   ENDDO
   ENDDO
   END IF

   i = ids
   IF( ( i .GE. its ) .AND. ( i .LE. MIN( ite , ide-1) ) ) THEN
   DO j =min(jte, jde-1), jts, -1
   DO k =kte, kts, -1
   a_moist_new(i,k,j,n) =(1.0 +(1.0)*sign(1.0,moist_new(i,k,j,n)-0.))*0.5*a_moist_new(i,k,j,n)
   ENDDO
   ENDDO
   END IF

   j = jde-1
   IF( ( j .GE. jts ) .AND. ( j .LE. MIN( jte , jde-1) ) ) THEN
   DO k =kte, kts, -1
   DO i =min(ite, ide-1), its, -1
   a_moist_new(i,k,j,n) =(1.0 +(1.0)*sign(1.0,moist_new(i,k,j,n)-0.))*0.5*a_moist_new(i,k,j,n)
   ENDDO
   ENDDO
   END IF

   j = jds
   IF( ( j .GE. jts ) .AND. ( j .LE. MIN( jte , jde-1) ) ) THEN
   DO k =kte, kts, -1
   DO i =min(ite, ide-1), its, -1
   a_moist_new(i,k,j,n) =(1.0 +(1.0)*sign(1.0,moist_new(i,k,j,n)-0.))*0.5*a_moist_new(i,k,j,n)
   ENDDO
   ENDDO
   END IF

  ENDDO
 END IF


   RETURN
   END SUBROUTINE a_microphysics_zero_outb

   END MODULE a_module_microphysics_zero_out

