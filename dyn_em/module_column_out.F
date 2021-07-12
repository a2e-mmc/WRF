!WRF
!
MODULE module_column_out


CONTAINS

  SUBROUTINE compute_column_output( ucol, vcol, wcol, tcol, pcol,     &
                                    tkecol, t11col, t22col, t33col,   &
                                    t12col, t13col, t23col, h3col,    &
                                    zcol, uscol, z0col, col_i, col_j, &
                                    u, v, w, t, p, pb, ph, phb, tke,  &
                                    t11, t22, t33, t12, t13, t23, h3, &
                                    ust, z0, ht,                      &
                                    ids, ide, jds, jde, kds, kde,     &
                                    ims, ime, jms, jme, kms, kme,     &
                                    its, ite, jts, jte, kts, kte      )

!=======================================================================
!
!  This subroutine prepares a 1D column of parameters for an auxiliary 
!  output stream. All variables are at cell center horizontal locations.
!
!  The approach is to find the patch upon which the desired column
!  of data reside, and compute the variables (cell center, z, etc.).
!  Next, transfer the column information from the patch upon which it
!  resides onto the master processor where it will be written to the
!  output file. This is accomplished by setting each local patch
!  column value to a large negative number. Then, for the patch
!  upon which the column exists, we compute the desired quantities.
!  These values will be greater than the large negative values
!  originally prescribed. Then we use the built-in function
!  wrf_dm_max_real to select the maximum value from all the patches,
!  which recovers the desired column.
!    
!  Perhaps not the most elegant approach, but it seems to work.
!    
!  Wrtitten by Jeff Mirocha (LLNL) May, 2021
!
!======================================================================= 

#ifdef DM_PARALLEL
    USE module_dm
#endif
    
   IMPLICIT NONE

   REAL, DIMENSION(kms:kme), INTENT( OUT ) :: ucol, vcol, wcol, tcol, pcol,      &
                                              tkecol, t11col, t22col, t33col,    &
                                              t12col, t13col, t23col, h3col, zcol 

   REAL, INTENT( OUT ) :: uscol, z0col

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT( IN ) :: u, v, w,           &
                                                             t, p, pb,          &
                                                             ph, phb, tke,      &
                                                             t11, t22, t33,     & 
                                                             t12, t13, t23, h3 
                                                         
   REAL, DIMENSION(ims:ime,jms:jme), INTENT( IN ) :: ust, z0, ht

   INTEGER, INTENT( IN ) :: col_i, col_j

   INTEGER, INTENT( IN  ) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

! local variable
   INTEGER :: k

!-----------------------------------------------------------------------
! Executable starts here
!
   
   DO k=kts,kte - 1
      ucol(k) = -999.0
      vcol(k) = -999.0
      wcol(k) = -999.0
      tcol(k) = -999.0
      pcol(k) = -999.0
      zcol(k) = -999.0
      tkecol(k) = -999.0
      t11col(k) = -999.0
      t22col(k) = -999.0
      t33col(k) = -999.0
      t12col(k) = -999.0
      t13col(k) = -999.0
      t23col(k) = -999.0
      h3col(k)  = -999.0
   ENDDO

   wcol(kte) = -999.0
   zcol(kte) = -999.0
   t13col(kte) = -999.0
   t23col(kte) = -999.0
   h3col(kte)  = -999.0

   uscol = -999.0
   z0col = -999.0

! If this is the processor containing the desired column, compute/assign the values.


   IF ( ( col_i .GE. its ) .AND. ( col_i .LE. ite ) .AND. ( col_j .GE. jts ) .AND. ( col_j .LE. jte ) ) THEN

      DO k = kts,kte - 1

         ucol(k) = 0.5*( u(col_i,k,col_j) + u(col_i+1,k,col_j) )
         vcol(k) = 0.5*( v(col_i,k,col_j) + v(col_i,k,col_j+1) )
         wcol(k) = w(col_i,k,col_j)
         tcol(k) = t(col_i,k,col_j)
         pcol(k) = pb(col_i,k,col_j) + p(col_i,k,col_j)
         zcol(k) = ( phb(col_i,k,col_j) + ph(col_i,k,col_j) )/9.81 - ht(col_i,col_j)
         tkecol(k) = tke(col_i,k,col_j)
         t11col(k) = t11(col_i,k,col_j)
         t22col(k) = t22(col_i,k,col_j)      
         t33col(k) = t33(col_i,k,col_j)   
         t12col(k) = 0.25*( t12(col_i,k,col_j) + t12(col_i+1,k,col_j) + t12(col_i+1,k,col_j+1) + t12(col_i,k,col_j+1) )
         t13col(k) = 0.5*( t13(col_i,k,col_j) + t13(col_i+1,k,col_j) )
         t23col(k) = 0.5*( t23(col_i,k,col_j) + t23(col_i,k,col_j+1) )
         h3col(k)  = h3(col_i,k,col_j)
      ENDDO
   
      wcol(kte) = w(col_i,kte,col_j)
      zcol(kte) = ( phb(col_i,kte,col_j) + ph(col_i,kte,col_j) )/9.81 - ht(col_i,col_j)
      t13col(kte) = 0.5*( t13(col_i,kte,col_j) + t13(col_i+1,kte,col_j) )
      t23col(kte) = 0.5*( t23(col_i,kte,col_j) + t23(col_i,kte,col_j+1) )
      h3col(kte)  = h3(col_i,kte,col_j)
    
      uscol = ust(col_i,col_j)  
      z0col = z0(col_i,col_j)  
   
   ENDIF

! Select maximum value from all processors
   
   DO k=kts,kte - 1
      ucol(k) = wrf_dm_max_real(ucol(k)) 
      vcol(k) = wrf_dm_max_real(vcol(k))
      wcol(k) = wrf_dm_max_real(wcol(k))
      tcol(k) = wrf_dm_max_real(tcol(k))
      pcol(k) = wrf_dm_max_real(pcol(k))
      zcol(k) = wrf_dm_max_real(zcol(k))
      tkecol(k) = wrf_dm_max_real(tkecol(k))
      t11col(k) = wrf_dm_max_real(t11col(k))
      t22col(k) = wrf_dm_max_real(t22col(k))
      t33col(k) = wrf_dm_max_real(t33col(k))
      t12col(k) = wrf_dm_max_real(t12col(k))
      t13col(k) = wrf_dm_max_real(t13col(k))
      t23col(k) = wrf_dm_max_real(t23col(k))
      h3col(k)  = wrf_dm_max_real(h3col(k))
   ENDDO

   wcol(kte) =  wrf_dm_max_real(wcol(kte))
   zcol(kte) =  wrf_dm_max_real(zcol(kte))
   t13col(kte) =  wrf_dm_max_real(t13col(kte))
   t23col(kte) =  wrf_dm_max_real(t23col(kte))
   h3col(kte)  =  wrf_dm_max_real(h3col(kte))

   uscol =  wrf_dm_max_real(uscol)
   z0col =  wrf_dm_max_real(z0col)
   

   RETURN

END SUBROUTINE compute_column_output

END MODULE module_column_out
