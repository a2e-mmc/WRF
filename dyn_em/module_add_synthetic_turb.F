!WRF prog
!module
!
MODULE module_add_synthetic_turb
#ifdef DM_PARALLEL
        USE MODULE_DM
#endif
!USE module_model_constants
CONTAINS

! Raj Rai wrote for MMC -- Apply synthetic turbulence along the inflow boundary
SUBROUTINE read_uvw_fluctuation(xtime,       & ! time-passed after simulation begins 
                                u_fluctu,    & ! u-velocity (synthetic) 
                                v_fluctu,    & ! v-velocity (synthetic)
                                w_fluctu,    & ! w-velocity (synthetic)
                                t_valu,      & ! time in minutes in 
                                nxy_Cnfg,    & ! nx or ny in inflow data
                                nz_Cnfg,     & ! nz in inflow  
                                nt_Cnfg,     & ! nt in inflow
                                stSynth,     & ! starttime, e.g. 120 after 2 hr simulation  
                                fname)         ! file name
implicit none
integer, intent(in)                       :: nxy_Cnfg, nz_Cnfg, nt_Cnfg
real, intent(in)                          :: xtime, stSynth
real, dimension(1:nxy_Cnfg,1:nz_Cnfg,1:nt_Cnfg),intent(inout) :: u_fluctu,v_fluctu, w_fluctu
real, dimension(1:nt_Cnfg), intent(inout)                     :: t_valu
character (28), intent(in) :: fname

! local variables
integer        :: nxyTemp, nzTemp, i, j, k, k3 
integer        :: k21, k22, k23, k11, k12, k13
real           :: t_min, add2xtime 
!character (28) :: fname

add2xtime = stSynth

open(21,file=trim(fname),form='formatted',status='old')

! reading file
read(21,*) nxyTemp, nzTemp
! check if the gridpoints in the file and namelists are different
if ((nxyTemp.ne.nxy_Cnfg).or.(nzTemp.ne.nz_Cnfg)) then
    print*,'Check the  inputfile for datasize'
end if

do k3=1,nt_Cnfg
   read(21,*) t_min
   !write(*,*) t_min
   t_valu(k3)=t_min+add2xtime
   do k21 = 1,nz_Cnfg
      read(21,*) (u_fluctu(k11,k21,k3),k11=1,nxy_Cnfg)
   end do
   do k22 =1,nz_Cnfg
      read(21,*) (v_fluctu(k12,k22,k3),k12=1,nxy_Cnfg)
   end do
   do k23 = 1,nz_Cnfg
      read(21,*) (w_fluctu(k13,k23,k3),k13=1,nxy_Cnfg)
   end do
end do
close(21)

END SUBROUTINE read_uvw_fluctuation 

!---------------------                ----------------
!-----------------------------------------------------
!---------------------                ----------------

SUBROUTINE add_synthetic_turb(muu, muv, mut,                & 
                              xtime,                        &
                              rdz,                          &
                              ru_tendf,rv_tendf, rw_tendf,  &
                              uFluctu,vFluctu, wFluctu,     & !syntheticturbulence
                              tvalu,                        & !time in minute
                              nxy_Cnfg, nz_Cnfg, nt_Cnfg,   & !dimofu/v/WFluctu
                              sideSynth,                    & !sidename, e.g.,'S','W',etc.
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte)

implicit none
integer, intent(in)                            :: nxy_Cnfg, nz_Cnfg, nt_Cnfg
integer, intent(in)                            :: ids, ide, jds, jde, kds, kde, &
                                                  ims, ime, jms, jme, kms, kme, &
                                                  its, ite, jts, jte, kts, kte
real, intent(in)                               :: xtime
character (50), intent(in)                     :: sideSynth
real, dimension(1:nt_Cnfg),intent(in)          :: tvalu
real, dimension(ims:ime,jms:jme), intent(in)               :: muu, muv, mut
real, dimension(ims:ime,kms:kme,jms:jme),intent(in)        :: rdz
real, dimension(ims:ime, kms:kme, jms:jme), intent(inout)  :: ru_tendf,rv_tendf, rw_tendf
real, dimension(1:nxy_Cnfg,1:nz_Cnfg,1:nt_Cnfg),intent(in) :: uFluctu,vFluctu, wFluctu

! local variables
integer                       :: i, ii, i1, k, k1, j, j1 
integer                       :: nt, tc, stpw, stpe, stps, stpn
integer                       :: k_synth, z_synth, z_synthVal
real                          :: fac_tc, fac_z,tms
real                          :: utend_mass, vtend_mass, wtend_mass
real, dimension(1:nz_Cnfg)                       :: zvertUVW
real, dimension(its-1:ite+1,kts:kte,jts-1:jte+1) :: zm

!write(*,*) 'size ... ', nxy_Cnfg, nz_Cnfg, nt_Cnfg

tms = 1.0
! finding current time
do nt = 1, nt_Cnfg-1
   if ((xtime.ge.tvalu(nt)).and.(xtime.lt.tvalu(nt+1))) exit
end do
tc = nt 
fac_tc = (xtime-tvalu(tc))/(tvalu(tc+1)-tvalu(tc))

!print*, 'xtime->', xtime, 'tvalu(tc)->', tvalu(tc+1), 'fac_tc->', fac_tc 
! vertical heights for sythetic data
do i =1,nz_Cnfg
   zvertUVW(i) = real(i)*10. !since vertical grid-spacing 10m
end do
  
! find z in each mass point
do j1 = jts-1, jte+1
   do i1 = its-1, ite+1
      zm(i1,kts,j1)=1.0/rdz(i1,kts,j1)
      do k1 = kts+1,min(kte,kde-1)
         zm(i1,k1,j1)=zm(i1,k1-1,j1)+1.0/rdz(i1,k1,j1)
      end do
   end do
end do
! adding perturbation along boundaries
if (index(trim(sideSynth),'W').ne.0) then !===============West boundary================
   if ((ims.le.1).and.(ime.ge.1)) then
      stpw = abs(ims)+1
      !print*,stpw
      do j = jts,jte
         ! find model heights within sythetic data
         do z_synth = kts, kte-1
            if ((zvertUVW(nz_Cnfg)-zm(stpw,z_synth,j)).le.0.) exit
         end do
         z_synthVal=z_synth-1
         do k = kts,z_synthVal
            if (zm(stpw,k,j).le.zvertUVW(1)) then !fluctuations below firstsynthetic data
               utend_mass=(1.-fac_tc)*uFluctu(j,k,tc)*tms+fac_tc*uFluctu(j,k,tc+1)*tms
               vtend_mass=(1.-fac_tc)*vFluctu(j,k,tc)*tms+fac_tc*vFluctu(j,k,tc+1)*tms
               wtend_mass=(1.-fac_tc)*wFluctu(j,k,tc)*tms+fac_tc*wFluctu(j,k,tc+1)*tms
            else if (zm(stpw,k,j).gt.zvertUVW(1)) then !interpolate databetween
               do k_synth=1,nz_Cnfg-1
                  if((zm(stpw,k,j).ge.zvertUVW(k_synth)).and.(zm(stpw,k,j).lt.zvertUVW(k_synth+1))) exit
               end do
               fac_z      = (zm(stpw,k,j)-zvertUVW(k_synth))/(zvertUVW(k_synth+1)-zvertUVW(k_synth))
               utend_mass = ((1.0-fac_z)*((1.-fac_tc)*uFluctu(j,k_synth,tc)*tms+fac_tc*uFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*uFluctu(j,k_synth+1,tc)*tms+fac_tc*uFluctu(j,k_synth+1,tc+1)*tms))*muu(stpw,j)
               vtend_mass = ((1.0-fac_z)*((1.-fac_tc)*vFluctu(j,k_synth,tc)*tms+fac_tc*vFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*vFluctu(j,k_synth+1,tc)*tms+fac_tc*vFluctu(j,k_synth+1,tc+1)*tms))*muv(stpw,j)
               wtend_mass = ((1.0-fac_z)*((1.-fac_tc)*wFluctu(j,k_synth,tc)*tms+fac_tc*wFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*wFluctu(j,k_synth+1,tc)*tms+fac_tc*wFluctu(j,k_synth+1,tc+1)*tms))*mut(stpw,j)
            else
               print*,'west- check vertical resolution'
            end if
            ru_tendf(stpw,k,j) = ru_tendf(stpw,k,j)+utend_mass
            rv_tendf(stpw,k,j) = rv_tendf(stpw,k,j)+vtend_mass
            rw_tendf(stpw,k,j) = rw_tendf(stpw,k,j)+wtend_mass
            !print*, 'West:ru_tendf->', ru_tendf(stpw,k,j),'West:utend_mass->',utend_mass 
         end do
      end do
   end if
end if


if (index(trim(sideSynth),'E').ne.0) then !==============East boundary=================
   if ((ims.le.ide).and.(ime.ge.ide)) then
      stpe = abs(ide)-4
      do j = jts,jte
         ! find model heights within sythetic data
         do z_synth = kts, kte-1
            if ((zvertUVW(nz_Cnfg)-zm(stpe,z_synth,j)).le.0.) exit
         end do
         z_synthVal=z_synth-1
         do k = kts,z_synthVal
            if (zm(stpe,k,j).le.zvertUVW(1)) then !fluctuations below firstsynthetic data
               utend_mass = (1.-fac_tc)*uFluctu(j,k,tc)*tms+fac_tc*uFluctu(j,k,tc+1)*tms
               vtend_mass = (1.-fac_tc)*vFluctu(j,k,tc)*tms+fac_tc*vFluctu(j,k,tc+1)*tms
               wtend_mass = (1.-fac_tc)*wFluctu(j,k,tc)*tms+fac_tc*wFluctu(j,k,tc+1)*tms
            else if (zm(stpe,k,j).gt.zvertUVW(1)) then !interpolate databetween
               do k_synth=1,nz_Cnfg-1
                  if((zm(stpe,k,j).ge.zvertUVW(k_synth)).and.(zm(stpe,k,j).lt.zvertUVW(k_synth+1))) exit
               end do
               fac_z      = (zm(stpe,k,j)-zvertUVW(k_synth))/(zvertUVW(k_synth+1)-zvertUVW(k_synth))
               utend_mass = ((1.0-fac_z)*((1.-fac_tc)*uFluctu(j,k_synth,tc)*tms+fac_tc*uFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*uFluctu(j,k_synth+1,tc)*tms+fac_tc*uFluctu(j,k_synth+1,tc+1)*tms))*muu(stpe,j)
               vtend_mass = ((1.0-fac_z)*((1.-fac_tc)*vFluctu(j,k_synth,tc)*tms+fac_tc*vFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*vFluctu(j,k_synth+1,tc)*tms+fac_tc*vFluctu(j,k_synth+1,tc+1)*tms))*muv(stpe,j)
               wtend_mass = ((1.0-fac_z)*((1.-fac_tc)*wFluctu(j,k_synth,tc)*tms+fac_tc*wFluctu(j,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*wFluctu(j,k_synth+1,tc)*tms+fac_tc*wFluctu(j,k_synth+1,tc+1)*tms))*mut(stpe,j)
            else
               print*,'east- check vertical resolution'
            end if
            ru_tendf(stpe,k,j) = ru_tendf(stpe,k,j)+utend_mass
            rv_tendf(stpe,k,j) = rv_tendf(stpe,k,j)+vtend_mass
            rw_tendf(stpe,k,j) = rw_tendf(stpe,k,j)+wtend_mass
         end do
      end do
   end if
end if

if (index(trim(sideSynth),'S').ne.0) then !===========South boundary==================
   if ((jms.le.1).and.(jme.ge.1)) then          
      stps = abs(jms)+1                         
      do i = its,ite
         ! find model heights within sythetic data
         do z_synth = kts, kte-1
           if ((zvertUVW(nz_Cnfg)-zm(i,z_synth,stps)).le.0.) exit
         end do
         z_synthVal=z_synth-1
         do k = kts,z_synthVal
            if (zm(i,k,stps).le.zvertUVW(1)) then !fluctuations belowfirst synthetic data
               utend_mass=(1.-fac_tc)*uFluctu(i,k,tc)*tms+fac_tc*uFluctu(i,k,tc+1)*tms
               vtend_mass=(1.-fac_tc)*vFluctu(i,k,tc)*tms+fac_tc*vFluctu(i,k,tc+1)*tms
               wtend_mass=(1.-fac_tc)*wFluctu(i,k,tc)*tms+fac_tc*wFluctu(i,k,tc+1)*tms
            else if (zm(i,k,stps).gt.zvertUVW(1)) then !interpolate databetween
               do k_synth=1,nz_Cnfg-1
                  if((zm(i,k,stps).ge.zvertUVW(k_synth)).and.(zm(i,k,stps).lt.zvertUVW(k_synth+1))) exit
               end do
               fac_z      = (zm(i,k,stps)-zvertUVW(k_synth))/(zvertUVW(k_synth+1)-zvertUVW(k_synth))
               utend_mass = ((1.0-fac_z)*((1.-fac_tc)*uFluctu(i,k_synth,tc)*tms+fac_tc*uFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*uFluctu(i,k_synth+1,tc)*tms+fac_tc*uFluctu(i,k_synth+1,tc+1)*tms))*muu(i,stps)
               vtend_mass = ((1.0-fac_z)*((1.-fac_tc)*vFluctu(i,k_synth,tc)*tms+fac_tc*vFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*vFluctu(i,k_synth+1,tc)*tms+fac_tc*vFluctu(i,k_synth+1,tc+1)*tms))*muv(i,stps)
               wtend_mass = ((1.0-fac_z)*((1.-fac_tc)*wFluctu(i,k_synth,tc)*tms+fac_tc*wFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*wFluctu(i,k_synth+1,tc)*tms+fac_tc*wFluctu(i,k_synth+1,tc+1)*tms))*mut(i,stps)
            else
               print*,'south- check vertical resolution'
            end if
            ru_tendf(i,k,stps) = ru_tendf(i,k,stps)+utend_mass
            rv_tendf(i,k,stps) = rv_tendf(i,k,stps)+vtend_mass
            rw_tendf(i,k,stps) = rw_tendf(i,k,stps)+wtend_mass
            !print*, 'South:rw_tendf,wtend_mass->',rw_tendf(i,k,stps),wtend_mass
         end do
      end do
   end if
end if

if (index(trim(sideSynth),'N').ne.0) then !==========North boundary====================
   if ((jms.le.jde).and.(jme.ge.jde)) then   
      stpn = abs(jde)-4                       
      do i = its,ite
         ! find model heights within sythetic data
         do z_synth = kts, kte-1
            if ((zvertUVW(nz_Cnfg)-zm(i,z_synth,stpn)).le.0.) exit
         end do
         z_synthVal=z_synth-1
         do k = kts,z_synthVal
            if (zm(i,k,stpn).le.zvertUVW(1)) then !fluctuations belowfirst synthetic data
               utend_mass=(1.-fac_tc)*uFluctu(i,k,tc)*tms+fac_tc*uFluctu(i,k,tc+1)*tms
               vtend_mass=(1.-fac_tc)*vFluctu(i,k,tc)*tms+fac_tc*vFluctu(i,k,tc+1)*tms
               wtend_mass=(1.-fac_tc)*wFluctu(i,k,tc)*tms+fac_tc*wFluctu(i,k,tc+1)*tms
            else if (zm(i,k,stpn).gt.zvertUVW(1)) then !interpolate databetween
               do k_synth=1,nz_Cnfg-1
                  if((zm(i,k,stpn).ge.zvertUVW(k_synth)).and.(zm(i,k,stpn).lt.zvertUVW(k_synth+1))) exit
               end do
               fac_z      = (zm(i,k,stpn)-zvertUVW(k_synth))/(zvertUVW(k_synth+1)-zvertUVW(k_synth))
               utend_mass = ((1.0-fac_z)*((1.-fac_tc)*uFluctu(i,k_synth,tc)*tms+fac_tc*uFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*uFluctu(i,k_synth+1,tc)*tms+fac_tc*uFluctu(i,k_synth+1,tc+1)*tms))*muu(i,stpn)
               vtend_mass = ((1.0-fac_z)*((1.-fac_tc)*vFluctu(i,k_synth,tc)*tms+fac_tc*vFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*vFluctu(i,k_synth+1,tc)*tms+fac_tc*vFluctu(i,k_synth+1,tc+1)*tms))*muv(i,stpn)
               wtend_mass = ((1.0-fac_z)*((1.-fac_tc)*wFluctu(i,k_synth,tc)*tms+fac_tc*wFluctu(i,k_synth,tc+1)*tms)+&
                            fac_z*((1.-fac_tc)*wFluctu(i,k_synth+1,tc)*tms+fac_tc*wFluctu(i,k_synth+1,tc+1)*tms))*mut(i,stpn)
            else
               print*,'north- check vertical resolution'
            end if
            ru_tendf(i,k,stpn) = ru_tendf(i,k,stpn)+utend_mass
            rv_tendf(i,k,stpn) = rv_tendf(i,k,stpn)+vtend_mass
            rw_tendf(i,k,stpn) = rw_tendf(i,k,stpn)+wtend_mass
         end do
      end do
   end if
endif

return
END SUBROUTINE add_synthetic_turb
END MODULE module_add_synthetic_turb


