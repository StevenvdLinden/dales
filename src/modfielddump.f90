!> \file modfielddump.f90
!!  Dumps 3D fields of several variables

!>
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myidx.myidy.expnr
!! If netcdf is true, this module leads the fielddump.myidx.myidy..expnr.nc output
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
module modfielddump

  use modglobal, only : longint, nsv

implicit none
private
PUBLIC :: initfielddump, fielddump,exitfielddump
save
!NetCDF variables
  integer :: nvar = 7
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.xxx.nc'
  character(80),dimension(:,:), allocatable :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, tmin, tmax
  integer(kind=longint) :: idtav,tnext,itmax,itmin
  integer :: klow,khigh,ncoarse=1
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !< switch for doing direct access writing (on/off)
  logical :: lu = .true.         !< switch for saving the u field
  logical :: lv = .true.         !< switch for saving the v field
  logical :: lw = .true.         !< switch for saving the w field
  logical :: lqt = .true.        !< switch for saving the qt field
  logical :: lql = .true.        !< switch for saving the ql field
  logical :: lthl = .true.       !< switch for saving the thl field
  logical :: lbuoy = .true.      !< switch for saving the buoy field
  logical :: lcli = .false.       !< switch for saving the cli field
  logical :: lclw = .false.       !< switch for saving the clw field
  logical :: lta = .false.        !< switch for saving the ta field
  logical :: lplw = .false.       !< switch for saving the plw field
  logical :: lpli = .false.       !< switch for saving the pli field
  logical :: lhus = .false.       !< switch for saving the hus field
  logical :: lhur = .false.       !< switch for saving the hur field
  logical :: ltntr = .false.      !< switch for saving the tntr field
  logical :: ltntrs = .false.     !< switch for saving the tntrs field
  logical :: ltntrl = .false.     !< switch for saving the tntrl field
  logical :: lsv(100) = .true.   !< switches for saving the sv fields

  ! indices for the variables in the netCDF vars array
  integer :: ind, ind_u=-1, ind_v=-1, ind_w=-1, ind_qt=-1, ind_ql=-1, ind_thl=-1, ind_buoy=-1, ind_sv(100)=-1
  integer :: ind_cli=-1, ind_clw=-1, ind_ta=-1, ind_plw=-1, ind_pli=-1, ind_hus=-1, ind_hur=-1, ind_tntr=-1, ind_tntrs=-1, ind_tntrl=-1

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,comm3d,mpi_logical,mpi_integer,myidx,myidy &
                       , D_MPI_BCAST
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime,tres,&
         checknamelisterror, output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modmicrodata, only : imicro, imicro_sice, imicro_sice2
    implicit none
    integer :: ierr, n
    character(3) :: csvname

    namelist/NAMFIELDDUMP/ &
         dtav,lfielddump,ldiracc,lbinary,klow,khigh,ncoarse, tmin, tmax,&
         lu, lv, lw, lqt, lql, lthl, lbuoy, lcli, lclw, lta, lplw, lpli, lhus, lhur, ltntr, ltntrs, ltntrl, lsv

    dtav=dtav_glob
    klow=1
    khigh=kmax
    tmin = 0.
    tmax = 1e8
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMFIELDDUMP')
      write(6 ,NAMFIELDDUMP)
      close(ifnamopt)

      if ((lcli .or. lclw) .and. .not. (imicro ==  imicro_sice .or. imicro == imicro_sice2)) then
         write (*,*) "FIELDDUMP: cli and clw output works only with simpleice microphysics. Turning off."
         lcli = .false.
         lclw = .false.
      end if
    end if
    call D_MPI_BCAST(ncoarse     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(klow        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(khigh       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(dtav        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmin        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(tmax        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lfielddump  ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ldiracc     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lbinary     ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lu          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lv          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lw          ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lqt         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lql         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lthl        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lbuoy       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lcli        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lclw        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lta         ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lplw        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lpli        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhus        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lhur        ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntr       ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrs      ,1,0,comm3d,ierr)
    call D_MPI_BCAST(ltntrl      ,1,0,comm3d,ierr)
    call D_MPI_BCAST(lsv       ,100,0,comm3d,ierr)

    idtav = dtav/tres
    itmin = tmin/tres
    itmax = tmax/tres

    tnext      = idtav   +btime
    if(.not.(lfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    if (lnetcdf) then
      write(fname,'(A,i3.3,A,i3.3,A)') 'fielddump.', myidx, '.', myidy, '.xxx.nc'
      fname(19:21) = cexpnr
      nvar = 17+nsv ! maximum number of variables
      allocate(ncname(nvar,4))
      call nctiminfo(tncname(1,:))
      ind = 1
      if (lu) then
         ind_u = ind
         ind = ind + 1
         call ncinfo(ncname(ind_u,:),'u','West-East velocity','m/s','mttt')
      end if
      if (lv) then
         ind_v = ind
         ind = ind + 1
         call ncinfo(ncname(ind_v,:),'v','South-North velocity','m/s','tmtt')
      end if
      if (lw) then
         ind_w = ind
         ind = ind + 1
         call ncinfo(ncname(ind_w,:),'w','Vertical velocity','m/s','ttmt')
      end if
      if (lqt) then
         ind_qt = ind
         ind = ind + 1
         call ncinfo(ncname(ind_qt,:),'qt','Total water specific humidity','kg/kg','tttt')
      end if
      if (lql) then
         ind_ql = ind
         ind = ind + 1
         call ncinfo(ncname(ind_ql,:),'ql','Liquid water specific humidity','kg/kg','tttt')
      end if
      if (lthl) then
         ind_thl = ind
         ind = ind + 1
         call ncinfo(ncname(ind_thl,:),'thl','Liquid water potential temperature','K','tttt')
      end if
      if (lbuoy) then
         ind_buoy = ind
         ind = ind + 1
         call ncinfo(ncname(ind_buoy,:),'buoy','Buoyancy','K','tttt')
      end if
      if (lcli) then
         ind_cli = ind
         ind = ind + 1
         call ncinfo(ncname(ind_cli,   :),'cli','mass fraction of cloud ice','kg/kg','tttt') ! new
      end if
      if (lclw) then
         ind_clw = ind
         ind = ind + 1
         call ncinfo(ncname(ind_clw,   :),'clw','mass fraction of cloud liquid water','kg/kg','tttt') ! new
      end if
      if (lta) then
         ind_ta = ind
         ind = ind + 1
         call ncinfo(ncname(ind_ta,    :),'ta','air temperature ','K','tttt') ! new
      end if
      if (lplw) then
         ind_plw = ind
         ind = ind + 1
         call ncinfo(ncname(ind_plw,   :),'plw','mass fraction of precipitating liquid water','kg/kg','tttt') !new
      end if
      if (lpli) then
         ind_pli = ind
         ind = ind + 1
         call ncinfo(ncname(ind_pli,   :),'pli','mass fraction of precipitating ice','kg/kg','tttt')   !new
      end if
      if (lhus) then
         ind_hus = ind
         ind = ind + 1
         call ncinfo(ncname(ind_hus,   :),'hus','specific humidity','kg/kg','tttt')   !new
      end if
      if (lhur) then
         ind_hur = ind
         ind = ind + 1
         call ncinfo(ncname(ind_hur,   :),'hur','relative humidity','%','tttt')      !new
      end if
      if (ltntr) then
         ind_tntr = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntr,  :),'tntr','tendency of air temperature due to radiative heating','K/s','tttt')   !new
      end if
      if (ltntrs) then
         ind_tntrs = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntrs, :),'tntrs','tendency of air temperature due to shortwave radiative heating','K/s','tttt')  !new
      end if
      if (ltntrl) then
         ind_tntrl = ind
         ind = ind + 1
         call ncinfo(ncname(ind_tntrl, :),'tntrl','tendency of air temperature due to longwave radiative heating','K/s','tttt')   !new
      end if

      do n=1,nsv
        if (lsv(n)) then
           ind_sv(n) = ind
           ind = ind + 1
           write (csvname(1:3),'(i3.3)') n
           call ncinfo(ncname(ind_sv(n),:),'sv'//csvname,'Scalar '//csvname//' specific concentration','(kg/kg)','tttt')
        end if
      end do
      nvar = ind - 1 ! total number of fields actually in use

      call open_nc(trim(output_prefix)//fname,  ncid,nrec,n1=ceiling(1.0*imax/ncoarse),n2=ceiling(1.0*jmax/ncoarse),n3=khigh-klow+1)

      if (nrec==0) then
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid, ncoarse, klow)
     end if
     call define_nc( ncid, NVar, ncname(1:NVar,:))
     ! must slice ncname here because define_nc expects first dimension to be NVar
     ! and NVar may have decreased if some fields are not saved
    end if

  end subroutine initfielddump

!> Do fielddump.
!> if lbinary, collect data to truncated (2 byte) integers, and write them to file
!> if lnetcdf, write to netCDF (as float32).
  subroutine fielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,thv0h,thvh,tmp0,rhof,exnf,presf
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,k1,rk3step,dzf, &
                          timee,dt_lim,cexpnr,ifoutput,rtimee,cp,tdn,tup
    use modmpi,    only : myid,cmyidx, cmyidy
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr, imicro, imicro_none, tuprsg, tdnrsg
    use modraddata, only   :lwu,lwd,swu,swd
    use modthermodynamics, only: qsat_tab
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:)
    real, allocatable :: vars(:,:,:,:)
    integer i,j,k,n,ii,jj
    integer :: writecounter = 1
    integer :: reclength


    if (.not. lfielddump) return
    if (rk3step/=3) return

    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if

    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    ! Only write fields if time is in the range (tmin, tmax)
    if (timee < itmin .or. timee > itmax) return

    if (lbinary) allocate(field(2-ih:i1+ih,2-jh:j1+jh,k1))
    if (lnetcdf) allocate(vars(ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1,nvar))

    reclength = ceiling(1.0*imax/ncoarse)*ceiling(1.0*jmax/ncoarse)*(khigh-klow+1)*2


    if (lnetcdf  .and. lu) vars(:,:,:,ind_u) = u0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E3*u0,2)
      if (ldiracc) then
        open (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbuu.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif


    if (lnetcdf .and. lv) vars(:,:,:,ind_v) = v0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E3*v0,2)
      if (ldiracc) then
        open (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbvv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif


    if (lnetcdf .and. lw) vars(:,:,:,ind_w) = w0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E3*w0,2)
      if (ldiracc) then
        open (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbww.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    if (lnetcdf .and. lqt) vars(:,:,:,ind_qt) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E5*qt0,2)
      if (ldiracc) then
        open (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbqt.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    if (lnetcdf .and. lql) vars(:,:,:,ind_ql) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E5*ql0,2)
      if (ldiracc) then
        open (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbql.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    if (lnetcdf .and. lthl) vars(:,:,:,ind_thl) = thl0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
    if (lbinary) then
      field = NINT(1.0E2*(thl0-300),2)
      if (ldiracc) then
        open (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbthl.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    end if

    if (lbinary) then
      if(imicro/=imicro_none) then
         do i=2-ih,i1+ih
            do j=2-jh,j1+jh
               do k=1,k1
                  field(i,j,k) = NINT(1.0E5*sv0(i,j,k,iqr),2)
               enddo
            enddo
         enddo
      else
         field = 0.
      endif
      if (ldiracc) then
        open (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbqr.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    ! buoyancy
    if (lnetcdf .and. lbuoy) then
      vars(:,:,:,ind_buoy) = thv0h(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      do k=klow,khigh
        vars(:,:,k-klow+1,ind_buoy) = vars(:,:,k-klow+1,ind_buoy) - thvh(k)
      end do
    end if

    if (lbinary) then
      do i=2-ih,i1+ih, ncoarse
         do j=2-jh,j1+jh, ncoarse
            do k=2,k1
               field(i,j,k) = NINT(1.0E2*(thv0h(i,j,k)-thvh(k)),2)
            enddo
         enddo
      enddo

      if (ldiracc) then
        open (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
        write (ifoutput, rec=writecounter) field(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)
      else
        open  (ifoutput,file='wbthv.'//cmyidx//'.'//cmyidy//'.'//cexpnr,form='unformatted',position='append')
        write (ifoutput) (((field(i,j,k),i=2,i1, ncoarse),j=2,j1, ncoarse),k=klow,khigh)
      end if
      close (ifoutput)
    endif

    if (lnetcdf .and. lcli) vars(:,:,:,ind_cli) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) * &
         (1 - max(0.,min(1.,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdn)/(tup-tdn))))

    if (lnetcdf .and. lclw) vars(:,:,:,ind_clw) = ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) * &
         max(0.,min(1.,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdn)/(tup-tdn)))

    if (lnetcdf .and. lta) vars(:,:,:,ind_ta) = tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    ! liquid and ice precip
    ! assuming simpleice is used
    !rsgratio_= max(0.,min(1.,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg))) ! rain vs snow/graupel partitioning   rsg = 1 if t > tuprsg
    if (lnetcdf .and. lplw) vars(:,:,:,ind_plw) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,iqr) * &
         max(0.,min(1.,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdnrsg)/(tuprsg-tdnrsg)))
    if (lnetcdf .and. lpli) vars(:,:,:,ind_pli) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,iqr)  * &
         (1 - max(0.,min(1.,(tmp0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)-tdnrsg)/(tuprsg-tdnrsg))))

    if (lnetcdf .and. lhus) vars(:,:,:,ind_hus) = qt0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh) - ql0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh)

    if (lnetcdf .and. lhur) then
       do k = klow, khigh
          do j = 1, ceiling(1.0*jmax/ncoarse)
             jj = 2 + (j-1) * ncoarse
             do i = 1, ceiling(1.0*imax/ncoarse)
                ii = 2 + (i-1) * ncoarse
                vars(i,j,k-klow+1,ind_hur) = 100 * (qt0(ii,jj,k) - ql0(ii,jj,k)) / &
                     qsat_tab(tmp0(ii,jj,k), presf(k))
             end do
          end do
       end do
    end if

    if (lnetcdf .and. ltntr) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntr) = ( &    !! need to loop over k here
               -swd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -swu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +swd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +swu(2:i1:ncoarse,2:j1:ncoarse,k) &
               -lwd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -lwu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +lwd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +lwu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    if (lnetcdf .and. ltntrs) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntrs) = ( & !! need to loop over k here
               -swd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -swu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +swd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +swu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    if (lnetcdf .and. ltntrl) then
       do k = klow,khigh
          vars(:,:,k-klow+1,ind_tntrl) = ( & !! need to loop over k here
               -lwd(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               -lwu(2:i1:ncoarse,2:j1:ncoarse,k+1) &
               +lwd(2:i1:ncoarse,2:j1:ncoarse,k) &
               +lwu(2:i1:ncoarse,2:j1:ncoarse,k) &
               )  / (rhof(k)*exnf(k)*cp*dzf(k))
       end do
    end if

    ! scalar variables
    if (lnetcdf) then
       do n=1,nsv
          if (lsv(n)) then
             vars(:,:,:,ind_sv(n)) = sv0(2:i1:ncoarse,2:j1:ncoarse,klow:khigh,n)
          end if
       end do
    end if


    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,ceiling(1.0*imax/ncoarse),ceiling(1.0*jmax/ncoarse),khigh-klow+1)
    end if

    if(lbinary) then
      if (myid==0) then
        open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
        write(ifoutput,'(F12.1 3F12.5)') timee,thls, qts,thvs
        close(ifoutput)
      end if
    endif

    writecounter=writecounter+1

    if (lbinary) deallocate(field)
    if (lnetcdf) deallocate(vars)

  end subroutine fielddump
!> Clean up when leaving the run
  subroutine exitfielddump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
