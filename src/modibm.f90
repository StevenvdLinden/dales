!> \file modibm.f90
!! SvdL, started 25 september 2024
!! By Michael Koene, email: mich4all@live.nl, TU Delft, section Atmospheric Physics, October 8, 2019
!! TODO: Write comments to output file
!!       clean up code (write statements)
!!       Test performance

module modibm
  use modprecision, only : field_r
  use modibmdata,   only : lapply_ibm, z0_wall, thlinit, thlwall, svwall, qtwall, e12wall, &
                           ibc_dirichlet, ibc_neumann, ibc_flux, ibc_thl, ibc_tke, ibc_qt, ibc_sv
                         !lpoislast, lreadfile_obstacles, lwallheat, &
                         !thlwall, thlroof,qtroof, thlibm,qtibm, &
                         !libm,Nair, ibm_adv_mask,kibm_max,&
  implicit none
  save
  public :: initibm, exitibm, applyibm

  integer :: nib_max, nfp_max
  integer :: n_idw_min = 4
  integer :: n_idw_pts = 8

  real(field_r) :: epsdis = 10**(-7.)               !< Small distance for which points are taken to "align"

  integer :: nfp_u                                  !< Number of u-points to be forced
  integer :: nfp_v                                  !< Number of v-points to be forced
  integer :: nfp_w                                  !< Number of w-points to be forced
  integer :: nfp_s                                  !< Number of scalar points to be forced
  integer :: nib_u                                  !< Number of internal u-points to be forced
  integer :: nib_v                                  !< Number of internal v-points to be forced
  integer :: nib_w                                  !< Number of internal w-points to be forced
  integer :: nib_s                                  !< Number of internal scalar points to be forced

  !< Fields, arrays and vectors for ibm (these "live" on each separate process)
  ! real(field_r), allocatable :: sdf_tmp(:,:,:)      !< Temporary field for full SDF     
  real(field_r), allocatable :: sdf_u  (:,:,:)      !< Signed Distance Function for the u-position      
  real(field_r), allocatable :: sdf_v  (:,:,:)      !< Signed Distance Function for the v-position
  real(field_r), allocatable :: sdf_w  (:,:,:)      !< Signed Distance Function for the w-position
  real(field_r), allocatable :: sdf_s  (:,:,:)      !< Signed Distance Function for the scalar position

  !< SvdL, 20240925: added allocatable arrays for x and y coordinates
  real(field_r), allocatable :: xh(:)               !<  x-position at center of grid cell [m]
  real(field_r), allocatable :: xf(:)               !<  x-position at face of grid cell [m]
  real(field_r), allocatable :: yh(:)               !<  y-position at center of grid cell [m]
  real(field_r), allocatable :: yf(:)               !<  y-position at face of grid cell [m]

  integer, allocatable :: fp_u(:,:)                 !< Indices of u-points to be forced
  integer, allocatable :: fp_v(:,:)                 !< Indices of v-points to be forced
  integer, allocatable :: fp_w(:,:)                 !< Indices of w-points to be forced
  integer, allocatable :: fp_s(:,:)                 !< Indices of scalar points to be forced

  integer, allocatable :: ib_u(:,:)                 !< Indices of internal u-points (= within boundary)
  integer, allocatable :: ib_v(:,:)                 !< Indices of internal v-points
  integer, allocatable :: ib_w(:,:)                 !< Indices of internal w-points
  integer, allocatable :: ib_s(:,:)                 !< Indices of internal scalar points

  integer, allocatable :: ipu_u(:,:,:)              !< Indices of u-interpolation points for forcing points of u >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
  integer, allocatable :: ipu_v(:,:,:)              !< v-interpolation points for forcing points of u
  integer, allocatable :: ipu_w(:,:,:)              !< w-interpolation points for forcing points of u
  integer, allocatable :: ipu_s(:,:,:)              !< scalar interpolation points for forcing points of u

  integer, allocatable :: ipv_u(:,:,:)              !< Indices of interpolation points for forcing points of v >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
  integer, allocatable :: ipv_v(:,:,:)              !< v-interpolation points for forcing points of v
  integer, allocatable :: ipv_w(:,:,:)              !< w-interpolation points for forcing points of v
  integer, allocatable :: ipv_s(:,:,:)              !< scalar interpolation points for forcing points of v

  integer, allocatable :: ipw_u(:,:,:)              !< Indices of interpolation points for forcing points of w >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
  integer, allocatable :: ipw_v(:,:,:)              !< v-interpolation points for forcing points of w
  integer, allocatable :: ipw_w(:,:,:)              !< w-interpolation points for forcing points of w
  integer, allocatable :: ipw_s(:,:,:)              !< scalar interpolation points for forcing points of w

  integer, allocatable :: ips_u(:,:,:)              !< Indices of interpolation points for forcing points of scalars >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
  integer, allocatable :: ips_v(:,:,:)              !< v-interpolation points for forcing points of scalars
  integer, allocatable :: ips_w(:,:,:)              !< w-interpolation points for forcing points of scalars
  integer, allocatable :: ips_s(:,:,:)              !< scalar interpolation points for forcing points of scalars

  real(field_r), allocatable :: loc_dist_u(:,:)     !< Nearest locations on boundary of forcing point + their distance, locations of interpolation points + their distance (xb,yb,zb,db,xi,yi,zi,di)
  real(field_r), allocatable :: loc_dist_v(:,:)     !< for v
  real(field_r), allocatable :: loc_dist_w(:,:)     !< for w
  real(field_r), allocatable :: loc_dist_s(:,:)     !< for scalars
  
  real(field_r), allocatable :: c_weightu_u(:,:)    !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of u. Dimensions: number of fp x 8 interpolation points
  real(field_r), allocatable :: c_weightu_v(:,:)    !< v-interpolation points for forcing points of u
  real(field_r), allocatable :: c_weightu_w(:,:)    !< w-interpolation points for forcing points of u
  real(field_r), allocatable :: c_weightu_s(:,:)    !< scalar interpolation points for forcing points of u

  real(field_r), allocatable :: c_weightv_u(:,:)    !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of v. Dimensions: number of fp x 8 interpolation points
  real(field_r), allocatable :: c_weightv_v(:,:)    !< v-interpolation points for forcing points of v
  real(field_r), allocatable :: c_weightv_w(:,:)    !< w-interpolation points for forcing points of v
  real(field_r), allocatable :: c_weightv_s(:,:)    !< scalar interpolation points for forcing points of v

  real(field_r), allocatable :: c_weightw_u(:,:)    !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of w. Dimensions: number of fp x 8 interpolation points
  real(field_r), allocatable :: c_weightw_v(:,:)    !< v-interpolation points for forcing points of w
  real(field_r), allocatable :: c_weightw_w(:,:)    !< w-interpolation points for forcing points of w
  real(field_r), allocatable :: c_weightw_s(:,:)    !< scalar interpolation points for forcing points of w

  real(field_r), allocatable :: c_weights_u(:,:)    !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of scalars. Dimensions: number of fp x 8 interpolation points
  real(field_r), allocatable :: c_weights_v(:,:)    !< v-interpolation points for forcing points of scalars
  real(field_r), allocatable :: c_weights_w(:,:)    !< w-interpolation points for forcing points of scalars
  real(field_r), allocatable :: c_weights_s(:,:)    !< scalar interpolation points for forcing points of scalars

  real(field_r), allocatable :: rot_u(:,:,:)        !< Elements of local terrain roration for forcing points of u. Dimensions: number of fp x 3 x 3 interpolation points
  real(field_r), allocatable :: rot_v(:,:,:)        !< for v
  real(field_r), allocatable :: rot_w(:,:,:)        !< for w
  real(field_r), allocatable :: rot_s(:,:,:)        !< for scalars

  real(field_r), allocatable :: thl_u(:)           !< Temperature of boundary point nearest to forcing points of u 
  real(field_r), allocatable :: thl_v(:)           !< Temperature of boundary point nearest to forcing points of v 
  real(field_r), allocatable :: thl_w(:)           !< Temperature of boundary point nearest to forcing points of w 
  real(field_r), allocatable :: thl_s(:)           !< Temperature of boundary point nearest to forcing points of scalars

  real(field_r), allocatable :: e12wall_l(:)
  real(field_r), allocatable :: qtwall_l(:)
  real(field_r), allocatable :: svwall_l(:,:)
contains
  subroutine initibm
    use modglobal,  only : i1, ih, j1, jh, k1, &
                           zh, zf, itot, jtot, imax, jmax, kmax, &
                           ifnamopt, ifinput, &
                           fname_options, &
                           dx, dy, cu, cv, nsv, &
                           iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
                           iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,iadv_upw,iadv_hybrid,&
                           iadv_hybrid_f, ibas_prf, output_prefix
    use modmpi,     only : myid, comm3d, mpierr, myidx, myidy, d_mpi_bcast, excjs, D_MPI_ALLREDUCE, &
                           mpi_max, mpi_sum, cmyid                       
    implicit none

    integer :: ierr, i,j,k, n, s
    integer :: bcvarr(3)
    integer :: advarr(4)
    character(100) :: name

    !< Read IBM settings from NAMOPTIONS
    namelist/NAMIBM/ lapply_ibm, z0_wall, thlinit, thlwall, e12wall, qtwall, svwall, ibc_thl, ibc_tke, ibc_qt, ibc_sv, n_idw_pts, n_idw_min

    if(myid==0) then !< If on main process, then read NAMOPTIONS settings 
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMIBM,iostat=ierr)

      if (ierr > 0) then
        print *, 'Problem in namoptions NAMIBM'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMIBM'
      endif

      bcvarr = (/ibc_thl,ibc_tke,ibc_qt/)
      if (any(bcvarr > 3).or.any(ibc_sv(1:nsv) > 3)) then
        print *, 'Illegal value for boundary conditions on walls'
        print *, 'Use 1 for Dirichlet, 2 for Neumann, 3 for flux'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMIBM'
      endif

      !< SvdL, 20241009: check use of kappa scheme, ideally tracers should work with this scheme. Also check if other schemes can easily be allowed/implemented
      advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
      if (any(advarr/=2).or.any(iadv_sv(1:nsv)/=2)) then
        print *, 'Current IBM implementation only works with 2nd order advection'
        print *, 'Proper check for kappa advection of scalars to be implemented'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMIBM'
      endif

      write(6 ,NAMIBM)
      close(ifnamopt)

    endif

    !< Broadcast settings to all processes
    call D_MPI_BCAST(lapply_ibm   ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(z0_wall      ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(thlinit      ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(thlwall      ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(qtwall       ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(e12wall      ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(svwall(1:nsv),  nsv,  0, comm3d, mpierr)
    call D_MPI_BCAST(ibc_thl      ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(ibc_qt       ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(ibc_sv(1:nsv),  nsv,  0, comm3d, mpierr)
    call D_MPI_BCAST(n_idw_pts    ,    1,  0, comm3d, mpierr)
    call D_MPI_BCAST(n_idw_min    ,    1,  0, comm3d, mpierr)

    !< Immediately step out of remaining subroutine if ibm is switched off
    if ( .not.(lapply_ibm) ) return

    !< Test for conlicting settings and set constant density (if not done so before)
    if ( lapply_ibm ) then 
      ibas_prf = 2                    !< SvdL, 20240926: actually not sure if required in new ibm setup, retain for now.
      if ( myid==0 ) then
         write (6,*) 'ibas_prf is set to 2 (Boussinesq constant density)'
         write (6,*) 'height dependent density gives problems with correction vertical advective'
         write (6,*) 'tendencies at the tops of obstacles'
      endif
    endif

    if ( abs(cu)>1e-15 .or. abs(cv)>1e-15 ) then
      if( myid==0 ) then
        write (6,*) 'Problem in namoptions'
        write (6,*) 'cu or cv cannot be nonzero while using IBM'
        write (6,*) 'The buildings would move in that case'
        write (6,*) 'Set cu and cv to 0. to solve this problem or simulate without buildings'
        stop 'ERROR: Problem in namoptions NAMIBM with cu and cv'
      endif  
    endif

    !< SvdL, 20240927: few lines below commented out.. will likely not be practical (no method yet to read in full field, 
    !< NOT broadcast to all process, yet still have fields per domain)
    !< REVERT to reading directly per process, WHICH requires that SDF fields are already domain-decomposed before the run.. 
    !< ALSO don't want all processes to access same file at the same time (so prevents having 1 field file)
    !< Allocate temporary field for reading of full SDF file; to prevent broadcasting huge fields to all subprocesses
    ! allocate(sdf_tmp  (itot+1,jtot+1,k1))   !< SvdL, 20240926: check later if size is really correct... don't want to create wrong offset

    ! stop 'ERROR: SvdL, still increase nr of ghostcells to minimum 3 for use with IBM, requires extra switch'

    !< Allocate total fields for the Signed Distance Functions per process
    allocate(sdf_u    (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sdf_v    (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sdf_w    (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sdf_s    (2-ih:i1+ih,2-jh:j1+jh,k1))

    !< Allocate arrays for local grid
    allocate(xh       (2-ih:i1+ih))
    allocate(xf       (2-ih:i1+ih))
    allocate(yh       (2-jh:j1+jh))
    allocate(yf       (2-jh:j1+jh))

    !< Allocate all necessary arrays, use as default size: int(imax*jmax*kmax/10).. could be too small, then throw error later
    nfp_max = INT(imax*jmax*kmax/10)
    nib_max = INT(imax*jmax*kmax/10)

    allocate(fp_u   (nfp_max,3))                 !< Indices of u-points to be forced
    allocate(fp_v   (nfp_max,3))                 !< Indices of v-points to be forced
    allocate(fp_w   (nfp_max,3))                 !< Indices of w-points to be forced
    allocate(fp_s   (nfp_max,3))                 !< Indices of scalar points to be forced
  
    allocate(ib_u   (nib_max,3))                 !< Indices of internal u-points (= within boundary)
    allocate(ib_v   (nib_max,3))                 !< Indices of internal v-points
    allocate(ib_w   (nib_max,3))                 !< Indices of internal w-points
    allocate(ib_s   (nib_max,3))                 !< Indices of internal scalar points
  
    allocate(ipu_u   (nfp_max,8,3))              !< Indices of u-interpolation points for forcing points of u >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
    allocate(ipu_v   (nfp_max,8,3))              !< v-interpolation points for forcing points of u
    allocate(ipu_w   (nfp_max,8,3))              !< w-interpolation points for forcing points of u
    allocate(ipu_s   (nfp_max,8,3))              !< scalar interpolation points for forcing points of u
  
    allocate(ipv_u   (nfp_max,8,3))              !< Indices of u-interpolation points for forcing points of v >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
    allocate(ipv_v   (nfp_max,8,3))              !< v-interpolation points for forcing points of v
    allocate(ipv_w   (nfp_max,8,3))              !< w-interpolation points for forcing points of v
    allocate(ipv_s   (nfp_max,8,3))              !< scalar interpolation points for forcing points of v

    allocate(ipw_u   (nfp_max,8,3))              !< Indices of u-interpolation points for forcing points of w >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
    allocate(ipw_v   (nfp_max,8,3))              !< v-interpolation points for forcing points of w
    allocate(ipw_w   (nfp_max,8,3))              !< w-interpolation points for forcing points of w
    allocate(ipw_s   (nfp_max,8,3))              !< scalar interpolation points for forcing points of w

    allocate(ips_u   (nfp_max,8,3))              !< Indices of u-interpolation points for forcing points of scalars >> Dimensions: number of fp x 8 interpolation points x 3 grid indices
    allocate(ips_v   (nfp_max,8,3))              !< v-interpolation points for forcing points of scalars
    allocate(ips_w   (nfp_max,8,3))              !< w-interpolation points for forcing points of scalars
    allocate(ips_s   (nfp_max,8,3))              !< scalar interpolation points for forcing points of scalars

    allocate(loc_dist_u   (nfp_max,8))     !< Nearest locations on boundary of forcing point + their distance, locations of interpolation points + their distance (xb,yb,zb,db,xi,yi,zi,di)
    allocate(loc_dist_v   (nfp_max,8))     !< for v
    allocate(loc_dist_w   (nfp_max,8))     !< for w
    allocate(loc_dist_s   (nfp_max,8))     !< for scalars

    !< SvdL, 20241009: naming scheme here is quite messy and unfortunate at the moment. See to change/improve later
    allocate(c_weightu_u   (nfp_max,8))     !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of u. Dimensions: number of fp x 8 interpolation points
    allocate(c_weightu_v   (nfp_max,8))     !< v-interpolation points for forcing points of u
    allocate(c_weightu_w   (nfp_max,8))     !< w-interpolation points for forcing points of u
    allocate(c_weightu_s   (nfp_max,8))     !< scalar interpolation points for forcing points of u
  
    allocate(c_weightv_u   (nfp_max,8))     !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of v. Dimensions: number of fp x 8 interpolation points
    allocate(c_weightv_v   (nfp_max,8))     !< v-interpolation points for forcing points of v
    allocate(c_weightv_w   (nfp_max,8))     !< w-interpolation points for forcing points of v
    allocate(c_weightv_s   (nfp_max,8))     !< scalar interpolation points for forcing points of v

    allocate(c_weightw_u   (nfp_max,8))     !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of w. Dimensions: number of fp x 8 interpolation points
    allocate(c_weightw_v   (nfp_max,8))     !< v-interpolation points for forcing points of w
    allocate(c_weightw_w   (nfp_max,8))     !< w-interpolation points for forcing points of w
    allocate(c_weightw_s   (nfp_max,8))     !< scalar interpolation points for forcing points of w

    allocate(c_weights_u   (nfp_max,8))     !< Interpolation coefficients (weights) of the u-interpolation points for forcing points of scalars. Dimensions: number of fp x 8 interpolation points
    allocate(c_weights_v   (nfp_max,8))     !< v-interpolation points for forcing points of scalars
    allocate(c_weights_w   (nfp_max,8))     !< w-interpolation points for forcing points of scalars
    allocate(c_weights_s   (nfp_max,8))     !< scalar interpolation points for forcing points of scalars

    allocate(rot_u   (nfp_max,3,3))        !< Elements of local terrain roration for forcing points of u. Dimensions: number of fp x 3 x 3 interpolation points
    allocate(rot_v   (nfp_max,3,3))        !< for v
    allocate(rot_w   (nfp_max,3,3))        !< for w
    allocate(rot_s   (nfp_max,3,3))        !< for scalars
  
    allocate(thl_u   (nfp_max))           !< Temperature of boundary point nearest to forcing points of u 
    allocate(thl_v   (nfp_max))           !< Temperature of boundary point nearest to forcing points of v 
    allocate(thl_w   (nfp_max))           !< Temperature of boundary point nearest to forcing points of w 
    allocate(thl_s   (nfp_max))           !< Temperature of boundary point nearest to forcing points of scalars

    allocate(qtwall_l   (nfp_max))           !< Temperature of boundary point nearest to forcing points of u 
    allocate(e12wall_l  (nfp_max))           !< Temperature of boundary point nearest to forcing points of v 
    allocate(svwall_l   (nfp_max,nsv))           !< Temperature of boundary point nearest to forcing points of w 

    ! thl_u = thlwall, thl_v = thlwall, thl_w = thlwall, thl_s = thlwall !< Set these allocated temperature arrays to thlinit (set to proper value later)
    c_weightu_u = 0. !0.0_field_r, c_weightu_v = 0.0_field_r, c_weightu_w = 0.0_field_r, c_weightu_s = 0.0_field_r
    !c_weightv_u = 0.0_field_r, c_weightv_v = 0.0_field_r, c_weightv_w = 0.0_field_r, c_weightv_s = 0.0_field_r
    !c_weightw_u = 0.0_field_r, c_weightw_v = 0.0_field_r, c_weightw_w = 0.0_field_r, c_weightw_s = 0.0_field_r
    !c_weights_u = 0.0_field_r, c_weights_v = 0.0_field_r, c_weights_w = 0.0_field_r, c_weights_s = 0.0_field_r  
    
    !< SvdL, 20241012: SILLY WAY TO SET THESE ARRAYS, DO NICER LATEr
    do s = 1,nfp_max
      qtwall_l(s) = qtwall
      e12wall_l(s) = e12wall
      do n = 1, nsv
        svwall_l(s, n) = svwall(n)
      enddo
    enddo

    !< Read in SDF of u
    name = 'sdf.V.XXXXYYYY'
    name(5:5) = 'u'
    name(7:14)=cmyid
    if (myid == 0) write(6,*) 'Loading ',name

    open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted', status='old', access='stream') !<< JUST A TEST
      read(ifinput)  (((sdf_u    (i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(ifinput)

    !< Read in SDF of v  
    name = 'sdf.V.XXXXYYYY'
    name(5:5) = 'v'
    name(7:14)=cmyid
    if (myid == 0) write(6,*) 'Loading ',name

    open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted', status='old', access='stream')
      read(ifinput)  (((sdf_v    (i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(ifinput)

    !< Read in SDF of w      
    name = 'sdf.V.XXXXYYYY'
    name(5:5) = 'w'
    name(7:14)=cmyid
    if (myid == 0) write(6,*) 'Loading ',name
  
    open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted', status='old', access='stream')
      read(ifinput)  (((sdf_w    (i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(ifinput)

    !< Read in SDF of s
    name = 'sdf.V.XXXXYYYY'
    name(5:5) = 's'
    name(7:14)=cmyid
    if (myid == 0) write(6,*) 'Loading ',name
  
    open(unit=ifinput,file=trim(output_prefix)//name,form='unformatted', status='old', access='stream')
      read(ifinput)  (((sdf_s    (i,j,k),i=2,i1),j=2,j1),k=1,kmax)
    close(ifinput)

    !< Set ghostcells for subdomain
    call excjs( sdf_u  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( sdf_v  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( sdf_w  , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( sdf_s  , 2,i1,2,j1,1,k1,ih,jh)

    !< Setting up interpolation scheme requires local grid coordinates to be known
    call constructlocalgrid

    !< Now all SDFs are loaded, proceed with identifying forcing points, ib points, interpolation points and boundary points
    call process_sdf(sdf_u, fp_u, ib_u, loc_dist_u, rot_u, nfp_u, nib_u, xh, yf, zf, 1)
    call process_sdf(sdf_v, fp_v, ib_v, loc_dist_v, rot_v, nfp_v, nib_v, xf, yh, zf, 1)
    call process_sdf(sdf_w, fp_w, ib_w, loc_dist_w, rot_w, nfp_w, nib_w, xf, yf, zh, 2)
    call process_sdf(sdf_s, fp_s, ib_s, loc_dist_s, rot_s, nfp_s, nib_s, xf, yf, zf, 1)

    !< SvdL, 20241009: yet to write, func that setups interpolation
    call setupinterpolation(sdf_u, sdf_v, sdf_w, sdf_s, fp_u, loc_dist_u, nfp_u, xf, yf, xh, yh, ipu_u, ipu_v, ipu_w, ipu_s, c_weightu_u, c_weightu_v, c_weightu_w, c_weightu_s)
    call setupinterpolation(sdf_u, sdf_v, sdf_w, sdf_s, fp_v, loc_dist_v, nfp_v, xf, yf, xh, yh, ipv_u, ipv_v, ipv_w, ipv_s, c_weightv_u, c_weightv_v, c_weightv_w, c_weightv_s)
    call setupinterpolation(sdf_u, sdf_v, sdf_w, sdf_s, fp_w, loc_dist_w, nfp_w, xf, yf, xh, yh, ipw_u, ipw_v, ipw_w, ipw_s, c_weightw_u, c_weightw_v, c_weightw_w, c_weightw_s)
    call setupinterpolation(sdf_u, sdf_v, sdf_w, sdf_s, fp_s, loc_dist_s, nfp_s, xf, yf, xh, yh, ips_u, ips_v, ips_w, ips_s, c_weights_u, c_weights_v, c_weights_w, c_weights_s)

    !< Deallocate SDF fields
    deallocate(sdf_u)
    deallocate(sdf_v)
    deallocate(sdf_w)
    deallocate(sdf_s)

    !< Deallocate local grid (not required anymore)
    deallocate(xh)
    deallocate(yh)
    deallocate(xf)
    deallocate(yf)
    
    return

  end subroutine initibm

  !< Routine to calculate and set tendencies for both the forcing points and the internal points on all active variables, except for eddy viscosities (these are set earlier)
  subroutine applyibm
    use modglobal, only : nsv, lmoist
    use modfields, only : up, vp, wp, thlp, e12p, qtp, svp, u0, v0, w0, thl0, e120, qt0, sv0 !< Make tendencies ("p") and fields at current substep ("0") known for use in ibm

    implicit none

    integer :: n

    if (.not. (lapply_ibm)) return

    !< Calculate and set tendencies to enforce IBM via forcing points
    call forcingpoints_u
    call forcingpoints_v
    call forcingpoints_w

    !< THESE MUST BE IMPROVED WHEN IMPLEMENTING NEUMANN AND FLUX CONDITIONS, TAKE VARIABLES INSIDE THE FUNCTION CALL (AS WE WANT TO AVOID RECALCULATING USTAR ALL THE TIME)
    call forcingpoints_s(thl0, thlp, thl_s  , ibc_thl)
    call forcingpoints_s(e120, e12p, e12wall_l, ibc_tke)

    if (lmoist) call forcingpoints_s(qt0, qtp, qtwall_l, ibc_qt)
    if ( nsv > 0 ) then
      do n = 1, nsv
        call forcingpoints_s(sv0(:,:,:,n), svp(:,:,:,n), svwall_l(:,n), ibc_sv(n))
      enddo
    endif

    !< Calculate and set internal points (inside boundary)
    !< these functions are the same and therefore require separate variable input
    call internalpoints(u0, up, ib_u, nib_u, 0.0_field_r)
    call internalpoints(v0, vp, ib_v, nib_v, 0.0_field_r)
    call internalpoints(w0, wp, ib_w, nib_w, 0.0_field_r)
    call internalpoints(thl0, thlp, ib_s, nib_s, thlinit) !< this is actually a reasonable fill value, it doesn't impact any flux (at least for now, 20241009)

    if (lmoist) call internalpoints(qt0, qtp, ib_s, nib_s, 0.0_field_r)
    if ( nsv > 0 ) then
      do n = 1, nsv
        call internalpoints(sv0(:,:,:,n), svp(:,:,:,n), ib_s, nib_s, 0.0_field_r)
      enddo
    endif

    return

  end subroutine applyibm

  subroutine exitibm
    implicit none

    if (.not. (lapply_ibm)) return

    deallocate(fp_u)
    deallocate(fp_v)
    deallocate(fp_w)
    deallocate(fp_s)
  
    deallocate(ib_u)
    deallocate(ib_v)
    deallocate(ib_w)
    deallocate(ib_s)
  
    deallocate(ipu_u)
    deallocate(ipu_v)
    deallocate(ipu_w)
    deallocate(ipu_s)
  
    deallocate(ipv_u)
    deallocate(ipv_v)
    deallocate(ipv_w)
    deallocate(ipv_s)

    deallocate(ipw_u)
    deallocate(ipw_v)
    deallocate(ipw_w)
    deallocate(ipw_s)

    deallocate(ips_u)
    deallocate(ips_v)
    deallocate(ips_w)
    deallocate(ips_s)

    deallocate(loc_dist_u)
    deallocate(loc_dist_v)
    deallocate(loc_dist_w)
    deallocate(loc_dist_s)
    
    deallocate(c_weightu_u)
    deallocate(c_weightu_v)
    deallocate(c_weightu_w)
    deallocate(c_weightu_s)

    deallocate(c_weightv_u)
    deallocate(c_weightv_v)
    deallocate(c_weightv_w)
    deallocate(c_weightv_s)

    deallocate(c_weightw_u)
    deallocate(c_weightw_v)
    deallocate(c_weightw_w)
    deallocate(c_weightw_s)

    deallocate(c_weights_u)
    deallocate(c_weights_v)
    deallocate(c_weights_w)
    deallocate(c_weights_s)
  
    deallocate(rot_u)
    deallocate(rot_v)
    deallocate(rot_w)
    deallocate(rot_s)
  
    deallocate(thl_u)
    deallocate(thl_v)
    deallocate(thl_w)
    deallocate(thl_s)           

    return

  end subroutine exitibm

  subroutine constructlocalgrid
    use modglobal, only: i1, j1, ih, jh, imax, jmax, dx, dy 
    use modmpi,    only: myidx, myidy
    implicit none

    integer :: i,j 

    !< Set "starting" values of subdomain at indices 2
    xh(2) = myidx*imax*dx !< SvdL, 2024109: if I am correct, myidx and myidy start at zero for first subdomains in both directions
    yh(2) = myidy*jmax*dy
    xf(2) = myidx*imax*dx + 0.5_field_r * dx
    yf(2) = myidy*jmax*dy + 0.5_field_r * dy

    do i = 3, i1+ih
      xh(i) = xh(i-1) + dx
      xf(i) = xf(i-1) + dx
    enddo
    do i = 1, 2-ih
      xh(i) = xh(i+1) - dx
      xf(i) = xf(i+1) - dx
    enddo

    do j = 3, j1+jh
      yh(j) = yh(j-1) + dy
      yf(j) = yf(j-1) + dy
    enddo
    do j = 1, 2-jh
      yh(j) = yh(j+1) - dy
      yf(j) = yf(j+1) - dy
    enddo

    return

  end subroutine constructlocalgrid

  subroutine process_sdf(sdf_in, fp_out, ib_out, loc_out, rot_out, nfp_out, nib_out, x_in, y_in, z_in, kstart)
    use modglobal,  only : i1, ih, j1, jh, k1, &
                           zh, zf, itot, jtot, imax, jmax, kmax, &
                           delta, &
                           dx, dy, cu, cv
    use modmpi,     only : myid, comm3d, mpierr,  myidx, myidy, d_mpi_bcast, excjs, D_MPI_ALLREDUCE, &
                           mpi_max, mpi_sum                       
    implicit none

    !< Inputs of process_sdf function
    real(field_r), intent(in)    :: sdf_in(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(in)    :: x_in(i1)
    real(field_r), intent(in)    :: y_in(j1)
    real(field_r), intent(in)    :: z_in(k1)
    integer, intent(in) :: kstart

    !< Outputs of process_sdf function
    integer, intent(inout) :: fp_out(nfp_max,3)
    integer, intent(inout) :: ib_out(nib_max,3)
    integer, intent(inout) :: nfp_out
    integer, intent(inout) :: nib_out
    real(field_r), intent(inout) :: loc_out(nfp_max,8)
    real(field_r), intent(inout) :: rot_out(nfp_max,3,3)

    !< Internally scoped subvariables
    real(field_r) :: dsdx, dsdy, dsdz, norm, faci, scale
    integer       :: i,j,k
    integer       :: nf = 0
    integer       :: ni = 0

    !< Loop over all values of the signed distance fields
    do k = kstart, kmax          !< Either starts at 1 (u,v,s) or 2 (w)
      do j = 2, j1
        do i = 2, i1

          !< Check if either internal (in boundary) or in forcing layer
          if ( sdf_in(i,j,k) < 0.0_field_r ) then

            ni = ni + 1 
            if ( ni > nib_max) then
              stop 'ERROR: too many ib points on this subdomain, increase nib_max OR use more processes'
            endif

            ib_out(ni,1) = i
            ib_out(ni,1) = j
            ib_out(ni,1) = k

          elseif ( (sdf_in(i,j,k) >= 0.0_field_r) .and. (sdf_in(i,j,k) < delta(k)) ) then

            nf = nf + 1
            if ( nf > nfp_max) then
              stop 'ERROR: too many fp points on this subdomain, increase nfp_max OR use more processes'
            endif    

            fp_out(nf,1) = i 
            fp_out(nf,2) = j 
            fp_out(nf,3) = k

            !< Calculate local normal vector (2nd order accuracy)
            dsdx = 0.5_field_r * ( sdf_in(i+1,j,k) - sdf_in(i-1,j,k) ) / dx
            dsdy = 0.5_field_r * ( sdf_in(i,j+1,k) - sdf_in(i,j-1,k) ) / dy
            dsdz = 0.5_field_r * ( sdf_in(i,j,k+1) - sdf_in(i,j,k-1) ) / delta(k)

            norm = (dsdx*dsdx + dsdy*dsdy + dsdz*dsdz) ** 0.5

            dsdx = dsdx/norm
            dsdy = dsdy/norm
            dsdz = dsdz/norm

            if ( sdf_in(i,j,k) < 0.5_field_r*delta(k)) then
              faci = 1.25_field_r * delta(k) / sdf_in(i,j,k)
            else
              faci = 1.75_field_r * delta(k) / sdf_in(i,j,k)
            endif
            
            !< Store corresponding point on boundary and interpolation point (xb,yb,zb,db,xi,yi,zi,di)
            loc_out(nf,1) = x_in(i) - sdf_in(i,j,k) * dsdx
            loc_out(nf,2) = y_in(j) - sdf_in(i,j,k) * dsdy
            loc_out(nf,3) = z_in(k) - sdf_in(i,j,k) * dsdz
            loc_out(nf,4) = sdf_in(i,j,k)
            loc_out(nf,5) = x_in(i) + faci * sdf_in(i,j,k) * dsdx
            loc_out(nf,6) = y_in(j) + faci * sdf_in(i,j,k) * dsdy
            loc_out(nf,7) = z_in(k) + faci * sdf_in(i,j,k) * dsdz
            loc_out(nf,8) = (1.0_field_r + faci) * sdf_in(i,j,k)

            !< Calculate and store values of rotation matrix 
            if ( ( dsdx == 1.0_field_r .and. dsdy == 0.0_field_r ) .and. dsdz == 0.0_field_r ) then
              rot_out(nf,1,1) = 0.0_field_r
              rot_out(nf,1,2) = 0.0_field_r
              rot_out(nf,1,3) = -1.0_field_r
              rot_out(nf,2,1) = 0.0_field_r
              rot_out(nf,2,2) = 1.0_field_r
              rot_out(nf,2,3) = 0.0_field_r
              rot_out(nf,3,1) = 1.0_field_r
              rot_out(nf,3,2) = 0.0_field_r
              rot_out(nf,3,3) = 0.0_field_r
            elseif ( ( dsdx == -1.0_field_r .and. dsdy == 0.0_field_r ) .and. dsdz == 0.0_field_r ) then
              rot_out(nf,1,1) = 0.0_field_r
              rot_out(nf,1,2) = 0.0_field_r
              rot_out(nf,1,3) = 1.0_field_r
              rot_out(nf,2,1) = 0.0_field_r
              rot_out(nf,2,2) = 1.0_field_r
              rot_out(nf,2,3) = 0.0_field_r
              rot_out(nf,3,1) = -1.0_field_r
              rot_out(nf,3,2) = 0.0_field_r
              rot_out(nf,3,3) = 0.0_field_r
            else
              scale = ( dsdy*dsdy + dsdz*dsdz )**0.5

              rot_out(nf,1,1) = ( dsdy*dsdy + dsdz*dsdz ) / scale
              rot_out(nf,1,2) = -dsdx*dsdy / scale
              rot_out(nf,1,3) = -dsdx*dsdz / scale
              rot_out(nf,2,1) = 0.0_field_r 
              rot_out(nf,2,2) = dsdz / scale 
              rot_out(nf,2,3) = - dsdy / scale 
              rot_out(nf,3,1) = dsdx
              rot_out(nf,3,2) = dsdy
              rot_out(nf,3,3) = dsdz
            endif
          endif 
          
        end do
      end do
    end do
    
    nfp_out = nf
    nib_out = ni

    return

  end subroutine process_sdf  

  subroutine setupinterpolation(sdfu_in, sdfv_in, sdfw_in, sdfs_in, fp_in, loc_dist_in, nfp_in, xf_in, yf_in, xh_in, yh_in, ip_u_out, ip_v_out, ip_w_out, ip_s_out, cweight_u_out, cweight_v_out, cweight_w_out, cweight_s_out)
    use modglobal, only : i1, j1, ih, jh, k1, delta, dx, dy, dzf, dzh, zf, zh
    !SvdL, 20241010: later think about including this instead of own crappy sorting implementation; use stdlib_sorting, only : sort_index !< external standard sorting library for modern fortran
    implicit none

    !< Variables going purely in of the function
    real(field_r), intent(in) :: sdfu_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(in) :: sdfv_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(in) :: sdfw_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(in) :: sdfs_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
    integer, intent(in)       :: fp_in      (nfp_max,3)
    real(field_r), intent(in) :: loc_dist_in(nfp_max,8)
    integer, intent(in)       :: nfp_in
    real(field_r), intent(in) :: xf_in(2-ih:i1+ih), xh_in(2-ih:i1+ih), yf_in(2-jh:j1+jh), yh_in(2-jh:j1+jh)

    !< Variables going in/out of the function
    integer, intent(inout)        :: ip_u_out(nfp_max,8,3), ip_w_out(nfp_max,8,3), ip_v_out(nfp_max,8,3), ip_s_out(nfp_max,8,3)
    real(field_r), intent(inout)  :: cweight_u_out(nfp_max,8), cweight_v_out(nfp_max,8), cweight_w_out(nfp_max,8), cweight_s_out(nfp_max,8)

    !< Internally scoped auxiliary variables
    integer             :: i, j, k, di, dj, dk, dk0, dk0h, nn
    integer             :: inf, jnf, inh, jnh, knf, knh, kf, valn, rr ! start of face levels should be one higher, as zh(1) = 0.0
    real(field_r)       :: xi, yi, zi, distance
    real(field_r)       :: cweight_sum = 0.0_field_r, distmax ! for now only use one maximum distance based on delta via dzf

    ! !< Internally scoped arrays to (temporarilty) store potential neighbours (125 is a HUGE search range, should think about shrinking this)
    ! integer, dimension(8,3)     :: u_neighbour_ind, v_neighbour_ind, w_neighbour_ind, s_neighbour_ind
    ! real(field_r), dimension(8) :: u_neighbour_dis, v_neighbour_dis, w_neighbour_dis, s_neighbour_dis

    !< Loop over all points to be forced
    do nn = 1, nfp_in

      xi = loc_dist_in(nn,5) !< coordinates of interpolation point
      yi = loc_dist_in(nn,6)
      zi = loc_dist_in(nn,7)

      kf = fp_in(nn,3)

      distmax = 2.0_field_r * delta(kf)

      !< SvdL, 20241011: Opt for different implementation than before... !!
      !< this works in combination with "faci" in subroutine process_sdf. Together they should kind of make sure any interpolation point is 
      !< approximately surrounded by 8 fluid points of each kind. In such case, finding these fluid locations from nearest integers inf, inh, jng, etc.
      !< becomes a rather trivial matter. Only thing left is to check/confirm that they are indeed not forcings points, and that in total enough (4,5,6?) remain
      !< Method starts by finding the "bottom-left-left" corner (k,j,i) of grid volume in which the point xi,yi,zi must lie.
      !< Subsequently going +1 index sideways in x and y, and +1 in the vertical should give 8 fluid points, here exclude if they are not far enough from boundary
      !< As long as minimum remains, 4 may be too small... better use 6 for now.. we should have enough free fluid points around the interpolation point. 

      !< Start by estimating grid positions of interpolation points
      inf = FLOOR((xi - xf_in(2))/dx) + 2 !< NINT = round to nearest integer; expression -> (xtest - x(istart)) / dx + istart; where istart=2 for DALES
      jnf = FLOOR((yi - yf_in(2))/dy) + 2

      inh = FLOOR((xi - xh_in(2))/dx) + 2
      jnh = FLOOR((yi - yh_in(2))/dy) + 2

      !< Do different procedure for k (in case of stretched grid)
      !< Even when zi < 0 (which should be allowed anyway, knh should default to 1)
      knf = 1
      knh = 1

      k = 1
      do while ( zf(k) < zi )
        knf = k
        k   = k+1
      enddo 

      k = 1
      do while ( zh(k) < zi )
        knh = k
        k   = k+1
      enddo 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !< Do for u-locations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cweight_sum = 0.0_field_r !< reset to zero
      valn = 0 !< local counter for valid neighbours

      do dk = 0, 1
        do dj = 0, 1
          do di = 0, 1
            if ( (sdfu_in(inh+di, jnf+dj, knf+dk) > delta(kf)).and.(sdfu_in(inh+di, jnf+dj, knf+dk) < 4*delta(kf) )) then
              valn = valn + 1 !< increase number of identified valid neighbours
              distance = ( (xi - xh_in(inh+di))**2 + (yi - yf_in(jnf+dj))**2 + (zi - zf(knf+dk))**2 )**(1./2.) !< calculate distance between gridpoint and interpolation point

              !< if distance is zero (to within some small precision), this point gets weight 1 and RESET all earlier ones to zero
              if ( distance < epsdis ) then 
                cweight_u_out(nn,1) = 1.0_field_r
                ip_u_out(nn,1,1) = inh+di
                ip_u_out(nn,1,2) = jnf+dj
                ip_u_out(nn,1,3) = knf+dk

                do rr = 2,8
                  cweight_u_out(nn,rr) = 0.0_field_r
                  ip_u_out(nn,rr,1) = inh+di    !< just a "FillValue", but with weight zero
                  ip_u_out(nn,rr,2) = jnf+dj    !< just a "FillValue"
                  ip_u_out(nn,rr,3) = knf+dk    !< just a "FillValue"
                enddo 
                valn = 100 !< also set valn to high number (condition unnessary in this case)
                goto 765   !< no need anymore to continue any do-loop
              endif
              
              765 continue 
              !< otherwhise start filling arrays with data for interpolation
              ip_u_out(nn,valn,1) = inh+di
              ip_u_out(nn,valn,2) = jnf+dj
              ip_u_out(nn,valn,3) = knf+dk

              !< SvdL, 20241012: CRITICALLY REVIEW DISTMAX I.C.W. TEST ON SDF VALUE ABOVE.. BOTH NEEDED? NEED TO CHANGE?
              cweight_u_out(nn,valn) = ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
              cweight_sum            = cweight_sum + ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
            endif
          enddo
        enddo
      enddo

      !< Let the code give ERROR when too litte neighbours are found (it also has to as ip_x_out structures are unfilled otherwhise)
      if ( valn < n_idw_min ) then
        write (6,* ) knf, knh
        write (6,*) 'too little valid neighbours identified for interpolation.'
        write (6,*) 'try increasing the grid resolution.'
        stop 'ERROR: setupinterpolation for IBM has failed.'
      endif

      !< Finally, normalize weights
      do rr = 1,8
        cweight_u_out(nn,rr) = cweight_u_out(nn,rr) / cweight_sum
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !< Do for v-locations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cweight_sum = 0.0_field_r !< reset to zero
      valn = 0 !< local counter for valid neighbours

      do dk = 0, 1
        do dj = 0, 1
          do di = 0, 1
            if ( (sdfv_in(inf+di, jnh+dj, knf+dk) > delta(kf)).and.(sdfv_in(inf+di, jnh+dj, knf+dk) < 4*delta(kf) )) then
              valn = valn + 1 !< increase number of identified valid neighbours
              distance = ( (xi - xf_in(inf+di))**2 + (yi - yh_in(jnh+dj))**2 + (zi - zf(knf+dk))**2 )**(1./2.) !< calculate distance between gridpoint and interpolation point

              !< if distance is zero (to within some small precision), this point gets weight 1 and RESET all earlier ones to zero
              if ( distance < epsdis ) then 
                cweight_u_out(nn,1) = 1.0_field_r
                ip_v_out(nn,1,1) = inf+di
                ip_v_out(nn,1,2) = jnh+dj
                ip_v_out(nn,1,3) = knf+dk

                do rr = 2,8
                  cweight_v_out(nn,rr) = 0.0_field_r
                  ip_v_out(nn,rr,1) = inf+di    !< just a "FillValue", but with weight zero
                  ip_v_out(nn,rr,2) = jnh+dj    !< just a "FillValue"
                  ip_v_out(nn,rr,3) = knf+dk    !< just a "FillValue"
                enddo 
                valn = 100 !< also set valn to high number (condition unnessary in this case)
                go to 821
              endif
 
              821 continue

              !< otherwhise start filling arrays with data for interpolation
              ip_v_out(nn,valn,1) = inf+di
              ip_v_out(nn,valn,2) = jnh+dj
              ip_v_out(nn,valn,3) = knf+dk

              !< SvdL, 20241012: CRITICALLY REVIEW DISTMAX I.C.W. TEST ON SDF VALUE ABOVE.. BOTH NEEDED? NEED TO CHANGE?
              cweight_v_out(nn,valn) = ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
              cweight_sum            = cweight_sum + ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
            endif
          enddo 
        enddo 
      enddo 

      !< Let the code give ERROR when too litte neighbours are found (it also has to as ip_x_out structures are unfilled otherwhise)
      if ( valn < n_idw_min ) then
        write (6,* ) knf, knh
        write (6,*) 'too little valid neighbours identified for interpolation.'
        write (6,*) 'try increasing the grid resolution.'
        stop 'ERROR: setupinterpolation for IBM has failed.'
      endif      

      !< Finally, normalize weights
      do rr = 1,8
        cweight_v_out(nn,rr) = cweight_v_out(nn,rr) / cweight_sum
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !< Do for w-locations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cweight_sum = 0.0_field_r !< reset to zero
      valn = 0 !< local counter for valid neighbours

      do dk = 0, 1
        do dj = 0, 1
          do di = 0, 1
            if ( (sdfw_in(inf+di, jnf+dj, knh+dk) > delta(kf)).and.(sdfw_in(inf+di, jnf+dj, knh+dk) < 4*delta(kf) )) then
              valn = valn + 1 !< increase number of identified valid neighbours
              distance = ( (xi - xf_in(inf+di))**2 + (yi - yf_in(jnf+dj))**2 + (zi - zh(knh+dk))**2 )**(1./2.) !< calculate distance between gridpoint and interpolation point

              !< if distance is zero (to within some small precision), this point gets weight 1 and RESET all earlier ones to zero
              if ( distance < epsdis ) then 
                cweight_u_out(nn,1) = 1.0_field_r
                ip_w_out(nn,1,1) = inf+di
                ip_w_out(nn,1,2) = jnf+dj
                ip_w_out(nn,1,3) = knh+dk

                do rr = 2,8
                  cweight_w_out(nn,rr) = 0.0_field_r
                  ip_w_out(nn,rr,1) = inf+di    !< just a "FillValue", but with weight zero
                  ip_w_out(nn,rr,2) = jnf+dj    !< just a "FillValue"
                  ip_w_out(nn,rr,3) = knh+dk    !< just a "FillValue"
                enddo 
                valn = 100 !< also set valn to high number (condition unnessary in this case)
                go to 878
              endif
 
              878 continue 

              !< otherwhise start filling arrays with data for interpolation
              ip_w_out(nn,valn,1) = inf+di
              ip_w_out(nn,valn,2) = jnf+dj
              ip_w_out(nn,valn,3) = knh+dk

              !< SvdL, 20241012: CRITICALLY REVIEW DISTMAX I.C.W. TEST ON SDF VALUE ABOVE.. BOTH NEEDED? NEED TO CHANGE?
              cweight_w_out(nn,valn) = ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
              cweight_sum            = cweight_sum + ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
            endif
          enddo 
        enddo 
      enddo 

      ! !< Let the code give ERROR when too litte neighbours are found (it also has to as ip_x_out structures are unfilled otherwhise)
      ! if ( valn < n_idw_min ) then
      !   write (6,* ) knf, knh
      !   write (6,*) 'too little valid neighbours identified for interpolation.'
      !   write (6,*) 'try increasing the grid resolution.'
      !   stop 'ERROR: setupinterpolation for IBM has failed.'
      ! endif      

      !< Finally, normalize weights
      do rr = 1,8
        cweight_w_out(nn,rr) = cweight_w_out(nn,rr) / cweight_sum
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !< Do for s-locations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      cweight_sum = 0.0_field_r !< reset to zero
      valn = 0 !< local counter for valid neighbours

      do dk = 0, 1
        do dj = 0, 1
          do di = 0, 1
            if ( (sdfs_in(inf+di, jnf+dj, knf+dk) > delta(kf)).and.(sdfs_in(inf+di, jnf+dj, knf+dk) < 4*delta(kf) )) then
              valn = valn + 1 !< increase number of identified valid neighbours
              distance = ( (xi - xf_in(inf+di))**2 + (yi - yf_in(jnf+dj))**2 + (zi - zf(knf+dk))**2 )**(1./2.) !< calculate distance between gridpoint and interpolation point

              !< if distance is zero (to within some small precision), this point gets weight 1 and RESET all earlier ones to zero
              if ( distance < epsdis ) then 
                cweight_s_out(nn,1) = 1.0_field_r
                ip_s_out(nn,1,1) = inf+di
                ip_s_out(nn,1,2) = jnf+dj
                ip_s_out(nn,1,3) = knf+dk

                do rr = 2,8
                  cweight_s_out(nn,rr) = 0.0_field_r
                  ip_s_out(nn,rr,1) = inf+di    !< just a "FillValue", but with weight zero
                  ip_s_out(nn,rr,2) = jnf+dj    !< just a "FillValue"
                  ip_s_out(nn,rr,3) = knf+dk    !< just a "FillValue"
                enddo 
                valn = 100 !< also set valn to high number (condition unnessary in this case)
                go to 935
              endif
 
              935 continue

              !< otherwhise start filling arrays with data for interpolation
              ip_s_out(nn,valn,1) = inf+di
              ip_s_out(nn,valn,2) = jnf+dj
              ip_s_out(nn,valn,3) = knf+dk

              !< SvdL, 20241012: CRITICALLY REVIEW DISTMAX I.C.W. TEST ON SDF VALUE ABOVE.. BOTH NEEDED? NEED TO CHANGE?
              cweight_s_out(nn,valn) = ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
              cweight_sum            = cweight_sum + ( max(0.0_field_r, (distmax - distance) ) / (distmax * distance) )**(1./2.)
            endif
          enddo 
        enddo 
      enddo 

      !< Let the code give ERROR when too litte neighbours are found (it also has to as ip_x_out structures are unfilled otherwhise)
      if ( valn < n_idw_min ) then
        write (6,* ) knf, knh
        write (6,*) 'too little valid neighbours identified for interpolation.'
        write (6,*) 'try increasing the grid resolution.'
        stop 'ERROR: setupinterpolation for IBM has failed.'
      endif      

      !< Finally, normalize weights
      do rr = 1,8
        cweight_s_out(nn,rr) = cweight_s_out(nn,rr) / cweight_sum
      enddo

    enddo

    return

  end subroutine setupinterpolation  

  subroutine forcingpoints_u
    use modglobal, only : ih, i1, jh, j1, k1, rk3step, rdt
    use modfields, only : up, vp, wp, u0, v0, w0
    implicit none

    real(field_r) :: uip, vip, wip, uip_la, vip_la, wip_la, ufp_la, vfp_la, wfp_la
    integer :: nn, pp, if, jf, kf, iip, jip, kip
    real    :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !< Loop over all points to be forced
    do nn = 1, nfp_u
      uip = 0.0_field_r
      vip = 0.0_field_r
      wip = 0.0_field_r

      !< 1. Interpolate surrounding neighbours to interpolation point (use the auxiliary velocity, i.e., the intermediate velocity at next substep without pressure forcing)
      do pp = 1, n_idw_pts
        uip = uip + c_weightu_u(nn,pp) * ( u0( ipu_u(nn,pp,1), ipu_u(nn,pp,2), ipu_u(nn,pp,3) ) + rk3coef * up( ipu_u(nn,pp,1), ipu_u(nn,pp,2), ipu_u(nn,pp,3) ) )
        vip = vip + c_weightu_v(nn,pp) * ( v0( ipu_v(nn,pp,1), ipu_v(nn,pp,2), ipu_v(nn,pp,3) ) + rk3coef * vp( ipu_v(nn,pp,1), ipu_v(nn,pp,2), ipu_v(nn,pp,3) ) )
        wip = wip + c_weightu_w(nn,pp) * ( w0( ipu_w(nn,pp,1), ipu_w(nn,pp,2), ipu_w(nn,pp,3) ) + rk3coef * wp( ipu_w(nn,pp,1), ipu_w(nn,pp,2), ipu_w(nn,pp,3) ) )
      enddo

      !< 2. Rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
      uip_la = rot_u(nn,1,1) * uip + rot_u(nn,1,2) * vip + rot_u(nn,1,3) * wip
      vip_la = rot_u(nn,2,1) * uip + rot_u(nn,2,2) * vip + rot_u(nn,2,3) * wip
      wip_la = rot_u(nn,3,1) * uip + rot_u(nn,3,2) * vip + rot_u(nn,3,3) * wip

      !< for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
      !< (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)

      if ( loc_dist_u(nn, 4) > z0_wall ) then
        !< 3. Calculate (locally-aligned) velocity at forcing point
        ufp_la = uip_la * log(loc_dist_u(nn, 4) / z0_wall) / log(loc_dist_u(nn, 8) / z0_wall);
        vfp_la = vip_la * log(loc_dist_u(nn, 4) / z0_wall) / log(loc_dist_u(nn, 8) / z0_wall);
        wfp_la = wip_la * (loc_dist_u(nn, 4) / loc_dist_u(nn, 8))**2

        !< 4. Rotate back to standard grid alginment (only one component is needed here), 
        !< AND overwrite old tendency at forcing point with new one to achieve this.
        up( fp_u(nn,1), fp_u(nn,2), fp_u(nn,3) ) = ( ( rot_u(nn,1,1)*ufp_la + rot_u(nn,2,1)*vfp_la + rot_u(nn,3,1)*wfp_la ) - u0( fp_u(nn,1), fp_u(nn,2), fp_u(nn,3) ) ) * rk3coefi
      else
        up( fp_u(nn,1), fp_u(nn,2), fp_u(nn,3) ) = ( 0.0_field_r - u0( fp_u(nn,1), fp_u(nn,2), fp_u(nn,3) ) ) * rk3coefi
      endif 
    enddo
    
    return

  end subroutine forcingpoints_u

  subroutine forcingpoints_v
    use modglobal, only : ih, i1, jh, j1, k1, rk3step, rdt
    use modfields, only : up, vp, wp, u0, v0, w0
    implicit none

    real(field_r) :: uip, vip, wip, uip_la, vip_la, wip_la, ufp_la, vfp_la, wfp_la
    integer :: nn, pp, if, jf, kf, iip, jip, kip
    real    :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !< Loop over all points to be forced
    do nn = 1, nfp_v
      uip = 0.0_field_r
      vip = 0.0_field_r
      wip = 0.0_field_r

      !< 1. Interpolate surrounding neighbours to interpolation point (use the auxiliary velocity, i.e., the intermediate velocity at next substep without pressure forcing)
      do pp = 1, n_idw_pts
        uip = uip + c_weightv_u(nn,pp) * ( u0( ipv_u(nn,pp,1), ipv_u(nn,pp,2), ipv_u(nn,pp,3) ) + rk3coef * up( ipv_u(nn,pp,1), ipv_u(nn,pp,2), ipv_u(nn,pp,3) ) )
        vip = vip + c_weightv_v(nn,pp) * ( v0( ipv_v(nn,pp,1), ipv_v(nn,pp,2), ipv_v(nn,pp,3) ) + rk3coef * vp( ipv_v(nn,pp,1), ipv_v(nn,pp,2), ipv_v(nn,pp,3) ) )
        wip = wip + c_weightv_w(nn,pp) * ( w0( ipv_w(nn,pp,1), ipv_w(nn,pp,2), ipv_w(nn,pp,3) ) + rk3coef * wp( ipv_w(nn,pp,1), ipv_w(nn,pp,2), ipv_w(nn,pp,3) ) )
      enddo

      !< 2. Rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
      uip_la = rot_v(nn,1,1) * uip + rot_v(nn,1,2) * vip + rot_v(nn,1,3) * wip
      vip_la = rot_v(nn,2,1) * uip + rot_v(nn,2,2) * vip + rot_v(nn,2,3) * wip
      wip_la = rot_v(nn,3,1) * uip + rot_v(nn,3,2) * vip + rot_v(nn,3,3) * wip

      !< for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
      !< (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)

      if ( loc_dist_v(nn, 4) > z0_wall ) then
        !< 3. Calculate (locally-aligned) velocity at forcing point
        ufp_la = uip_la * log(loc_dist_v(nn, 4) / z0_wall) / log(loc_dist_v(nn, 8) / z0_wall);
        vfp_la = vip_la * log(loc_dist_v(nn, 4) / z0_wall) / log(loc_dist_v(nn, 8) / z0_wall);
        wfp_la = wip_la * (loc_dist_v(nn, 4) / loc_dist_v(nn, 8))**2

        !< 4. Rotate back to standard grid alginment (only one component is needed here), 
        !< AND overwrite old tendency at forcing point with new one to achieve this.
        vp( fp_v(nn,1), fp_v(nn,2), fp_v(nn,3) ) = ( ( rot_v(nn,1,2)*ufp_la + rot_v(nn,2,2)*vfp_la + rot_v(nn,3,2)*wfp_la ) - v0( fp_v(nn,1), fp_v(nn,2), fp_v(nn,3) ) ) * rk3coefi
      else
        vp( fp_v(nn,1), fp_v(nn,2), fp_v(nn,3) ) = ( 0.0_field_r - v0( fp_v(nn,1), fp_v(nn,2), fp_v(nn,3) ) ) * rk3coefi
      endif 
    enddo
    
    return

  end subroutine forcingpoints_v

  subroutine forcingpoints_w
    use modglobal, only : ih, i1, jh, j1, k1, rk3step, rdt
    use modfields, only : up, vp, wp, u0, v0, w0
    implicit none

    real(field_r) :: uip, vip, wip, uip_la, vip_la, wip_la, ufp_la, vfp_la, wfp_la
    integer :: nn, pp, if, jf, kf, iip, jip, kip
    real    :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !< Loop over all points to be forced
    do nn = 1, nfp_w
      uip = 0.0_field_r
      vip = 0.0_field_r
      wip = 0.0_field_r

      !< 1. Interpolate surrounding neighbours to interpolation point (use the auxiliary velocity, i.e., the intermediate velocity at next substep without pressure forcing)
      do pp = 1, n_idw_pts
        uip = uip + c_weightw_u(nn,pp) * ( u0( ipw_u(nn,pp,1), ipw_u(nn,pp,2), ipw_u(nn,pp,3) ) + rk3coef * up( ipw_u(nn,pp,1), ipw_u(nn,pp,2), ipw_u(nn,pp,3) ) )
        vip = vip + c_weightw_v(nn,pp) * ( v0( ipw_v(nn,pp,1), ipw_v(nn,pp,2), ipw_v(nn,pp,3) ) + rk3coef * vp( ipw_v(nn,pp,1), ipw_v(nn,pp,2), ipw_v(nn,pp,3) ) )
        wip = wip + c_weightw_w(nn,pp) * ( w0( ipw_w(nn,pp,1), ipw_w(nn,pp,2), ipw_w(nn,pp,3) ) + rk3coef * wp( ipw_w(nn,pp,1), ipw_w(nn,pp,2), ipw_w(nn,pp,3) ) )
      enddo

      !< 2. Rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
      uip_la = rot_w(nn,1,1) * uip + rot_w(nn,1,2) * vip + rot_w(nn,1,3) * wip
      vip_la = rot_w(nn,2,1) * uip + rot_w(nn,2,2) * vip + rot_w(nn,2,3) * wip
      wip_la = rot_w(nn,3,1) * uip + rot_w(nn,3,2) * vip + rot_w(nn,3,3) * wip

      !< for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
      !< (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)

      if ( loc_dist_w(nn, 4) > z0_wall ) then
        !< 3. Calculate (locally-aligned) velocity at forcing point
        ufp_la = uip_la * log(loc_dist_w(nn, 4) / z0_wall) / log(loc_dist_w(nn, 8) / z0_wall);
        vfp_la = vip_la * log(loc_dist_w(nn, 4) / z0_wall) / log(loc_dist_w(nn, 8) / z0_wall);
        wfp_la = wip_la * (loc_dist_w(nn, 4) / loc_dist_w(nn, 8))**2

        !< 4. Rotate back to standard grid alginment (only one component is needed here), 
        !< AND overwrite old tendency at forcing point with new one to achieve this.
        wp( fp_w(nn,1), fp_w(nn,2), fp_w(nn,3) ) = ( ( rot_w(nn,1,3)*ufp_la + rot_w(nn,2,3)*vfp_la + rot_w(nn,3,3)*wfp_la ) - w0( fp_w(nn,1), fp_w(nn,2), fp_w(nn,3) ) ) * rk3coefi
      else
        wp( fp_w(nn,1), fp_w(nn,2), fp_w(nn,3) ) = ( 0.0_field_r - w0( fp_w(nn,1), fp_w(nn,2), fp_w(nn,3) ) ) * rk3coefi
      endif 
    enddo
    
    return

  end subroutine forcingpoints_w

  subroutine forcingpoints_s(a_in, ap_out, bcval, ibc)
    use modglobal, only : ih, i1, jh, j1, k1, rk3step, rdt
    use modfields, only : up, vp, wp, u0, v0, w0
    implicit none

    real(field_r), intent(in)     :: a_in  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(inout)  :: ap_out(2-ih:i1+ih,2-jh:j1+jh,k1)

    real(field_r), intent(in)     :: bcval (nfp_max)
    integer,       intent(in)     :: ibc

    real(field_r) :: sip, uip, vip, wip, uip_la, vip_la, wip_la, ufp_la, vfp_la, wfp_la
    integer :: nn, pp, if, jf, kf, iip, jip, kip
    real    :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    !< Select appropriate boundary condition
    if     ( ibc == ibc_dirichlet ) then
      !< Loop over all forcing points
      do nn = 1, nfp_s
        sip = 0.0_field_r

        !< 1. Interpolate surrounding neighbours to interpolation point
        do pp = 1, n_idw_pts
          sip = sip + c_weights_s(nn,pp) * (a_in( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) + rk3coef * ap_out( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) )
        enddo

        !< 2. Calculate scalar at forcing point and set tendency to achieve this
        if ( loc_dist_s(nn, 4) > z0_wall ) then
          ap_out( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = ( ( ( sip - bcval(nn) ) * log(loc_dist_s(nn, 4) / z0_wall) / log(loc_dist_s(nn, 8) / z0_wall) + bcval(nn) ) -  a_in( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) )  * rk3coefi
        else
          ap_out( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = ( bcval(nn) -  a_in( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) ) * rk3coefi
        endif
      enddo
    elseif ( ibc == ibc_neumann )   then
      stop 'ERROR: Neumann boundary conditions for walls not yet implemented'
    elseif ( ibc == ibc_flux )      then
      stop 'ERROR: flux boundary conditions for walls not yet implemented'
    endif 

    return

  end subroutine forcingpoints_s

  !< SvdL, 20241009: notes for now:
  ! (1) still implement variable z0_wall(n) and pass these
  ! (2) thlwall moet overschreven worden door temp_s (zijn namelijk op zelfde punt gedefinieerd),
  ! denk nog na over temp_s (hernoemen) en veranderd maken als gevolg van bijv. een process
  ! (3) BELANGRIJK: bij gebruik Deardorff schema: als we e120 forceren nabij de wand, dan zou ekm en ek automatisch ook goed moeten gaan, 
  ! in ieder geval zou het inconsistent zijn als we hier andere waardes forcering die niet compatibel zijn. 
  ! Alleen bij gebruik Smagorinsky schema zou dus ekm/ekh aangepast moeten worden, lokaal...

  !< This one is special: it does not "integrate" to auxiliary velocity, that is because eddy viscosities ARE PASSIVE SCALARS. Also, they are set prior to calculation of all other tendencies.
  !< SvdL, 20241009: STILL THINK ABOUT HOW TO DO THIS IN DALES, AS REGULAR VISCOSITY UPDATES ARE PERFORMED ONLY LATER......
  subroutine forcingpoints_ek
    use modglobal, only : ih, i1, jh, j1, k1, fkar, ekmin
    use modfields, only :  up, vp, wp, u0, v0, w0
    use modsubgriddata, only : ekm, ekh, Prandtl, lsmagorinsky 
    implicit none

    ! if (lsmagorinsky .and. lapply_ibm) then

    real(field_r) :: uip, vip, wip, uip_la, vip_la, wip_la, ufp_la, vfp_la, wfp_la, uhor
    integer :: nn, pp, if, jf, kf, iip, jip, kip

    if (lsmagorinsky .and. lapply_ibm) then
      !< Loop over all points to be forced
      do nn = 1, nfp_s
        uip = 0.0_field_r
        vip = 0.0_field_r
        wip = 0.0_field_r

        !< 1. Interpolate surrounding neighbours to interpolation point
        do pp = 1, n_idw_pts
          uip = uip + c_weights_u(nn,pp) * ( u0( ips_u(nn,pp,1), ips_u(nn,pp,2), ips_u(nn,pp,3) ) ) 
          vip = vip + c_weights_v(nn,pp) * ( v0( ips_v(nn,pp,1), ips_v(nn,pp,2), ips_v(nn,pp,3) ) )
          wip = wip + c_weights_w(nn,pp) * ( w0( ips_w(nn,pp,1), ips_w(nn,pp,2), ips_w(nn,pp,3) ) )
        enddo

        !< 2. Rotate velocities to locally align with surface tangent (under the assumption that flow at second layer still aligns)
        uip_la = rot_s(nn,1,1) * uip + rot_s(nn,1,2) * vip + rot_s(nn,1,3) * wip
        vip_la = rot_s(nn,2,1) * uip + rot_s(nn,2,2) * vip + rot_s(nn,2,3) * wip
        wip_la = rot_s(nn,3,1) * uip + rot_s(nn,3,2) * vip + rot_s(nn,3,3) * wip

        uhor   = ( uip_la*uip_la + vip_la*vip_la )**0.5
      
        if ( loc_dist_u(nn, 4) > z0_wall ) then
          !< 3. calculate local shear velocity (ustar) and thereby eddy viscosity at forcing point
          !< for now, (1) neglect flow rotation over height, (2) neglect stability effects (requires "fine enough" grid),
          !< (3) assume both points are in logarithmic layer, and (4) assume zero-valued Dirichlet conditions for momentum (i.e. no-slip condition)
          ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = fkar**2 * uhor * loc_dist_u(nn, 4) / log( loc_dist_u(nn, 8) / z0_wall )
          ekh( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) / Prandtl

        else
          !< SvdL, 20241009: this is probably a "wrong" value but "best guess" we can do at the moment
          ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = fkar**2 * uhor * z0_wall / log( loc_dist_u(nn, 8) / z0_wall )
          ekh( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) / Prandtl
          
        endif 

        ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = max(ekm( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ), ekmin)
        ekh( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ) = max(ekh( fp_s(nn,1), fp_s(nn,2), fp_s(nn,3) ), ekmin)
      enddo
    endif

    return

  end subroutine forcingpoints_ek

  subroutine internalpoints(a_in, ap_out, ibpoints, nibpoints, val)
    use modglobal, only : ih, i1, jh, j1, k1, rk3step, rdt

    implicit none

    real(field_r), intent(in)     :: a_in  (2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(inout)  :: ap_out(2-ih:i1+ih,2-jh:j1+jh,k1)

    integer, intent(in)           :: ibpoints(nib_max,3)
    integer, intent(in)           :: nibpoints
    real(field_r), intent(in)     :: val

    integer   :: iib,jib,kib,nib
    real      :: rk3coef,rk3coefi

    rk3coef = rdt / (4. - dble(rk3step))
    rk3coefi = 1. / rk3coef

    do nib = 1, nibpoints
      iib = ibpoints(nib,1) 
      jib = ibpoints(nib,2) 
      kib = ibpoints(nib,3) 
 
      ap_out(iib,jib,kib) = ( val - a_in(iib,jib,kib) ) * rk3coefi
    enddo

    return

  end subroutine internalpoints 

  !< SvdL, 20241009 Zou dus alleen bij Smagorinsky relevant moeten zijn, anders niet.. bij Deardorff automatisch goed omdat e120 op 0 wordt gezet.
  subroutine internalpoints_ek
    use modsubgriddata, only : ekm, ekh, Prandtl, lsmagorinsky
    implicit none

    integer :: nn

    if (lsmagorinsky .and. lapply_ibm) then

      do nn = 1, nib_s
        ekm( ib_s(nn,1), ib_s(nn,2), ib_s(nn,3) ) = 0.0_field_r
        ekh( ib_s(nn,1), ib_s(nn,2), ib_s(nn,3) ) = 0.0_field_r
      enddo
    endif

    return

  end subroutine internalpoints_ek

  ! subroutine initibm_Stephan
  !   use modglobal,  only : zh, zf, itot, jtot, ih, i1, jh, j1, k1, imax, jmax, kmax, cexpnr, ifnamopt, ifinput, &
  !                         fname_options, nsv, cu, cv, ijtot, &
  !                         iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
  !                         iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,iadv_upw,iadv_hybrid,&
  !                         iadv_hybrid_f,ibas_prf,&
  !                         dx,dy,fkar
  !   use modmpi,     only : myid, comm3d, mpierr,  myidx, myidy, d_mpi_bcast, excjs, D_MPI_ALLREDUCE, &
  !                         mpi_max, mpi_sum
  !   use modfields,   only : ksfc
  !   !use mpi
  !   implicit none

  !   !< Field for the immersed boundary height
  !   real(field_r), allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
  !   integer, allocatable :: Nairl(:)
  !   integer       :: i, j, k, ierr,ii,jj,kk,n  !cstep , kmin
  !   integer       :: advarr(4)
  !   integer       :: ibm_adv_mask_imin, ibm_adv_mask_imax, & !< use 2nd order advection near obstacles
  !                   ibm_adv_mask_kmin, ibm_adv_mask_kmax
  !   integer          :: kibm_maxl
  !   character(100) :: readstring

  !   namelist/NAMIBM/ lapply_ibm, lreadfile_obstacles, &
  !                             lwallheat, &
  !                             thlwall, thlibm, thlroof, qtibm,lpoislast, z0_wall
  !   !if(.not.(myid==0)) return

  !   dx_half = 0.5 * dx
  !   dy_half = 0.5 * dy
  !   Cm_xwall = (fkar/(log(dx_half/z0_wall)))**2  !cstep similar to neutral boundary layer with a log wind profile
  !   Cm_ywall = (fkar/(log(dy_half/z0_wall)))**2
  !   if (lwallheat) then
  !     Cd_xwall = Cm_xwall   !cstep offer possibility to set it to zero (zero heat flux BC)
  !     Cd_ywall = Cm_ywall
  !   else
  !     Cd_xwall = 0.
  !     Cd_ywall = 0.
  !   endif

  !   write(6,*) 'allocating fields in modibm'

  !   allocate(bc_height (itot+1,jtot+1))
  !   allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))
  !   allocate(ibm_adv_mask(2-ih:i1+ih,2-jh:j1+jh,k1))
  !   allocate(Nair(k1))
  !   allocate(Nairl(k1))

  !   allocate(tempsvp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
  !   allocate(tempthlp(2-ih:i1+ih,2-jh:j1+jh,k1))
  !   allocate(tempqtp(2-ih:i1+ih,2-jh:j1+jh,k1))

  !   write(6,*) 'succesfully allocated fields in modibm'

  !   ibm_adv_mask (:,:,:) = 0.
  !   bc_height(:,:) = 0
  !   libm (:,:,:) = .false.
  !   call D_MPI_BCAST(bc_height,(itot+1)*(jtot+1),0,comm3d,mpierr)

  !   kibm_maxl = 0. !cstep  the index of the highest obstacle in the entire domain
  !   do i=2,i1
  !     do j=2,j1
  !     ! do k=1,kmax  !skip zero values which are ground surface points
  !     !   if (zh(k+1).LE.bc_height(i+myidx*imax,j+myidy*jmax)) then  !cstep read in heights rather than indices
  !             !if zh is 1 mm above building height, height is set to dz below (maybe better use zf as a criterion for nicer rounding off)

  !     !     libm (i,j,k) = .true.
  !     !     ksfc (i,j)   = k + 1

  !       do k=1,kmax
  !           if (zf(k).LE.bc_height(i+myidx*imax,j+myidy*jmax)) then  !obstacle height is above mid point of vertical grid
  !             libm (i,j,k) = .true.
  !             ksfc (i,j)   = k + 1   !half (flux) level
  !             if (ksfc(i,j).gt.kibm_maxl) then
  !               kibm_maxl = k
  !             endif
  !             write (6,*) 'libm',i+myidx*imax,j+myidy*jmax,i,j,k,libm(i,j,k),bc_height(i+myidx*imax,j+myidy*jmax),zh(ksfc(i,j))
  !         endif
  !       end do

  !     end do
  !   end do

  !   call excjs( libm  , 2,i1,2,j1,1,k1,ih,jh)
  !   call D_MPI_ALLREDUCE(kibm_maxl,kibm_max,1,MPI_MAX,comm3d,mpierr)


  !   Nair(:) = 0.
  !   Nairl(:) = 0
  !   do i=2,i1
  !     do j=2,j1
  !       do k=ksfc(i,j),k1
  !         Nairl(k) = Nairl(k)+1
  !       enddo
  !     enddo
  !   enddo
  !   !cstep write (6,*) 'Nairl ',Nairl
  !   call D_MPI_ALLREDUCE(Nairl, Nair, k1, MPI_SUM, comm3d,mpierr)
  !   !cstep write (6,*) 'Nair',Nair

  !   call constructboundarytypes


  !   advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
  !   if     (any(advarr==iadv_cd6).or.any(iadv_sv(1:nsv)==iadv_cd6)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 3
  !     ibm_adv_mask_kmin = 3
  !     ibm_adv_mask_kmax = 3
  !   elseif (any(advarr==iadv_62).or.any(iadv_sv(1:nsv)==iadv_62)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 3
  !     ibm_adv_mask_kmin = 1
  !     ibm_adv_mask_kmax = 1
  !   elseif (any(advarr==iadv_5th).or.any(iadv_sv(1:nsv)==iadv_5th)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 3
  !     ibm_adv_mask_kmin = 2
  !     ibm_adv_mask_kmax = 2
  !   elseif (any(advarr==iadv_52).or.any(iadv_sv(1:nsv)==iadv_52)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 3
  !     ibm_adv_mask_kmin = 1
  !     ibm_adv_mask_kmax = 1
  !   elseif (any(advarr==iadv_hybrid).or.any(iadv_sv(1:nsv)==iadv_hybrid)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 2
  !     ibm_adv_mask_kmin = 3
  !     ibm_adv_mask_kmax = 2
  !   elseif (any(advarr==iadv_hybrid_f).or.any(iadv_sv(1:nsv)==iadv_hybrid_f)) then
  !     ibm_adv_mask_imin = 3
  !     ibm_adv_mask_imax = 2
  !     ibm_adv_mask_kmin = 3
  !     ibm_adv_mask_kmax = 2
  !   elseif (any(advarr==iadv_kappa).or.any(iadv_sv(1:nsv)==iadv_kappa)) then
  !     ibm_adv_mask_imin = 2
  !     ibm_adv_mask_imax = 1
  !     ibm_adv_mask_kmin = 2
  !     ibm_adv_mask_kmax = 1
  !   elseif (any(advarr==iadv_cd2).or.any(iadv_sv(1:nsv)==iadv_cd2)) then
  !     ibm_adv_mask_imin = 1
  !     ibm_adv_mask_imax = 1
  !     ibm_adv_mask_kmin = 1
  !     ibm_adv_mask_kmax = 1
  !   end if

  !   do k=1,kmax-ibm_adv_mask_kmax
  !     do i=2,imax
  !       do j=2,jmax
  !         if (libm(i,j,k)) then
  !           do ii=i-ibm_adv_mask_imax,i+ibm_adv_mask_imin
  !           ibm_adv_mask(ii,j,k) = 1
  !           enddo
  !           do jj=j-ibm_adv_mask_imax,j+ibm_adv_mask_imin
  !           ibm_adv_mask(i,jj,k) = 1
  !           enddo
  !           do kk=k,k+ibm_adv_mask_kmin
  !           ibm_adv_mask(i,j,kk) = 1
  !           enddo
  !         endif
  !       enddo
  !     enddo
  !   enddo

  !   write(6,* ) 'start deallocate'
  !   deallocate (bc_height,Nairl)
  !   write(6,*) 'exit initibm'

  !   return
  ! end subroutine initibm_Stephan

  !!!! BACKUP MATERIAL

  ! subroutine setupinterpolation(sdfu_in, sdfv_in, sdfw_in, sdfs_in, fp_in, loc_dist_in, npf_in, xf_in, yf_in, zf_in, xh_in, yh_in, zh_in)
  !   use modglobal, only : i1, j1, ih, jh, k1, delta
  !   !SvdL, 20241010: later think about including this instead of own crappy sorting implementation; use stdlib_sorting, only : sort_index !< external standard sorting library for modern fortran
  !   implicit none

  !   !< Variables going in/out of the function
  !   real(field_r), intent(in) :: sdfu_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
  !   real(field_r), intent(in) :: sdfv_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
  !   real(field_r), intent(in) :: sdfw_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
  !   real(field_r), intent(in) :: sdfs_in    (2-ih:i1+ih,2-jh:j1+jh,k1)
  !   integer, intent(in)       :: fp_in      (nfp_max,3)
  !   real(field_r), intent(in) :: loc_dist_in(nfp_max,8)
  !   integer, intent(in)       :: nfp_in
  !   real(field_r), intent(in) :: xf_in(2-ih:i1+ih), xh_in(2-ih:i1+ih), yf_in(2-jh:j1+jh), yh_in(2-jh:j1+jh), zf_in(k1), zh_in(k1)

  !   !< Internally scoped auxiliary variables
  !   integer             :: i, j, k, di, dj, dk, dk0, dk0h
  !   integer             :: inf, jnf, inh, jnh, knf = 1 , knh = 2 ! start of face levels should be one higher, as zh(1) = 0.0
  !   real(field_r)       :: xi, yi, zi
  !   real(field_r)       :: c_weight_sum = 0.0_field_r, distmax ! for now only use one maximum distance based on delta via dzf

  !   !< Internally scoped arrays to (temporarilty) store potential neighbours (125 is a HUGE search range, should think about shrinking this)
  !   integer, dimension(125,3)     :: u_neighbour_ind, v_neighbour_ind, w_neighbour_ind, s_neighbour_ind
  !   real(field_r), dimension(125) :: u_neighbour_dis, v_neighbour_dis, w_neighbour_dis, s_neighbour_dis

  !   !< Loop over all points to be forced
  !   do nn = 1, nfp_in

  !     xi = loc_dist_in(nn,5) !< coordinates of interpolation point
  !     yi = loc_dist_in(nn,6)
  !     zi = loc_dist_in(nn,7)

  !     kf = fp_in      (nn,3)

  !     distmax = 2.0_field_r * delta(kf)

  !     !< Start by estimating grid positions of interpolation points
  !     inf = NINT((xi - xf_in(2))/dx) + 2 !< NINT = round to nearest integer; expression -> (xtest - x(istart)) / dx + istart; where istart=2 for DALES
  !     jnf = NINT((yi - yf_in(2))/dy) + 2

  !     inh = NINT((xi - xh_in(2))/dx) + 2
  !     jnh = NINT((yi - yh_in(2))/dy) + 2

  !     !< Do different procedure for k (in case of stretched grid)
  !     knf = 1
  !     knh = 1

  !     do k = 1,(k1-1) !< .. still check why I let this run only to k1-1
  !       if ( zf(k) > zi ) then
  !         if ( (zf(k) - zi ) < (zf(k+1) - zi) ) then
  !           knf = k
  !         else
  !           knf = k+1
  !         endif
  !         exit
  !       endif
  !     enddo

  !     do k = 1,(k1-1) !< .. still check why I let this run only to k1-1
  !       if ( zh(k) > zi ) then
  !         if ( (zh(k) - zi ) < (zh(k+1) - zi) ) then
  !           knh = k
  !         else
  !           knh = k+1
  !         endif
  !         exit
  !       endif
  !     enddo
      
  !     !< Limit vertical stencil near surface (interpolation locations may not be below surface)
  !     dk0  = max(-2, 1-kn);
  !     dk0h = max(-2, 1-knh);

  !     !< Do for u-locations
  !     c_weight_sum = 0.0_field_r !< reset to zero

  !     !< SvdL, 20241010: NOTE ZELF: kijk kritisch naar deze search range (ook VISUEEL/TEKENEN), mij lijkt dat deze search range best kleiner kan, waarmee je standaard al genoeg punten bereikt...
  !     !< SvdL, 20241011: Opt for different implementation than before... 
  !     !< this works in combination with "faci" in subroutine process_sdf. Together they should kind of make sure any interpolation point is 
  !     !< approximately surrounded by 8 fluid points of each kind. In such case, finding these fluid locations from nearest integers inf, inh, jng, etc.
  !     !< becomes a rather trivial matter. Only thing left is to check/confirm that they are indeed not forcings points, and that in total enough (4,5,6?) remain

  !     do dk = dk0, 2
  !       do dj = -2, 2
  !         do di -2, 2
  !           if ( (sdf_in(inh+di, jnf+dj, knf+dk) > delta(kf)).and. (sdf_in(inh+di, jnf+dj, knf+dk) < 4*delta(kf) )) then
  !           endif
  !         enddo
  !       enddo
  !     enddo

  !     if (sdfu[ijk_test] < Delta || sdfu[ijk_test] > TF(4.)*Delta)
  !     continue;

  !     const TF distance = absolute_distance(xi, yi, zi, xh[inh+di], y[jn+dj], z[kn+dk]);
  !     Neighbour<TF> tmp_neighbour = {inh+di, jn+dj, kn+dk, distance};
  !     u_neighbours.push_back(tmp_neighbour);

  !   enddo

  !   return

  ! end subroutine setupinterpolation 

end module modibm
