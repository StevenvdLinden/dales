!> \file modibmdata.f90
!! SvdL, started 25 september 2024
!! By Michael Koene (MK), TU Delft, section Atmospheric Physics, 28 January 2019
!! cstep: subroutine airslabsum  moved from modmpi to here to avoid mutual dependencies

module modibmdata
  use modprecision, only : field_r
  implicit none
  save

  logical :: lapply_ibm     = .false.         !< Switch to enable immersed boundary method 
  
  integer :: ibc_thl = 0, ibc_tke = 0, ibc_qt = 0, ibc_sv(100) = 0
  integer, parameter :: ibc_dirichlet = 0
  integer, parameter :: ibc_neumann   = 1
  integer, parameter :: ibc_flux      = 2

  real(field_r)    :: thlinit        = 293.            !< Initial temperature at the boundaries
  real(field_r)    :: thlwall        = 293.            !< Temperature at the boundaries
  real(field_r)    :: e12wall        = 0.              !< Value of subgrid tke at the boundaries
  real(field_r)    :: qtwall         = 0.              !< Value of moisture at the boundaries
  real(field_r)    :: svwall(100)    = 0.              !< Scalar value of moisture at the boundaries

  real(field_r)    :: z0_wall       = 0.03            !< compare with 0.03 m for open flat terrain, grass, few isolated obstacles

  ! logical :: lreadfile_sdf  = .false.         !< Switch to read positions and model level height from a file
  ! logical :: lwallheat      = .false.         !< Switch to apply lateral heat flux from buildings
  ! logical :: lpoislast      = .true.          !< Switch to use the Poisson solver after the Immersed boundary method
  !                                             !  .false. will set the order to: ZeroVelocity -> PoissonSolver -> IBM

  ! real    :: thlroof        = 293.            !< Obstacle roof (top) temperature
  ! real    :: qtroof         = 0.              !< Obstacle roof specific humidity
  ! real    :: thlibm         = 293             !< Interior potential temperature of obstacle
  ! real    :: qtibm          = 0.              !< Wall specific humidity for the latent heat flux at the sides and top of the buildings
  !                                             !< In modsurface it will be set to the saturation value (but this needs to be adapted)

  ! !< Boolean for applying IBM
  ! logical, allocatable :: libm (:,:,:)       !< Is true inside a building point
  ! !< Mask field for advection. IBM works with 2nd order advection scheme
  ! integer, allocatable :: ibm_adv_mask (:,:,:)
  ! !< Number of grid points in a slab excluding obstacles, and the number of obstacle points
  ! integer, allocatable :: Nair (:)
  ! integer :: kibm_max                        !< index of vertical layer that contains the highest obstacle        


end module modibmdata
