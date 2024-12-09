!> \file modnudge.f90
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>

!>
!!  Nudges theta_l and q_t profiles to the initial profiles on a timescale tnudgeT
!>
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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

module modnudge
  use modprecision, only: field_r

  implicit none

  private

  save

  ! Switches for enabling/disabling nudging
  logical :: lnudge = .false.
  logical :: lunudge = .false.
  logical :: lvnudge = .false.
  logical :: lwnudge = .false.
  logical :: lthlnudge = .false.
  logical :: lqtnudge = .false.

  ! Nudging profiles
  real(field_r), allocatable :: tnudge(:,:)
  real(field_r), allocatable :: unudge(:,:)
  real(field_r), allocatable :: vnudge(:,:)
  real(field_r), allocatable :: wnudge(:,:)
  real(field_r), allocatable :: thlnudge(:,:)
  real(field_r), allocatable :: qtnudge(:,:)

  ! List of nudging time steps
  real(field_r), allocatable :: timenudge(:)

  ! Nudging timescale
  real(field_r) :: tnudgefac = 1.
  ! Number of nudging time steps
  integer       :: ntnudge = 10000

  public :: initnudge
  public :: nudge
  public :: exitnudge

contains

  subroutine initnudge
    use modmpi,    only: myid, mpierr, comm3d, D_MPI_BCAST
    use modglobal, only: ifnamopt, fname_options, runtime, cexpnr, ifinput, &
                         k1, kmax, checknamelisterror

    integer      :: ierr, k, t
    character(1) :: chmess1
    real, allocatable, dimension(:) :: height

    namelist /NAMNUDGE/ lnudge, tnudgefac

    call D_MPI_BCAST(lnudge, 1, 0, comm3d, mpierr)

    if (.not. lnudge) return

    allocate(tnudge(k1,ntnudge), unudge(k1,ntnudge), vnudge(k1,ntnudge), &
              wnudge(k1,ntnudge), thlnudge(k1,ntnudge), qtnudge(k1,ntnudge))
    allocate(timenudge(0:ntnudge), height(k1))

    tnudge = 0
    unudge = 0
    vnudge = 0
    wnudge = 0
    thlnudge = 0
    qtnudge = 0
    timenudge = 0

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read(ifnamopt, NAMNUDGE, iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMNUDGE')
      write(6, NAMNUDGE)
      close(ifnamopt)
    end if

    if (myid == 0) then
      t = 0
      open(ifinput, file='nudge.inp.'//cexpnr)

      do while (timenudge(t) < runtime)
        t = t + 1
        if (t > ntnudge) then
          write(*, *) "Too many time points in file ", 'nudge.inp.'//cexpnr, &
          &", the limit is ntnudge = ", ntnudge
          stop
        end if

        chmess1 = "#"
        ierr = 1 ! not zero
        !search for the next line consisting of "# time", from there onwards the profiles will be read
        do while (.not. (chmess1 == "#" .and. ierr == 0))
          read(ifinput, *, iostat=ierr) chmess1, timenudge(t)
          if (ierr < 0) then
            stop 'STOP: No time dependend nudging data for end of run'
          end if

        end do
        write(6, *) 'time', timenudge(t)
        write(6, *) ' height    t_nudge    u_nudge    v_nudge    w_nudge    &
        &thl_nudge    qt_nudge'
        do k = 1, kmax
          read(ifinput, *) &
            height(k), &
            tnudge(k,t), &
            unudge(k,t), &
            vnudge(k,t), &
            wnudge(k,t), &
            thlnudge(k,t), &
            qtnudge(k,t)
        end do

        do k = kmax, 1, -1
          write(6, '(f7.1,6e12.4)') &
            height(k), &
            tnudge(k,t), &
            unudge(k,t), &
            vnudge(k,t), &
            wnudge(k,t), &
            thlnudge(k,t), &
            qtnudge(k,t)
        end do
      end do
      close (ifinput)
      tnudge = tnudgefac * tnudge
    end if
    call D_MPI_BCAST(timenudge, ntnudge + 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(tnudge, k1 * ntnudge, 0, comm3d, mpierr)
    call D_MPI_BCAST(unudge, k1 * ntnudge, 0, comm3d, mpierr)
    call D_MPI_BCAST(vnudge, k1 * ntnudge, 0, comm3d, mpierr)
    call D_MPI_BCAST(wnudge, k1 * ntnudge, 0, comm3d, mpierr)
    call D_MPI_BCAST(thlnudge, k1 * ntnudge, 0, comm3d, mpierr)
    call D_MPI_BCAST(qtnudge, k1 * ntnudge, 0, comm3d, mpierr)

    lunudge = any(abs(unudge) > 1e-8)
    lvnudge = any(abs(vnudge) > 1e-8)
    lwnudge = any(abs(wnudge) > 1e-8)
    lthlnudge = any(abs(thlnudge) > 1e-8)
    lqtnudge = any(abs(qtnudge) > 1e-8)

    !$acc enter data copyin(timenudge, tnudge, unudge, vnudge, wnudge, &
                            thlnudge, qtnudge)
  end subroutine initnudge

  subroutine nudge
    use modglobal, only: timee, rtimee, i1, j1, kmax, rdt
    use modfields, only: up, vp, wp, thlp, qtp, u0av, v0av, qt0av, thl0av

    integer :: i, j, k, t
    real    :: dtm, dtp, currtnudge

    if (.not. (lnudge)) return

    if (timee == 0) return

    t = 1
    do while (rtimee > timenudge(t))
      t = t + 1
    end do
    if (rtimee > timenudge(1)) then
      t = t - 1
    end if

    dtm = (rtimee - timenudge(t)) / (timenudge(t + 1) - timenudge(t))
    dtp = (timenudge(t + 1) - rtimee) / (timenudge(t + 1) - timenudge(t))

    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          currtnudge = max(1.0_field_r * rdt, &
                           tnudge(k,t) * dtp + tnudge(k,t + 1) * dtm)
          if (lunudge) up(i,j,k) = up(i,j,k) - &
                (u0av(k) - (unudge(k, t)*dtp + unudge(k,t + 1)*dtm))/currtnudge
          if (lvnudge) vp(i,j,k) = vp(i,j,k) - &
                (v0av(k) - (vnudge(k,t)*dtp + vnudge(k,t + 1)*dtm))/currtnudge
          if (lwnudge) wp(i,j,k) = wp(i,j,k) - &
                         (-(wnudge(k, t)*dtp + wnudge(k,t + 1)*dtm))/currtnudge
          if (lthlnudge) thlp(i,j,k) = thlp(i, j, k) - &
          (thl0av(k) - (thlnudge(k,t)*dtp + thlnudge(k,t + 1)*dtm))/currtnudge
          if (lqtnudge) qtp(i,j,k) = qtp(i,j,k) - &
             (qt0av(k) - (qtnudge(k,t)*dtp + qtnudge(k,t + 1)*dtm))/currtnudge
        end do
      end do
    end do
  end subroutine nudge

  subroutine exitnudge
    deallocate(timenudge)
  end subroutine exitnudge

end module
