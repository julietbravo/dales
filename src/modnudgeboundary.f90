!> \file modrelaxboundary.f90
!>
!! Nudge the lateral boundaries
!>
!! \author Bart van Stratum, KNMI
!
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modnudgeboundary
implicit none

public  :: initnudgeboundary, nudgeboundary, exitnudgeboundary
save
    logical :: lnudge_boundary = .false.
    real, dimension(:), allocatable :: nudgefac_west,  nudgefac_east
    real, dimension(:), allocatable :: nudgefac_south, nudgefac_north
    real :: nudge_offset=-1, nudge_width=-1, tau=-1

contains
    subroutine initnudgeboundary
        use modmpi,    only : myid, mpierr, comm3d, mpi_logical, my_real, myidx, myidy, nprocx, nprocy
        use modglobal, only : ifnamopt, fname_options, i1, j1, dx, dy, xsize, ysize
        implicit none

        integer :: ierr, i, j
        real :: x, y

        !
        ! Read namelist settings
        !
        namelist /NAMNUDGEBOUNDARY/ lnudge_boundary, nudge_offset, nudge_width, tau

        if (myid==0) then
            open(ifnamopt, file=fname_options, status='old', iostat=ierr)
            read (ifnamopt, NAMNUDGEBOUNDARY, iostat=ierr)
            if (ierr > 0) then
                stop 'ERROR: Problem in namoptions NAMNUDGEBOUNDARY'
            endif
            write(6, NAMNUDGEBOUNDARY)
            close(ifnamopt)
        end if

        call MPI_BCAST(lnudge_boundary, 1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(nudge_offset,    1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_width,     1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(tau,             1, my_real,     0, comm3d, mpierr)

        if (lnudge_boundary) then
            !
            ! Init and calculate nudge factors
            !
            allocate( nudgefac_west(2:i1),  nudgefac_east(2:i1) )
            allocate( nudgefac_south(2:j1), nudgefac_north(2:j1) )

            do i=2,i1
                x = myidx * (xsize / nprocx) + (i-1.5)*dx
                nudgefac_west(i) = exp(-0.5*((x-       nudge_offset )/nudge_width)**2)
                nudgefac_east(i) = exp(-0.5*((x-(xsize-nudge_offset))/nudge_width)**2)
            end do

            do j=2,j1
                y = myidy * (ysize / nprocy) + (j-1.5)*dy
                nudgefac_south(j) = exp(-0.5*((y-       nudge_offset )/nudge_width)**2)
                nudgefac_north(j) = exp(-0.5*((y-(ysize-nudge_offset))/nudge_width)**2)
            end do
        end if ! lnudge_boundary
    end subroutine initnudgeboundary

    subroutine nudgeboundary
        use modglobal, only : i1, j1, kmax, rdt, cu, cv, eps1
        use modfields, only : u0, up, v0, vp, w0, wp, thl0, thlp, qt0, qtp, uprof, vprof, thlprof, qtprof
        implicit none

        integer :: i,j,k
        real :: tau_i

        if (lnudge_boundary) then
            if (tau <= eps1) then
                tau_i = 1. / rdt  ! Nudge on time scale equal to current time step
            else
                tau_i = 1. / tau  ! Nudge on specified time scale
            end if

            do k=1,kmax
                do j=2,j1
                    do i=2,i1
                        up(i,j,k)   = up(i,j,k)   + nudgefac_west(i) * tau_i * (uprof(k)   - (u0(i,j,k)+cu)) + &
                                                & + nudgefac_east(i) * tau_i * (uprof(k)   - (u0(i,j,k)+cu))

                        vp(i,j,k)   = vp(i,j,k)   + nudgefac_west(i) * tau_i * (vprof(k)   - (v0(i,j,k)+cv)) + &
                                                & + nudgefac_east(i) * tau_i * (vprof(k)   - (v0(i,j,k)+cv))

                        wp(i,j,k)   = wp(i,j,k)   + nudgefac_west(i) * tau_i * (0.         - w0(i,j,k)) + &
                                                & + nudgefac_east(i) * tau_i * (0.         - w0(i,j,k))

                        thlp(i,j,k) = thlp(i,j,k) + nudgefac_west(i) * tau_i * (thlprof(k) - thl0(i,j,k)) + &
                                                & + nudgefac_east(i) * tau_i * (thlprof(k) - thl0(i,j,k))

                        qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_west(i) * tau_i * (qtprof(k)  - qt0(i,j,k) ) + &
                                                & + nudgefac_east(i) * tau_i * (qtprof(k)  - qt0(i,j,k) )
                    end do
                end do
            end do

        end if ! lnudge_boundary
    end subroutine nudgeboundary

    subroutine exitnudgeboundary
        implicit none
        if (lnudge_boundary) then
            deallocate( nudgefac_west, nudgefac_east, nudgefac_south, nudgefac_north )
        end if
    end subroutine exitnudgeboundary

end module modnudgeboundary
