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

public  :: initnudgeboundary, nudgeboundary, exitnudgeboundary, write_fields
save
    logical :: lnudge_boundary = .false.    ! Switch boundary nudging
    logical :: lperturb_boundary = .false.  ! Switch perturbation of thl near boundary
    logical :: lnudge_w = .true.            ! Nudge w to zero or not
    logical :: lnudge_spectral = .false.    ! Spectral nudging to large scale field

    real, dimension(:,:), allocatable :: nudge_factor    ! Nudging factor (0-1) near lateral boundaries
    real, dimension(:,:), allocatable :: perturb_factor  ! Perturbing factor (0-1) near lateral boundaries

    real, dimension(:,:,:,:), allocatable :: lbc_u, lbc_v, lbc_thl, lbc_qt  ! Input for nudging to external field
    real, dimension(:,:,:),   allocatable :: inc_u, inc_v, inc_thl, inc_qt  ! Spectral nudging increment

    ! Boundary nudging settings:
    real :: nudge_offset=0, nudge_width=0, nudge_radius=0, nudge_tau=-1

    ! Inflow boundary perturbations settings:
    real :: perturb_offset=0, perturb_width=0, perturb_radius=0, perturb_ampl=0, perturb_zmax=100
    integer :: perturb_blocksize=1, kmax_perturb=0

    ! Spectral nudging
    real :: dt_nudge_spectral=-1

    ! Misc
    real :: dt_input_lbc=-1
    integer :: lbc_index=1

contains


    subroutine initnudgeboundary

        use modmpi,      only : myid, mpierr, comm3d, mpi_logical, mpi_int, my_real
        use modglobal,   only : ifnamopt, fname_options, i1, j1, k1, ih, jh, lwarmstart, &
            & kmax, zf, rtimee
        use modboundary, only : boundary
        implicit none

        integer :: ierr, k

        ! Read namelist settings
        namelist /NAMNUDGEBOUNDARY/ lnudge_boundary, lperturb_boundary, lnudge_w, &
            & nudge_offset, nudge_width, nudge_radius, nudge_tau, &
            & perturb_offset, perturb_width, perturb_radius, perturb_ampl, perturb_zmax, &
            & dt_input_lbc, perturb_blocksize, lnudge_spectral, dt_nudge_spectral

        if (myid==0) then
            open(ifnamopt, file=fname_options, status='old', iostat=ierr)
            read (ifnamopt, NAMNUDGEBOUNDARY, iostat=ierr)
            if (ierr > 0) then
                stop 'ERROR: Problem in namoptions NAMNUDGEBOUNDARY'
            endif
            write(6, NAMNUDGEBOUNDARY)
            close(ifnamopt)
        end if

        call MPI_BCAST(lnudge_boundary,   1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(lperturb_boundary, 1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(lnudge_w,          1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(lnudge_spectral,   1, mpi_logical, 0, comm3d, mpierr)

        call MPI_BCAST(perturb_blocksize, 1, mpi_int,     0, comm3d, mpierr)

        call MPI_BCAST(nudge_offset,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_width,       1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_radius,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_tau,         1, my_real,     0, comm3d, mpierr)

        call MPI_BCAST(perturb_offset,    1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_width,     1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_radius,    1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_ampl,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_zmax,      1, my_real,     0, comm3d, mpierr)

        call MPI_BCAST(dt_input_lbc,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(dt_nudge_spectral, 1, my_real,     0, comm3d, mpierr)

        if (lnudge_boundary) then

            ! Init and calculate nudge and perturb factors
            allocate( nudge_factor(2:i1, 2:j1) )

            ! Two time steps (last dim) are kept in memory,
            ! and are linearly interpolated in time
            allocate( lbc_u  (2-ih:i1+ih, 2-jh:j1+jh, k1, 2) )
            allocate( lbc_v  (2-ih:i1+ih, 2-jh:j1+jh, k1, 2) )
            allocate( lbc_thl(2-ih:i1+ih, 2-jh:j1+jh, k1, 2) )
            allocate( lbc_qt (2-ih:i1+ih, 2-jh:j1+jh, k1, 2) )

            ! Read the first two input times
            lbc_index = rtimee/dt_input_lbc+1
            call read_new_LBCs((lbc_index-1)*dt_input_lbc)
            call read_new_LBCs(lbc_index*dt_input_lbc)

            ! Hack - read full initial 3D field
            if (.not. lwarmstart) then
                call read_fields('lbc', .true.)
            end if

            ! Make sure the ghost cells are set correctly..
            call boundary

            ! Calculate the nudging factor
            call calc_weighting_factor(nudge_factor, nudge_offset, nudge_width, nudge_radius)

            if (lperturb_boundary) then
                allocate( perturb_factor(2:i1, 2:j1) )

                call calc_weighting_factor(perturb_factor, perturb_offset, perturb_width, perturb_radius)

                ! Find maximum grid level to which the perturbations are applied
                do k=1,kmax
                    if (zf(k) > perturb_zmax) then
                        kmax_perturb = k-1
                        exit
                    end if
                end do
            end if ! lperturb_boundary

        end if ! lnudge_boundary

        if (lnudge_spectral .and. lwarmstart) then
            ! Allocate fields which contain the increments
            ! from the spectral nudging
            allocate( inc_u  (2-ih:i1+ih, 2-jh:j1+jh, k1) )
            allocate( inc_v  (2-ih:i1+ih, 2-jh:j1+jh, k1) )
            allocate( inc_thl(2-ih:i1+ih, 2-jh:j1+jh, k1) )
            allocate( inc_qt (2-ih:i1+ih, 2-jh:j1+jh, k1) )

            call read_increments
        end if

    end subroutine initnudgeboundary


    subroutine nudgeboundary

        use modglobal, only : i1, j1, imax, jmax, kmax, rdt, cu, cv, &
            & eps1, rtimee, btime, tres, lwarmstart
        use modfields, only : u0, up, v0, vp, w0, wp, thl0, thlp, qt0, qtp
        use modmpi, only    : myidx, myidy, nprocx, nprocy

!#ifdef __INTEL_COMPILER
!        use ifport
!#endif
        implicit none

        integer :: i, j, k, blocki, blockj, subi, subj
        real :: tau_i, perturbation, t0, t1, tfac, runtime
        real :: lbc_u_int, lbc_v_int, lbc_w_int, lbc_t_int, lbc_q_int

        if (lnudge_boundary) then

            if (nudge_tau <= eps1) then
                tau_i = 1. / rdt  ! Nudge on time scale equal to current time step
            else
                tau_i = 1. / nudge_tau  ! Nudge on specified time scale
            end if

            ! Read new LBC (if required)
            if (rtimee > lbc_index*dt_input_lbc) then
                lbc_index = lbc_index + 1
                call read_new_LBCs(lbc_index*dt_input_lbc)
            end if

            ! Calculate time interpolation factor
            t0   = (lbc_index - 1) * dt_input_lbc    ! Time of previous boundary
            t1   = (lbc_index    ) * dt_input_lbc    ! Time of next boundary
            tfac = 1.-(rtimee - t0) / (t1 - t0)      ! Interpolation factor

            if (myidx == 0 .or. myidx == nprocx-1 .or. myidy ==0 .or. myidy == nprocy-1) then
                do k=1,kmax
                    do j=2,j1
                        do i=2,i1

                            ! Interpolate LBC in time
                            lbc_u_int = tfac * lbc_u  (i,j,k,1) + (1.-tfac) * lbc_u  (i,j,k,2)
                            lbc_v_int = tfac * lbc_v  (i,j,k,1) + (1.-tfac) * lbc_v  (i,j,k,2)
                            lbc_t_int = tfac * lbc_thl(i,j,k,1) + (1.-tfac) * lbc_thl(i,j,k,2)
                            lbc_q_int = tfac * lbc_qt (i,j,k,1) + (1.-tfac) * lbc_qt (i,j,k,2)
                            lbc_w_int = 0.

                            ! Nudge the boundaries
                            up(i,j,k)   = up(i,j,k)   + nudge_factor(i,j) * tau_i * (lbc_u_int - (u0(i,j,k)+cu))
                            vp(i,j,k)   = vp(i,j,k)   + nudge_factor(i,j) * tau_i * (lbc_v_int - (v0(i,j,k)+cv))
                            thlp(i,j,k) = thlp(i,j,k) + nudge_factor(i,j) * tau_i * (lbc_t_int - thl0(i,j,k))
                            qtp(i,j,k)  = qtp(i,j,k)  + nudge_factor(i,j) * tau_i * (lbc_q_int - qt0(i,j,k) )

                            if (lnudge_w) then
                                wp(i,j,k)   = wp(i,j,k)   + nudge_factor(i,j) * tau_i * (lbc_w_int - w0(i,j,k))
                            end if

                        end do
                    end do
                end do
            end if


            ! BvS; quick-and-dirty test with perturbing the inflow boundary.
            if (lperturb_boundary) then

                do k=1,kmax_perturb

                    do blockj=0, jmax/perturb_blocksize-1
                        do blocki=0, imax/perturb_blocksize-1
                            perturbation = perturb_ampl*(rand(0)-0.5)

                            do subj=0, perturb_blocksize-1
                                do subi=0, perturb_blocksize-1
                                    i = blocki*perturb_blocksize + subi + 2
                                    j = blockj*perturb_blocksize + subj + 2

                                    thlp(i,j,k) = thlp(i,j,k) + perturb_factor(i,j) * perturbation / rdt
                                end do
                            end do

                        end do
                    end do

                end do
            end if

        end if ! lnudge_boundary


        if (lnudge_spectral .and. lwarmstart) then
            runtime = rtimee - (btime*tres)

            if (runtime <= dt_nudge_spectral) then
                tau_i = 1./dt_nudge_spectral
                
                do k=1,kmax
                    do j=2,j1
                        do i=2,i1
                            up(i,j,k)   = up(i,j,k)   + inc_u  (i,j,k) * tau_i
                            vp(i,j,k)   = vp(i,j,k)   + inc_v  (i,j,k) * tau_i
                            thlp(i,j,k) = thlp(i,j,k) + inc_thl(i,j,k) * tau_i
                            qtp(i,j,k)  = qtp(i,j,k)  + inc_qt (i,j,k) * tau_i
                        end do
                    end do
                end do

            end if
        end if ! lnudge_spectral

    end subroutine nudgeboundary


    subroutine exitnudgeboundary

        implicit none
        if (lnudge_boundary) then
            deallocate( nudge_factor, lbc_u, lbc_v, lbc_thl, lbc_qt )
            if (lperturb_boundary) deallocate( perturb_factor )
            if (lnudge_spectral) deallocate( inc_u, inc_v, inc_thl, inc_qt )
        end if

    end subroutine exitnudgeboundary


    function corner_factor(x, y, x_center, y_center, radius, width) result(f)

        implicit none
        real, intent(in) :: x, y, x_center, y_center, radius, width
        real :: D, f

        D = sqrt((x-x_center)**2 + (y-y_center)**2) - radius
        f = exp(-0.5*(D/width)**2)

    end function


    subroutine calc_weighting_factor(factor, offset, width, radius)

        use modglobal, only : xsize, ysize, j1, i1, dx, dy
        use modmpi,    only : myidx, myidy, nprocx, nprocy
        implicit none

        real, dimension(2:i1,2:j1), intent(inout) :: factor
        real, intent(in) :: offset, width, radius

        real :: dc, x, y
        integer :: i, j

        ! Total size of corners
        dc = offset + radius

        ! Calculate weighting factor
        do j=2, j1
            do i=2,i1
                ! Location in domain, accounting for MPI
                x = myidx * (xsize / nprocx) + (i-1.5)*dx
                y = myidy * (ysize / nprocy) + (j-1.5)*dy

                ! Smooth corners of nudging zone
                if (y < dc .and. x < dc) then    ! SW-corner
                    factor(i,j) = corner_factor(x, y, dc, dc, radius, width)

                else if (y < dc .and. x > xsize-dc) then    ! SE-corner
                    factor(i,j) = corner_factor(x, y, xsize-dc, dc, radius, width)

                else if (y > ysize-dc .and. x < dc) then    ! NW-corner
                    factor(i,j) = corner_factor(x, y, dc, ysize-dc, radius, width)

                else if (y > ysize-dc .and. x > xsize-dc) then   ! NE-corner
                    factor(i,j) = corner_factor(x, y, xsize-dc, ysize-dc, radius, width)

                else
                    factor(i,j) = min( 1., exp(-0.5*((x-       offset )/width)**2) + &
                                         & exp(-0.5*((x-(xsize-offset))/width)**2) + &
                                         & exp(-0.5*((y-       offset )/width)**2) + &
                                         & exp(-0.5*((y-(ysize-offset))/width)**2) )
                end if
            end do
        end do

    end subroutine calc_weighting_factor


    subroutine read_fields(fld_name, lsurface)

        use modfields, only   : u0, v0, um, vm, thlm, thl0, qtm, qt0
        use modsurfdata, only : tskin
        use modglobal, only   : i1, j1, iexpnr, kmax, rtimee
        use modmpi,    only   : myidx, myidy

        implicit none
        character(3), intent(in) :: fld_name
        logical, intent(in) :: lsurface
        integer :: iminute, ihour
        character(80) :: input_file = '___HHHhMMm_x___y___.___'

        ihour = floor(rtimee/3600.)
        iminute = mod(rtimee, 3600.)/60.

        write(input_file(1:3),   '(A3)')   fld_name
        write(input_file(4:6),   '(i3.3)') ihour
        write(input_file(8:9),   '(i2.2)') iminute
        write(input_file(13:15), '(i3.3)') myidx
        write(input_file(17:19), '(i3.3)') myidy
        write(input_file(21:23), '(i3.3)') iexpnr

        print*,'Reading field: ', input_file

        open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')
        read(666) u0   (2:i1,2:j1,1:kmax)
        read(666) v0   (2:i1,2:j1,1:kmax)
        read(666) thl0 (2:i1,2:j1,1:kmax)
        read(666) qt0  (2:i1,2:j1,1:kmax)
        if (lsurface) then 
            read(666) tskin(2:i1,2:j1)
        endif
        close(666)

        um     (2:i1,2:j1,1:kmax) = u0   (2:i1,2:j1,1:kmax)
        vm     (2:i1,2:j1,1:kmax) = v0   (2:i1,2:j1,1:kmax)
        thlm   (2:i1,2:j1,1:kmax) = thl0 (2:i1,2:j1,1:kmax)
        qtm    (2:i1,2:j1,1:kmax) = qt0  (2:i1,2:j1,1:kmax)

    end subroutine read_fields


    subroutine read_increments

        use modglobal, only   : i1, j1, iexpnr, kmax, rtimee
        use modmpi,    only   : myidx, myidy

        implicit none
        integer :: iminute, ihour
        character(80) :: input_file = 'incHHHhMMm_x___y___.___'

        ihour = floor(rtimee/3600.)
        iminute = mod(rtimee, 3600.)/60.

        write(input_file(4:6),   '(i3.3)') ihour
        write(input_file(8:9),   '(i2.2)') iminute
        write(input_file(13:15), '(i3.3)') myidx
        write(input_file(17:19), '(i3.3)') myidy
        write(input_file(21:23), '(i3.3)') iexpnr

        print*,'Reading increments: ', input_file

        open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')
        read(666) inc_u  (2:i1,2:j1,1:kmax)
        read(666) inc_v  (2:i1,2:j1,1:kmax)
        read(666) inc_thl(2:i1,2:j1,1:kmax)
        read(666) inc_qt (2:i1,2:j1,1:kmax)
        close(666)

    end subroutine read_increments


    subroutine read_new_LBCs(time)

        use modglobal, only : i1, j1, iexpnr, kmax
        use modmpi,    only : myidx, myidy, nprocx, nprocy

        implicit none
        real, intent(in) :: time !< Input: time to read (seconds)
        integer :: ihour, iminute
        character(80) :: input_file = 'lbc___h__m_x___y___.___'

        ! Only the MPI tasks at the domain edges read the LBCs:
        if (myidx == 0 .or. myidx == nprocx-1 .or. myidy == 0 .or. myidy == nprocy-1) then

            ! File name to read
            ihour = floor(time/3600.)
            iminute = mod(time, 3600.)/60.
            write(input_file( 4: 6), '(i3.3)') ihour
            write(input_file( 8: 9), '(i2.2)') iminute
            write(input_file(13:15), '(i3.3)') myidx
            write(input_file(17:19), '(i3.3)') myidy
            write(input_file(21:23), '(i3.3)') iexpnr

            print*,'Processing LBC: ', input_file

            ! Copy old second time to new first time
            lbc_u  (:,:,:,1) = lbc_u  (:,:,:,2)
            lbc_v  (:,:,:,1) = lbc_v  (:,:,:,2)
            lbc_thl(:,:,:,1) = lbc_thl(:,:,:,2)
            lbc_qt (:,:,:,1) = lbc_qt (:,:,:,2)

            ! Read new LBC for next time
            open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')
            read(666) lbc_u  (2:i1,2:j1,1:kmax,2)
            read(666) lbc_v  (2:i1,2:j1,1:kmax,2)
            read(666) lbc_thl(2:i1,2:j1,1:kmax,2)
            read(666) lbc_qt (2:i1,2:j1,1:kmax,2)
            close(666)

        end if

    end subroutine read_new_LBCs


    subroutine write_fields(fld_name)

        ! BvS - this should really go somewhere else, probably modstartup...
        use modfields, only   : u0, v0, um, vm, thlm, thl0, qtm, qt0
        use modsurfdata, only : tskin
        use modglobal, only   : i1, j1, iexpnr, kmax, timeleft, rk3step, rtimee
        use modmpi,    only   : myidx, myidy

        implicit none
        character(3), intent(in) :: fld_name
        integer :: ihour, iminute
        character(80) :: input_file = '___HHHhMMm_x___y___.___'

        if (lnudge_boundary .and. lnudge_spectral .and. rk3step == 3 .and. timeleft==0) then

            ihour = floor(rtimee/3600.)
            iminute = mod(rtimee, 3600.)/60.

            write(input_file(1:3),   '(A3)')   fld_name
            write(input_file(4:6),   '(i3.3)') ihour
            write(input_file(8:9),   '(i2.2)') iminute
            write(input_file(13:15), '(i3.3)') myidx
            write(input_file(17:19), '(i3.3)') myidy
            write(input_file(21:23), '(i3.3)') iexpnr
            print*,'Writing field: ', input_file

            open(666, file=input_file, form='unformatted', status='unknown', action='write', access='stream')
            write(666) u0   (2:i1,2:j1,1:kmax)
            write(666) v0   (2:i1,2:j1,1:kmax)
            write(666) thl0 (2:i1,2:j1,1:kmax)
            write(666) qt0  (2:i1,2:j1,1:kmax)
            close(666)
        end if

    end subroutine write_fields

end module modnudgeboundary
