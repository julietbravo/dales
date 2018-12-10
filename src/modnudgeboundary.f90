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
    logical :: lnudge_boundary = .false., lperturb_boundary = .false., lsorbjan=.false.
    integer :: nudge_mode = 2  ! 1=to initial profile, 2=to mean profile, 3=input
    real, dimension(:), allocatable :: nudgefac_west,  nudgefac_east, nudgefac_south,  nudgefac_north
    real, dimension(:), allocatable :: perturbfac_west,  perturbfac_east, perturbfac_south,  perturbfac_north
    real, dimension(:), allocatable :: unudge, vnudge, thlnudge, qtnudge                        ! For nudging to profile
    real, dimension(:,:,:,:), allocatable :: unudge_inp, vnudge_inp, thlnudge_inp, qtnudge_inp  ! Input for nudging to external field
    real :: nudge_offset=-1, nudge_width=-1, tau=-1, perturb_ampl=0, zmax_perturb=0, dt_input_lbc=-1
    integer :: blocksize=1, kmax_perturb=0, lbc_index=1


contains
    subroutine initnudgeboundary
        use modmpi,    only : myid, mpierr, comm3d, mpi_logical, mpi_int, my_real, myidx, myidy, nprocx, nprocy
        use modglobal, only : ifnamopt, fname_options, i1, j1, k1, ih, jh, dx, dy, xsize, ysize, zf, kmax
        use modboundary, only : boundary
        implicit none

        integer :: ierr, i, j, k
        real :: x, y

        ! Read namelist settings
        namelist /NAMNUDGEBOUNDARY/ lnudge_boundary, nudge_offset, nudge_width, tau, nudge_mode, &
            & lperturb_boundary, perturb_ampl, blocksize, zmax_perturb, lsorbjan, dt_input_lbc

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
        call MPI_BCAST(lsorbjan,          1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(nudge_mode,        1, mpi_int,     0, comm3d, mpierr)
        call MPI_BCAST(blocksize,         1, mpi_int,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_offset,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_width,       1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(tau,               1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_ampl,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(zmax_perturb,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(dt_input_lbc,      1, my_real,     0, comm3d, mpierr)

        if (lnudge_boundary) then
            ! Init and calculate nudge and perturb factors
            allocate( nudgefac_west(2:i1),  nudgefac_east(2:i1) )
            allocate( nudgefac_south(2:j1),  nudgefac_north(2:j1) )
            allocate( perturbfac_west(2:i1),  perturbfac_east(2:i1) )
            allocate( perturbfac_south(2:i1),  perturbfac_north(2:i1) )

            ! Nudge boundaries to single profile
            if (nudge_mode < 3) then
                allocate( unudge(k1), vnudge(k1), thlnudge(k1), qtnudge(k1) )
            end if

            ! Nudge boundaries to input (3D) data
            if (nudge_mode == 3) then
                ! Two time steps (last dim) are kept in memory,
                ! and are linearly interpolated in time
                allocate(unudge_inp  (2-ih:i1+ih, 2-jh:j1+jh, k1, 2))
                allocate(vnudge_inp  (2-ih:i1+ih, 2-jh:j1+jh, k1, 2))
                allocate(thlnudge_inp(2-ih:i1+ih, 2-jh:j1+jh, k1, 2))
                allocate(qtnudge_inp (2-ih:i1+ih, 2-jh:j1+jh, k1, 2))

                ! Read the first two input times
                call read_new_LBCs(0.)
                call read_new_LBCs(dt_input_lbc)
                lbc_index = 1

                ! Hack - read full initial 3D field
                call read_initial_fields

                ! Make sure the ghost cells are set correctly..
                call boundary
            end if

            ! Calculate the nudge factors
            do i=2,i1
                x = myidx * (xsize / nprocx) + (i-1.5)*dx
                nudgefac_west(i)   = exp(-0.5*((x-       nudge_offset )/nudge_width)**2)
                nudgefac_east(i)   = exp(-0.5*((x-(xsize-nudge_offset))/nudge_width)**2)
            end do
            do j=2,j1
                y = myidy * (ysize / nprocy) + (j-1.5)*dy
                nudgefac_south(j)  = exp(-0.5*((y-       nudge_offset )/nudge_width)**2)
                nudgefac_north(j)  = exp(-0.5*((y-(ysize-nudge_offset))/nudge_width)**2)
            end do

            if (lperturb_boundary) then
                ! Find maximum grid level to which the perturbations are applied
                do k=1,kmax
                    if (zf(k) > zmax_perturb) then
                        kmax_perturb = k-1
                        exit
                    end if
                end do

                ! Calculate the perturbation factors
                do i=2,i1
                    x = myidx * (xsize / nprocx) + (i-1.5)*dx
                    perturbfac_west(i) = exp(-0.5*((x-       2*nudge_offset )/nudge_width)**2)
                    perturbfac_east(i) = exp(-0.5*((x-(xsize-2*nudge_offset))/nudge_width)**2)
                end do
                do j=2,j1
                    y = myidy * (ysize / nprocy) + (j-1.5)*dy
                    perturbfac_west(j) = exp(-0.5*((y-       2*nudge_offset )/nudge_width)**2)
                    perturbfac_east(j) = exp(-0.5*((y-(ysize-2*nudge_offset))/nudge_width)**2)
                end do
            end if

        end if ! lnudge_boundary
    end subroutine initnudgeboundary


    subroutine read_initial_fields
        ! BvS - this should really go somewhere else, probably modstartup...
        use modfields, only : u0, v0, um, vm, thlm, thl0, qtm, qt0
        use modglobal, only : i1, j1, k1, iexpnr, itot, jtot, kmax
        use modmpi, only : myidx, myidy, nprocx, nprocy

        implicit none

        character(80) :: input_file = 'lbc000h00m_x___y___.___'
        write(input_file(13:15), '(i3.3)') myidx
        write(input_file(17:19), '(i3.3)') myidy
        write(input_file(21:23), '(i3.3)') iexpnr

        print*,'Reading initial field: ', input_file

        open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')
        read(666) u0  (2:i1,2:j1,1:kmax)
        read(666) v0  (2:i1,2:j1,1:kmax)
        read(666) thl0(2:i1,2:j1,1:kmax)
        read(666) qt0 (2:i1,2:j1,1:kmax)
        close(666)

        um  (2:i1,2:j1,1:kmax) = u0  (2:i1,2:j1,1:kmax)
        vm  (2:i1,2:j1,1:kmax) = v0  (2:i1,2:j1,1:kmax)
        thlm(2:i1,2:j1,1:kmax) = thl0(2:i1,2:j1,1:kmax)
        qtm (2:i1,2:j1,1:kmax) = qt0 (2:i1,2:j1,1:kmax)

    end subroutine read_initial_fields


    subroutine read_new_LBCs(time)
        use modglobal, only : i1, j1, k1, iexpnr, itot, jtot, kmax
        use modmpi, only : myidx, myidy, nprocx, nprocy

        implicit none
        real, intent(in) :: time !< Input: time to read (seconds)
        integer :: ihour, imin, k
        character(80) :: input_file = 'lbc___h__m_x___y___.___'

        ! Only the MPI tasks at the domain edges read the LBCs:
        if (myidx == 0 .or. myidx == nprocx-1 .or. myidy == 0 .or. myidy == nprocy-1) then

            ! File name to read
            ihour = floor(time/3600)
            imin  = floor((time-ihour*3600)/3600.*60.)
            write(input_file( 4: 6), '(i3.3)') ihour
            write(input_file( 8: 9), '(i2.2)') imin
            write(input_file(13:15), '(i3.3)') myidx
            write(input_file(17:19), '(i3.3)') myidy
            write(input_file(21:23), '(i3.3)') iexpnr
            
            print*,'Processing LBC: ', input_file

            ! Copy t+1 to t
            unudge_inp  (:,:,:,1) = unudge_inp  (:,:,:,2)
            vnudge_inp  (:,:,:,1) = vnudge_inp  (:,:,:,2)
            thlnudge_inp(:,:,:,1) = thlnudge_inp(:,:,:,2)
            qtnudge_inp (:,:,:,1) = qtnudge_inp (:,:,:,2)

            ! Read new LBC for t+1
            open(666, file=input_file, form='unformatted', status='unknown', action='read', access='stream')
            read(666) unudge_inp  (2:i1,2:j1,1:kmax,2)
            read(666) vnudge_inp  (2:i1,2:j1,1:kmax,2)
            read(666) thlnudge_inp(2:i1,2:j1,1:kmax,2)
            read(666) qtnudge_inp (2:i1,2:j1,1:kmax,2)
            close(666)

        end if

    end subroutine read_new_LBCs


    subroutine nudgeboundary
        use modglobal, only : i1, j1, imax, jmax, kmax, rdt, cu, cv, eps1, zf, rtimee
        use modfields, only : u0, up, v0, vp, w0, wp, thl0, thlp, qt0, qtp, &
                            & uprof, vprof, thlprof, qtprof, &
                            & u0av,  v0av,  thl0av,  qt0av
        use modmpi, only    : myidx, myidy, nprocx, nprocy
        implicit none

        integer :: i, j, k, blocki, blockj, subi, subj
        real :: tau_i, perturbation, zi, thetastr, t0, t1, tfac
        real :: unudge_int, vnudge_int, wnudge_int, tnudge_int, qnudge_int

        if (lnudge_boundary) then

            if (tau <= eps1) then
                tau_i = 1. / rdt  ! Nudge on time scale equal to current time step
            else
                tau_i = 1. / tau  ! Nudge on specified time scale
            end if

            ! Select which profile to nudge to
            if (nudge_mode == 1) then
                unudge   = uprof
                vnudge   = vprof
                thlnudge = thlprof
                qtnudge  = qtprof
            else if (nudge_mode == 2) then
                unudge   = u0av
                vnudge   = v0av
                thlnudge = thl0av
                qtnudge  = qt0av
            end if

            ! Nudge boundary to single profile
            !if (nudge_mode < 3) then
            !    do k=1,kmax
            !        do j=2,j1
            !            do i=2,i1
            !                !up(i,j,k)   = up(i,j,k)   + nudgefac_west(i) * tau_i * (unudge(k)   - (u0(i,j,k)+cu)) + &
            !                !                        & + nudgefac_east(i) * tau_i * (unudge(k)   - (u0(i,j,k)+cu))

            !                !vp(i,j,k)   = vp(i,j,k)   + nudgefac_west(i) * tau_i * (vnudge(k)   - (v0(i,j,k)+cv)) + &
            !                !                        & + nudgefac_east(i) * tau_i * (vnudge(k)   - (v0(i,j,k)+cv))

            !                !wp(i,j,k)   = wp(i,j,k)   + nudgefac_west(i) * tau_i * (0.          - w0(i,j,k)) + &
            !                !                        & + nudgefac_east(i) * tau_i * (0.          - w0(i,j,k))

            !                !thlp(i,j,k) = thlp(i,j,k) + nudgefac_west(i) * tau_i * (thlnudge_inp(i,j,k,1) - thl0(i,j,k)) + &
            !                !                        & + nudgefac_east(i) * tau_i * (thlnudge_inp(i,j,k,1) - thl0(i,j,k))

            !                !qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_west(i) * tau_i * (qtnudge(k)  - qt0(i,j,k) ) + &
            !                !                        & + nudgefac_east(i) * tau_i * (qtnudge(k)  - qt0(i,j,k) )
            !            end do
            !        end do
            !    end do
            !end if

            ! Nudge boundary to input data
            if (nudge_mode == 3) then

                ! Read new LBC (if required)
                if (rtimee > lbc_index * dt_input_lbc) then
                    call read_new_LBCs(rtimee)
                    lbc_index = lbc_index + 1
                end if

                ! Calculate time interpolation factor
                t0   = (lbc_index - 1) * dt_input_lbc    ! Time of previous boundary
                t1   = (lbc_index    ) * dt_input_lbc    ! Time of next boundary
                tfac = 1.-(rtimee - t0) / (t1 - t0)      ! Interpolation factor

                ! Zonal nudging
                if (myidx == 0 .or. myidx == nprocx-1) then
                    do k=1,kmax
                        do j=2,j1
                            do i=2,i1

                                ! Interpolate LBC in time
                                unudge_int = tfac * unudge_inp  (i,j,k,1) + (1.-tfac) * unudge_inp  (i,j,k,2)
                                vnudge_int = tfac * vnudge_inp  (i,j,k,1) + (1.-tfac) * vnudge_inp  (i,j,k,2)
                                tnudge_int = tfac * thlnudge_inp(i,j,k,1) + (1.-tfac) * thlnudge_inp(i,j,k,2)
                                qnudge_int = tfac * qtnudge_inp (i,j,k,1) + (1.-tfac) * qtnudge_inp (i,j,k,2)
                                wnudge_int = 0.


                                up(i,j,k)   = up(i,j,k)   + nudgefac_west(i) * tau_i * (unudge_int - (u0(i,j,k)+cu)) + &
                                                        & + nudgefac_east(i) * tau_i * (unudge_int - (u0(i,j,k)+cu))

                                vp(i,j,k)   = vp(i,j,k)   + nudgefac_west(i) * tau_i * (vnudge_int - (v0(i,j,k)+cv)) + &
                                                        & + nudgefac_east(i) * tau_i * (vnudge_int - (v0(i,j,k)+cv))

                                wp(i,j,k)   = wp(i,j,k)   + nudgefac_west(i) * tau_i * (wnudge_int - w0(i,j,k)) + &
                                                        & + nudgefac_east(i) * tau_i * (wnudge_int - w0(i,j,k))

                                thlp(i,j,k) = thlp(i,j,k) + nudgefac_west(i) * tau_i * (tnudge_int - thl0(i,j,k)) + &
                                                        & + nudgefac_east(i) * tau_i * (tnudge_int - thl0(i,j,k))

                                qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_west(i) * tau_i * (qnudge_int - qt0(i,j,k) ) + &
                                                        & + nudgefac_east(i) * tau_i * (qnudge_int - qt0(i,j,k) )
                            end do
                        end do
                    end do
                end if

                if (myidy == 0 .or. myidy == nprocy-1) then
                    do k=1,kmax
                        do j=2,j1
                            do i=2,i1

                                ! Interpolate LBC in time
                                unudge_int = tfac * unudge_inp(i,j,k,1)   + (1.-tfac) * unudge_inp(i,j,k,2)
                                vnudge_int = tfac * vnudge_inp(i,j,k,1)   + (1.-tfac) * vnudge_inp(i,j,k,2)
                                tnudge_int = tfac * thlnudge_inp(i,j,k,1) + (1.-tfac) * thlnudge_inp(i,j,k,2)
                                qnudge_int = tfac * qtnudge_inp(i,j,k,1)  + (1.-tfac) * qtnudge_inp(i,j,k,2)
                                wnudge_int = 0.

                                up(i,j,k)   = up(i,j,k)   + nudgefac_south(j) * tau_i * (unudge_int - (u0(i,j,k)+cu)) + &
                                                        & + nudgefac_north(j) * tau_i * (unudge_int - (u0(i,j,k)+cu))

                                vp(i,j,k)   = vp(i,j,k)   + nudgefac_south(j) * tau_i * (vnudge_int - (v0(i,j,k)+cv)) + &
                                                        & + nudgefac_north(j) * tau_i * (vnudge_int - (v0(i,j,k)+cv))

                                wp(i,j,k)   = wp(i,j,k)   + nudgefac_south(j) * tau_i * (wnudge_int - w0(i,j,k)) + &
                                                        & + nudgefac_north(j) * tau_i * (wnudge_int - w0(i,j,k))

                                thlp(i,j,k) = thlp(i,j,k) + nudgefac_south(j) * tau_i * (tnudge_int - thl0(i,j,k)) + &
                                                        & + nudgefac_north(j) * tau_i * (tnudge_int - thl0(i,j,k))

                                qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_south(j) * tau_i * (qnudge_int  - qt0(i,j,k) ) + &
                                                        & + nudgefac_north(j) * tau_i * (qnudge_int  - qt0(i,j,k) )
                            end do
                        end do
                    end do
                end if
            end if

            ! BvS; quick-and-dirty test with perturbing the inflow boundary.
            if (lperturb_boundary) then

                do k=1,kmax_perturb
                    if (lsorbjan) then
                        ! BvS; even quicker-and-dirtier test using variance scaling from Sorbjan (1989)
                        thetastr = -8e-3 / 0.28        ! BOMEX
                        zi       = zmax_perturb

                        perturb_ampl = 2*((2*(zf(k)/zi)**(-2/3.) * (1-(zf(k)/zi))**(4/3.) + 0.94*(zf(k)/zi)**(4/3.) * &
                            & (1-(zf(k)/zi))**(-2/3.)) * thetastr**2.)**0.5
                    end if

                    do blockj=0, jmax/blocksize-1
                        do blocki=0, imax/blocksize-1
                            perturbation = perturb_ampl*(rand(0)-0.5)

                            do subj=0, blocksize-1
                                do subi=0, blocksize-1
                                    i = blocki*blocksize + subi + 2
                                    j = blockj*blocksize + subj + 2

                                    thlp(i,j,k) = thlp(i,j,k) + perturbfac_west(i) * perturbation / rdt
                                end do
                            end do

                        end do
                    end do
                end do
           end if

        end if ! lnudge_boundary
    end subroutine nudgeboundary


    subroutine exitnudgeboundary
        implicit none
        if (lnudge_boundary) then
            deallocate( nudgefac_west, nudgefac_east, nudgefac_south, nudgefac_north )
            deallocate( perturbfac_west, perturbfac_east, perturbfac_south, perturbfac_north )

            if (nudge_mode  < 3) deallocate( unudge, vnudge, thlnudge, qtnudge )
            if (nudge_mode == 3) deallocate( unudge_inp, vnudge_inp, thlnudge_inp, qtnudge_inp )
        end if
    end subroutine exitnudgeboundary

end module modnudgeboundary
