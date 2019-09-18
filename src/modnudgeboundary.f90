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
    logical :: lnudge_boundary = .false., lnudge_mean=.false., lperturb_boundary = .false.
    integer :: nudge_mode = 2  ! 1=to initial profile, 2=to mean profile
    real, dimension(:), allocatable :: nudgefac_west, perturbfac_west
    real, dimension(:), allocatable :: unudge, vnudge, thlnudge, qtnudge
    real, dimension(:), allocatable :: ub_westl, vb_westl, thlb_westl, qtb_westl
    real, dimension(:), allocatable :: ub_westg, vb_westg, thlb_westg, qtb_westg
    real, dimension(:), allocatable :: ut_west, vt_west, thlt_west, qtt_west
    real :: nudge_offset=-1, nudge_width=-1, tau=-1
    real :: perturb_offset=-1, perturb_width=-1, perturb_ampl=0, zmax_perturb=0
    integer :: blocksize=1, kmax_perturb=0

contains
    subroutine initnudgeboundary
        use modmpi,    only : myid, mpierr, comm3d, mpi_logical, mpi_int, my_real, myidx, myidy, nprocx, nprocy
        use modglobal, only : ifnamopt, fname_options, i1, j1, k1, dx, dy, xsize, ysize, zf, kmax
        implicit none

        integer :: ierr, i, j, k
        real :: x, y

        !
        ! Read namelist settings
        !
        namelist /NAMNUDGEBOUNDARY/ lnudge_boundary, nudge_offset, nudge_width, tau, nudge_mode, &
            & lperturb_boundary, perturb_offset, perturb_width, perturb_ampl, blocksize, zmax_perturb, &
            & lnudge_mean

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
        call MPI_BCAST(lnudge_mean,       1, mpi_logical, 0, comm3d, mpierr)
        call MPI_BCAST(nudge_mode,        1, mpi_int,     0, comm3d, mpierr)
        call MPI_BCAST(blocksize,         1, mpi_int,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_offset,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(nudge_width,       1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(tau,               1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_offset,    1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_width,     1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(perturb_ampl,      1, my_real,     0, comm3d, mpierr)
        call MPI_BCAST(zmax_perturb,      1, my_real,     0, comm3d, mpierr)

        if (lnudge_boundary) then
            !
            ! Init and calculate nudge factors
            !
            allocate( nudgefac_west(2:i1),  perturbfac_west(2:i1) )
            allocate( unudge(k1), vnudge(k1), thlnudge(k1), qtnudge(k1) )
            allocate( ub_westl(k1), vb_westl(k1), thlb_westl(k1), qtb_westl(k1) )
            allocate( ub_westg(k1), vb_westg(k1), thlb_westg(k1), qtb_westg(k1) )
            allocate( ut_west(k1), vt_west(k1), thlt_west(k1), qtt_west(k1) )

            do i=2,i1
                x = myidx * (xsize / nprocx) + (i-1.5)*dx
                nudgefac_west(i)   = exp(-0.5*((x-nudge_offset )/nudge_width)**2)
            end do

            if (lperturb_boundary) then
                !
                ! Calculate perturbation factor
                !
                do i=2,i1
                    x = myidx * (xsize / nprocx) + (i-1.5)*dx
                    perturbfac_west(i) = exp(-0.5*((x-perturb_offset )/perturb_width)**2)
                end do

                !
                ! Find maximum grid level to which the perturbations are applied
                !
                do k=1,kmax
                    if (zf(k) > zmax_perturb) then
                        kmax_perturb = k-1
                        exit
                    end if
                end do
            end if

        end if ! lnudge_boundary
    end subroutine initnudgeboundary

    subroutine nudgeboundary
        use modglobal, only : i1, j1, imax, jmax, kmax, rdt, cu, cv, eps1, zf
        use modfields, only : u0, up, v0, vp, w0, wp, thl0, thlp, qt0, qtp, &
                            & uprof, vprof, thlprof, qtprof, &
                            & u0av,  v0av,  thl0av,  qt0av
        use modmpi, only    : nprocy, my_real, commcol, mpierr, mpi_sum
        implicit none

        integer :: i, j, k, blocki, blockj, subi, subj, n
        real :: tau_i, perturbation, zi, thetastr

        if (lnudge_boundary) then

            if (tau <= eps1) then
                tau_i = 1. / rdt  ! Nudge on time scale equal to current time step
            else
                tau_i = 1. / tau  ! Nudge on specified time scale
            end if

            !
            ! Switch between different quantities to nudge to
            !
            if (nudge_mode == 1) then       ! Nudge to initial profiles
                unudge   = uprof
                vnudge   = vprof
                thlnudge = thlprof
                qtnudge  = qtprof
            else if (nudge_mode == 2) then  ! Nudge to mean profiles
                unudge   = u0av
                vnudge   = v0av
                thlnudge = thl0av
                qtnudge  = qt0av
            else
                stop "unsupported nudge_mode"
            end if

            if (lnudge_mean) then
                !
                ! Nudge only the mean state
                !

                ub_westl   = 0
                vb_westl   = 0
                thlb_westl = 0
                qtb_westl  = 0

                n = 0
                do i=2,i1
                    if (nudgefac_west(i) > 0.01) then
                        n = n+1
                        do k=1,kmax
                            do j=2,j1
                                ub_westl  (k) = ub_westl  (k) + u0  (i,j,k)+cu
                                vb_westl  (k) = vb_westl  (k) + v0  (i,j,k)+cv
                                thlb_westl(k) = thlb_westl(k) + thl0(i,j,k)
                                qtb_westl (k) = qtb_westl (k) + qt0 (i,j,k)
                            end do
                        end do
                    end if
                end do

                if (n > 0) then

                    ub_westl   = ub_westl   / (n*jmax)
                    vb_westl   = vb_westl   / (n*jmax)
                    thlb_westl = thlb_westl / (n*jmax)
                    qtb_westl  = qtb_westl  / (n*jmax)

                    !
                    ! Calculate MPI mean in y-direction
                    ! TO-DO: write routine in modmpi
                    !
                    call MPI_ALLREDUCE(ub_westl,   ub_westg,   kmax, my_real, mpi_sum, commcol, mpierr)
                    call MPI_ALLREDUCE(vb_westl,   vb_westg,   kmax, my_real, mpi_sum, commcol, mpierr)
                    call MPI_ALLREDUCE(thlb_westl, thlb_westg, kmax, my_real, mpi_sum, commcol, mpierr)
                    call MPI_ALLREDUCE(qtb_westl,  qtb_westg,  kmax, my_real, mpi_sum, commcol, mpierr)

                    ub_westg   = ub_westg   / nprocy
                    vb_westg   = vb_westg   / nprocy
                    thlb_westg = thlb_westg / nprocy
                    qtb_westg  = qtb_westg  / nprocy

                    !
                    ! Calculate mean tendency
                    !
                    ut_west   = tau_i * (unudge   - ub_westg)
                    vt_west   = tau_i * (vnudge   - vb_westg)
                    thlt_west = tau_i * (thlnudge - thlb_westg)
                    qtt_west  = tau_i * (qtnudge  - qtb_westg)

                    !
                    ! Apply tendency
                    !
                    do k=1,kmax
                        do j=2,j1
                            do i=2,i1
                                up(i,j,k)   = up(i,j,k)   + nudgefac_west(i) * ut_west  (k)
                                vp(i,j,k)   = vp(i,j,k)   + nudgefac_west(i) * vt_west  (k)
                                thlp(i,j,k) = thlp(i,j,k) + nudgefac_west(i) * thlt_west(k)
                                qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_west(i) * qtt_west (k)
                            end do
                        end do
                    end do

                end if ! n>0

            else

                do k=1,kmax
                    do j=2,j1
                        do i=2,i1
                            up(i,j,k)   = up(i,j,k)   + nudgefac_west(i) * tau_i * (unudge(k)   - (u0(i,j,k)+cu))
                            vp(i,j,k)   = vp(i,j,k)   + nudgefac_west(i) * tau_i * (vnudge(k)   - (v0(i,j,k)+cv))
                            wp(i,j,k)   = wp(i,j,k)   + nudgefac_west(i) * tau_i * (0.          - w0(i,j,k))
                            thlp(i,j,k) = thlp(i,j,k) + nudgefac_west(i) * tau_i * (thlnudge(k) - thl0(i,j,k))
                            qtp(i,j,k)  = qtp(i,j,k)  + nudgefac_west(i) * tau_i * (qtnudge(k)  - qt0(i,j,k) )
                        end do
                    end do
                end do

            end if

            ! BvS; quick-and-dirty test with perturbing the inflow boundary.
            if (lperturb_boundary) then

                do k=1,kmax_perturb
                    do blockj=0, jmax/blocksize-1
                        do blocki=0, imax/blocksize-1
                            perturbation = perturb_ampl*(rand(0)-0.5)

                            do subj=0, blocksize-1
                                do subi=0, blocksize-1
                                    i = blocki*blocksize + subi + 2
                                    j = blockj*blocksize + subj + 2

                                    if (perturbfac_west(i) > 0.01) then
                                        !thlp(i,j,k) = thlp(i,j,k) + perturbfac_west(i) * perturbation / rdt
                                        thlp(i,j,k) = thlp(i,j,k) + perturbation / rdt
                                    end if
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
            deallocate( nudgefac_west, perturbfac_west )
            deallocate( unudge, vnudge, thlnudge, qtnudge)
        end if
    end subroutine exitnudgeboundary

end module modnudgeboundary
