!> \file modstat_nc.f90
! Write NetCDF statistics for individual columns
!
!  \author Bart van Stratum, KNMI/WUR
!
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
! Copyright 1993-2019 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
module modcolumnstat
use modglobal, only : longint
implicit none

public  :: initcolumnstat, exitcolumnstat, columnstat
save
    logical :: lcolumnstat = .false.    ! Switch column statistics

    ! NetCDF output
    integer, parameter :: nvar = 13
    character(80) :: fname = 'column.i_____j_____.___.nc'
    character(80), dimension(nvar,4) :: ncname
    character(80), dimension(1,4)    :: tncname
    integer, allocatable :: nrec(:), ncid(:)

    ! Statistics timing
    real :: dtav=-1, dtsamp=-1
    integer(kind=longint) :: idtav, idtsamp, itnextav, itnextsamp

    ! Number and location (indices) columns
    integer :: ncolumns_loc
    integer, allocatable :: column_i(:), column_j(:)

    ! Array to sample/average statistics
    real, allocatable :: column_data(:,:,:)
    integer :: nsamples = 0

contains

    subroutine columnstat
        !
        ! Main subroutine (called from program.f90), which
        ! takes care of the time book keeping, and calls the routines
        ! which sample and write the column statistics
        !
        use modglobal, only    : rk3step, timee, dt_lim
        implicit none

        if (.not. lcolumnstat) return
        if (rk3step /= 3)      return

        if (timee == itnextsamp) then
            call sample_columns
            itnextsamp = itnextsamp + idtsamp
        else
            dt_lim = min(dt_lim, itnextsamp-timee)
        end if

        if (timee == itnextav) then
            call write_columns
            itnextav = itnextav + idtav
        else
            dt_lim = min(dt_lim, itnextav-timee)
        end if

    end subroutine columnstat


    subroutine sample_columns
        !
        ! Subroutine which gathers the statistics
        !
        use modglobal, only    : cu, cv, kmax, nsv
        use modfields, only    : um, vm, wm, thlm, qtm, ql0, svm
        use modmicrodata, only : imicro, imicro_bulk, precep
        use modraddata, only   : swd, swu, lwd, lwu
        implicit none

        integer :: n, ii, jj

        do n=1, ncolumns_loc
            ii = column_i(n) + 1
            jj = column_j(n) + 1

            column_data(:,1,n) = column_data(:,1,n) + um  (ii,jj,1:kmax)+cu
            column_data(:,2,n) = column_data(:,2,n) + vm  (ii,jj,1:kmax)+cv
            column_data(:,3,n) = column_data(:,3,n) + wm  (ii,jj,1:kmax)
            column_data(:,4,n) = column_data(:,4,n) + thlm(ii,jj,1:kmax)
            column_data(:,5,n) = column_data(:,5,n) + qtm (ii,jj,1:kmax)
            column_data(:,6,n) = column_data(:,6,n) + ql0 (ii,jj,1:kmax)
            if (nsv>1) then
                column_data(:,7,n) = column_data(:,7,n) + svm(ii,jj,1:kmax,1)
                column_data(:,8,n) = column_data(:,8,n) + svm(ii,jj,1:kmax,2)
            end if
            if (imicro == imicro_bulk) then ! This is very conservative...
                column_data(:,9,n) = column_data(:,9,n) + precep(ii,jj,1:kmax)
            end if
            column_data(:,10,n) = column_data(:,10,n) + swd(ii,jj,1:kmax)
            column_data(:,11,n) = column_data(:,11,n) + swu(ii,jj,1:kmax)
            column_data(:,12,n) = column_data(:,12,n) + lwd(ii,jj,1:kmax)
            column_data(:,13,n) = column_data(:,13,n) + lwu(ii,jj,1:kmax)
        end do

        nsamples = nsamples + 1

    end subroutine sample_columns


    subroutine write_columns
        !
        ! Subroutine which writes the statistics to the NetCDF files
        !
        use modglobal, only  : rtimee, kmax
        use modstat_nc, only : writestat_nc
        implicit none

        integer :: n

        ! Calculate mean
        column_data(:,:,:) = column_data(:,:,:) / nsamples

        ! Write to NetCDF files
        do n=1, ncolumns_loc
            call writestat_nc(ncid(n), 1, tncname, (/rtimee/), nrec(n), .true.)
            call writestat_nc(ncid(n), nvar, ncname, column_data(1:kmax,:,n), nrec(n), kmax)
        end do

        ! Reset field and counter
        column_data(:,:,:) = 0.
        nsamples = 0

    end subroutine write_columns

    subroutine initcolumnstat
        !
        ! Subroutine which initialises the column statistics
        !
        use modglobal, only  : ifnamopt, fname_options
        use modmpi,    only  : mpierr, mpi_logical, my_real, comm3d, myid,&
                               myidx, myidy, nprocx, nprocy, mpi_integer
        use modglobal, only  : xsize, ysize, dx, dy, itot, jtot, cexpnr, kmax,&
                               tres, dt_lim, btime
        use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc
        implicit none

        integer :: ierr, n, i, ncolumns, ii, jj
        integer :: maxcolumns = 100
        real    :: x0_mpi, y0_mpi, x1_mpi, y1_mpi, xsize_mpi, ysize_mpi, itot_mpi, jtot_mpi
        character(5) :: ci, cj

        real, allocatable :: column_x(:), column_y(:)
        integer, allocatable :: column_i_tmp(:), column_j_tmp(:)

        namelist/NAMCOLUMNSTAT/ &
        lcolumnstat, column_x, column_y, dtav, dtsamp

        allocate(column_x(maxcolumns), column_y(maxcolumns))
        column_x(:) = -1
        column_y(:) = -1

        if (myid == 0) then
            open(ifnamopt, file=fname_options, status='old', iostat=ierr)
            read (ifnamopt, NAMCOLUMNSTAT, iostat=ierr)
            if (ierr > 0) then
                print *, 'Problem in namoptions NAMCOLUMNSTAT'
                print *, 'iostat error: ', ierr
                stop 'ERROR: Problem in namoptions NAMCOLUMNSTAT'
            endif
            write(6, NAMCOLUMNSTAT)
            close(ifnamopt)

            ! Check how many columns we actually have
            do n=1, maxcolumns
                if (column_x(n) < 0 .or. column_y(n) < 0) then
                    ncolumns = n-1
                    exit
                end if
            end do
        end if

        ! Broadcast namelist input and some other settings
        call mpi_bcast(lcolumnstat, 1,          mpi_logical, 0, comm3d, mpierr)
        call mpi_bcast(column_x,    maxcolumns, my_real,     0, comm3d, mpierr)
        call mpi_bcast(column_y,    maxcolumns, my_real,     0, comm3d, mpierr)
        call mpi_bcast(ncolumns,    1,          mpi_integer, 0, comm3d, mpierr)
        call mpi_bcast(dtav,        1,          my_real,     0, comm3d, mpierr)
        call mpi_bcast(dtsamp,      1,          my_real,     0, comm3d, mpierr)

        ! Statistics/output timing
        idtav   = dtav  /tres
        idtsamp = dtsamp/tres

        itnextav   = btime + idtav
        itnextsamp = btime + idtsamp

        if (.not. lcolumnstat) return

        dt_lim = min(min(dt_lim, idtav), idtsamp)

        ! Check how many columns live on this MPI task
        ! MPI subdomain properties:
        xsize_mpi = xsize / nprocx
        ysize_mpi = ysize / nprocy

        itot_mpi = itot / nprocx
        jtot_mpi = jtot / nprocy

        x0_mpi = myidx     * xsize_mpi
        x1_mpi = (myidx+1) * xsize_mpi
        y0_mpi = myidy     * ysize_mpi
        y1_mpi = (myidy+1) * ysize_mpi

        allocate(column_i_tmp(ncolumns), column_j_tmp(ncolumns))
        column_i_tmp(:) = -1
        column_j_tmp(:) = -1

        ! Count columns on this MPI task, and check for duplicates
        i = 1
        do n=1, ncolumns
            if (column_x(n) >= x0_mpi .and. column_x(n) < x1_mpi .and. column_y(n) >= y0_mpi .and. column_y(n) < y1_mpi) then
                ! Local index in field (without accounting for ghost cells)
                ii = int(column_x(n)/dx)+1 - myidx*itot_mpi
                jj = int(column_y(n)/dy)+1 - myidy*jtot_mpi

                ! Check for duplicates
                if ( any(column_i_tmp == ii) .and. any(column_j_tmp == jj) ) then
                    print*,'WARNING: duplicate column: x=', column_x(n), ', y=', column_y(n), ', skipping..!'
                else
                    column_i_tmp(i) = ii
                    column_j_tmp(i) = jj
                    i = i+1
                end if
            end if
        end do

        ncolumns_loc = i-1

        allocate(column_i(ncolumns_loc), column_j(ncolumns_loc))
        do n=1, ncolumns_loc
            column_i(n) = column_i_tmp(n)
            column_j(n) = column_j_tmp(n)
        end do

        ! Cleanup!
        deallocate(column_i_tmp, column_j_tmp, column_x, column_y)

        ! Create NetCDF files
        allocate(nrec(ncolumns_loc), ncid(ncolumns_loc))

        ! Define variables and their names / units / ..
        ! Time only:
        call ncinfo(tncname( 1,:), 'time', 'Time', 's', 'time')
        ! Vertical profiles:
        call ncinfo(ncname ( 1,:), 'u', 'West-east velocity', 'm/s', 'tt')
        call ncinfo(ncname ( 2,:), 'v', 'South-North velocity', 'm/s', 'tt')
        call ncinfo(ncname ( 3,:), 'w', 'Vertical velocity', 'm/s', 'mt')
        call ncinfo(ncname ( 4,:), 'thl','Liquid water potential temperature','K', 'tt')
        call ncinfo(ncname ( 5,:), 'qt', 'Specific humidity','kg/kg', 'tt')
        call ncinfo(ncname ( 6,:), 'ql', 'Liquid water specific humidity','kg/kg', 'tt')
        call ncinfo(ncname ( 7,:), 'sv001', 'Scalar 001 specific mixing ratio','kg/kg', 'tt')
        call ncinfo(ncname ( 8,:), 'sv002', 'Scalar 002 specific mixing ratio','kg/kg', 'tt')
        call ncinfo(ncname ( 9,:), 'rainrate', 'Rain rate','kg/kg m/s', 'mt')
        call ncinfo(ncname (10,:), 'swd', 'Short wave downward radiative flux','W/m2', 'tt')
        call ncinfo(ncname (11,:), 'swu', 'Short wave upward radiative flux','W/m2', 'tt')
        call ncinfo(ncname (12,:), 'lwd', 'Long wave downward radiative flux','W/m2', 'tt')
        call ncinfo(ncname (13,:), 'lwu', 'Long wave upward radiative flux','W/m2', 'tt')

        ! One NetCDF file per column (for now..)
        do n=1, ncolumns_loc
            ! Global i,j indices & file name
            ii = column_i(n) + itot_mpi*myidx
            jj = column_j(n) + jtot_mpi*myidy
            write(ci,'(i5.5)') ii
            write(cj,'(i5.5)') jj
            fname(21:23) = cexpnr
            fname( 9:13) = ci
            fname(15:19) = cj

            ! Define or open NetCDF file
            ncid(n) = 666 + n
            call open_nc(fname, ncid(n), nrec(n), n3=kmax)
            if (nrec(n) == 0) then
                call define_nc(ncid(n), 1, tncname)
                call writestat_dims_nc(ncid(n))
            end if
            call define_nc(ncid(n), nvar, ncname)
        end do

        ! Allocate array which gathers all statistics
        allocate( column_data(1:kmax, nvar, ncolumns_loc) )
        column_data(:,:,:) = 0.

    end subroutine initcolumnstat


    subroutine exitcolumnstat
        !
        ! Cleanup of the column statistics at end of model run
        !
        use modstat_nc, only : exitstat_nc
        implicit none
        integer :: n

        if (.not. lcolumnstat) return

        ! Close all NetCDF files
        do n=1, ncolumns_loc
            call exitstat_nc( ncid(n) )
        end do

         ! Clean memory
        deallocate(nrec, ncid, column_i, column_j, column_data)

    end subroutine exitcolumnstat

end module modcolumnstat
