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
implicit none

public  :: initcolumnstat!, columnstat, exitcolumnstat
save
    logical :: lcolumnstat = .false.    ! Switch column statistics
    integer :: ncolumns_loc             ! Number of columns on this MPI task

    ! NetCDF output
    integer, parameter :: nvar = 13
    character(80) :: fname = 'column.i_____j_____.___.nc'
    character(80), dimension(nvar,4) :: ncname
    character(80), dimension(1,4)    :: tncname
    integer, allocatable :: nrec(:), ncid(:)

    ! Location columns
    integer, allocatable :: column_i(:), column_j(:)

contains

subroutine initcolumnstat
    use modglobal, only  : ifnamopt, fname_options
    use modmpi,    only  : mpierr, mpi_logical, my_real, comm3d, myid,&
                           myidx, myidy, nprocx, nprocy, mpi_integer,&
                           cmyidx, cmyidy
    use modglobal, only  : xsize, ysize, dx, dy, itot, jtot, cexpnr, kmax
    use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc
    implicit none

    integer :: ierr, n, i, ncolumns, ii, jj
    integer :: maxcolumns = 100
    real    :: x0_mpi, y0_mpi, x1_mpi, y1_mpi, xsize_mpi, ysize_mpi, itot_mpi, jtot_mpi
    character(5) :: ci, cj

    real, allocatable :: column_x(:), column_y(:)
    integer, allocatable :: column_i_tmp(:), column_j_tmp(:)

    namelist/NAMCOLUMNSTAT/ &
    lcolumnstat, column_x, column_y

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

    if (lcolumnstat) then

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
                ! Local index in field
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
        call ncinfo(tncname( 1,:), 'time', 'Time', 's', 'time')
        call ncinfo(ncname ( 1,:), 'u', 'West-east velocity', 'm/s', 'tt')
        call ncinfo(ncname ( 2,:), 'v', 'South-North velocity', 'm/s', 'tt')
        call ncinfo(ncname ( 3,:), 'w', 'Vertical velocity', 'm/s', 'mt')
        call ncinfo(ncname ( 4,:), 'thl','Liquid water potential temperature','K', 'tt')
        call ncinfo(ncname ( 5,:), 'qt', 'Specific humidity','kg/kg', 'tt')
        call ncinfo(ncname ( 6,:), 'ql', 'Liquid water specific humidity','kg/kg', 'tt')
        call ncinfo(ncname ( 7,:), 'sv001', 'Scalar 001 specific mixing ratio','kg/kg', 'tt')
        call ncinfo(ncname ( 8,:), 'sv002', 'Scalar 002 specific mixing ratio','kg/kg', 'tt')
        call ncinfo(ncname ( 9,:), 'rainrate', 'Rain rate','W/m2', 'tt')
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

    end if ! lcolumnstat

end subroutine initcolumnstat

end module modcolumnstat
