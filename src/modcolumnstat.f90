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
    integer, parameter :: nvar = 5
    integer :: nrec = 0, ncid = 666
    character(80) :: fname = 'column.x___y___.___.nc'
    character(80), dimension(nvar,4) :: ncname
    character(80), dimension(1,4)    :: tncname

    ! Location columns (in x,y and index space)
    real, allocatable    :: column_x(:), column_y(:)
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

    integer :: ierr, n, i, ncolumns
    integer :: maxcolumns = 5
    real    :: x0_mpi, y0_mpi, x1_mpi, y1_mpi, xsize_mpi, ysize_mpi, itot_mpi, jtot_mpi

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
        xsize_mpi = xsize / nprocx
        ysize_mpi = ysize / nprocy

        itot_mpi = itot / nprocx
        jtot_mpi = jtot / nprocy

        x0_mpi = myidx     * xsize_mpi
        x1_mpi = (myidx+1) * xsize_mpi
        y0_mpi = myidy     * ysize_mpi
        y1_mpi = (myidy+1) * ysize_mpi

        ncolumns_loc = 0
        do n=1, ncolumns
            if (column_x(n) >= x0_mpi .and. column_x(n) < x1_mpi .and. column_y(n) > y0_mpi .and. column_y(n) < y1_mpi) then
                ncolumns_loc = ncolumns_loc + 1
            end if
        end do

        ! Find indices (nearest grid center) of columns
        allocate(column_i(ncolumns_loc), column_j(ncolumns_loc))

        i = 1
        do n=1, ncolumns
            if (column_x(n) >= x0_mpi .and. column_x(n) < x1_mpi .and. column_y(n) > y0_mpi .and. column_y(n) < y1_mpi) then
                column_i(i) = int(column_x(n)/dx)+1 - myidx*itot_mpi
                column_j(i) = int(column_y(n)/dy)+1 - myidy*jtot_mpi
                i = i+1
            end if
        end do

        ! Create NetCDF file (one per MPI tasks for now...)
        fname(17:19) = cexpnr
        fname(9:11)  = cmyidx
        fname(13:15) = cmyidy

        ! Define variables and their names / units / ..
        call ncinfo(tncname(1,:), 'time', 'Time', 's', 'time')
        call ncinfo(ncname (1,:), 'u', 'West-east velocity', 'm/s', 'tt')
        call ncinfo(ncname (2,:), 'v', 'South-North velocity', 'm/s', 'tt')
        call ncinfo(ncname (3,:), 'w', 'Vertical velocity', 'm/s', 'mt')
        call ncinfo(ncname (4,:), 'thl','Liquid water potential temperature','K', 'tt')
        call ncinfo(ncname (5,:), 'qt', 'Specific humidity','kg/kg', 'tt')

        call open_nc(fname, ncid, nrec, n3=kmax)
        if (nrec == 0) then
            call define_nc(ncid, 1, tncname)
            call writestat_dims_nc(ncid)
        end if
        call define_nc(ncid, nvar, ncname)

    end if ! lcolumnstat

end subroutine initcolumnstat

end module modcolumnstat
