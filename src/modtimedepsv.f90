!> \file modtimedepsv.f90
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars

!>
!!  Prescribes surface values, fluxes and LS forcings at certain times for scalars
!>
!!  \author Roel Neggers, KNMI
!!  \author Thijs Heus,MPI-M
!!  \author Stephan de Roode, TU Delft
!!  \author Simon Axelsen, UU
!!  \par Revision list
!! \todo documentation
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



module modtimedepsv

use netcdf
implicit none
private
public :: inittimedepsv, timedepsv,ltimedepsv,exittimedepsv
save
! switches for timedependent surface fluxes and large scale forcings
  logical       :: ltimedepsv     = .false. !< Overall switch, input in namoptions
  logical       :: ltimedepsvz    = .false. !< Switch for large scale forcings
  logical       :: ltimedepsvsurf = .true.  !< Switch for surface fluxes

  integer, parameter    :: kflux = 100
  integer, parameter    :: kls   = 100
  real, allocatable     :: timesvsurf (:)
  real, allocatable     :: svst     (:,:) !< Time dependent surface scalar concentration

  real, allocatable     :: timesvz  (:)
  real, allocatable     :: svzt(:,:,:) !< Time dependent, height dependent scalar concentrations



contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedepsv
    use modmpi,   only :myid,my_real,mpi_logical,mpierr,comm3d
    use modglobal,only :cexpnr,kmax,k1,ifinput,runtime,nsv,&
                         nc_input  !cstep

    implicit none

    character (80):: chmess
    character (1) :: chmess1
    character (16) :: outputfmt !format used to write the input read to stdout
    integer :: k,t,n, ierr
    real :: dummyr
    real, allocatable, dimension (:) :: height

    !cstep : NCDF vars
    integer :: ncid,varID,sts ! checks for presence ncdf input file
    integer :: nsv_in         ! number of passive scalars present in ncdf input file. 
                              !if nsv_in < nsv, profiles are assumed to be zero
    integer :: dimIDkls,dimIDkflux,dimIDnsv
    integer :: kflux_nc, kls_nc
    character(len=nf90_max_name) :: tmpName


    if (nsv==0 .or. .not.ltimedepsv ) return  !cstep : only do if ltimedepsv = .true.

    allocate(height(k1))
    allocate(timesvsurf (0:kflux))
    allocate(svst  (kflux,nsv))
    allocate(timesvz  (0:kls))

    allocate(svzt(k1,kls,nsv))
    timesvsurf = 0
    timesvz   = 0
    svst       = 0
    svzt       = 0

    if (myid==0) then

       outputfmt = '(f10.3,100e10.3)'
       write(outputfmt(8:10),'(I3)') nsv

!    --- load lsforcings---

      if (nc_input) then
         sts        = nf90_open('case_setup.nc',nf90_nowrite,ncid)

         sts          = nf90_inq_dimid(ncid,"Ntimesvz",dimIDkls)
         if (sts.ne.nf90_noerr) stop  'ERROR:  time dependent scalars not present'
         sts          = nf90_inquire_dimension(ncid, dimIDkls, tmpName,kls_nc)

         sts          = nf90_inq_dimid(ncid,"Ntimesvsurf",dimIDkflux)
         if (sts.ne.nf90_noerr) stop  'ERROR:  time dependent scalar fluxes not present'
         sts          = nf90_inquire_dimension(ncid, dimIDkflux, tmpName,kflux_nc)

         sts          = nf90_inq_dimid(ncid,"nsv",dimIDnsv)
         if (sts.ne.nf90_noerr) stop 'no scalars provided in ncdf file'

         sts          = nf90_inquire_dimension(ncid, dimIDnsv, tmpName,nsv_in)
         if (nsv_in.ne.nsv) then 
            stop 'nr of scalars in ncdf input file does not match nsv in namoptions'
         else
             sts          = nf90_inq_varid(ncid,"timesvz",varID)
             sts          = nf90_get_var(ncid, varID, timesvz(1:kls_nc))
 
             sts          = nf90_inq_varid(ncid,"timesvsurf",varID)
             sts          = nf90_get_var(ncid, varID, timesvsurf(1:kflux_nc))

         endif
    
         sts          = nf90_inq_varid(ncid,"zf",varID)
         sts          = nf90_get_var(ncid, varID, height(1:kmax))

         sts          = nf90_inq_varid(ncid,"svst",varID)
         sts          = nf90_get_var(ncid, varID, svst(1:kflux_nc,1:nsv))

         sts          = nf90_inq_varid(ncid,"svzt",varID)
         sts          = nf90_get_var(ncid, varID, svzt(1:kmax,1:kls_nc,1:nsv))


      else !cstep input from ascii file

         open(ifinput,file='ls_fluxsv.inp.'//cexpnr)
         read(ifinput,'(a80)') chmess
         write(6,*) chmess
         read(ifinput,'(a80)') chmess
         write(6,*) chmess
         read(ifinput,'(a80)') chmess
         write(6,*) chmess


!      --- load fluxes---
         t    = 0
         ierr = 0
         do while (timesvsurf(t)< runtime.and.t.lt.kflux)
            t=t+1
            read(ifinput,*, iostat = ierr) timesvsurf(t), (svst(t,n),n=1,nsv)
            if (ierr < 0) then
                stop 'STOP: No time dependend data for end of run (surface fluxes of scalar)'
            end if
         end do

         ! flush to the end of fluxlist
         do while (ierr ==0)
            read (ifinput,*,iostat=ierr) dummyr
         end do
!       ---load large scale forcings----
         t = 0

         do while (timesvz(t) < runtime.and.t.lt.kls)
            t = t + 1
            chmess1 = "#"
            ierr = 1 ! not zero
            !search for the next line consisting of "# time", from there onwards the profiles will be read
            do while (.not.(chmess1 == "#" .and. ierr ==0))
               read(ifinput,*,iostat=ierr) chmess1,timesvz(t)
               if (ierr < 0) then
                  stop 'STOP: No time dependend data (scalars) for end of run'
               end if
            end do
            write (*,*) 'timesvz = ',timesvz(t)
            do k=1,kmax
               read (ifinput,*) height(k), (svzt(k,t,n),n=1,nsv)
            end do
         end do

         close(ifinput)
      endif  !nc_input

      t    = 0
      do while (timesvsurf(t)< runtime)
         t=t+1
         write(*,'(f9.1,4e12.4)') timesvsurf(t), (svst(t,n),n=1,nsv)
      end do

      t = 0
      do while (timesvz(t) < runtime)
          t = t + 1
          write (*,*) 'timesvz = ',timesvz(t)
          do k=kmax,1,-1
             write (6,outputfmt) height(k),(svzt(k,t,n),n=1,nsv)
          end do
      enddo

      if(timesvsurf(1)>runtime) then
         write(6,*) 'Time dependent surface variables do not change before end of'
         write(6,*) 'simulation. --> only large scale changes in scalars'
         ltimedepsvsurf=.false.
      endif

      if ((timesvz(1) > runtime) .or. (timesvsurf(1) > runtime)) then
         write(6,*) 'Time dependent large scale forcings sets in after end of simulation -->'
         write(6,*) '--> only time dependent surface variables (scalars)'
         ltimedepsvz=.false.
      end if
   endif   !myid

   
    call MPI_BCAST(timesvsurf(1:kflux),kflux,MY_REAL,0,comm3d,mpierr)
    !cstep  call MPI_BCAST(timesvz(1:kflux),kflux,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(svst             ,kflux*nsv,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(timesvz(1:kls)    ,kls,MY_REAL  ,0,comm3d,mpierr)
    call MPI_BCAST(ltimedepsvsurf ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltimedepsvz    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    do n=1,nsv
         call MPI_BCAST(svzt(1:k1,1:kls,n),kmax*kls,MY_REAL,0,comm3d,mpierr)
    enddo
    call timedepsv

    deallocate(height)

  end subroutine inittimedepsv

  subroutine timedepsv
    use modglobal, only : nsv
    implicit none

    if(nsv==0 .or. .not.ltimedepsv) return
    call timedepsvz
    call timedepsvsurf

  end subroutine timedepsv

  subroutine timedepsvz
  implicit none

    if(.not.(ltimedepsvz)) return
    stop 'Modtimedepsv: time dependent scalars at all levels not programmed'

    return
  end subroutine timedepsvz

  subroutine timedepsvsurf
    use modglobal,   only : rtimee,nsv
    use modsurfdata,  only : svs
    implicit none
    integer t,n
    real fac

    if(.not.(ltimedepsvsurf)) return

  !     --- interpolate! ----
    t=1
    do while(rtimee>timesvsurf(t))
      t=t+1
    end do
    if (rtimee>timesvsurf(t)) then
      t=t-1
    end if

    fac = ( rtimee-timesvsurf(t) ) / ( timesvsurf(t+1)-timesvsurf(t))
    do n=1,nsv
       svs(n) = svst(t,n) + fac * (svst(t+1,n) - svst(t,n))
    enddo
    return
  end subroutine timedepsvsurf


  subroutine exittimedepsv
    use modglobal, only : nsv
    implicit none
    if (nsv==0 .or. .not.ltimedepsv) return
    deallocate(timesvz,svzt,timesvsurf)
  end subroutine exittimedepsv

end module modtimedepsv
