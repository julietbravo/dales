!> \file modcrosssection.f90
!!   Dumps an instantenous crosssection of the field

!>
!! Dumps an instantenous crosssection of the field.
!>
!! Crosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

!!  \par Revision list
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
module modcrosssection


  use modglobal, only : longint,kmax

implicit none
private
PUBLIC :: initcrosssection, crosssection,exitcrosssection
save
!NetCDF variables
  integer,parameter :: nvar = 12       !< number of variables in xy, xz and yz crossections
  integer,parameter :: nvar_path = 2   !< number of variables in the xy-path crosssections
  integer,parameter :: nvar_span = 6   !< number of variables in the spanwise averaged crosssections

  integer :: ncid1 = 0                 !< xz
  integer,allocatable :: ncid2(:)      !< xy
  integer :: ncid3 = 1                 !< yz
  integer :: ncid4 = 2                 !< xy-path
  integer :: ncid5 = 3                 !< spanwise xz

  integer :: nrec1 = 0
  integer,allocatable :: nrec2(:)
  integer :: nrec3 = 0
  integer :: nrec4 = 0
  integer :: nrec5 = 0

  integer :: crossheight(100)
  integer :: nxy = 0
  integer :: cross

  character(4) :: cheight
  character(80) :: fname1 = 'crossxz.xxxxyxxx.xxx.nc'
  character(80) :: fname2 = 'crossxy.xxxx.xxxxyxxx.xxx.nc'
  character(80) :: fname3 = 'crossyz.xxxxyxxx.xxx.nc'
  character(80) :: fname4 = 'crossxy.xxxxyxxx.xxx.nc'
  character(80) :: fname5 = 'crossxzspan.xxxxyxxx.xxx.nc'

  character(80),dimension(nvar,4) :: ncname1
  character(80),dimension(1,4) :: tncname1
  character(80),dimension(nvar,4) :: ncname2
  character(80),dimension(1,4) :: tncname2
  character(80),dimension(nvar,4) :: ncname3
  character(80),dimension(1,4) :: tncname3
  character(80),dimension(nvar_path,4) :: ncname4
  character(80),dimension(1,4) :: tncname4
  character(80),dimension(nvar_span,4) :: ncname5
  character(80),dimension(1,4) :: tncname5

  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcross = .false. !< switch for doing the crosssection (on/off)
  logical :: lbinary = .false. !< switch for doing the crosssection (on/off)
  integer :: crossplane = 2 !< Location of the xz crosssection
  integer :: crossortho = 2 !< Location of the yz crosssection

  logical :: lxy = .true.   !< switch for doing xy crosssections
  logical :: lxz = .true.   !< switch for doing xz crosssections
  logical :: lyz = .true.   !< switch for doing yz crosssections
  logical :: lpath = .true. !< switch for doing xy path crosssections
  logical :: lspan = .true. !< switch for doing spanwise avearged xz crosssections

contains
!> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid,myidx,myidy
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,kmax,i1,dt_lim,cexpnr,tres,btime
    use modstat_nc,only : lnetcdf,open_nc, define_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr,k

    namelist/NAMCROSSSECTION/ &
    lcross, lbinary, dtav, crossheight, crossplane, crossortho, lxy, lxz, lyz, lpath, lspan

    allocate(ncid2(kmax),nrec2(kmax))
    crossheight(1)=2
    crossheight(2:100)=-999
    ncid2(1)=4
    ncid2(2:kmax)=0
    nrec2(1:kmax)=0

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCROSSSECTION,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCROSSSECTION'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCROSSSECTION'
      endif
      write(6 ,NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcross     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lxy        ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lxz        ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lyz        ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lpath      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lspan      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call mpi_bcast(lbinary    ,1,mpi_logical,0,comm3d,mpierr)
    call MPI_BCAST(crossheight(1:100),100,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossplane ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossortho ,1,MPI_INTEGER,0,comm3d,mpierr)

    nxy=0
    k=1
    do while (crossheight(k) > 0)
    nxy=nxy+1
    ncid2(k)=k+1
    nrec2(k)=0
    k=k+1
    end do

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcross)) return
    dt_lim = min(dt_lim,tnext)

    if(any((crossheight(1:100).gt.kmax)) .or. crossplane>j1 .or. crossortho> i1 ) then
      stop 'CROSSSECTION: crosssection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'CROSSSECTION: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then

      if (myidy==0 .and. lxz) then
        fname1(9:16) = cmyid
        fname1(18:20) = cexpnr
        call ncinfo(tncname1(1,:),'time','Time','s','time')
        call ncinfo(ncname1( 1,:),'uxz', 'xz crosssection of the West-East velocity','m/s','m0tt')
        call ncinfo(ncname1( 2,:),'vxz', 'xz crosssection of the South-North velocity','m/s','t0tt')
        call ncinfo(ncname1( 3,:),'wxz', 'xz crosssection of the Vertical velocity','m/s','t0mt')
        call ncinfo(ncname1( 4,:),'thlxz','xz crosssection of the Liquid water potential temperature','K','t0tt')
        call ncinfo(ncname1( 5,:),'thvxz','xz crosssection of the Virtual potential temperature','K','t0tt')
        call ncinfo(ncname1( 6,:),'qtxz','xz crosssection of the Total water specific humidity','kg/kg','t0tt')
        call ncinfo(ncname1( 7,:),'qlxz','xz crosssection of the Liquid water specific humidity','kg/kg','t0tt')
        call ncinfo(ncname1( 8,:),'buoyxz','xz crosssection of the Buoyancy','K','t0tt')
        call ncinfo(ncname1( 9,:),'qrxz','xz crosssection of the Rain water specific humidity','kg/kg','t0tt')
        call ncinfo(ncname1( 10,:),'nrxz','xz crosssection of the Number concentration','-','t0tt')
        call ncinfo(ncname1( 11,:),'cloudnrxz','xz crosssection of the cloud number','-','t0tt')
        call ncinfo(ncname1( 12,:),'e120xz','xz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','t0tt')
        call open_nc(fname1,  ncid1,nrec1,n1=imax,n3=kmax)
        if (nrec1 == 0) then
          call define_nc( ncid1, 1, tncname1)
          call writestat_dims_nc(ncid1)
        end if
        call define_nc( ncid1, NVar, ncname1)
      end if

      if (lxy) then
        do cross=1,nxy
          write(cheight,'(i4.4)') crossheight(cross)
          fname2(9:12) = cheight
          fname2(14:21) = cmyid
          fname2(23:25) = cexpnr
          call ncinfo(tncname2(1,:),'time','Time','s','time')
          call ncinfo(ncname2( 1,:),'uxy','xy crosssections of the West-East velocity','m/s','mt0t')
          call ncinfo(ncname2( 2,:),'vxy','xy crosssections of the South-North velocity','m/s','tm0t')
          call ncinfo(ncname2( 3,:),'wxy','xy crosssections of the Vertical velocity','m/s','tt0t')
          call ncinfo(ncname2( 4,:),'thlxy','xy crosssections of the Liquid water potential temperature','K','tt0t')
          call ncinfo(ncname2( 5,:),'thvxy','xy crosssections of the Virtual potential temperature','K','tt0t')
          call ncinfo(ncname2( 6,:),'qtxy','xy crosssections of the Total water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2( 7,:),'qlxy','xy crosssections of the Liquid water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2( 8,:),'buoyxy','xy crosssection of the Buoyancy','K','tt0t')
          call ncinfo(ncname2( 9,:),'qrxy','xy crosssection of the Rain water specific humidity','kg/kg','tt0t')
          call ncinfo(ncname2(10,:),'nrxy','xy crosssection of the rain droplet number concentration','-','tt0t')
          call ncinfo(ncname2(11,:),'cloudnrxy','xy crosssection of the cloud number','-','tt0t')
          call ncinfo(ncname2(12,:),'e120xy','xy crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','tt0t')
          call open_nc(fname2,  ncid2(cross),nrec2(cross),n1=imax,n2=jmax)
          if (nrec2(cross)==0) then
            call define_nc( ncid2(cross), 1, tncname2)
            call writestat_dims_nc(ncid2(cross))
          end if
          call define_nc( ncid2(cross), NVar, ncname2)
        end do
      end if

      if (myidx==0 .and. lyz) then
        fname3(9:16) = cmyid
        fname3(18:20) = cexpnr
        call ncinfo(tncname3(1,:),'time','Time','s','time')
        call ncinfo(ncname3( 1,:),'uyz','yz crosssection of the West-East velocity','m/s','0ttt')
        call ncinfo(ncname3( 2,:),'vyz','yz crosssection of the South-North velocity','m/s','0mtt')
        call ncinfo(ncname3( 3,:),'wyz','yz crosssection of the Vertical velocity','m/s','0tmt')
        call ncinfo(ncname3( 4,:),'thlyz','yz crosssection of the Liquid water potential temperature','K','0ttt')
        call ncinfo(ncname3( 5,:),'thvyz','yz crosssection of the Virtual potential temperature','K','0ttt')
        call ncinfo(ncname3( 6,:),'qtyz','yz crosssection of the Total water specific humidity','kg/kg','0ttt')
        call ncinfo(ncname3( 7,:),'qlyz','yz crosssection of the Liquid water specific humidity','kg/kg','0ttt')
        call ncinfo(ncname3( 8,:),'buoyyz','yz crosssection of the Buoyancy','K','0ttt')
        call ncinfo(ncname3( 9,:),'qryz','yz crosssection of the Rain water specific humidity','kg/kg','0ttt')
        call ncinfo(ncname3(10,:),'nryz','yz crosssection of the Number concentration','-','0ttt')
        call ncinfo(ncname3(11,:),'cloudnryz','yz crosssection of the cloud number','-','0ttt')
        call ncinfo(ncname3(12,:),'e120yz','yz crosssection of sqrt(turbulent kinetic energy)','m^2/s^2','0ttt')
        call open_nc(fname3, ncid3, nrec3, n2=jmax, n3=kmax)
        if (nrec3==0) then
          call define_nc(ncid3, 1, tncname3)
          call writestat_dims_nc(ncid3)
        end if
        call define_nc( ncid3, NVar, ncname3)
      end if

      if (lpath) then
        fname4(9:16) = cmyid
        fname4(18:20) = cexpnr
        call ncinfo(tncname4(1,:),'time','Time','s','time')
        call ncinfo(ncname4( 1,:),'lwpxy','Cloud liquid water path','kg/m2','tt0t')
        call ncinfo(ncname4( 2,:),'rwpxy','Rain water path','kg/m2','tt0t')
        call open_nc(fname4, ncid4, nrec4, n1=imax, n2=jmax)
        if (nrec4==0) then
          call define_nc( ncid4, 1, tncname4)
          call writestat_dims_nc(ncid4)
        end if
        call define_nc( ncid4, nvar_path, ncname4)
      end if

      if (myidy==0 .and. lspan) then
        fname5(13:20) = cmyid
        fname5(22:24) = cexpnr
        call ncinfo(tncname5(1,:), 'time', 'Time', 's', 'time')
        call ncinfo(ncname5 (1,:), 'u2rxz',   'spanwise avaraged xz crosssection of the resolved u variance',   'm2/s2',   'm0tt')
        call ncinfo(ncname5 (2,:), 'v2rxz',   'spanwise avaraged xz crosssection of the resolved v variance',   'm2/s2',   't0tt')
        call ncinfo(ncname5 (3,:), 'w2rxz',   'spanwise avaraged xz crosssection of the resolved w variance',   'm2/s2',   't0mt')
        call ncinfo(ncname5 (4,:), 'thl2rxz', 'spanwise avaraged xz crosssection of the resolved thl variance', 'K2',      't0tt')
        call ncinfo(ncname5 (5,:), 'qt2rxz',  'spanwise avaraged xz crosssection of the resolved qt variance',  'kg2/kg2', 't0tt')
        call ncinfo(ncname5 (6,:), 'ql2rxz',  'spanwise avaraged xz crosssection of the resolved ql variance',  'kg2/kg2', 't0tt')
        call open_nc(fname5, ncid5, nrec5, n1=imax, n3=kmax)
        if (nrec5 == 0) then
          call define_nc( ncid5, 1, tncname5)
          call writestat_dims_nc(ncid5)
        end if
        call define_nc(ncid5, nvar_span, ncname5)
      end if

    end if


  end subroutine initcrosssection
!>Run crosssection. Mainly timekeeping
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lxz)   call wrtvert
    if (lxy)   call wrthorz
    if (lyz)   call wrtorth
    if (lpath) call wrtpath
    if (lspan) call wrtspan

  end subroutine crosssection


!> Do the xz crosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k,n
  character(20) :: name

  real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

  if( myidy /= 0 ) return 

  allocate(thv0(2:i1,1:kmax),buoy(2:i1,1:kmax))


    do  i=2,i1
    do  k=1,kmax
      thv0(i,k) = (thl0(i,crossplane,k)+rlv*ql0(i,crossplane,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(i,crossplane,k)-rv/rd*ql0(i,crossplane,k))
      buoy(i,k) = thv0(i,k)-thvf(k)
    enddo
    enddo

    if(lbinary) then
      open(ifoutput,file='movv_u.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(i,crossplane,k)+cu,i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_v.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(i,crossplane,k)+cv,i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_w.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_thl.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_thv.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(i,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_buoy.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(i,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_qt.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movv_ql.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(i,crossplane,k),i=2,i1),k=1,kmax)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_tnn.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(i,crossplane,k,n),i=2,i1),k=1,kmax)
        close(ifoutput)
      end do
    end if

    if (lnetcdf) then
      allocate(vars(1:imax,1:kmax,nvar))
      vars(:,:,1) = um(2:i1,crossplane,1:kmax)+cu
      vars(:,:,2) = vm(2:i1,crossplane,1:kmax)+cv
      vars(:,:,3) = wm(2:i1,crossplane,1:kmax)
      vars(:,:,4) = thlm(2:i1,crossplane,1:kmax)
      vars(:,:,5) = thv0(2:i1,1:kmax)
      vars(:,:,6) = qtm(2:i1,crossplane,1:kmax)
      vars(:,:,7) = ql0(2:i1,crossplane,1:kmax)
      vars(:,:,8) = buoy(2:i1,1:kmax)
      if(nsv>1) then
      vars(:,:,9) = svm(2:i1,crossplane,1:kmax,2)
      vars(:,:,10) = svm(2:i1,crossplane,1:kmax,1)
      else
      vars(:,:,9) = 0.
      vars(:,:,10) = 0.
      end if
      vars(:,:,11) = cloudnr(2:i1,crossplane,1:kmax)
      vars(:,:,12) = e120(2:i1,crossplane,1:kmax)
      call writestat_nc(ncid1,1,tncname1,(/rtimee/),nrec1,.true.)
      call writestat_nc(ncid1,nvar,ncname1(1:nvar,:),vars,nrec1,imax,kmax)
      deallocate(vars)
    end if
    deallocate(thv0,buoy)

  end subroutine wrtvert

!> Do the xy crosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr,inr
    implicit none


    ! LOCAL
    integer i,j,n
    character(40) :: name
    real, allocatable :: thv0(:,:,:),vars(:,:,:),buoy(:,:,:)

    allocate(thv0(2:i1,2:j1,nxy),buoy(2:i1,2:j1,nxy))

    do  cross=1,nxy
    do  j=2,j1
    do  i=2,i1
      thv0(i,j,cross) =&
       (thl0(i,j,crossheight(cross))+&
       rlv*ql0(i,j,crossheight(cross))/&
       (cp*exnf(crossheight(cross)))) &
                    *(1+(rv/rd-1)*qt0(i,j,crossheight(cross))&
                    -rv/rd*ql0(i,j,crossheight(cross)))
      buoy(i,j,cross) =thv0(i,j,cross)-thvf(crossheight(cross))
    enddo
    enddo
    enddo

    if(lbinary) then
      do  cross=1,nxy
      write(cheight,'(i4.4)') crossheight(cross)
      open(ifoutput,file='movh_u.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(i,j,crossheight(cross))+cu,i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_v.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(i,j,crossheight(cross))+cv,i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_w.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_thl.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_thv.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(i,j,cross),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_buoy.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(i,j,cross),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_qt.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      open(ifoutput,file='movh_ql.'//cheight//'.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(i,j,crossheight(cross)),i=2,i1),j=2,j1)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_snn.'//trim(cheight)//'.'//cmyid//'.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(i,j,crossheight(cross),n),i=2,i1),j=2,j1)
        close(ifoutput)
      end do
      end do
    endif

    if (lnetcdf) then
      do cross=1,nxy
      allocate(vars(1:imax,1:jmax,nvar))
      vars=0.
      vars(:,:,1) = um(2:i1,2:j1,crossheight(cross))+cu
      vars(:,:,2) = vm(2:i1,2:j1,crossheight(cross))+cv
      vars(:,:,3) = wm(2:i1,2:j1,crossheight(cross))
      vars(:,:,4) = thlm(2:i1,2:j1,crossheight(cross))
      vars(:,:,5) = thv0(2:i1,2:j1,cross)
      vars(:,:,6) = qtm(2:i1,2:j1,crossheight(cross))
      vars(:,:,7) = ql0(2:i1,2:j1,crossheight(cross))
      vars(:,:,8) = buoy(2:i1,2:j1,cross)
      if(nsv>1) then
      vars(:,:,9) = svm(2:i1,2:j1,crossheight(cross),iqr)
      vars(:,:,10) = svm(2:i1,2:j1,crossheight(cross),inr)
      else
      vars(:,:,9) = 0.
      vars(:,:,10) = 0.
      end if
      vars(:,:,11) = cloudnr(2:i1,2:j1,crossheight(cross))
      vars(:,:,12) = e120(2:i1,2:j1,crossheight(cross))
      call writestat_nc(ncid2(cross),1,tncname2,(/rtimee/),nrec2(cross),.true.)
      call writestat_nc(ncid2(cross),nvar,ncname2(1:nvar,:),vars,nrec2(cross),imax,jmax)
      deallocate(vars)
      end do
    end if

    deallocate(thv0,buoy)

  end subroutine wrthorz

!> Do the xy-path crosssections and dump them to file
  subroutine wrtpath
    use modglobal, only : imax,jmax,i1,j1,nsv,rtimee,dzf
    use modfields, only : svm,ql0,rhof
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmicrodata, only : iqr
    implicit none

    integer i,j,k
    real, allocatable :: lwp(:,:), rwp(:,:), vars(:,:,:)

    allocate(lwp(2:i1,2:j1),rwp(2:i1,2:j1))

    ! Calculate cloud and rain water paths
    do j=2,j1
      do i=2,i1
        lwp(i,j) = 0
        rwp(i,j) = 0
        do k=1,kmax
          lwp(i,j) = lwp(i,j) + rhof(k) * ql0(i,j,k) * dzf(k)

          if(nsv>1) then
            rwp(i,j) = rwp(i,j) + rhof(k) * svm(i,j,k,iqr) * dzf(k)
          end if
        end do
      end do
    end do

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,nvar_path))
      vars=0.
      vars(:,:,1) = lwp(2:i1,2:j1)
      vars(:,:,2) = rwp(2:i1,2:j1)
      call writestat_nc(ncid4,1,tncname4,(/rtimee/),nrec4,.true.)
      call writestat_nc(ncid4,nvar_path,ncname4,vars,nrec4,imax,jmax)
      deallocate(vars)
    end if

    deallocate(lwp,rwp)

  end subroutine wrtpath

  ! yz cross section
  subroutine wrtorth
    use modglobal, only : jmax,kmax,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput,rtimee
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,e120,exnf,thvf,cloudnr
    use modmpi,    only : cmyid, myidx
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    integer j,k,n
    character(21) :: name

    real, allocatable :: thv0(:,:),vars(:,:,:),buoy(:,:)

    if( myidx /= 0 ) return 


    allocate(thv0(1:j1,1:kmax),buoy(1:j1,1:kmax))

    do  j=1,j1
    do  k=1,kmax
      thv0(j,k) =&
       (thl0(crossortho,j,k)+&
       rlv*ql0(crossortho,j,k)/&
       (cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(crossortho,j,k)&
                    -rv/rd*ql0(crossortho,j,k))
       buoy(j,k) =thv0(j,k)-thvf(k)
    enddo
    enddo

    if(lbinary) then
      open(ifoutput,file='movo_u.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((um(crossortho,j,k)+cu,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_v.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((vm(crossortho,j,k)+cv,j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_w.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((wm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thlm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((thv0(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_buoy.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((buoy(j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*qtm(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      open(ifoutput,file='movo_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
      write(ifoutput,'(es12.5)') ((1.e3*ql0(crossortho,j,k),j=2,j1),k=1,kmax)
      close(ifoutput)

      do n = 1,nsv
        name = 'movh_tnn.'//cmyid//'.'//cexpnr
        write(name(7:8),'(i2.2)') n
        open(ifoutput,file=name,position='append',action='write')
        write(ifoutput,'(es12.5)') ((svm(crossortho,j,k,n),j=2,j1),k=1,kmax)
        close(ifoutput)
      end do
    end if

    if (lnetcdf) then
      allocate(vars(1:jmax,1:kmax,nvar))
      vars(:,:,1) = um(crossortho,2:j1,1:kmax)+cu
      vars(:,:,2) = vm(crossortho,2:j1,1:kmax)+cv
      vars(:,:,3) = wm(crossortho,2:j1,1:kmax)
      vars(:,:,4) = thlm(crossortho,2:j1,1:kmax)
      vars(:,:,5) = thv0(2:j1,1:kmax)
      vars(:,:,6) = qtm(crossortho,2:j1,1:kmax)
      vars(:,:,7) = ql0(crossortho,2:j1,1:kmax)
      vars(:,:,8) = buoy(2:j1,1:kmax)
      if(nsv>1) then
      vars(:,:,9) = svm(crossortho,2:j1,1:kmax,2)
      vars(:,:,10) = svm(crossortho,2:j1,1:kmax,1)
      else
      vars(:,:,9) = 0.
      vars(:,:,10) = 0.
      end if
      vars(:,:,11) = cloudnr(crossortho,2:j1,1:kmax)
      vars(:,:,12) = e120(crossortho,2:j1,1:kmax)
      call writestat_nc(ncid3,1,tncname3,(/rtimee/),nrec3,.true.)
      call writestat_nc(ncid3,nvar,ncname3(1:nvar,:),vars,nrec3,jmax,kmax)
      deallocate(vars)
    end if

    deallocate(thv0,buoy)

  end subroutine wrtorth

!> Do the spanwise averaged xz crosssections and dump them to file
  subroutine wrtspan
    use modglobal, only : imax, kmax, rtimee
    use modfields, only : um, vm, wm, thlm, qtm, ql0
    use modmpi,    only : myidy
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    real, allocatable :: var_u(:,:), var_v(:,:), var_w(:,:), var_thl(:,:), var_qt(:,:), var_ql(:,:)
    real, allocatable :: vars(:,:,:)

    allocate(var_u  (imax, kmax), var_v (imax, kmax), var_w (imax, kmax), &
           & var_thl(imax, kmax), var_qt(imax, kmax), var_ql(imax, kmax))

    call calc_spanwise_variance(um,   var_u)
    call calc_spanwise_variance(vm,   var_v)
    call calc_spanwise_variance(wm,   var_w)
    call calc_spanwise_variance(thlm, var_thl)
    call calc_spanwise_variance(qtm,  var_qt)
    call calc_spanwise_variance(ql0,  var_ql)

    if (lnetcdf .and. myidy==0) then
      allocate(vars(1:imax, 1:kmax, nvar_span))
      vars(:,:,1) = var_u
      vars(:,:,2) = var_v
      vars(:,:,3) = var_w
      vars(:,:,4) = var_thl
      vars(:,:,5) = var_qt
      vars(:,:,6) = var_ql

      call writestat_nc(ncid5, 1, tncname5, (/rtimee/), nrec5, .true.)
      call writestat_nc(ncid5, nvar_span, ncname5(1:nvar_span,:), vars, nrec5, imax, kmax)
      deallocate(vars)
    end if

    deallocate(var_u, var_v, var_w, var_thl, var_qt, var_ql)

  end subroutine wrtspan

  ! Help routine to calculate spanwise averaged variance
  subroutine calc_spanwise_variance(field, variance)
    use modglobal, only : imax, jmax, kmax, i1, j1, k1, ih, jh
    use modmpi,    only : sumcol, nprocy
    implicit none

    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: field
    real, dimension(imax, kmax),               intent(out) :: variance
    integer :: i,j,k

    real, allocatable :: mean_local(:,:), mean_global(:,:)
    allocate(mean_local(imax,kmax), mean_global(imax,kmax))

    mean_local = 0

    ! Calculate spanwise averaged quantities
    ! Strip ghost cells for MPI routines...
    do k=1,kmax
      do j=2,j1
        do i=1,imax
          mean_local(i,k) = mean_local(i,k) + field(i+1,j,k)
        end do
      end do
    end do

    call sumcol(mean_local, mean_global, imax, kmax)
    mean_global = mean_global / (jmax * nprocy)

    ! Calculate variance
    mean_local=0

    do k=1,kmax
      do j=2,j1
        do i=1,imax
          mean_local(i,k) = mean_local(i,k) + (field(i+1,j,k) - mean_global(i,k))**2.
        end do
      end do
    end do

    call sumcol(mean_local, variance, imax, kmax)
    variance = variance / (jmax * nprocy)

    deallocate(mean_local, mean_global)

  end subroutine calc_spanwise_variance

!> Clean up when leaving the run
  subroutine exitcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myidx, myidy
    implicit none

    if(lcross .and. lnetcdf) then
      if (myidy ==0 .and. lxz)   call exitstat_nc(ncid1)

      if (lxy) then
        do cross=1,nxy
          call exitstat_nc(ncid2(cross))
        end do
      end if

      if (myidx==0 .and. lyz)   call exitstat_nc(ncid3)
      if (lpath)                call exitstat_nc(ncid4)
      if (myidy==0 .and. lspan) call exitstat_nc(ncid5)
    end if

  end subroutine exitcrosssection

end module modcrosssection
