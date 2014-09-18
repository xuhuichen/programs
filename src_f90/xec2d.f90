      subroutine xec
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
!     MPI variables

      integer status(MPI_STATUS_SIZE)

      integer n_lc, nunit_temp, nunit_el, i
      integer wallm
!
      double precision drmin
      real cputime1,cputime2
!
      character *43 localeventfile, evlfilename
!
      real etime, elapse(2), et0
!
!      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
!
!     wall time of the program in minutes. 
!     When time is approaching this wall, the program will write out most variables for a restart,
!     and stop the current run.
      wallm=60*24
!
!
!
!     master node part
      if(myid.eq.master) then
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       write(*,*) 'compton2d_volpara'
      etotal = 0
      do 50 i = 1, num_nt
 50   nelectron(i) = 0
      nunit_evt = myid+1000
      nunit_temp = 15
      nunit_el = 16
      open(unit=nunit_evt, file=eventfile, access='stream')
      open(unit=nunit_temp, file=temp_file, status='unknown')
!
 100  continue

      etotal_old = etotal
      etotal=etime(elapse)
      write(*,*) 'elapsed time:',etotal
      write(*,*) 'myid=', myid, ' calling xec_bcast,'
      call xec_bcast
      write(*,*) 'myid=', myid, ' called xec_bcast,'

      if(etotal.gt.(0.95*wallm-3)*60)then
        write(*,*) 'keeping a record of the current run...'
        call write_record
        write(*,*) 'myid=',myid,'finished recording. exit.'
        goto 900
      endif

      call xec_add
      write(*,*) 'myid=', myid, ' called xec_add'

!
      write(4,110) ncycle,time
      write(*,110) ncycle,time
 110  format(i4,'. time step: t = ',e14.7,' s')
!

      write(*,*) 'Calling imcgen ... '
      call imcgen2d
      et0=etime(elapse)
      write(*,*) 'myid=', myid, 'finished imcgen at(s):',et0,&
     & ' Calling imcfield  ... ncycle=', ncycle
      call imcfield2d 
      et0=etime(elapse)
      write(*,*)'myid=',myid, 'finished imcfield at(s):',et0,&
     &  'Calling imcvol ...'
      call imcvol2d
      et0=etime(elapse)
      write(*,*) 'myid=', myid,'finished imcvol at(s):',et0,&
     & ' Calling imcsurf ...'
      call imcsurf2d
      et0=etime(elapse)
      write(*,*) 'myid=', myid,'finished imcsurf at(s):',et0,&
     & ' Calling imcredist ...'
      call imcredist
      write(*,*) 'myid=',myid, 'finished imcredist at(s):',et0,&
     & ' Calling update ...'
      if ((fp_sw.ne.0).and.(ncycle.gt.0))&
     & call update
      et0=etime(elapse)
      write(*,*) 'myid=', myid, ' called update at:',et0
      write(*,*) 'Calling graphics ...'
      call graphics_collect
      write(*,*) 'myid=',myid,' af graphics_collect,'
      if (ncycle.gt.0) call graphics
      write(*,*) 'myid=',myid,' af graphics'
!

!     advance in time step
      dt(2) = dt(1)
      if (ncycle.gt.0) then
         time = time + dt(1)
      else if ((fp_sw.eq.0).and.(pair_switch.eq.1)) then
         drmin = dmin1((r(nr) - rmin), (z(nz) - zmin))
         dt(1) = 3.d-12*drmin
      endif
      ncycle = ncycle + 1
!
       if ((time-dt(2)).lt.tstop) then
        if((time-dt(2)) < currenttimestop) then
         goto 100
        else
        write(*,*) "next mhdtimestep, currenttimestop=",currenttimestop,&
     &             "time=",time-dt(2)
         call reader(1)
         call setup(1)
         goto 100
        endif
       endif 

!
      close(nunit_evt)
      close(nunit_temp)
!
      open(unit=nunit_el, file='esp.dat', status='unknown')
      do 200 i = 1, num_nt
         if (nelectron(i).gt.0) then
            write(nunit_el, 210) gnt(i), nelectron(i)
         else
            write(nunit_el, 215) gnt(i), 1.d-2
         endif
 200  continue
      close(nunit_el)
!
 210  format (e14.7,1x,i7)
 215  format (e14.7,1x,e14.7)
!
!
!
!     slave node part
      else
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      localeventfile = evlfilename(eventfile, myid)

      nunit_evt = myid+1000
      open(unit=nunit_evt, file=localeventfile, access='stream')
!
 101  continue
      lkcount = 0
      cmcount1 = 0
      cmcount2 = 0 
!      write(*,*) 'myid=', myid, ' calling xec_bcast'
      call xec_bcast
!      write(*,*) 'myid=', myid, ' called xec_bcast'

      if(etotal.gt.(0.95*wallm-3)*60.and.etotal_old.le.wallm*60)then
        write(*,*) 'keeping a record of the current run...'
        call write_record
        write(*,*) 'myid=',myid,'finished recording. exit'
        goto 900
      endif

      call xec_add
!      write(*,*) 'myid=', myid, ' called xec_add'

!
!
!
!      write(*,*) 'myid=', myid, 'Calling imcgen ... '
!      write(*,*) 'edep=', edep(1,1)
      call imcgen2d
!
      call imcfield2d
!      write(*,*)'myid=',myid, 'finished imcfield,',
!     1  'Calling imcvol ...'
      call imcvol2d
!      write(*,*) 'myid=', myid,'finished imcvol,',
!     1 ' Calling imcsurf ...'
      call imcsurf2d
!      write(*,*) 'myid=', myid,'finished imcsurf,',
!     1 ' Calling imcredist ...'
      call imcredist
!      write(*,*) 'myid=',myid, 'finished imcredist:',
!     1 ' Calling update ...'
!
      if ((fp_sw.ne.0).and.(ncycle.gt.0))&
     & call update
!      write(*,*) 'myid=', myid, ' called update'
      call graphics_collect
!      write(*,*) 'myid=',myid,' af graphics_collect'
!
!
!
      if(cmcount1.gt.0)write(*,*) 'myid=',myid,' time=',time,&
     &     ' cmcount1=',cmcount1,'cmcount2=',cmcount2 ! Xuhui
      if (time.lt.tstop) then 
       if(time < currenttimestop) then
        goto 101
       else
!        write(*,*) "myid=",myid,"currenttimestop=",currenttimestop,
!     1             "time=",time
        call setup_bcast
        goto 101
       endif
      endif
      close(nunit_evt)
!
!
!
      endif
 900  return
      end
!
!
!
!
!
      character*43 function evlfilename(eventfile, num)
      implicit none
      character*30 eventfile
      integer num
!
      character*43 localeventfile
!
      evlfilename(14:43) = eventfile(1:30)
         evlfilename(9:9) = 'p'
      if(num.lt.10) then
         evlfilename(10:10) = '0'
         evlfilename(11:11) = '0'
         evlfilename(12:12) = char(num+48)

      else if(num.lt.100) then
         evlfilename(10:10) = '0'
         evlfilename(11:11) = char(num/10+48)
         evlfilename(12:12) = char(num-num/10*10+48)
      else
         evlfilename(10:10) = char(num/100+48)
         evlfilename(11:11) = char((num-num/100*100)/10+48)
         evlfilename(12:12) = char(num-num/100*100&
     &        -(num-num/100*100)/10*10+48)
      endif
         evlfilename(13:13) = '_'
         evlfilename(1:8)='photons/'
!
      return
      end

!     Returns name of the local record file on node num
      character*43 function recordfilename(recordfile, num)
      implicit none
      character*30 recordfile
      integer num
!
      character*43 localrecordfile
!
      recordfilename(14:43) = recordfile(1:30)
         recordfilename(9:9) = 'p'
      if(num.lt.10) then
         recordfilename(10:10) = '0'
         recordfilename(11:11) = '0'
         recordfilename(12:12) = char(num+48)

      else if(num.lt.100) then
         recordfilename(10:10) = '0'
         recordfilename(11:11) = char(num/10+48)
         recordfilename(12:12) = char(num-num/10*10+48)
      else
         recordfilename(10:10) = char(num/100+48)
         recordfilename(11:11) = char((num-num/100*100)/10+48)
         recordfilename(12:12) = char(num-num/100*100&
     &        -(num-num/100*100)/10*10+48)
      endif
         recordfilename(13:13) = '_'
         recordfilename(1:8) = 'restart/'
!
      return
      end

!
!
!
!
!
!
!
!
!     The subroutine xec_bcast broadcasts ncycle, time, 
!     other variables to the slaves at the beginning of 
!     every time step.  xec_add gathers Ed_ref, Ed_in, and ecens
!     from the slaves and sumes them for the master process.
!     xec_bcast is called in xec2d.
!     J. Finke, 2 May 2005.
      subroutine xec_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer status(MPI_STATUS_SIZE)

      real etime, elapse(2)
!

!
!     broadcast nc
      call MPI_BCAST(ncycle, 1, MPI_INTEGER, master, &
     &      MPI_COMM_WORLD, ierr)
!     broadcast times
      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
!     broadcast elapsed times
      call MPI_BCAST(etotal, 1, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(etotal_old, 1, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
!     broadcast fpair
      call MPI_BCAST(f_pair, jmax*kmax, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!     broadcast mubins
      call MPI_BCAST(mu, nmumax, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!     broadcast gg_abs
      call MPI_BCAST(E_gg, nmumax, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!     broadcast reflection
      call MPI_BCAST(P_ref, n_ref*n_ref, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(E_ref, n_ref, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(W_abs, n_ref*n_ref, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!
!
!
      return
      end

      subroutine xec_add
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer int_temp
      integer status(MPI_STATUS_SIZE)
      integer j, k
      double precision Ed_temp

      do 10 k=1, nr
!     collect ddh from slave nodes to master node
         call MPI_REDUCE(Ed_ref(k), Ed_temp, 1, MPI_DOUBLE_PRECISION,&
     &        MPI_SUM, master, MPI_COMM_WORLD, ierr)
         if(myid.eq.master) Ed_ref(k) = Ed_temp
         call MPI_REDUCE(Ed_in(k), Ed_temp, 1, MPI_DOUBLE_PRECISION,&
     &        MPI_SUM, master, MPI_COMM_WORLD, ierr)
         if(myid.eq.master) Ed_in(k) = Ed_temp

         do 20 j=1, nz
!           collect ecens to master node (previously done in imcgen)
            if((ncycle.le.1).or.(fp_sw.eq.0)) then
            call MPI_REDUCE(ecens(j,k), Ed_temp, 1, MPI_DOUBLE_PRECISION&
     &           , MPI_SUM, master, MPI_COMM_WORLD, ierr)
            if(myid.eq.master) ecens(j,k) = Ed_temp
            endif
            call MPI_REDUCE(npcen(j,k), int_temp, 1, MPI_INTEGER, &
     &           MPI_SUM, master, MPI_COMM_WORLD, ierr)
            if(myid.eq.master)npcen(j,k) = int_temp
!           npcen for each node is now meaningless since MC particles changed
!           nodes in imcredist.
!
 20      continue
 10   continue

      return
      end
!
!
!
!
!     This subroutine collects variables from the slaves needed in 
!     the graphics subroutine
!     J. Finke, 19 July 2005
      subroutine graphics_collect
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer n, lgp
      double precision temporary
!
!
      do 10 n = 1, nmu
         do 20 lgp=1,nph_lc
!            fout(n,lgp) = 0.d0 
            call MPI_REDUCE(edout(n,lgp), temporary, 1, &
     &           MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, &
     &           ierr)
            if(myid.eq.master) edout(n,lgp) = temporary
 20      continue
         do 30 lgp=1, nphomax
            if(myid.eq.master) fout(n,lgp) = 0.d0
            call MPI_REDUCE(fout(n,lgp), temporary, 1, &
     &           MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, &
     &           ierr)
            if(myid.eq.master) fout(n,lgp) = temporary
 30      continue
 10   continue
!
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:19:37 EDT 2006
! version: 2
! Name: J. Finke
! Fixed bug with fp_sw.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun 16 12:33:24 EDT 2006
! version: 3
! Name: J. Finke
! xec now calls imcfield.     
