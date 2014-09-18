      program main
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!     MPI variables
      logical continued
      integer status(MPI_STATUS_SIZE)
!
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      write(*,*) 'myid=',myid
      master = 0      
      inquire(file='p000_misc.dat',exist=continued)

      if(continued)then
        if(myid.eq.master)open(unit=4, file='log.txt', access='append' )
        call read_record

      else
        if(myid.eq.master) then
!
          open(unit=4, file='log.txt')
          write(4,*) 'Number of Processors:  ', numprocs
!
          call reader(0)
!
        endif
        call setup(0)
      endif
!
      call xec
      if(myid.eq.master) then
!
      close(4)
!
      endif
      call MPI_FINALIZE(ierr)
      write(*,*) 'myid=',myid,' end of compton2d'
      end
!
!

!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun 16 16:29:29 EDT 2006
! version: 2
! Name: J. Finke
! Prints number of processors to the log file 
! at beginning of program.     
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! 2013.2nd quater
! Xuhui
! changed the first splitting in imctrk to a 
! deterministic description, i.e. increased
! scattering probability of 1/spl1 of the ew, and
! the other ew has no chance of being scattered
!cccccccccccccccccc
! 2013.4th quater
! Xuhui Chen
! added the option to use helical magnetic field,
! which gives anisotropic synchrotron emissivity.
! But synchrotron self-absorption still does not
! take the anisotropy into account.
! ------------------
! Mean while
! Haocheng Zhang 
! developed another code that can
! use the electron and B information of this MCFP code
! and calculate time-dependent synchrotron polarization.
!ccccccccccccccccccc
! 2014.02.05
! Xuhui Chen
! added particle diffusion. It is implemented in update
! subroutine.
!ccccccccccccccccccc
! 2014.02.07
! Xuhui Chen
! volume_em is moved from imcgen to imcvol_para so that
! this part of the code is also parallelized.
