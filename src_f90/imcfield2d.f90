!     Each processor reads in large amount of existing photons 
!     from the previous time step and tracks them.
!     comment by Xuhui Chen 2012-07-15
      subroutine imcfield2d 
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer l, sender
!
!
!
!     master node
      if(myid.eq.master) then
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      write(*,*) 'myid=', myid, ' calling field_send_synch'
      call field_send_synch
!      write(*,*) 'myid=', myid, ' calling field_calc'
!      call field_calc
!      write(*,*) 'myid=', myid, ' called field_calc'
!
      do 10 l = 1, (numprocs-1)
         write(*,*) 'l=', l, ' calling field_recv, numprocs=', numprocs
         call field_recv(sender)
         write(*,*) 'l=', l, ' sender=', sender, ' called field_recv'
 10   continue
!
!
!     slave node
      else
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      write(*,*) 'myid=', myid, ' calling field_recv_synch'
      call field_recv_synch
!
      write(*,*) 'myid=', myid, ' calling field_calc'
      call field_calc
      write(*,*) 'myid=', myid, ' calling field_send'
      call field_send
      write(*,*) 'myid=', myid, ' called field_sendq'
!
!
!     end of imcfield2d
      endif
!      write(*,*) 'myid=', myid, ' in imcfield ecens=', ecens(1,1)
      return
      end
!
!
!
!
!     Read in photons from census file (previous time step) and track them.
!     J. Finke, 2 May 2005
      subroutine field_calc
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
       integer ncrem, nread, lwad, lwai, ndxin
!
       double precision cdt
!
!
!
!      
       ncrem = ndxout
       ndxout = 0 ! The re-initialization of ndxout. used to be in imcgen Xuhui 08/15/11
!       
!
       if (ncrem.eq.0) goto 900
!
       if (ncrem.ge.ucens) then
          write(*,*)'too many photons needed to read'
          stop
          nread = ucens
          ncrem = ncrem - ucens
       else
          nread = ncrem
          ncrem = 0
       endif
!
     
       dbufin(:)=dbufout(:)
       ibufin(:)=ibufout(:)
       ndxin = 0
!  
 110   continue
!      
         lwad = 6*ndxin
         lwai = 6*ndxin
         rpre = dbufin(lwad+1)
         zpre = dbufin(lwad+2)
         wmu = dbufin(lwad+3)
         phi = dbufin(lwad+4)
         ew = dbufin(lwad+5)
         xnu = dbufin(lwad+6)
         jgpsp = ibufin(lwai+1)
         jgplc = ibufin(lwai+2)
         jgpmu = ibufin(lwai+3)
         jph = ibufin(lwai+4)
         kph = ibufin(lwai+5)
         emitype = ibufin(lwai+6)
         call initialize_rand(jph,kph) ! immediately returned when rand_switch is 2
         dcen = c_light*dt(1)
!
         if (abs(wmu).gt.1.d0+1.d-15)write(*,'(e24.17)')&
     &    'What?! wmu too extreme! Eta=',wmu
         if (wmu.gt.1.d0) wmu = 1.d0
         if (wmu.lt.-1.d0) wmu = -1.d0
!         if (wmu.gt.0.99999999d0) wmu = 0.99999999d0
!         if (wmu.lt.-0.99999999d0) wmu = -0.99999999d0
!
         ndxin = ndxin+1
!
         call imctrk2d(-1)
!         
         if(ndxin .ge. nread) then
            if (ncrem.le.0) goto 900
!            if (ncrem.ge.ucens) then
!               nread = ucens
!               ncrem = ncrem - ucens
!            else
               nread = ncrem
               ncrem = 0
!            endif
!           call read_cens(nread)
            ndxin = 0
         endif
      
         goto 110
!
 900    continue
!
        return
        end
!
!
!
!
!     This subroutine is here for synchronization purposes.  In the
!     future, if something is needed for sending in imcfield,
!     it may be added here.
!     J. Finke, 9 May 2005
      subroutine field_send_synch
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer node
!
      do 10 node=1, (numprocs-1)
         call MPI_SEND(MPI_BOTTOM, 0, MPI_INTEGER, node, myid,&
     &        MPI_COMM_WORLD, ierr)
 10   continue
!
      return
      end
!
!
!
!     This subroutine is here for synchronization purposes.  In the
!     future, if something is needed for recieving in imcfield,
!     it may be added here.
!     J. Finke, 9 May 2005
      subroutine field_recv_synch
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer dummy
      integer status(MPI_STATUS_SIZE)
!
!
      call MPI_RECV(dummy, 1, MPI_INTEGER, master, MPI_ANY_TAG,&
     &     MPI_COMM_WORLD, status, ierr)
!
!
      return
      end
!
!
!
!
!     Sends number of the node that completed its field 
!     calculation.  Paired with field_recv.
!     J. Finke, 2 May 2005
      subroutine field_send
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      call MPI_SEND(master, 1, MPI_DOUBLE_PRECISION,&
     &     master, myid, MPI_COMM_WORLD, ierr)
!
      return
      end
!
!
!
!     Recieves number of the node that completed its field 
!     calculation.  Paired with field_send.
!     J. Finke, 2 May 2005
      subroutine field_recv(sender)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer sender
!
      integer status(MPI_STATUS_SIZE)
      double precision dummy
!
!
      call MPI_RECV(dummy, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
     &     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      sender = status(MPI_TAG)
!
      return
      end
!
!
!
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun 16 12:26:36 EDT 2006
! version: 2
! Name: J. Finke
! fibran seed is now stored in census files. 
! imcfield is now determinable.     
