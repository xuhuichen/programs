!     Broadcasts to the slaves data needed for the volume calculation.
!     J. Finke, 18 May 2005
      subroutine vol_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer status(MPI_STATUS_SIZE)
      double precision esurf
!
!
!     broadcast integers
!
!     common block timeindex
!      call MPI_BCAST(t, 1, MPI_INTEGER, master, 
!     1      MPI_COMM_WORLD, ierr)
!     common block izones
!      call MPI_BCAST(nz, 1, MPI_INTEGER, master, 
!     1      MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(nr, 1, MPI_INTEGER, master, 
!     1      MPI_COMM_WORLD, ierr)
!
!     broadcast doubles
!
!     common block zones
!      call MPI_BCAST(z, jmax, MPI_DOUBLE_PRECISION, master, 
!     1      MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(r, kmax, MPI_DOUBLE_PRECISION, master, 
!     1      MPI_COMM_WORLD, ierr)
!     common block times
!      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION, master, 
!     1      MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(dt, 2, MPI_DOUBLE_PRECISION, master, 
!     1      MPI_COMM_WORLD, ierr)
!
!     common block eweights
!      write(*,*) 'myid=',myid,' in vol_bcast'
      call MPI_BCAST(ewsv, jmax*kmax, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)   
!     common block vol_em
      call MPI_BCAST(E_ph, n_vol, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(eps_th, n_vol*jmax*kmax, MPI_DOUBLE_PRECISION, 
!     1      master, MPI_COMM_WORLD, ierr)
!
!
      return
      end
!
!
!
!     Creates struct type for sending jobs to the slaves.
!     J. Finke, 18 May 2005
      subroutine vol_create_job_type
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!

      integer lengths(2), displacements(2), types(2), address(2)

!

!
      lengths(1) = int_ans_size
      lengths(2) = double_ans_size
      call MPI_ADDRESS( integer_ans(1), address(1), ierr)
      call MPI_ADDRESS( double_ans(1), address(2), ierr )
      displacements(1) = 0
      displacements(2) = address(2) - address(1)
      types(1) = MPI_INTEGER
      types(2) = MPI_DOUBLE_PRECISION

      call MPI_TYPE_STRUCT(2, lengths, displacements, types, &
     &     vol_job_type, ierr)
      call MPI_TYPE_COMMIT(vol_job_type, ierr)
!
      return
      end
!
!
!
!     Sends volume job for zone to available_proc.  
!     Paired with vol_recv_job.
!     J. Finke, 18 May 2005
      subroutine vol_send_job(available_proc, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer available_proc, zone
!
!
      integer j, k, i, i_a
!
      integer status(MPI_STATUS_SIZE)
!
!
      call get_j_k(zone, j, k, nr)
!
!     put items to send in arrays
!
!     common block random
      integer_ans(1) = seeds(j,k)
!     common block ph_numbers
      integer_ans(2) = nsv(j,k)
!     common block zones
      double_ans(1) = zsurf(j,k)
!     common block vol_em
      double_ans(2) = Eloss_th(j,k)
      double_ans(3) = Eloss_tot(j,k)
      double_ans(4) = B_field(j,k)
      double_ans(5) = theta_local(j,k)
!      do 10 i = 1, n_vol
!         double_ans(3+i) = eps_tot(i,j,k)
! 10   continue
!      do i=1,n_vol
!         do i_a=1,num_a
!           double_ans(3+n_vol+(i-1)*num_a+i_a) = eps_a(i_a,i,j,k)
!         enddo
!      enddo
!
      call MPI_SEND(integer_ans, 1, vol_job_type, available_proc, &
     &     zone, MPI_COMM_WORLD, ierr)
!
!
      return
      end
!

      subroutine vol_send_end_signal(available_proc,zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer node
!
!
      integer end_signal, zone, available_proc
      integer status(MPI_STATUS_SIZE)
!
!
      end_signal = jmax*kmax+1
!      write(*,*) 'myid=', myid, ' in send_end_signal, end_signal=', 
!     1     end_signal
      write(*,*)'sender=',available_proc,&
     & 'send end. end_signal=',end_signal
      call MPI_SEND(integer_ans, 1, vol_job_type, available_proc, &
     &     end_signal, MPI_COMM_WORLD, ierr)
            write(*,*)'sender=',available_proc,'end successfully sent'
!
      return
      end
!
!     Recieves volume job for zone from master node.  
!     Paired with vol_send_job.
!     J. Finke, 18 May 2005
      subroutine vol_recv_job(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer zone
!
!
      integer j, k, i, i_a
!
      integer status(MPI_STATUS_SIZE)
      integer end_signal
!
      integer lengths(2), displacements(2), types(2), address(2)
!
!
      integer v
      logical flag
!
!
      end_signal = jmax*kmax+1
!
!
!      write(*,*)'myid=',myid,'will recieve'
      call MPI_RECV(integer_ans, 1, vol_job_type, master,MPI_ANY_TAG,&
     &     MPI_COMM_WORLD, status, ierr)
!      write(*,*)'myid=',myid,'recieved'
      zone = status(MPI_TAG)
!      write(*,*)'zone tag=',zone
      if(zone.eq.end_signal) return
      call get_j_k(zone, j, k, nr)
!
!     common block random
      seeds(j,k) = integer_ans(1)
!     common block ph_numbers
      nsv(j,k) = integer_ans(2)
!     common block zones
      zsurf(j,k) = double_ans(1)
!     common block vol_em
      Eloss_th(j,k) = double_ans(2)
      Eloss_tot(j,k) = double_ans(3)
      B_field(j,k) = double_ans(4)
      theta_local(j,k) = double_ans(5)

!      do 10 i = 1, n_vol
!         eps_tot(i,j,k) = double_ans(3+i)
! 10   continue
!      do i=1,n_vol
!         do i_a=1,num_a
!           eps_a(i_a,i,j,k) = double_ans(3+n_vol+(i-1)*num_a+i_a)
!         enddo
!      enddo
!
      return
      end
!
!
!
!     Sends result of volume calculation for zone to master node.
!     The results include kappa_tot calculated by volume_em
!     Paired with vol_recv_result.
!     J. Finke, 18 May 2005
      subroutine vol_send_result(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer zone
      integer j, k
      double precision ans(n_vol)
!
      call get_j_k(zone, j, k, nr)
      ans(:) = kappa_tot(:,j,k)
!
!      write(*,*)'myid=',myid,'will send result'
      call MPI_SEND(ans, n_vol, MPI_DOUBLE_PRECISION, master, zone,&
     &     MPI_COMM_WORLD, ierr)
!      write(*,*)'myid=',myid,'result sent'
!
!
      return
      end
!
!
!
!     Recieves result of volume calculation for zone from sender.
!     paired with vol_send_result.
!     J. Finke, 18 May 2005
      subroutine vol_recv_result(sender, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer sender, zone
      integer j, k
      double precision ans(n_vol)
!
      integer junk
      integer status(MPI_STATUS_SIZE)
!

!
!
!      write(*,*)'myid=',myid,'will recieve result'
      call MPI_RECV(ans, n_vol, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
     &     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
!      write(*,*)'myid=',myid,'result recieved'
!
      sender = status(MPI_SOURCE)
      zone = status(MPI_TAG)

      call get_j_k(zone, j, k, nr)
      kappa_tot(:,j,k) = ans(:)
!
!
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:36:26 EDT 2006
! version: 2
! Name: J. Finke
! Changed common block 'random'.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Thu Oct  5 13:22:38 EDT 2006
! version: 3
! Name: J. Finke
! Changed routines vol_create_job_type, vol_send_job, and vol_recv_job 
! so that 
! eps_tot is sent to the slave nodes.  
