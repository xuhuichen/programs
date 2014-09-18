!     redistribute the photons in the field, so that each nodes handles
!     similar number of field MC particles.
!     maximum nodes number 1000
!     Xuhui Chen, August 2011
      subroutine imcredist
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'

      integer nctot, nc_temp, nc_mod, nc_uneven(1000), nc_even(1000),&
     &        i, islave, overload, oltot, ol6, tot6, nd6
      integer status(MPI_STATUS_SIZE)
!      integer ibuf_r(iucens), icollected(iucens)
!      double precision dbuf_r(ducens), dcollected(ducens)


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!cc   master node
      if(myid.eq.master)then
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cc   collect the number of photons each slave node has
      nctot = 0
      do i = 1, numprocs-1
        call MPI_RECV(nc_temp, 1, MPI_INTEGER, MPI_ANY_SOURCE,&
     &       MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        islave = status(MPI_SOURCE)
        nc_uneven(islave) = nc_temp
        nctot = nctot + nc_temp
      enddo
       write(*,*) 'in imcredist nctot=', nctot
!cc   calculate the number of photons each slave node should have
      nc_mod = mod(nctot,numprocs-1)
      nc_even(:)=nctot/(numprocs-1)
      do islave = 1, nc_mod
         nc_even(islave)=nc_even(islave)+1
      enddo
!
      do islave = 1, numprocs-1
         overload = nc_uneven(islave)-nc_even(islave)
         call MPI_SEND(overload, 1, MPI_INTEGER, islave, 2,&
     &        MPI_COMM_WORLD, ierr)
      enddo
!cc   collect all the overload photons
      oltot = 0
      do islave = 1, numprocs-1
        overload = nc_uneven(islave)-nc_even(islave)

        if(overload.gt.0)then
          ol6 = overload*6
          tot6 = oltot*6

!cc     bufout in master load is used to collect overload photons
          call MPI_RECV(ibufout(tot6+1), ol6, MPI_INTEGER, islave,&
     &         101, MPI_COMM_WORLD, status, ierr)
!          icollected((tot6+1):(tot6+ol6))=ibuf_r(1:ol6)
          call MPI_RECV(dbufout(tot6+1), ol6, MPI_DOUBLE_PRECISION,&
     &         islave, 102, MPI_COMM_WORLD, status, ierr)
!          dcollected((tot6+1):(tot6+ol6))=dbuf_r(1:ol6)

          oltot = oltot+overload
          if(oltot.gt.ucens)then
            write(*,*)'too many photons collected in master node'
            stop
          endif

        endif
      enddo
!cc   distribute all the overload photons to underloaded nodes
      do islave = 1, numprocs-1
        overload = nc_uneven(islave)-nc_even(islave)

        if(overload.lt.0)then 
          overload = -overload   ! overload now means underload
          ol6 = overload*6 
          tot6 = oltot*6

!          ibuf_r(1:ol6) = icollected((tot6-ol6+1):tot6)
          call MPI_SEND(ibufout(tot6-ol6+1), ol6, MPI_INTEGER, islave,&
     &         201, MPI_COMM_WORLD, ierr)
!          dbuf_r(1:ol6) = dcollected((tot6-ol6+1):tot6)
          call MPI_SEND(dbufout(tot6-ol6+1), ol6, MPI_DOUBLE_PRECISION,&
     &         islave, 202, MPI_COMM_WORLD,ierr)

          oltot = oltot-overload
        endif
      enddo
      write(*,*)'photons left in master node:',oltot
!cc   slave node
      else
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!cc   send the number of photons each node has
      call MPI_SEND(ndxout, 1, MPI_INTEGER, master, myid,&
     &     MPI_COMM_WORLD, ierr)
!cc   recieve the overload of photons each node has
      call MPI_RECV(overload, 1, MPI_INTEGER, master, 2,&
     &     MPI_COMM_WORLD, status, ierr)
!cc   send the overloaded photons to the master node
      if(overload.gt.0)then
        ol6 = overload*6
        nd6 = ndxout*6

!        ibuf_r(1:ol6) = ibufout((nd6-ol6+1):nd6) 
        call MPI_SEND(ibufout(nd6-ol6+1), ol6, MPI_INTEGER, master,&
     &       101, MPI_COMM_WORLD, ierr)
!        dbuf_r(1:ol6) = dbufout((nd6-ol6+1):nd6) 
        call MPI_SEND(dbufout(nd6-ol6+1), ol6, MPI_DOUBLE_PRECISION,&
     &       master, 102, MPI_COMM_WORLD, ierr)

        ndxout = ndxout-overload
      endif
!cc   underloaded nodes recieve the redistributed photons from the master node
      if(overload.lt.0)then
        overload = - overload ! overload now means underload
        ol6 = overload*6
        nd6 = ndxout*6

        call MPI_RECV(ibufout(nd6+1), ol6, MPI_INTEGER, master,&
     &       201, MPI_COMM_WORLD, status, ierr)
!        ibufout((nd6+1):(nd6+ol6))=ibuf_r(1:ol6)
        call MPI_RECV(dbufout(nd6+1), ol6, MPI_DOUBLE_PRECISION, &
     &       master, 202, MPI_COMM_WORLD, status, ierr)
!        dbufout((nd6+1):(nd6+ol6))=dbuf_r(1:ol6)

        ndxout = ndxout+overload
      endif
!      write(*,*)'myid=',myid,'ndxout=',ndxout
      endif

      return
      end
