        program main
        implicit none
	!include 'mpif.h'
        double precision array(0:1,1:2), ran1, b(1000000),x
        integer a(2)
        integer i,j,k,seed, myid, numprocs, ierr, mid
        integer, parameter  :: dp = selected_real_kind(8)
	real etime, elapse(2), et0
        character*30 fmt, str
        !call MPI_INIT(ierr)
        !call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
        !call MPI_COMM_RANK(MPI_COMM_WORLD, i, ierr)
        !write(*,*) 'myid=',myid,'i=',i
	str = 'abcdefg'
	str(2:3)='rt'
        seed=10020
        array(0,1)=11
        array(0,2)=-12
        array(1,1)=21
        array(1,2)=5
        do i=1,1000000
          b(i)= -i*0.1    
        enddo

 
        do j =1,100
        x = -3447.0
        a(1) = 1
        a(2) = 1000000
        k = 0
        do while (x.le.b(a(1)).and.x.ge.b(a(2)))
          k = k +1
          mid =(a(1)+a(2))/2
          if(x.lt.b(mid))then
            a(1)=mid
          else
            a(2)=mid
            if(x.eq.b(mid))exit
          endif
          if(a(2)-a(1).eq.1)exit
        enddo
        enddo
        et0=etime(elapse)
        write(*,*)'a(2)=',a(2),'et0=',et0, 'k=',k


        do j =1,100
        x = -3447.3
        a(1) = 1
        a(2) = 1000000
        do while (x.le.b(a(1)))
          a(1)=a(1)+1
        enddo
        enddo
        et0=etime(elapse)
        write(*,*)'a(1)=',a(1),'et0=',et0
        end

     
        double precision FUNCTION ran1(idum)
        implicit none
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        integer counter
        REAL AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     @          NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
c        PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
c     @          NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=AM,RNMX=1.d0-EPS)
        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        save counter
        DATA iv /NTAB*0/, iy /0/
        data counter /0/
c      write(4,*) 'iy=',iy
        counter = counter + 1
        if (idum.le.0.or.iy.eq.0) then
                idum=max(-idum,1)
                do j=NTAB+8,1,-1
                        k=idum/IQ
                        idum=IA*(idum-k*IQ)-IR*k
                        if (idum.lt.0) idum=idum+IM
                        if (j.le.NTAB) iv(j)=idum
                end do
                iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran1=int(min(AM*iy,RNMX)*1.e7)/1.d7
c        write(4,*) 'ran1=', ran1
c        write(4,*) 'idum=',idum
c        write(4,*) 'iv1=', iv(1)
        return
        end

