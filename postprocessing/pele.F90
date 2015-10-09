!	postprocessing of the electron distribution from the Monte Carlo/Fokker Planck code.
      program main
        implicit none
	integer jmax,kmax,num_nt, nunit_fnt
	parameter (jmax = 99)
	parameter (kmax = 99)
	parameter (num_nt = 200)
	double precision f_ele(jmax,kmax,num_nt), n_e(jmax,kmax), E_e(jmax,kmax), Eloss_sy(jmax,kmax)
	double precision gnt(num_nt), B_field(jmax,kmax), theta(jmax,kmax),fnt_tot(num_nt),fnt_esc(num_nt), &
               spec_e(num_nt), f_escape(num_nt)
	double precision time, emax, nmax2date, emax2date, smax2date, gmax2date(10), ztot, rtot, dz, dr, &
               u_tot, u_tot_old, u_increase, sum_g_1, ne
        integer i,j,k, ncycle,lggam, gi(10)
        integer dum, nz, nr 
	real etime, elapse(2), et0
	character*4 cstep
        character*30 fmtele, felename, emapname
        character*80 buffer
	logical f_exist

	nunit_fnt=22
	nmax2date = 0.d0
	emax2date = 0.d0
	smax2date = 0.d0
	gmax2date(:) = 0.d0
!       read in E_e.dat
	open(nunit_fnt, file='E_e.dat', status='unknown')
         read(nunit_fnt,'("# Volume size: Z=",e14.7," R=",e14.7)')ztot,rtot
         read(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
         read(nunit_fnt,'("# Energy grid for electron spectrum")')
         do i = 1, num_nt
           read(nunit_fnt, '(e14.7)') gnt(i)
         enddo
        close(nunit_fnt)
	dz=ztot/nz
	dr=rtot/nr

!       determine the sample energies
	lggam = 1
	do i=1,num_nt-1
	  if(gnt(i)+1.lt.10**lggam.and.gnt(i+1)+1.gt.10**lggam)then
	    gi(lggam)=i
	    lggam=lggam+1
	  endif
	enddo

	felename = 'fnt_0001.dat'
        cstep(1:1)='_'
	emapname = 'g5map_0001.dat'
	do ncycle= 1, 9999
          cstep(1:1) = char(48 + int(ncycle/1000))
          cstep(2:2) = char(48 + int(ncycle/100)-10*int(ncycle/1000))
          cstep(3:3) = char(48 + int(ncycle/10)-10*int(ncycle/100))
          cstep(4:4) = char(48 + ncycle-10*int(ncycle/10))

!       read in fnt_***.dat
	  felename(1:3) = 'fnt'
	  felename(5:8) = cstep
	  inquire( file=felename,exist=f_exist)
	  if(.not.f_exist)exit
	  emapname(7:10) = cstep
	  inquire( file=emapname,exist=f_exist)
          if(f_exist)cycle

          open(nunit_fnt, file=felename, status='unknown')
           read(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') dum,time
           write(fmtele,'("("a9","i4"e14.7)")'),'"#B=    "',nz*nr
           read(nunit_fnt, fmtele)(((B_field(j,k)), k=1,nr), j=1,nz)
           write(fmtele,'("("a9","i4"e14.7)")'),'"#theta="',nz*nr
           read(nunit_fnt, fmtele)(((theta(j,k)), k=1,nr), j=1,nz)
           write(fmtele,'("("i4"e14.7)")'),nz*nr+2
           do i = 1, num_nt
             read(nunit_fnt,fmtele)(((f_ele(j,k,i)), k=1,nr),j=1,nz),fnt_tot(i),f_escape(i)
           enddo
	  close(nunit_fnt)
!        ne =0.d0
!        do i =1, num_nt-1
!          ne = ne + fnt_tot(i)*(gnt(i+1)-gnt(i))
!        enddo
!        write(*,*)'ne=',ne

!       calculate escaped (from the center) particle spectrum
	  fnt_esc(:)=0.d0 !fnt_tot(:)*ztot*rtot**2!*pi
	  do j=1,nz
	    do k=1,nr
	      if(k.gt.1.and.j.ne.12.and.j.ne.13)fnt_esc(:)=fnt_esc(:)+f_ele(j,k,:)*dz*((k*dr)**2-((k-1)*dr)**2)!*pi
	    enddo
	  enddo
	  fnt_esc(:)=fnt_esc(:)/ztot/rtot**2!/pi  They all have pi factor. Canceled.

!       calculate spectral index of the total electron spectrum
         do i = 1, num_nt-1
           spec_e(i) = dlog10(fnt_tot(i+1)/fnt_tot(i))/dlog10((gnt(i+1)+1.d0)/(gnt(i)+1.d0))
         enddo

!	calculate density n_e, energy density E_e, and synchrotron energy loss (erg/s/cm^3)
          n_e(:,:) = 0.d0
          E_e(:,:) = 0.d0
	  do j=1,nz
	  do k=1,nr
	    sum_g_1 = 0.d0
            do i = 1, num_nt-1
              n_e(j,k) = n_e(j,k) + (gnt(i+1) - gnt(i))*f_ele(j,k,i)
              E_e(j,k) = E_e(j,k) + (gnt(i+1) - gnt(i))*(gnt(i)+1.d0)*f_ele(j,k,i)
              sum_g_1 = sum_g_1 + ((gnt(i)+1.d0)**2 -1.d0)*f_ele(j,k,i)*(gnt(i+1)-gnt(i))
            enddo
            Eloss_sy(j,k) = 1.058d-15*B_field(j,k)**2.d0*sum_g_1
	  enddo
          enddo
          u_tot_old = u_tot
          u_tot = 0.d0
          do i = 1,num_nt-1
            u_tot = u_tot + (gnt(i+1) - gnt(i))*(gnt(i)+1.d0)*fnt_tot(i)
          enddo
          u_increase = (u_tot-u_tot_old)/u_tot_old

          write(fmtele,'("("i3"e14.7)")'),2*nr

!       output sample electron files
	  felename(1:3) = 'smp'
          open(nunit_fnt, file=felename, status='unknown')
	   write(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') dum,time
           write(nunit_fnt,'("# total energy density (erg/cm^3):",e14.7," increased by ",e10.3)') u_tot,u_increase
           write(fmtele,'("("i4"e14.7,e14.6)")'),nz*3+2
	   do i = 1, num_nt
	     write(nunit_fnt,fmtele)(((f_ele(j,min(nr,k*nr/2+1),i)), k=0,2), j=1,nz),&
                  fnt_tot(i),f_escape(i),spec_e(i)
	   enddo
	  close(nunit_fnt)

!	output particle density
	  emapname(1:2) = 'ne'
	  emapname(7:10) = cstep
	  emax=maxval(n_e(:,:))
	  if(nmax2date.lt.emax)nmax2date = emax
	  open(nunit_fnt, file=emapname, status='unknown')
           write(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') ncycle,time
           write(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
	   write(nunit_fnt, '("# maximum value = ",e14.7," (log:",f5.2,") max to date =",e14.7)')emax,log10(emax),nmax2date
           do j=1,nz
             write(nunit_fnt,fmtele)(n_e(j,k), k=nr,1,-1),(n_e(j,k), k=1,nr)
	   enddo
	  close(nunit_fnt)

!	output energy density
	  emapname(1:2) = 'ee'
	  emapname(7:10) = cstep
	  emax=maxval(E_e(:,:))
	  if(emax2date.lt.emax)emax2date = emax
	  open(nunit_fnt, file=emapname, status='unknown')
           write(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') ncycle,time
           write(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
           write(nunit_fnt, '("# maximum value = ",e14.7," (log:",f5.2,") max to date =",e14.7)')emax,log10(emax),emax2date
           do j=1,nz
             write(nunit_fnt,fmtele)(E_e(j,k), k=nr,1,-1),(E_e(j,k), k=1,nr)
	   enddo
	  close(nunit_fnt)

!	output synchrotron loss rate
	  emapname(1:2) = 'se'
	  emapname(7:10) = cstep
	  emax=maxval(Eloss_sy(:,:))
	  if(smax2date.lt.emax)smax2date = emax
	  open(nunit_fnt, file=emapname, status='unknown')
           write(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') ncycle,time
           write(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
           write(nunit_fnt, '("# maximum value = ",e14.7," (log:",f5.2,") max to date =",e14.7)')emax,log10(emax),smax2date
           do j=1,nz
             !write(nunit_fnt,fmtele)(dmax1(Eloss_sy(j,k),1.d-20), k=nr,1,-1),(dmax1(Eloss_sy(j,k),1.d-20), k=1,nr)
             write(nunit_fnt,fmtele)(Eloss_sy(j,k), k=nr,1,-1),(Eloss_sy(j,k), k=1,nr)

	   enddo
	  close(nunit_fnt)


	  emapname(1:1) = 'g'
!	out put differential energy density at certain sample energies
	  do lggam=1,5
	    emapname(2:2) = char(48 + lggam)
	    emax=maxval(f_ele(:,:,gi(lggam)))
	    if(gmax2date(lggam).lt.emax)gmax2date(lggam) = emax
	    open(nunit_fnt, file=emapname, status='unknown')
             write(nunit_fnt,'("# ",i4,". time step: t = ",e14.7," s")') ncycle,time
             write(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
             write(nunit_fnt, '("# maximum value = ",e14.7," (log:",f5.2,") max to date =",e14.7)')emax,log10(emax),gmax2date(lggam)
             do j=1,nz
               write(nunit_fnt,fmtele)(f_ele(j,k,gi(lggam)), k=nr,1,-1), (f_ele(j,k,gi(lggam)), k=1,nr)
	     enddo
	    close(nunit_fnt)
	  enddo
	enddo

        end

     
