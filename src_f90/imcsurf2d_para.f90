      subroutine imcsurf2d
!     This subroutine creates new photons based on the input energy at four surfaces 
!     (inner, outer, upper, lower) from external source
!     and begins to track them. !comment by Xuhui Chen 7/8/08
!     

  
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer i, n, js, ks, nxs
      integer npsurf, isurfmu
!     MPI variables
      integer num_sent, end_signal, sender
      integer zone
!
      double precision esurf
      double precision psi, wv1, wv2, tpl
!
      double precision fibran, t_average, x1, x2
!
!

!
!
!     MPI Initialization 
!      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      num_sent = 0
      end_signal = jmax*kmax+1
!      write(*,*) 'myid=', myid, 'end_signal=', end_signal
!
!
!
!     master node part
      if(myid.eq.master) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!
!      write(*,*) 'beginning of imcsurf2d'
      isurfmu = 1
      npsurf = 0
      esurf = 0.d0
!
!      do 120 n = 1, nphtotal
! 120  F(n) = 0.d0
!
!
!    Scan through time bins to find 
!      the correct time index t 
!
       if (ncycle.eq.0) then
          ti = 1
          goto 135
       endif
       t_average = time + 5.d-1*dt(1)
!
      do 130 ti = 1, ntime
         if (t1(ti).gt.t_average) goto 135
 130  continue
 135  continue
!
!
!
!     Vertical boundaries
!
      call z_surf_bcast
!
!     This sends the first round of zones to the slaves for processing.
         do 300 js = 1, min(nz, numprocs-1)
            zone = js
            call z_surf_send_job(js, zone)
!            write(*,*) 'myid=', myid, ' sent zone=', zone
            num_sent = num_sent + 1
 300     continue
!
!     As slaves complete processing a zone, this recieves the results
!     and sends the next zone to the slaves.
         do 302 js=1, nz
!            write(*,*) 'js=', js
            call z_surf_recv_result(sender, zone)
            npsurf = npsurf + npsurfz(js)
            if(num_sent.lt.nz) then
!           if there are still more zones to send, send next zone
!           to the available node.
               zone = num_sent + 1
               call z_surf_send_job(sender, zone)
               num_sent = num_sent + 1
            else
               call surf_send_end_signal(sender)
            endif
 302     continue
!
!
!
!
!     Radial boundaries
!
         num_sent = 0
         call r_surf_bcast
!
!     This sends the first round of zones to the slaves for processing.
         do 400 ks = 1, min(nr, numprocs-1)
            zone = ks
!            write(*,*) 'myid=', myid, ' before send zone=', zone,
!     1           ' nsurfl(zone)=', nsurfl(zone)
            call r_surf_send_job(ks, zone)
!            write(*,*) 'myid=', myid, ' sent zone=', zone
            num_sent = num_sent + 1
 400     continue
!
!     As slaves complete processing a zone, this recieves the results
!     and sends the next zone to the slaves.
         do 402 ks=1, nr
            call r_surf_recv_result(sender, zone)
!            write(*,*) 'master recieved sender=', sender, ' zone=', zone
!     1           , ' ks=', ks
            npsurf = npsurf + npsurfr(ks)
            if(num_sent.lt.nr) then
!           if there are still more zones to send, send next zone
!           to the available node.
               zone = num_sent + 1
!            write(*,*) 'myid=', myid, ' before send zone=', zone,
!     1           ' nsurfl(zone)=', nsurfl(zone)
               call r_surf_send_job(sender, zone)
!            write(*,*) 'myid=', myid, ' sent zone=', zone
               num_sent = num_sent + 1
            else
               call surf_send_end_signal(sender)
            endif
 402     continue
!
!
!     End of master part.
!
!      close(10)             
!
 800  format(' rpre = ',e14.7,', zpre = ',e14.7)
 801  format(' xnu = ',e14.7)
 802  format(' phi = ',e14.7,',  wmu = ',e14.7)
 803  format(' jgpsp = ',i3,',  jgplc = ',i2)
 804  format(' jgpmu = ',i2)
 805  format(' jgpsp = ',i3,', ew = ',e14.7)
 806  format('     F = ',e14.7,', xnu = ',e14.7)
 807  format('     F = ',e14.7,', Delta (F) = ',e14.7)
!
!
!
!
!     slave nodes part.
      else if(myid.ne.master) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!
!     vertical zones part
!
      call z_surf_bcast
!
!     if there are more nodes than work skip this node.
      if(myid.gt.nz) goto 990
!
!     as long as the node doesn't recieve the end_signal, it
!     will keep tracking photons for the zones.
!      write(*,*) 'myid=', myid, ' before recv zone=', zone
 991  call z_surf_recv_job(zone)
      if(zone.eq.end_signal) goto 990
      call z_surf_calc(zone)
!      write(*,*) 'myid=', myid, ' js=', zone, ' after z_surf_calc'
      call z_surf_send_result(zone)
!      write(*,*) 'myid=', myid, ' zone=', zone, ' result was sent'
      goto 991
 990  continue
!      write(*,*) 'myid=', myid, ' z recieved end signal'
!      stop
!
!
!     radial zones part.
!
      call r_surf_bcast
!
!     if there are more nodes than work skip this node.
      if(myid.gt.nr) goto 995
!
!     as long as the node doesn't recieve the end_signal, it
!     will keep tracking photons for the zones.
!      write(*,*) 'myid=', myid, ' before recv zone=', zone
 996  call r_surf_recv_job(zone)
!      write(*,*) 'myid=', myid, ' zone=', zone, ' recved job', 
!     1     ' nsurfl(zone)=', nsurfl(zone)
      if(zone.eq.end_signal) goto 995
!      write(*,*) 'myid=', myid, ' before r_surf_calc, ncycle=', ncycle
      call r_surf_calc(zone)
!      if( (myid.eq.2).or.(myid.eq.3) ) write(*,*) 'myid=', myid, 
!     1     ' after r_surf_calc, zone=', zone
!      write(*,*) 'myid=', myid, ' ks=', zone, ' after r_surf_calc'
      call r_surf_send_result(zone)
!      if(myid.eq.2) write(*,*) 'myid=', myid, ' zone=', 
!     1     zone, ' result was sent'
      goto 996
 995  continue
!      write(*,*) 'myid=', myid, ' r recieved end signal'
!
!
!
!     end of slave part
      endif
!
!
!
!     end of imcsurf2d
!      write(*,*) 'myid=', myid, ' end imsurf2d ecens=', ecens(1,1)
 900  return
      end
!
!
!
!     This subroutine calculates and tracks the photons produced 
!     at the inner and outer surface boundaries of the zones.
!     J. Finke 4 April 2005
      subroutine z_surf_calc(js)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer js
!
      integer nxs
      integer npsurf
!
      double precision fibran, t_average, x1, x2
      double precision esurf
      double precision tpl, F(nphomax)
!
!
      call initialize_zrand(js)
!      write(*,*) 'js=',js,' rseed=',rseed
!      write(*,*) 'js=',js,' rand=',fibran()
!      write(*,*) 'js=',js,' rand_switch=',rand_switch
!
!      rseed = 8535
!         added time greater than t0 to make sure the external radiation begins after a certain time ! Xuhui 1/28/10
              if (tbbi(js, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))&
     &            call file_sp(i_fname(js, ti))

              if (nsurfi(js).le.0) goto 201
              npsurfz(js) = 0
              do 200 nxs = 1, nsurfi(js)
!                 write(*,*) 'js=', js, ' nxs=', nxs
                 jph = js
                 npsurfz(js) = npsurfz(js) + 1
!
!                 wmu = 2.d0*drand(0) - 1.d0
                 wmu = 2.d0*fibran() - 1.d0
                 if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                 if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
!
!                 phi = -1.1d1/7.d0 + 2.2d1/7.d0*drand(0)
                 phi = -1.1d1/7.d0 + 2.2d1/7.d0*fibran()
                 if (phi.lt.-1.5707963d0) phi = -1.57079063d0
                 if (phi.gt.1.5707963d0) phi = 1.5707963d0
                 rpre = rmin
                 if (jph.eq.1) then          
                      zpre = z(1)*fibran()
                 else
                      zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                 endif
                 ew = ewsurfi(jph)
                 esurf = esurf+ew
!                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
!
                 if (tbbi(jph,ti).gt.0.d0) then
                    tpl = tbbi(jph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
!
!                 F(jgpsp) = F(jgpsp) + ew
!                
                 kph = 1
                 emitype = 10
!
                 call imctrk2d(-1)
!
 200          continue
 201          continue
!
              if (tbbo(js, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))&
     &           call file_sp(o_fname(js, ti))
              do 250 nxs = 1, nsurfo(js)
                  jph = js
                  npsurfz(js) = npsurfz(js) + 1
!                  wmu = 2.d0*drand(0) - 1.d0
                  wmu = 2.d0*fibran() - 1.d0
                  if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                  if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
!
                  rpre = r(nr)
                  if (jph.eq.1) then          
                      zpre = z(1)*fibran()
                  else
                      zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                  endif
!                  x1 = drand(0)
!                  x2 = drand(0)
                  x1 = fibran()
                  x2 = fibran()
                  if (x1.lt.0.5) then
                     phi = 1.1d1/7.d0 + (1.1d1/7.d0)*x2
                     if (phi.lt.1.57079638d0) phi = 1.57079638d0
                  else 
                     phi = -1.1d1/7.d0 - (1.1d1/7.d0)*x2
                     if (phi.gt.-1.57079638d0) phi = -1.57079638d0
                  endif 
                  ew = ewsurfo(jph)  
                  esurf = esurf+ew
!                  dcen = drand(0)*c_light*dt(1)
                  dcen = fibran()*c_light*dt(1)
!
                  if (tbbo(jph,ti).gt.0.d0) then
                     tpl = tbbo(jph,ti)
                     call planck(tpl)
                  else
                     call file_sample()
                  endif
!
!                  F(jgpsp) = F(jgpsp) + ew
!                 
                  kph = nr
                  emitype = 10
!
                  call imctrk2d(-1)
!
 250           continue
!
!
!
!             write(*,*) 'end js=',js,' rand=',fibran()
               return
               end
!
!
!
!     This subroutine calculates and tracks the photons produced 
!     at the upper and lower surface boundaries of the zones.
!     J. Finke 4 April 2005
      subroutine r_surf_calc(ks)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer ks
!

      integer npsurf, nxs
!
      double precision psi
      double precision fibran, t_average, x1, x2
      double precision esurf
      double precision tpl, F(nphomax)

!
!

!
!
      call initialize_rrand(ks)

             if (tbbu(ks, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti) )&
     &             call file_sp(u_fname(ks, ti))

             npsurfr(ks) = 0
             do 330 nxs = 1, nsurfu(ks)
                 kph = ks
                 npsurfr(ks) = npsurfr(ks) + 1
!                 
!                 wmu = -drand(0)
!                 phi = 4.4d1/7.d0*drand(0)
!                 
                 wmu = -fibran()
                 if (wmu.gt.0.9999999999d0) wmu = 0.9999999999d0
                 if (wmu.lt.-0.9999999999d0) wmu = -0.9999999999d0
!
                 phi = 2.d0*pi*fibran()
!                 psi = drand(0)
                 psi = fibran()
                 Zpre= z(nz)    
                 if (kph.eq.1) then
                       rpre = dsqrt((rmin**2) + psi*(r(kph)**2&
     &                                             - rmin**2))
                 else
                       rpre = dsqrt((r(kph-1)**2) + psi*(r(kph)**2 &
     &                                               - r(kph-1)**2))
                 endif            
                 ew = ewsurfu(kph)
                 esurf = esurf+ew
!
!                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
!
                 if (tbbu(kph,ti).gt.0.d0) then
                    tpl = tbbu(kph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
!
!                 F(jgpsp) = F(jgpsp) + ew
!
                 jph = nz
                 emitype = 10
!
                 call imctrk2d(-1)
!
 330          continue
!
!              write(4,*) 'Lower boundary, zone ',ks,': ',nsurfl(ks),
!     1                   ' photons'
!              write(*,*) 'Lower boundary, zone ',ks,': ',nsurfl(ks),
!     1                   ' photons'
!
              if (tbbl(ks, ti).lt.0.d0 .and. time+0.5d0*dt(1).ge.t0(ti))&
     &                call file_sp(l_fname(ks, ti))

              do 360 nxs = 1, nsurfl(ks)
!                 write(*,*) 'nxs=', nxs
                 kph = ks
                 npsurfr(ks) = npsurfr(ks) + 1
!
!                  wmu = drand(0)
                 if(star_switch.eq.1) then
                    wmu = 9.9999999d-1
                 else
!                   all the beamed external radiation are in the up direction
!                    wmu = fibran()
                     wmu = 9.9999999d-1 ! Xuhui 2/2/10
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                 endif
                  if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                  if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
!
!
!                 wv1 = drand(0)
!                 wv2 = drand(0)
!                 wmu = dmax1(wv1, wv2)
!
!                 phi = (4.4d1/7.d0)*drand(0)
                 phi = 2.d0*pi*fibran()  
                 psi = fibran()
                 Zpre = zmin
                 if (kph.eq.1) then
                       rpre = dsqrt((rmin**2)+ psi*(r(kph)**2&
     &                                            - rmin**2))
                 else
                       rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2&
     &                                              - r(kph-1)**2))
                 endif             
                 ew = ewsurfl(kph)
                 esurf = esurf+ew
!
                 if (tbbl(kph,ti).gt.0.d-20) then
                    tpl = tbbl(kph,ti)
                    call planck(tpl)
                 else
                    call file_sample()
                 endif
!                 dcen = drand(0)*c_light*dt(1)
                 dcen = fibran()*c_light*dt(1)
!
!                 F(jgpsp) = F(jgpsp) + ew
!
                 jph = 1
                 emitype = 10
!
                 call imctrk2d(-1)
!              
 360       continue
!
!
!
!
      return
      end
!
!
!
!
!     Read external spectrum file;
!     calculate flux normalizations
!     and probability distributions
!
!
      subroutine file_sp(fname)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      character *30 fname
!
      integer,parameter :: ndisk = 10
       
!
      integer i, n_input, j, ijump
      double precision mu_disk(ndisk), dmu(ndisk)
      double precision r_d(ndisk),&
     &  Ltot_disk, Ftot_blr, Ftot_ir, Ftot_blr_norm, Ftot_ir_norm
      double precision Id_file(nfmax,ndisk), F_como(nfmax,ndisk),&
     &  L_disk(nfmax), F_blr(nfmax), F_ir(nfmax)
!
      double precision Isum
      double precision beta, dopp(ndisk)
      character *10 junk
!
!

!
      n_input = 25

      open(n_input, file=fname, status='unknown')
      i = 0
      read(unit=n_input, err=120, fmt=*)junk
 100  i = i + 1
      read(unit=n_input, err=120, fmt=*) &
     & E_file(i), L_disk(i), F_blr(i), F_ir(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     calculation of radiative intensity from the accretion disk.
!     followed the calculation in Ghisellini & Tavecchio 2009.
!
!     dopp is the doppler factor of the photons travel from the
!     accretion disk to the jet.
      do j=1,ndisk
        r_d(j)=R_disk*(j-0.5)/ndisk
      enddo
      beta = sqrt(1.d0-1.d0/g_bulk**2)
      mu_disk(:) = d_jet/sqrt(d_jet**2+r_d(:)**2)
      do j=1,ndisk
         dmu(j) = d_jet/sqrt(d_jet**2+(R_disk*(j-1)/ndisk)**2)-&
     &   d_jet/sqrt(d_jet**2+(R_disk*j/ndisk)**2)
      enddo
      dopp(:) = g_bulk*(1.d0-beta*mu_disk(:))
      Id_file(i,:)=L_disk(i)/(4*pi**2*R_disk**2*mu_disk(:))
      F_como(i,:)=2*pi*dopp(:)*Id_file(i,:)*dmu(:)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if ((E_file(i).gt.0.d0).and.(i.le.(nfmax-1))) goto 100
 120  close(n_input)
      nfile = i - 1

!     The following procedure has assumed E_file to be evenly spaced
!     in logarithm. It moves the new spectrum to energy grids of its own
!      F_file(:)=1.d-29
!      do j=1,ndisk
!       if(dopp(j).ge.1.d0)then
!         do ijump=0,nfmax-1
!            if(E_file(ijump+1).ge.E_file(1)*dopp(j))exit
!         enddo
!         do i=1,nfmax-ijump
!            F_file(i+ijump)=F_file(i+ijump)+F_como(i,j)
!         enddo
!       else
!         do ijump=0,nfmax-1
!            if(E_file(nfmax-ijump).le.E_file(nfmax)*dopp(j))exit
!         enddo
!         do i=1,nfmax-ijump
!            F_file(i)=F_file(i)+F_como(i+ijump,j)
!         enddo
!       endif
!      enddo

!     Ltot_disk is the total power of the disk radiation.
!     E_file is the median of every energy grid, so sqrt() is divided for
!     the actual power.
      Ltot_disk = 0.d0
      Ftot_blr = 0.d0
      Ftot_ir = 0.d0
      do i=1,nfmax-1
        Ltot_disk = Ltot_disk+L_disk(i)*(E_file(i+1)-E_file(i))
        Ftot_blr = Ftot_blr+F_blr(i)*(E_file(i+1)-E_file(i))
        Ftot_ir = Ftot_ir+F_ir(i)*(E_file(i+1)-E_file(i))
      enddo
      Ltot_disk = Ltot_disk/sqrt(E_file(2)/E_file(1))
!     Ftot_blr is the total power of the BLR radiation of the input file.
!     Ftot_ir is the total power of the torus radiation of the input file.
      Ftot_blr = Ftot_blr/sqrt(E_file(2)/E_file(1))
      Ftot_ir = Ftot_ir/sqrt(E_file(2)/E_file(1))
!     Ftot_blr_norm and Ftot_ir_norm are the total power the BLR and torus 
!     radiation should have according
!     to Ghisellini & Madau 1996 (timed c for flux)
      Ftot_blr_norm = 17.d0/48.d0/pi*g_bulk**2*fr_blr*Ltot_disk/R_blr**2
      Ftot_ir_norm =  1.d0/4.d0/pi*g_bulk**2*fr_ir*Ltot_disk/R_ir**2
!     F_blr in the end is the normalized BLR and torus flux.
      F_blr(:) = F_blr(:)/Ftot_blr*Ftot_blr_norm
      F_ir(:) = F_ir(:)/Ftot_ir*Ftot_ir_norm



!************** Here it decides which of the above calculation counts*******
!      F_file(:)=F_file(:)+F_blr(:)
      F_file(:)=F_blr(:)+F_ir(:)

!     flare by the change of external photons
!     the duration of the change is obtained by assuming the blob crossing
!     a photon rich/depleted region. The sigma_z is translated and rounded
!     in units of dz
!     inj_switch 125 has theta and external flares at the same time
      if((inj_switch.eq.2.or.inj_switch.eq.125).and.(time-inj_t).gt.0.d0.and.&
     &  (time-inj_t).lt.int(sigma_z/dz+0.5d0)*dz/inj_v)then
         F_file(:) = F_file(:)*flare_amp
         write(*,*)'external flux increased.'
      endif   
!******************************c
      if (nfile.lt.2) then
         write(*,*) 'Error in input spectrum file:'
         write(*,*) 'Less than 2 lines of input read!'
         write(*,*)
         stop
      endif
!
!      do 130 i = 1, nfile
!         write(*,*) i, E_file(i), F_file(i)
! 130  continue
!
      Isum = 0.d0
      do 150 i = 1,nfile-1
         alpha(i) = log(F_file(i+1)/F_file(i))&
     &             /log(E_file(i+1)/E_file(i))
         a1(i) = alpha(i) + 1.d0
         if (a1(i).gt.2.d1) a1(i) = 2.d1
         if (a1(i).lt.-2.d1) a1(i) = -2.d1
         if (dabs(a1(i)).lt.1.d-3) then
            I_file(i) = F_file(i)*E_file(i)&
     &                 *log(E_file(i+1)/E_file(i))
         else
            I_file(i) = F_file(i)*E_file(i)&
     &                *((E_file(i+1)/E_file(i))**a1(i)&
     &                  - 1.d0)/a1(i)
         endif
         Isum = Isum + I_file(i)
         P_file(i) = Isum
 150  continue
!
      do 200 i = 1, nfile-1
         P_file(i) = P_file(i)/Isum
 200  continue

      int_file = Isum
!
      return
      end
!
!
!
!
!     Subroutine to sample a photon energy
!     from a user-specified input spectrum
!
!
      subroutine file_sample
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer i
      integer jbot, jtop, jmid, m, n
!
      double precision hubot, hutop, Isum
      double precision fibran, x1, x2
!
!      write(*,*) 'in file_sample'
!      write(*,*) 'xnu=', xnu
!      write(*,*) 'nfile=', nfile
!
!      x1 = drand(0)
      x1 = fibran()
      do 210 i = 1, nfile-1
 210  if (P_file(i).gt.x1) goto 220
 220  continue
!
!      x2 = drand(0)
      x2 = fibran()
      xnu = E_file(i)*((a1(i)*I_file(i)*x2&
     &     /(F_file(i)*E_file(i)) + 1.d0)**(1./a1(i)))
!
!
!
!     Determine photon group of new photon energy
!     in photon energy bin structure for spectrum
!
      jbot = 1
      jtop = nphtotal + 1
      hubot = 1.000001d0*hu(1)
      hutop =.999999d0*hu(jtop)
      if( xnu .ge. hutop ) then
!         xnu = hutop
!         jgp = numgps
         jgpsp = 0
         go to 300
      endif
      if( xnu .le. hubot ) then
!         xnu = hubot
!         jgp = 1
         jgpsp = 0
         go to 300
      endif
!
  140 continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 160
      if(xnu .eq. hu(jmid)) go to 160
      if(xnu .lt. hu(jmid)) then
         jtop = jmid
      else
         jbot = jmid
      endif
      go to 140
!
  160 continue
      jgpsp = jmid
!
  300 continue
!
!
!     Determine photon group of new photon energy
!   in photon energy bin structure for light curves
!
      jgplc = 0
      do 400 m = 1, nph_lc
         if ((xnu.gt.Elcmin(m)).and.&
     &       (xnu.le.Elcmax(m))) then
            jgplc = m
            goto 600
         endif
  400 continue
!
  600 continue
!
!
!    Determine angular bin of photon
!
      do 700 n = 1, nmu
         if (wmu.le.mu(n)) then
            jgpmu = n
            goto 900
         endif
  700 continue
      jgpmu = nmu
!
  900 continue      
!
      return
      end

!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:25:40 EDT 2006
! version: 2
! Name: J. Finke
! Changed call to 'initialize_rand' to 'initialize_zrand' and 'initialize_rrand'. 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Sun Aug  6 18:51:03 EDT 2006
! version: 3
! Name: J. Finke
! Added possibility of using upper surface to represent 
! a star. Upper surface photons can come in 
! parallel to each other, representing the star, instead 
! of in random directions.     
