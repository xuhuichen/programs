!      write the photon to the output as the photon hits the volume surface and becomes part 
!      of the emission that can be observed.
!      comment by Xuhui Chen 2013-07-15
      subroutine imcleak
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
!     NAME OF ROUTINE...imcleak
!
!     AUTHOR OF ROUTINE...E. H. Canfield July 27, 1983
!
!     PURPOSE...
!
!        Called by imctrk when particle hits inner or outer surface, or
!        by imcfield if hydro motion sweeps either surface past particle.
!
!        default mode:
!        Particles at the outer surface escape the system, and all leakage
!        and pinex/pipe data are updated.  Particles at the inner surface
!        are kept for tracking across the void (except for slabs).
!
!        Options allow specular reflection at surfaces, and killing particles
!        entering an interior void.
!
!        For reflection at outer surface, no energy escapes, but leakage
!        and pinex/pipe arrays are updated to show internal spectral data.
!
!     GLOSSARY OF IMPORTANT VARIABLES...
!
!        izbnd = 1  particle is at xq(1)
!              = jm particle is at xq(jm).
!        ibrad = 1  xq(1 ) is a reflecting surface.
!              = 2  xq(jm) is a reflecting surface.
!              = 3  both are reflecting surfaces.
!        lkin  = 1  particle is to be killed on entering interior void.
!        idead = 1  set if particle leaked from system.
!              = 3  set if particle was reflected.
!              = 4  set if particle is to be tracked across void.
!
!     LIST CHANGES HERE...
!
!ccccc
! 09/02/88 spectral "enhancement" code for dermer.
!      accumulate data for <xnu> and <xnu**2>
! 11/23/04 implicit none added at begining.  J.D. Finke
!ccccc
!/c
!
!
      integer ef_switch, n_in, n_out, i_0, i_1,&
     &        jbot, jtop, jmid
      integer n, m
!     To calculate the spectrum incident on the upper and lower
!     boundaries, such as in photon bubble simulations, set
!     spec_switch to 1.  Otherwise, to calculate the spectrum
!     that leaves the region, set it to 0.
!     J. Finke 7/22/06
!
      double precision rnum, ewref, ewnew
      double precision hubot, hutop
      double precision f
      double precision xnu_sv
      double precision fibran
!


!
!
      xnu_sv = xnu
!
      if (kph.eq.0) then
!
!        Particle is at inner r-boundary:
!           If rmin > 0, absorb photon,
!         else continue on the other side
!
         if (rmin.gt.1.d-10) then
            idead = 1
            erlki(jph) = erlki(jph) + ew
!            incounter = incounter + 1
            goto 900
         else
            idead = 0
            phi = 1.d-6
            kph = 1
            goto 900
         endif
      endif
!
      if(jph .gt. 0) go to 100
!
!     particle is at lower surface:
!         Test for reflection.
!
      if (tbbl(kph,ti).gt.0.d0) then
         Ed_in(kph) = Ed_in(kph) + ew
         erlkl(kph) = erlkl(kph) + ew
         if (kph.eq.2) then
!         write(*,*) 'ew=',ew,' k=',k,' Ed_in=',Ed_in(k)
!         stop
         endif
      endif
!
!       Compton reflection at lower boundary
!       (cr_sent = 1 or cr_sent = 3)
!
      if ((cr_sent.eq.1).or.(cr_sent.eq.3).or.(cr_sent.eq.4)) then
!         dE_in = dE_in + ew
!
!      write(*,*) 'Compton reflection at lower boundary.'
!      write(*,*) 'xnu = ',xnu
!
         if ((tbbl(kph,ti).le.0.d0).or.(cr_sent.eq.4)) then
!            write(*,*) 'b4 refl wmu=', wmu
            wmu = -wmu
            if ( (jgpsp.gt.0).and.(spec_switch.eq.1) )&
     &         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
!            write(*,*) 'j=', j, ' k=', k, ' refl.'
!            write(*,*) ' wmu=', wmu
!            stop
            goto 60
         endif
!
         n_in = 1
 10      if ((E_ref(n_in).lt.xnu).and.(n_in.lt.n_ref)) then
            n_in = n_in + 1
            goto 10
         endif
!
!         write(*,*) 'n_in = ',n_in
!
         n_out = 1
!         rnum = drand(0)
         rnum = fibran()
 20      if ((P_ref(n_out, n_in).lt.rnum).and.(n_out.lt.n_ref)) then
            n_out = n_out + 1
            goto 20
         endif
!
!         write(*,*) 'n_out = ',n_out
!
         if (n_out.gt.1) then
!            xnu = E_ref(n_out-1) 
!     1          + drand(0)*(E_ref(n_out) - E_ref(n_out-1))
            xnu = E_ref(n_out-1) &
     &          + fibran()*(E_ref(n_out) - E_ref(n_out-1))
         else
            xnu = E_ref(n_out)
         endif
!
!         write(*,*) 'xnu_new = ',xnu
!         write(*,*) 'W_abs = ',W_abs(n_out, n_in)
!
         ewnew = ew*W_abs(n_out, n_in)*xnu/xnu_sv
         Ed_ref(kph) = Ed_ref(kph) + ewnew
         wmu = dabs(wmu)
         ew = ewnew
!
!
 60      continue
         call get_bin
!
         idead = 0
         jph = 1
         goto 900
!
      else
!     no reflection, photon exits.
!     Added by J. Finke, 7 August 2006.
         if (ncycle.gt.0) then
!            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu,
            write(nunit_evt) t_bound, xnu, ew, rpre, zpre, wmu,&
     &                            phi, emitype
            if (jgplc.gt.0) &
     &           edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )&
     &           fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
         idead = 1
         goto 900
      endif
!
 105  format(7(e14.7,1x),i4)
!
!
!         particle is at outer surface.
!    Update all spectral and light curve data
!
  100 continue
!
      if (jph.eq.nz+1) goto 500
!
!       Particle is at outer r-boundary:
!         Test for Compton reflection
!
      erlko(jph) = erlko(jph) + ew
!      outcounter = outcounter + 1
!
      if ((wmu.gt.0.).or.(cr_sent.eq.0).or.&
     &    (cr_sent.eq.1).or.(cr_sent.eq.4)) then
!
!             Write photon 
!       to spectrum and event file
!
         t_bound = time + dt(1) - rad_cp*dcen
         if (ncycle.gt.0) then
!            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu, 
            write(nunit_evt) t_bound, xnu, ew, rpre, zpre, wmu,&
     &                            phi, emitype
            if (jgplc.gt.0) &
     &         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )&
     &         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
         idead = 1
         goto 900
!
!
!
      else
!
!      Compton reflection at outer disk 
!
!         dE_in = dE_in + ew*f_ref
!
         n_in = 1
 110     if ((E_ref(n_in).lt.xnu).and.(n_in.lt.n_ref)) then
            n_in = n_in + 1
            goto 110
         endif
         n_out = 1
!         rnum = drand(0)
         rnum = fibran()
 120     if ((P_ref(n_out, n_in).lt.rnum).and.(n_out.lt.n_ref)) then
            n_out = n_out + 1
            goto 120
         endif
!
         if (n_out.gt.1) then
!            xnu = E_ref(n_out-1) 
!     1          + drand(0)*(E_ref(n_out) - E_ref(n_out-1))
            xnu = E_ref(n_out-1) &
     &          + fibran()*(E_ref(n_out) - E_ref(n_out-1))
         else
            xnu = E_ref(n_out)
         endif
!
!
         ewref = ew*W_abs(n_out, n_in)*xnu/xnu_sv
         ew = ewref
         if (dabs(wmu).gt.1.d-6) then 
            t_bound = time + dt(1) &
     &              - rad_cp*(dcen + zpre/wmu)
            zpre = 0.d0
            f = zpre*dsqrt(1. - wmu**2)/dabs(wmu)
            rpre = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*cos(phi))
         else 
            t_bound = 1.d20
         endif
!         wmu = drand(0)
         wmu = fibran()
!
!
 180     continue
         call get_bin
!
         if (ncycle.gt.0) then
!            write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu,
            write(nunit_evt) t_bound, xnu, ew, rpre, zpre, wmu,&
     &                            phi, emitype
            if (jgplc.gt.0) &
     &         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
            if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )&
     &         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
         endif
!
         idead = 1
         goto 900
!
!
      endif
!
!
!      Photon is at upper z surface.
!       Write photon to event file
!
!
 500  continue

!     Reflection at upper boundary added
!     J. Finke, 28 July 2005
!      if(upper_sent.eq.1) then
!         wmu = -wmu
!         if ( (jgpsp.gt.0).and.(spec_switch.eq.1) )
!     1        fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
!         do 270 n = 1, nmu
!            if(wmu.le.mu(n)) then
!               jgpmu = n
!               goto 180
!            endif
! 270     continue
!         jgpmu = nmu
!      endif
!      Xuhui deleted this because the goto command is getting warning from the
!      compiler

      erlku(kph) = erlku(kph) + ew
      t_bound = time + dt(1) - rad_cp*dcen
      if (ncycle.gt.0.and.wmu.lt.0.98) then
!          write(nunit_evt, 105) t_bound, xnu, ew, rpre, zpre, wmu,
          write(nunit_evt) t_bound, xnu, ew, rpre, zpre, wmu,&
     &                          phi, emitype
         ! wmu<0.98 to remove upper boudary record of large amount of external
         ! radiation coming out there. Xuhui 5/28/11
          if (jgplc.gt.0) &
     &         edout(jgpmu, jgplc) = edout(jgpmu, jgplc) + ew/dt(1)
          if ( (jgpsp.gt.0).and.(spec_switch.eq.0) )&
     &         fout(jgpmu, jgpsp)  = fout(jgpmu, jgpsp) + ew
      endif
!
      idead = 1
!
!
!
  900 continue
!
      return
      end
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
!
!
!     This subroutine determines the spectral, lightcurve, and
!     anglular bin a photon should be in.
!     J. Finke 2 January 2007
      subroutine get_bin
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
      integer ef_switch, n_in, n_out, i_0, i_1,&
     &        jbot, jtop, jmid
      integer m, n
      double precision hubot, hutop
!
!

!
!     Determine photon group of new photon energy
!     in photon energy bin structure for spectrum
!
      jbot = 1
      jtop = nphtotal + 1
      hubot = 1.000001*hu(1)
      hutop =.999999*hu(jtop)
      if( xnu .ge. hutop ) then
         jgpsp = 0
         go to 30
      endif
      if( xnu .le. hubot ) then
         jgpsp = 0
         go to 30
      endif
!
 24   continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 26
      if(xnu .eq. hu(jmid)) go to 26
      if(xnu .lt. hu(jmid)) then
         jtop = jmid
      else
         jbot = jmid
      endif
      go to 24
!     
 26   continue
      jgpsp = jmid
!     
 30   continue
!
!
!     Determine photon group of new photon energy
!   in photon energy bin structure for light curves
!
      jgplc = 0
      do 40 m = 1, nph_lc
         if ((xnu.gt.Elcmin(m)).and.&
     &        (xnu.le.Elcmax(m))) then
            jgplc = m
            goto 60
         endif
 40   continue
!
 60   continue
!
!
!    Determine angular bin of photon
!
      do 70 n = 1, nmu
         if (wmu.le.mu(n)) then
            jgpmu = n
            goto 80
         endif
 70   continue
      jgpmu = nmu
!
 80   continue
!     
!     
!
      return
      end
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun  2 15:14:30 EDT 2006
! version: 1
! Name: J. Finke
! Replaced ran1 with fibran and ran1. Added reflection 
! at lower boundary.      
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun  2 15:22:22 EDT 2006
! version: 2
! Name: J. Finke
! Added fibran as well and ran1 for random 
! number generation. Added total reflection at lower boundary. 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Sun Jun  4 19:37:16 EDT 2006
! version: 3
! Name: J. Finke
! added zseeds and rseeds     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Thu Jul 27 13:07:43 EDT 2006
! version: 4
! Name: J. Finke
! Added spec_switch, so quick looks spectra can be 
! spectra incident on top and bottom boundaries.  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Aug  8 12:12:06 EDT 2006
! version: 5
! Name: J. Finke
! Fixed bug at end, where photons at upper 
! z surface are read to event and spectra 
! files twice. Also changed leak to lower z 
! surface, so that photons there are read to 
! event and spectra files.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jan  2 20:15:39 EST 2007
! version: 5
! Name: J. Finke
! Created subroutine get_bin for determining bins of a 
! photon.        
