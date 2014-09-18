      subroutine planck(tpl)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!     called by genimc, imcsurf, imcvol
!
!    AUTHOR...Gene Canfield
!
!     samples normalized planckian (or wien) spectrum
!
!    GLOSSARY...
!
!     input: tpl(j) = temperature array
!             iwien = 0  select from planck spectrum
!                   = 1  select from   wien spectrum
!             hu(l) = photon group boundary array [keV]
!     output:
!             xnu = photon energy [keV]
!           jgpsp = photon group (spectrum binning)
!           jgplc = light curve photon group
!           jgpmu = angular bin no.
!
!
      double precision tpl
!
      integer iwien
      integer jbot, jtop, jmid, m, n
!
      double precision hubot, hutop
      double precision u4, ap0, ap1, ap2, ap3, rn1
      double precision fibran
!
!
!
!
      iwien = 0
!
!   99 u4 = drand(0)
!      u4 = u4*drand(0)
!      u4 = u4*drand(0)
!      u4 = u4*drand(0)
!
   99 u4 = fibran()
      u4 = u4*fibran()
      u4 = u4*fibran()
      u4 = u4*fibran()
      if (u4.le.1.d-200) goto 99
      ap0 = -log(u4)
      ap1 = 1.d0
      ap2 = 1.d0
      ap3 = 1.d0
      if (iwien .eq. 1) go to 120
!      rn1 = 1.08232d0*drand(0)
      rn1 = 1.08232d0*fibran()
!
  100 continue
      if( rn1 .le. ap1 ) go to 120
      ap2 = ap2 + 1.d0
      ap3 = 1.d0/ap2
      ap1 = ap1 + ap3**4
      go to 100
!
  120 continue
      xnu = ap0*ap3*tpl
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
!      write(10, 800) rpre, zpre
!      write(10, 801) xnu
!      write(10, 802) phi, wmu
!      write(10, 803) jgpsp, jgplc
!      write(10, 804) jgpmu
!
!
      return
!
      end
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:32:12 EDT 2006
! version: 2
! Name: J. Finke
! Changed common block 'random'.     
