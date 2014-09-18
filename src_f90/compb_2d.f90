!        actual Compton scattering of an individual photon.
!        for details refer to Podznyakov et al. (1983)
        subroutine compb2d(i_gam)
        implicit none
        include 'general.pa'
        include 'commonblock.f90'
!
        
        integer i_gam, m, n, icoms
        integer isurfmu
        integer jbot, jtop, bdotr, ncoll, jmid
       
        double precision wa, wb, swa, cazes, omege,  fuzz, cazs, caz    
        double precision phis, phie, rhos, alphas, theta, omeg
        double precision znue, znue3, betz, gamz, xxx, sz, &
     &                   games, gams, phat, omegs, omeges
        double precision costheta, cosphis, cos_alphas, tl,&
     &                   tr, znues, znus, xknot, sintheta, cosphi
        double precision thetas, cosrhos, cosalphas, ucaz, ucazs,&
     &                   uca, coscaz, znu
        double precision xnus, wmus, hubot, hutop
        double precision tt
        double precision emasskev, gamm, betb
        double precision cosdphi, dphi
        double precision fibran
!
!
!
!
        data fuzz / 1.d-10 /
!
        icoms = 1
        emasskev = 5.11d2      
        znu = xnu / emasskev
        tt = tea(jph,kph) / emasskev
!
!        
 100  continue        
!
!
!     Draw new particle energy and velocity
!      from electron population according
!       to FP solution (MB, 13/May/1999)
!
!      write(*,*) 'xnu = ',xnu
!      write(*,*) 'wmu = ',wmu
!      write(*,*) 'phi = ',phi
!
      call nth2d(jph, kph, gamm, betb, i_gam)
!
!      write(*,*) 'gamma = ',gamm
!      write(*,*) 'beta = ', betb
!
!
!     betb = sample of speed from isotropic relativistic
!     Maxwellian + non-thermal power-law distribution.
!     This prospective target may end up rejected.  next select
!     omeg = lab cosine of angle between electron and incident photon
!
!         omeg = 2.d0*drand(0) - 1.d0
         omeg = 2.d0*fibran() - 1.d0
         if (omeg.gt.9.9999999d-1) omeg = 9.9999999d-1
         if (omeg.lt.-9.9999999d-1) omeg = -9.9999999d-1
!         tl = drand(0)
         tl = fibran()
         tr = 0.5d0*(1.0d0 - betb*omeg)
         if(tl .gt. tr) omeg = -omeg

         if (omeg.gt.9.9999999d-1) omeg = 9.9999999d-1
         if (omeg.lt.-9.9999999d-1) omeg = -9.9999999d-1
!
!     znue = doppler shifted photon frequency.
!     accept target electron with prob. = xknot = total K-N/Thomson.
!     for small znue, use series to prevent bad roundoff.
!     (for znue = .01, error .lt. (544/7)*znue**5 = e-08)
!
      znue = (1.-betb*omeg)*znu*gamm
      if  (znue.lt.1.d-10) goto 100
      if(znue .le. 1.d-2) then
        xknot = 1.0d0-znue&
     &         *(2.0d0-znue*(5.2d0-znue*(13.3d0-1.144d3*znue/3.5d1)))
      else
        znue3 = znue*znue*znue
        betz   = 1.0d0 + 2.0d0*znue
        gamz   = znue*(znue-2.0d0) - 2.0d0
        xxx = 4.0d0*znue + 2.0d0*znue3*(1.0d0+znue)/betz**2 &
     &      + gamz*dlog(betz)
        xknot = 3.75d-1*xxx/znue3
      endif
!
!      write(*,*) 'xknot = ',xknot
!      write(*,*) 'znue = ',znue
!        
!      if(drand(0) .gt. xknot) go to 100
       if(fibran() .gt. xknot) go to 100
!     ncoll = ncoll + 1
!       
      betz = 1.0d0 + 2.0d0*znue
!
 200  continue
 202  sz = (1.0d0 + 2.0d0*znue*fibran())/betz
!     sz = (1.0 + 2.0*znue*drand(0))/betz
      games = 1.0d0 + (1.0d0-1.0d0/sz)/znue
      if ((1.d0 - games**2).lt.0.) goto 202
      tr = games**2 - 1.0d0 + sz + 1.0d0/sz
      phat = betz + 1.d0/betz
!     if(drand(0)*phat .gt. tr) go to 200
      if(fibran()*phat .gt. tr) go to 200
      znues = znue*sz
      
!      write(*,*) 'znues = ',znues
!
  210 continue
!      wa = drand(0)
!      wb = 2.*drand(0) - 1.
      wa = fibran()
      wb = 2.d0*fibran() - 1.d0
      swa = wa*wa+wb*wb
      if((swa.ge.1.0d0) .or. (swa.le.1.d-20)) go to 210
!
      cazes = (wa*wa-wb*wb)/swa
      
!      write(*,*) 'cazes = ',cazes
!        
!     omege = cosine of angle between photon and electron in rest frame.
!     omeges = cosine between elec. direction and scattered photon.
!
      omege = (omeg-betb)/(1.-betb*omeg)
      if (omege.gt.9.9999999d-1) omege = 9.9999999d-1
      if (omege.lt.-9.9999999d-1) omege = -9.9999999d-1
      omeges = games*omege + &
     &         cazes*dsqrt((1.-omege**2+fuzz)*(1.-games**2))
      if (omeges.gt.9.9999999d-1) omeges = 9.9999999d-1
      if (omeges.lt.-9.9999999d-1) omeges = -9.9999999d-1
!
!      write(*,*) 'omeges = ',omeges
!      write(*,*) 'omege = ',omege
!
!     xform back to lab frame.
!
!     omegs = cosine between electron and scattered photon in lab frame.
!     znus = scattered photon frequency in lab.
!     gams = cosine of scattering angle in lab.
!
      omegs = (omeges+betb)/(1.0d0+omeges*betb)
      if (omegs.gt.9.9999999d-1) omegs = 9.9999999d-1
      if (omegs.lt.-9.9999999d-1) omegs = -9.9999999d-1
      znus = (1.d0 + betb*omeges)*gamm*znues
      gams = 1.d0 - (znue - znues)/(znu*znus)
      if (gams.gt.9.9999999d-1) gams = 9.9999999d-1
      if (gams.lt.-9.9999999d-1) gams = -9.9999999d-1
!
!      write(*,*) 'omegs = ',omegs
!      write(*,*) 'znus = ',znus
!      write(*,*) 'gams = ',gams
!
!     calculate lab scattered direction vector relative to radius.
!     cazs = cos. of azimuthal scattering angle (uniform from symmetry.)
!     wmus = lab cosine between scattered photon and radius.
!
 220  continue
!      wa = drand(0)
!      wb = 2.*drand(0) - 1.d0
      wa = fibran()
      wb = 2.d0*fibran() - 1.d0
      swa = wa*wa + wb*wb
      if((swa.ge.1.0d0) .or. (swa.le.1.d-20)) go to 220
      cazs = (wa*wa-wb*wb)/swa
      if (cazs.gt.9.9999999d-1) cazs = 9.9999999d-1
      if (cazs.lt.-9.9999999d-1) cazs = -9.9999999d-1
      wmus = wmu*gams + cazs*dsqrt((1.0d0-gams**2)&
     &      *(1.0d0-wmu**2 + fuzz))
      if (wmus.gt.9.9999999d-1) wmus = 9.9999999d-1
      if (wmus.lt.-9.9999999d-1) wmus = -9.9999999d-1
      xnus = znus*emasskev
!
!      write(*,*) 'cazs = ',cazs
!      write(*,*) 'wmus = ',wmus
!      write(*,*) 'xnus = ',xnus
!
!  Introducing the two dimensional scenario in compton scattering
!  Now using two sets of axes 'z' and 'r'.
!  phie =  Azimuthal scattering angle of the electron - randomly distributed.
!  phis =  Azimuthal scattering angle of the scattered photon with the 'r' axis.
!  phi =  
!  theta =  Angle between the 'z' axis and the electron.
!  rhos = Angle between the scattered photon and 'r' axis.
!
!
!      phie = 4.4d1*drand(0)/7.d0
!
!      caz = (omegs- gams *omeg)/ dsqrt((1.-omegs**2)*(1.-gams**2 ))
!      if (caz.gt.0.99999999d0) caz = 0.99999999d0
!      if (caz.lt.-0.99999999d0) caz = -0.99999999d0
!      ucaz  = acos(caz)
!      ucazs = acos(cazs) 
!      uca   = ucaz + ucazs
!      coscaz = cos(uca)
!      
!      write(*,*) 'caz = ',caz
!      write(*,*) 'uca = ',uca
!      write(*,*) 'coscaz = ',coscaz
!      
!      costheta = omeg*wmu + (dsqrt((1.-omeg**2)*(1.-wmu**2))*coscaz)
!      if (costheta.gt.0.99999999d0) costheta = 0.99999999d0
!      if (costheta.lt.-0.99999999d0) costheta = -0.99999999d0
!      theta = acos(costheta)
!      sintheta = sin(theta)
!      if (dabs(sintheta).lt.1.d-10) sintheta = 1.d-10
!      cosphi = (omegs - wmus*costheta)/(sintheta*dsqrt(1.-wmus**2))
!      if (cosphi.gt.0.99999999d0) cosphi = 0.99999999d0
!      if (cosphi.lt.-0.99999999d0) cosphi = -0.99999999d0
!      
!      write(*,*) 'costheta = ',costheta
!      write(*,*) 'sintheta = ',sintheta
!      write(*,*) 'cosphi = ',cosphi
!
!
!      phi = acos(cosphi)
!      phis = phie-phi
!      thetas = acos(wmus)
!
!
!      alphas = ( 2.2d1/1.4d1 - thetas)
!      cosphis= cos(phis)
!      cos_alphas= cos(alphas)
!      cosrhos = cosphis*cos_alphas
!      if (cosrhos.gt.0.99999999d0) cosrhos = 0.99999999d0
!      if (cosrhos.lt.-0.99999999d0) cosrhos = -0.99999999d0
!      rhos = acos(cosrhos)
!
       cosdphi = (gams - wmu*wmus)&
     &          /dsqrt((1.0d0 - wmu**2)*(1.0d0 - wmus**2))
       if (cosdphi.gt.9.9999999d-1) cosdphi = 9.9999999d-1
       if (cosdphi.lt.-9.9999999d-1) cosdphi = -9.9999999d-1
       dphi = acos(cosdphi)
       phis = phi + dphi
!
!       write(*,*) 'cosdphi = ',cosdphi
!       write(*,*) 'dphi = ',dphi
!       write(*,*) 'phis = ',phis
!
!
  250 continue
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
      if( xnus .ge. hutop ) then
         jgpsp = 0
         go to 300
      endif
      if( xnus .le. hubot ) then
         jgpsp = 0
         go to 300
      endif
!
  140 continue
      jmid = (jbot+jtop)/2
      if( jmid .eq. jbot ) go to 160
      if(xnus .eq. hu(jmid)) go to 160
      if(xnus .lt. hu(jmid)) then
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
         if ((xnus.gt.Elcmin(m)).and.&
     &       (xnus.le.Elcmax(m))) then
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
         if (wmus.le.mu(n)) then
            jgpmu = n
            goto 900
         endif
  700 continue
      jgpmu = nmu
!
!
!
 900  continue
      ew = ew*xnus/xnu
      xnu = xnus
      wmu = wmus
      phi = phis

!  
!      write(*,*) 'xnus = ',xnu
!      write(*,*) 'wmus = ',wmu
!      write(*,*) 'phis = ',phis
!
      return
      end
     
      
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:24:02 EDT 2006
! version: 2
! Name: J. Finke
! Fixed bug in common block 'random'.   
