!     Tracks the evolution of individula photons for one time step.
!     This subroutine is intensively called by imcfield, imcsurf and imcvol.
!     It calculate three distances: 
!     1. the distance the photon can travel if it moves freely for one time step;
!     2. the distance the photon can travel until it hits a boundary;
!     3. the distance the photon can travel before it gets Compton scattered.
!     The shortest one will be what happens to the photon.
!     The photon also gets absorbed by Synchrotron Self-absorption and pair production processes.
!     If the ew (statistical energy weight) of the photon dropped below wtmin, the photon is killed.
!
!     scat_flag = -1, normal tracking including splitting of photons,
!        for the portion of photon that has a chance to be scattered,
!            For (1) and (2), the photon is not written to the memory, but combined to the other portion,
!            For (3), imctrk(-1) is called again after scattering.
!        for the portion of photon that is not scattered, imctrk(0) is called.
!            For (1), the photon is written to the memory,
!            For (2), imcleak is called if the boudary is a surface.
!     scat_flag = 0, track of unscattered portion of the photons,
!        the photon will not be scattered, but only written to memory when it comes to rest, 
!        or leaked from the surfaces.
!     scat_flag = 1, track and write photons when they are scattered. obsolete
!
!     wmu is the sine of the elevation angle
!     phi is the azimuth angle, with the radial direction as angle 0
      recursive subroutine imctrk2d(scat_flag)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!, ucens, iucens, ducens
      integer scat_flag
!
!
      integer i_gam
      integer i 
!___________________________________________________________
      integer iii,ii,ii2, nscat, iitrk
      integer jcsv,kcsv,jgpspcsv,&
     &        jgplccsv,jgpmucsv,emitypecsv
      double precision ewcsv, xnucsv, wmucsv,rprecsv,&
     &                 zprecsv, phicsv, dcencsv
      integer jtsv,ktsv,jgpsptsv,&
     &        jgplctsv,jgpmutsv,emitypetsv
      double precision ewtsv, xnutsv, wmutsv,rpretsv,&
     &                 zpretsv, phitsv, dcentsv
      
      integer nptrks, knew, jnew, kbnd, jbnd, npkill
      integer lwad, lwai
      integer i_gg, inout, npbnd
      integer isurfmu
      integer ikind
      integer Eta_switch
!
      double precision ewold
!      double precision comac_ar(jmax,kmax), enexc_ar(jmax,kmax)
      double precision Zr, Rr, eoutr(jmax, kmax)
      double precision disp, psq, disbr, dpbsq, Zbnd, rbnd,&
     &                 xqsqleft, trld, trldb, dcol       
      double precision wtmin, wkth, wmue, ewsv2
      double precision Egg_min
      double precision kgg, kgg0, kgg1
      double precision pair_enhance
      double precision theta, Eta, comac, enexc
      double precision sigabs, velfact, facomp
      double precision colmfp, mb_ran, denom, f
      double precision sigcom, sigcomex, sigsc, deleabs, delecomp
      double precision ewnew, xabs, ekill, ekillt, ewpl, sstar
      double precision delpr, exan, wmusv, xnusv, wmustar
      double precision kappa_ic


      double precision fibran, rnew, znew
      real etime, elapse(2), et0
!
!
       idead = 0
       sigabs = 1.d-40
       velfact = 1.d0
       facomp = 1.d0
!       wkth = 1.d-2
       wkth = 1.d-10
!
!      idead = -1 from imcvol if rwlk turned on.
!              =  0 particle returned from rwlk alive.  track it.
!              =  1 killed in rwlk or imcleak.
!              =  2 retire rwlk particle to census.
!              =  3 particle was reflected from surface in imcleak.
!              =  4 particle leaked into void for tracking.
!
!     
         wtmin = wkth*ew
!
!         
!        Sum up squared photon energies
!        and squares for recoil  
!         
!          sum_e = sum_e + 1.957d-3*ew*xnu
!          sum_e2 = sum_e2+ew*((1.957d-3*xnu)**2)
!
!      Arrays to store the IC cross section and fractional energy change for one MC particle.
!      Save computer time.
!       comac_ar = 0.d0
!       enexc_ar = 0.d0

!      This is the first split of MC particles
       nscat = 0 ! the number of scattering happened in the first splitting
!      store all the MC particle informations before it is splitted.
       if(scat_flag.eq.-1)then
         ewtsv = ew/split1 ! first splitting.
       else
         ewtsv = ew
       endif
       dcentsv = dcen
       xnutsv = xnu
       zpretsv = zpre
       rpretsv = rpre
       wmutsv = wmu
       phitsv = phi
       jtsv = jph
       ktsv = kph
       jgpsptsv = jgpsp
       jgplctsv = jgplc
       jgpmutsv = jgpmu
       emitypetsv = emitype
!################################################################################################
       ew = ewtsv

  100  continue
!       sigabs = 1.d-40 ! Xuhui 5/24/09
!         mb_ran = drand(0)
        if(scat_flag .eq. 0)then
          mb_ran = 1.d-10
        else
          mb_ran = fibran()
! Xuhui 5/1/09
         if(rand_switch.eq.2)mb_ran=int(mb_ran*1.d6)/1.d6+1.d-6*fibran()
! added to combine two 7 digits random number into one 13 digits random number.
! the end point of the highest order of fibran() may appear more frequent
! than other points. So we choose to use the second highest order.
         mb_ran = 1.d0-(1.d0-mb_ran)/split1 ! Xuhui risky
        endif
        colmfp = -dlog(mb_ran)
!
! Prepare IC scattering probability. it will not be zeroed when the 
! particle crossed the boundary.
       kappa_ic=0.d0
  110  continue
       sigabs = 1.d-40 ! moved from after 100 to after 110 on 10/23/12
!
       if (ew.lt.1.d-40) goto 900
      if (abs(wmu).gt.1.d0+1.d-15)&
     &    write(*,'("wmu too extreme! wmu=",e24.17)') wmu
      if (wmu.gt.1.d0) wmu = 1.d0
      if (wmu.lt.-1.d0) wmu = -1.d0
!       if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
!       if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
!          
         if (pair_switch.eq.1) then
             pair_enhance = 1.d0 + 2.d0*f_pair(jph, kph)
         else
             pair_enhance = 1.d0
         endif
!  
!        if( scat_flag .eq. -1 )then
!          if(comac_ar(jph,kph).eq.0.d0)then

              call comtot(jph, kph, xnu, 1,comac,enexc)

!              comac_ar(jph,kph) = comac
!              enexc_ar(jph,kph) = enexc
!          else
!              comac = comac_ar(jph,kph)
!              enexc = enexc_ar(jph,kph)
!          endif ! Xuhui 3/4/09
!        elseif( scat_flag .eq. 1 )then
!             call comtot(jph, kph, xnu, 1,comac,enexc)
!        elseif( scat_flag.eq.0 )then
!             comac = 0.
!             enexc = 0.
!        else
!             write(*,*)'wrong flag'
!             stop
!        endif

!         if((wmu.gt.0.48107d0).and.(wmu.lt.0.48109d0)) then
!            write(*,*) 'af comtot xnu=',xnu
!         endif
!
!         if (xnu.gt.4.3d2) then
!            write(*,*) 'comac = ',comac
!            write(*,*) 'enexc = ',enexc
!         endif
!
         sigcom = velfact*comac*pair_enhance
         sigcomex = velfact*enexc*pair_enhance
         sigsc = sigcom
!
         if (kph.eq.1) then
            xqsqleft = rmin**2
         else
            xqsqleft = r(kph-1)**2
         endif
!         
!         
        if(scat_flag .ne. 0)then
           dcol = colmfp/sigsc
        else
           dcol = 100*dmax1(r(nr),z(nz))
        endif

        if ( dcen .le. dcol ) then
           trld = dcen
           ikind = 2
        else
          trld = dcol
          ikind = 3
        endif
!
!        if (xnu.gt.4.3d2) write(*,*) 'trld(col/cen) = ',trld
!
!______________________________________________________________________________________
!      Geometrical calculation begins:
       Eta = cos(phi)
!     Eta_switch is +1 if phi is in quadrants I or II, and -1
!     if it is in quadrants III or IV.  This way
!     the quadrant of phi is not lost when taking its cosine.  
!     J. Finke, 6 Sept. 2005
       if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
          Eta_switch = 1
       else
          Eta_switch = -1
       endif
!       if(ncycle.gt.0) 
!         if((wmu.gt.0.48107d0).and.(wmu.lt.0.48109d0)) then
!            write(*,*) 'af phi=', phi, ' eta_switch=', eta_switch
!         endif
!       
!       if (Eta.gt.9.9999999d-1) Eta = 9.9999999d-1
!       if (Eta.lt.-9.9999999d-1) Eta = -9.9999999d-1
!
       disp = Eta*rpre
       psq = rpre*rpre*(1. - Eta**2)
!
       if ((Eta.lt.0.d0) .and. (psq.lt.xqsqleft))  then
!
          incounter = incounter + 1
          kbnd = kph-1
          inout = -1
          if (kph.gt.1) then
             rbnd = r(kph-1)
          else
             rbnd = rmin
          endif
!
       else
          outcounter = outcounter + 1
          kbnd = kph
          inout = 1
          rbnd = r(kph)

       endif     
!
       dpbsq = rbnd**2 - psq
       if (dpbsq.lt.1.d-6) dpbsq = 1.d-6
!
       disbr = dble(inout)*sqrt(dpbsq) - disp
       trldb = disbr/dsqrt(1.d0 - (wmu**2))
       f     = disbr
       Zr = Zpre + wmu*trldb
!
        if (jph.eq.1) then
!        
           if  ( Zr.gt.Z(jph)) then
!
!         nearest boundary = upward z boundary
!          
              Zbnd = Z(jph)
              knew = kph
              jnew = jph+1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt (rpre**2 + f**2 + 2.d0*rpre*f*Eta)
              rbnd =  Rr     
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
                               
           else if (Zr.lt.zmin) then
!
!        nearest boundary = downward z boundary
!       
              Zbnd = zmin
              knew = kph
              jnew = jph - 1
              f = (zmin - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt (rpre**2 + f**2 + 2.d0* rpre*f*Eta)           
              rbnd = Rr
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
     
           else 
!
!        nearest boundary = r boundary
!
               knew = kph + inout
               jnew = jph
               if (kbnd.gt.0) then
                  Rr = r(kbnd)
               else
                  Rr = rmin
               endif
               rbnd = Rr
               Zbnd = Zr
              
           endif
!
        else
!       j is greater than 1   cccccccccccccccccccccc
!
           if  ( Zr.gt.Z(jph)) then
!
!         nearest boundary = upward z boundary
!          
              Zbnd = Z(jph)
              knew = kph
              jnew = jph+1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*Eta)
              rbnd =  Rr     
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
                           
           else if (Zr.lt.Z(jph-1)) then
!
!        nearest boundary = downward z boundary
!       
              Zbnd = Z(jph-1)
              knew = kph
              jnew = jph - 1
              f = (zbnd - zpre)*dsqrt(1.d0 - wmu**2)/wmu
              Rr = dsqrt(rpre**2 + f**2 + 2.d0*rpre*f*Eta)            
              rbnd = Rr
              trldb = dsqrt(f**2 + (Zbnd - Zpre)**2)
!
           else
!
!        nearest boundary = r boundary
!
               knew = kph + inout
               jnew = jph
               if (kbnd.gt.0) then
                  Rr = r(kbnd)
               else
                  Rr = rmin
               endif
               rbnd = Rr
               Zbnd = Zr

           endif
         endif
!      end of the Geometrical calculation
!      got the trldb, the travel length to the boundary.
!______________________________________________________________________________
       kappa_ic=kappa_ic+dmin1(trldb,dcen)*sigsc
       if(ncycle.eq.410.and.jph.eq.10.and.kph.eq.1.and.&
     &  xnu.gt.1e-3.and.xnu.lt.3e-3)write(*,*)'kappa_ic=',kappa_ic
       if(kappa_ic*split1.gt.0.5d0)then
         write(*,*)'optically not thin enough, reduce split1'
         write(*,*)'kappa_ic',kappa_ic,'xnu=',xnu
         stop
       endif
!_______________________________________________________________________________
!
       if (trldb.lt.trld) then
!  
           ikind = 1
           trld = trldb
           rnew = rbnd
           znew = Zbnd
!
       else 
           jnew = jph
           knew = kph
           f = trld*dsqrt(1.d0 - wmu**2)
           rnew = dsqrt(f**2 + rpre**2 + 2.d0*f*rpre*Eta)
           znew = zpre + trld*wmu
!
       endif
!
!
       do i = 1, n_vol-1
          if (xnu.lt.E_ph(i+1)) exit
       enddo
!
       if (pair_switch.eq.1)then
         do i_gg = 1, n_gg
           if (xnu.lt.E_gg(i_gg)) exit
         enddo
       endif
!
         sigabs = sigabs + pair_enhance*kappa_tot(i, jph, kph)
!
!
         if (pair_switch.eq.1) then
            if (xnu.lt.E_gg(1)) then
               kgg = (xnu/E_gg(1))*k_gg(1, jph, kph)
!            else if (i_gg.eq.n_gg) then
!               kgg = k_gg(n_gg, jph, kph)
            else if (i_gg.gt.n_gg) then
                kgg = k_gg(n_gg, jph, kph)
            else
!               kgg0 = k_gg(i_gg, jph, kph)
!               kgg1 = k_gg(i_gg+1, jph, kph)
!               kgg = kgg0 + (xnu - E_gg(i_gg))*(kgg1 - kgg0)
!     1                     /(E_gg(i_gg+1) - E_gg(i_gg))   
               kgg0 = k_gg(i_gg-1, jph, kph)
               kgg1 = k_gg(i_gg, jph, kph)
               kgg = kgg0 + (xnu - E_gg(i_gg-1))*(kgg1 - kgg0)&
     &                     /(E_gg(i_gg) - E_gg(i_gg-1))  ! Xuhui 2013/6/28
            endif
            sigabs = sigabs + kgg
         endif
         if (sigabs.lt.1.d-40) sigabs = 1.d-40
!
         xabs = sigabs*trld
!
         if (xabs.lt.100.) then
            ewnew = ew*dexp(-xabs)
        else
            ewnew = 0.
        endif
!        if ( ewnew .le. wtmin ) then
!        if(ewnew .le. wtmin/split1)then
!           npkill = npkill + 1
!           ekill  = ekill + ewnew
!           ekillt = ekillt + ewnew
!           ewnew  = 0.
!       endif
!
!       if (xnu.gt.4.d2) write(*,*) 'ewnew = ', ewnew
!
!       if ((xnu.gt.47.d0).and.(sigabs.gt.1.d-40)) then
       if((xnu.gt.47.d0).and.(sigabs.gt.1.d-40).and.&
     &       (pair_switch.eq.1))then ! Xuhui 5/19/09
!          deleabs = (ew - ewnew)*(sigabs - k_gg(i_gg, jph, kph))/sigabs
          deleabs = (ew - ewnew)*(sigabs - kgg)/sigabs ! Xuhui 2013/6/28
       else
          deleabs = ew - ewnew
       endif
       if (deleabs.lt.1.d-50) deleabs = 1.d-50
!
!
!  ewpl = energy weight path length
!
       if ( xabs .le. .00001d0 ) then
         ewpl = ew*trld*(1.d0-.5d0*xabs)
         wmustar = wmu
      else
         ewpl = deleabs/sigabs
!  119    mb_ran = drand(0)
  119    mb_ran = fibran()
         if (mb_ran.lt.(ew/deleabs)) then
            sstar = -dlog(1. - mb_ran*deleabs/ew)/sigabs
         else
            goto 119
         endif
         
         denom = dsqrt (rpre**2 + 2.*wmu*rpre*sstar + sstar**2)
         wmustar =  (wmu*rpre + sstar)/denom
      endif
!
      delpr = deleabs*wmustar*c_light
      delecomp = facomp*sigcomex*ewpl
      if(scat_flag .eq. 0)then
        edep(jph, kph) = edep(jph, kph) + deleabs ! + delecomp Xuhui 2/15/11
        prdep(jph, kph) = prdep(jph, kph) + delpr
      endif !scat_flag=0 is non-scattering tracking.
!
!      if (ewnew .le. 1.d-40) go to 900
      if(ewnew.le.wtmin) go to 900 ! Xuhui 5/26/09
!
      ew = ewnew
      dcen = dcen - trld
!      if (rnew.gt.1.d-10) then
!         Eta = (trld - Eta*rpre)/rnew
!      else      
!         Eta = (trld + Eta*rpre)/rnew ! 17/10/2012 Xuhui
          Eta = (f + Eta*rpre)/rnew

!      endif
      if (abs(Eta).gt.1.d0+1.d-15)&
     &    write(*,'("What?! Eta too extreme! Eta=",e24.17)') Eta
      if (Eta.gt.1.d0) Eta = 1.d0
      if (Eta.lt.-1.d0) Eta = -1.d0
      phi = acos(Eta)
      if (eta_switch.eq.-1) phi = 2.d0*pi - phi
       if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
          Eta_switch = 1
       else
          Eta_switch = -1
       endif 
      rpre = rnew
      Zpre = znew
!
!      call bin_add(phi)  
!       
!        if(ncycle .eq. 0 .and. myid .eq. 1)then
!          iii = int(log10(dcol))
!          write(unit=112,fmt='(1X,e14.6,1X,I1)')dcol,ikind
!        endif
!___________________________________________________
      if ( ikind .eq. 1 ) then
!     particle has reached boundary, r(kbnd) (ikind = 1)
        colmfp = colmfp - sigsc*trld  ! decrease the travel length for scattering
      
        if((jnew .eq. nz+1) .or. (jnew.eq.0) .or.&
     &     (knew .eq. nr+1) .or. (knew.eq.0)) then
!
           jph = jnew
           kph = knew
!
           if(scat_flag.eq.-1)then
              goto 900
           else
              call imcleak
              lkcount = lkcount + 1
           endif

         if (idead .eq. 1) go to 900
!         if (idead .eq. 4) then
!            kph = 0
!            go to 100
!         endif
!
          go to 110
!          go to 100 ! Xuhui 3/4/09
        else
!     
          kph = knew
          jph = jnew
          go to 110
!          go to 100 ! Xuhui 3/4/09
        endif
!___________________________________________________
      elseif (ikind .eq. 2 .and. scat_flag .ne. -1) then
!
!       case ikind = 2: write particle to census
!
!
      npcen(jph,kph) = npcen(jph,kph) + 1
      ecens(jph, kph) = ecens(jph, kph) + ew
!
!
      if (pair_switch.eq.1) then
         do i_gg = 1, n_gg-1
            if (xnu.lt.E_gg(i_gg)) exit
         enddo
         Egg_min = (E_gg(1)**2)/E_gg(2)
         if (xnu.gt.Egg_min) then
             n_ph(i_gg, jph, kph) = n_ph(i_gg, jph, kph)+(ew*6.25d8)/xnu
         endif
      endif
!
      do i_gg = 1, nphfield-1
        if (xnu.lt.E_field(i_gg+1)) exit
      enddo
      Egg_min = (E_field(1)**2)/E_field(2)
!
!
      if (xnu.gt.Egg_min) then
         n_field(i_gg, jph, kph) = n_field(i_gg, jph, kph) &
     &                       + 6.25d8*ew/xnu
      endif
!
       lwad = 6*ndxout
       lwai = 6*ndxout
       dbufout(lwad+1) = rpre
       dbufout(lwad+2) = Zpre
       dbufout(lwad+3) = wmu
       dbufout(lwad+4) = phi
       dbufout(lwad+5) = ew
       dbufout(lwad+6) = xnu
       ibufout(lwai+1) = jgpsp
       ibufout(lwai+2) = jgplc
       ibufout(lwai+3) = jgpmu
       ibufout(lwai+4) = jph
       ibufout(lwai+5) = kph
       ibufout(lwai+6) = emitype !int( fibran()*1.d5 )
       ndxout = ndxout + 1
       if(ndxout .ge. ucens) then
          write(*,*) 'too many photons'
          stop ! Xuhui cens
       endif
       go to 900
!_______________________________________________
      elseif (ikind .eq. 3) then
!
!     case ikind = 3: Compton scattering
!
       nscat = nscat + 1 ! Xuhui
!       

!*****************XUHUI*******************************************************
!      This is the second splitting of the MC particles
!      store all the MC particle information before it is splitted.
       ewcsv = ew
       dcencsv = dcen
       xnucsv = xnu
       zprecsv = zpre
       rprecsv = rpre
       wmucsv = wmu
       phicsv = phi
       jcsv = jph
       kcsv = kph
       jgpspcsv = jgpsp
       jgplccsv = jgplc
       jgpmucsv = jgpmu
       emitypecsv = emitype

       if(scat_flag.eq.1)write(*,*)'multiple scattering'
       cmcount1 = cmcount1+1 ! count the number of 2nd split scattering happend in this time step, in this node.
!       if(cmcount1.eq.1)write(*,*)'first scattering ew=',ew
       ewcsv = ewcsv/split2
       
       do 220, ii=1,split2
       ew = ewcsv
       dcen = dcencsv
       xnu = xnucsv
       zpre = zprecsv
       rpre = rprecsv
       wmu = wmucsv
       phi = phicsv
       jph = jcsv
       kph = kcsv
       jgpsp = jgpspcsv
       jgplc = jgplccsv
       jgpmu = jgpmucsv
       emitype = emitypecsv
       
       
       ewold = ew
       call compb2d(i_gam)
!_______________________________________
!      This is the 3rd splitting of the MC particles.
!      This happens when the MC particle is scattered to a very high energy.
         if(ew.gt.ewold*split2*split1*spl3_trg)then
           ewold = ewcsv/split3
           do 217, ii2=1,split3
215          continue
             ew = ewcsv/split3
             dcen = dcencsv
             xnu = xnucsv
             zpre = zprecsv
             rpre = rprecsv
             wmu = wmucsv
             phi = phicsv
             jph = jcsv
             kph = kcsv
             jgpsp = jgpspcsv
             jgplc = jgplccsv
             jgpmu = jgpmucsv
             emitype = emitypecsv
             call compb2d(i_gam)
            if(ew.le.ewold*split2*split1*spl3_trg)goto 215
            edep(jph,kph) = edep(jph,kph) + ew - ewold + deleabs ! Xuhui 2013-07-15
            prdep(jph, kph) = prdep(jph, kph) + delpr
            E_IC(i_gam) = E_IC(i_gam) + ew - ewold
!            cmener = cmener + ew - ewold
            cmcount2 = cmcount2 +1
            if(phi.gt.2.d0*pi) phi = phi - 2.d0*pi
            if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
               Eta_switch = 1
            else
               Eta_switch = -1
            endif
            emitype = emitype+1
            call imctrk2d(-1)
217       continue
!___________________________________________________________________
         else
           edep(jph,kph) = edep(jph,kph) + ew - ewold + deleabs ! Xuhui 2013-07-15
           prdep(jph, kph) = prdep(jph, kph) + delpr
           E_IC(i_gam) = E_IC(i_gam) + ew - ewold
!         cmener = cmener + ew - ewold
!       if (ew .lt. wtmin)goto 220 ! Xuhui 5/18/09
!      if the energy weight is below wtmin, then kill the photon.
           cmcount2 = cmcount2 +1 ! Xuhui

           if(phi.gt.2.d0*pi) phi = phi - 2.d0*pi
           if(( phi.le.pi).and.(phi.ge.1.d-10) ) then
            Eta_switch = 1
           else
            Eta_switch = -1
           endif
           emitype = emitype+1
           call imctrk2d(-1)
         endif
       
220    continue 
!
!      goto 100
      endif
 900  continue
!############################################################################################
!     combine the photons that are not scattered in the first splitting.
      if((split1-nscat).gt.0.and.scat_flag.eq.-1) then
       ew = (split1-nscat)*ewtsv
       dcen = dcentsv
       xnu = xnutsv
       zpre = zpretsv
       rpre = rpretsv
       wmu = wmutsv
       phi = phitsv
       jph = jtsv
       kph = ktsv
       jgpsp = jgpsptsv
       jgplc = jgplctsv
       jgpmu = jgpmutsv
       emitype = emitypetsv
       call imctrk2d(0)
      endif

!
      return
      end    
!

       
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:26:40 EDT 2006
! version: 2
! Name: J. Finke
! Changed common block 'random'.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun 16 12:31:21 EDT 2006
! version: 3
! Name: J. Finke
! seed is now stored in census files. imcfield 
! is now determinable.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc
! Mon Sep 22 21:33  CDT 2008
! Name : Xuhui Chen
! changed the code to split the photon before it is
! determined whether the photon will be compton
! scattered, and recombine the unscattered ones, and
! further split the scattered one (twice). This is to avoid
! the lack of the number of high energy photon.
      
