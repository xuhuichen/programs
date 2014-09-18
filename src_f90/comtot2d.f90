!      calculate the Compton cross section and hence the inverse of the free travel length of a photon 
!      in the medium of electrons in one particular zone.
      subroutine comtot(j, k, xnuc, ienexc,comac,enexc)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
!      name of routine...comtot
!      purpose...
!
!        input:
!           xnuc   = photon frequency (kev)
!           tcomp  = electron temperature (kev)
!           j      = zone -1 (for zfree, abar, etc.)
!           ienexc = 0 calculate comac
!                  = 1 calculate comac and enexc.
!
!        output:
!           comac = hot compton sigma (1/cm)
!           enexc = expected value of comac*fractional energy change
!                   (1/cm)*(kev/kev)
!
!        these are hot,macroscopic values interpolated from the
!        microscopic data-loaded comp0,enxtab.
!
!        comp0(i) = cold compton cross-section (cm**2).
!             i=1,201   xnuc = 5.*(i-1) kev.   0. .le. xnu .le. 1000.
!
!        enxtab(i,k) = hot microscopic compton cross section times
!        expected value of fractional enery loss.  (cm**2*kev/kev)
!        positive value means photon is losing energy.
!        i=1,13    tea = 10.*(i-1) kev.  tea = 0. ( 10 kev ) 120.
!        k=1,52 has equal log spacing in xnuc.
!        xnuc(1)=0. xnuc(2)=2.3142 kev.  ......xnuc(52)=emasskev.
!        xnuc(k+1) - xnuc(k) = 1.05 * ( xnuc(k) - xnuc(k-1) )
!        k=53 is for xnuc = 1000. kev
!
!
!
!        Inclusion of hot cross sections for nonthermal distributions:
!                         M. Boettcher, 07/03/97
!
!        For the total cross section, a fit to the correction factor
!                    <sigma> = (1 - f)*sigma_cold 
!                           is used where
!
!                                          (-1.9 p)
!        f = 0.48*(0.66*log10<gamma> - 15 e         )*
!
!              (-60/xnu^1.5)   0.11      -(xnu/<gamma>)^(1/2)
!           *(e            *xnu       + e                    )
!
!  <gamma> = gbar_nth = average Lorentz factor of nonthermal distribution;
!        p =    p_nth = Spectral index of nonthermal distribution.
!
!
!         Tables for hot cross section times expected value of
!         fractional energy loss (subroutine imcdate) generated
!         by numerical integration. (M. Boettcher, 07/02/97)
!
!         exn_nth(i, j, k, l) in units of (cm^2 keV/keV)
!
!    Index i: Power-law index: p_ind = 1.6 + 0.2*(i-1)
!                                    = 1.6, ..., 3.6
!                                      (i = 1, ..., 11)
!
!    Index j: low-energy cut-off: gmin = 1.2^(j-1)
!                                      = 1, 1.2, ..., 2.48832
!                                        (j = 1, ..., 6)
!
!    Index k: high-energy cut-off: gmax = 100*2^(k-1)
!                                       = 100, 200, ..., 1.28e4
!                                        (k = 1, ..., 8)
!
!    Index l: Photon energy: xnu(0) = 0.
!                            xnu(1) = 2.314177408 keV
!                               .
!                               .
!                               .
!
!                            xnu(l+1) = xnu(l) + 1.05*(xnu(l) - xnu(l-1))
!                                      (l = 1, ..., 66)
!                               .
!                               .
!                               .
!
!                            xnu(66) = 1057.1116 keV
!
!
!
!/c
! 07/20/82 ehc-deleted prexc and pxtab1.
! 03/11/83 ehc-trapped cosig .le. 0.  changed first fraci to frai.
! 05/22/85 ehc-added input arguments. renamed some local variables.
! 02/19/86 ehc-put wienke fit for sig(nu,te), if icoms=4, for howard-liang.
!              corrected cosig .lt. 0 fix.
! 03/04/86 ehc-added icoms=5 for liang (use wienke fit for 3-d maxwell).
! 08/22/88 ehc-fixed 5/85 (probably) error. missing parentheses shifted
!         enexc to low end of freq. interpolation interval.
!        for ordinary compton, energy deposit in matter was too small.
!         for inverse compton, error in energy removed from matter varied.
!
      integer j, k, ienexc
      double precision xnuc
!
      integer icoms, icmp, itea, kxnu, i_p, i_min, i_max
!
      double precision cosig, enexc_th, frac_min, enexc_nt
      double precision abar, zfree, exan, emasskev, alf, alf21
      double precision alf22, alf1, alf2, cosigth, sigthom
      double precision xbin, frai, cosignth, T_corr, compfac
      double precision tbin, fraci, frack, wva, wvb, enex_th
      double precision x_p, fracp, x_min, x_max, frac_max
      double precision p0j0k0, p0j0k1, p0j1k0, p0j1k1, p0j0
      double precision p0j1, p1j0k0, p1j0k1, p1j1k0, p1j1k1
      double precision p1j0, p1j1, wp0, wp1, enex_nt, tcomp
      double precision comac, enexc
      integer i
      double precision x,csigma
      double precision sigma_E,intg_v,gamma0,betta
!       c           c            c            c
!
!
      tcomp = tea(j,k)
      icoms = 6
      abar = 1.d0
      zfree = 1.d0
      exan = 1.0d-72
      emasskev = 5.11d2
      sigthom = 6.6516d-25
      compfac = n_e(j,k)*zfree/abar
!
      if(icoms .eq. 5)then
!___________________________________________________________________________
!
!     b. wienke (pvt. communication, early 1985) gave hot sigma fit,
!     good to 3% for all xnu, te in (0, 5000 kev).
!     note: fit gives cosig in barns, not cm**2.  use avogad24.
!
      alf = xnuc/emasskev
      if(xnuc .le. 1.0d2) then
         alf21 = .032507 + 59.9954*alf - 186.572*alf**2
      else
         alf21 = .14912*(89.6552*alf**.1324 - 40.9714)
      endif
      alf22 = 2.07114*(1.-dexp(-10.161*alf))**1.0885
      alf1 = 1. - 12.986*alf + 14.2996*alf**.9809
      alf2 = 1. + alf21*(tcomp/5000.) - alf22*(tcomp/5000.)**2
      cosigth = sigthom/(alf1*alf2)
!
      xbin = dmin1((2.d-1*xnuc + 1.d0), 2.d2)
      icmp = xbin
      frai = xbin - icmp
      T_corr = 1. - .48*(.66*log10(gbar_nth(j,k)) - 15.*&
     &        exp(-19.*p_nth(j,k)))*(exp(-60./(xnuc**1.5))*xnuc**.11&
     &        + exp(-sqrt(xnuc/gbar_nth(j,k))))
!      cosignth = T_corr*(comp0(icmp) + frai*(comp0(icmp+1) 
!     1                 - comp0(icmp)))
      x = xnuc/emasskev
      if(x .lt. 1.d-3)then
        csigma = 6.65d-25*(1.d0-2.d0*x+26.d0*x*x/5.d0)
      else
        csigma = 6.65d-25 *3.d0/4.d0* (  (1.d0+x)/x**3.d0 &
     &        *(2.d0*x*(1.d0+x)/(1.d0+2.d0*x)-dlog(1.d0+2.d0*x))&
     &        + 1.d0/2.d0/x*dlog(1.d0+2.d0*x)&
     &        - (1.d0+3.d0*x)/(1.d0+2.d0*x)**2.d0 )  
      endif
      cosignth = T_corr* csigma     ! Xuhui 4/13/09
      cosig = amxwl(j,k)*cosigth + (1. - amxwl(j,k))*1.d24*cosignth
!
!      compfac = n_e(j,k)*zfree/abar
      if(cosig .le. 0.d0) then
         comac = 1.d-99
      else
         comac = compfac*cosig
      endif
      enexc = 0.d0
      go to 900
!
!        interpolate in xnuc.  assumes comp0 constant for xnuc .ge. 995kev.
!
  100 continue
!___________________________________________________________________________
      else if(icoms.lt.4)then
      xbin = dmin1((2.d-1*xnuc + 1.d0), 2.d2)
      icmp = xbin
      frai = xbin - icmp
!
!        apply tea correction to cold data.
!
!      compfac = n_e(j,k)*zfree/abar
      T_corr = 1.d0 - tcomp*xnuc/(4.7703d+04 + xnuc*6.3769d+02)
      cosigth = T_corr*(comp0(icmp) + frai*(comp0(icmp+1) &
     &                - comp0(icmp)))
     
      T_corr = 1.d0 - .48d0*(.66d0*log10(gbar_nth(j,k)) - 1.5d1*&
     &   exp(-1.9d0*p_nth(j,k)))*(exp(-6.d1/(xnuc**1.5))*xnuc**.11&
     &   + exp(-sqrt(xnuc/gbar_nth(j,k))))
!      cosignth = T_corr*(comp0(icmp) + frai*(comp0(icmp+1) 
!     1                 - comp0(icmp)))
      x = xnuc/551.d0
      if(x .lt. 1.d-3)then
        csigma = 6.65d-25*(1.d0-2.d0*x+26.d0*x*x/5.d0)
      else
        csigma = 6.65d-25 *3.d0/4.d0* (  (1.d0+x)/x**3.d0 &
     &        *(2.d0*x*(1.d0+x)/(1.d0+2.d0*x)-dlog(1.d0+2.d0*x))&
     &        + 1.d0/2.d0/x*dlog(1.d0+2.d0*x)&
     &        - (1.d0+3.d0*x)/(1.d0+2.d0*x)**2.d0 )  
      endif
      cosignth = T_corr*csigma
      cosig = amxwl(j,k)*cosigth + (1.d0 - amxwl(j,k))*cosignth
      if(cosig .le. 1.d-30) then
         comac = 1.d-30
      else
         comac = n_e(j,k)*cosig
      endif
!_____________________________________________________________________________
      else if(icoms.eq.6)then
       cosig = 0.d0
       x =xnuc/emasskev
       do i=1,num_nt-1
         gamma0 = gnt(i)+1.d0
         betta = dsqrt(1.d0-1.d0/gamma0**2)
!       D_1 = 1.d0 - betta*mu
!       I_0 =4.d0*Pi*gamma0**2
!       I_1 =4.d0*Pi*gamma0**2/(1.d0+2.d0*x*gamma0*D_1)
!       I_2 =2.d0*Pi/(alpha_c*D_1)*dlog(1.d0+2.d0*x*gamma0*D_1)
!       I_3 =4.d0*Pi*gamma0**4*alpha_c*D_1/(1.d0+2.d0*x*gamma0*D_1)**2
!       sigma_v = 0.625d0/Pi*sigthom*gamma0**2*D_1*
!     1       (  0.5d0*I_0-(I_2-I_1)/(alpha_c*D_1*gamma0**2)+
!     1          (I_0+I_1-2.d0*I_2)/(2.d0*alpha_c**2*D_1**2*gamma0**4)+
!     1          0.5d0*(I_2-I_3)  )
          if(x*gamma0*(1+betta).lt.1.d-2)then
             sigma_E = sigthom*(1.d0-2.d0*x*gamma0)
          else
             sigma_E = 9.375d-2*sigthom/gamma0**2/betta/x**2*&
     &  (intg_v(2*gamma0*(1+betta)*x)-intg_v(2*gamma0*(1-betta)*x))
          endif
          cosig = cosig + sigma_E*f_nt(j,k,i)*(gnt(i+1)-gnt(i))
       enddo
       if(cosig .lt. 1.d-40) then
         comac = 1.d-40
       else
         comac = n_e(j,k)*cosig
       endif
      endif
!
!``````````````````````````````````````````````````````````````````````````
      if(ienexc .ne. 0) then
!
!        interpolate in tea and xnuc for expected values.
!        for tea .gt. 120. and/or xnuc .gt. 1000., extrapolate.
!
       tbin = 1.d-1*tcomp + 1.d0
       itea = tbin
       if( itea .gt. 12 ) itea = 12
       fraci = tbin - itea
       if( xnuc .ge. emasskev ) then
         kxnu = 52
         frack = (xnuc-emasskev)/4.88994d2
       else
!
!        solve for frequency -1, xbin .lt. 52.
!        xnuc = 2.3142*( 1.05**(xbin-1.) - 1. )/.05
!
          xbin = 1.d0 + 2.04959d1*log(1.d0 + 2.16059d-2*xnuc)
          kxnu = xbin
          frack = xbin - kxnu
       endif
       wva = fraci*( enxtab(itea+1,kxnu  ) - enxtab(itea,kxnu  ) )&
     &     +enxtab(itea,kxnu  )
       wvb = fraci*( enxtab(itea+1,kxnu+1) - enxtab(itea,kxnu+1) )&
     &     +enxtab(itea,kxnu+1)
       enexc_th = compfac*(wva+frack*(wvb-wva))
!
       x_p = 5.d0*(p_nth(j,k) - 1.6d0) + 1.d0
       i_p = x_p
       if (i_p.gt.25) i_p = 25
       if (i_p.lt.1) i_p = 1
       fracp = x_p - i_p
       x_min = 5.484815d0*dlog(gmin(j,k)) + 1.d0
       i_min = x_min
       if (i_min.gt.5) i_min = 5
       if (i_min.lt.1) i_min = 1
       frac_min = x_min - i_min
       x_max = 1.442695d0*dlog(gmax(j,k)*.01) + 1.d0
       i_max = x_max
       if (i_max.gt.7) i_max = 7
       if (i_max.lt.1) i_max = 1
       frac_max = x_max - i_max
!
       p0j0k0 = enx_nth(i_p, i_min, i_max, kxnu) + frack*&
     &           (enx_nth(i_p, i_min, i_max, kxnu+1) &
     &          - enx_nth(i_p, i_min, i_max, kxnu))
       p0j0k1 = enx_nth(i_p, i_min, i_max+1, kxnu) + frack*&
     &           (enx_nth(i_p, i_min, i_max+1, kxnu+1) &
     &          - enx_nth(i_p, i_min, i_max+1, kxnu))
       p0j0 = p0j0k0 + frac_max*(p0j0k1 - p0j0k0)
!
       p0j1k0 = enx_nth(i_p, i_min+1, i_max, kxnu) + frack*&
     &           (enx_nth(i_p, i_min+1, i_max, kxnu+1) &
     &          - enx_nth(i_p, i_min+1, i_max, kxnu))
       p0j1k1 = enx_nth(i_p, i_min+1, i_max+1, kxnu) + frack*&
     &           (enx_nth(i_p, i_min+1, i_max+1, kxnu+1) &
     &          - enx_nth(i_p, i_min+1, i_max+1, kxnu))
       p0j1 = p0j1k0 + frac_max*(p0j1k1 - p0j1k0)
       wp0 = p0j0 + frac_min*(p0j1 - p0j0)
!
       p1j0k0 = enx_nth(i_p+1, i_min, i_max, kxnu) + frack*&
     &           (enx_nth(i_p+1, i_min, i_max, kxnu+1) &
     &          - enx_nth(i_p+1, i_min, i_max, kxnu))
       p1j0k1 = enx_nth(i_p+1, i_min, i_max+1, kxnu) + frack*&
     &           (enx_nth(i_p+1, i_min, i_max+1, kxnu+1) &
     &          - enx_nth(i_p+1, i_min, i_max+1, kxnu))
       p1j0 = p1j0k0 + frac_max*(p1j0k1 - p1j0k0)
!
       p1j1k0 = enx_nth(i_p+1, i_min+1, i_max, kxnu) + frack*&
     &           (enx_nth(i_p+1, i_min+1, i_max, kxnu+1) &
     &          - enx_nth(i_p+1, i_min+1, i_max, kxnu))
       p1j1k1 = enx_nth(i_p+1, i_min+1, i_max+1, kxnu) + frack*&
     &           (enx_nth(i_p+1, i_min+1, i_max+1, kxnu+1) &
     &          - enx_nth(i_p+1, i_min+1, i_max+1, kxnu))
       p1j1 = p1j1k0 + frac_max*(p1j1k1 - p1j1k0)
       wp1 = p1j0 + frac_min*(p1j1 - p1j0)
       enexc_nt = compfac*(wp0 + fracp*(wp1 - wp0))
!
       enexc = amxwl(j,k)*enexc_th + (1.d0 - amxwl(j,k))*enexc_nt
!
      endif
!``````````````````````````````````````````````````````````````````````````````
  900 continue
      return
      end


      double precision function intg_v(x)
!
!     This function calculate the integral of eqn(2.3) in Coppi et al. 1990
!     Xuhui Chen 5/4/2009
!
      implicit none
      double precision dilog, x
      double precision intg_1, intg_2, intg_3
!
      intg_1 = -x/2.d0+0.5d0/(1.d0+x)
      intg_2 = 4.d0*dilog(-x)
      intg_3 = (9.d0+x+8.d0/x)*dlog(1.d0+x)
!
      intg_v = (intg_1 + intg_2 + intg_3)
      return 
      end



      double precision function dilog(x)
!
!     The DiLogarithm function
!     Code translated by R.Brun from CERNLIB DILOG function C332
!     This fortran version is transplanted from a C++ version in ROOT
!     Xuhui Chen 5/4/2009
!
      implicit none
      include 'general.pa'
      double precision x, HF, PI2, PI3, PI6, PI12, C(20)
      parameter(HF  = 0.5d0)
      parameter(PI2 = pi*pi)
      parameter(PI3 = PI2/3)
      parameter(PI6 = PI2/6)
      parameter(PI12 = PI2/12)
      data C /0.42996693560813697d0, 0.40975987533077105d0,&
     &     -0.01858843665014592d0, 0.00145751084062268d0,&
     &     -0.00014304184442340d0,&
     &      0.1588415541880d-4,-0.190784959387d-5, 0.024195180854d-5,&
     &     -0.003193341274d-5, 0.000434545063d-5,-0.000060578480d-5,&
     &      0.000008612098d-5,-0.000001244332d-5, 0.000000182256d-5,&
     &     -0.000000027007d-5, 0.000000004042d-5,-0.000000000610d-5,&
     &      0.000000000093d-5,-0.000000000014d-5, 0.000000000002d-5/
      double precision T,H,Y,S,A,ALFA,B1,B2,B0
      integer i

!
      if(x.eq.1)then
         H = PI6
      else if(x.eq.-1)then
         H = -PI12
      else
         T = -x
         if(T.le.-2)then
            Y = -1/(1+T)
            S = 1
            B1 = dlog(-T)
            B2 = dlog(1+1/T)
            A = -PI3+HF*(B1*B1-B2*B2)
          else if(T.lt.-1)then
            Y = -1-T
            S = -1
            A = dlog(-T)
            A = -PI6+A*(A+dlog(1+1/T))
          else if(T.le.-0.5d0)then
            Y = -(1+T)/T
            S = 1
            A = dlog(-T)
            A = -PI6+A*(-HF*A+dlog(1+T))
          else if (T .lt. 0)then
                 Y = -T/(1+T)
                 S = -1
                 B1= dlog(1+T)
                 A = HF*B1*B1
          else if (T .le. 1)then
                 Y = T
                 S = 1
                 A = 0
          else
                 Y = 1/T
                 S = -1
                 B1= dlog(T)
                 A = PI6+HF*B1*B1
          endif
             H    = Y+Y-1
             ALFA = H+H
             B1   = 0
             B2   = 0
             do i=20,1,-1 
                B0 = C(i) + ALFA*B1-B2
                B2 = B1
                B1 = B0
             enddo
             H = -(S*(B0-H*B2)+A)
      endif
      dilog = H
      return
      end

!ccccccccccccccccccccccccccccccccccccc
! Wed Sep 24, 2008 CDT
! Xuhui Chen
! Aborted use of common area 'sigmas'. Used the
! passing parameter in stead.
!ccccccccccccccccccccccccccccccccccccc
! Sat May 2, 2009 CDT
! Xuhui Chen
! added icoms=6 case for computing the hot Compton
! cross section for the nonthermal electron
! distribution, according to Hamada et al. 1978 PTP.
