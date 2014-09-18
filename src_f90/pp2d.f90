!
!      These subroutines calculate the
!     differential pair production rate
!    
!
      subroutine pairprod(j, k)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
!
!
      integer i_p1, i_p2, i_g
!
      double precision gamma, sum, E1, eps1, E2, eps2, sum1, sd
      double precision f_pprod, emin
!
!
      if (j.eq.1.and.k.eq.1) &
     &   open(28, file='dn_e1.dat', status='unknown')
      if (j.eq.2.and.k.eq.1) &
     &   open(29, file='dn_e2.dat', status='unknown')
!
      do 500 i_g = 1, num_nt
         gamma = gnt(i_g) + 1.0d0
         sum = 0.d0
!
         do 400 i_p1 = 1, n_gg-1
            if (n_ph(i_p1, j, k).lt.1.d-5) goto 400
            E1 = E_gg(i_p1)
            eps1 = 1.957d-3*E1
            sum1 = 0.d0
            emin = dmax1((1.d0/eps1), (gamma + 1.d0 - eps1))
!
            do 300 i_p2 = 1, n_gg-1
               E2 = E_gg(i_p2)
               eps2 = 1.957d-3*E2
               if (eps2.lt.emin) goto 300
               if (n_ph(i_p2, j, k).lt.1.d-5) goto 300
!
               sd = f_pprod(eps1, eps2, gamma)
               sum1 = sum1 + sd*(E_gg(i_p2+1) - E2)&
     &                      *n_ph(i_p2, j, k)/(eps2**2)
  300       continue
            sum = sum + sum1*(E_gg(i_p1+1) - E1)*n_ph(i_p1, j, k)&
     &                 /(eps1**2)
  400    continue
!
         dn_pp(j, k, i_g) = 1.496d-14*sum
!
         if (j.eq.1.and.k.eq.1) then
           write(28, 600) gamma, dmax1(1.d-20, dn_pp(j, k, i_g))
         else if (j.eq.2.and.k.eq.1) then
           write(29, 600) gamma, dmax1(1.d-20, dn_pp(j, k, i_g))
         endif
!
  500 continue
!
  600 format (e10.3,1x,e10.3)
      if (j.eq.1.and.k.eq.1) then
         close(28)
      else if (j.eq.2.and.k.eq.1) then
         close(29)
      endif
!
      return
      end
!
!
!
      double precision function f_pprod(eps1, eps2, gamma)
      implicit none
      double precision eps1, eps2, gamma
!
      double precision f_inner, ecm_U, ecm_L, x
      double precision E, edagger, estar, Det
      double precision Det2, estar2, edagger2
!
      E = eps1 + eps2
      x = gamma*(E - gamma)
      Det2 = (x + 1.d0)**2 - E**2
      if (Det2.lt.0.d0) then
         f_pprod = 0.d0
         goto 900
      endif
      Det = dsqrt(Det2)
!
      estar2 = 5.d-1*(x + 1.d0 + Det)
      edagger2 = 5.d-1*(x + 1.d0 - Det)
      if ((estar2.lt.0.d0).or.(edagger2.lt.0.d0)) then
         f_pprod = 0.d0
         goto 900
      endif
      estar = dsqrt(estar2)
      edagger = dsqrt(edagger2)
!        
      ecm_U = dmin1(dsqrt(eps1*eps2), estar)
      ecm_L = dmax1(1.d0, edagger)
!     
      f_pprod = f_inner(ecm_U, eps1, eps2, gamma)&
     &        - f_inner(ecm_L, eps1, eps2, gamma)
!
 900  return
      end
!
!
!
      double precision function f_inner(ecm, eps1, eps2, gamma)
      implicit none
      double precision ecm, eps1, eps2, gamma
!
      double precision E, f1, f12, H
!
      E = eps1 + eps2
!
      f12 = E**2 - 4.d0*(ecm**2)
      if (f12.lt.0.d0) then
        f_inner = 0.d0
        goto 900
      endif
      f1 = 2.5d-1*dsqrt(f12)
!
      f_inner = f1 + H(ecm, eps1, eps2, gamma)&
     &             + H(ecm, eps2, eps1, gamma)
!
 900  return
      end
!
!
!
      double precision function H(ecm, eps1, eps2, gamma)
      implicit none
      double precision ecm, eps1, eps2, gamma
!
      double precision c, d, d2, I_pm, ee
!
      ee = eps1*eps2
      c = (eps1 - gamma)**2 - 1.d0
      d = eps1**2 + ee + gamma*(eps2 - eps1)
      d2 = ee + c*(ecm**2)
      if (d2.le.0.d0) then
         H = 0.d0
         goto 900
      endif
!
      if (dabs(c).gt.1.d-10) then
         H = -1.25d-1*ecm*(d/ee + 2.d0/c)/dsqrt(d2)&
     &     + 2.5d-1*(2.d0 - (ee - 1.d0)/c)&
     &             *I_pm(ecm, eps1, eps2, c)&
     &     + 2.5d-1*dsqrt(d2)*(ecm/c + 1.d0/(ecm*ee))
      else
         H = ((ecm**3)/1.2d1 - 1.25d-1*ecm*d)/(ee**1.5)&
     &     + ((ecm**3)/6.d0 + 5.d-1*ecm + 2.5d-1/ecm)&
     &       /dsqrt(ee)
      endif
!
 900  return
      end
!
!
!
!
      double precision function I_pm(ecm, eps1, eps2, c)
      implicit none
      double precision ecm, eps1, eps2, c
!
      double precision d2
!
      if (c.gt.1.d-40) then
         d2 = eps1*eps2 + c*(ecm**2)
         I_pm = dlog(ecm*dsqrt(c) + dsqrt(d2))/dsqrt(c)
      else if (c.lt.-1.d-40) then
         d2 = -c/(eps1*eps2)
         I_pm = asin(ecm*dsqrt(d2))/dsqrt(-c)
      else
         I_pm = 0.d0
      endif
!
 900  return
      end
!
!
!      This subroutine calculates
!      the pair annihilation rates
!
!
      subroutine pa_calc(j, k, ne_local)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision ne_local
!
!
      integer i
!
      double precision gamma, pa_el, pa_pos,&
     &                 Delta_ne, Delta_np, Delta_g
!
!
      if (k.eq.1) then
         if (j.eq.1) open(27, file='pa_1.dat', status='unknown')
         if (j.eq.2) open(28, file='pa_2.dat', status='unknown')
      endif
!
      Delta_ne = 0.d0
      Delta_np = 0.d0
      do 100 i = 1, num_nt
         gamma = 1.d0 + gnt(i)
         if (f_nt(j, k, i).gt.0.d0) then
            dne_pa(j, k, i) = -ne_local*f_nt(j, k, i)*pa_el(j, k, gamma)
         else
            dne_pa(j, k, i) = 0.d0
         endif
         if (n_pos(j, k, i).gt.0.d0) then
            dnp_pa(j, k, i) = -n_pos(j, k, i)*ne_local*pa_pos(j,k,gamma)
         else
            dnp_pa(j, k, i) = 0.d0
         endif
!
         if (k.gt.1) goto 50
         if (j.eq.1) then
            write(27, 200) gamma, dmax1(1.d-20, -dne_pa(j, k, i)), &
     &                     dmax1(1.d-20, -dnp_pa(j, k, i))
         else if (j.eq.2) then
            write(28, 200) gamma, dmax1(1.d-20, -dne_pa(j, k, i)), &
     &                     dmax1(1.d-20, -dnp_pa(j, k, i))
         endif
 50      continue
!
         if (i.lt.num_nt) then
            Delta_g = gnt(i+1) - gnt(i)
            Delta_ne = Delta_ne + dne_pa(j, k, i)*Delta_g
            Delta_np = Delta_np + dnp_pa(j, k, i)*Delta_g
         endif
 100  continue
!
 200  format (e12.3,1x,e12.3,1x,e12.3)
!
      if (k.eq.1) then
         if (j.eq.1) then
            close(27)
         else if (j.eq.2) then
            close(28)
         endif
      endif
!
      return
      end
!
!
!
!
      double precision function pa_el(j, k, ge)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision ge
!
!
      integer i
!
      double precision gp
      double precision Delta_g, sum, vsigma
!
!
      sum = 0.d0
      do 100 i = 1, num_nt-1
         if (n_pos(j, k, i).gt.0.d0) then
            gp = gnt(i) + 1.d0
            Delta_g = gnt(i+1) - gnt(i)
            sum = sum + n_pos(j, k, i)*Delta_g*vsigma(ge, gp)
         endif
 100  continue
      pa_el = sum
!
      return
      end
!
!
!
!
      double precision function pa_pos(j, k, gp)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision gp
!
!
      integer i
!
      double precision ge
      double precision Delta_g, sum, vsigma
!
!
      sum = 0.d0
      do 100 i = 1, num_nt-1
         if (f_nt(j, k, i).gt.0.d0) then
            ge = gnt(i) + 1.d0
            Delta_g = gnt(i+1) - gnt(i)
            sum = sum + f_nt(j, k, i)*Delta_g*vsigma(ge, gp)
         endif
 100  continue
      pa_pos = sum
!
      return
      end
!
!
!
      double precision function vsigma(ge, gp)
      implicit none
      double precision ge, gp
!
      double precision gcm_min, gcm_max, f_vs, be, bp
      double precision gmin2, gmax2
!
      be = dsqrt(1.d0 - 1.d0/(ge**2))
      bp = dsqrt(1.d0 - 1.d0/(gp**2))
      gmin2 = 5.d-1*(1.d0 + ge*gp*(1.d0 - be*bp))
      if (gmin2.gt.1.00002d0) then
         gcm_min = dsqrt(gmin2)
      else
         gcm_min = 1.00001d0
      endif
      gmax2 = 5.d-1*(1.d0 + ge*gp*(1.d0 + be*bp))
      if (gmax2.gt.1.00002d0) then
         gcm_max = dsqrt(gmax2)
      else
         gcm_max = 1.00001d0
      endif
!
      if (gcm_max.gt.gcm_min) then
         vsigma = 7.48d-15*(f_vs(gcm_max) - f_vs(gcm_min))&
     &           /(be*bp*((ge*gp)**2))
      else
         vsigma = 0.d0
      endif
!
      return
      end
!
!
!
      double precision function f_vs(gcm)
      double precision gcm
!      
      double precision bcm, L
!
      bcm = dsqrt(1.d0 - 1.d0/(gcm**2))
      L = dlog((1.d0 + bcm)/(1.d0 - bcm))
!
      f_vs = (bcm**3)*(gcm**2)*L - 2.d0*(gcm**2) + 7.5d-1*(L**2)
!
      return
      end
!
!
!      This subroutine smoothes the calculated
!       internal photon spectrum by fitting a
!     spectral shape n(E) = N E^(-alpha) e^(-E/E_0)
!
!
      subroutine nph_smooth(j, k, Te)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision Te

!
      integer i, k1, l, m, n1, n2
!
      double precision f_s, N, a, E0, chi2, chi2min, y
      double precision Nmin, E0min, amin, N0, E00, a0, dchi2
!
!
!
      n1 = 2
      n2 = 10
!
 100  if ((n_ph(n2, j, k).le.1.d0).or.(n_ph(n1, j, k).le.1.d0)) &
     &     goto 900
!
      a0 = dlog(n_ph(n1, j, k)/n_ph(n2, j, k))&
     &    /dlog(E_gg(n2)/E_gg(n1))
      if (((a0.lt.0.d0).or.(a0.gt.4.0d0)).and.(n2.lt.15)) then
         n2 = n2 + 1
         goto 100
      endif
      if (n2.eq.15) goto 900
      if (a0.lt.1.d-2) a0 = 1.d-2
      if (a0.gt.4.0d0) a0 = 4.0d0
      N0 = n_ph(3, j, k)
      E00 = Te
!
!      write(*,*) 'a0 = ',a0
!      write(*,*) 'N0 = ',N0
!      write(*,*) 'E00 = ',E00
!
      chi2min = 1.d50
!
      N = 5.d-1*N0
      do 450 k1 = 1, 21
         a = a0 - 5.d-1
         do 400 l = 1, 13
            E0 = 3.5d-1*E00
            do 350 m = 1, 16
               chi2 = 0.d0
               do 300 i = 1, n_gg
                  y = E_gg(i)/E0
                  if (y.lt.2.d1) then
                     f_s = N/(((E_gg(i)/E_gg(3))**a)*dexp(y))
                  else
                     f_s = 0.d0
                  endif
                  if ((f_s.gt.1.d0).and.&
     &                (n_ph(i, j, k).gt.1.d0)) then
                     dchi2 = ((n_ph(i, j, k) - f_s)**2)/f_s
                     chi2 = chi2 + dchi2
                  endif
 300           continue
!
               if (chi2.le.chi2min) then
                  E0min = E0
                  amin = a
                  Nmin = N
                  chi2min = chi2
               endif
               E0 = E0*1.15
 350        continue
            a = a + 5.d-2
 400     continue        
         N = N*1.075d0
 450  continue
!
      N = Nmin
      a = amin
      E0 = E0min
!
      do 500 i = 1, n_gg
         y = E_gg(i)/E0
         if (y.lt.2.d1) then
            n_ph(i, j, k) = N/(((E_gg(i)/E_gg(3))**a)*dexp(y))
         else
            n_ph(i, j, k) = 0.d0
         endif
 500   continue
!
!      write(*,*) 'Zone ',j,' N = ',N
!      write(*,*) ' a = ',a,'; E0 = ',E0
!
 900  return
      end
!
