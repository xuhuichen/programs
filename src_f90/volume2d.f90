!
!       This subroutine calculates the emissivities
!       (cumulative probability distributions)
!       and absorption opacities due to thermal
!       bremsstrahlung, cyclotron, non-thermal
!       synchrotron, and pair annihilation
!       emission. (MB/03/Dec/99)
!
!       cyclotron and thermal bremsstrahlung have been deactivated for 
!       the purpose of blazar (XC 2014.02.07)
!	the results of this subroutine includes eps_a, eps_tot, eps_th and kappa_tot
!
!
      subroutine volume_em(zone)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer zone, j, k
      double precision T_keV, ne_local, B, l_min
!
      integer i, m, n_harmonics, i_el, i_pos, i_a
!
      double precision mm, f_m
      double precision P(n_vol), P_sum, Theta, kappa_C
      double precision dE, T, nu, E, x, y, G_ff, sum_sd, kappa_br
      double precision kappa_sy, kappa_cy, j_sy,j_br, j_cy
      double precision sum_k, xgg, g, dg1, sd, sd_k, p_1
      double precision N_e_nt, f_cy
      double precision sum_th, P_th(n_vol), tau_tot, j_th
      double precision E_m, nu_m, D_m, x_br, g_av, gamma_r, f_rz
      double precision K2, nu_c, v, nu_min, B_nu, f_rel_br
      double precision ge, gp, vdsigma, j_pa, sum1_pa, eps
      double precision nu_p, gamma_bar, McDonald
      double precision j_sy_a(num_a), P_a(num_a,n_vol), P_sum_a(num_a), &
     &                 a_eb, sina
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i2, i3
      double precision F_sync, &
     &       gamma0(num_nt), nu_b, face, Ub, facg, tt, &
     &       eq43, eq13, expk43, expk13, ff, es, gamp(num_nt)
      real(8), parameter ::  sigmaT = 6.6524616d-25
      real(8), parameter ::  hplanck= 6.626075d-27
      real(8), parameter ::  emc2   = 8.187111059d-7
      real(8), parameter ::  ee     = 4.803d-10
      real(8), parameter ::  em     = 9.109d-28
!
!
      call get_j_k(zone, j, k, nr)
      if ((j.eq.1).and.(k.eq.1)) then
         l_min = dmin1(z(1), (r(1) - rmin))
      else if (j.eq.1) then
         l_min = dmin1(z(1), (r(k) - r(k-1)))
      else if (k.eq.1) then
         l_min = dmin1((z(j) - z(j-1)), (r(1) - rmin))
      else
         l_min = dmin1((z(j) - z(j-1)), (r(k) - r(k-1)))
      endif
      T_keV = tea(j,k)
      ne_local = n_e(j,k)
      B = B_field(j,k)
!
      n_harmonics = 5
!
      nu_b = ee*B/(2*pi*em*c_light)
      Ub = B**2.d0/(8.d0*pi)
      face = 3.d0**1.5d0*sigmaT*c_light*Ub/(pi*nu_b)
      !    = 2*dsqrt(3)*B*ee**3/(em*c_light**2)
      ! gamma*p (p==momentum), needed for the self-abs. frequency
      gamma0(:) = gnt(:)+1.d0
      gamp(:) = gamma0(:)*sqrt(gamma0(:)**2.d0-1.d0)

      dE = dexp(dlog(1.d20)/dble(n_vol))
      T = 1.16d7*T_keV
      Theta = T_keV/5.11d2
      kappa_C = 6.65d-25*ne_local
      if (Theta.lt.2.d-1) then
         K2 = 1.2533d0*dsqrt(Theta)*(1. + 1.875d0*Theta  &
     &       + 8.2031d-1*(Theta**2) - 2.03d-1*(Theta**3))/dexp(Theta)
!      else
!         K2 = 2.*(Theta**2)
!      endif
      else
         K2 = McDonald(2.d0, (1./Theta))
      endif
      nu_c = 2.8d6*B
      nu_min = dble(n_harmonics)*nu_c
      nu_p = 9.d3*dsqrt(ne_local)
!
      if (Theta.gt.0.1) then
         f_rel_br = 1.41d0*dsqrt(Theta)*(dlog(2.d0*Theta) + 9.228d-1) &
     &            - 1.d0
         f_rel_br = 1.d0 + (Theta**2)*f_rel_br/(1.d0 + (Theta**2))
         if (f_rel_br.lt.1.d0) f_rel_br = 1.d0
      else
         f_rel_br = 1.d0
      endif
!
      dg1 = 1.05d0
      p_1 = 1. - p_nth(j,k)
      if ((amxwl(j,k).gt.9.9999d-1).or. &
     &    (gmax(j,k).lt.1.01*gmin(j,k))) then
         N_e_nt = 0.d0
      else if (dabs(p_1).gt.1.d-3) then
         N_e_nt = ne_local*(1. - amxwl(j,k))*p_1 &
     &           /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
      else
         N_e_nt = ne_local*(1. - amxwl(j,k))/log(gmax(j,k)/gmin(j,k))
      endif
!
      P_sum = 0.
      P_sum_a(:) = 0.
      sum_th = 0.
      Eloss_cy(j,k) = 0.
      Eloss_th(j,k) = 0.
!
      E = 1.d-10/dE

      x_br = 7.353d1*T_keV
!
      g_av = gamma_bar(Theta)
!
      gamma_R = 2.1d-3*sqrt(ne_local)/(B*dsqrt(g_av))
      y = gamma_R/g_av
      if (y.lt.1.d2) then
         f_rz = dexp(-y)
      else
         f_rz = 0.
      endif
!
      do 50 i = 1, n_vol
         E = E*dE
         E_ph(i) = E
         nu = 2.41487d17*E
!
!           Thermal bremsstrahlung
!
         y = E/T_keV
!
         if ((x_br.lt.1.d8).and.(x_br.gt.1.d-8).and. &
     &       (y.gt.1.d-8).and.(y.lt.1.d8)) then
            if (y.ge.1.d0) then
               if (x_br.gt.1.d0) then
                  G_ff = 9.772d-1/dsqrt(y)
               else
                  if (x_br*y.lt.1.d0) then
                     G_ff = 1.d0
                  else
                     G_ff = 3.464d0/dsqrt(x_br*y)
                  endif
               endif
            else
               if (x_br.lt.1.) then
                  if (y.gt.dsqrt(x_br)) then
                     G_ff = 1.d0
                  else
                     G_ff = 2.7566d0*dlog(1.5802d1*x_br/y)
                  endif
               else
                  G_ff = 5.513d-1*dlog(6.9298d0/y)
               endif
            endif
         else
            G_ff = 0.
         endif
!
         if (y.lt.150.) then
            j_br = G_ff*dexp(-y)*1.6416d-20*(ne_local**2)*amxwl(j,k) &
     &            *f_rel_br/dsqrt(T)
         else
            j_br = 0.
         endif
!
         kappa_br = 2.6d-44*(ne_local**2)*G_ff*amxwl(j,k) &
     &             *(1. - dexp(-y))*f_rel_br/((T**1.5)*(E**3))
!
!           Non-thermal synchrotron
!
!         if ((nu.le.nu_p).or.(amxwl(j,k).gt.9.999d-1)) then
         if (nu.le.nu_p)then
            j_sy = 0.
            kappa_sy = 0.
            goto 30
         endif

!cccccccc          c            c               c              c
!       pitch angle dependent synchrotron emissivity
!       i_a is the sequence for the angle between electron and magnetic field
        do i_a = 1,num_a
          sum_sd = 0.
          a_eb=pi/num_a*(i_a-0.5d0)
          sina=dsin(a_eb)
          do i2=1,num_nt-1
            xgg = (2*pi)*nu*em*c_light/(1.5d0*ee*B*sina)
            !xgg =5.75d10*E/(B*sina) !xgg/g**2 is the 'x' in R&L book.
            gamma0(i2) = gnt(i2) + 1.d0
            y = xgg/(gamma0(i2)*gamma0(i2))
            if(y.le.0.001d0)then
              F_sync = 2.7083d0*(y/2)**0.3333333
            else if(y.ge.10.d0)then
              F_sync = 1.2533d0*dexp(-y)*dsqrt(y)
            else
                do i3=1,36
                  if(y.lt.F_sync_x(i3)) exit
                enddo
                F_sync =(y-F_sync_x(i3-1))/(F_sync_x(i3)-F_sync_x(i3-1)) &
     &                 *(F_sync_t(i3)-F_sync_t(i3-1))+F_sync_t(i3-1)
            endif
            sd = f_nt(j,k,i2)*F_sync
            sum_sd = sum_sd + (gnt(i2+1)-gnt(i2))*sd
          enddo
          j_sy_a(i_a) = 0.5d0*face*sum_sd*ne_local/(4.0*pi)*sina
          !j_sy_a(i_a) = 9.01d-6*sum_sd*ne_local
        enddo

!          j_sy_a(i_a) = sum_sd*ne_local/(4.d0*pi)*sina
!        integration by a from 0 to pi is: sina*2*pi*sina*d(a)=2*pi*d(a)*sina**2, then average by 4pi
!          j_sy = j_sy + sum_sd*ne_local/(4.d0*pi)*2*pi*(pi/num_a)*sina**2/(4.d0*pi)=
!          j_sy = j_sy + sum_sd*ne_local/(8.d0*num_a)*sina**2
!          kappa_sy = kappa_sy + sum_k*ne_local/(8.d0*pi*em*nu**2)*2*pi*(pi/num_a)*sina**2/(4.d0*pi)=
!          kappa_sy = kappa_sy + 
!     1               sum_k*ne_local/(16.d0*em*nu**2)/num_a*sina**2
!cccccccc          c            c               c              c
!         pitch angle averaged synchrotron emissivity
          sum_sd = 0.
          sum_k = 0.
          do i2=1,num_nt-1
              facg = 3.d0*(gamma0(i2)**2.d0)*nu_b
              tt = nu/facg
              if( tt < 1.0d4 ) then
                 eq43 = expk43(tt)
                 eq13 = expk13(tt)
                 ff = tt*tt*(eq43*eq13 - 0.6*tt*(eq43-eq13)*(eq43+eq13))
                 es = face*ff*dexp(-2.d0*tt)
              else
                 es = 0.0
              end if

            sd = f_nt(j,k,i2)*es
            sd_k = gamp(i2)*es
            sum_sd = sum_sd + (gnt(i2+1)-gnt(i2))*sd
            sum_k = sum_k +(f_nt(j,k,i2)/gamp(i2)-f_nt(j,k,i2+1)/ &
     &              gamp(i2+1))*sd_k
          enddo
          j_sy = sum_sd*ne_local/(4.d0*pi)
          kappa_sy = sum_k*ne_local/(8.d0*pi*em*nu**2)
!
!
!
!           Thermal cyclotron
!
 30      f_m = 1.d0
         j_cy = 0.
         kappa_cy = 0.
!
!        Sum first n harmonics
!
         if (nu.le.nu_p) then
            j_cy = 0.d0
            kappa_cy = 0.d0
            goto 26
         endif
!
         do 25 m = 1, n_harmonics
            mm = dble(m)
            f_m = f_m/(4.*mm)
            nu_m = mm*nu_c
            E_m = 4.14d-18*nu_m
            D_m = 7.07d-1*Theta*E_m
            x = ((E - E_m)/D_m)**2
            y = E_m/T_keV
!
            if (x.lt.50.) then
               f_cy = f_rz*dexp(-x)*ne_local*(B**2)*(Theta**(mm-1.5d0)) &
     &               *(mm + 1.d0)*f_m*(mm**(2.d0*mm + 1.d0))
               j_cy = j_cy + 8.46d-14*f_cy*(E**2)/(E_m**3)
               if (y.lt.150.) then
                  kappa_cy = kappa_cy + 5.705d33*(dexp(y) - 1.d0) &
     &                                 *f_cy/(nu*(nu_m**3))
               else
                  x = y - x
                  if (x.gt.150.) then 
                     kappa_cy = 1.d70
                  else if (x.gt.-100.) then
                     kappa_cy = kappa_cy+f_rz*5.705d33*dexp(x)*ne_local &
     &                         *(B**2)*(Theta**(mm - 1.5d0))*f_m &
     &                         *(mm + 1.d0)*(mm**(2.d0*mm + 1.d0)) &
     &                         /(nu*(nu_m**3))
                  endif
               endif
            endif
  25     continue
!
!        Use Mahadevan, Narayan & Yi (1996) formula
!          for the higher harmonics (nu > n*nu_c)
!
         if (nu.gt.nu_min) then
            v = nu/(nu_c*(Theta**2))
            y = 4.5*v
            if (y.lt.1.d6) then
               j_cy = j_cy + 4.652d-12*ne_local*nu/(K2*(v**1.6666667d-1) &
     &                      *dexp(y**3.33333d-1))
            endif
            y = E/T_keV
            if (y.gt.100.) then
               B_nu = 1.d-50
            else if (y.lt.1.d-6) then
               B_nu = 3.56d-30*(nu**3)/y
            else
               B_nu = 3.56d-30*(nu**3)/(dexp(E/T_keV) - 1.d0)
            endif
            if (y.lt.100.) kappa_cy = kappa_cy + j_cy/B_nu
         endif
  26     continue
!
!
!             Pair annihilation radiation
!
!
         j_pa = 0.d0
         if ((pair_switch.eq.0).or.(f_pair(j,k).lt.1.d-10)) goto 47
         eps = 1.957d-3*E
!
         do 45 i_el = 1, num_nt-1
            if (f_nt(j, k, i_el).lt.1.d-10) goto 45
            ge = gnt(i_el) + 1.d0
            sum1_pa = 0.d0
            do 40 i_pos = 1, num_nt-1
               if (n_pos(j, k, i_pos).lt.1.d-10) goto 40
               gp = gnt(i_pos) + 1.d0
!
               sum1_pa = sum1_pa + (gnt(i_pos+1) - gnt(i_pos)) &
     &                      *n_pos(j, k, i_pos)*vdsigma(eps, ge, gp)
  40        continue
            j_pa = j_pa + (gnt(i_el+1) - gnt(i_el))  &
     &                   *ne_local*f_nt(j, k, i_el)*sum1_pa
  45     continue
         j_pa = j_pa*eps*1.6d-9
!
!
!  47     kappa_tot(i, j, k) = kappa_br + kappa_sy + kappa_cy
  47     kappa_tot(i, j, k) = kappa_sy
!        Xuhui 6/1/09 deactivated br and cy absorption

!
!         if (kappa_tot(i, j, k).lt.(1.d1*kappa_C)) then
         if (kappa_tot(i, j, k).lt.dmax1(1.d0/l_min, &
     &       1d1*kappa_C))then ! Xuhui 5/8/09
!            P_sum = P_sum + (j_br + j_sy + j_cy + j_pa)*E*(dE - 1.d0)
            P_sum_a(:) = P_sum_a(:) + j_sy_a(:)*E*(dE - 1.d0)
            P_sum = P_sum + j_sy*E*(dE - 1.d0)
!         Xuhui 6/1/09; deactivated any spectrum except synchrotron
!            Eloss_cy(j, k) = Eloss_cy(j, k) + j_cy*E*(dE - 1.d0)
         else
            x = E/T_keV
            tau_tot = kappa_tot(i, j, k)*l_min
            if (x.lt.1.d2) then
               j_th = 1.47d-47*(nu**3)/(dexp(x) - 1.d0)
            else
               j_th = 1.d-50
            endif
            if (tau_tot.lt.5.d1)  &
     &         j_th = j_th*(1.d0 - dexp(-tau_tot))
            sum_th = sum_th + j_th*E*(dE - 1.d0)
!            Eloss_th(j, k) = Eloss_th(j, k) + j_th*E*(dE - 1.d0)
         endif
         P_a(:,i) = P_sum_a(:)
         P(i) = P_sum
         P_th(i) = sum_th
!
  50  continue

!
      do i_a=1,num_a
         if (P_sum_a(i_a).gt.1.d-50) then
           eps_a(i_a, :, j, k) = P_a(i_a,:)/P_sum_a(i_a)
         else
!           write(*,*)'P_sum_a =',P_sum_a(i_a),'angle bin=',i_a
           eps_a(i_a, :, j, k) = 0.
         endif
      enddo

      if (P_sum.gt.1.d-50) then
         eps_tot(:, j, k) = P(:)/P_sum
      else
         eps_tot(:, j, k) = 0.
      endif

      if (sum_th.gt.1.d-50) then
         eps_th(:, j, k) = P_th(:)/sum_th
      else
         eps_th(:, j, k) = 0.
      endif
!
!
!
      return
      end
!
!
!        This subroutine calculates the gamma-gamma
!        pair production opacity at photon energy E
!
!
      double precision function kgg_calc(E, j, k)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      double precision E
      integer j, k
!

      integer i
      double precision sum_mu, mu_thr, mu_local, dmu, eps1, eps2, beta,f
!
      kgg_calc = 0.d0
      eps2 = 1.957d-3*E
!
      do 200 i = 1, n_gg-1
         eps1 = 1.957d-3*E_gg(i)
         if ((n_ph(i, j, k).gt.1.d-20).and.(eps1.gt.1.d0/eps2)) then
            sum_mu = 0.
            mu_thr = 1.d0 - 2.d0/(eps1*eps2)
            if (mu_thr.gt.1.d0) mu_thr = 1.d0
            dmu = .01*(1. + mu_thr)
            mu_local = 5.d-1*dmu - 1.d0
!
 100        continue
            beta = dsqrt(1.d0 - 2.d0/(eps1*eps2*(1.d0 - mu_local)))
            if ((beta.ge.1.d0).or.(beta.le.0.d0)) goto 150
            f = (1.d0 - beta**2)*((3.d0 - beta**4) &
     &          *dlog((1.d0 + beta)/(1.d0 - beta)) &
     &          - 2.d0*beta*(2.d0 - beta**2))
            sum_mu = sum_mu + (1.d0 - mu_local)*f*dmu
 150        mu_local = mu_local + dmu
            if (mu_local.lt.mu_thr) goto 100
!
            kgg_calc = kgg_calc  &
     &               + sum_mu*n_ph(i, j, k)*(E_gg(i+1) - E_gg(i))
         endif
 200  continue
!
      kgg_calc = kgg_calc*6.234d-26
      return
      end
!
!
!       This subroutine calculates the photon-energy-
!      differential cross section for pair annihilation
!
!
      double precision function vdsigma(eps, ge, gp)
      double precision eps, ge, gp
!
      double precision be, bp, gcm_u, gcm_l, f_vds
      double precision gcm2, gcmstar, gcmmax, eps_u, eps_l
!
      if ((ge.lt.1.000001d0).or.(gp.lt.1.0000001d0)) then
         vdsigma = 0.d0
         goto 900
      endif
      be = dsqrt(1.d0 - 1.d0/(ge**2)) + 1.d-10
      bp = dsqrt(1.d0 - 1.d0/(gp**2)) + 1.d-10
      eps_u = 5.d-1*(gp*(1.d0 + bp) + ge*(1.d0 + be))
      eps_l = 5.d-1*(gp*(1.d0 - bp) + ge*(1.d0 - be))
      if ((eps.le.eps_l).or.(eps.ge.eps_u)) then
         vdsigma = 0.d0
         goto 900
      endif
!
      gcm2 = 5.d-1*(1.d0 + ge*gp*(1.d0 - be*bp))
      if (gcm2.le.1.00001d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcm_l = dsqrt(gcm2)
!
      gcm2 = 5.d-1*(1.d0 + ge*gp*(1.d0 + be*bp))
      if (gcm2.le.1.d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcmmax = dsqrt(gcm2)
      gcm2 = eps*(ge + gp - eps)
      if (gcm2.le.1.00001d0) then
         vdsigma = 0.d0
         goto 900
      endif
      gcmstar = dsqrt(gcm2)
      gcm_u = dmin1(gcmstar, gcmmax)
!
      if (gcm_u.gt.(1.0001d0*gcm_l)) then
         vdsigma = 7.48d-15*(f_vds(gcm_u, ge, gp, eps)  &
     &                     - f_vds(gcm_l, ge, gp, eps)) &
     &                     /(be*bp*((ge*gp)**2))
      else
         vdsigma = 0.d0
      endif
!
 900  return
      end
!
!
!
      double precision function f_vds(gcm, ge, gp, eps)
      double precision gcm, ge, gp, eps
!
      double precision H_pa, D
!
      D = (ge + gp)**2 - 4.*(gcm**2)
      if (D.le.1.d-20) then
         f_vds = 0.d0
         goto 900
      endif
!
      f_vds = dsqrt(D) + H_pa(gcm, ge, gp, eps)  &
     &                 + H_pa(gcm, gp, ge, eps)
!
 900  return
      end
!
!
!
      double precision function H_pa(gcm, ge, gp, eps)
      double precision gcm, ge, gp, eps
!
      double precision c, d, u, I_pa, u2, gcms2, gstar
!
!
      c = (ge - eps)**2 - 1.d0
      d = ge*(gp + ge) + eps*(gp - ge)
!
      gcms2 = eps*(ge + gp - eps)
      if (gcms2.lt.1.00001d0) then
         H_pa = 0.d0
         goto 900
      endif
      gstar = dsqrt(gcms2)
      u2 = c*(gcm**2) + gcms2
      if (u2.lt.1.d-20) then
         H_pa = 0.d0
         goto 900
      endif
      u = dsqrt(u2)
!
      if (dabs(c).gt.1.d-8) then
         H_pa = (2.d0 + (1.d0 - gcms2)/c)*I_pa(c, gcm, gstar, u) &
     &        + (1.d0/gcm - gcm/c + 5.d-1*gcm*(2.d0*c - d)/gcms2)/u &
     &        + gcm*u/c
      else
         H_pa = (2.d0*(gcm**3)/3.d0 + 2.d0*gcm + 1.d0/gcm)/gstar &
     &        + 5.d-1*(2.d0*(gcm**3)/3.d0 - d*gcm)/(gstar**3)
      endif
!
 900  return
      end
!
!
!
      double precision function I_pa(c, gcm, gcmstar, u)
      double precision c, gcm, gcmstar, u
!
      if (c.ge.1.d-8) then
         I_pa = dlog(gcm*dsqrt(c) + u)/dsqrt(c)
      else if (c.le.-1.d-8) then
         I_pa = asin(gcm*dsqrt(-c)/gcmstar)/dsqrt(-c)
      else
         I_pa = 0.d0
      endif
!
 900  return
      end
!
!
!
      double precision function gamma_bar(Theta)
      double precision Theta
!
      double precision sum, sd, t, ts, s, dt, K2, K3
      double precision McDonald
!
      if (Theta.lt.0.2) then
         gamma_bar = (1. + 4.375*Theta + 7.383*(Theta**2)  &
     &              + 3.384*(Theta**3))/(1. + 1.875*Theta  &
     &              + .8203*(Theta**2)) - Theta
         goto 900
      endif
!
      K2 = McDonald(2.d0, (1.d0/Theta))
      K3 = McDonald(3.d0, (1.d0/Theta))
!
      gamma_bar = K3/K2 - Theta
!
 900  continue
      if (gamma_bar.lt.1.d0) gamma_bar = 1.d0
!
      return
      end
!
!
!
      double precision function McDonald(nu, z)
      double precision nu, z
!
      double precision sum, sd, t, ts, s, dt, d, y, a, &
     &                 GammaF
!
      sum = 0.d0
      t = 1.d0
      dt = 1.001d0
      d = dt - 1.d0
      s = 5.d-1*(1.d0 + dt)
      a = nu - 5.d-1
!
 100  ts = t*s
      y = z*ts
      if (y.lt.2.25d2) then
         sd = ((ts*ts - 1.d0)**a)/dexp(y)
      else
         sd = 0.d0
      endif
      sum = sum + d*t*sd
      t = t*dt
      if ((t.lt.2.d0).or.(sd.gt.1.d-8)) goto 100
!
      McDonald = dsqrt(3.14159265d0)*((5.d-1*z)**nu)*sum &
     &          /GammaF(5.d-1 + nu)
!
      return
      end
!
!
!
!     This calculates the gamma function of xx.  It uses
!     the function gammln from "Numerical Recipies".  It is
!     called by McDonald.
      Function GammaF(xx)
      implicit none
      double precision xx, GammaF, gammln
      GammaF = exp(gammln(xx))
      return
      end
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccc
!  This function is from "Numerical Recipies in FORTRAN",
!  second edition, page207.  It calculates the ln of the
!  gamma function of xx.

      FUNCTION gammln(xx)
      implicit none
      DOUBLE PRECISION gammln, xx
      INTEGER j
      DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
      SAVE cof, stp
      DATA cof, stp/76.18009172947146d0,-86.50532032941677d0, &
     &    24.01409824083091d0,-1.231739572450155d0,  &
     &    .1208650973866179d-2,-.5395239384953d-5, &
     &    2.5066282746310005d0/
       x=xx
       y=x
       tmp=x+5.5d0
       tmp=(x+0.5d0)*log(tmp)-tmp
       ser=1.000000000190015d0
       do 11 j=1,6
          y=y+1.d0
          ser=ser+cof(j)/y
 11    continue
       gammln=tmp+log(stp*ser/x)
       return
       END
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        double precision function expk13(t)

!-------------------------------------------------------------------------------
!       Returns exp(t)*(Modified Bessel function order 1/3)(t)
!-------------------------------------------------------------------------------

        implicit none
        include 'general.pa'

        !! save c1,c2
        double precision  z3, zs, z, t, poly
        double precision  c1, c2
        double precision  f1, f2

        data c1,c2 / 0.35502805 , 0.25881940 /

        if( t .le. 1.0 ) then 
           !--------------------------------------------------
           ! Small argument use Airy function expansion A&S p446
           !
           z3 = 1.5d0*t
           zs = z3**0.3333333333333333d0
           z  = zs*zs
           z3 = z3*z3

           f1 = 1.0d0 + z3/6.0d0*(1.0d0 + z3/30.0d0*(1.0d0 + z3/56.0d0))
           f2 = z*(1.0d0 + z3/12.0d0*(1.0d0  &
     &          + z3/42.0d0*(1.0d0 + z3/90.0d0)))
           expk13 = exp(t)*pi*1.7320508d0/zs*(c1*f1 - c2*f2)

        else 
           !--------------------------------------------------
           ! Large arguments use asymptotic expansion
           !
           z   = 1.0d0/(72.0d0*t)
           poly= 1.0d0 - 5.0d0*z*(1.0d0 - 38.5d0*z)
           expk13 = dsqrt(0.5d0*pi/t)*poly &
     &              /(1.0d0 + 1.0d0/(1.0d0 + 58.0d0*t*t))

        end if 

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        double precision function expk43(t)

!-------------------------------------------------------------------------------
!      Returns exp(t)**(modified Bessel function order 4/3)(t)
!-------------------------------------------------------------------------------
        
        implicit none
        include 'general.pa'

        real(8) :: z, poly, t

        if( t .le. 1.0 ) then 
           !--------------------------------------------------
           ! Small argument use fit.
           !
           poly   = 1.0d0 + t*(0.9757317d0 - 7.6790616d-2*t)
           expk43 = 0.44648975d0*(2.0d0/t)**1.333333333d0*poly

        else 
           !--------------------------------------------------
           ! Large arguments use asymptotic argument.
           !
           z    = 1.0/(72.0*t)
           poly = 1.0+55.0*z*(1.0-8.5*z)
           expk43 = dsqrt(0.5*pi/t)*poly*(1.0 + 1.0/(1.0 + 50.0*t*t))

        end if

        return
        end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun  6 23:56:12 EDT 2006
! version: 2
! Name: J. Finke
! Changed a few double constants with 'e' to 
! 'd'. Changed variable 'm' to an integer from 
! double.        
