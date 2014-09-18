!
!    This subroutine calculates the pure Compton
!       reflection matrix according to White,
!    Lightman \& Zdziarski (1988; ApJ, 331, 939)
!
!
      subroutine Pref_calc()
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
      integer n_in, n_out
      double precision sum, x0, x1, y0, y1, dy, dyc, dE, Gx, Gy,&
     &                 A, B, alpha_pref, beta, dx, E_trans
!
!
      dE = dexp(dlog(1.d3)/dble(n_ref))
      E_ref(1) = 1.d0
      do 100 n_in = 2, n_ref
 100  E_ref(n_in) = E_ref(n_in-1)*dE
      E_trans = 20.
!
      do 500 n_in = 1, n_ref
         x0 = 1.957d-3*E_ref(n_in)
         y0 = 1.d0/x0
         sum = 0.
         do 450 n_out = 1, n_ref
            if (n_out.gt.n_in) then
               P_ref(n_out, n_in) = 1.d0
               goto 450
            endif
!
            if (E_ref(n_in).le.E_trans) then
               if (n_out.lt.n_in) then
                  P_ref(n_out, n_in) = 0.d0
               else
                  P_ref(n_out, n_in) = 1.d0
               endif
               goto 450
            endif
!
            x1 = 1.957d-3*E_ref(n_out)
            y1 = 1.d0/x1
            dy = y1 - y0
            dyc = 1.d3 - y0
            A = 5.6d-1 + 1.12/(y0**0.785) - 3.4d-1/(y0**1.04)
            alpha_pref = -3.d-1/(y0**.51) + 6.d-2/(y0**.824)
            beta = 3.7d-1 - (y0**.85)
            if (dabs(alpha_pref + .5).lt.1.d-4) then
               B = (1.d0 - A*(2.d0 + dlog(.5*dyc))/dsqrt(dyc))&
     &            /((y0**(1.d0 - beta))*((y0 + 2.d0)**beta)&
     &              *((1.d0 + 2.d0/y0)**(1.d0 - beta) - 1.d0))&
     &            *(1.d0 - beta)
            else
               B = (1.d0 - A*(2.d0 + (((.5*dyc)**(alpha_pref + .5))&
     &              - 1.d0)/(alpha_pref + .5))/dsqrt(dyc))&
     &            /((y0**(1.d0 - beta))*((y0 + 2.d0)**beta)&
     &              *((1.d0 + 2.d0/y0)**(1.d0 - beta) - 1.d0))&
     &            *(1.d0 - beta)
            endif
!
            if (dy.lt.2.d0) then
               Gy = B*(((y0 + 2.d0)/(y0 + dy))**beta)
            else if (dy.lt.dyc) then
               Gy = A*((dyc/dy)**alpha_pref)/(dy**1.5)
            else
               Gy = A/(dy**1.5)
            endif
            Gx = Gy/(x1**2)
            dx = dE*x1
            sum = sum + Gx*dx
!            if (sum.gt.1.d0) sum = 1.d0
            P_ref(n_out, n_in) = sum
!
 450     continue 
         do 480 n_out = 1, n_ref
            if ((E_ref(n_out).le.E_ref(n_in)).and.&
     &          (E_ref(n_in).gt.E_trans))&
     &         P_ref(n_out, n_in) = P_ref(n_out, n_in)/sum
 480     continue
 500  continue
!
!
      return
      end
!
!
!
!    This subroutine calculates the disk ionization
!     parameter, finds the respective ionization
!    state fractions and photoelectric absorption
!     opacities, and defines the absorption matrix
!          for Compton reflection calculations
!
!
      subroutine Wref_calc()
      implicit none
      include 'general.pa'
      include 'commonblock.f90'

!
      integer n_in, n_out, i, n
      double precision x0, x1, y, xi, ne_local, k_nu, kappa_c,eps(n_ref)
      double precision sigma_ions
      double precision ab0_he, ab0_c, ab0_n, ab0_o, ab0_ne, ab0_mg,&
     &                 ab0_si, ab0_s, ab0_ar, ab0_ca, ab0_fe, ab0_ni
      double precision ionf_he(2), ionf_c(6), ionf_n(7), ionf_o(8),&
     &                 ionf_ne(10), ionf_mg(12), ionf_si(14),&
     &                 ionf_s(16), ionf_ar(18), ionf_ca(20), &
     &                 ionf_fe(26), ionf_ni(28)
      double precision ab_he(2), ab_c(6), ab_n(7), ab_o(8),&
     &                 ab_ne(10), ab_mg(12), ab_si(14), ab_s(16),&
     &                 ab_ar(18), ab_ca(20), ab_fe(26), ab_ni(28)
      double precision sigma_he(2), sigma_c(6), sigma_n(7), sigma_o(8),&
     &                 sigma_ne(10), sigma_mg(12), sigma_si(14), &
     &                 sigma_s(16), sigma_ar(18), sigma_ca(20),&
     &                 sigma_fe(26), sigma_ni(28)
      double precision sigma_c2(6), sigma_n2(7), sigma_o2(8), &
     &                 sigma_ne2(10), sigma_mg2(12), sigma_si2(14), &
     &                 sigma_s2(16), sigma_ar2(18), sigma_ca2(20), &
     &                 sigma_fe2(26), sigma_ni2(28)
!
      double precision he_edge(2), c_edge(6), n_edge(7), o_edge(8), &
     &                 ne_edge(10), mg_edge(12), si_edge(14), &
     &                 s_edge(16), ar_edge(18), ca_edge(20), &
     &                 fe_edge(26), ni_edge(28)
      double precision c_edge2(6), n_edge2(7), o_edge2(8), ne_edge2(10), &
     &                 mg_edge2(12), si_edge2(14), s_edge2(16), &
     &                 ar_edge2(18), ca_edge2(20), fe_edge2(26), &
     &                 ni_edge2(28)
!
      data ionf_he / 1.d0, 0.d0 /
      data ionf_c  / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_n  / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_o  / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_ne / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0 /
      data ionf_mg / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_si / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_s  / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_ar / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0 /
      data ionf_ca / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0 /
      data ionf_fe / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0 /
      data ionf_ni / 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
     &               0.d0, 0.d0, 0.d0, 0.d0 /
!
      data ab0_he, ab0_c, ab0_n,  ab0_o,  ab0_ne, ab0_mg, &
     &     ab0_si, ab0_s, ab0_ar, ab0_ca, ab0_fe, ab0_ni /&
     &     6.33d-2, 3.90d-4, 8.12d-5, 6.47d-4, 9.14d-5, 3.73d-5,&
     &     3.52d-5, 1.76d-5, 3.73d-6, 2.20d-6, 3.16d-5, 1.68e-6 /

!
      data sigma_he / 9.0d-18, 1.4d-18 /
      data sigma_c / 1.d-18, 9.5d-19, 8.d-19, 6.d-19, 4.5d-19, &
     &               1.6d-19 /
      data sigma_n / 9.d-19, 7.0d-19, 6.d-19, 5.d-19, 4.5d-19, &
     &               3.d-19, 1.1d-19 /
      data sigma_o / 6.d-19, 5.0d-19, 4.d-19, 4.d-19, 3.1d-19, &
     &               3.d-19, 2.1d-19, 1.d-19 /
      data sigma_ne/ 4.d-19, 4.0d-19, 3.d-19, 3.d-19, 2.4d-19, &
     &              2.2d-19, 2.1d-19, 2.d-19,1.5d-19, 6.0d-20 /
      data sigma_mg/ 2.d-19, 2.0d-19, 2.d-19, 2.d-19, 1.9d-19, &
     &              1.5d-19, 1.5d-19,1.5d-19,1.2d-19, 1.1d-19,&
     &               1.d-19, 4.0d-20 /
      data sigma_si/1.3d-19, 1.3d-19,1.3d-19,1.3d-19, 1.3d-19,&
     &              1.2d-19, 1.2d-19,1.1d-19, 1.d-19, 1.0d-19,&
     &               9.d-20, 7.5d-20, 7.d-20, 3.d-20 /
      data sigma_s / 1.d-19, 1.0d-19, 1.d-19, 1.d-19, 1.d-19,&
     &               1.d-19, 1.0d-19, 1.d-19, 8.d-20, 8.d-20,&
     &               8.d-20, 7.0d-20, 6.d-20, 5.d-20, 5.d-20,&
     &              2.5d-20 /
      data sigma_ar/ 8.d-20, 8.0d-20, 8.d-20, 8.d-20, 8.d-20,&
     &               8.d-20, 7.5d-20, 7.d-20, 7.d-20, 7.d-20, &
     &               6.d-20, 6.0d-20, 5.d-20, 5.d-20, 5.d-20,&
     &               4.d-20, 4.0d-20, 2.d-20 /
      data sigma_ca/ 7.d-20, 6.5d-20, 6.d-20, 6.d-20, 6.d-20,&
     &               6.d-20, 6.0d-20, 6.d-20, 6.d-20, 5.d-20,&
     &               5.d-20, 5.0d-20, 5.d-20, 5.d-20, 4.d-20,&
     &               4.d-20, 4.0d-20,3.5d-20, 3.d-20,1.5d-20 /
      data sigma_fe/ 3.d-20, 3.0d-20,3.5d-20, 3.d-20, 3.d-20,&
     &               3.d-20, 3.0d-20, 3.d-20, 3.d-20, 3.d-20,&
     &               3.d-20, 2.7d-20, 3.d-20, 3.d-20, 3.d-20,&
     &               3.d-20, 3.0d-20, 3.d-20,2.5d-20,2.5d-20,&
     &              2.5d-20, 2.5d-20, 2.d-20, 2.d-20, 2.d-20,&
     &               1.d-20 /
      data sigma_ni/ 3.d-20, 3.0d-20, 3.d-20, 3.d-20, 3.d-20, &
     &               3.d-20, 2.8d-20, 3.d-20, 3.d-20, 3.d-20,&
     &               3.d-20, 2.5d-20,2.5d-20,2.5d-20,2.5d-20,&
     &              2.2d-20, 2.1d-20, 2.d-20, 2.d-20, 2.d-20,&
     &               2.d-20, 2.0d-20, 2.d-20, 2.d-20,1.7d-20, &
     &              1.6d-20, 1.6d-20, 8.d-21 /
!
      data he_edge/ .024, .052 /
      data c_edge / .30, .31, .31, .33, .40, .49 /
      data n_edge / .40, .41, .42, .48, .51, .56, .69 /
      data o_edge / .52, .57, .59, .61, .64, .68, .74, .88 /
      data ne_edge/ .88, .89, .91, .95, 1.0, 1.1, 1.1, 1.1, 1.2,&
     &              1.3 /
      data mg_edge/ 1.2, 1.3, 1.3, 1.4, 1.4, 1.4, 1.4, 1.5, 1.6, &
     &              1.7, 1.7, 2.0 /
      data si_edge/ 1.8, 1.8, 1.9, 1.9, 1.9, 2.0, 2.0, 2.1, 2.1,&
     &              2.1, 2.2, 2.3, 2.4, 2.6 /
      data s_edge / 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.5, 2.6, 2.8,&
     &              2.8, 2.9, 2.9, 3.1, 3.1, 3.2, 3.4 /
      data ar_edge/ 3.1, 3.1, 3.1, 3.1, 3.2, 3.2, 3.3, 3.3, 3.3,&
     &              3.3, 3.5, 3.5, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2 /
      data ca_edge/ 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1,&
     &              4.2, 4.2, 4.2, 4.5, 4.6, 4.7, 4.8, 4.8, 5.0,&
     &              5.1, 5.4 /
      data fe_edge/ 7.1, 7.1, 7.1, 7.3, 7.3, 7.4, 7.4, 7.4, 7.4, &
     &              7.4, 7.4, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9, 7.9,&
     &              8.1, 8.1, 8.2, 8.2, 8.5, 8.6, 8.6, 9.2 /
      data ni_edge/ 8.2, 8.2, 8.3, 8.3, 8.3, 8.4, 8.5, 8.5, 8.6, &
     &              8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 9.0, 9.0, 9.0, &
     &              9.0, 9.0, 9.5, 9.5, 9.7, 10., 10., 10., 10.,&
     &              11. /
!
      data sigma_c2 / 3.d-16, 2.d-17, 2.d-18, 5.d-19, 0.d0, 0.d0 /
      data sigma_n2 / 3.d-16, 6.d-16, 1.d-17, 9.d-19, 4.d-19, 0.d0, &
     &                0.d0 /
      data sigma_o2 / 5.d-16, 3.d-17, 7.d-18, 4.d-18, 1.d-18, 3.d-19,&
     &                0.d0, 0.d0 /
      data sigma_ne2/ 1.d-15, 7.d-17, 3.d-17, 5.d-18, 3.d-18, 2.d-18,&
     &                6.d-19, 2.d-19, 0.d0,   0.d0 /
      data sigma_mg2/ 4.d-17, 5.d-17, 2.d-17, 1.d-17, 5.d-18, 3.d-18, &
     &                2.d-18, 1.d-18, 3.d-19,1.2d-19, 0.d0,   0.d0 /
      data sigma_si2/ 4.d-17, 2.d-17, 1.d-17, 7.d-18, 5.d-18, 4.d-18,&
     &                3.d-18, 2.d-18, 1.d-18, 7.d-19, 2.d-19,1.5d-19,&
     &                  0.d0, 0.d0 /
      data sigma_s2 / 1.d-17, 7.d-18, 6.d-18, 5.d-18, 4.d-18, 4.d-18,&
     &                2.d-18, 1.d-18, 1.d-18, 1.d-18, 6.d-19, 4.d-19, &
     &                3.d-19,1.3d-19,   0.d0, 0.d0 /
      data sigma_ar2/ 7.d-18, 6.d-18, 5.d-18, 5.d-18, 4.d-18, 3.d-18,&
     &                3.d-18, 2.d-18, 2.d-18, 1.d-18, 1.d-18, 8.d-19,&
     &                4.d-19, 3.d-19, 1.d-19, 5.d-20,   0.d0,  0.d0 /
      data sigma_ca2/ 4.d-18, 4.d-18, 4.d-18, 3.d-18, 3.d-18, 3.d-18,&
     &                2.d-18, 2.d-18, 2.d-18,1.5d-18,1.5d-18, 1.d-18,&
     &                6.d-19, 4.d-19, 3.d-19, 2.d-19, 1.d-19, 5.d-20,&
     &                  0.d0,  0.d0 /
      data sigma_fe2/2.5d-18, 2.d-18, 2.d-18, 2.d-18, 2.d-18,1.5d-18,&
     &               1.5d-18,1.2d-18, 1.d-18, 1.d-18, 9.d-19, 7.d-19,&
     &                7.d-19, 7.d-19, 6.d-19, 5.d-19, 4.d-19, 4.d-19,&
     &                3.d-19, 2.d-19,1.5d-19, 1.d-19, 5.d-20, 3.d-20, &
     &                  0.d0,  0.d0 /
      data sigma_ni2/ 2.d-18, 2.d-18,1.5d-18,1.5d-18,1.5d-18, 1.d-18,&
     &                1.d-18, 1.d-18, 9.d-19, 9.d-19, 7.d-19, 7.d-19,&
     &                6.d-19, 6.d-19, 5.d-19, 4.d-19, 3.d-19, 3.d-19,&
     &               2.5d-19,2.5d-19, 2.d-19,1.5d-19, 1.d-19, 7.d-20,&
     &                5.d-20, 2.d-20,   0.d0, 0.d0 /
!
      data c_edge2 / .011, .024,.049, .065, 1.d6, 1.d6 /
      data n_edge2 / .014, .030,.046, .075, .098, 1.d6, 1.d6 /
      data o_edge2 / .013, .033,.053, .075, .11, .13, 1.d6, 1.d6 /
      data ne_edge2/ .021, .040,.062, .097, .14, .16, .21, .24, &
     &               1.d6,1.d6 /
      data mg_edge2/.054,.064,.080, .12, .13, .18, .22, .27, .32,&
     &               .35,1.d6,1.d6 /
      data si_edge2/ .11, .11, .12, .13, .15, .20, .23, .30, .33,&
     &               .38, .49, .52,1.d6,1.d6 /
      data s_edge2 / .16, .18, .20, .22, .24, .24, .28, .32, .38,&
     &               .42, .51, .55, .64, .71,1.d6,1.d6 /
      data ar_edge2/ .23, .25, .28, .30, .31, .33, .35, .38, .41, &
     &               .49, .53, .62, .69, .74, .88, .91,1.d6,1.d6 /
      data ca_edge2/ .35, .35, .38, .40, .41, .42, .45, .51, .52,&
     &               .55, .59, .63, .71, .82, .91, .93, 1.1, 1.1,&
     &              1.d6, 1.d6 /
      data fe_edge2/ .71, .73, .74, .74, .79, .82, .88, .88, .91,&
     &               .95, 1.0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.1, 1.2,&
     &               1.3, 1.5, 1.7, 1.8, 2.0, 2.1,1.d6, 1.d6 /
      data ni_edge2/ .89, .89, .91, .94, .94, 1.0, 1.0, 1.0, 1.0,&
     &               1.0, 1.1, 1.1, 1.2, 1.2, 1.3, 1.3, 1.4, 1.4,&
     &               1.5, 1.6, 1.6, 1.8, 1.9, 2.1, 2.2, 2.3,1.d6,&
     &               1.d6 /
!
!
!      xi = F_ib/n_disk
!
      xi = 0.d0
      n_disk = 1.d18
      ne_local = n_disk
      kappa_c = 6.65d-25*ne_local
!
!
!     Algorithm or table for the ionization
!      fractions needs to be inserted here!
!
!
      do 10 i = 1, 2
 10   ab_he(i) = ab0_he*ionf_he(i)
!
      do 15 i = 1, 6
 15   ab_c(i) = ab0_c*ionf_c(i)
!
      do 20 i = 1, 7
 20   ab_n(i) = ab0_n*ionf_n(i)
!
      do 25 i = 1, 8
 25   ab_o(i) = ab0_o*ionf_o(i)
!
      do 30 i = 1, 10
 30   ab_ne(i) = ab0_ne*ionf_ne(i)
!
      do 35 i = 1, 12
 35   ab_mg(i) = ab0_mg*ionf_mg(i)
!
      do 40 i = 1, 14
 40   ab_si(i) = ab0_si*ionf_si(i)
!
      do 45 i = 1, 16
 45   ab_s(i) = ab0_s*ionf_s(i)
!
      do 50 i = 1, 18
 50   ab_ar(i) = ab0_ar*ionf_ar(i)
!
      do 55 i = 1, 20
 55   ab_ca(i) = ab0_ca*ionf_ca(i)
!
      do 60 i = 1, 26
 60   ab_fe(i) = ab0_fe*ionf_fe(i)
!
      do 65 i = 1, 28
 65   ab_ni(i) = ab0_ni*ionf_ni(i)
!
!
      do 400 n = 1, n_ref
         sigma_ions = 0.
!
         do 112 i = 1, 2
            if (ab_he(i).gt.0.) then
               sigma_ions = sigma_ions + sigma_he(i)*ab_he(i) &
     &                                  /((E_ref(n)/he_edge(i))**3)
            endif
  112    continue
!
         do 113 i = 1, 6
            if ((ab_c(i).gt.0.).and.(E_ref(n).gt.c_edge(i))) then
               sigma_ions = sigma_ions + sigma_c(i)*ab_c(i) &
     &                                  /((E_ref(n)/c_edge(i))**3)
             endif
             if ((ab_c(i).gt.0.).and.(E_ref(n).gt.c_edge2(i))) then
               sigma_ions = sigma_ions + sigma_c2(i)*ab_c(i) &
     &                                  /((E_ref(n)/c_edge2(i))**3)
            endif
  113    continue
!
         do 114 i = 1, 7
            if ((ab_n(i).gt.0.).and.(E_ref(n).gt.n_edge(i))) then
               sigma_ions = sigma_ions + sigma_n(i)*ab_n(i) &
     &                                  /((E_ref(n)/n_edge(i))**3)
            endif
            if ((ab_n(i).gt.0.).and.(E_ref(n).gt.n_edge2(i))) then
               sigma_ions = sigma_ions + sigma_n2(i)*ab_n(i) &
     &                                  /((E_ref(n)/n_edge2(i))**3)
            endif
  114    continue
!
         do 115 i = 1, 8
            if ((ab_o(i).gt.0.).and.(E_ref(n).gt.o_edge(i))) then
               sigma_ions = sigma_ions + sigma_o(i)*ab_o(i) &
     &                                  /((E_ref(n)/o_edge(i))**3)
            endif
            if ((ab_o(i).gt.0.).and.(E_ref(n).gt.o_edge2(i))) then
               sigma_ions = sigma_ions + sigma_o2(i)*ab_o(i) &
     &                                  /((E_ref(n)/o_edge2(i))**3)
            endif
  115    continue
!
         do 116 i = 1, 10
            if ((ab_ne(i).gt.0.).and.(E_ref(n).gt.ne_edge(i))) then
               sigma_ions = sigma_ions + sigma_ne(i)*ab_ne(i) &
     &                                  /((E_ref(n)/ne_edge(i))**3)
            endif
            if ((ab_ne(i).gt.0.).and.(E_ref(n).gt.ne_edge2(i))) then
               sigma_ions = sigma_ions + sigma_ne2(i)*ab_ne(i) &
     &                                  /((E_ref(n)/ne_edge2(i))**3)
           endif
  116    continue
!
         do 117 i = 1, 12
            if ((ab_mg(i).gt.0.).and.(E_ref(n).gt.mg_edge(i))) then
               sigma_ions = sigma_ions + sigma_mg(i)*ab_mg(i) &
     &                                  /((E_ref(n)/mg_edge(i))**3)
            endif
            if ((ab_mg(i).gt.0.).and.(E_ref(n).gt.mg_edge2(i))) then
               sigma_ions = sigma_ions + sigma_mg2(i)*ab_mg(i) &
     &                                  /((E_ref(n)/mg_edge2(i))**3)
            endif
  117    continue
!
         do 118 i = 1, 14
            if ((ab_si(i).gt.0.).and.(E_ref(n).gt.si_edge(i))) then
               sigma_ions = sigma_ions + sigma_si(i)*ab_si(i) &
     &                                  /((E_ref(n)/si_edge(i))**3)
            endif
            if ((ab_si(i).gt.0.).and.(E_ref(n).gt.si_edge2(i))) then
               sigma_ions = sigma_ions + sigma_si2(i)*ab_si(i) &
     &                                 /((E_ref(n)/si_edge2(i))**3)
            endif
  118    continue
!
         do 119 i = 1, 16
            if ((ab_s(i).gt.0.).and.(E_ref(n).gt.s_edge(i))) then
               sigma_ions = sigma_ions + sigma_s(i)*ab_s(i) &
     &                                  /((E_ref(n)/s_edge(i))**3)
            endif
            if ((ab_s(i).gt.0.).and.(E_ref(n).gt.s_edge2(i))) then
               sigma_ions = sigma_ions + sigma_s2(i)*ab_s(i) &
     &                                  /((E_ref(n)/s_edge2(i))**3)
            endif
  119    continue
!
         do 120 i = 1, 18
            if ((ab_ar(i).gt.0.).and.(E_ref(n).gt.ar_edge(i))) then
               sigma_ions = sigma_ions + sigma_ar(i)*ab_ar(i) &
     &                               /((E_ref(n)/ar_edge(i))**3)
            endif
            if ((ab_ar(i).gt.0.).and.(E_ref(n).gt.ar_edge2(i))) then
               sigma_ions = sigma_ions + sigma_ar2(i)*ab_ar(i) &
     &                                  /((E_ref(n)/ar_edge2(i))**3)
            endif
  120    continue
!
         do 121 i = 1, 20
            if ((ab_ca(i).gt.0.).and.(E_ref(n).gt.ca_edge(i))) then
               sigma_ions = sigma_ions + sigma_ca(i)*ab_ca(i) &
     &                                  /((E_ref(n)/ca_edge(i))**3)
            endif
            if ((ab_ca(i).gt.0.).and.(E_ref(n).gt.ca_edge2(i))) then
               sigma_ions = sigma_ions + sigma_ca2(i)*ab_ca(i) &
     &                                  /((E_ref(n)/ca_edge2(i))**3)
            endif
  121    continue
!
         do 122 i = 1, 26
            if ((ab_fe(i).gt.0.).and.(E_ref(n).gt.fe_edge(i))) then
               sigma_ions = sigma_ions + sigma_fe(i)*ab_fe(i) &
     &                                  /((E_ref(n)/fe_edge(i))**3)
            endif
            if ((ab_fe(i).gt.0.).and.(E_ref(n).gt.fe_edge2(i))) then
               sigma_ions = sigma_ions + sigma_fe2(i)*ab_fe(i) &
     &                                  /((E_ref(n)/fe_edge2(i))**3)
            endif
  122    continue
!
         do 123 i = 1, 28
            if ((ab_ni(i).gt.0.).and.(E_ref(n).gt.ni_edge(i))) then
               sigma_ions = sigma_ions + sigma_ni(i)*ab_ni(i) &
     &                                   /((E_ref(n)/ni_edge(i))**3)
            endif
            if ((ab_ni(i).gt.0.).and.(E_ref(n).gt.ni_edge2(i))) then
               sigma_ions = sigma_ions + sigma_ni2(i)*ab_ni(i) &
     &                                   /((E_ref(n)/ni_edge2(i))**3)
            endif
  123    continue
!
         k_nu = sigma_ions*n_disk
         eps(n) = k_nu/(k_nu + kappa_c)
!
  400 continue
!
      do 500 n_in = 1, n_ref
         x0 = 1.957d-3*E_ref(n_in)
         do 450 n_out = 1, n_ref
            if (n_out.gt.n_in) then
               W_abs(n_out, n_in) = 0.d0
               goto 450
            endif
            x1 = 1.957d-3*E_ref(n_out)
!           
            if (E_ref(n_in).gt.20.) then
               y = 2.5d-6*(1.d0/(x0**4) - 1.d0/(x1**4))
               if (y.ge.-50.) then
                  W_abs(n_out, n_in) = dmin1(1.d0, dexp(y))
               else
                  W_abs(n_out, n_in) = 0.d0
               endif
            else
               W_abs(n_out, n_in) = (1.d0 - dsqrt(eps(n_in)))&
     &                             /(1.d0 + dsqrt(eps(n_in)))
            endif
!
  450    continue
  500 continue
!
      return
      end
!
