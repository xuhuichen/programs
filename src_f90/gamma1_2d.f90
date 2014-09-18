!
      subroutine gam_min(j, k, e_temp)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision e_temp
!
      double precision Ntherm, Nnth, Theta, p_1, fth, fnth, beta
      double precision p_2
      double precision Inttherm, gamma_bar
!
      Theta = e_temp/511.
      if (amxwl(j,k).gt.0.99999999d0) then
         gmin(j,k) = 1.d0 + 1.d2*Theta
         N_nth(j,k) = 1.d-20         
         goto 15
      endif
      p_1 = 1.d0 - p_nth(j,k)
      p_2 = 2.d0 - p_nth(j,k)
      if (dabs(p_1).lt.1.d-4) p_1 = p_1 + 2.d-4
      if (dabs(p_2).lt.1.d-4) p_2 = p_2 + 2.d-4

      if (amxwl(j,k).lt.1.d-4 .or. e_temp.lt.1.d0) then
        if (gmin(j,k).lt.1.01d0) gmin(j,k) = 1.01d0
        N_nth(j,k) = p_1/(gmax(j,k)**p_1 - gmin(j,k)**p_1)
        goto 15
      endif

!     gmin(j,k) = dmin1(.9*gmax(j,k), 1.0d1)

  10  gmin(j,k) = gmin(j,k)/1.01
!      if (gmin(j,k).le.1.01) goto 15 ! Xuhui 05/25/09
      beta = sqrt(1. - 1./(gmin(j,k)*gmin(j,k)))
      N_nth(j,k) = (1. - amxwl(j,k))*p_1&
     &            /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
      Ntherm = amxwl(j,k)/Inttherm(gmin(j,k), Theta)
      fth = Ntherm*beta*gmin(j,k)*gmin(j,k)*exp(-gmin(j,k)/Theta)
      fnth = N_nth(j,k)/(gmin(j,k)**p_nth(j,k))

  15  continue
  
      if (amxwl(j,k).gt.1.d-4) then
         gbar_nth(j,k) = gamma_bar(Theta)
      else
         gbar_nth(j,k) = N_nth(j,k)&
     &                 *(gmax(j,k)**p_2 - gmin(j,k)**p_2)/p_2
      endif

      return
      end
!
!
      double precision function Inttherm(g1, Theta)
      double precision g1, Theta
      double precision t, sum, dt, sd, sw, s, y
!
      t = 1.
      sum = 0.
      dt = 1.005
      s = .5*(1. + dt)
      sw = dt - 1.
!
  20  ts = t*s
      y = ts/Theta
      if (y.lt.2.d2) then
         sd = ts*sqrt(ts*ts - 1.)*exp(-y)
      else
         sd = 1.d-200
      endif
      sum = sum + sd*t*sw
      t = t*dt
      if (t.lt.g1) goto 20
!
      Inttherm = sum
      return
      end
!
!
!
      double precision function Ith_new(g1, Theta)
      double precision g1, Theta
      double precision t, sum, dt, sd, sw, s, y
!
      t = 1.
      sum = 0.
      dt = 1.005
      s = .5*(1. + dt)
      sw = dt - 1.
!
  20  ts = t*s
      y = (ts - 1.)/Theta
      if (y.lt.2.d2) then
         sd = ts*sqrt(ts*ts - 1.)*exp(-y)
      else
         sd = 1.d-200
      endif
      sum = sum + sd*t*sw
      t = t*dt
      if (t.lt.g1) goto 20
!
      Ith_new = sum
      return
      end
!
!
      double precision function K_two(x)
      double precision x
      double precision t, sum, dt, sd, s0, sw

      t = 1.0001
      sum = 0.
      dt = 1.005
      s0 = 0.
      sw = dt - 1.

  30  sd = ((t*t - 1.)**1.5)*exp(-x*t)
      sum = sum + sd*t*sw
      t = t*dt
      if (sd.ge.s0) s0 = sd
      if (abs(sd/s0).gt.1.e-8) goto 30

      K_two = .3333*x*x*sum
      return
      end



