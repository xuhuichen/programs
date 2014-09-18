!
      subroutine P_nontherm(j, k)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
!
!
      integer nunit_fnt
      integer i, wr_switch
      double precision x, y, sum, g, Nnth, Ntherm, Theta, p_1, &
     &                 beta, Ith, n_eplus, f_pos(num_nt)
      double precision Inttherm, Ith_new
      double precision g_1, Delta_g, g_trans
!
!
      character*20 pho_file, spec_file, lc_file, es_file
      common/out_names/pho_file, spec_file, lc_file, es_file
!
      Theta = tea(j,k)/5.11d2
      p_1 = 1. - p_nth(j,k)
      if (amxwl(j,k).lt.1.d-4) then
         Ntherm = 0.d0
      else
         Ith = Ith_new(gmin(j,k), Theta)
         Ntherm = amxwl(j,k)/Ith
      endif
      if (amxwl(j,k).lt.0.99999999) then
         Nnth = (1. - amxwl(j,k))*p_1/(gmax(j,k)**p_1 - gmin(j,k)**p_1)
      else
         Nnth = 0.d0
      endif
!
!       write(*,*) 'Ntherm = ',Ntherm,'; Nnth = ',Nnth
!      if (es_freq*((ncycle-1)/es_freq).eq.(ncycle-1)) then
!         wr_switch = 1
!      else
!         wr_switch = 0
!      endif
!
      wr_switch = 1
      es_file = 'output/seb.dat'
!
      nunit_fnt=22
      if (wr_switch.eq.1) &
     &    open(unit=nunit_fnt, file=es_file, status='unknown')
!
      if (ncycle.gt.0) goto 300
!
      sum = 0.d0
      dg = 1.1d0
      g_1 = 2.d-1/dg

!
      do 100 i = 1, num_nt
         g = 1.d0 + g_1
         beta = sqrt(1. - 1./(g**2))
         if (g.lt.gmin(j,k)) then
           if(amxwl(j,k).gt.1.d-4)then
             y = (g - 1.d0)/Theta
             if (y.lt.1.d2) then
                f_nt(j, k, i) = Ntherm*(g**2)*beta/exp(y)
             else
                f_nt(j, k, i) = 0.d0
             endif
           else
!          Add the broken power-law for the cases without thermal component.
             if(g.gt.gmin(j,k))then
              y = gmin(j,k)/gmax(j,k)
              f_nt(j,k,i) = Nnth/((gmin(j,k)**p_nth(j,k))*dexp(y))&
     &         /(g/gmin(j,k))**1.5!*dexp(sqrt(g/gmin(j,k))-1) !
             else
              f_nt(j,k,i) =0.d0
             endif
           endif
         else
            y = g/gmax(j,k)
            if (y.lt.1.d2) then
               f_nt(j, k, i) = Nnth/((g**p_nth(j, k))*dexp(y))
            else 
               f_nt(j, k, i) = 0.d0
            endif
         endif
         gnt(i) = g_1
!
         if (i.gt.1) then
!            sum = sum + f_nt(j, k, i)*(gnt(i) - gnt(i-1))
            sum = sum + f_nt(j, k, i-1)*(gnt(i) - gnt(i-1)) ! Xuhui 11/14/08
            Pnt(j, k, i-1) = sum
         endif
!         Pnt(j, k, i-1) = sum
         if (i.eq.1) then
            g_1 = 2.d-1 
         else
            g_1 = g_1*dg
         endif
 100  continue
      Pnt(j,k,num_nt) = sum
!
      n_eplus = 0.d0
      if (pair_switch.eq.1) then
         do 120 i = 1, num_nt-1
         n_eplus = n_eplus + n_pos(j, k, i)*(gnt(i+1) - gnt(i))
 120     continue
      endif
      do 140 i = 1, num_nt-1
         if (n_eplus.gt.1.d-10) then
            f_pos(i) = n_pos(j, k, i)/n_eplus
         else
            f_pos(i) = 0.d0
         endif
 140  continue
      f_pos(num_nt) = 0.d0
!
      do 200 i = 1, num_nt
        Pnt(j, k, i) = Pnt(j, k, i)/sum
        f_nt(j,k,i) = f_nt(j,k,i)/sum ! Xuhui 11/14/08
        if (wr_switch.eq.1)&
     &     write(nunit_fnt, 50) gnt(i), dmax1(1.d-30, f_nt(j, k, i)), &
     &                   dmax1(1.d-30, f_pos(i))
 200  continue
  50  format(e14.7, 1x, e14.7, 1x, e14.7)
      goto 500
!
 300  n_eplus = 0.d0
      if (pair_switch.eq.1) then
         do 320 i = 1, num_nt-1
         n_eplus = n_eplus + n_pos(j, k, i)*(gnt(i+1) - gnt(i))
 320     continue
      endif
      do 340 i = 1, num_nt-1
         if (n_eplus.gt.1.d-10) then
            f_pos(i) = n_pos(j, k, i)/n_eplus
         else
            f_pos(i) = 0.d0
         endif
 340  continue
      f_pos(num_nt) = 0.d0
!
      if (wr_switch.eq.1) then
         do 400 i = 1, num_nt
           write(nunit_fnt, 50) gnt(i), dmax1(1.d-30, f_nt(j, k, i)), &
     &                   dmax1(1.d-30, f_pos(i))
 400     continue
      endif
!
 500  if (wr_switch.eq.1) close(nunit_fnt)
      return
      end
!
!
!
!
      subroutine nth2d(j, k, gamm, betb, i)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer j, k
      double precision gamm, betb
!
!
      integer i
      double precision rnum, g
      double precision fibran
!
!
      rnum = fibran()
!     Xuhui added this line to get higher precision random number. 6/13/09
      rnum=int(rnum*1.d6)/1.d6+1.d-6*fibran()
      do 100 i = 2, num_nt
 100  if (Pnt(j, k, i).gt.rnum) goto 200
 200  continue
      gamm = sqrt(gnt(i)*gnt(i-1)) + 1.d0
      betb = dsqrt(1.d0 - 1.d0/(gamm*gamm))
      nelectron(i) = nelectron(i) + 1
!
      return
      end

!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:31:19 EDT 2006
! version: 2
! Name: J. Finke
! Changed 'random' common block.     
