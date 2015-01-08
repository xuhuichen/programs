      subroutine imcgen2d
!     this subroutine calculate some general stuff, like the radiative energy   
!     loss in volume, the total energy of the system, and photon numbers to be
!     produced in the current time step. 
!     All works are done in the master node
!     !comment by Xuhui Chen 7/8/08

      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      double precision sigma
      parameter(sigma = 1.0267d24)
!
      integer i, j, k, npwant, nst_local, nt, t, n
      integer n_new
      integer n_input
      integer syc_switch, jk_sample, dh_sentinel
!
      double precision sum_g_1, t_average
      double precision fac, bingo, fcens, gamma_R, fbias, bias
      double precision fas(jmax, kmax)
      double precision p_1, isy, Th, uB
      double precision g_av, Th_e, N_e_nt, f_rel_br
      double precision kgg_calc, int_sy, McDonald
      double precision gamma_bar
      double precision Eloss_pa, tau_integ
      real etime, elapse(2), et0

!
!
!
!        Formats
!
  10  format('Zone no. ',i2,',',i2)
  15  format('Eloss_sy = ',e14.7,' erg')
  20  format('Eloss_cy = ',e14.7,' erg')
  25  format('Eloss_th = ',e14.7,' erg')
  30  format('Eloss_br = ',e14.7,' erg (f_rel = ',e14.7,')')
  35  format('Eloss_pa = ',e14.7,' erg')
  40  format('Ingoing energy at inner boundary, zone ',i2,': ',e14.7,&
     &       ' erg')
  41  format('Ingoing energy at outer boundary, zone ',i2,': ',e14.7,&
     &       ' erg')
  42  format('Ingoing energy at upper boundary, zone ',i2,': ',e14.7,&
     &       ' erg')
  43  format('Ingoing energy at lower boundary, zone ',i2,': ',e14.7, &
     &       ' erg')
  50  format('Eloss_total = ',e14.7,' erg')
  55  format('Total ingoing energy: bingo = ',e14.7,' erg')
  60  format('Photons produced at inner boundary, zone ',i2,': ',i7)
  61  format('Photons produced at outer boundary, zone ',i2,': ',i7)
  62  format('Photons produced at upper boundary, zone ',i2,': ',i7)
  63  format('Photons produced at lower boundary, zone ',i2,': ',i7)
  64  format('Photons produced in zone (',i2,',',i2,'): ',i7)
  70  format(e14.7,1x,e14.7)
  71  format(e14.7,1x,e14.7,1x,e14.7)
! 
  90  format('Boundary temp. at upper boundary, zone ',i2,': ',&
     &        e14.7,' keV')
  91  format('Boundary temp. at lower boundary, zone ',i2,': ',&
     &        e14.7,' keV')
  92  format('Boundary temp. at inner boundary, zone ',i2,': ',&
     &        e14.7,' keV')
  93  format('Boundary temp. at outer boundary, zone ',i2,': ',&
     &        e14.7,' keV')
  94  format('Ed_abs = ',e14.7,' ergs; Ed_ref = ',e14.7,' ergs.')
!
!
!      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
!
!     master node
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(myid.eq.master) then
!
!
!
!
!
      n_input = 19
!
!
!       Initilaizing arrays
!
!      ndxout = 0
      do 100 j = 1, nz
         do 110  k = 1, nr
            npcen(j,k) = 0
            fas(j,k) = 0.d0
            edep(j,k) = 0.d0
            if (ncycle.eq.0) then
               ec_old(j,k) = 0.d0
            else
               ec_old(j,k) = ecens(j,k)
            endif
            prdep(j,k) = 0.d0
 110     continue
 100  continue 
!
      do 120 n = 1, nmu
         do 115 i = 1, nph_lc
 115        edout(n, i) = 0.d0
 120  continue
!
      bingo = 0.d0 
!
!
!     Calculate (time-dependent) surface emission
!
      if (ncycle.eq.0) then
         t = 1
         goto 135
      endif
      t_average = time + 5.d-1*dt(1)
!
      do 130 t = 1, ntime
         if (t1(t).gt.t_average) goto 135
 130  continue
 135  continue
      if(t.gt.ntime)then
        tbbu(:,t)=0.0d0
        tbbl(:,t)=0.0d0
        tbbi(:,t)=0.0d0
        tbbo(:,t)=0.0d0
      endif
 
 
!
!
      do 140 j = 1, nz
!         added 'time greater than t0' to make sure the external radiation begins after a certain time ! Xuhui 1/28/10
          if (tbbi(j,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t))then
             call file_sp(i_fname(j,t))
             erini(j) = dt(1)*Asurfi(j)*int_file
             write(4, 40) j, erini(j)
          else
             erini(j) = dt(1)*Asurfi(j)*sigma*(tbbi(j, t)**4.d0)
             if (tbbi(j,t).gt.0.d0) then
                write(4, 92) j, tbbi(j, t)
                write(4, 40) j, erini(j)
             endif
          endif
!
          if (tbbo(j,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(o_fname(j,t))
             erino(j) = dt(1)*Asurfo(j)*int_file
             write(4, 41) j, erino(j)
          else
             erino(j) = dt(1)*Asurfo(j)*sigma*(tbbo(j, t)**4.d0)
             if (tbbo(j,t).gt.0.d0) then
                write(4, 93) j, tbbo(j, t)
                write(4, 41) j, erino(j)
             endif
          endif
          erlko(j) = 0.d0
          erlki(j) = 0.d0
 140  continue    
!
      do 150 k = 1, nr          
          Ed_abs(k) = Ed_in(k) - Ed_ref(k)
          if (tbbu(k,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(u_fname(k,t))
             erinu(k) = dt(1)*Asurfu(k)*int_file
             write(4, 42) k, erinu(k)
          else
             erinu(k) = dt(1)*Asurfu(k)*sigma*(tbbu(k, t)**4.d0)
             if (star_switch.eq.1) &
     &            erinu(k) = erinu(k)*(Rstar/dist_star)**2.d0
!            if weight the luminosity by the radius of the star
!            and the distance to the star.
!            J. Finke 20 July 2006.
             if (tbbu(k,t).gt.0.d0) then
                write(4, 90) k, tbbu(k,t)
                write(4, 42) k, erinu(k)
             endif

          endif
!
          dh_sentinel = 0
          if (tbbl(k,t).lt.0.d0 .and. (time+0.5d0*dt(1)).ge.t0(t)) then
             call file_sp(l_fname(k,t))
             erinl(k) = dt(1)*Asurfl(k)*int_file
             write(4, 43) k, erinl(k)
          else
             erinl(k) = dt(1)*Asurfl(k)*sigma*(tbbl(k, t)**4.d0)
             if (dh_sentinel.eq.1) then
                erinl(k) = erinl(k) + Ed_abs(k)*dt(1)/dt(2)
                if (tbbl(k,t).lt.1.d-20) erinl(k) = 0.d0
                write(4, 94) Ed_abs(k), Ed_ref(k)
             endif
             if (tbbl(k,t).gt.0.d0) then
                write(4, 91) k, tbbl(k,t)
                write(4, 43) k, erinl(k)
             endif
          endif
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
!
 150  continue
      et0=etime(elapse)
      write(*,*)'before vol emiss in imcgen, et0',et0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate volume emissivities etc.
!
!
      open(16, file='n_ph1.dat', status='unknown')
      open(17, file='nph1_smooth.dat', status='unknown')
      open(18, file='n_ph2.dat', status='unknown')
      open(19, file='nph2_smooth.dat', status='unknown')
!
      Eloss_cy(:,:) = 0.d0
      Eloss_sy(:,:) = 0.d0
      Eloss_br(:,:) = 0.d0
      Eloss_th(:,:) = 0.d0
      Eloss_tot(:,:) = 0.d0
      do 400 j = 1, nz
         do 390 k = 1, nr
            Eloss_pa = 0.d0
            jk_sample = 0
            if((j.eq.1.or.mod(j,nz/3).eq.0).and.(k.eq.1.or.mod(k,nr/3).eq.0))jk_sample=1
!
           if(jk_sample)write(4, 10) j,k
!
!          Calculate the Magnetic Field 
            if (ep_switch(j,k).eq.1) then
               Th = 1.957d-3*tea(j,k)
               if (Th.lt.1.d-2) then
                  uB = 1.5d0*Th + 7.5d0*(Th**2.d0)
               else
                  uB = McDonald(3.d0, (1.d0/Th))&
     &                /McDonald(2.d0, (1.d0/Th)) - Th - 1.d0
               endif
               uB = uB*n_e(j,k)*8.176d-7*(1. + 2.d0*f_pair(j,k))
               B_field(j,k) = dsqrt(2.513d1*uB)
            else if (ep_switch(j,k).eq.2) then
               Th = 1.066d-6*tna(j,k)
               if (Th.lt.1.d-2) then
                  uB = 1.5d0*Th + 7.5d0*(Th**2.d0)
               else
                  uB = McDonald(3.d0, (1.d0/Th))&
     &                /McDonald(2.d0, (1.d0/Th)) - Th - 1.d0
               endif
               uB = uB*n_e(j,k)*1.5d-3
               B_field(j,k) = dsqrt(2.513d1*uB)
            endif            
!          flare caused by increase of magnetic field. Magnetic turbulence extend 5*dz distance
!          the thickness of the flare is translated, and rounded to number of vertical zones.
            B_field(j,k) = B_input(j,k)
            theta_local(j,k)=theta_b(j,k)
            acc_prob_local(j,k)=acc_prob
!------------------------------------------------------------------
            if(inj_switch.ne.0.and.k.le.int(sigma_r/dr+0.5d0).and. &
     &          (time-inj_t).gt.(j-1)*dz/inj_v.and.&
     &          (time-inj_t).lt.(j-1+int(sigma_z/dz+0.5d0))*dz/inj_v)then
!          inj_switch 134 added for a special case with both B and
!          acceleration changes
              if(inj_switch.eq.3)then
                B_field(j,k) = B_input(j,k)*dsqrt(flare_amp)
                write(*,*)'Amplified B_field=',B_field(j,k),'j=',j,'k=',k
              elseif(inj_switch.eq.134)then
                B_field(j,k) = B_input(j,k)*dsqrt(5.d0)
                write(*,*)'Amplified B_field=',B_field(j,k),'j=',j,'k=',k
!          inj_switch 5 changes theta while keeping B_z. Therefore B
!          increases.
!          inj_switch 6 changes theta while keeping B_phi. Therefore B
!          decreases.
              elseif(inj_switch.eq.5 .or. inj_switch.eq.115)then
                theta_local(j,k)=flare_amp/180*pi
                write(*,*)'During flares, theta=',theta_local(j,k),'j=',j
                B_field(j,k) = B_input(j,k)*dcos(theta_b(j,k))/dcos(theta_local(j,k))
              elseif(inj_switch.eq.6 .or. inj_switch.eq.116)then
                theta_local(j,k)=flare_amp/180*pi
                write(*,*)'During flares, theta=',theta_local(j,k),'j=',j
                B_field(j,k) = B_input(j,k)*dsin(theta_b(j,k))/dsin(theta_local(j,k))
!          inj_switch 145,146 and 125,126 added for a special case with both theta (to 75 degree) and
!          acceleration or external field changes
              elseif(inj_switch.eq.145.or.inj_switch.eq.125)then
                theta_local(j,k)=75.d0/180*pi
                write(*,*)'During flares, theta=',theta_local(j,k),'j=',j
                B_field(j,k) = B_input(j,k)*dcos(theta_b(j,k))/dcos(theta_local(j,k))
              elseif(inj_switch.eq.146.or.inj_switch.eq.126)then
                theta_local(j,k)=75.d0/180*pi
                write(*,*)'During flares, theta=',theta_local(j,k),'j=',j
                B_field(j,k) = B_input(j,k)*dsin(theta_b(j,k))/dsin(theta_local(j,k))
              elseif(inj_switch.eq.7)then
                acc_prob_local(j,k) = dmin1(acc_prob_local(j,k)*flare_amp,1.d0)
                write(*,*)'During flares, acc_prob=',acc_prob_local(j,k),'j=',j
              endif
            endif
!--------------------------------------------------------------------------
            if(inj_switch.eq.13)then
              if((time-inj_t)*inj_v.gt.((j-1)*dz+(k-1)*dr*(sigma_z/2)/r(nr)).and.&
     &        (time-inj_t)*inj_v.lt.(sigma_z+(j-1)*dz-(k-1)*dr*(sigma_z/2)/r(nr)))then
                B_field(j,k) = B_input(j,k)*dsqrt(flare_amp)
                write(*,*)'Amplified B_field=',B_field(j,k),'j=',j,'k=',k
              endif
            endif

            if(inj_switch.eq.15)then
              if((time-inj_t)*inj_v.gt.((j-1)*dz+(k-1)*dr*(sigma_z/2)/r(nr)).and.&
     &          (time-inj_t)*inj_v.lt.(sigma_z+(j-1)*dz-(k-1)*dr*(sigma_z/2)/r(nr)))then
                theta_local(j,k)=flare_amp/180*pi
                write(*,*)'During flares, theta=',theta_local(j,k),'j=',j,'k=',k
                B_field(j,k) = B_input(j,k)*dcos(theta_b(j,k))/dcos(theta_local(j,k))
              endif
            endif
!

!
!            call volume_em(zone)
            Th_e = tea(j,k)/5.11d2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        syc_switch = 0
!        syc_switch 1 means syc using power-law fit. 0 means without power-law fit.
        if(syc_switch.eq.1)then
!           use power law fit to calculate synchrotron loss
            if (amxwl(j,k).gt.9.9999d-1) then
               Eloss_sy(j,k) = 0.d0
               goto 200
            endif

         p_1 = 1. - p_nth(j,k)
         if (dabs(p_1).gt.1.d-3) then
            N_e_nt = (1.d0 - amxwl(j,k))*p_1*n_e(j,k)&
     &              /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
         else
            N_e_nt = (1.d0 - amxwl(j,k))*n_e(j,k)&
     &              /log(gmax(j,k)/gmin(j,k))
         endif

         g_av = gamma_bar(Th_e)
         gamma_R = 2.1d-3*sqrt(n_e(j,k))/(B_field(j,k)*dsqrt(g_av))

         isy = int_sy(gmin(j,k), gmax(j,k), p_nth(j,k), gamma_R)
         Eloss_sy(j,k) = 1.05838d-15*dt(1)*N_e_nt*(B_field(j,k)**2.d0)&
     &      *vol(j,k)*isy
         write(*,*) 'j,k=',j,k,' Eloss_sy=',Eloss_sy(j,k)
         write(*,*) 'dt1=',dt(1),' N_e_nt=',N_e_nt
       else
!      use actual electron distribution to calculate synchrotron loss
          sum_g_1 = 0.d0
            do 170 i=1, num_nt-1
  170       sum_g_1 = sum_g_1 + ((gnt(i)+1.d0)**2.d0 - 1.d0)*f_nt(j,k,i)&
     &                         *(gnt(i+1) - gnt(i))
          Eloss_sy(j,k) = 1.058d-15*n_e(j,k)*dt(1)*B_field(j,k)**2.d0&
     &                    *sum_g_1*vol(j,k)
       endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Xuhui syc
!
 200     continue
!       write(*,*)'gmax=',gmax(j,k),'gmin=',gmin(j,k),'slope=',p_nth(j,k)

!         Eloss_cy(j,k) = dt(1)*vol(j,k)*Eloss_cy(j,k)

!
!         Eloss_th(j,k) = dt(1)*zsurf(j,k)*Eloss_th(j,k)

!
         if (Th_e.gt.0.1d0) then
            f_rel_br = 1.41d0*dsqrt(Th_e)*(dlog(2.d0*Th_e) + 9.228d-1)&
     &               - 1.d0
            f_rel_br = 1.d0 + (Th_e**2.d0)*f_rel_br/(1.d0 + &
     &           (Th_e**2.d0))
            if (f_rel_br.lt.1.d0) f_rel_br = 1.d0
         else
            f_rel_br = 1.d0
         endif
!
         Eloss_br(j,k) = 5.34d-24*vol(j,k)*dt(1)*amxwl(j,k)&
     &                *sqrt(tea(j,k))*f_rel_br*(2.828d0*f_pair(j,k) &
     &                + 1.d0)*(n_e(j,k)**2.d0)         
!
!
         if (f_pair(j,k).lt.1.d-20) then
            Eloss_pa = 1.223d-20*vol(j,k)*dt(1)*f_pair(j,k)&
     &                *(1.d0 + f_pair(j,k))*(n_e(j,k)**2.d0)&
     &                /(1.d0/(1.d0 + 6.d0*Th_e)&
     &              + Th_e/(dlog(1.123d0*Th_e + 1.d0) + 2.5d-1))
         endif
!

!
!
!         Eloss_tot(j,k) = Eloss_sy(j,k) + Eloss_br(j,k) + Eloss_cy(j,k) 
!     1                  + Eloss_th(j,k) + Eloss_pa
         Eloss_tot(j,k) = Eloss_sy(j,k)
!        Xuhui 6/1/09 deactivated any emission except synchrotron.
!
         fas(j,k) = Eloss_tot(j,k)
!
         if(jk_sample)then
           if (Eloss_sy(j,k).gt.1.d-20) write(4,15) Eloss_sy(j,k)
           if (Eloss_cy(j,k).gt.1.d-20) write(4,20) Eloss_cy(j,k)
           if (Eloss_th(j,k).gt.1.d-20) write(4,25) Eloss_th(j,k)
           if (Eloss_br(j,k).gt.1.d-20) write(4,30) Eloss_br(j,k),f_rel_br
           if (Eloss_pa.gt.1.d-20) write(4,35) Eloss_pa
           write(4,50) Eloss_tot(j,k)
         endif
!
!          Check if blob is inside current zone;
!              if so, add blob emission
!
!         if ((r_blob.lt.xq(j,k)).and.(r_blob.gt.xq(j-1,k))) then
!            kT_bl = bc(3)*((xq(jm)/r_blob)**pT_b)
!            R_bl = h_r*xq(jm)
!            L_bl = 1.29d25*(kT_bl**4)*(R_bl**2)
!            E_bl = L_bl*1.d-8*dt(1)
!            fas(j) = Eloss_tot(j) + E_bl
!            p_blob(j) = E_bl/fas(j)
!
!
!         else
!            p_blob(j) = 0.d0
!            fas(j) = Eloss_tot(j)
!         endif
!
!
         if (pair_switch.eq.0) goto 385
!
!         write(*,*) 'Normalizing photon density spectra ...'
!
         do 370 i = 1, n_gg-1
             n_ph(i, j, k) = n_ph(i, j, k)&
     &                      /(vol(j, k)*(E_gg(i+1) - E_gg(i)))
             if (j.eq.1.and.k.eq.1) then
                write(16,70) E_gg(i), dmax1(1.d-10, n_ph(i, j, k))
             else if (j.eq.2.and.k.eq.1) then
                write(18,70) E_gg(i), dmax1(1.d-10, n_ph(i, j, k))
             endif
  370    continue
!
!         write(*,*) 'Smoothing photon density spectra ...'
!
         call nph_smooth(j, k, tea(j, k))
!
!
!       Calculate pair production opacities and rates
!
!         write(*,*) 'Calculating gamma-gamma opacities ...'
!
         do 380 i = 1, n_gg-1
             k_gg(i, j, k) = kgg_calc(E_gg(i), j, k)
             if (k_gg(i, j, k).lt.1.d-50) k_gg(i, j, k) = 0.d0
!
             if (j.eq.1.and.k.eq.1) then
                write(17, 71) E_gg(i), dmax1(1.d-10, n_ph(i, j, k)),&
     &                        dmax1(1.d-10, k_gg(i, j, k))
             else if (j.eq.2.and.k.eq.1) then
                write(19, 71) E_gg(i), dmax1(1.d-10, n_ph(i, j, k)),&
     &                        dmax1(1.d-10, k_gg(i, j, k))
             endif
  380     continue
!
!         write(*,*) 'Calculating pair production ...'
!
         call pairprod(j,k)        
!
!         write(*,*) 'Pair production rates done.'
!
  385     continue
  390     continue
  400 continue
!
      close(16)
      close(17)
      close(18)
      close(19)

!       calculate the optical depth across the region from bottom center
!       to top center.
        i = 41
          tau_integ = kappa_tot(i,1,1)*z(1)
          do j =2 ,nz
            tau_integ = tau_integ + kappa_tot(i,j,1)*(z(j)-z(j-1))
          enddo
          if(tau_integ.gt.1.)write(*,*)'tau_integ = ', tau_integ,'at i=',i
        
!
!
!      Calculate total input energy
!
      Emiss_old = Emiss_tot
      if(ncycle.eq.0)Emiss_old = 0.d0
      Emiss_tot = 0.d0
      do 420 j = 1, nz
         do 410 k = 1, nr
            bias = 1.d0
            bingo = bingo + bias*(ecens(j,k) + fas(j,k))
            Emiss_tot = Emiss_tot + fas(j,k)
  410    continue
  420 continue
!
      do 430 j = 1, nz
  430 bingo = bingo + erini(j) + erino(j)
      do 440 k = 1, nr
  440 bingo = bingo + erinu(k) + erinl(k)
      write(4,55) bingo
!
!     less photons are needed for the electron preparation phase
      if(inj_switch.ne.0.and.time.lt.(inj_t-5.d0*dmax1(z(nz),2*r(nr))/c_light))then
        nst_local = nst/30
        write(*,*)'Electron preparation, nst_local =',nst_local
      else
        nst_local = nst
      endif
!
      n_new = 0
!
      do 700 j = 1, nz
          nsurfi(j) = 0
          nsurfo(j) = 0
          if(tbbi(j,t).lt..0)nsurfi(j) = nst_local/nz
          if(tbbo(j,t).lt..0)nsurfo(j) = nst_local/nz
          n_new = n_new + nsurfi(j) + nsurfo(j)
 700  continue
       
      do 710 k = 1, nr
          nsurfu(k) = 0
          nsurfl(k) = 0
          if(k.eq.1)then
           if(tbbu(k,t).lt..0)nsurfu(k)=nst_local*(r(k)**2-rmin**2)/r(nr)**2
           if(tbbl(k,t).lt..0)nsurfl(k)=nst_local*(r(k)**2-rmin**2)/r(nr)**2
          else
           if(tbbu(k,t).lt..0)nsurfu(k)=nst_local*(r(k)**2-r(k-1)**2)/r(nr)**2
           if(tbbl(k,t).lt..0)nsurfl(k)=nst_local*(r(k)**2-r(k-1)**2)/r(nr)**2
          endif
          n_new = n_new + nsurfu(k) + nsurfl(k)
 710   continue      
       
       do 720 j = 1, nz
           do 730 k = 1, nr
             nsv(j,k) = 0.5d0*nst_local*fas(j,k)/Emiss_tot !(r(k)**2-r(k-1)**2)/r(nr)**2/nz
                         ! Xuhui new way to distribute MC particles 4/1/11
             n_new = n_new + nsv(j, k)
             if (nsv(j,k).gt.0) then
                ewsv(j,k) = fas(j,k)/dble(nsv(j,k))
             else
                ewsv(j,k) = 0.d0
             endif
 730      continue
 720  continue
!
! Calculating the surface energy weights
!
!
       do 800 j = 1,nz
          if (nsurfi(j).gt.0) then
             ewsurfi(j) = erini(j)/dble(nsurfi(j))
          else
             ewsurfi(j) = 0.d0
          endif
          if (nsurfo(j).gt.0) then
             ewsurfo(j) = erino(j)/dble(nsurfo(j)) 
          else
             ewsurfo(j) = 0.d0
          endif
  800  continue                    
                      
       do 810 k =1,nr                   
          if (nsurfu(k).gt.0) then
             ewsurfu(k) = erinu(k)/dble(nsurfu(k))
          else
             ewsurfu(k) = 0.d0
          endif
          if (nsurfl(k).gt.0) then
             ewsurfl(k) = erinl(k)/dble(nsurfl(k))
          else 
             ewsurfl(k) = 0.d0
          endif
  810  continue
!
!    Introduce "bias" correction if more than
!     10*nst_local particles would be produced
!
      if (n_new.gt.(10*nst_local)) then
         fbias = dble(10*nst_local)/dble(n_new)
         do 815 j = 1, nz
            nsurfi(j) = int(dble(nsurfi(j))*fbias)
            nsurfo(j) = int(dble(nsurfo(j))*fbias)
            ewsurfi(j) = ewsurfi(j)/fbias
            ewsurfo(j) = ewsurfo(j)/fbias
 815     continue
         do 820 k = 1, nr
            nsurfu(k) = int(dble(nsurfu(k))*fbias)
            nsurfl(k) = int(dble(nsurfl(k))*fbias)
            ewsurfu(k) = ewsurfu(k)/fbias
            ewsurfl(k) = ewsurfl(k)/fbias
 820     continue
        do 830 j = 1, nz
           do 825 k = 1, nr
!              nsv(j, k) = int(dble(nsv(j, k)/1000)*fbias)*1000
               nsv(j,k) = int(nsv(j,k)*fbias)
              ewsv(j, k) = ewsv(j, k)/fbias
 825       continue
 830    continue
!         
      endif
      do 835 j = 1, nz, nz/3
         if (nsurfi(j).gt.10) write(4,60) j, nsurfi(j)
         if (nsurfo(j).gt.10) write(4,61) j, nsurfo(j)
 835  continue
      do 836 k = 1, nr, nr/3
         if (nsurfu(k).gt.10) write(4,62) k, nsurfu(k)
         if (nsurfl(k).gt.10) write(4,63) k, nsurfl(k)
 836  continue  
!    
      do 838 j = 1, nz, nz/3
         do 837 k = 1, nr, nr/3
            if (nsv(j,k).gt.0) write(4,64) j, k, nsv(j,k)
 837     continue
 838  continue
!
!
      ecens(:,:) = 0.d0
      n_field(:,:,:) = 0.d0
!     k=1,nr
      Ed_abs(:) = 0.d0
      Ed_in(:) = 0.d0
      Ed_ref(:) = 0.d0
!     i=1,num_nt
      E_IC(:)=0.d0
!
      call seed_zone(rseed)


!
!     slave node
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else
!
      do 900 j=1, nz
          erlko(j) = 0.d0
          erlki(j) = 0.d0
 900  continue
      do 901 k=1, nr
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
 901  continue
!
!      ndxout = 0
      do 851 j = 1, nz
        do 841 k = 1, nr
            npcen(j,k) = 0
            fas(j,k) = 0.d0
            edep(j,k) = 0.d0
            ecens(j,k) = 0.d0
           do 842 i = 1, nphfield
 842          n_field(i, j, k) = 0.d0
 841       continue
 851    continue
!
      do 716 n = 1, nmu
         do 715 i = 1, nph_lc
 715        edout(n, i) = 0.d0
 716     continue
!
      do 861 k = 1, nr
         Ed_abs(k) = 0.d0
         Ed_in(k) = 0.d0
         Ed_ref(k) = 0.d0
 861  continue
!
      do 210 j=1, nz
          erlko(j) = 0.d0
          erlki(j) = 0.d0      
 210   continue
      do 220 k=1, nr
          erlku(k) = 0.d0
          erlkl(k) = 0.d0
 220   continue
      do i=1,num_nt
        E_IC(i) =0.d0
      enddo
!
!
!
!     end of imcgen2d
      endif
      return 
      end
!
!
!
!
!       Subroutine to integrate the synchrotron
!           emissivity over photon energy
!
!
      double precision function int_sy(gmin, gmax, p, g_R)
      double precision gmin, gmax, p, g_R
!
      double precision sum, g, gs, dg, d, s, sd, x, y
!
      dg = 1.d-2
      d = 1.d0 + dg
      s = 1.d0 + 5.d-1*d
      sum = 0.d0
      g = gmin
      y = 2.d0 - p
!
 100  gs = g*s
      x = -g_R/gs
      if (x.gt.-1.d2) then
         sd = (gs**y)*dexp(x)
         sum = sum + g*dg*sd
      endif
      g = g*d
      if (g.lt.gmax) goto 100
!
      int_sy = sum
      return
      end
!
! 
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Sun Aug  6 18:31:12 EDT 2006
! version: 3
! Name: J. Finke
! Added possibility of using upper surface to represent 
! a star. Weighs the calculation of erinu by 
! star parameters.       
