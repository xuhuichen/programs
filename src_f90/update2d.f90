!     This subroutines solves the Fokker-Planck
!      equation for the electrons and updates
!          electron spectral parameters
!     Processes includes radiative cooling, stochastic acceleration, particle escape, particle injection.
!
!
!
      subroutine update
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!

      double precision temp_min, temp_max


      parameter (temp_min = 5.d0)
      parameter (temp_max = 1.d3)

!
!
!   rad_cp = parameter limiting the minimum time step 
!            to ensure efficient radiative coupling
!            between zones. rad_cp = 3.33d-11 corresponds
!            to (Delta t)_min = (Delta R)/c for 1. zone
!
!   df_implicit = maximum relative temperature change for
!                 implicit Fokker-Planck time step
!
!
!    Tmin = minimum electron temperature [keV]
!    temp_max = maximum electron temperature [keV]
!
!    df_T = Maximum allowed relative temperature
!           increment/decrement per time step
!
!
!     variables used only in update2d and its subroutines.
      integer i, j, k, l, fp_steps, i_nt, zone
      integer nzacc_min, nzacc_max, nzacc_move, j_new, j_acc, nz_acc, nr_acc
      integer nunit_fnt
      character*30 fntname, emapname, nmapname, fmtele
      double precision t_fp, n_p
      double precision Te_old(jmax, kmax)
      double precision sum_dt, Delta_T, sum_min
      double precision Th_p, Th_e, g_av, gamma_bar, gamma_R
      double precision h_T, dT_coulp, y, volume, d_t, df_time
      double precision dT_sy, dT_c, dT_br, dT_A, dT_total(jmax,kmax)
      double precision E_el, E_pos, ne, ne_new, n_positron, n_lept
      double precision gamma1(num_nt), Delta_g
      double precision f_old(num_nt), npos_old(num_nt)
      double precision f_new(num_nt), npos_new(num_nt)
      double precision a_i(num_nt), b_i(num_nt), c_i(num_nt)
      double precision f_th, hr_th_Coul, hr_th_c, hr_th_sy
      double precision hr_th_br, hr_th_A, Omega, Om_p, v_a2, v_a, t_A
      double precision rhoh_wp, nu_A, k_min, k_max, vth_p, vth_e
      double precision te_mo, hr_th_total, sum_g11, sum_g_1
      double precision hr_nt_mo, hr_nt_C, hr_nt_br, hr_nt_sy
      double precision hr_nt_Coul, hr_nt_A, hr_st_A, heating_nt
      double precision hrmo_old, f_br, f_sy, fdisp_A, fdg_A
      double precision g_thr, Th_K2, McDonald, g_read, grid_temp
      double precision dg_cp(num_nt), disp_cp(num_nt)
      double precision dg_ce(num_nt), disp_ce(num_nt)
      double precision beta, The_mo, dg_sy(num_nt), dg_br(num_nt)
      double precision dg_ic(num_nt), dgcp_old, Intdgcp, dg_mo
      double precision p_g, k_res, om_R, Gamma_k, tau_k, xx
      double precision disp_mo, dg_A(num_nt), Intd2cp
      double precision disp_A(num_nt), fcorr_turb, dgA_original
      double precision dte_mo, hr_max, heat_total, rho, q_nm
      double precision D_g2, dgdt(num_nt), disp(num_nt)
      double precision sum_p, sum_old
      double precision pp_rate, pa_rate
      double precision dfmax, Delta_ne, Delta_np
      double precision d_temp, gbar, curv, curv_old, curv_2
      double precision sum_nt, sum_th, sump_old, dg2, p_1
      double precision sum_g, sumg_old, N_nt, f_pl, l_fraction
      double precision dE_fraction, dt_min, temp0, temp1
      double precision The_new
      double precision rmid, zmid, tl_flare, tlev, Tp_flare
      double precision fcorr_coul, f_disp_corr
!     particle diffusion
      double precision diff_z(0:jmax,kmax,num_nt),diff_r(jmax,kmax,num_nt), escape_total,var_scale,&
     &                 grad_ne(num_nt), fnt_tot(num_nt), ne_tot, ne_ave
      double precision r_acc_copy(jmax,kmax), fibran
!

!     variables common from outside modified by update.
!     gbar_nth and N_nth are not modified or used in update, but
!     are in the same common block as the others which are modified.
!     tea(jmax, kmax) also probably belongs here.

!     MPI variables
      integer status(MPI_STATUS_SIZE), num_sent, end_signal, sender,&
     &        num_zones
      double precision ans
!
!
!     variables in trid are shared between FP_calc, tridag and trid_p in one 
!     step
      common / trid_update / a_i, b_i, c_i, f_old, f_new, npos_old,&
     &                npos_new
!     Below are libraries from the subroutne coulomb.


!
!
!        Output formats
!
  5   format('Evolution of thermal population in zone ',i2,&
     &       ',',i2,':')
 10   format ('   Coulomb heating/cooling rate: ',e14.7,' erg/s')
 15   format ('       Synchrotron cooling rate: ',e14.7,' erg/s')
 20   format ('           Compton cooling rate: ',e14.7,' erg/s')
 25   format ('    Bremsstrahlung cooling rate: ',e14.7,' erg/s')
 30   format ('Hydromagnetic acceleration rate: ',e14.7,' erg/s')
 35   format ('     Total heating/cooling rate: ',e14.7,' keV/s')
 45   format('Te_new(',i2,',',i2,') = ',e14.7,' keV')
 50   format(' Adjusted time step: ',e14.7,' s')
 70   format ('Te = ',e14.7,'; file name: ',a20)
 75   format ('Tp = ',e14.7,'; file name: ',a20)
 95   format (' Acceleration rate of suprathermal particles: ',&
     &        e14.7,' erg/s.')
 1000 format ('     Total heat input (Coulomb + hydromagn.): ',&
     &        e14.7,' erg/s.')
 1020 format ('Total energy (old): ',e14.7,&
     &        ' ergs; (new): ',e14.7,' ergs.')
 1025 format ('(Delta E)/E = ',e14.7)
 1035 format ('  dT_max = ',e14.7)
 1040 format('           dt_new = ',e14.7,' s')
 1045 format('current time step = ',e14.7,' s')
 1050 format('New time step = ',e14.7,' s')
!
!
!     MPI Initialization 
!      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      num_sent = 0
      end_signal = jmax*kmax+1
      if(ncycle.lt.2) then
         call make_FP_bcast_type
      endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(myid.eq.master) then

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Calculate particle diffussion between zone boundaries
!        flow from j to j+1 or k to k+1 is positive 
!        unit: s^-1, i.e. flow particle #
!
!      Get central 2x2 particle density
      ne_tot = 0.d0
      do j=nz/2,nz/2+1
        do k=1,2
          ne_tot = ne_tot + n_e(j,k)*vol(j,k)
        enddo
      enddo
      ne_ave = ne_tot/pi/(z(nz/2+1)-z(nz/2-1))/r(2)**2
      write(*,*)'central 2x2 density:',ne_ave
!      Get average particle density
      ne_tot = 0.d0
      do j=1,nz
        do k=1,nr
          ne_tot = ne_tot + n_e(j,k)*vol(j,k)
        enddo
      enddo
      ne_ave = ne_tot/pi/z(nz)/r(nr)**2
      write(*,*)'average density:',ne_ave
!      Vertical boundaries, 0 ~ nz
      n_diff(:,:,:) = 0.d0
      if(esc_sw.ge.2)then
        dz = z(2)-z(1)
        dr = r(2)-r(1)
!       var_scale = half_scale*2, where half_scale is the minimum length scale 
!       it takes for the density to decrease to 1/2 of the high value
!        var_scale = 0.1d0*z(nz)
      if(fp_sw.eq.4)then
          do i = 1, num_nt
            diff_cf(i) = dz*c_light/r_esc *dmin1((gnt(i)+1),1.d5)/1.d5
          enddo
      else
          do i = 1, num_nt
            diff_cf(i) = dz*c_light/r_esc
          enddo
      endif
        write(*,*)'diff_cf=',diff_cf(1),diff_cf(50),diff_cf(200)
!       N is proportional to ln(R_max/r) from any non-source point to the 0 density boundary R_max 
!       (dN/dr*r=const)
        do k = 1, nr
          diff_z(0,k,:) = diff_cf(:)*-1.d0*f_nt(1,k,:)*n_e(1,k)/dz*vol(1,k)/dz
          do j = 1, nz-1
            grad_ne(:) = (f_nt(j,k,:)*n_e(j,k)-f_nt(j+1,k,:)*n_e(j+1,k))/dz
!            grad_ne(:) = dsign(
!     1      dmin1(abs(grad_ne(:)),dmax1(f_nt(j,k,:)*n_e(j,k),f_nt(j+1,k,:)*n_e(j+1,k))/var_scale),grad_ne(:))
            diff_z(j,k,:) = diff_cf(:)*grad_ne(:)*vol(1,k)/dz
          enddo
          diff_z(nz,k,:) = diff_cf(:)*f_nt(nz,k,:)*n_e(nz,k)/dz*vol(1,k)/dz
        enddo
        if(esc_sw.eq.3.or.esc_sw.eq.5)then
           diff_z(0,:,:) = 0.d0
           diff_z(nz,:,:) = 0.d0
        endif
!       Radial boundaries, 1 ~ nr (no inner boundary diffusion)
        do j= 1, nz
          do k = 1, nr-1
            grad_ne(:) = (f_nt(j,k,:)*n_e(j,k)-f_nt(j,k+1,:)*n_e(j,k+1))/dr
!            grad_ne(:) = dsign(
!     1      dmin1(abs(grad_ne(:)),dmax1(f_nt(j,k,:)*n_e(j,k),f_nt(j,k+1,:)*n_e(j,k+1))/var_scale),grad_ne(:))
            diff_r(j,k,:) = diff_cf(:)*grad_ne(:)*2*pi*r(k)*dz
          enddo
          diff_r(j,nr,:) = diff_cf(:)*f_nt(j,nr,:)*n_e(j,nr)/dr*2*pi*r(nr)*dz
        enddo
        if(esc_sw.eq.4.or.esc_sw.eq.5)diff_r(:,nr,:) = 0.d0
	write(*,*)'diff_r(1,nr,100)=',diff_r(1,nr,100)
!       change of particle spectra in each zone
!       unit: s^-1 cm^-3, i.e. density change, energy dependent
        do j=1,nz
          n_diff(j,1,:) = (diff_z(j-1,1,:)-diff_z(j,1,:)-diff_r(j,1,:))/vol(j,1)
          do k=2,nr
            n_diff(j,k,:) = (diff_z(j-1,k,:)+diff_r(j,k-1,:)-diff_z(j,k,:)-diff_r(j,k,:))/vol(j,k)
          enddo
        enddo
!       sweeping of ambient matter at the j=1 bottom boundary
!       the sweep compensates the particle number lost by escape
!       unit: s^-1 cm^-3
        escape_total = 0.d0
         do i = 1, num_nt-1
          do k = 1, nr
            escape_total = escape_total+(diff_z(nz,k,i)-diff_z(0,k,i))*(gnt(i+1)-gnt(i))
          enddo
          do j = 1, nz
            escape_total = escape_total+diff_r(j,nr,i)*(gnt(i+1)-gnt(i))
          enddo
         enddo   ! sum up of boundary particle escape
         write(*,*)'total escape number/s:',escape_total
!       escape from the central 2x2 region
         escape_total = 0.d0
         do i = 1, num_nt-1
          do k = 1,2 ! r numbers, for z diffusion
            escape_total = escape_total+(diff_z(nz/2+1,k,i)-diff_z(nz/2-1,k,i))*(gnt(i+1)-gnt(i))
          enddo
          do j = nz/2,nz/2+1 ! z numbers, for r diffusion
            escape_total = escape_total+diff_r(j,2,i)*(gnt(i+1)-gnt(i))
          enddo
         enddo
         write(*,*)'central 2x2 escape number/s:',escape_total
        if(ncycle.eq.1)then
          if(pick_sw.eq.2)sweep=3.2d-4*ne_ave/r_esc
          if(pick_sw.eq.3)sweep=2.6d-4*ne_ave/r_esc
          if(pick_sw.eq.4)sweep=4.d-3*ne_ave/r_esc ! pick up in central zones
          if(pick_sw.eq.5)sweep=4.d-3*ne_ave/r_esc ! pick up in acceleration zones
          write(*,*)'sweep=',sweep
        endif
      endif
!     end of calculating diffusion
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(turb_sw.eq.1)then
!       single zones in the center have random chance to accelerate
        r_acc(:,:)=1.d40
        nr_acc=2
        nz_acc=2
        do j=1,nz-nz_acc+1
          if(fibran().lt.acc_prob)then
            do j_acc=j,j+nz_acc-1
              do k=1,nr_acc
                r_acc(j_acc,k)=1/(1/r_acc(j_acc,k)+1/r_acc_peak)
                write(*,*)'r_acc(j,k)=',j_acc,k,r_acc(j_acc,k)
              enddo
            enddo
          endif
        enddo
!     begin random jumping of the accelerating zones, 1D
      else if(turb_sw.eq.2)then
        do j=1,nz
           if(r_acc(j,1).lt.1.d30)then
              nzacc_min=j
              exit
           endif
        enddo
        do j=nzacc_min,nz
           if(r_acc(j,1).lt.1.d30)nzacc_max=j
        enddo

        nzacc_move = int(nz*fibran())
        do while(nzacc_min+nzacc_move.le.nz.and.nzacc_max+nzacc_move.gt.nz)
           nzacc_move = int(nz*fibran())
        enddo
        write(*,*)'nzacc_move=',nzacc_move 
        r_acc_copy(:,:)=r_acc(:,:)
        do j=1,nz
          j_new=j+nzacc_move
          if(j_new.gt.nz)j_new=j_new-nz
          r_acc(j_new,:)=r_acc_copy(j,:)
        enddo
        write(*,*)'r_acc(:,1)=',r_acc(:,1)
      endif
!     end of random jumping of acceleration zone, 1D
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     Fraction of photon time step to be used
!       as elementary implicit FP time step
!
      f_t_implicit = 2.d-1
      lnL = 20.d0
      dT_max = 0.d0
!
      do 110 j = 1, nz
         do 100 k = 1, nr
            Te_new(j,k) = tea(j,k)
            Te_old(j,k) = tea(j,k)
 100     continue
 110  continue
!
      write(*,*) 'myid=', myid, ' before cens_add_up'
      call cens_add_up
      write(*,*) 'myid=', myid, ' after cens_add_up'
!
!
      if (ncycle.gt.1) goto 160
!     The initial time step lets the simulation volume fill
!     up with photons, then resets the timer to zero.
!
!     Only the master nodes does the initial time
!     step calculation.  JF, 2 Feb. 2005
      call photon_fill
!   After the first time step, start "update" here.
!
  160 continue
!
!
!       Output the current electron distribution

!
!     output the electron E_e file that contains the energy grid used for
!     all the electron spectrum files.
         if(ncycle.eq.1)then
            fntname = 'electrons/fnt_0001.dat'
            nunit_fnt=22
            open(nunit_fnt, file='electrons/E_e.dat', status='unknown')
            write(nunit_fnt,'("# Volume size: Z=",e14.7," R=",e14.7)')z(nz),r(nr)
            write(nunit_fnt, '("# nz= ",i4," nr= ",i4)') nz, nr
            write(nunit_fnt,'("# Energy grid for electron spectrum")')
            do i = 1, num_nt
               write(nunit_fnt, '(e14.7)') gnt(i)
            enddo
            close(nunit_fnt)
         endif
         fntname(15:15) = char(48 + int(ncycle/1000))
         fntname(16:16) = char(48 + int(ncycle/100)-10*int(ncycle/1000))
         fntname(17:17) = char(48 + int(ncycle/10)-10*int(ncycle/100))
         fntname(18:18) = char(48 + ncycle-10*int(ncycle/10))

!        calculate average electron spectrum across the whole region, and output it
!        at the end of the electron file
         fnt_tot(:)=0
         do i = 1, num_nt
           do j = 1, nz
             do k = 1, nr 
               fnt_tot(i) = fnt_tot(i) + f_nt(j,k,i)*n_e(j,k)*vol(j,k)
             enddo
           enddo
         enddo
         fnt_tot(:) = fnt_tot(:)/pi/z(nz)/r(nr)**2

         open(nunit_fnt, file=fntname, status='unknown')
         write(nunit_fnt,'("# ",i4,". time step: t = ",&
     &         e14.7," s")') ncycle,time
!         write(nunit_fnt,'("# j=1,k=1;...   j=1,k=5;...    j=1,k=9;
!     1 ...   j=15,k=1;...  j=15,k=5;...  j=15,k=9;
!     1 ... j=30,k=1;...  j=30,k=5;...  j=30,k=9")')

         write(fmtele,'("("a9","i4"e14.7)")'),'"#B=    "',nz*nr
         write(nunit_fnt, fmtele)&
     &        (((B_field(j,k)), k=1,nr), j=1,nz)
         write(fmtele,'("("a9","i4"e14.7)")'),'"#theta="',nz*nr
         write(nunit_fnt, fmtele)&
     &        (((theta_local(j,k)), k=1,nr), j=1,nz)
         write(fmtele,'("("i4"e14.7)")'),nz*nr+1
         do i = 1, num_nt
           write(nunit_fnt,fmtele)(((dmax1(1.d-30,f_nt(j,k,i)*n_e(j,k))), k=1,nr), j=1,nz),&
     &           dmax1(1.d-30,fnt_tot(i))
         enddo
         close(nunit_fnt)

!
!     Solve Fokker-Planck equation for thermal + nonthermal
!            particle population (MB, 20/July/2001)
!
      hr_st_total = 0.d0
      hr_total = 0.d0
      E_tot_old = 0.d0 
      E_tot_new = 0.d0
!
!
!
!
!     Broadcast data to all the processes needed in the
!     Fokker-Planck routine.
      write(*,*) 'myid=',myid,' b4 fp_bcast dt1=',dt(1),' t0=',t0(2)
      call FP_bcast
      write(*,*) 'myid=',myid,' af fp_bcast dt1=',dt(1),' t0=',t0(2)
      num_zones = nr*nz
!
!     This sends the first round of zones to the slaves for processing.
      do 900 l = 1, min(num_zones, (numprocs-1))
         zone = l   
         call FP_send_job(l, zone)
      write(*,*) 'myid=',myid, ' sending zone=', zone
         num_sent = num_sent + 1
  900 continue
!
!     As slaves complete processing a zone, this recieves the results
!     and sends the next zone to the slaves.
      do 902 l = 1, num_zones
!         write(*,*) 'myid=', myid, ' recieve result'
            call FP_recv_result(sender, zone)
!          write(*,*) 'myid=', myid, ' recieves zone=', zone, ' from node=', sender
          if(num_sent.lt.num_zones) then
               zone = num_sent+1
               call FP_send_job(sender, zone)
               num_sent = num_sent + 1
            else
               write(*,*) 'myid=', myid, ' sending end signal'
               call FP_send_end_signal(sender)
               write(*,*) 'myid=', myid, ' end signal sent to node=', &
     &              sender
            endif
 902  continue
!
      write(*,*) 'myid=',myid,' b4 old,new=',E_tot_old, E_tot_new
      call E_add_up
      write(*,*) 'myid=',myid,' af E_add_up'
      write(4, *)
      write(4, *)
      write(4,1000) hr_total
      write(4,95) hr_st_total
      dE_fraction = (E_tot_new - E_tot_old)/E_tot_old
      write(4, *) 'Energy Check:'
      write(4, 1020) E_tot_old, E_tot_new
      write(*,*) 'myid=',myid,' old, new:  ',E_tot_old, E_tot_new
      write(4, 1025) dE_fraction
      write(4, *)
!
 910  continue

      write(4, 1035) dT_max
      if (dT_max.lt.(0.2*df_T)) then ! Xuhui
         dt_new = 3.d0*dt(1)
      else if (dT_max.lt.(.75*df_T)) then
         dt_new = 1.1d0*dt(1)
      else if (dT_max.gt.(5.*df_T)) then ! Xuhui
         dt_new = 0.33d0*dt(1)
      else if (dT_max.gt.(1.25*df_T)) then
         dt_new = 7.5d-1*dt(1)
      else
         dt_new = dt(1)
      endif

      write(4,1040) dt_new
      dt(2) = dt(1)
      write(4,1045) dt(1)
!      if (dabs(dE_fraction).lt.1.d-2) then
!         dt(1) = dmin1((1.1d0*dt(1)), dt_new)
!      else if (dabs(dE_fraction).gt.2.d-2) then
!         dt(1) = dmin1((7.5d-1*dt(1)), dt_new)
!      endif ! Xuhui del 5/11/09

!      Set dt_min so that efficient radiative coupling
!         between adjacent zones is guaranteed
!
      dt_min = rad_cp*dmin1(dr, dz)

!
      write(4,1050) dt(1)
      write(4,*)
      write(4,*)
!
!     store results for next MC time step
      do 930 j = 1, nz
         do 920 k = 1, nr     
            if (tna(j,k).gt.1.) then
               tea(j,k) = Te_new(j,k)
               temp0 = temp_min
               temp1 = temp_max
               tea(j,k) = dmin1(temp1, tea(j,k))
               tea(j,k) = dmax1(temp0, tea(j,k))
            endif
 920     continue
 930  continue
      call FP_end_bcast  ! Xuhui
!
!
!     end master part
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      else if(myid.ne.master) then
!     beginning of slave part
!
      lnL = 20.d0
      write(*,*) 'myid=',myid,' b4 cens_add_up'
      call cens_add_up
      write(*,*) 'myid=',myid,' af cens_add_up'

      hr_st_total = 0.d0
      hr_total = 0.d0
      E_tot_old = 0.d0 
      E_tot_new = 0.d0
!
!     recieve broadcast of parameters used in FP_calc.
      write(*,*) 'myid=',myid,' b4 fp_bcast dt1=',dt(1),' t0=',t0(2)
      call FP_bcast
      write(*,*) 'myid=',myid,' af fp_bcast dt1=',dt(1),' t0=',t0(2)
!
      num_zones = nr*nz
!
!     if there are more nodes than work skip this node.
      if(myid.gt.num_zones) goto 990
!
!     as long as the node doesn't recieve the end_signal, it
!     will keep performing the FP calcuation for zones.
 991  call FP_recv_job(zone)
      call get_j_k(zone, j, k, nr)
!      write(*,*) 'myid=',myid, ' recieved, zone=', zone, ' j=',j,' k=',k
      if(zone.eq.end_signal) goto 990
      call FP_calc(zone)
!      write(*,*) 'myid=', myid, ' zone=', zone, ' done with fp_calc'
      call FP_send_result(zone)
      goto 991
 990  continue
!
      write(*,*) 'myid=',myid,' b4 old,new=',E_tot_old, E_tot_new
      call E_add_up
      write(*,*) 'myid=',myid,' af E_add_up'
      write(*,*) 'myid=',myid,' old, new:  ',E_tot_old, E_tot_new
      call FP_end_bcast  ! Xuhui 4/5/09
      endif
!     end of slave part
!
!
 950  return
      end
!
!
!
!
!========================================================================================
!========================================================================================
!     This performs the Fokker-Planck calculation.
!     Based on Change & Cooper (1970).
!     Xuhui Chen, 2009
      subroutine FP_calc(zone)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      
      integer zone
!
      double precision temp_min, temp_max
      parameter (temp_min = 5.d0)
      parameter (temp_max = 1.d3)
!
      logical ex
!
!     variables used only in FP_calc
      integer i, j, k, fp_steps,ii, inj_dis
      integer dgr_p, dgr_e, dge_cycle, dte_stop, temp_i, i_ph
      integer i_nt
      integer i_100000, i_10000, i_1000, i_100, i_10, i_1,i_5,&
     &        id_100,id_10,id_1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision sum_gg
      double precision n_inject, inject_ne(num_nt), inj_dur, inj_y,inj_sum,&
     & inj_E,inj_rho, t_esc, t_acc, inj_rate, inj_g2var, diff_ne(num_nt)
      double precision D_gminus,D_gplus,smw(num_nt)&
     &                 ,bigW(num_nt),bigB,bigC(num_nt)
      character *10  dateup, timeup
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision t_fp
      double precision Delta_T
      double precision Th_p, Th_e, g_av, gamma_bar, gamma_R
      double precision h_T, y, volume, d_t, d_dT_At
      double precision E_el, E_pos, ne, ne_new, n_positron, n_lept
      double precision gamma1(num_nt), Delta_g
      double precision f_old(num_nt), npos_old(num_nt)
      double precision f_new(num_nt), npos_new(num_nt)
      double precision a_i(num_nt), b_i(num_nt), c_i(num_nt)
      double precision f_th, hr_th_Coul, hr_th_c, hr_th_sy
      double precision hr_th_br, hr_th_A, Omega, Om_p, v_a2, v_a, t_A
      double precision rhoh_wp, nu_A, k_min, k_max, vth_p, vth_e
      double precision te_mo, hr_th_total, sum_g11, sum_g_1
      double precision hr_nt_mo, hr_nt_C, hr_nt_br, hr_nt_sy
      double precision hr_nt_Coul, hr_nt_A, hr_st_A, heating_nt
      double precision hrmo_old, f_br, f_sy, fdisp_A, fdg_A
      double precision g_thr, Th_K2, McDonald, g_read, grid_temp
      double precision dg_cp(num_nt), disp_cp(num_nt)
      double precision dg_ce(num_nt), disp_ce(num_nt)
      double precision beta, The_mo, dg_sy(num_nt), dg_br(num_nt)

      double precision dg_ic(num_nt), dgcp_old, Intdgcp, dg_mo
      double precision p_g, k_res, om_R, Gamma_k, tau_k, xx
      double precision disp_mo, dg_A(num_nt), Intd2cp
      double precision disp_A(num_nt), fcorr_turb, dgA_original
      double precision dte_mo, hr_max, heat_total, rho, q_nm
      double precision D_g2, dgdt(num_nt), disp(num_nt)
      double precision sum_p, sum_old
      double precision pp_rate, pa_rate
      double precision dfmax, Delta_ne, Delta_np
      double precision d_temp, gbar, curv, curv_old, curv_2
      double precision sum_nt, sum_th, dg2, p_1
      double precision sum_g, sumg_old, N_nt, f_pl, l_fraction
      double precision The_new
      double precision rmid, zmid, tl_flare, tlev, Tp_flare
      double precision fcorr_coul, f_disp_corr
      double precision n_p
      double precision dT_total(jmax,kmax)
!
!     variables common from outside update used within update and fp_calc.
!     the above Asurfs are not used in fp_calc but are in the common
!     block with z, r etc.
!     Below are libraries from the subroutine coulomb.
!
!     variables common between FP_calc, update, and photon_fill
!     variables common between update and FP_calc.
      double precision dg_icold
!
      character *30 name_dge, name_dgp, dgdtfile, dispfile
! 
      common / trid_update / a_i, b_i, c_i, f_old, f_new, npos_old,&
     &                npos_new
!     Below are libraries from the subroutine coulomb.
!     d_update has variables that are common between photon_fill, FP_calc, 
!     and update.
!     to_fp_calc has variables common between fp_calc and update
!____________________________________________________________________
!
      data name_dge/'rates/p001_dge0005.dat'/
      data name_dgp/'rates/p001_dgp0005.dat'/
!
 10   format ('   Coulomb heating/cooling rate: ',e14.7,' erg/s')
 15   format ('       Synchrotron cooling rate: ',e14.7,' erg/s')
 20   format ('           Compton cooling rate: ',e14.7,' erg/s')
 25   format ('    Bremsstrahlung cooling rate: ',e14.7,' erg/s')
 30   format ('Hydromagnetic acceleration rate: ',e14.7,' erg/s')
 40   format ('   Moeller heating/cooling rate: ',e14.7,&
     &        ' erg/s (kT_e = ',e14.7,' keV)')
! 55   format(e14.7,1x,e14.7)
 60   format ('Thermal temperature estimate: ',e14.7,' keV')
 65   format(e14.7,1x,e14.7,1x,e14.7)
 80   format(e12.5, 1x, e12.5, 1x, e12.5, 1x, e12.5, 1x,&
     &       e12.5, 1x, e12.5, 1x, e12.5, 1x, e12.5)
 85   format(e12.5, 1x, e12.5, 1x, e12.5, 1x, e12.5, 1x, e12.5)
 90   format ('   Total heating/cooling rate (FP electrons): ',&
     &        e14.7,' erg/s.')
 95   format (' Acceleration rate of suprathermal particles: ',&
     &        e14.7,' erg/s.')
 1000 format ('     Total heat input (Coulomb + hydromagn.): ',&
     &        e14.7,' erg/s.')
 1005 format('FP electron temperature in zone ',i2,&
     &       ',',i2,': ',e11.4,' keV (',e11.4,')')
 1010 format ('Pair fraction: ',e14.7)
 1015 format ('PP-rate: ',e14.7,' cm^(-3)/s; PA-rate: ',&
     &         e14.7,' cm^(-3)/s')
 1030 format('delta^2 = ',e14.7,'; l_wp/l_Coul = ',e14.7)
!
!
!
      call get_j_k(zone, j, k, nr)
      if(j.eq.1.and.k.eq.1)then
        write(*,*) 'zone=', zone, ' E_tot_new=', E_tot_new, ' j,k=',j,k
        write(*,*) 'zone=', zone, ' ecens=',ecens(j,k)
      endif
!
      if(esc_sw.eq.1)then
         t_esc = r_esc*z(nz)/c_light ! Xuhui escape
      else
         t_esc = 1.d40*z(nz)/c_light
      endif
      t_acc = r_acc(j,k)*z(nz)/c_light ! Xuhui acceleration
!===============================================================
!       Flare related to acceleration change
!----------------------------------------------------------
!          inj_switch 134 added for a special case with both B and
!          acceleration changes
!          inj_switch 145, 146 added for special cases with both theta and
!          acceleration changes
      if((inj_switch.eq.4.or.inj_switch.eq.134.or.inj_switch.eq.145.or.inj_switch.eq.146) &
     &        .and.k.le.int(sigma_r/dr+0.5d0))then
        if((time-inj_t).gt.(j-1)*dz/inj_v.and.&
     &     (time-inj_t).lt.(j-1+int(sigma_z/dz+0.5d0))*dz/inj_v)then
          t_acc = t_acc/flare_amp
          write(*,*)'increased acceleration: r_acc=',t_acc/z(nz)*c_light,'j=',j,'k=',k
        endif
!----------------------------------------------------------
      elseif(inj_switch.eq.14)then
          if((time-inj_t)*inj_v.gt.((j-1)*dz+(k-1)*dr*(sigma_z/2)/r(nr)).and.&
     &      (time-inj_t)*inj_v.lt.(sigma_z+(j-1)*dz-(k-1)*dr*(sigma_z/2)/r(nr)))then
           t_acc = t_acc/flare_amp
           write(*,*)'increased acceleration: r_acc=',t_acc/z(nz)*c_light,'j=',j,'k=',k
        endif
      endif
!        End of flare related to acceleration change
!=====================================================================
      t_fp = 0.d0
      fp_steps = 0
      Te_new(j,k) = tea(j,k)
!           if (tna(j,k).lt.1.d0) return 
!
!
!         Retrieve particle distribution
!
       volume = vol(j,k)
       E_el = 0.d0
       E_pos = 0.d0
       n_p = n_e(j,k)  ! proton density
       ne = n_p*(1. +f_pair(j,k)) ! electron density
       ne_new = ne
       n_positron = n_p*f_pair(j,k) ! positron density
       n_lept = ne + n_positron ! lepton density
       if (n_lept.lt.1.d-11) return ! Xuhui 11/20/08
       if(f_pair(j,k).gt.1.d-1)write(*,*)'f_pair=',f_pair(j,k)
       if(j.eq.1.and.k.eq.1)write(*,*)'n_e 1,1 =',ne
       if(j.eq.15.and.k.eq.9)write(*,*)'n_e 15,9 =',ne
       if(j.eq.30.and.k.eq.9)write(*,*)'n_e 30,9 =',ne
!
       do 180 i = 1, num_nt
          gamma1(i) = gnt(i) + 1.d0
          if (i.gt.1) then
             Delta_g = gnt(i) - gnt(i-1)
             if (pair_switch.eq.1) then
                n_positron = n_positron + Delta_g*n_pos(j,k,i)
                E_pos = E_pos + Delta_g*gamma1(i)*n_pos(j,k,i)
             endif
             E_el = E_el + Delta_g*gamma1(i)*f_nt(j, k, i)
          endif
  180  continue
!
       E_el = E_el*ne*8.176d-7*volume
       if (pair_switch.eq.1) E_pos = E_pos*8.176d-7*volume
       E_tot_old = E_tot_old + E_el + E_pos + ec_old(j, k)
       E_tot_new = E_tot_new + ecens(j, k)
!
      sum_p = 0.
      do 185 i = 1, num_nt-1
  185   sum_p = sum_p + (gnt(i+1) - gnt(i))*f_nt(j, k, i)
!
      do 190 i = 1, num_nt
          if (pair_switch.eq.1) npos_old(i) = n_pos(j, k, i)
          f_nt(j, k, i) = f_nt(j, k, i)/sum_p
          f_old(i) = f_nt(j, k, i)
  190  continue
       f_old(num_nt) = 0.d0
       if (pair_switch.eq.1) then
          npos_old(num_nt) = 0.d0
       endif

!       particle diffusion in one MC time step
        if(esc_sw.ge.2)then
           diff_ne(:) = n_diff(j,k,:)*dt(1)
           do i = 1, num_nt
             if(abs(diff_ne(i)/ne).lt.1.d-20)diff_ne(i)=0.d0
           enddo
           f_old(:) = f_old(:) + diff_ne(:)/ne
           n_inject = 0.d0
           do i = 1, num_nt-1
             if(f_old(i).lt.0.d0)then
               write(*,*)'f_nt negative',f_old(i),i
               stop
             endif
             n_inject = n_inject + diff_ne(i)*(gnt(i+1)-gnt(i))
           enddo
           ne_new = ne + n_inject
           ne = ne_new          ! electron density increases
           n_p = n_p + n_inject !  protons are also injected
           n_e(j,k) = n_p
           n_lept = n_lept+n_inject ! total lepton density increases

           sum_p = 0.
           do i = 1, num_nt-1
             sum_p = sum_p + (gnt(i+1) - gnt(i))*f_old(i)
           enddo
           do i = 1, num_nt
               f_old(i) = f_old(i)/sum_p
           enddo
           f_old(num_nt) = 0.d0
        endif



!
!            if ((ncycle.eq.2).and.(j.eq.1).and.(k.eq.1)) then
!               open(21, file='f_old.dat', status='unknown')
!               do 195 i = 1, num_nt
!  195          write(21, 55) gnt(i), f_old(i)
!               close(21)
!            endif
!
!         Calculate pair annihilation rates
!
       if (pair_switch.eq.1) call pa_calc(j, k, n_e)
!
!
!         Calculate energy loss, acceleration, 
!                and dispersion rates
!                  (MB, 20/July/2001)
!
!            Th_p = tna(j,k)/5.382d5
!
       if (k.gt.1) then
          rmid = 5.d-1*(r(k) + r(k-1))
       else
          rmid = 5.d-1*(r(k) + rmin)
       endif
       if (j.gt.1) then
          zmid = 5.d-1*(z(j) + z(j-1))
       else
          zmid = 5.d-1*(z(j) + zmin)
       endif
!
       if (cf_sentinel.eq.1) then
          y = 5.d-1*(((rmid - r_flare)/sigma_r)**2.d0&
     &             + ((zmid - z_flare)/sigma_z)**2.d0&
     &             + ((time - t_flare)/sigma_t)**2.d0)
!
         if (y.lt.1.d2) then
            tl_flare = flare_amp/dexp(y)
         else
            tl_flare = 0.d0
         endif
       else
         tl_flare = 0.d0
       endif
!
       tlev = turb_lev(j,k) + tl_flare
       Tp_flare = tna(j,k)*(1.d0 + tl_flare)
!
       Th_p = Tp_flare/9.382d5
       Th_e = tea(j,k)/5.11d2
       f_th = 1.5d0*volume*n_lept
!
!      Thermal cooling rates in erg/s
!
       hr_th_br = -Eloss_br(j,k)/dt(1)
!            hr_th_c = edep(j,k)/dt(1)
           do i=1,num_nt-1
              dg_ic(i) = 0.d0
              do i_ph = 1, nphfield
                    dg_ic(i) = dg_ic(i) &
     &                       - n_field(i_ph, j, k)*F_IC(i, i_ph)/volume
              enddo
           enddo
!
!
  200       g_av = gamma_bar(Th_e)

            hr_th_c =0.d0
            do i=1,num_nt-1
               hr_th_c = hr_th_c - 8.176d-7*dg_ic(i)*&
     &             f_old(i)*(gnt(i+1)-gnt(i))*volume*n_lept
            enddo

            if(fp_steps.gt.1000000) then 
               write(*,*) 'zone=',zone, &
     &              ' fp_steps=', fp_steps, ' t_fp=', &
     &              t_fp, ' dt(1)=', dt(1)
               write(*,*) ' d_t=',d_t, ' f_th=',f_th,&
     &              ' hr_th_total=',hr_th_total
               write(*,*) 'Te_new=',Te_new(j,k),&
     &              ' df_implicit=',df_implicit
               write(*,*) 'hr_th_Coul=',hr_th_Coul
               write(*,*) 'hr_th_sy=',hr_th_sy
               write(*,*) 'hr_th_br=',hr_th_br
               write(*,*) 'hr_th_C=',hr_th_C
               write(*,*) 'hr_th_A=',hr_th_A
               stop
            endif
            gamma_R = 2.1d-3*sqrt(n_lept)/(B_field(j,k)*sqrt(g_av))
!
            h_T = .79788*(2.*((Th_e + Th_p)**2.d0) + 2.d0*(Th_e + Th_p) &
     &             + 1.d0)/(((Th_e + Th_p)**1.5d0)*(1.d0 + 1.875d0*Th_e &
     &             + .8203d0*(Th_e**2.d0)))
            hr_th_Coul = f_th*1.7386d-26*n_p*lnL*h_T&
     &                  *(Tp_flare - Te_new(j,k))
!
            y = gamma_R/g_av
            if (y.lt.100.d0) then
!              hr_th_sy = -(Eloss_cy(j,k) + Eloss_sy(j,k) 
!     1                     + Eloss_th(j,k))/(dt(1)*dexp(y))
            hr_th_sy = -Eloss_sy(j,k)/(dt(1)*dexp(y))
            else
               hr_th_sy = 0.d0
            endif
!
            Omega = 1.76d7*B_field(j,k)
            Om_p = Omega/1.836d3
            v_a2 = 4.765d22*(B_field(j,k)**2.d0)/ne
            if (v_a2.gt.9.d20) v_a2 = 9.d20
            v_a = dsqrt(v_a2)
!
            if (k.eq.1) then
               dr = r(k) - rmin
            else
               dr = r(k) - r(k-1)
            endif
            if (j.eq.1) then
               dz = z(j) - zmin
            else
               dz = z(j) - z(j-1)
            endif
            t_A = dmin1(dr, dz)/v_a
!
!            if (ft_turb(j,k).gt.1.d-20) then
!               rhoh_wp = ft_turb(j,k)*hr_th_Coul/volume
!            else
!               rhoh_wp = turb_lev(j,k)*hr_th_Coul/volume
!            endif
!
            hr_th_A = tlev*hr_th_Coul
            if (hr_th_A.lt.1.d-20) hr_th_A = 1.d-20
            nu_A = .5d0*(q_turb(j,k) + 3.d0)
            k_min = 2.d0*pi/dmin1(dr, dz)
            k_max = Omega/dsqrt(v_a2)
            vth_p = c_light*dsqrt(3.*Th_p + 2.25d0*(Th_p**2.d0))&
     &             /(1.d0 + 1.5d0*Th_p)
            vth_e = c_light*dsqrt(3.*Th_e + 2.25d0*(Th_e**2.d0))&
     &             /(1.d0 + 1.5d0*Th_e)
            if (vth_p.gt.c_light) vth_p = c_light
            if (vth_e.gt.c_light) vth_e = c_light
!
!
!            hr_th_total = hr_th_Coul + hr_th_sy + hr_th_br + hr_th_c
!     1                  + hr_th_A
            hr_th_total = hr_th_sy + hr_th_c + hr_th_A ! Xuhui

!
!          Determine estimated (thermal) electron temperature
!      after current time step for implicit Moeller FP coefficients
!
            dT_total(j,k) = 6.25d8*dt(1)*hr_th_total/f_th
            f_t_implicit = df_implicit*Te_new(j,k)&
     &                    /dabs(dT_total(j,k))
            if (f_t_implicit.gt.df_T) f_t_implicit = df_T
            te_mo = Te_new(j,k) + f_t_implicit*dT_total(j,k)
!
            if (te_mo.lt.temp_min) then
               te_mo = temp_min
            else if (te_mo.gt.temp_max) then
               te_mo = temp_max
            endif
!
            sum_g11 = 0.d0
            do 210 i = 1, num_nt-1
  210       sum_g11 = sum_g11 + (gamma1(i)**1.1d0)*f_old(i)&
     &                         *(gnt(i+1) - gnt(i))
!
            sum_g_1 = 0.d0
            do 220 i=1, num_nt-1
  220       sum_g_1 = sum_g_1 + (gamma1(i)**2.d0 - 1.d0)*f_old(i)&
     &                         *(gnt(i+1) - gnt(i))
!
!
!
!            Create file names for Coulomb and Moeller  
!                energy loss and dispersion rates     
!
             hr_nt_mo = 1.d80
             dge_cycle = 0
             dte_stop = 0
!
             dgr_p = 0
 230         dgr_e = 0
             The_mo = 1.9569d-3*te_mo
             dge_cycle = dge_cycle + 1
             temp_i = 0
             grid_temp = 0.d0
 240         temp_i = temp_i + 1
             grid_temp = grid_temp + 1.d0
             if (dabs(te_mo - grid_temp).le.5.d-1) then
                i_1000 = temp_i/1000
                name_dge(15:15) = char(48+i_1000)
                i_100 = (temp_i - 1000*i_1000)/100
                name_dge(16:16) = char(48+i_100)
                i_10 = (temp_i - 1000*i_1000 - 100*i_100)/10
                name_dge(17:17) = char(48+i_10)
                i_1 = temp_i - 1000*i_1000 - 100*i_100 - 10*i_10
                name_dge(18:18) = char(48+i_1)
                id_100 =myid/100
                name_dge(8:8) = char(48+id_100)
                id_10 = myid/10-id_100*10
                name_dge(9:9) = char(48+id_10)
                id_1 = myid--id_100*100-id_10*10
                name_dge(10:10) = char(48+id_1)
!
             endif
!
             if ((te_mo.ge.(grid_temp + 5.d-1))) goto 240 
!
             grid_temp = 0.d0
             temp_i = 0
 245         temp_i = temp_i + 1
             grid_temp = grid_temp + 1.d1
!
         if ((dgr_p.eq.0).and.&
     &            (dabs(Tp_flare - grid_temp).le.5.d0)) then
         i_5 = temp_i/100000
         if (i_5.gt.0) then
            name_dgp(13:13) = char(48+i_5)
         else
            name_dgp(13:13) = 'g'
         endif
         i_100000 = (temp_i - 100000*i_5)/10000
         if((i_100000).gt.0.or.(i_5.gt.0)) then
            name_dgp(14:14) = char(48+i_100000)
         else
            name_dgp(14:14) = 'p'
         endif
         i_10000 = (temp_i - 100000*i_5 - 10000*i_100000)/1000
         name_dgp(15:15) = char(48+i_10000)
         i_1000 = (temp_i - 100000*i_5 - 10000*i_100000 -&
     &        1000*i_10000)/100
         name_dgp(16:16) = char(48+i_1000)
         i_100 = (temp_i - 100000*i_5 - 10000*i_100000&
     &        - 1000*i_10000 - 100*i_1000)/10
         name_dgp(17:17) = char(48+i_100)
         i_10 = temp_i -100000*i_5 - 10000*i_100000 - 1000*i_10000 -&
     &        100*i_1000 - 10*i_100
         name_dgp(18:18) = char(48+i_10) 
         id_100 =myid/100
         name_dgp(8:8) = char(48+id_100)
         id_10 = myid/10-id_100*10
         name_dgp(9:9) = char(48+id_10)
         id_1 = myid--id_100*100-id_10*10
         name_dgp(10:10) = char(48+id_1)
!     Below is the old way of getting name_dgp.  Changed to above.
!     J. Finke 20 June 2006
!             if ((dgr_p.eq.0).and.
!     1           (dabs(Tp_flare - grid_temp).le.5.d0)) then
!                i_5 = temp_i/10000
!                if (i_5.gt.0) name_dgp(9:9) = char(48+i_5)
!                i_10000 = (temp_i - 10000*i_5)/1000
!                name_dgp(10:10) = char(48+i_10000)
!                i_1000 = (temp_i - 10000*i_5 
!     1                 - 1000*i_10000)/100
!                name_dgp(11:11) = char(48+i_1000)
!                i_100 = (temp_i - 10000*i_5 - 1000*i_10000 
!     1                - 100*i_1000)/10
!                name_dgp(12:12) = char(48+i_100)
!                i_10 = temp_i - 10000*i_5 - 1000*i_10000 
!     1               - 100*i_1000 - 10*i_100
!                name_dgp(13:13) = char(48+i_10)
!
!                write(n3, 75) Tp_flare, name_dgp
!            if((j.eq.1).and.(k.eq.1)) write(*,*) 'name_dgp:  ', name_dgp
!
             endif
!          
             if (Tp_flare.ge.(grid_temp + 5.d0)) goto 245 
!
!          Inquire if energy loss / dispersion rate files exist;
!                   if files exist, attempt to read
!                   (check for electron energy grid)
!

            inquire(file=name_dge, EXIST=ex)
            if (ex) then
               open(unit=17, file=name_dge, status='unknown')
               do 250 i = 1, num_nt
                  read(17, 65) g_read, dg_ce(i), disp_ce(i)
                  if ((gamma1(i).gt.1.002).and.&
     &                (dabs((g_read - gamma1(i))/gamma1(i)).gt.1.d-2))&
     &            then
                      close(17)
                   goto 260
                  endif  
                  dg_ce(i) = dg_ce(i)*n_lept
                  disp_ce(i) = disp_ce(i)*n_lept
 250          continue
               dgr_e = 1
               close(17)
!
            endif
!
 260        continue
!
            inquire(file=name_dgp, EXIST=ex)
            if (ex) then
               open(unit=17, file=name_dgp, status='unknown')
               do 270 i = 1, num_nt
                  read(17, 65) g_read, dg_cp(i), disp_cp(i)
                  if ((gamma1(i).gt.1.002).and.&
     &                (dabs((g_read - gamma1(i))/gamma1(i)).gt.1.d-2))&
     &            then
                      close(17)
                      goto 280
                  endif
                  dg_cp(i) = dg_cp(i)*n_p
                  disp_cp(i) = disp_cp(i)*n_p
 270           continue
               dgr_p = 1
               close(17)
!
            endif
!
 280        continue
!
            if (k.gt.1) then
               rmid = 5.d-1*(r(k) + r(k-1))
            else
               rmid = 5.d-1*(r(k) + rmin)
            endif
            if (j.gt.1) then
               zmid = 5.d-1*(z(j) + z(j-1))
            else
               zmid = 5.d-1*(z(j) + zmin)
            endif
!
            if (cf_sentinel.eq.1) then
               y = 5.d-1*(((rmid - r_flare)/sigma_r)**2.d0&
     &                  + ((zmid - z_flare)/sigma_z)**2.d0&
     &                  + ((time - t_flare)/sigma_t)**2.d0)
!
               if (y.lt.1.d2) then
                  tl_flare = flare_amp/dexp(y)
               else
                  tl_flare = 0.d0
               endif
            else
               tl_flare = 0.d0
            endif
!
            tlev = turb_lev(j,k) + tl_flare
            fdg_A = pi*(q_turb(j,k) - 1.d0)*(Omega**(2.d0 &
     &             - q_turb(j,k)))&
     &           *(k_min**(q_turb(j,k) - 1.d0))*tlev*v_a2&
     &           *(c_light**(q_turb(j,k) - 3.d0))/q_turb(j,k)
            fdisp_A = 2.d0*pi*(q_turb(j,k) - 1.d0)*(Omega**(2.d0 &
     &                - q_turb(j,k)))*(k_min**(q_turb(j,k) - 1.d0))&
     &               *tlev*v_a2*(C_light**(q_turb(j,k) - 3.d0))&
     &               /(q_turb(j,k)*(q_turb(j,k) + 2.d0))
!            f_sy = (Eloss_sy(j,k) + Eloss_cy(j,k) + Eloss_th(j,k))
!     1             /(8.176d-7*volume*dt(1)*n_lept*sum_g_1)
!            f_sy = Eloss_sy(j,k)/(8.176d-7*volume*dt(1)*n_lept*sum_g_1)
            f_sy = 1.058d-15*B_field(j,k)**2/8.176d-7 ! Xuhui 2/18/11
!
            f_br = Eloss_br(j,k)&
     &            /(8.176d-7*volume*dt(1)*n_lept*sum_g11)
!
!
            hrmo_old = dabs(hr_nt_mo)
            heating_nt = 0.d0
            hr_st_A = 0.d0
            hr_nt_A = 0.d0
            hr_nt_C = 0.d0
            hr_nt_sy = 0.d0
            hr_nt_Coul = 0.d0
            hr_nt_br = 0.d0
            hr_nt_mo = 0.d0
            g_thr = 1.d0 + 4.d0*Th_e
            Th_K2 = The_mo*McDonald(2.d0, (1.d0/The_mo))
!
            do 300 i = 1, num_nt
               beta = dsqrt(1.d0 - 1.d0/(gamma1(i)**2.d0))
               y = gamma_R/gamma1(i)
               if (y.lt.100.d0) then
                  dg_sy(i) = -f_sy*(gamma1(i)**2.d0 - 1.d0)/dexp(y)
               else
                  dg_sy(i) = -1.d-50
               endif
!
               dg_br(i) = -f_br*(gamma1(i)**1.1d0)
!
!               dg_ic(i) = 0.d0
!               do 290 i_ph = 1, nphfield
!                  dg_ic(i) = dg_ic(i) 
!     1                     - n_field(i_ph, j, k)*F_IC(i, i_ph)/volume
!
! 290           continue ! Xuhui moved this to before the start of the FP loop 3/9/11
!
               if (dgr_p.eq.0) then
                  if (gamma1(i).lt.3.) then
                     dg_cp(i) = 1.194d-14*n_p*lnL&
     &                   *Intdgcp(gamma1(i), beta, Tp_flare)&
     &                  /((1. + 1.875d0*Th_p + .8203d0*(Th_p**2.d0))&
     &                   *sqrt(Th_p)*(gamma1(i)**2.d0)*beta)
                  else
                     dg_cp(i) = dgcp_old
                  endif
               endif
!
               dgcp_old = dg_cp(i)
!
               if (dgr_e.eq.0) then
                  dg_ce(i) = 1.496d-14*lnL*(n_lept/Th_K2)&
     &                      *dg_mo(gamma1(i), beta, The_mo)&
     &                      /(gamma1(i)*gamma1(i)*beta)
               endif
!
!
!         Account for Landau damping by reducing the
!          according to the average of the absorbed
!           wave energy through the current region
!                    (MB, 22/Sept/2000)
!
               p_g = gamma1(i)*beta
               k_res = Omega/(c_light*p_g)
               om_R = k_res*v_a
!
               if ((k_res.lt.k_min).or.(k_res.gt.k_max)) then
                  tau_k = 1.d30
               else
                  y = om_R - 2.d0*Om_p
                  if (dabs(y).lt.1.d-20) then
                     Gamma_k = 1.d50
                  else
                     Gamma_k = 1.77245d0*Om_p*((Om_p - om_R)**2.d0)&
     &                        /(om_R - 2.d0*Om_p)
                  endif
                  y = ((om_R - Omega)/(k_res*vth_e))**2.d0
                  if (y.lt.2.d2) then
                     xx = 1.836d3/(dexp(y)*k_res*vth_e)
                  else
                     xx = 0.d0
                  endif
                  y = ((om_R - Om_p)/(k_res*vth_p))**2.d0
                  if (y.lt.2.d2)&
     &               xx = xx + 1.d0/(dexp(y)*k_res*vth_p)
                  Gamma_k = Gamma_k*xx
                  tau_k = Gamma_k*t_A
!
               endif
               if (tau_k.lt.1.d-6) then
                  dg_A(i) = fdg_A*(p_g**(q_turb(j,k) - 1.d0))
               else if (tau_k.lt.2.d2) then
                  dg_A(i) = fdg_A*(p_g**(q_turb(j,k) - 1.d0))&
     &                     *(1.d0 - dexp(-tau_k))/tau_k
               else
                  dg_A(i) = fdg_A*(p_g**(q_turb(j,k) - 1.d0))/tau_k
               endif
!
               if (i.lt.num_nt) then
                  hr_nt_Coul = hr_nt_Coul + dg_cp(i)*f_nt(j, k, i)&
     &                                     *(gnt(i+1) - gnt(i))
                  hr_nt_A = hr_nt_A + dg_A(i)*f_nt(j, k, i)&
     &                                     *(gnt(i+1) - gnt(i))
               endif
!
               if (dgr_e.eq.0) then
!                   if (gamma1(i).lt.1.d3) then
!                   if (gamma1(i).lt.1.d6) then ! Xuhui 10/30/08
                      f_disp_corr = 2.5d-1
                      disp_ce(i) = f_disp_corr*2.99d-14*lnL&
     &                            *(n_lept/Th_K2)&
     &                            *disp_mo(gamma1(i), beta, The_mo)&
     &                            /(gamma1(i)*gamma1(i)*beta)
!                   else
!                      disp_ce(i) = 0.
!                   endif
               endif
!
               if (dgr_p.eq.0) then
                  if (gamma1(i).lt.3.) then
                     disp_cp(i) = 1.194d-14*n_p&
     &                   *Intd2cp(gamma1(i), beta, Tp_flare)&
     &                   /((Th_p**1.5d0)*(1.d0 + 1.875d0*Th_p &
     &                    + .8203d0*(Th_p**2.d0))*(gamma1(i)**2.d0)*beta)
                  else
                     disp_cp(i) = 0.d0
                  endif
               endif
!
               if (tau_k.lt.1.d-6) then
                  disp_A(i) = fdisp_A*(gamma1(i)**q_turb(j,k))&
     &                       *(beta**(q_turb(j,k) + 1.d0))
               else if (tau_k.lt.2.d2) then
                  disp_A(i) = fdisp_A*(gamma1(i)**q_turb(j,k))&
     &                       *(beta**(q_turb(j,k) + 1.d0))&
     &                       *(1.d0 - dexp(-tau_k))/tau_k
               else       
                  disp_A(i) = fdisp_A*(gamma1(i)**q_turb(j,k))&
     &                    *(beta**(q_turb(j,k) + 1.d0))/tau_k
               endif
!
  300       continue
!
!
            if (dgr_e.eq.0) then
               open(unit=17, file=name_dge, status='unknown')
               do 310 i = 1, num_nt
  310          write(17, 65) gamma1(i), (dg_ce(i)/n_lept), &
     &                        (disp_ce(i)/n_lept)
               g_read = 0.d0
               write(17, 65) g_read, g_read, g_read
               close(17)
            endif
!
            if (dgr_p.eq.0) then
               open(unit=17, file=name_dgp, status='unknown')
               do 320 i = 1, num_nt
 320           write(17, 65) gamma1(i), (dg_cp(i)/n_p), (disp_cp(i)/n_p)
               g_read = 0.d0
               write(17, 65) g_read, g_read, g_read
               close(17)
            endif
!
            do 330 i = 1, num_nt
 330        if (gamma1(i).gt.9.d0) goto 340
!
 340        p_g = dsqrt(gamma1(i)**2.d0 - 1.d0)
            dgA_original = fdg_A*(p_g**(q_turb(j,k) - 1.d0))
            fcorr_turb = dgA_original/dg_A(i)
            hr_nt_A = 0.d0
!
            hr_nt_Coul = hr_nt_Coul*8.176d-7*n_lept*volume
            fcorr_coul = hr_th_Coul/hr_nt_Coul
            hr_nt_Coul = hr_th_Coul

            if(fp_sw.eq.2)then
               dg_A(:) = gamma1(:)**1/t_acc*1.d5/gamma1(:)
               disp_A(:) = gamma1(:)**2/t_acc/2.d0*1.d5/gamma1(:)
            else if(fp_sw.eq.3)then ! First order Fermi
               dg_A(:) = gamma1(:)**1/t_acc
               disp_A(:) = 1.d-40
            else
               dg_A(:) = gamma1(:)**1/t_acc
               disp_A(:) = gamma1(:)**2/t_acc/2.d0
            endif
            do 350 i = 1, num_nt
!               dg_A(i) = dg_A(i)*fcorr_turb
!               disp_A(i) = disp_A(i)*fcorr_turb
               dg_cp(i) = dg_cp(i)*fcorr_coul
!               disp(i) = disp_ce(i) + disp_cp(i) + disp_A(i)
!               dgdt(i) = dg_sy(i) + dg_br(i) + dg_ic(i) + dg_cp(i) 
!     1                 + dg_ce(i) + dg_A(i)
!               disp(i) = disp_ce(i) + disp_A(i)
!               dgdt(i) = dg_sy(i) + dg_br(i) + dg_ic(i) 
!     1                 + dg_ce(i) + dg_A(i)
               disp(i) = disp_A(i)
               dgdt(i) = dg_sy(i) + dg_ic(i) + dg_A(i) ! Xuhui 6/11/08 deactivate some processes

!
               if (i.lt.num_nt) then
                  heating_nt = heating_nt + dgdt(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  hr_nt_sy = hr_nt_sy + dg_sy(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  hr_nt_br = hr_nt_br + dg_br(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  hr_nt_C = hr_nt_C + dg_ic(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  hr_nt_mo = hr_nt_mo + dg_ce(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  hr_nt_A = hr_nt_A + dg_A(i)*f_old(i)&
     &                                     *(gamma1(i+1) - gamma1(i))
                  if (gamma1(i).gt.g_thr) &
     &                hr_st_A = hr_st_A + dg_A(i)*f_old(i)&
     &                                   *(gamma1(i+1) - gamma1(i))
               endif
 350        continue
            dgdt(1) = 0.d0
            dgdt(num_nt) = 0.d0
            disp(1) = 0.d0
            disp(2) = 0.d0
            disp(num_nt-1) = 0.d0
            disp(num_nt) = 0.d0 ! momentum boundary flux is 0
!
            heating_nt = heating_nt*8.176d-7*n_lept*volume
            hr_nt_mo = hr_nt_mo*8.176d-7*n_lept*volume
            hr_st_A = hr_st_A*8.176d-7*n_lept*volume
            hr_nt_A = hr_nt_A*8.176d-7*n_lept*volume
            hr_nt_C = hr_nt_C*8.176d-7*n_lept*volume
            hr_nt_sy = hr_nt_sy*8.176d-7*n_lept*volume
            hr_nt_br = hr_nt_br*8.176d-7*n_lept*volume
!
            heat_total = hr_nt_Coul + hr_nt_A
            l_fraction = hr_nt_A/hr_nt_Coul
!
            hr_max = dmax1(dabs(hr_nt_Coul), dabs(hr_nt_sy))
            hr_max = dmax1(hr_max, dabs(hr_nt_C))
            hr_max = dmax1(hr_max, dabs(hr_nt_br))
            hr_max = dmax1(hr_max, dabs(hr_nt_A))
            if (dte_stop.eq.1) goto 360
!
            if (dabs(hr_nt_mo).gt.(1.d-1*hr_max)) then
               if (dabs(hr_nt_mo).gt.hrmo_old) then
                  te_mo = te_mo - dte_mo
                  dte_stop = 1
                  goto 230
               endif
               if (hr_nt_mo.gt.0.d0) then
                  dte_mo = -1.d0
               else
                  dte_mo = 1.d0
               endif
               te_mo = te_mo + dte_mo
               if ((te_mo.gt.temp_max).or.(te_mo.lt.temp_min) &
     &             .or.(dge_cycle.gt.50)) goto 360
               goto 230
            endif
!
 360        E_tot_old = E_tot_old + heat_total*f_t_implicit*dt(1)
!
            if (fp_steps.eq.0) then
               hr_total = hr_total + heat_total
               hr_st_total = hr_st_total + hr_st_A
            else
               goto 380
            endif
!
           if(myid.eq.1)then
            dgdtfile='rates/dgdt000.dat'
            dispfile='rates/disp000.dat'
            dgdtfile(11:11) = char(48 + int(ncycle/100))
            dgdtfile(12:12) = char(48 + int(ncycle/10)-&
     &                 10*int(ncycle/100))
            dgdtfile(13:13) = char(48 + ncycle-&
     &                 10*int(ncycle/10))
            dispfile(11:11) = char(48 + int(ncycle/100))
            dispfile(12:12) = char(48 + int(ncycle/10)-&
     &                 10*int(ncycle/100))
            dispfile(13:13) = char(48 + ncycle-&
     &                 10*int(ncycle/10))

            open(unit=24, file=dgdtfile, status='unknown')
            open(unit=25, file=dispfile, status='unknown')
            
            do 370 i = 1, num_nt
               write(24, 80) gnt(i), dg_sy(i), dg_br(i), dg_ic(i), &
     &                       dg_cp(i), dg_ce(i), dg_A(i), dgdt(i)
               write(25, 85) gnt(i), disp_ce(i), disp_cp(i), &
     &                       disp_A(i), disp(i)
 370        continue
            close(24)
            close(25)
           endif
!
  380       continue
            d_t = f_t_implicit*dt(1)
!
!           make the FP time step stops right at the end of the MC time step.
!           Xuhui 3/10/10
            if(d_t.gt.(dt(1)-t_fp))d_t=1.00001d0*(dt(1)-t_fp)
!
  390       continue
!
! Solve the Fokker-Planck Equation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!===================================================================
!           Implicit version of the
!       scheme of Nayakshin & Melia (1998)
!       **********************************
!       used the method in Chang and Cooper(1970) instead of N&M98
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           If pairs are present, 
!      re-normalize positron distribution
!
            if (pair_switch.eq.0) then
               f_pair(j,k) = 0.d0
               n_positron = 0.d0
               goto 465
            endif
!
            n_positron = f_pair(j,k)*n_p
            if (n_positron.gt.1.d-20) then
               sum_p = 0.d0
               do 430 i = 1, num_nt-1
 430           sum_p = sum_p + npos_new(i)*(gamma1(i+1) - gamma1(i))
!
               sum_p = sum_p/n_positron  ! normalize npos back to particle density
               if (sum_p.gt.1.d-20) then
                  do 440 i = 1, num_nt
                     npos_new(i) = npos_new(i)/sum_p
 440              continue
               endif
            endif
!
!
!         Add pair production and annihilation
!
 450        ne_new = ne
            pp_rate = 0.d0
            pa_rate = 0.d0
            n_positron = 0.d0
            do 460 i = 1, num_nt-1
              Delta_ne = (dn_pp(j,k,i) + dne_pa(j,k,i))*d_t
              Delta_np = (dn_pp(j,k,i) + dnp_pa(j,k,i))*d_t
!
              f_old(i) = f_old(i) + Delta_ne/ne
              npos_new(i) = npos_new(i) + Delta_np
              pp_rate = pp_rate +&
     &                  dn_pp(j,k,i)*(gamma1(i+1) - gamma1(i))
              pa_rate = pa_rate +&
     &                  dnp_pa(j,k,i)*(gamma1(i+1) - gamma1(i))
!
              if (f_old(i).lt.1.d-50) f_old(i) = 0.d0
              if (npos_new(i).lt.1.d-50) npos_new(i) = 0.d0
!

              if ((Delta_ne.gt.ne).or.(Delta_np.gt.ne)) then
                 write(*,*) 'j = ',j
                 write(*,*) 'gamma = ',gamma1(i)
                 write(*,*) 'Delta_ne = ',Delta_ne
                 write(*,*) 'Delta_np = ',Delta_np
                 write(*,*) 'dne_pp = ',dn_pp(j,k,i)
                 write(*,*) 'dne_pa = ',dne_pa(j,k,i)
                 write(*,*) 'dnp_pa = ',dnp_pa(j,k,i)
              endif
              n_positron = n_positron &
     &                   + npos_new(i)*(gamma1(i+1) - gamma1(i))
 460        continue
            if (n_positron.lt.1.d-50) n_positron = 0.d0
 465        ne_new = n_p + n_positron
            ne = ne_new
            f_pair(j,k) = dmax1(0.d0, (n_positron/n_p))

           
!____________________________________________________________________
!         Electron injection 
           n_inject = 0.d0
!----------------------------------------------------------------------
!         constant electron pick up at medium energy
          if(pick_sw.ne.0)then
           inj_sum = 0.d0
           inj_E = 0.d0
           do i=1, num_nt-1
              inject_ne(i) = pick_dis(i)
!             if(pick_sw.eq.1)then
!              inject_ne(i) = 1.d2*dexp(-(gamma1(i)-inj_gg)**2/
!     1                 2/inj_sigma**2)/(inj_sigma*dsqrt(2*pi))
!             else if(pick_sw.eq.2)then
!               inj_y = gamma1(i)/2.d0
!               if(gamma1(i).gt.1.d0)then
!                 inject_ne(i) = 1.d2/((gamma1(i)**0.d0)*dexp(inj_y))
!               else
!                 inject_ne(i) = 0.d0
!               endif
!             endif
                      inj_sum = inj_sum + inject_ne(i)*(gnt(i+1)-gnt(i))
                    inj_E = inj_E + inject_ne(i)*(gnt(i+1)-gnt(i))*&
     &                   gamma1(i)
           enddo

           if(pick_sw.eq.1)then
             inj_rho = pick_rate(j,k)*d_t
           elseif((pick_sw.eq.2.and.j.eq.1) &
     &         .or.(pick_sw.eq.3.and.j.eq.(nz/6+1).and.k.le.nr*2/3) &
!     1         .or.(pick_sw.eq.3.and.j.eq.(nz/6+1))&
     &         .or.(pick_sw.eq.4.and.(j.eq.nz/2+1.or.j.eq.nz/2).and.k.le.2)&
     &         .or.(pick_sw.eq.5.and.r_acc(j,k).lt.1.e35))then
!          particle sweep at the lower boundary
               inj_rho = sweep*d_t
           else
             inj_rho = 0.d0
           endif
           if(inj_switch.eq.6.and.k.le.int(sigma_r/dr+0.5d0))then
             if((time+t_fp-inj_t).gt.(j-1)*dz/inj_v.and.&
     &       (time+t_fp-inj_t).lt.(j-1+int(sigma_z/dz+0.5d0))*dz/inj_v)then
               inj_rho = inj_rho*flare_amp
               write(*,*)'increased pickup: ',inj_rho,'j=',j,'k=',k
             endif
           endif

           
           inject_ne(:) = inj_rho*inject_ne(:)/inj_sum
           f_old(:) = f_old(:) + inject_ne(:)/ne
           do i = 1, num_nt-1
             n_inject = n_inject + inject_ne(i)*(gnt(i+1)-gnt(i))
           enddo
          endif
!----------------------------------------------------------------------
!         Injection by a shock! Xuhui inj 11/17/08
!         inj_switch 115, 116 for cases with both injection theta change
          inj_dis=2
         if(inj_switch .eq. 1 .or. inj_switch .eq. 115 .or. inj_switch.eq.116) then
           inj_sum = 0.d0
           inj_E = 0.d0
           if((time+t_fp-inj_t).gt.dz/inj_v*(j-1) &
     &        .and.(time+t_fp-inj_t).lt.dz/inj_v*j &
     &        .and.k.le.int(sigma_r/dr+0.5d0))then ! small region injection
                  do i = 1, num_nt-1
!          determine the electron distribution      ccccccccccccccccccccc
                    if(inj_dis.eq.1) then
!          inject electrons of a Gaussian distribution
                     inject_ne(i) = 1.d2*dexp(-(gamma1(i)-inj_gg)**2/&
     &                    2/inj_sigma**2)/(inj_sigma*dsqrt(2*pi))
                    else if(inj_dis.eq.2) then
!          inject electrons of a Power-law distribution
                      inj_g2var = inj_g2*10**&
     &                     ((time+t_fp-inj_t)*inj_v/z(nz))
                      if(gamma1(i).gt.inj_g1)then
                          if(g2var_switch.eq.1)then
                            inj_y = gamma1(i)/inj_g2var
                          else
                            inj_y = gamma1(i)/inj_g2
                          endif
                          if(inj_y.lt.1.d2)then
                      inject_ne(i) =1.d2/((gamma1(i)**inj_p)*dexp(inj_y))
                          else
                             inject_ne(i) = 0.d0
                          endif
                      else
                          inject_ne(i) = 0.d0
                      endif
!          c        c       c       c       c        ccccccccccccccccccccc
                    endif
                      inj_sum = inj_sum + inject_ne(i)*(gnt(i+1)-gnt(i))
                    inj_E = inj_E + inject_ne(i)*(gnt(i+1)-gnt(i))*&
     &                   gamma1(i)
                  enddo
               inj_E = inj_E/inj_sum
               inj_rate = inj_L/8.186d-7/inj_E/(pi*r(nr)**2*dz)
               inj_rho = inj_rate*d_t
               write(*,*)'injection rate:(cm^-3*s^-1)',inj_rate
                  do i = 1, num_nt-1
                    inject_ne(i) = inj_rho*inject_ne(i)/inj_sum
                    f_old(i) = f_old(i) + inject_ne(i)/ne
                    n_inject = n_inject + inject_ne(i)*(gnt(i+1)-gnt(i))
                  enddo
                  inj_E= inj_E*8.186d-7*(pi*r(nr)**2*dz)*inj_rho/d_t
                  write(*,*)'injected L (erg/s) = , n_inject=',inj_E, n_inject

           endif

!          The t_esc/(t_esc+d_t) term is caused by particle escape.
!          This is calculated when using N^j+1 instead N^j with t_escape
         endif
           ne_new = ne + n_inject
           ne = ne_new          ! electron density increases
           n_p = n_p + n_inject !  protons are also injected
           n_e(j,k) = n_p
           n_lept = n_lept+n_inject ! total lepton density increases
!___________________________________________________________________
!          Particle escape
           ne_new = ne_new*t_esc/(t_esc+d_t)
           ne = ne_new
           n_p = n_p*t_esc/(t_esc+d_t)
           n_e(j,k) = n_p
           n_lept = n_lept*t_esc/(t_esc+d_t)

!          Calculate the coefficients
!
!
            a_i(1) = 0.d0
            b_i(1) = 1.d0
            c_i(1) = 0.d0
            a_i(num_nt) = 0.d0
            b_i(num_nt) = 1.d0
            c_i(num_nt) = 0.d0
!            do 410 i = 2, num_nt-1
!               Delta_g = (gnt(i+1) - gnt(i-1))
!               D_g2 = (gnt(i+1) - gnt(i-1))
!     1               *(gnt(i+1) - gnt(i))
!
!               q_nm = gnt(i+1)/gnt(i)
!
!               alpha = 2.d0/(1. + q_nm)
!               rho = 2.d0*q_nm/(1. + q_nm)
!               if (i.eq.3) rho = q_nm
!
!               if (i.eq.num_nt-1) then
!                  c_i(i) = 0.d0
!               else 
!                  c_i(i) = d_t*(dgdt(i+1)/Delta_g - disp(i+1)/D_g2)
!               endif
!
!               if (i.eq.2) then
!                  b_i(i) = 1.d0 + d_t*(dgdt(i)/Delta_g
!     1                               + disp(i)/(alpha*D_g2))
!               else if (i.eq.num_nt-1) then
!                  b_i(i) = 1.d0 - d_t*(dgdt(i)/Delta_g
!     1                           - rho*disp(i)/(alpha*D_g2))
!               else
!                  b_i(i) = 1.d0 + 2.d0*d_t*disp(i)/(alpha*D_g2)
!               endif
!
!               if (i.eq.2) then
!                  a_i(i) = 0.d0
!               else
!                  a_i(i) = -d_t*(dgdt(i-1)/Delta_g 
!     1                      + rho*disp(i-1)/(alpha*D_g2))
!               endif


!  410       continue

!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! The current solving of the equation is based on the equation used in the X-ray time lags paper
            do 475 i = 2, num_nt-1
               D_gminus = gnt(i) - gnt(i-1)
               D_gplus = gnt(i+1) - gnt(i)
               Delta_g = sqrt(gnt(i)/gnt(i-1))*D_gminus
               
               if(i.eq.2)then
                  bigB = -(dgdt(1)+dgdt(2))/2.d0 !+(disp(2)-disp(1))
!     1                  /D_gminus/2.d0
                  ! factor 1/2 in bigB added as considered a previous err. 2013/1/9
                  bigC(1) = (disp(1)+disp(2))/2.d0 !/4.d0
                  smw(1) = D_gminus*bigB/bigC(1)
                  bigW(1) = smw(1)/(exp(smw(1))-1.d0)
               endif

               bigB = -(dgdt(i)+dgdt(i+1))/2.d0 !+(disp(i+1)-disp(i))
!     1                /2.d0/D_gplus  ! Xuhui 4/24/11
               bigC(i) = (disp(i)+disp(i+1))/2.d0 !/4.d0 
               smw(i) = D_gplus*bigB/bigC(i)
               bigW(i) = smw(i)/(exp(smw(i))-1.d0)
!
               c_i(i) = -d_t*(bigC(i)*smw(i)/&
     &                     (1.d0-exp(-smw(i)))/Delta_g/D_gplus)

                b_i(i) = 1.d0 + d_t/Delta_g*(bigC(i)*bigW(i)/&
     &                    D_gplus+bigC(i-1)*smw(i-1)/&
     &                    (1.d0-exp(-smw(i-1)))/D_gminus) +d_t/t_esc

               a_i(i) = -d_t/Delta_g*bigC(i-1)*bigW(i-1)/D_gminus

  475       continue  ! Xuhui
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!
!

!         Now, solve the tridiagonal matrix of the 
!       implicitly discretized Fokker-Planck equation
!
            call tridag
            if (pair_switch.eq.1) call trid_p
!
            f_new(num_nt) = 0.d0
            if (pair_switch.eq.1) then
               npos_new(num_nt) = 0.d0
               npos_new(1) = 0.d0
            endif
            f_new(1) = 0.d0
!
!==================================================================
!


            dfmax = 0.
            sum_p = 0.
            sum_E(j,k) = 0.
            do 480 i = 1, num_nt-1
               sum_p = sum_p + (gnt(i+1) - gnt(i))*f_new(i)
               sum_E(j,k) = sum_E(j,k) + &
     &            (gnt(i+1) - gnt(i))*gamma1(i)*f_new(i)
               Pnt(j, k, i) = sum_p
 480        continue
            sum_E(j,k) = sum_E(j,k)/sum_p ! Xuhui inj
            sum_old = sum_p
            if(fp_steps.eq.1.and.j.eq.1.and.k.eq.1) write(*,*)'sum_p after FP is:',sum_p, 'j,k',j,k
!
!
!            Advance in time
!
            t_fp = t_fp + d_t
            fp_steps = fp_steps + 1
!
            do 485 i = 1, num_nt
               f_new(i) = f_new(i)/sum_p
               f_old(i) = f_new(i)
               if (pair_switch.eq.1) npos_old(i) = npos_new(i)
 485        continue
!
!
!_____________________________________________________________________
!            Determine new electron temperature
!
!
             gbar = 0.d0
             do 490 i = 1, num_nt-1
                gbar = gbar + gamma1(i)*f_new(i)*(gnt(i+1) - gnt(i))
 490         continue
!
             The_new = Th_e
             d_temp = .005*Th_e
             if(gbar.gt.g_av)then
                do 500, while(gbar.gt.g_av)
!                   The_new = The_new + d_temp
                   The_new = The_new*1.005 ! Xuhui
                   g_av = gamma_bar(The_new)
!                   if (The_new.gt.1.5) goto 520 ! Xuhui 10/30/08
 500            continue
             else
                do 510, while(gbar.lt.g_av)
!                  The_new = The_new - d_temp
                  The_new = The_new/1.005 ! Xuhui
                  g_av = gamma_bar(The_new)
                  if (The_new.lt.1.d-2) goto 520
 510            continue
             endif
             
!
 520        continue
!            if (The_new.gt.1.5d0) The_new = 1.5d0 ! Xuhui 10/30/08
            Te_new(j,k) = 5.11d2*The_new
            Th_e = The_new
!_______________________________________________________________________
!
!            Next time step
!
            if (t_fp.lt.dt(1)) goto 200
            write(*,*)'myid=',myid,'fp_steps=',fp_steps
            write(*,*)'tea=',Te_new(j,k)
!
!
!
!      Normalize new electron and positron distributions
!
            E_el = 0.d0
            E_pos = 0.d0
            do 540 i = 1, num_nt
               f_nt(j, k, i) = f_new(i)
               if (pair_switch.eq.1) n_pos(j,k,i) = npos_new(i)
               Pnt(j, k, i) = Pnt(j, k, i)/sum_p
               if (i.gt.1) then
                  E_el = E_el + f_nt(j, k, i)*gamma1(i)&
     &                         *(gnt(i) - gnt(i-1))
                  if (pair_switch.eq.1)&
     &               E_pos = E_pos + n_pos(j,k,i)*gamma1(i)&
     &                              *(gnt(i) - gnt(i-1))
               endif
 540        continue
            E_el = E_el*ne*8.176d-7*volume
            if (pair_switch.eq.1) E_pos = E_pos*8.176d-7*volume
            E_tot_new = E_tot_new + E_el + E_pos
!
            Delta_T = dabs(Te_new(j,k) - tea(j,k))/Te_new(j,k)
            if (Delta_T.gt.dT_max) dT_max = Delta_T
! 
            if(j.eq.1.and.k.eq.1)write(*, 1005) j, k, Te_new(j,k), tea(j,k)
!_______________________________________________________________________
!
!          Determine approx. parameters of 
!               nonthermal component:
!      gamma_1 = point where curvature of logarithm.
!                distrib. fct. becomes positive;
!           fit spectrum beyond gamma_1 with PL*Exp
!                     (MB, 14/May/1999)
            do 550 i = 5, num_nt-5
 550             if(f_new(i).gt.1.d-10)goto 555
 555        continue
            gmin(j,k) = gamma1(i)
            i_nt = i
            do 560 i = num_nt-5,5,-1
 560             if(f_new(i).gt.1.d-15)goto 565
 565        continue
            gmax(j,k)=gamma1(i)
            sum_nt = 0.d0
            sum_th = 0.d0
            do 570 i = 1, num_nt-1
               if (i.lt.i_nt) then
                  sum_th = sum_th + (gamma1(i+1) - gamma1(i))*f_new(i)
               else
                  sum_nt = sum_nt + (gamma1(i+1) - gamma1(i))*f_new(i)
               endif
  570       continue
            amxwl(j,k) = sum_th/(sum_nt + sum_th)
!
!
            if (amxwl(j,k).gt.9.999d-1) then
               amxwl(j,k) = 1.d0
               goto 620
            endif
!
!
! then
!               amxwl(j,k) = 1.d0
!               goto 620
!            else
!               p_nth(j,k) = dlog(f_new(i_nt+5)/f_new(i_nt+2))
!     1                     /dlog(gamma1(i_nt+2)/gamma1(i_nt+5))
!            endif 
!
!            gmax(j,k) = 2.d0*gmin(j,k)
!            dg2 = 5.d-2*gmax(j,k)
!            p_1 = 1.d0 - p_nth(j,k)
            p_nth(j,k) = 0.1
            sum_g = 1.d50
  580       sumg_old = sum_g
            sum_g = 0.d0
!            dg2 = 5.d-2*gmax(j,k)
            sum_gg = 0.d0
            p_1 = 1.d0 - p_nth(j,k) 
            if (dabs(p_1).gt.1.d-4) then
               N_nt = (1. - amxwl(j,k))*p_1&
     &               /(gmax(j,k)**p_1 - gmin(j,k)**p_1)
            else
               N_nt = (1.d0 - amxwl(j,k))/log(gmax(j,k)/gmin(j,k))
            endif
!            do 590 i = i_nt+2, num_nt-2
            do 590 i = i_nt, num_nt-2
              y = gamma1(i)/gmax(j,k)
              if (y.lt.100.d0) then
                 f_pl = N_nt/((gamma1(i)**p_nth(j,k))*exp(y))
                 sum_g = sum_g + f_pl*gamma1(i)*(gnt(i+1)-gnt(i))
                 sum_gg = sum_gg +f_pl*(gnt(i+1)-gnt(i))
!                 if (f_pl.gt.1.d-50)
!     1              sum_g = sum_g + ((f_new(i) - f_pl)**2.d0)/f_pl

              else
                 goto 600
              endif
  590      continue
  600      continue
           sum_g = sum_g/sum_gg 
           sum_g = abs(sum_g - sum_E(j,k))
           if (sum_g.lt.sumg_old) then
!              if (gmax(j,k).lt.(2.d0*gamma1(num_nt))) then
!                 gmax(j,k) = gmax(j,k) + dg2
             if(p_nth(j,k).lt.10.)then
               p_nth(j,k) = p_nth(j,k) + 0.5d-1
             else
                 goto 610
             endif
              goto 580
!
           endif
!
 610        continue
!
 620        continue
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            return
            end
!
!
!
!====================================================================================================
!     This performs the first time step.  It lets the simulation
!     volume fill up with photons.
!     J. Finke, 5 May 2005

      subroutine photon_fill
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!       Calculate thermal cooling/heating rates
!      and new temperature of thermal population
!                 (MB, 11/May/1999)
!
!     The initial time step lets the simulation volume fill
!     up with photons, then resets the timer to zero.
!
!     variables used only in update2d and its subroutines.
      integer i, j, k, i_ph, jk_sample
      double precision t_fp, n_p
      double precision sum_dt, sum_min, Delta_T
      double precision Th_p, Th_e, g_av, gamma_bar, gamma_R
      double precision h_T, dT_coulp, y, volume, d_t, df_time
      double precision dT_sy, dT_c, dT_br, dT_A, dT_total(jmax,kmax)
      double precision n_positron
      double precision rmid, zmid, tl_flare, tlev, Tp_flare
!
!     variables common from outside update used within update.
!     The Elosses are modified by photon_fill to adjust for the
!     size of the time step.  Te_new is also modified by 
!     photon_fill.
      double precision dg_ic(num_nt)
!
!     d_update has variables that are common between photon_fill, FP_calc, 
!     and update.
!
!     Output formats
  5   format('Evolution of thermal population in zone ',i2,&
     &       ',',i2,':')
 10   format ('   Coulomb heating/cooling rate: ',e14.7,' erg/s')
 15   format ('       Synchrotron cooling rate: ',e14.7,' erg/s')
 20   format ('           Compton cooling rate: ',e14.7,' erg/s')
 25   format ('    Bremsstrahlung cooling rate: ',e14.7,' erg/s')
 30   format ('Hydromagnetic acceleration rate: ',e14.7,' erg/s')
 35   format ('     Total heating/cooling rate: ',e14.7,' keV/s')
 45   format('Te_new(',i2,',',i2,') = ',e14.7,' keV')
 50   format(' Adjusted time step: ',e14.7,' s')
 55   format('maximum relative temperature change for implicit Fokker-Planck &
     &       time step. df_implicit=',e14.7)
!
!
!
      dt_new = dt(1)
      write(4,*) 'Starting thermal evolution calculation.'
!
      do 130 j = 1, nz
         do 129 k = 1, nr
            sum_dt = 0.d0
            Delta_T = 0.d0
            sum_min = -1.d-10
            jk_sample = 1
            if((j.eq.1.or.mod(j,nz/3).eq.0).and.(k.eq.1.or.mod(k,nr/3).eq.0))jk_sample=1
!
            n_p = n_e(j,k)
!            if (n_e(j,k).lt.1.d-2) goto 129
            if (n_e(j,k).lt.1.d-11) goto 129
            if (tna(j,k).lt.1.d0) goto 129
!
            write(4, *)
            if(jk_sample)write(4,5) j, k
!
!            Th_p = tna(j,k)/5.382d5
!
            if (k.gt.1) then
               rmid = 5.d-1*(r(k) + r(k-1))
            else
               rmid = 5.d-1*(r(k) + rmin)
            endif
            if (j.gt.1) then
               zmid = 5.d-1*(z(j) + z(j-1))
            else
               zmid = 5.d-1*(z(j) + zmin)
            endif
!
            if (cf_sentinel.eq.1) then
               y = 5.d-1*(((rmid - r_flare)/sigma_r)**2.d0&
     &                  + ((zmid - z_flare)/sigma_z)**2.d0&
     &                  + ((time - t_flare)/sigma_t)**2.d0)
!
               if (y.lt. 1.d2) then
                  tl_flare = flare_amp/dexp(y)
               else
                  tl_flare = 0.d0
               endif
            else
               tl_flare = 0.d0
            endif
            tlev = turb_lev(j,k) + tl_flare
!
            Tp_flare = tna(j,k)*(1.d0 + tl_flare)
            Th_p = Tp_flare/9.382d5
 125        Th_e = tea(j,k)/5.11d2
            g_av = gamma_bar(Th_e)
            gamma_R = 2.1d-3*sqrt(n_e(j,k))/(B_field(j,k)*sqrt(g_av))
!
!             Thermal cooling rates in erg/s
!
            h_T = .79788d0*(2.d0*((Th_e + Th_p)**2.d0) + 2.d0*&
     &           (Th_e + Th_p) + 1.d0)/(((Th_e + Th_p)**1.5d0)* &
     &           (1.d0 + 1.875d0*Th_e+ .8203d0*(Th_e**2.d0)))
            dT_coulp = 2.608d-26*n_p*lnL*(Tp_flare - tea(j,k))*h_T
!
            volume = vol(j,k)
!
            y = gamma_R/g_av
            if (y.lt.100.d0) then
!               dT_sy = -6.66667d-1*(Eloss_cy(j,k) + Eloss_sy(j,k) 
!     1            + Eloss_th(j,k))/(volume*n_e(j,k)*dt(1)*dexp(y))
               dT_sy = -6.66667d-1*Eloss_sy(j,k)/&
     &                 (volume*n_e(j,k)*dt(1)*dexp(y))
            else
               dT_sy = 0.d0
            endif
!
            n_positron = n_p*f_pair(j,k)
            dT_br = -6.66667d-1*Eloss_br(j,k)/(volume*n_e(j,k)*dt(1))
!
!            dT_c = 6.666667d-1*edep(j,k)/(dt(1)*volume*n_e(j,k))

         dT_c =0.d0
         do i=1,num_nt-1
            dg_ic(i) = 0.d0
            do i_ph = 1, nphfield
                  dg_ic(i) = dg_ic(i) &
     &                     - n_field(i_ph, j, k)*F_IC(i, i_ph)/volume
            enddo
            dT_c = dT_c - 6.666667d-1*8.176d-7*dg_ic(i)*&
     &             f_nt(j,k,i)*(gnt(i+1)-gnt(i))
         enddo  ! Xuhui 3/9/11
!
            dT_A = tlev*dT_coulp
!
!           Total thermal heating/cooling rate in keV/s
!
            dT_total(j,k) = (dT_coulp + dT_sy + dT_br + dT_c &
     &                     + dT_A)/1.6d-9
          if(jk_sample)then
            write(4,10) dT_coulp
            write(4,15) dT_sy
            write(4,20) dT_c
            write(4,25) dT_br
            write(4,30) dT_A
            write(4,35) dT_total(j,k)
            write(4, *)
          endif
!
!         Determine optimum time step, d_t, so that
!   electron temperature changes by a factor of df_T.
!
            d_t = df_T*tea(j,k)/dabs(dT_total(j,k))
!            if (d_t.lt.dt_new) dt_new = d_t ! Xuhui 2/17/11
!
  129    continue
  130 continue
!      if(dt_new.gt.dt(1)) dt_new = dt(1)
!
      df_time = dt_new/dt(1)
      do 150 j = 1, nz
         do 149 k = 1, nr
            edep(j,k) = edep(j,k)*df_time
            Eloss_tot(j,k) = Eloss_tot(j,k)*df_time
            Eloss_br(j,k) = Eloss_br(j,k)*df_time
            Eloss_cy(j,k) = Eloss_cy(j,k)*df_time
            Eloss_sy(j,k) = Eloss_sy(j,k)*df_time
            Eloss_th(j,k) = Eloss_th(j,k)*df_time
            Te_new(j,k) = tea(j,k) + dt_new*dT_total(j,k)
            dT_max = df_T
  149   continue
  150 continue       
!     The above loop renormalizes Eloss's to the new time step.
!      time = time - dt(1) + dt_new
!      dt(1) = dt_new ! Xuhui del 5/11/09
      write(4, 50) dt(1)
      write(4, 55) df_implicit
!
      return
      end
!
!
!
!============================================================================================
!     cens_add_up sums up necessary items from slave nodes:
!     edep, ecens, and n_field.
!     J. Finke, 22 June 2006
      subroutine cens_add_up
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer i, j, k
      double precision E_temp(jmax,kmax), ecens_temp(jmax,kmax), &
     &                 n_temp(nphfield,jmax,kmax)
      character*30 nfld_name
!
!
      Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      do 10 j=1, nz
!         do 20 k=1, nr
!            call MPI_REDUCE(edep(j,k), E_temp, 1, MPI_DOUBLE_PRECISION,
!     1           MPI_SUM, master, MPI_COMM_WORLD, ierr)
!            if(myid.eq.master) edep(j,k) = E_temp
!            call MPI_REDUCE(ecens(j,k), ecens_temp, 1, 
!     1           MPI_DOUBLE_PRECISION, MPI_SUM, master, 
!     2           MPI_COMM_WORLD, ierr)
!            if(myid.eq.master) ecens(j,k) = ecens_temp
!            do 30 i=1,nphfield
!               call MPI_REDUCE(n_field(i,j,k), E_temp, 1,
!     1              MPI_DOUBLE_PRECISION, MPI_SUM, master,
!     2              MPI_COMM_WORLD, ierr)
!               if(myid.eq.master) n_field(i,j,k) = E_temp
! 30         continue
! 20      continue
! 10   continue
       call MPI_REDUCE(edep, E_temp,jmax*kmax,&
     &        MPI_DOUBLE_PRECISION, MPI_SUM, master,&
     &        MPI_COMM_WORLD, ierr)
       if(myid.eq.master) edep=E_temp

       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_REDUCE(ecens, ecens_temp,jmax*kmax,&
     &        MPI_DOUBLE_PRECISION, MPI_SUM, master,&
     &        MPI_COMM_WORLD, ierr)
       if(myid.eq.master) ecens=ecens_temp

       Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_REDUCE(n_field, n_temp,nphfield*jmax*kmax,&
     &        MPI_DOUBLE_PRECISION, MPI_SUM, master,&
     &        MPI_COMM_WORLD, ierr)
       if(myid.eq.master) n_field=n_temp

      if(myid.eq.master)then
           nfld_name='temp/nfld_01_01.dat'
           do k=1,nr
            nfld_name(14:14)=char(48+int(k/10))
            nfld_name(15:15)=char(48+k-int(k/10))
            open(23, file=nfld_name, status='unknown')
            do i=1,nphfield
              write(23, *) E_field(i),dmax1(1.d-20,&
     &                     n_field(i,1,k)/vol(1,k))
            enddo
            close(23)
           enddo
      endif
!
      return
      end
!
!
!
!==================================================================================
!     E_add_up sums the erlk's from the different nodes and
!     zones and puts the sum in E_tot_old.  It also sums
!     the erin's and puts the sum in E_tot_new.
!     J. Finke, 2 May 2005
      subroutine E_add_up
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer i
      integer j, k
      double precision E_temp, E_tot_out
      double precision erlku_tot(kmax), erlkl_tot(kmax), &
     &                 erlko_tot(jmax), erlki_tot(jmax)
!
!     d_update has variables that are common between photon_fill, FP_calc, 
!     and update.
!     to_fp_calc has variables common between fp_calc and update
!
!      E_tot_old = 0.d0 
!      E_tot_new = 0.d0
       E_tot_out = 0.d0
!
      call MPI_REDUCE(erlki, erlki_tot, nz, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(erlko, erlko_tot, nz, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(erlku, erlku_tot, nr, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(erlkl, erlkl_tot, nr, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
!
      if(myid.eq.master) then
         write(*,*) 'erini(5)=',erinl(5),' erlkl5=',erlkl_tot(5)
      do 170 j = 1, nz
         E_tot_old = E_tot_old + erini(j) + erino(j)
         E_tot_new = E_tot_new + erlki_tot(j) + erlko_tot(j)
         E_tot_out = E_tot_out + erlki_tot(j) + erlko_tot(j) ! Xuhui Chen 
 170  continue
      do 175 k = 1, nr
         E_tot_old = E_tot_old + erinu(k) + erinl(k)
         E_tot_new = E_tot_new + erlku_tot(k) + erlkl_tot(k)
         E_tot_out = E_tot_out + erlku_tot(k) + erlkl_tot(k) ! Xuhui Chen

 175  continue
      endif
!
      call MPI_REDUCE(dT_max, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_MAX, master, MPI_COMM_WORLD, ierr)
      if(myid.eq.master) dT_max = E_temp
!
      call MPI_REDUCE(E_tot_new, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      if(myid.eq.master) E_tot_new = E_temp
      call MPI_REDUCE(E_tot_old, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)         
      if(myid.eq.master) E_tot_old = E_temp

      do i = 1, num_nt
        call MPI_REDUCE(E_IC(i), E_temp, 1, MPI_DOUBLE_PRECISION,&
     &       MPI_SUM, master, MPI_COMM_WORLD, ierr)
        if(myid.eq.master)E_IC(i) = E_temp
      enddo
      if(myid.eq.master)then
           open(23, file='output/eic.dat', status='unknown')
           do i=1,num_nt
             write(23, *) gnt(i), dmax1(1.d-20,E_IC(i))
           enddo
           close(23)
      endif
 
      call MPI_REDUCE(E_tot_out, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      if(myid.eq.master)then
         E_tot_out = E_temp ! Xuhui Chen
         write(4,*)'Luminosity =',E_tot_out/dt(1)
      endif

      call MPI_REDUCE(hr_total, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      if(myid.eq.master) hr_total = E_temp
      call MPI_REDUCE(hr_st_total, E_temp, 1, MPI_DOUBLE_PRECISION,&
     &     MPI_SUM, master, MPI_COMM_WORLD, ierr)
      if(myid.eq.master) hr_st_total = E_temp
!
!     end of E_add_up
      return
      end
!
!
!
!
      double precision function Intdgcp(g, b, kTp)
      implicit none
      double precision g, b, kTp
!
      double precision sum, gr, br, gcp, gcm, Om_m
      double precision Om_p, q, s, dgr, sgr, grs
      double precision sd, s0, gs, bs, d, me, mp
      double precision E10, E1s, p10, p1s, xm, xp
      double precision Om1, Om2
!
      sum = 0.
      s0 = 0.
      me = 5.11d2
      mp = 9.38d5
!
      gr = 1.d0
      dgr = 1.001d0
      d = dgr - 1.
      sgr = .5*(1. + dgr)
!
 100  grs = gr*sgr
      br = dsqrt(1. - 1./(grs**2.d0))
      s = mp**2 + me**2.d0 + 2.d0*mp*me*grs
      q = dsqrt(s)/kTp
      gs = (mp*grs + me)/dsqrt(s)
      bs = dsqrt(1.d0 - 1.d0/(gs**2.d0))
      E10 = me*g
      E1s = me*gs
      p10 = me*g*b
      p1s = me*mp*grs*br/dsqrt(s)
      gcp = (E10*E1s + p10*p1s)/(me**2.d0)
      gcm = (E10*E1s - p10*p1s)/(me**2.d0)
      xm = (mp + g*me)/kTp - q*gcm
      xp = (mp + g*me)/kTp - q*gcp
      if (xm.gt.-200.d0) then
         Om1 = dexp(xm)
      else
         Om1 = 0.d0
      endif
      if (xp.gt.-200.d0) then
         Om2 = dexp(xp)
      else
         Om2 = 0.d0
      endif
      Om_p = Om1 + Om2
      Om_m = Om1 - Om2
      sd = (Om_m*(g*((bs*gs)**2.d0) + gs/q) - Om_p*b*g*bs*(gs**2.d0))&
     &     /(grs*(br**3.d0))
      sum = sum + sd*gr*d
      gr = gr*dgr
      if (dabs(sd).gt.dabs(s0)) s0 = sd
      if ((dabs(s0).lt.1.d-100).or.(dabs(sd/s0).gt.1.d-8)) goto 100
!
      Intdgcp = sum
      return
      end
!
!
!
!
      double precision function Intd2cp(g, b, kTp)
      implicit none
      double precision g, b, kTp
!
      double precision sum, gr, br, gcp, gcm, eta0
      double precision eta1, eta2, q, s, dgr, sgr
      double precision grs, sd, s0, gs, bs, d, me, lnL
      double precision mp, E10, E1s, p10, p1s, const_A
      double precision const_B, Inteta, tau, p0, p1, p2
!
      sum = 0.
      s0 = 0.
      me = 5.11d2
      mp = 9.38d5
      lnL = 2.d1
!
      gr = 1.d0
      dgr = 1.001d0
      d = dgr - 1.
      sgr = .5*(1. + dgr)
!
 100  grs = gr*sgr
      br = dsqrt(1. - 1./(grs**2.d0))
      const_A = lnL - .25d0*(1.d0 + br**2.d0)
      const_B = lnL - .25d0*(6.d0 + br**2.d0)
      s = mp**2.d0 + me**2.d0 + 2.d0*mp*me*grs
      gs = (mp*grs + me)/dsqrt(s)
      bs = dsqrt(1.d0 - 1.d0/(gs**2.d0))
      E10 = me*g
      E1s = me*gs
      p10 = me*g*b
      p1s = me*mp*grs*br/dsqrt(s)
      gcp = (E10*E1s + p10*p1s)/(me**2.d0)
      gcm = (E10*E1s - p10*p1s)/(me**2.d0)
      q = dsqrt(s)/kTp
      tau = (mp + g*me)/kTp
      p0 = 0.d0
      p1 = 1.d0
      p2 = 2.d0
      eta0 = Inteta(gcm, gcp, p0, q, tau)
      eta1 = Inteta(gcm, gcp, p1, q, tau)
      eta2 = Inteta(gcm, gcp, p2, q, tau)
      sd = (-eta0*(const_A*((bs*gs)**2.d0) + const_B*(g**2.d0)) &
     &     + 2.d0*eta1*const_B*g*gs + eta2*(const_A*((bs*gs)**2.d0) &
     &     - const_B*(gs**2.d0)))/(gs*bs*(br**2.d0))
      sum = sum + sd*gr*d
      gr = gr*dgr
      if (dabs(sd).gt.dabs(s0)) s0 = sd
      if ((dabs(s0).lt.1.d-100).or.(dabs(sd/s0).gt.1.d-8)) goto 100
!
      Intd2cp = sum
      return
      end
!
!
!
!
      double precision function Intdgmo(g, b, Th_e)
      implicit none
      double precision g, b, Th_e
!
      double precision sum, gr, br, gcp, gcm, Om_m
      double precision Om_p, q, dgr, sgr, grs, sd
      double precision s0, gs, bs, d, Y, const_A
      double precision const_B, const_C, lnL, xm, xp
      double precision Om1, Om2
!
      sum = 0.d0
      s0 = 0.d0
      lnL = 2.d1
!
      gr = 1.d0
      dgr = 1.002d0
      if (Th_e.lt.2.d-1) dgr = 1.001d0
      d = dgr - 1.d0
      sgr = 5.d-1*(1.d0 + dgr)
!
 100  grs = gr*sgr
      br = dsqrt(1.d0 - 1.d0/(grs**2.d0))
      q = dsqrt(2.*(grs + 1.d0))/Th_e
      gs = dsqrt(.5d0*(grs + 1.d0))
      bs = dsqrt(1.d0 - 1.d0/(gs**2))
      gcp = g*gs*(1.d0 + b*bs)
      gcm = g*gs*(1.d0 - b*bs)
      xm = (1.d0 + g)/Th_e - q*gcm
      xp = (1.d0 + g)/Th_e - q*gcp
      if (xm.gt.-200.d0) then
         Om1 = dexp(xm)
      else
         Om1 = 0.d0
      endif
      if (xp.gt.-200.d0) then
         Om2 = dexp(xp)
      else
         Om2 = 0.d0
      endif
      Om_p = Om1 + Om2
      Om_m = Om1 - Om2
      const_A = (2.d0*(gs**2.d0) - 1.d0)**2.d0
      const_B = 2.*(gs**4.d0) - gs**2.d0 - .25d0
      const_C = (gs**2.d0 - 1.d0)**2.d0
      Y = 4.98885d-25*(.5d0*const_A*(lnL + 8.465736d-1) &
     &         - 6.9314718d-1*const_B + .125d0*const_C)&
     &                /((gs*(gs**2.d0 - 1.d0))**2.d0)
      sd = grs*br*Y*(Om_m*(g*((bs*gs)**2.d0) + gs/q) &
     &             - Om_p*b*g*bs*(gs**2.d0))
      sum = sum + sd*gr*d
      gr = gr*dgr
      if (dabs(sd).gt.dabs(s0)) s0 = sd
      if ((dabs(s0).lt.1.d-100).or.(dabs(sd/s0).gt.1.d-10)) goto 100
!
      Intdgmo = sum
      return
      end
!
!
!
!
      double precision function Intd2mo(g, b, Th_e)
      implicit none
      double precision g, b, Th_e
!
      double precision sum, gr, br, gcp, gcm, I1
      double precision I2, q, dgr, sgr, grs, sd
      double precision s0, gs, bs, d, const_A, p0
      double precision const_B, const_C, lnL, Inteta
      double precision eta0, eta1, eta2, tau, p1, p2
!
      sum = 0.
      s0 = 0.
      lnL = 20.
!
      gr = 1.d0
      dgr = 1.002d0
      if (Th_e.lt.2.d-1) dgr = 1.001d0
      d = dgr - 1.d0
      sgr = 5.d-1*(1. + dgr)
!
 100  grs = gr*sgr
      br = dsqrt(1.d0 - 1.d0/(grs**2.d0))
      q = dsqrt(2.d0*(grs + 1.d0))/Th_e
      gs = dsqrt(.5d0*(grs + 1.d0))
      bs = dsqrt(1.d0 - 1.d0/(gs**2.d0))
      gcp = g*gs*(1.d0 + b*bs)
      gcm = g*gs*(1.d0 - b*bs)
      const_A = (2.d0*(gs**2.d0) - 1.d0)**2.d0
      const_B = 2.d0*(gs**4.d0) - gs**2.d0 - .25d0
      const_C = (gs**2.d0 - 1.d0)**2.d0
!
      I1 = 4.98885d-25*(.5d0*const_A - 3.8629436d-1*const_B &
     &    + 8.3333333d-2*const_C)/((gs*(gs**2.d0 - 1.d0))**2.d0)
      I2 = 4.98885d-25*(const_A*(lnL + .34657d0) - const_B &
     &     + 1.6666667d-1*const_C)/((gs*(gs**2.d0 - 1.d0))**2.d0)
!
!       I1 = 0.d0
!       I2 = 4.98885d-25*const_A*(lnL + .34657)
!     1      /((gs*(gs**2 - 1.))**2)
!
      tau = (1.d0 + g)/Th_e
      p0 = 0.d0
      p1 = 1.d0
      p2 = 2.d0
      eta0 = Inteta(gcm, gcp, p0, q, tau)
      eta1 = Inteta(gcm, gcp, p1, q, tau)
      eta2 = Inteta(gcm, gcp, p2, q, tau)
!
      sd = ((grs**2.d0 - 1.d0)/bs*gs)*(eta0*(I1*(g**2.d0)&
     &    - .5d0*I2*(g**2.d0 + (bs*gs)**2.d0)) + 2.d0*eta1*g*gs&
     &     *(.5d0*I2 - I1) - eta2*(.5d0*I2 - I1*(gs**2.d0)))
      sum = sum + sd*gr*d
      gr = gr*dgr
      if (dabs(sd).gt.dabs(s0)) s0 = sd
      if ((dabs(s0).lt.1.d-100).or.(dabs(sd/s0).gt.1.d-10)) goto 100
!
      Intd2mo = sum
      return
      end
!
!
!
!
!        Moeller Energie exchange coefficient 
!         neglecting large-angle scatterings:
!              (Nayakshin & Melia 1998)
!
!
      double precision function dg_mo(g, b, Th)
      implicit none
      double precision g, b, Th
!
      double precision sum, sd, x, d, y, chi, gs, bs
      double precision gplus, gminus, ch_f
!
      sum = 0.d0
      x = 1.d0
      d = 1.d-3*Th
!
 100  gs = x + 5.d-1*d
      bs = dsqrt(1.d0 - 1.d0/(gs**2.d0))
      y = gs/Th
      if (y.lt.5.d2) then
         gplus = g*gs*(1.d0 + b*bs)
         gminus = g*gs*(1.d0 - b*bs)
         if (gplus.gt.1.0001d0*gminus) then
            chi = ch_f(gplus) - ch_f(gminus)
            sd = 5.d-1*(gs - g)*chi/dexp(y)
            sum = sum + d*sd
         endif
      endif
      x = x + d
      if (x.lt.(1.d0 + 1.d1*Th)) goto 100
!
      dg_mo = sum
      return
      end
!
!
!
!       Moeller Energie dispersion coefficient 
!         neglecting large-angle scatterings:
!              (Nayakshin & Melia 1998)
!
!
      double precision function disp_mo(g, b, Th)
      implicit none
      double precision g, b, Th
!
      double precision sum, sd, x, d, y, chi, zeta, gs, bs
      double precision gplus, gminus, z_f, ch_f
!
      sum = 0.d0
      x = 1.d0
      d = 1.d-3*Th
!
 100  gs = x + 5.d-1*d
      bs = dsqrt(1.d0 - 1.d0/(gs**2.d0))
      y = gs/Th
      if (y.lt.5.d2) then
         gplus = g*gs*(1.d0 + b*bs)
         gminus = g*gs*(1.d0 - b*bs)
         if (gplus.gt.1.0001d0*gminus) then
            chi = ch_f(gplus) - ch_f(gminus)
            zeta = z_f(g, gs, gplus) - z_f(g, gs, gminus)
            sd = (-5.d-1*((g - gs)**2.d0)*chi + zeta)/dexp(y)
            sum = sum + d*sd
         endif
      endif
      x = x + d
      if (x.lt.(1.d0 + 1.d1*Th)) goto 100
!
      disp_mo = sum
      return
      end
!
!
!
      double precision function ch_f(x)
      implicit none
      double precision x
!
      double precision z, x1, x2, x3, z_f
!
      if (x.lt.1.00000001d0) then
         z_f = 0.d0
         ch_f = 0.d0 ! Xuhui 5/19/09
         goto 100
      endif
      z = dsqrt(5.d-1*(x - 1.d0))
      x1 = 2.d0*dlog(z + dsqrt(z**2.d0 + 1.d0))
      x2 = dsqrt(x**2.d0 - 1.d0)
      x3 = dsqrt((x + 1.d0)/(x - 1.d0))
!
      ch_f = x1 + x2 - x3
 100  return
      end
!
!
!
      double precision function z_f(g, g1, x)
      implicit none
      double precision g, g1, x
!
      double precision z, y, I1, I2
!
      if (x.lt.1.00000001d0) then
         z_f = 0.d0
         goto 100
      endif
      y = x**2.d0 - 1.d0
!
      I1 = dsqrt(y) - dlog(x + dsqrt(y)) + dsqrt((x - 1.d0)/(x + 1.d0))
      I2 = 5.d-1*(x*dsqrt(y) + dlog(x + dsqrt(y)))
!
      z_f = 5.d-1*((g + g1)**2.d0)*I1 - I2
 100  return
      end
!
!
      double precision function Inteta(x0, x1, p, q, tau)
      implicit none
      double precision x0, x1, p, q, tau
!
      double precision x, xs, dx, s, d, sum, sd, y
!
      sum = 0.d0
      x = x0
      d = dmin1((5.d-3*(x1 - x0)/x0), 5.d-3)
      s = 1.d0 + 5.d-1*d
      dx = 1.d0 + d
!
  100 xs = x*s
      y = tau - q*xs
      if (y.gt.-200.d0) then
         if (p.lt.0.1d0) then
           sd = dexp(y)
         else 
           sd = (xs**p)*dexp(y)
         endif
      else
         sd = 0.d0
      endif
      sum = sum + x*d*sd
      x = x*dx
      if ((x*s).lt.x1) goto 100
!
      Inteta = sum
      return
      end
!
!
!
!
!========================================================================================
      subroutine tridag
      implicit none
      include 'general.pa'
!
!
      integer i, n
      double precision bet, gam(num_nt)
      double precision a_i(num_nt), b_i(num_nt), c_i(num_nt),&
     &                 f_old(num_nt), f_new(num_nt), &
     &                 npos_old(num_nt), npos_new(num_nt)
      common / trid_update / a_i, b_i, c_i, f_old, f_new, npos_old,&
     &                npos_new
!
      if (dabs(b_i(1)).le.1.d-100) then
          write(*,*) 'Error: b(1) = 0.'
         return
      endif
!
      bet = b_i(1)
      f_new(1) = f_old(1)/bet
!
      do 100 i = 2, num_nt 
         gam(i) = c_i(i-1)/bet
         bet = b_i(i) - a_i(i)*gam(i)
         if (dabs(bet).le.1.d-100) then
             write(*,*) 'Error: bet = 0.'
            do 50 n = 1, num_nt
  50        f_new(n) = 0.d0
            return
         endif
!         if(f_old(i).lt.1.d-10)f_old(i)=1.d-10 ! Xuhui
         f_new(i) = (f_old(i) - a_i(i)*f_new(i-1))/bet 
!         if (f_new(i).lt.0.d0) f_new(i) = 0.d0  ! Xuhui
 100  continue
!
      do 200 i = num_nt-1, 1, -1
         f_new(i) = f_new(i) - gam(i+1)*f_new(i+1)
!         if (abs(f_new(i+1)).lt.1.d-8) f_new(i+1) = 0.d0 ! Xuhui
         if (f_new(i+1).lt.0.d0) f_new(i+1) = 0.d0 ! Xuhui
 200  continue
!
      return
      end
!
!
!
!
!============================================================================================
      subroutine trid_p
      implicit none
      include 'general.pa'
!
!
      integer i, n
      double precision bet, gam(num_nt)
      double precision a_i(num_nt), b_i(num_nt), c_i(num_nt),&
     &                 f_old(num_nt), f_new(num_nt), &
     &                 npos_new(num_nt), npos_old(num_nt)
      common / trid_update / a_i, b_i, c_i, f_old, f_new, npos_old,&
     &                npos_new
!
      if (dabs(b_i(1)).le.1.d-100) then
          write(*,*) 'Error: b(1) = 0.'
         return
      endif
!
      bet = b_i(1)
      npos_new(1) = npos_old(1)/bet
!
      do 100 i = 2, num_nt
         gam(i) = c_i(i-1)/bet
         bet = b_i(i) - a_i(i)*gam(i)    
         if (dabs(bet).le.1.d-100) then
             write(*,*) 'Error: bet = 0.'
            do 50 n = 1, num_nt
  50        npos_new(n) = 0.
            return
         endif
         npos_new(i) = (npos_old(i) - a_i(i)*npos_new(i-1))/bet
 100  continue
!
      do 200 i = num_nt-1, 1, -1
         npos_new(i) = npos_new(i) - gam(i+1)*npos_new(i+1)
         if (npos_new(i).lt.0.d0) npos_new(i) = 0.d0
         if (npos_new(i+1).lt.0.d0) npos_new(i+1) = 0.d0 ! Xuhui
 200  continue
!
      return
      end


!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:35:43 EDT 2006
! version: 2
! Name: J. Finke
! Added MPI_REDUCE of 'ecens'. Changed various 'write' statements. 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Thu Jun 22 11:23:38 EDT 2006
! version: 3
! Name: J. Finke
! Added summing of 'n_field' to subroutine 'E_add_up'  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Wed Jul  5 12:23:31 EDT 2006
! version: 4
! Name: J. Finke
! Create routine cens_add_up, which adds up items related 
! to the census files needed for the FP 
! calculation. It takes the place of E_add_up. E_add_up 
! was moved to after the FP_calculation. E_add_up now 
! also determines the maximum 'dT_max', as well as 
! adding up 'E_tot_new' and 'E_tot_old'.   
