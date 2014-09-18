!     common blocks used in updated.f, imctrk2d.f, imcgen.f, xec2d.f, record.f
      integer ndxout
      integer cf_sentinel, nr, nz, ncycle, pair_switch
      integer ierr, myid, numprocs, master
      integer jph, kph ,jgpsp, jgplc, jgpmu, emitype
      integer nsurfu(kmax), nsurfl(kmax), nsurfo(jmax), nsurfi(jmax),&
     &        nsv(jmax,kmax), npcen(jmax, kmax)
      integer ibufin(iucens), ibufout(iucens)
      integer rseed, seeds(jmax,kmax), zseeds(jmax), &
     &     rseeds(kmax), rand_switch
      integer fp_sw
      integer incounter, outcounter,cmcount1,cmcount2,lkcount
      integer nunit_evt
      integer nmu
      integer nelectron(num_nt)
      integer nphreg, nphbins(nregmax)
      integer nphtotal, nph_lc
      integer cr_sent
      integer nst
      integer ntime
      integer ep_switch(jmax,kmax)
      integer turb_sw
      integer star_switch
      integer spec_switch
      integer inj_switch, g2var_switch, esc_sw, pick_sw
      integer FP_bcast_type, FP_Bcast_type2
      integer ti
      integer npsurfr(kmax)
      integer vol_job_type, integer_ans(int_ans_size)
      integer nfile
      integer npsurfz(jmax)
      integer randcounter
      integer idead
      integer split2, split3
      double precision split1, spl3_trg
      double precision r_flare, z_flare, t_flare, sigma_r, &
     &                 sigma_z, sigma_t, flare_amp
      double precision acc_prob
      double precision amxwl(jmax, kmax), gmin(jmax, kmax),&
     &                 gmax(jmax, kmax), p_nth(jmax, kmax),&
     &                 gbar_nth(jmax, kmax), N_nth(jmax, kmax)
      double precision dne_pa(jmax, kmax, num_nt), &
     &                 dnp_pa(jmax, kmax, num_nt)
      double precision dn_pp(jmax, kmax, num_nt)
      double precision q_turb(jmax, kmax), turb_lev(jmax,kmax)
      double precision ec_old(jmax,kmax)
      double precision n_pos(jmax,kmax,num_nt)
      double precision erini(jmax), erino(jmax), &
     &                 erinu(kmax), erinl(kmax)
      double precision erlki(jmax), erlko(jmax), &
     &                 erlku(kmax), erlkl(kmax)
      double precision ecens(jmax,kmax)
      double precision edep(jmax,kmax), prdep(jmax,kmax)
      double precision kappa_tot(n_vol, jmax, kmax),&
     &                 eps_tot(n_vol, jmax, kmax),&
     &                 eps_a(num_a, n_vol, jmax, kmax),&
     &                 eps_th(n_vol, jmax, kmax), E_ph(n_vol)
      double precision Eloss_tot(jmax,kmax) ,Eloss_br(jmax,kmax),  &
     &                 Eloss_cy(jmax,kmax), &
     &                 Eloss_sy(jmax, kmax), Eloss_th(jmax, kmax)
      double precision tea(jmax,kmax),tna(jmax,kmax),B_input(jmax,kmax),&
     &                 B_field(jmax,kmax), n_e(jmax, kmax),&
     &                 sum_E(jmax,kmax)
      double precision z(jmax), r(kmax), rmin, zmin, vol(jmax,kmax),&
     &                 zsurf(jmax,kmax)
      double precision Asurfu(kmax), Asurfl(kmax), &
     &                 Asurfi(jmax), Asurfo(jmax)
      double precision time, dt(2), tstop, mcdt, t0(ntmax), t1(ntmax)
      double precision f_pair(jmax,kmax)
      double precision Pnt(jmax,kmax,num_nt), gnt(num_nt)
      double precision f_nt(jmax,kmax,num_nt), dg
      double precision E_IC(num_nt)
      double precision n_field(nphfield,jmax,kmax), E_field(nphfield)
      double precision xnu, wmu, phi, rpre, zpre, dcen, ew
      double precision dbufin(ducens), dbufout(ducens)
      double precision n_ph(n_gg, jmax, kmax), k_gg(n_gg, jmax, kmax), &
     &                 E_gg(n_gg)    
      double precision mu(nmumax) 
      double precision P_ref(n_ref, n_ref), W_abs(n_ref,n_ref),&
     &                 E_ref(n_ref), F_ib, dE_in, ref_delay, &
     &                 f_ref, A_in, n_disk
      double precision Ed_abs(kmax), Ed_ref(kmax), Ed_in(kmax)
      double precision Ephmin(nregmax), Ephmax(nregmax)
      double precision hu(nphomax+1), Elcmin(nphlcmax), Elcmax(nphlcmax)
      double precision tbbu(kmax, ntmax), tbbl(kmax, ntmax), &
     &                 tbbi(jmax, ntmax), tbbo(jmax, ntmax)
      double precision Rstar, dist_star
      double precision inj_g1, inj_g2, inj_p, inj_t, inj_gg,&
     &                 inj_sigma, inj_L, inj_v, g_bulk,&
     &                 pick_rate(jmax,kmax), pick_dis(num_nt)
      double precision theta_b(jmax,kmax), theta_local(jmax,kmax)
      double precision r_esc, r_acc_peak, r_acc(jmax,kmax)
      double precision R_blr, fr_blr, R_ir, fr_ir, R_disk, d_jet 
      double precision F_IC(num_nt, nphfield)
      double precision Te_new(jmax, kmax),dt_new, lnL, dT_max
      double precision E_tot_old, E_tot_new, hr_total,&
     &                 hr_st_total, f_t_implicit, dr, dz
      double precision lib_dg_ce(num_nt,num_temp_max),&
     &                 lib_dg_cp(num_nt, num_temp_max),&
     &                 lib_disp_ce(num_nt, num_temp_max),&
     &                 lib_disp_cp(num_nt, num_temp_max),&
     &                 rate_gm1(num_nt), tea_min, tea_max,&
     &                 tna_min, tna_max
      double precision ewsurfu(kmax), ewsurfl(kmax), ewsurfo(jmax), &
     &                 ewsurfi(jmax), ewsv(jmax, kmax)
      double precision double_ans(double_ans_size)
      double precision E_file(nfmax), I_file(nfmax), P_file(nfmax),&
     &                 F_file(nfmax), alpha(nfmax), a1(nfmax), int_file
      double precision comp0(201), enxtab(13,66),&
     &                 enx_nth(26,6,8,66)
      double precision randlist(randmax)
      double precision T_sum(jmax,kmax), time_sum
      double precision edout(nmumax, nphlcmax),&
     &                 fout(nmumax, nphomax)
      double precision fac_old, Emiss_tot, Emiss_old
      double precision t_bound
      double precision F_sync_x(36), F_sync_t(36)
!     particle diffusion
      double precision n_diff(jmax,kmax,num_nt),diff_cf(num_nt),sweep
      double precision currenttimestop
      character *30 temp_file
      character *30 spname, phname, eventfile, lcname(nmumax)
      character *30 i_fname(jmax, ntmax), o_fname(jmax, ntmax),&
     &              u_fname(kmax, ntmax), l_fname(kmax, ntmax)
      real etotal, etotal_old

      common / ndx / ndxout
      common / i_flare / cf_sentinel
      common / nc / ncycle
      common / izones / nz, nr
      common / ps / pair_switch
      ! MPI block not to be recorded
      common / MPI / ierr, myid, numprocs, master
      ! iphoton block not to be recorded
      common / iphoton / jph, kph, jgpsp, jgplc, jgpmu, emitype
      common / ph_numbers / nsurfu, nsurfl, nsurfo, nsurfi, nsv, &
     &                      npcen
      common / ibuffer / ibufin, ibufout
      common / random / rseed, seeds, zseeds, rseeds, rand_switch
      common / tc / fp_sw
!     counter block not to be recorded
      common / counter / incounter,outcounter,cmcount1,cmcount2,lkcount
      common / event_pointer / nunit_evt
      common / imubins / nmu
      common / nel / nelectron
      common / ieb_setup / nphreg, nphbins
      common / ienergy / nphtotal, nph_lc
      common / icr / cr_sent
      common / phnumber / nst
      common / itimes / ntime
      common / izq / ep_switch
      common / turb / turb_sw
      common / star_sw / star_switch
      common / spec_sw / spec_switch
      common / injswi / inj_switch, g2var_switch, esc_sw, pick_sw
      common / Bcast / FP_bcast_type, FP_bcast_type2
      common / timeindex / ti
      common / imcsurfr / npsurfr
      common / vol_var_type  / vol_job_type, integer_ans
      common / inspi / nfile
      common / imcsurfz / npsurfz
      common / rcount / randcounter
      common / itrk / idead
      common / split / split2, split3
      common / split_trg / split1, spl3_trg

      common / d_flare / r_flare, z_flare, t_flare, sigma_r, &
     &                   sigma_z, sigma_t, flare_amp
      common / nontherm / amxwl, gmin, gmax, p_nth, gbar_nth, N_nth
      common / annihil / dne_pa, dnp_pa
      common / dnpp / dn_pp
      common / turbulence / q_turb, turb_lev
      common / ecold / ec_old
      common / n_pair / n_pos
      common / s_energies / erini, erino, erinu, erinl
      common / s_leakage / erlki, erlko, erlku, erlkl
      common / cens_energy / ecens
      common / deposition / edep, prdep
      common / vol_em / kappa_tot, eps_tot, eps_a, eps_th, E_ph, &
     &                 Eloss_tot, Eloss_br, Eloss_cy, Eloss_sy, Eloss_th
      common / zone_quantities / tea, tna, B_input, B_field, n_e, sum_E
      common / zones / z, r, rmin, zmin, vol, Asurfu, Asurfl,&
     &                 Asurfi, Asurfo, zsurf
      common / times / time, dt, tstop,currenttimestop, mcdt, t0, t1
      common / fpair / f_pair
      common / Pnth / Pnt, gnt, f_nt, dg
      common / IC_track / E_IC 
      common / photonfield / n_field, E_field
!     7 photon variables do not need to be recorded
      common / photon / xnu, wmu, phi, rpre, zpre, dcen, ew
      common / dbuffer / dbufin, dbufout
      common / gg_abs / n_ph, k_gg, E_gg
      common / mubins / mu
      common / reflection / P_ref, W_abs, E_ref, F_ib, dE_in,&
     &                      ref_delay, f_ref, A_in, n_disk
      common / ddh / Ed_abs, Ed_ref, Ed_in
      common / eb_setup / Ephmin, Ephmax
      common / energy / hu, Elcmin, Elcmax
      common / bdtemp / tbbu, tbbl, tbbi, tbbo
      common / star / Rstar, dist_star
      common / inj / inj_g1, inj_g2, inj_p, inj_t, inj_gg,&
     &               inj_sigma, inj_L, inj_v, g_bulk, pick_rate, pick_dis
      common / anis / theta_b, theta_local
      common / escacc / r_esc, r_acc_peak, r_acc, acc_prob
      common / extrad / R_blr, fr_blr, R_ir, fr_ir, R_disk, d_jet
      common / fic / F_IC
      common / d_update / Te_new, dt_new, lnL, dT_max
      common / to_fp_calc / E_tot_old, E_tot_new, hr_total, &
     &                      hr_st_total, f_t_implicit, dr, dz
!     these coulomb_sc not actually used
      common / coulomb_sc / lib_dg_ce, lib_disp_ce, lib_dg_cp, &
     &                      lib_disp_cp, rate_gm1, tea_min, tea_max,&
     &                      tna_min, tna_max
      common / eweights / ewsurfu, ewsurfl, ewsurfi, ewsurfo, ewsv
      common / vol_var_type_dbl / double_ans
      common / insp / E_file, I_file, P_file, F_file, alpha, a1,int_file
      common / ctot / comp0, enxtab, enx_nth
      common / rlist / randlist
      common / T_sums / T_sum, time_sum
      common / outputs / edout, fout
      common / fac_o / fac_old, Emiss_tot, Emiss_old
      common / t_esc / t_bound
      common / sync_table / F_sync_x, F_sync_t
      common / diffusion / n_diff,diff_cf,sweep

      common / tf / temp_file
      common / filenames / spname, phname, eventfile, lcname
      common / s_fnames / i_fname, o_fname, u_fname, l_fname

!     does not record elapse time. It is process dependent
      common / elapse_time / etotal, etotal_old

















