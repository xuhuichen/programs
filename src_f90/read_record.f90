!     eps_a, double_ans are not recorded
      subroutine read_record
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer nunit_misc, nunit_census
      integer i, j, k, l
      character *43 localrecord, recordfilename

        nunit_census=myid+20000
        nunit_misc = myid+10000

!     master node part
!      if(myid.eq.master)then
        localrecord = recordfilename('misc.dat',myid)
        open(unit=nunit_misc,file=localrecord)
!
!       integral variables
        read(nunit_misc,'(i10)')ndxout
        read(nunit_misc,'(5i10)')cf_sentinel,nr,nz,ncycle,pair_switch

        do j=1,nz
          read(nunit_misc,'(2i10)')nsurfo(j), nsurfi(j)
        enddo
        do k=1,nr
          read(nunit_misc,'(2i10)')nsurfu(k), nsurfl(k)
        enddo
       do j=1,nz
         do k=1,nr
           read(nunit_misc,'(3i10)')nsv(j,k), npcen(j, k)
         enddo
       enddo

       read(nunit_misc,'(2i10)')rseed, rand_switch
       do j=1,nz
         do k=1,nr
           read(nunit_misc,'(i10)')seeds(j,k)
         enddo
       enddo
       do j=1,nz
         read(nunit_misc,'(i10)')zseeds(j)
       enddo
       do k=1,nr
         read(nunit_misc,'(i10)')rseeds(k)
       enddo

       read(nunit_misc,'(i10)')fp_sw
       read(nunit_misc,'(2i10)')nunit_census
       read(nunit_misc,'(i10)')nunit_evt
       read(nunit_misc,'(i10)')nmu
       do i=1,num_nt
         read(nunit_misc,'(i10)')nelectron(i)
       enddo

       read(nunit_misc,'(i10)')nphreg 
       do i=1,nregmax
         read(nunit_misc,'(i10)')nphbins(i)
       enddo

       read(nunit_misc,'(2i10)')nphtotal, nph_lc
       read(nunit_misc,'(i10)')cr_sent
       read(nunit_misc,'(i10)')nst
       read(nunit_misc,'(i10)')ntime
       do j=1,nz
         do k=1,nr
           read(nunit_misc,'(i10)')ep_switch(j,k)
         enddo
       enddo
       read(nunit_misc,'(i10)')turb_sw
       read(nunit_misc,'(i10)')star_switch
       read(nunit_misc,'(i10)')spec_switch
       read(nunit_misc,'(4i10)')inj_switch,g2var_switch,esc_sw,pick_sw
       read(nunit_misc,'(2i10)')FP_bcast_type, FP_Bcast_type2
       read(nunit_misc,'(i10)')ti
       do k=1,nr
         read(nunit_misc,'(i10)')npsurfr(k)
       enddo

       read(nunit_misc,'(i10)')vol_job_type
       do i=1,int_ans_size
         read(nunit_misc,'(i10)')integer_ans(i)
       enddo

       read(nunit_misc,'(i10)')nfile

       do j=1,nz
         read(nunit_misc,'(i10)')npsurfz(j)
       enddo

       read(nunit_misc,'(i10)')randcounter
       read(nunit_misc,'(i10)')idead
       read(nunit_misc,'(2i10)')split2, split3
       read(nunit_misc,'(2e14.7)')split1, spl3_trg
!
!       double precision variables
      read(nunit_misc,'(8e14.7)') r_flare, z_flare, t_flare, sigma_r, &
     &                 sigma_z, sigma_t, flare_amp, acc_prob
      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(6e14.7)') amxwl(j,k), gmin(j,k),&
     &                 gmax(j,k), p_nth(j,k),&
     &                 gbar_nth(j,k), N_nth(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
          do i=1,num_nt
          read(nunit_misc,'(2e14.7)') dne_pa(j,k,i),dnp_pa(j,k,i)
          enddo
        enddo
      enddo

      do j=1,nz
        do k=1,nr
          do i=1,num_nt
          read(nunit_misc,'(e14.7)') dn_pp(j,k,i)
          enddo
        enddo
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(3e14.7)')q_turb(j,k),turb_lev(j,k),&
     &       r_acc(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') ec_old(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
          do i=1,num_nt
          read(nunit_misc,'(e14.7)') n_pos(j,k,i)
          enddo
        enddo
      enddo

      do j=1,nz
        read(nunit_misc,'(2e14.7)') erini(j), erino(j)
      enddo
      do k=1,nr
        read(nunit_misc,'(2e14.7)') erinu(k), erinl(k)
      enddo

      do j=1,nz
        read(nunit_misc,'(2e14.7)') erlki(j), erlko(j)
      enddo
      do k=1,nr
        read(nunit_misc,'(2e14.7)') erlku(k), erlkl(k)
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') ecens(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(2e14.7)') edep(j,k), prdep(j,k)
        enddo
      enddo

      do i=1,n_vol
        do j=1,nz
          do k=1,nr
          read(nunit_misc,'(3e14.7)') kappa_tot(i,j,k), eps_tot(i,j,k), &
     &                             eps_th(i,j,k)
          enddo
        enddo
      enddo
      do i=1,n_vol
       read(nunit_misc,'(e14.7)') E_ph(i)
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(5e14.7)') Eloss_tot(j,k),&
     &                 Eloss_br(j,k), Eloss_cy(j,k), &
     &                 Eloss_sy(j,k), Eloss_th(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(4e14.7)') tea(j,k), tna(j,k),&
     &                            B_input(j,k), n_e(j,k), sum_E(j,k)
        enddo
      enddo

      do j=1,nz
        read(nunit_misc,'(e14.7)') z(j)
      enddo 
      do k=1,nr
        read(nunit_misc,'(e14.7)') r(k)
      enddo
      read(nunit_misc,'(2e14.7)')rmin, zmin
      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(2e14.7)')vol(j,k), zsurf(j,k)
        enddo
      enddo

      do k=1,nr
        read(nunit_misc,'(2e14.7)') Asurfu(k), Asurfl(k)
      enddo
      do j=1,nz
        read(nunit_misc,'(2e14.7)')Asurfi(j), Asurfo(j)
      enddo

      read(nunit_misc,'(5e14.7)') time, dt(1), dt(2), tstop, mcdt
      do i=1,ntime
        read(nunit_misc,'(2e14.7)') t0(i), t1(i)
      enddo

      do j=1,nz
        do k=1,nr      
        read(nunit_misc,'(e14.7)') f_pair(j,k)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
          do i=1,num_nt
          read(nunit_misc,'(e14.7)') Pnt(j,k,i)
          enddo
        enddo
      enddo
      do i=1,num_nt
         read(nunit_misc,'(e14.7)') gnt(i)
      enddo

      do j=1,nz
        do k=1,nr
          do i=1,num_nt
          read(nunit_misc,'(e14.7)') f_nt(j,k,i)
          enddo
        enddo
      enddo
      read(nunit_misc,'(e14.7)') dg

      do i=1,num_nt
        read(nunit_misc,'(e14.7)') E_IC(i)
      enddo

      do i=1,nphfield
        do j=1,nz
          do k=1,nr
          read(nunit_misc,'(e14.7)') n_field(i,j,k) 
          enddo
        enddo
      enddo
      do i=1,nphfield
        read(nunit_misc,'(e14.7)') E_field(i)
      enddo

      do i=1,n_gg
        do j=1,nz
          do k=1,nr
          read(nunit_misc,'(2e14.7)') n_ph(i,j,k),k_gg(i,j,k)
          enddo
        enddo
      enddo
      do i=1,n_gg
        read(nunit_misc,'(e14.7)') E_gg(i)
      enddo

      do i=1,nmu
        read(nunit_misc,'(e14.7)') mu(i)
      enddo

      do i=1,n_ref
        do j=1,n_ref
        read(nunit_misc,'(2e14.7)') P_ref(i, j), W_abs(i,j)
        enddo
      enddo
      do i=1,n_ref
        read(nunit_misc,'(e14.7)') E_ref(i)
      enddo
      read(nunit_misc,'(6e14.7)') F_ib, dE_in, ref_delay,&
     &                             f_ref, A_in, n_disk

      do k=1,nr
        read(nunit_misc,'(3e14.7)') Ed_abs(k), Ed_ref(k), Ed_in(k)
      enddo

      do i=1,nregmax
        read(nunit_misc,'(2e14.7)') Ephmin(i), Ephmax(i)
      enddo
      do i=1,nphomax+1
        read(nunit_misc,'(e14.7)') hu(i)
      enddo
      do i=1,nphlcmax
        read(nunit_misc,'(2e14.7)')Elcmin(i), Elcmax(i)
      enddo

      do k=1,nr
        do i=1,ntime
        read(nunit_misc,'(2e14.7)') tbbu(k,i), tbbl(k,i)
        enddo
      enddo
      do j=1,nz
        do i=1,ntime
        read(nunit_misc,'(2e14.7)')tbbi(j, i), tbbo(j,i)
        enddo
      enddo
      read(nunit_misc,'(2e14.7)') Rstar, dist_star
      read(nunit_misc,'(9e14.7)') inj_g1, inj_g2, inj_p, inj_t,inj_gg,&
     &                 inj_sigma, inj_L, inj_v, g_bulk
      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') pick_rate(j, k)
        enddo
      enddo
      do i=1,num_nt
        read(nunit_misc,'(e14.7)') pick_dis(i)
      enddo
      read(nunit_misc,'(e14.7)') r_esc
      read(nunit_misc,'(4e14.7)') R_blr, fr_blr, R_ir, fr_ir, &
     &                 R_disk, d_jet

      do i=1,num_nt
        do j=1,nphfield
        read(nunit_misc,'(e14.7)') F_IC(i, j)
        enddo
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') Te_new(j, k)
        enddo
      enddo
      read(nunit_misc,'(3e14.7)') dt_new, lnL, dT_max

      read(nunit_misc,'(7e14.7)') E_tot_old, E_tot_new, hr_total,&
     &                 hr_st_total, f_t_implicit, dr, dz

      do i=1,num_nt
        do j=1,num_temp_max
        read(nunit_misc,'(4e14.7)') lib_dg_ce(i,j), lib_dg_cp(i,j),&
     &                 lib_disp_ce(i,j), lib_disp_cp(i,j)
        enddo
      enddo
      do i=1,num_nt
        read(nunit_misc,'(e14.7)') rate_gm1(i) 
      enddo
      read(nunit_misc,'(4e14.7)')tea_min, tea_max, tna_min, tna_max

      do k=1,nr
        read(nunit_misc,'(2e14.7)') ewsurfu(k), ewsurfl(k)
      enddo
      do j=1,nz
        read(nunit_misc,'(2e14.7)') ewsurfo(j), ewsurfi(j)
      enddo
      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') ewsv(j,k)
        enddo
      enddo

      do i=1,double_ans_size
        read(nunit_misc,'(e14.7)') double_ans(i)
      enddo

      do i=1,nfmax
        read(nunit_misc,'(6e14.7)') E_file(i), I_file(i), P_file(i),&
     &                 F_file(i), alpha(i), a1(i), int_file
      enddo

      do i=1,201
        read(nunit_misc,'(e14.7)') comp0(i) 
      enddo
      do i=1,13
        do j=1,66
        read(nunit_misc,'(e14.7)')enxtab(i,j)
        enddo
      enddo
      do i=1,26
        do j=1,6
          do k=1,8
            do l=1,66
            read(nunit_misc,'(e14.7)') enx_nth(i,j,k,l)
            enddo
          enddo
        enddo
      enddo

      do i=1,randmax
        read(nunit_misc,'(e14.7)') randlist(i)
      enddo

      do j=1,nz
        do k=1,nr
        read(nunit_misc,'(e14.7)') T_sum(j,k)
        enddo
      enddo
      read(nunit_misc,'(e14.7)') time_sum

      do i=1,nmu
        do j=1,nphlcmax
        read(nunit_misc,'(e14.7)') edout(i,j)
        enddo
      enddo
      do i=1,nmu
        do j=1,nphomax
        read(nunit_misc,'(e14.7)') fout(i,j)
        enddo
      enddo

      read(nunit_misc,'(3e14.7)') fac_old, Emiss_tot, Emiss_old
      read(nunit_misc,'(e14.7)') t_bound
!
!     characters
      read(nunit_misc,'(a30)') temp_file
      read(nunit_misc,'(3a30)') spname, phname, eventfile
      do i=1,nmu
        read(nunit_misc,'(a30)')lcname(i)
      enddo
      do j=1,nz
        do i=1,ntime
        read(nunit_misc,'(2a30)') i_fname(j,i), o_fname(j,i)
        enddo
      enddo
      do k=1,nr
        do i=1,ntime
        read(nunit_misc,'(2a30)')u_fname(k,i), l_fname(k,i)
        enddo
      enddo

        close(nunit_misc)

!     slave node part
!      else
        localrecord = recordfilename('census.dat',myid)
        open(unit=nunit_census, file=localrecord)
        call read_cens(ndxout,nunit_census)
        close(nunit_census)
!      endif

      call make_FP_bcast_type

      end


