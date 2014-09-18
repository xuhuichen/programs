!      
!
!
!
!
!     This broadcasts the names of the event file (eventfile)
!     to the slaves, as well as a few other things.
!     J. Finke, 19 July 2005
!     updated again, name_bcast (in xec2d) and setup_bcast combined to make
!     new setup_bcast.
!     J. Finke, 7 February 2006.
      subroutine setup_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer i
      double precision junk, ran1
      integer start_address, address
      integer displacements(64), types(64), lengths(64)
      integer setup_bcast_type
!
!
!
!
!     broadcast filenames names
      call MPI_BCAST(eventfile, 30, MPI_CHARACTER, master,&
     &     MPI_COMM_WORLD, ierr)
!
!     broadcast tc
      call MPI_BCAST(fp_sw, 1, MPI_INTEGER, master, &
     &      MPI_COMM_WORLD, ierr)
!
!     broadcast times
      call MPI_BCAST(tstop, 1, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(currenttimestop, 1, MPI_DOUBLE_PRECISION, master, &
     &      MPI_COMM_WORLD, ierr)
!
!     imubins
      call MPI_BCAST(nmu, 1, MPI_INTEGER, master,&
     &     MPI_COMM_WORLD, ierr)
!
!     zone_quantities
      call MPI_BCAST(n_e, jmax*kmax, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(B_input,jmax*kmax,MPI_DOUBLE_PRECISION,&
     &     master,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea,jmax*kmax,MPI_DOUBLE_PRECISION,&
     &     master,MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna,jmax*kmax,MPI_DOUBLE_PRECISION,&
     &     master,MPI_COMM_WORLD, ierr)
!
!     icr
      call MPI_BCAST(cr_sent, 1, MPI_INTEGER, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(acc_prob, 1, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!
!     star_switch
      call MPI_BCAST(star_switch, 1, MPI_INTEGER, master, &
     &     MPI_COMM_WORLD, ierr)
!
!     common block random
      call MPI_BCAST(rand_switch, 1, MPI_INTEGER, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rseed, 1, MPI_INTEGER, master,&
     &     MPI_COMM_WORLD, ierr)
      rseed = rseed + myid*84725
      if(rand_switch.eq.2) then
         do 10 i=1,myid+10
            junk = ran1(rseed)
!            write(*,*) i,' myid=',myid,' junk=',junk,' rseed=',rseed
 10      continue
      endif
!
!     common block photonfield
      call MPI_BCAST(E_field, nphfield, MPI_DOUBLE_PRECISION, master, &
     &     MPI_COMM_WORLD, ierr)
!
!     common block Pnth and npos
      call MPI_BCAST(Pnt,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gnt, num_nt, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_nt,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(n_pos,jmax*kmax*num_nt,MPI_DOUBLE_PRECISION,master,&
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dg, 1, MPI_DOUBLE_PRECISION, master,&
     &     MPI_COMM_WORLD, ierr)
!
!     integers
!
!     common block i_flare
      call MPI_ADDRESS(cf_sentinel, start_address, ierr)
      displacements(1) = 0
      lengths(1) = 1
      types(1) = MPI_INTEGER

!     common block izones
      call MPI_ADDRESS(nr, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = 1
      types(2) = MPI_INTEGER

      call MPI_ADDRESS(nz, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_INTEGER

!     common block ps
      call MPI_ADDRESS(pair_switch, address, ierr)
      displacements(4) = address - start_address
      lengths(4) = 1
      types(4) = MPI_INTEGER
!
!     doubles
!
!     common block zones.
      call MPI_ADDRESS(z, address, ierr)
      displacements(5) = address - start_address
      lengths(5) = jmax
      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(r, address, ierr)
      displacements(6) = address - start_address
      lengths(6) = kmax
      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(rmin, address, ierr)
      displacements(7) = address - start_address
      lengths(7) = 1
      types(7) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(zmin, address, ierr)
      displacements(8) = address - start_address
      lengths(8) = 1
      types(8) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(vol, address, ierr)
      displacements(9) = address - start_address
      lengths(9) = kmax*jmax
      types(9) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(zsurf, address, ierr)
      displacements(10) = address - start_address
      lengths(10) = kmax*jmax
      types(10) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfu, address, ierr)
      displacements(11) = address - start_address
      lengths(11) = kmax
      types(11) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfl, address, ierr)
      displacements(12) = address - start_address
      lengths(12) = kmax
      types(12) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfo, address, ierr)
      displacements(13) = address - start_address
      lengths(13) = jmax
      types(13) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(Asurfi, address, ierr)
      displacements(14) = address - start_address
      lengths(14) = kmax
      types(14) = MPI_DOUBLE_PRECISION

!     common block d_flare
      call MPI_ADDRESS(r_flare, address, ierr)
      displacements(15) = address - start_address
      lengths(15) = 1
      types(15) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(z_flare, address, ierr)
      displacements(16) = address - start_address
      lengths(16) = 1
      types(16) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t_flare, address, ierr)
      displacements(17) = address - start_address
      lengths(17) = 1
      types(17) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_r, address, ierr)
      displacements(18) = address - start_address
      lengths(18) = 1
      types(18) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_z, address, ierr)
      displacements(19) = address - start_address
      lengths(19) = 1
      types(19) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sigma_t, address, ierr)
      displacements(20) = address - start_address
      lengths(20) = 1
      types(20) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(flare_amp, address, ierr)
      displacements(21) = address - start_address
      lengths(21) = 1
      types(21) = MPI_DOUBLE_PRECISION
!
!     common block inj
      call MPI_ADDRESS(inj_g1, address, ierr)
      displacements(22) = address - start_address
      lengths(22) = 1
      types(22) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_g2, address, ierr)
      displacements(23) = address - start_address
      lengths(23) = 1
      types(23) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_p, address, ierr)
      displacements(24) = address - start_address
      lengths(24) = 1
      types(24) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_t, address, ierr)
      displacements(25) = address - start_address
      lengths(25) = 1
      types(25) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_v, address, ierr)
      displacements(26) = address - start_address
      lengths(26) = 1
      types(26) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_gg, address, ierr)
      displacements(27) = address - start_address
      lengths(27) = 1
      types(27) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_sigma, address, ierr)
      displacements(28) = address - start_address
      lengths(28) = 1
      types(28) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_switch, address, ierr)
      displacements(29) = address - start_address
      lengths(29) = 1
      types(29) = MPI_INTEGER
      call MPI_ADDRESS(esc_sw, address, ierr)
      displacements(30) = address - start_address
      lengths(30) = 1
      types(30) = MPI_INTEGER
      call MPI_ADDRESS(g2var_switch, address, ierr)
      displacements(31) = address - start_address
      lengths(31) = 1
      types(31) = MPI_INTEGER
      call MPI_ADDRESS(r_esc, address, ierr)
      displacements(32) = address - start_address
      lengths(32) = 1
      types(32) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(r_acc, address, ierr)
      displacements(33) = address - start_address
      lengths(33) = 1
      types(33) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(inj_L, address, ierr)
      displacements(34) = address - start_address
      lengths(34) = 1
      types(34) = MPI_DOUBLE_PRECISION  
      call MPI_ADDRESS(pick_sw, address, ierr)
      displacements(35) = address - start_address
      lengths(35) = 1
      types(35) = MPI_INTEGER
      call MPI_ADDRESS(pick_rate, address, ierr)
      displacements(36) = address - start_address
      lengths(36) = jmax*kmax
      types(36) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(g_bulk, address, ierr)
      displacements(37) = address - start_address
      lengths(37) = 1
      types(37) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(R_blr, address, ierr)
      displacements(38) = address - start_address
      lengths(38) = 1
      types(38) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(fr_blr, address, ierr)
      displacements(39) = address - start_address
      lengths(39) = 1
      types(39) = MPI_DOUBLE_PRECISION  
      call MPI_ADDRESS(R_ir, address, ierr)
      displacements(40) = address - start_address
      lengths(40) = 1
      types(40) = MPI_DOUBLE_PRECISION   
      call MPI_ADDRESS(fr_ir, address, ierr)
      displacements(41) = address - start_address
      lengths(41) = 1
      types(41) = MPI_DOUBLE_PRECISION    
      call MPI_ADDRESS(R_disk, address, ierr)
      displacements(42) = address - start_address
      lengths(42) = 1
      types(42) = MPI_DOUBLE_PRECISION 
      call MPI_ADDRESS(d_jet, address, ierr)
      displacements(43) = address - start_address
      lengths(43) = 1
      types(43) = MPI_DOUBLE_PRECISION     
      call MPI_ADDRESS(split1, address, ierr)
      displacements(44) = address - start_address
      lengths(44) = 1
      types(44) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(split2, address, ierr)
      displacements(45) = address - start_address
      lengths(45) = 1
      types(45) = MPI_INTEGER
      call MPI_ADDRESS(split3, address, ierr)
      displacements(46) = address - start_address
      lengths(46) = 1
      types(46) = MPI_INTEGER
      call MPI_ADDRESS(spl3_trg, address, ierr)
      displacements(47) = address - start_address
      lengths(47) = 1
      types(47) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(theta_b, address, ierr)
      displacements(48) = address - start_address
      lengths(48) = kmax*jmax
      types(48) = MPI_DOUBLE_PRECISION
      call MPI_ADDRESS(pick_dis, address, ierr)
      displacements(49) = address - start_address
      lengths(49) = num_nt
      types(49) = MPI_DOUBLE_PRECISION

!
      call MPI_TYPE_STRUCT(49, lengths, displacements, types, &
     &     setup_bcast_type, ierr)
      call MPI_TYPE_COMMIT(setup_bcast_type, ierr)
      call MPI_BCAST(cf_sentinel, 1, setup_bcast_type,&
     &      master, MPI_COMM_WORLD, ierr)
!
!
!
      return
      end
!
!
!
!
!
!
!     This subroutine calculates the j and k (vertical and radial
!     zone numbers) from the zone number, zone.
      subroutine get_j_k(zone, j, k, r)
      implicit none
      integer zone, j, k, r
!
      j = (zone-1)/r+1
      k = zone - r*(j-1)
      return 
      end
!
!
!     
!     This calculates the zone number from j and k, the radial
!     and vertical zone numbers.
      integer function get_zone(j, k, r)
      implicit none
      integer j, k, r
!
      get_zone = k+r*(j-1)
      return 
      end
!
!
!
!
!
!     This routine makes the types used in FP_bcast to efficiently
!     broadcast to the processes.
!     J. Finke, 7 February 2006
      subroutine make_FP_bcast_type
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer status(MPI_STATUS_SIZE)

!     the above Asurfs are not used in fp_calc but are in the common
!     block with z, r etc.


!     Below are libraries from the subroutine coulomb.
!     variables for building data type FP_bcast_type
      integer start_address, address
      integer displacements(28), types(28), lengths(28)
      integer FP_Bcast_type3
!

!
!
!     variables common from outside update.
!
!      write(*,*) 'myid=', myid, ' in make_fp_bcast_type'
!

!     doubles from outside update.
!
!     common block times
!      call MPI_ADDRESS(time, start_address, ierr)
!      displacements(1) = 0
!      lengths(1) = 1
!      types(1) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(dt, address, ierr)
!      displacements(2) = address - start_address
!      lengths(2) = 2
!      types(2) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t0, start_address, ierr)
      displacements(1) = 0
      lengths(1) = ntmax
      types(1) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(t1, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = ntmax
      types(2) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(tstop, address, ierr)
!      displacements(5) = address - start_address
!      lengths(5) = 1
!      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(mcdt, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(sweep, address, ierr)
      displacements(4) = address - start_address
      lengths(4) = 1
      types(4) = MPI_DOUBLE_PRECISION

!     common block fic.
      call MPI_ADDRESS(F_IC, address, ierr)
      displacements(5) = address - start_address
      lengths(5) = num_nt*nphfield
      types(5) = MPI_DOUBLE_PRECISION

!      call MPI_TYPE_INDEXED(2, lengths, displacements, 
!     1     MPI_DOUBLE_PRECISION, FP_bcast_type, ierr)
      call MPI_TYPE_STRUCT(5, lengths, displacements, &
     &     types, FP_bcast_type, ierr)
      call MPI_TYPE_COMMIT(FP_bcast_type, ierr)

!     common block photonfield.
!      call MPI_ADDRESS(n_field, start_address, ierr)
!      displacements(1) = 0
!      lengths(1) = nphfield*jmax*kmax
!      types(1) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(E_field, start_address, ierr)
      displacements(1) = 0
      lengths(1) = nphfield
      types(1) = MPI_DOUBLE_PRECISION

!      call MPI_TYPE_STRUCT(2, lengths, displacements, 
!     1     types, FP_bcast_type2, ierr)
!      call MPI_TYPE_COMMIT(FP_bcast_type2, ierr)

!     common block n_pair
!      call MPI_ADDRESS(n_pos, start_address, ierr)
!      displacements(1) = 0
!      lengths(1) =  kmax*jmax*num_nt
!      types(1) = MPI_DOUBLE_PRECISION

!     common block Pnth
!      call MPI_ADDRESS(Pnt, address, ierr)
!      displacements(2) = address - start_address
!      lengths(2) = kmax*jmax*num_nt
!      types(2) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(gnt, address, ierr)
!      displacements(3) = address - start_address
!      lengths(3) = num_nt
!      types(3) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(f_nt, address, ierr)
!      displacements(4) = address - start_address
!      lengths(6) =  kmax*jmax*num_nt
!      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dg, address, ierr)
      displacements(2) = address - start_address
      lengths(2) = 1
      types(2) = MPI_DOUBLE_PRECISION

!     common block dnpp
!      call MPI_ADDRESS(dn_pp, address, ierr)
!      displacements(6) = address - start_address
!      lengths(8) =  kmax*jmax*num_nt
!      lengths(8) = 1
!      types(8) = MPI_DOUBLE_PRECISION

!     common block annihil
!      call MPI_ADDRESS(dne_pa, address, ierr)
!      displacements(7) = address - start_address
!      lengths(9) =  kmax*jmax*num_nt
!      types(9) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(dnp_pa, address, ierr)
!      displacements(8) = address - start_address
!      lengths(10) = kmax*jmax*num_nt
!      types(10) = MPI_DOUBLE_PRECISION
!
!     variables common only between update, FP_calc, and photon_fill
!
!     common block d_update
!      call MPI_ADDRESS(Te_new, address, ierr)
!      displacements(9) = address - start_address
!      lengths(11) = kmax*jmax
!      types(11) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dt_new, address, ierr)
      displacements(3) = address - start_address
      lengths(3) = 1
      types(3) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(lnL, address, ierr)
!      displacements(11) = address - start_address
!      lengths(13) = 1
!      types(13) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(dT_max, address, ierr)
      displacements(4) = address - start_address
      lengths(4) = 1
      types(4) = MPI_DOUBLE_PRECISION

!     common block to_fp_calc
      call MPI_ADDRESS(E_tot_old, address, ierr)
      displacements(5) = address - start_address
      lengths(5) = 1
      types(5) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(E_tot_new, address, ierr)
      displacements(6) = address - start_address
      lengths(6) = 1
      types(6) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(hr_st_total, address, ierr)
      displacements(7) = address - start_address
      lengths(7) = 1
      types(7) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(hr_total, address, ierr)
      displacements(8) = address - start_address
      lengths(8) = 1
      types(8) = MPI_DOUBLE_PRECISION

      call MPI_ADDRESS(f_t_implicit, address, ierr)
      displacements(9) = address - start_address
      lengths(9) = 1
      types(9) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(dr, address, ierr)
!      displacements(20) = address - start_address
!      lengths(20) = 1
!      types(20) = MPI_DOUBLE_PRECISION

!      call MPI_ADDRESS(dz, address, ierr)
!      displacements(21) = address - start_address
!      lengths(21) = 1
!      types(21) = MPI_DOUBLE_PRECISION

      call MPI_TYPE_STRUCT(9, lengths, displacements, &
     &     types, FP_bcast_type2, ierr)
      call MPI_TYPE_COMMIT(FP_bcast_type2, ierr)
!
!
!

!      call MPI_BCAST(time, 1, MPI_DOUBLE_PRECISION,
!     1      master, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(time, 1, FP_bcast_type, 
!     1      master, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(n_field, 1, FP_bcast_type2, 
!     1      master, MPI_COMM_WORLD, ierr)
!
!
      return
      end
!
!
!
!
!
!
!     This subroutine will broadcast to all of the processes
!     data needed for the Fokker-Planck routine.
!     J. Finke, 7 March 2005
      subroutine FP_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer status(MPI_STATUS_SIZE)
!     the above Asurfs are not used in fp_calc but are in the common
!     block with z, r etc.
!
!
!      write(*,*) 'myid=', myid, ' in fp_bcast'
!
!
!      call MPI_BCAST(time, 1, FP_bcast_type, 
!     1      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(t0, 1, FP_bcast_type, &
     &      master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(E_field, 1, FP_bcast_type2, &
     &      master, MPI_COMM_WORLD, ierr)
!
!
!     
      return
      end
!
!
!
!     FP_send_job sends the zone number 
!     from the master node to the available slave.  It
!     must be paired with FP_recv_job by the slave node.
!     J. Finke,  7 March 2005
      subroutine FP_send_job(available_proc, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer available_proc, zone
!
      integer ans_size, last
      parameter (last = 23)
      parameter (ans_size = num_nt*8+nphfield+last)
!
      integer j, k, i
      integer status(MPI_STATUS_SIZE)
!
      double precision ans(ans_size)
!
!      write(*,*) 'send job myid=', myid, ' zone=', zone
      call get_j_k(zone, j, k, nr)
!
!     common block zone_quantities
      ans(1) = tea(j,k)
      ans(2) = tna(j,k)
      ans(3) = n_e(j,k)
      ans(4) = B_field(j,k)
!     common block vol_em
      ans(5) = Eloss_sy(j,k)
      ans(6) = Eloss_cy(j,k)
      ans(7) = Eloss_br(j,k)
      ans(8) = Eloss_th(j,k)
      ans(9) = Eloss_tot(j,k)
!     common block deposition
      ans(10) = edep(j,k)
      ans(11) = prdep(j,k)
!     common block cens_energy
      ans(12) = ecens(j,k)
!     common block ecold
      ans(13) = ec_old(j,k)
!     common block fpair
      ans(14) = f_pair(j,k)
!     common block turbulence.
      ans(15) = q_turb(j,k)
      ans(16) = turb_lev(j,k)
!     common block nontherm.  These variables are modified by
!     FP_calc and update and returned to xec2d.
      ans(17) = amxwl(j,k)
      ans(18) = gmin(j,k)
      ans(19) = gmax(j,k)
      ans(20) = gbar_nth(j,k)
      ans(21) = N_nth(j,k)
      ans(22) = p_nth(j,k)
      ans(last) = r_acc(j,k)
!     common block Pnth
      do 10 i=1, num_nt
         ans(last+i) = Pnt(j, k, i)
         ans(num_nt+i+last) = gnt(i)
         ans(num_nt+num_nt+i+last) = f_nt(j,k,i)
 10   continue

!      write(*,*) 'myid=',myid,' n_field=',n_field(50,1,1)
!      stop
!     common block n_pair
      do 20 i=1, num_nt
         ans(num_nt*3+last+i) = n_pos(j,k,i)
 20   continue
!     common block dnpp
      do 30 i=1, num_nt
         ans(num_nt*4+last+i) = dn_pp(j,k,i)
 30   continue
!     common block annihil
      do 40 i=1, num_nt
         ans(num_nt*5+last+i) = dne_pa(j,k,i)
         ans(num_nt*6+last+i) = dnp_pa(j,k,i)
 40   continue
      do 50 i=1, num_nt
         ans(num_nt*7+last+i) = n_diff(j,k,i)
 50   continue
!     common block photonfield
      do 60 i=1, nphfield
         ans(num_nt*8+last+i) = n_field(i,j,k)
 60   continue

!      ans(num_nt*3+23) = dg
!     common block to_fp_calc
!      ans(num_nt*3+24) = E_tot_old
!      ans(num_nt*3+25) = E_tot_new
!      ans(num_nt*3+26) = hr_st_total
!      ans(num_nt*3+27) = hr_total
!      ans(num_nt*3+28) = f_t_implicit
!      ans(num_nt*3+29) = dr
!      ans(num_nt*3+30) = dz
!
      call MPI_SEND(ans, ans_size, MPI_DOUBLE_PRECISION, &
     &     available_proc, zone, MPI_COMM_WORLD, ierr)
!      write(*,*) 'myid=', myid, ' end of fp_send_job'
!
!
      return
      end
!
!
!     
!
!     This subroutine send the signal to a node that there are
!     no more calculations for it to perform.
!     J. Finke 18 March 2005
      subroutine FP_send_end_signal(node)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer node
!
      integer end_signal
      integer status(MPI_STATUS_SIZE)
!
!
      end_signal = jmax*kmax+1
!      write(*,*) 'myid=', myid, ' in send_end_signal, end_signal=', 
!     1     end_signal
      call MPI_SEND(0.d0, 1, MPI_DOUBLE_PRECISION, node, &
     &     end_signal, MPI_COMM_WORLD, ierr)
!
      return
      end
!
!
!
!
!
!     FP_recv_job recieves the zone number 
!     from the master node by the slave node.  It must 
!     be paired with FP_send_job or FP_send_end_signal 
!     in the master node.
!     J. Finke, 7 March 2005
      subroutine FP_recv_job(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer zone, get_zone, last
!
      integer ans_size
      parameter (last = 23)
      parameter (ans_size = num_nt*8+nphfield+last)
!
      integer j, k, i, end_signal
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
!
!      write(*,*) 'recv job myid=', myid, ' zone=', zone
      end_signal = jmax*kmax+1
!      write(*,*) 'myid=', myid, ' in fp_recv_job'
!
      call MPI_RECV(ans, ans_size, MPI_DOUBLE_PRECISION, master, &
     &     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      zone = status(MPI_TAG)
      if(zone.eq.end_signal) then
         return
      else
         call get_j_k(zone,j,k,nr)
      endif
!
!     common block zone_quantities
      tea(j,k) = ans(1)
      tna(j,k) = ans(2)
      n_e(j,k) = ans(3)
      B_field(j,k) = ans(4)
!      write(*,*) 'myid=', myid, ' in fp_recv_job, j=', j      
!      write(*,*) 'myid=', myid, ' in fp_recv_job, k=', k
!      write(*,*) 'myid=', myid, ' in fp_recv_job, zone=', zone
!      tea(j,k) = ans
!      write(*,*) 'myid=', myid, ' tea(j,k) =', tea(j,k)

!     common block vol_em
      Eloss_sy(j,k) = ans(5)
      Eloss_cy(j,k) = ans(6)
      Eloss_br(j,k) = ans(7)
      Eloss_th(j,k) = ans(8)
      Eloss_tot(j,k) = ans(9)
!     common block deposition
      edep(j,k) = ans(10)
      prdep(j,k) = ans(11)
!     common block cens_energy
      ecens(j,k) = ans(12)
!     common block ecold
      ec_old(j,k) = ans(13)
!     common block fpair
      f_pair(j,k) = ans(14)
!     common block turbulence.
      q_turb(j,k) = ans(15)
      turb_lev(j,k) = ans(16)
!     common block nontherm.  These variables are modified by
!     FP_calc and update and returned to xec2d.
      amxwl(j,k) = ans(17)
      gmin(j,k) = ans(18)
      gmax(j,k) = ans(19)
      gbar_nth(j,k) = ans(20)
      N_nth(j,k) = ans(21)
      p_nth(j,k) = ans(22)
      r_acc(j,k) = ans(last)
!     common block Pnth
      do 10 i=1, num_nt
         Pnt(j,k,i) = ans(last+i)
         gnt(i) = ans(num_nt+i+last)
         f_nt(j,k,i) = ans(num_nt+num_nt+i+last)
 10   continue
!     common block n_pair
      do 20 i=1, num_nt
         n_pos(j,k,i) = ans(num_nt*3+last+i) 
 20   continue
!     common block dnpp
      do 30 i=1, num_nt
         dn_pp(j,k,i) = ans(num_nt*4+last+i)
 30   continue
!     common block annihil
      do 40 i=1, num_nt
         dne_pa(j,k,i) = ans(num_nt*5+last+i) 
         dnp_pa(j,k,i) = ans(num_nt*6+last+i) 
 40   continue
!     common black diffusion
      do 50 i=1, num_nt
         n_diff(j,k,i) = ans(num_nt*7+last+i)
 50   continue
!     common block photonfield
      do 60 i=1, nphfield
         n_field(i,j,k) = ans(num_nt*8+last+i) 
 60   continue
!      write(*,*) 'recv j,k=',j,k,' dne_pa=',dne_pa(j,k,num_nt),
!     1     ' ans=', ans(num_nt*5+nphfield+22+num_nt) 
!      write(*,*) 'recv j,k=',j,k,' dnp_pa=',dnp_pa(j,k,num_nt),
!     1     ' ans=', ans(num_nt*6+nphfield+22+num_nt) 
!      write(*,*) 'myid=',myid,' j,k=',j,k,' dn_pp=',dn_pp(j,k,num_nt),
!     1     ' dne_pa=',dne_pa(j,k,num_nt),' dnp_pa=',dnp_pa(j,k,num_nt)
!      dg = ans(num_nt*3+23)
!     common block to_fp_calc
!      E_tot_old = ans(num_nt*3+24)
!      E_tot_new = ans(num_nt*3+25)
!      hr_st_total = ans(num_nt*3+26)
!      hr_total = ans(num_nt*3+27)
!      f_t_implicit = ans(num_nt*3+28)
!      dr = ans(num_nt*3+29)
!      dz = ans(num_nt*3+30)
!
!
!
!      write(*,*) 'myid=', myid, ' end of recv_job'
      return
      end
!
!
!
!     FP_send_result sends the results from FP_calc from the
!     slave node to the master node.  It must be paired with
!     FP_recv_result.
!     J. Finke, 7 March 2005
      subroutine FP_send_result(zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer zone
!
      integer ans_size
      parameter (ans_size = num_nt*4+30)
!
      integer i, j, k
!
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
!
!     nontherm variables, except for gbar_nth and N_nth, are modified
!     by FP_calc and send back to the master node by this subroutine.
!     tea is calculated in FP_calc and returned to the master node.
!     to_fp_calc has variables modified by FP_calc and returned
!     to the master node with this subroutine.
!      write(*,*) 'send result myid=', myid, ' zone=', zone
      call get_j_k(zone, j, k, nr)
!
!
!      common block nontherm
      ans(1) = amxwl(j,k)
      ans(2) = gmax(j,k)
      ans(3) = gmin(j,k)
      ans(4) = p_nth(j,k)
!     electron temp.
      ans(5) = tea(j,k)
!
!     common block Pnth
      do 10 i=1, num_nt
         ans(i+5) = Pnt(j, k, i)
         ans(i+5+num_nt) = gnt(i)
         ans(i+5+num_nt+num_nt) = f_nt(j,k,i)
         ans(i+5+num_nt*3) = n_pos(j,k,i) 
 10   continue
      ans(6+num_nt*4) = dg
!     common block to_fp_calc
!      ans(7+num_nt*3) = E_tot_old
!      ans(8+num_nt*3) = E_tot_new
!      ans(9+num_nt*3) = hr_total
!      ans(10+num_nt*3) = hr_st_total
      ans(11+num_nt*4) = f_t_implicit
      ans(12+num_nt*4) = dr
      ans(13+num_nt*4) = dz
!     common block d_update.
      ans(14+num_nt*4) = Te_new(j,k)
!     common block zone quantities and npos.
      ans(15+num_nt*4) = n_e(j,k)
      ans(16+num_nt*4) = B_field(j,k)
      ans(17+num_nt*4) = tea(j,k)
      ans(18+num_nt*4) = tna(j,k)
      ans(19+num_nt*4) = sum_E(j,k)
!      write(*,*) 'myid=', myid, ' ans=', ans
!      write(*,*) 'myid=', myid, ' Te_new(j,k)=', Te_new(j,k)
!
!      write(*,*) 'Pnt(1,1,50) =', Pnt(1,1,50)
!      write(*,*) 'Pnt(1,1,70) =', Pnt(1,1,70)
!      write(*,*) 'f_nt(1,1,50) =', f_nt(1,1,50)
!      write(*,*) 'f_nt(1,1,70) =', f_nt(1,1,70)
!      write(*,*) 'Te_new=', Te_new(j,k)
      call MPI_SEND(ans, ans_size, MPI_DOUBLE_PRECISION, &
     &     master, zone, MPI_COMM_WORLD, ierr)
!
!
      return
      end
!
!
!
!     FP_recv_result recieves the results of FP_calc from the
!     slave node.  It must be paired with FP_send_job.
!     J. Finke, 7 March 2005
      subroutine FP_recv_result(sender, zone)
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
      integer zone, get_zone, sender
!
      integer ans_size
      parameter (ans_size = num_nt*4+30)
!
      integer i, j, k
!
      integer status(MPI_STATUS_SIZE)
      double precision ans(ans_size)
!
!     nontherm variables, except for gbar_nth and N_nth, are modified
!     by FP_calc and recieved by the master node by this subroutine.
!     tea is calculated in FP_calc and recieved by the master node.
!     to_fp_calc has variables modified by FP_calc and recieved
!     by the master node with this subroutine.
!
!      write(*,*) 'recv result myid=', myid, ' zone=', zone
!    recieving common block nontherm
      call MPI_RECV(ans, ans_size, MPI_DOUBLE_PRECISION, &
     &             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
     &             status, ierr)
!     the zone numbers are stored in the tag.
      sender = status(MPI_SOURCE)
      zone = status(MPI_TAG)

      call get_j_k(zone,j,k,nr)

      amxwl(j,k) = ans(1)
      gmax(j,k) = ans(2)
      gmin(j,k) = ans(3)
      p_nth(j,k) = ans(4)
!     electron temp.
      tea(j,k) = ans(5)
!     common block Pnth
      do 10 i=1, num_nt
         Pnt(j,k,i) = ans(i+5)
         gnt(i) = ans(i+num_nt+5)
         f_nt(j,k,i) = ans(i+num_nt+num_nt+5)
         n_pos(j,k,i) = ans(i+3*num_nt +5) 
 10   continue
      dg = ans(4*num_nt+6)
!     common block to_fp_calc
!      E_tot_old = ans(3*num_nt+7)
!      E_tot_new = ans(3*num_nt+8)
!      hr_total = ans(3*num_nt+9)
!      write(*,*) 'myid=', myid, ' recv hr_total=', hr_total
!      hr_st_total = ans(3*num_nt+10)
      f_t_implicit = ans(4*num_nt+11)
      dr = ans(4*num_nt+12)
      dz = ans(4*num_nt+13)
!      write(*,*) 'myid=', myid, ' recv dz=', dz
!     common block d_update.
!      write(*,*) 'myid=', myid, ' recv ans=', ans
      Te_new(j,k) = ans(4*num_nt+14)
!     common block zone quantities
      n_e(j,k) = ans(4*num_nt+15)
      B_field(j,k) = ans(4*num_nt+16)
      tea(j,k) = ans(4*num_nt+17)
      tna(j,k) = ans(4*num_nt+18)
      sum_E(j,k) = ans(4*num_nt+19)
!      write(*,*) 'myid=', myid, ' recv Te_new=', Te_new(j,k)
!
!
!      write(*,*) 'Pnt(1,1,50) =', Pnt(1,1,50)
!      write(*,*) 'Pnt(1,1,70) =', Pnt(1,1,70)
!      write(*,*) 'f_nt(1,1,50) =', f_nt(1,1,50)
!      write(*,*) 'f_nt(1,1,70) =', f_nt(1,1,70)
!      write(*,*) 'Te_new=', Te_new(j,k)
      return
      end

      
!     This routine broadcast the results after update to every slave node.
      subroutine FP_end_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'

      integer zone
!
!
      integer i, j, k
!
      integer status(MPI_STATUS_SIZE)
!
!     nontherm variables, except for gbar_nth and N_nth, are modified
!     by FP_calc and send back to the master node by this subroutine.
!     tea is calculated in FP_calc and returned to the master node.
!     to_fp_calc has variables modified by FP_calc and returned
!     to the master node with this subroutine.
!c
!     common block nontherm
      write(*,*)'myid, in FP_end',myid
      call MPI_BCAST(amxwl, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmax, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gmin, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(P_nth, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!     electron temp.
      call MPI_BCAST(amxwl, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!     common block Pnth and n_pair
      call MPI_BCAST(Pnt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(gnt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(f_nt, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dg, 1, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(n_pos, jmax*kmax*num_nt, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!     common block to_fp_calc
!      call MPI_BCAST(f_t_implicit, 1, MPI_DOUBLE_PRECISION,
!     1               master, MPI_COMM_WORLD, ierr)
!      write(*,*)'myid, in the mid of FP_end',myid
      call MPI_BCAST(dr, 1, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(dz, 1, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!     common block d_update.
      call MPI_BCAST(Te_new, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!     common block zone quantities and npos
      call MPI_BCAST(n_e, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(B_field, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna, jmax*kmax, MPI_DOUBLE_PRECISION,&
     &               master, MPI_COMM_WORLD, ierr)
!
!
      return
      end

!
!
!
!
!
!
      subroutine coulomb_bcast
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer status(MPI_STATUS_SIZE)
!
!
!
!     broadcast number of zones to all processes
      call MPI_BCAST(nz, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nr, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
!
!     broadcast coulomb scattering libraries to all processes
      call MPI_BCAST(lib_dg_ce, num_nt+num_temp_max, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_dg_cp, num_nt+num_temp_max, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_disp_ce, num_nt+num_temp_max, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(lib_disp_cp, num_nt+num_temp_max, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(rate_gm1, num_nt, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea_min, 1, &
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tea_max, 1,&
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna_min, 1,&
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tna_max, 1,&
     &     MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
!
!
      return
      end
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:22:41 EDT 2006
! version: 2
! Name: J. Finke
! Commented out 'write' statements.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Wed Jul  5 12:10:18 EDT 2006
! version: 3
! Name: J. Finke
! Get rid of sending and recieving 'hr_total' and 
! 'hr_st_total' in 'FP_send_result' and 'FP_recv_result'. 
! This facilitates the 
! adding up of 'hr_total' and 'hr_st_total' in routine 
! 'E_add_up' in update2d.f      
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Sun Aug  6 18:55:16 EDT 2006
! version: 4
! Name: J. Finke
! In setup_bcast, it now broadcasts star_switch.   
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Wed Nov  1 15:18:29 EST 2006
! version: 5
! Name: J. Finke
! Changed setup_bcast so that it now broadcasts rand_switch. 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Wed Nov  1 16:49:51 EST 2006
! version: 6
! Name: J. Finke
! Changed setup_bcast so now it broadcasts rseed. 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! Mon Apr 06 14:13 CDT 2009
! Name: Xuhui Chen
! setup_bcast now broadcasts Pnt, f_nt, n_pos, dg, B_field
! and the injection parameters.
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
! Mon Apr 06 14:36 CDT 2009
! Name: Xuhui Chen
! Established the FP_end_bcast subroutine to broadcast the
! updated electron informations to every slave node.
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
