      subroutine reader(readerflag)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer readerflag
      integer j, k, l, n, m, t
      integer i
      integer npht_old, errors, lstart
      integer nin
      integer mhdtimestep

!     If the upper surface represents a star:
!     Rstar is the star's radius, and dist_star is the star's distance.
!     J. Finke, 7/20/06.
!
      double precision fibran, dummy, ran1
      real cputime
!
      character *30 fname
      character *81 buffer, trash
      character *3 nj, nk
      character *3 xmhdtimestep

!

!__________________________________________________________________
!
      data fname,mhdtimestep / 'input/input_01_01_001.dat',1 /
      nin=2
!
!       Input formats
!
 100  format(a80,i2)
 101  format(a80,i3)
 102  format(a80,i10) 
 110  format(a80,e14.7)
 112  format(a80,a30)
!
!       Output formats
!
  75  format('Number of vertical zones: ',i2)
  80  format('Number of radial zones:   ',i2)
  81  format(' Radius of the star:  ',e14.7,' cm')
  82  format(' Distance to the star:  ', e14.7,' cm')  
  85  format(' z_max = ',e14.7,' cm')
  90  format(' r_min = ',e14.7,' cm')
  95  format(' r_max = ',e14.7,' cm')
  96  format(' star_switch = ',i2)
  97  format(' Time to stop simulation: tstop = ',e14.7,' s')
  98  format(' Maximum time step in (dr/c): ',e14.7,' s')
 115  format('Number of time steps: ntime = ',i2)
 116  format('Number of time steps: ntime = ',i4)
 120  format('Boundary temp. at upper boundary, zone ',i2,': ',e14.7)
 125  format('Boundary temp. at lower boundary, zone ',i2,': ',e14.7)
 130  format('Boundary temp. at inner boundary, zone ',i2,': ',e14.7)
 135  format('Boundary temp. at outer boundary, zone ',i2,': ',e14.7)
 140  format('Zone (',i2,',',i2,'): File name: ',a30)
 145  format('Electron temperature: ',e14.7,' keV')
 150  format('Proton temperature:   ',e14.7,' keV')
 155  format('Number density:       ',e14.7,' cm^(-3)')
 160  format('Equipartition switch: ',i1)
 165  format('Magnetic field:       ',e14.7,' G')
 166  format('the angle between the helical magnetic field and z&
     & (positive means counter clockwise with increasing z)', e14.7)
!       right hand magnetic field, traditional polar coodinates
 170  format('Maxwellian fraction:  ',e14.7)
 172  format('Low-energy cut-off:   ',e14.7)
 175  format('High-energy cut-off:  ',e14.7)
 180  format('Power-law index:      ',e14.7)
 181  format('Turbulence sp. index  ',e14.7)
 182  format('Turbulence level:     ',e14.7)
 183  format('Acceleration time scale (Z/c): ',e14.7)
 184  format('Number of photon energy regions: ',i2)
 185  format('Energy Minimum of photon in ',i2,'. region:     ',e14.7)
 190  format('Energy Maximum of photon in ',i2,'. region:     ',e14.7) 
 191  format('Number of photon energy bins in ',i2,'. region: ',i3)
 195  format('Beginning of time step ',i2,': ',e14.7,' s')
 196  format('      End of time step ',i2,': ',e14.7,' s')
 197  format('Number of mu-bins: ',i2)
 198  format('Energy minimum of ',i2,'. light curve energy bin  : ',&
     &       e14.7)
 199  format('Energy maximum of ',i2,'. light curve energy bin  : ',&
     &       e14.7)
 201  format('Number of light curve energy bins:  ',i2)
 202  format('File name for time-integrated energy spectra: ',a30)
 203  format('File name for time-integrated photon spectra: ',a30)
 204  format('File name for light curve in ',i2,'. angular bin: ',a30)
 206  format('Number of particles per photon cycle: ',i10)
 207  format('No Compton reflection')
 208  format('Compton reflection at lower boundary only')
 209  format('Compton reflection at outer disk only')
 211  format('Compton reflection at lower boundary and outer disk')
 212  format('Pair processes on')
 213  format('Pair processes off')
 214  format('Number of initial calls to random no. generator: ',i6)
 216  format('Census file name 1: ',a30)
 217  format('Census file name 2: ',a30)
 218  format('Disk heating off')
 219  format('Disk heating on')
 221  format('File name for event file: ',a30)
 222  format('Constant temperature (electron spectrum) selected.')
 223  format('Electron dynamics calculation option selected, fp_sw:', i2)
 224  format('No coronal flare.')
 226  format('Coronal flare parameters:')
 227  format('Center of coronal heating flare, radius [cm]: ', e14.7)
 228  format('Center of coronal heating flare, height [cm]: ', e14.7)
 229  format('Center of coronal heating flare, time [s]: ', e14.7)
 231  format('Width of coronal heating flare, radius [cm]: ', e14.7)
 232  format('Width of coronal heating flare, height [cm]: ', e14.7)
 233  format('Width of coronal heating flare, time [s]: ', e14.7)
 234  format('Flare amplitude: (B_max/B_0)^2; (delta T_p)/T_p: ', e14.7)
 240  format('Upper boundary input spectrum from file: ', a30)
 241  format('Lower boundary input spectrum from file: ', a30)
 242  format('Inner boundary input spectrum from file: ', a30)
 243  format('Outer boundary input spectrum from file: ', a30)
 244  format('File name for temperature history: ', a30)
!
!           Error messages
!
 205  format('nz reset to ',i3)
 210  format('nr reset to ',i3)
 215  format('nphbins(',i2,') reset to ',i3)
 220  format('Ephmin(',i2,') reset to ',e14.7,' keV')
 225  format('Ephmax must be larger than Ephmin in photon region '&
     &       ,i2,'.')
 230  format('nmu reset to ',i2)
 235  format('Ephmax must be larger than Ephmin in lc energy bin '&
     &       ,i2,'.')
 245  format('flare switch no(0), injection(1), external(2), B(3):',i2)
 246  format('particle injection Energy lower limit',e14.7)
 247  format('particle injection Energy upper limit',e14.7)
 248  format('particle injection power-law index',e14.7)
 249  format('time to start the particle injection',e14.7)
 250  format('injection rate (erg/s)',e14.7)
 251  format(' Particle escape no(0)/esc(1)/diff(2)/reflect(z:3/r:4/all:5) esc_sw',i2)
 252  format('injection Gaussian partcle energy:',e14.7)
 253  format('injection Gaussian sigma:',e14.7)
 254  format('lorentz factor of the bulk motion:',e14.7)
 255  format('particle escape time scale(z/c):',e14.7)
 257  format('injection gamma max variation switch on/off (1/0):',i2)
 258  format('medium energy Gaussian electron pick up on/off (1/0):',i2)
 259  format('electron pick up rate (1/s/cm^3)',e14.7)
 260  format('Radius of the broad line region (cm):',e14.7)
 261  format('Fraction of disk luminosity reprocessed by BLR:',e14.7)
 262  format('Radius of the infrared torus (cm):',e14.7)
 263  format('Fraction of disk luminosity reprocessed by torus:',e14.7)
 264  format('Radius of the accretion disk (cm):',e14.7)
 265  format('Distance from the jet to the central engine (cm):',e14.7)
 266  format('3 split numbers used, and the 3rd split trigger:',&
     &       e14.7,2i10,e14.7)
 267  format('the angle between the helical magnetic field and z &
     & (positive means counter clockwise with increasing z)', e14.7)

!
      if(readerflag/=0) goto 789
!      open(unit=nin, file='input.dat', status='unknown')
      open(unit=nin, file='input/input.dat', status='unknown')
      open(unit=3, file='errors.txt', status='unknown')
!
      errors = 0
!
      read(nin, 101) buffer, nz
      read(nin, 101) buffer, nr
      read(nin, 110) buffer, z(nz)
      read(nin, 110) buffer, rmin
      read(nin, 110) buffer, r(nr)
      read(nin, 101) buffer, star_switch
      if(star_switch.eq.1) then
         read(nin, 110) buffer, Rstar
         read(nin, 110) buffer, dist_star
      endif
      read(nin, 110) buffer, tstop
      read(nin, 110) buffer, mcdt
      read(nin, 101) buffer, ntime
      if (nz.gt.jmax) then
         nz = jmax
         write(*,*) 'Maximum number of vertical zones exceeded!'
         write(*,205) nz
         write(3,*) 'Maximum number of vertical zones exceeded!'
         write(3,205) nz 
         errors = errors + 1
      endif
      if (nr.gt.kmax) then
         nr = kmax
         write(*,*) 'Maximum number of radial zones exceeded!'
         write(*,210) nr
         write(3,*) 'Maximum number of radial zones exceeded!'
         write(3,210) nr
         errors = errors + 1
      endif
      if (ntime.gt.ntmax) then
         write(*,116) ntime
         write(*,*) 'Fatal error:'
         write(*,*) 'Maximum number of time steps exceeded!'
         write(*,*) 'Re-define input!'
         write(*,*)
         write(3,116) ntime
         write(3,*) 'Fatal error:'
         write(3,*) 'Maximum number of time steps exceeded!'
         write(3,*) 'Re-define input!'
         errors = errors + 1
         close(nin)
         close(3)
         close(4)
         stop
      endif
      write(4,*)
      write(4,75) nz
      write(4,80) nr
      write(4,85) z(nz)
      write(4,90) rmin
      write(4,95) r(nr)
      write(4,96) star_switch
      if(star_switch.eq.1) then
         write(4,81) Rstar
         write(4,82) dist_star
      endif
      write(4,97) tstop
      write(4,98) mcdt
      write(4,115) ntime
      
      
!
!   ... Read ntime;
!     iterate from here ntime times
!
      do 350 t = 1, ntime
         read(nin, 110) buffer, t0(t)
         read(nin, 110) buffer, t1(t)
         write(4,*)
         write(4, 195) t, t0(t)
         write(4, 196) t, t1(t)
!
         read(nin, 110) buffer, tbbu(1,1)
         read(nin, 112) buffer, u_fname(1, 1)
         read(nin, 110) buffer, tbbl(1,1)
         read(nin, 112) buffer, l_fname(1, 1)
         do 300 k = 2, nr
            tbbu(k,t)=tbbu(1,1)
            u_fname(k,t)=u_fname(1,1)
            tbbl(k,t)=tbbl(1,1)
            l_fname(k,t)=l_fname(1,1)
 300     continue
!  
         do 310 j = 1, nz
             tbbi(j,t) = 0.d0
             tbbo(j,t) = 0.d0 ! Xuhui 11/19/08
!            read(nin, 110) buffer, tbbi(j,t)
!            if (tbbi(j,t).lt.0.d0) read(nin, 112) buffer, i_fname(j, t)
!            read(nin, 110) buffer, tbbo(j,t)
!            if (tbbo(j,t).lt.0.d0) read(nin, 112) buffer, o_fname(j, t)
 310     continue
!  
            if (tbbu(nr,t).lt.0.d0) then
               write(4, 240) u_fname(nr,t)
            else
               write(4, 120) nr, tbbu(nr,t)
            endif
            if (tbbl(nr,t).lt.0.d0) then
               write(4, 241) l_fname(nr,t)
            else
               write(4, 125) nr, tbbl(nr,t)
            endif
!  
            if (tbbi(nz,t).lt.0.d0) then
               write(4, 242) i_fname(nz,t)
            else
               write(4, 130) nz, tbbi(nz,t)
            endif
            if (tbbo(nz,t).lt.0.d0) then
               write(4, 243) o_fname(nz,t)
            else
               write(4, 135) nz, tbbo(nz,t)
            endif
 350  continue
!
!      End iteration (ntime)
!
!      Read energy grid ...
!
      write(4,*) 'before spec_switch'
      read(nin,100) buffer, spec_switch
      read(nin,100) buffer, nphreg
      write(4,*)
      write(4,184) nphreg
      if (nphreg.gt.nregmax) then
         write(*,*) 'Fatal error:'
         write(*,*) 'Maximum number of photon regions exceeded!'
         write(*,*) 'Re-define input!'
         write(*,*)
         write(3,*) 'Fatal error:'
         write(3,*) 'Maximum number of photon regions exceeded!'
         write(3,*) 'Re-define input!'
         errors = errors + 1
         close(nin)
         close(3)
         close(4)
         stop
      endif
!
      nphtotal = 0
      Npht_old = 0
      do 370 m = 1, nphreg
         read(nin,110) buffer, Ephmin(m)
         read(nin,110) buffer, Ephmax(m)
         read(nin,101) buffer, nphbins(m)
         npht_old = nphtotal
         nphtotal = nphtotal + nphbins(m)
         if (nphtotal.gt.nphomax) then
            nphbins(m) = nphomax - npht_old
            nphtotal = nphomax
            write(*,*) 'Maximum number of photon bins exceeded!'
            write(*,215) m, nphbins(m)
            write(3,*) 'Maximum number of photon bins exceeded!'
            write(3,215) m, nphbins(m)
            errors = errors + 1
         endif
         if (m.gt.1) then
           if (dabs((Ephmin(m)-Ephmax(m-1))/Ephmin(m)).gt.1.d-6) then
              Ephmin(m) = Ephmax(m-1)
              write(*,*) 'Energy regions must be adjacent.'
              write(*,220) m,Ephmin(m)
              write(3,*) 'Energy regions must be adjacent.'
              write(3,220) m,Ephmin(m)
              errors = errors + 1
           endif
         endif
         if (Ephmin(m).ge.Ephmax(m)) then
            write(*,*)
            write(*,*) 'Fatal error:'
            write(*,225) m
            write(*,*) 'Re-define input!'
            write(3,*)
            write(3,*) 'Fatal error:'
            write(3,225) m
            write(3,*) 'Re-define input!'
            errors = errors + 1
            close(nin)
            close(3)
            close(4)
            stop
         endif
 370  continue
         
      do 380 m = 1, nphreg
         write(4,185) m, Ephmin(m)
         write(4,190) m, Ephmax(m)
         write(4,191) m, nphbins(m)
 380   continue
!
      write(4,*)
      read(nin,101) buffer, nmu
      if (nmu.gt.nmumax) then
         nmu = nmumax
         write(*,*)
         write(*,*) 'Maximum number of angular bins exceeded!'
         write(*,230) nmu
         write(3,*) 'Maximum number of angular bins exceeded!'
         write(3,230) nmu
         errors = errors + 1
      endif
      write(4,197) nmu
!      
!  ENERGY BINNING FOR LIGHT CURVES
!
      read(nin,101) buffer, nph_lc
      if (nph_lc.gt.nphlcmax) then
         write(*,*) 'Fatal error:'
         write(*,*) 'Maximum number of lc photon bins exceeded!'
         write(*,*) 'Re-define input!'
         write(*,*)
         write(3,*) 'Fatal error:'
         write(3,*) 'Maximum number of lc photon bins exceeded!'
         write(3,*) 'Re-define input!'
         errors = errors + 1
         close(nin)
         close(3)
         close(4)
         stop
      endif
      write(4,*)
      write(4,201) nph_lc
!                    
         do 385 m = 1, nph_lc
!        
          read(nin,110) buffer, Elcmin(m)
          read(nin,110) buffer, Elcmax(m)
          if (Elcmin(m).ge.Elcmax(m)) then
            write(*,*)
            write(*,*) 'Fatal error:'
            write(*,235) m
            write(*,*) 'Re-define input!'
            write(3,*)
            write(3,*) 'Fatal error:'
            write(3,235) m
            write(3,*) 'Re-define input!'
            errors = errors + 1
            close(nin)
            close(3)
            close(4)
            stop
          endif
  385   continue
!         
! 
         do 386 m = 1, nph_lc
!                   
         write(4, 198) m, Elcmin(m)
         write(4, 199) m, Elcmax(m)
  386   continue       
!
!         Read output file names for spectra
!               and light curves
!
         read(nin,112) buffer, spname
         read(nin,112) buffer, phname
         read(nin,112) buffer, lcname(1)
         read(nin,112) buffer, eventfile
         read(nin,112) buffer, temp_file
         write(4,202) spname
         write(4,203) phname
         write(4,204) 1, lcname(1)
!
         do 390 l = 1, 30
            if (lcname(1)(l:l).eq.'_') goto 395
 390     continue
 395     continue
         if (l.gt.25) then
            write(*,*)
            write(*,*) 'Fatal error:'
            write(*,*) 'File name for light curve too long!'
            write(*,*) 'Re-define input!'
            write(3,*)
            write(3,*) 'Fatal error:'
            write(3,*) 'File name for light curve too long!'
            write(3,*) 'Re-define input!'
            errors = errors + 1
            close(nin)
            close(3)
            close(4)
            stop
         endif
!         
         do 397 n = 2, nmu
            lcname(n) = lcname(1)
            lcname(n)(l+1:l+1) = char(48 + int(n/10))
            lcname(n)(l+2:l+2) = char(48 + n - 10*(int(n/10)))
            lcname(n)(l+3:l+6) = '.dat'
  397    continue
!
         do 398 n = 2, nmu
  398    write(4, 204) n, lcname(n)
         write(4, 221) eventfile
         write(4, 244) temp_file
!
      read(nin,102) buffer, nst
      write(4, 206) nst
      read(nin,102) buffer, rseed
      call cpu_time(cputime)
      write(4, 214) rseed
      read(nin,102) buffer, rand_switch
!
!
      if(rand_switch.eq.2) dummy = ran1(-rseed)
!
      read(nin,100) buffer, cr_sent
      if (cr_sent.eq.0) then
         write(4,207)
      else if (cr_sent.eq.1) then
         write(4,208)
      else if (cr_sent.eq.2) then
         write(4,209)
      else
         write(4,211)
      endif
!

      read(nin,100) buffer, turb_sw
      write(4,'("turbulent switch:",i2)')turb_sw

      read(nin,110) buffer, acc_prob
      write(4,'("acceleration probability:",e14.7)')acc_prob
!
      read(nin, 100) buffer, pair_switch
      if (pair_switch.eq.0) then
         write(4, 213)
      else
         write(4, 212)
      endif
!
      read(nin, 100) buffer, fp_sw
      if (fp_sw.eq.0) then
         write(4, 222)
      else
         write(4, 223) fp_sw
      endif
!
      read(nin, 100) buffer, cf_sentinel
!      if (cf_sentinel.eq.1) then
!         write(4, *)
!         write(4, 226)
         read(nin,110) buffer, r_flare
         read(nin,110) buffer, z_flare
         read(nin,110) buffer, t_flare
         read(nin,110) buffer, sigma_r
         read(nin,110) buffer, sigma_z
         read(nin,110) buffer, sigma_t
         read(nin,110) buffer, flare_amp
      if (cf_sentinel.eq.1) then
         write(4, *)
         write(4, 226)
         write(4, 227) r_flare
         write(4, 228) z_flare
         write(4, 229) t_flare
         write(4, 233) sigma_t
      else
         write(4, 224)
         r_flare = 0.d0
         z_flare = 0.d0
         t_flare = 0.d0
!         sigma_r = 1.d0
!         sigma_z = 1.d0
         sigma_t = 1.d0
!         flare_amp = 0.d0
      endif
         write(4, 231) sigma_r
         write(4, 232) sigma_z
         write(4, 234) flare_amp
!
      read(nin,101)buffer,inj_switch
      read(nin,100)buffer,g2var_switch
      read(nin,110)buffer,inj_g1
      read(nin,110)buffer,inj_g2
      read(nin,110)buffer,inj_p
      read(nin,110)buffer,inj_t
      read(nin,110)buffer,inj_L
      read(nin,100)buffer,esc_sw
      read(nin,110)buffer,r_esc
      read(nin,100)buffer,pick_sw
      read(nin,110)buffer,dummy
      read(nin,110)buffer,inj_gg
      read(nin,110)buffer,inj_sigma
      read(nin,110)buffer,g_bulk
      inj_v = dsqrt(1.d0-1.d0/g_bulk**2)*c_light
      if(fp_sw.ne.0)then
        write(4, *)
        write(4, 245)inj_switch
        write(4, 257)g2var_switch
        write(4, 246)inj_g1
        write(4, 247)inj_g2
        write(4, 248)inj_p
        write(4, 249)inj_t
        write(4, 250)inj_L
        write(4, 251)esc_sw
        write(4, 255)r_esc
        write(4, 258)pick_sw
        if(pick_sw.ge.1)then
          write(4, 252)inj_gg
          write(4, 253)inj_sigma
        endif
!        write(4, 259)dummy
        write(4, 254)g_bulk
      endif
      read(nin,110)buffer, R_blr
      read(nin,110)buffer, fr_blr
      read(nin,110)buffer, R_ir
      read(nin,110)buffer, fr_ir
      read(nin,110)buffer, R_disk
      read(nin,110)buffer, d_jet
      read(nin,110)buffer, split1
      read(nin,102)buffer, split2
      read(nin,102)buffer, split3
      read(nin,110)buffer, spl3_trg
        write(4, 260) R_blr
        write(4, 261) fr_blr
        write(4, 262) R_ir
        write(4, 263) fr_ir
        write(4, 264) R_disk
        write(4, 265) d_jet
        write(4, 266) split1,split2,split3,spl3_trg

      if (errors.eq.0) then
         write(3,*) 'No input errors detected.'
      endif
!
      close(nin)
      close(3)


 789  continue
!
!        Read zone quantities:
!
      write(xmhdtimestep,"(i3.3)") mhdtimestep
      fname(19:21)=trim(xmhdtimestep)
      do 500 j = 1, nz
!       
         nj(1:1) = char(48 + int(j/10))
         nj(2:2) = char(48 + j - 10*(int(j/10)))
!
         do 490 k = 1, nr
!
            nk(1:1) = char(48 + int(k/10))
            nk(2:2) = char(48 + k -10*(int(k/10)))
            fname(13:14) = nj(1:2)
            fname(16:17) = nk(1:2)
!
            write(4,*)
            write(4,140) j, k, fname      
!
            open(unit=nin, file=fname, status='unknown')
!
            gmin(j,k) = 1.d0
            read(nin,'(e14.7)') currenttimestop
            read(nin,'(e14.7)') tea(j,k)
            read(nin,'(e14.7)') tna(j,k)
            read(nin,'(e14.7)') n_e(j,k)
            read(nin,'(i2)') ep_switch(j,k)
            read(nin,'(e14.7)') B_input(j,k)
            read(nin,'(e14.7)') theta_b(j,k)
            read(nin,'(e14.7)') amxwl(j,k)
            read(nin,'(e14.7)') gmin(j,k)
            read(nin,'(e14.7)') gmax(j,k)
            read(nin,'(e14.7)') p_nth(j,k)
            read(nin,'(e14.7)') q_turb(j,k)
            read(nin,'(e14.7)') turb_lev(j,k)
            read(nin,'(e14.7)') r_acc(j,k)
!
            write(4,145) tea(j,k)
            write(4,150) tna(j,k)
            write(4,155) n_e(j,k)
            write(4,160) ep_switch(j,k)
            if (ep_switch(j,k).eq.0)then
               write(4,165) B_input(j,k)
               write(4,166) theta_b(j,k)
            endif
!     convert angle from degree to radian
            theta_b(j,k) = theta_b(j,k)/1.8d2*pi
            write(4,170) amxwl(j,k)
            if (amxwl(j,k).lt.9.999999d-1) then
               write(4,172) gmin(j,k)
               write(4,175) gmax(j,k)
               write(4,180) p_nth(j,k)
            endif
            write(4,181) q_turb(j,k) 
            write(4,182) turb_lev(j,k)
            write(4,183) r_acc(j,k)
!
            close(nin)
!
 490     continue
 500  continue
      mhdtimestep=mhdtimestep+1
      write(4,*) "Zone data complete, next mhdtimestep=",mhdtimestep

! exit
      return
      end                
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:32:54 EDT 2006
! version: 2
! Name: J. Finke
! Changed common block 'random'.     
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Sun Aug  6 18:14:36 EDT 2006
! version: 3
! Name: J. Finke
! Added spec_switch, so quick look spectra can be 
! spectra incident on top and bottom boundaries. This 
! is used in photon bubble lowdensity simulations. Also 
! added possibility of using upper surface to represent 
! a star. Added reading in of parameters star_switch, 
! dist_star, and Rstar.      
