       subroutine imcvol2d
!      This subroutine creates new photons based on the radiative energy lose of the zone in the current time
!      time step, and begins to track them. The exact position, direction, energy and creation time of each 
!      photon is decided in a stochastic way based on its appropriate probability.
!      ! comment by Xuhui Chen 2013-07-15

       implicit none
       include 'mpif.h'
       include 'general.pa'
       include 'commonblock.f90'
!
!     MPI variables
      integer status(MPI_STATUS_SIZE), &
     &        num_sent, end_signal, sender,&
     &        num_zones       
      integer l, zone
!
!
!     MPI Initialization 
!      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      num_sent = 0
      end_signal = jmax*kmax+1
!
!
!     master part
      if(myid.eq.master) then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      call vol_bcast
      num_zones = nr*nz
      call vol_create_job_type
!
!
!     This sends the first round of zones to the slaves for processing.
      do 500 l= 1, min(num_zones, (numprocs-1))
         zone = l
         call vol_send_job(l, zone)
         num_sent = num_sent + 1
 500  continue  
!
!     As slaves complete processing a zone, this recieves the results
!     and sends the next zone to the slaves.
      do 501 l = 1, num_zones
         call vol_recv_result(sender, zone)
         if(num_sent.lt.num_zones) then
            zone = num_sent+1
            call vol_send_job(sender, zone)
            num_sent = num_sent + 1
         else
            zone = num_sent
            call vol_send_end_signal(sender,zone)
         endif
 501  continue
!
!
!
      else if(myid.ne.master) then
!     beginning of slave part
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      call vol_bcast
      num_zones = nr*nz
      call vol_create_job_type
!
!     if there are more nodes than work skip this node.
      if(myid.gt.num_zones) goto 990
!
 991  continue
      call vol_recv_job(zone)
      if(zone.eq.end_signal) goto 990
      call volume_em(zone)
!     eps_* are used immediately after being calculated, by the local slave node. No need to
!     end them to the master node any more. XC 2014.02.07
      call vol_calc(zone)
      call vol_send_result(zone)
      goto 991
 990  continue
!
      endif
!     end of slave part
!
!
!
!     end of imcvol2d
      return
      end           
!
!
!
!
      subroutine vol_calc(zone)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer zone
!    
      integer i, jv, kv
      integer npsurf, isurfmu
      integer npvol, npick, ivol, more
      integer jmid, jbot, jtop
      integer m, n, i_a
!
      double precision fas(jmax, kmax)
      double precision esurf, ew_save
      double precision ewcd, f_thermal
      double precision rnum0, rnum
      double precision psi, x1, x2
      double precision hubot, hutop, evol, norwlk, xxmore, rnpick
      double precision f_outer, f_inn, f_lower, f_upper, delz
      double precision fibran
      double precision cosa_eb, sina_eb, sina_eb2

!
!

!
      call get_j_k(zone, jv, kv, nr)
      call initialize_rand(jv,kv)
!
            ewcd = ewsv(jv,kv)
            ew_save = ewcd
            more = nsv(jv,kv)
                        
            f_thermal = Eloss_th(jv,kv)/Eloss_tot(jv,kv)
         
            if (jv.eq.1) then
               delz = z(1)
            else
               delz = z(jv) - z(jv-1)
            endif

            if (kv.eq.1) then
               f_inn  = (4.4d1/7.d0*rmin*delz)/zsurf(jv, kv)
               f_outer  = (4.4d1/7.d0*r(kv)*delz)/zsurf(jv, kv)
               f_upper = (2.2d1/7.d0 *(r(kv)**2&
     &                        - rmin**2))/zsurf(jv, kv)
               f_lower = (2.2d1/7.d0 *(r(kv)**2&
     &                        - rmin**2))/zsurf(jv, kv)
            else
               f_inn  = (4.4d1/7.d0*r(kv-1)*delz)/zsurf(jv, kv)
               f_outer  =  (4.4d1/7.d0*r(kv)*delz)/zsurf(jv, kv)
               f_upper = (2.2d1/7.d0 *(r(kv)**2&
     &                        - r(kv-1)**2))/zsurf(jv, kv)
               f_lower = (2.2d1/7.d0 *(r(kv)**2&
     &                        - r(kv-1)**2))/zsurf(jv, kv)
            endif
            
            f_outer = f_inn + f_outer
            f_upper = f_outer + f_upper
            f_lower = f_upper + f_lower
!
  50        format (' f_lower = ',e14.7)
            if (dabs(f_lower - 1.d0).gt.1.e-2) then
               write(4,50) f_lower
            endif
!                  
!      rseed = 8535           
            do 450 ivol = 1, more
               
               jph = jv
               kph = kv
               ew = ew_save
!               dcen = c_light*dt(1)*drand(0)
               dcen = c_light*dt(1)*fibran()
!               rnum = drand(0)
               rnum = fibran()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       thermal emission
!
        if (rnum.lt.f_thermal) then
                   i = 0
!                   rnum = drand(0)    
                   rnum = fibran()     
 120               i = i+1
                   if (eps_th(i, jph, kph).lt.rnum .and. i.lt.n_vol) &
     &                 go to 120
!
                   if (i.lt.n_vol) then
!                      xnu = E_ph(i)
!     1                    + drand(0)*(E_ph(i+1) - E_ph(i))
                      xnu = E_ph(i)&
     &                    + fibran()*(E_ph(i+1) - E_ph(i))
                   else
                      xnu = E_ph(i) 
                   endif
!
!                   rnum0 = drand(0)
                   rnum0 = fibran()
                   if (rnum0.lt.f_inn) then
!                       wmu = 2.d0*drand(0) - 1.d0
                       wmu = 2.d0*fibran() - 1.d0
                       if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                       if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
!
!                       x1 = drand(0)
!                       x2 = drand(0)
                       x1 = fibran()
                       x2 = fibran()
                       if (x1.lt.0.5) then
                          phi = 1.1d1/7.d0 + (1.1d1/7.d0)*x2
                          if (phi.lt.1.57079638d0) phi = 1.57079638d0
                       else 
                          phi = -1.1d1/7.d0 - (1.1d1/7.d0)*x2
                          if (phi.gt.-1.57079638d0) phi = -1.57079638d0
                       endif 

                       if (kph.eq.1) then
                          rpre = 1.00001*rmin
                       else
                          rpre = 1.00001*r(kph-1)
                       endif
                       if  (jph.eq.1) then          
                          zpre = z(1)*fibran()
                       else
                          zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                       endif

                   else if (rnum0.lt.f_outer) then
!
!                        wmu = 2.d0*drand(0) - 1.d0
!
                        wmu = 2.d0*fibran() - 1.d0
                        if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                        if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                        rpre = 0.999999*r(kph)
!
                        if (jph.eq.1) then          
                            zpre = z(1)*fibran()
                        else
                            zpre = z(jph-1)+fibran()*(z(jph)-z(jph-1))
                        endif
!                        phi = -1.1d1/7.d0 + 2.2d1/7.d0*drand(0)
                        phi = -1.1d1/7.d0 + 2.2d1/7.d0*fibran()
                        if (phi.lt.-1.5707963d0) phi = -1.57079063d0
                        if (phi.gt.1.5707963d0) phi = 1.5707963d0
!
!                        
                    else if (rnum0.lt.f_upper) then
!
                        wmu = fibran()
                        if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                        if (wmu.lt.0.d0) wmu = 0.d0
                        phi = 4.4d1/7.d0*fibran()
                        if (phi.gt.2.d0*pi) phi = 2.d0*pi
!
!                        wmu = drand(0)
!                        phi = 4.4d1/7.d0*drand(0)
!
!                        psi = drand(0)
                        psi = fibran()
                        Zpre= 0.999999*z(jph)    
                        if (kph.eq.1) then
                           rpre = dsqrt((rmin**2) + psi*(r(kph)**2&
     &                                            - rmin**2))
                        else
                           rpre = dsqrt((r(kph-1)**2) + psi*(r(kph)**2 &
     &                                               - r(kph-1)**2))
                        endif
!
                    else

!                        wmu = -drand(0)
!                        phi = 4.4d1/7.d0*drand(0)
!
                        wmu = -fibran()
                        phi = 4.4d1/7.d0*fibran()
                        if (wmu.gt.0.d0) wmu = 0.d0
                        if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                        if (phi.gt.2.d0*pi) phi = 2.d0*pi
!
                        if (jph.eq.1) then
                           Zpre = 1.000001*zmin
                           if (zpre.le.zmin) zpre = zmin + 1.d-6
                        else
                           Zpre = 1.000001*z(jph-1)
                           if (zpre.le.z(jph-1)) zpre = z(jph-1) + 1.d-6
                        endif
                          
!                        psi = drand(0)
                        psi = fibran()
                        if (kph.eq.1) then
                           rpre = dsqrt((rmin**2)+ psi*(r(kph)**2&
     &                                            - rmin**2))
                        else
                           rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2&
     &                                            - r(kph-1)**2))  
                        endif
                    endif             
              
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!         synchrotron emission
!
          else
!---------------------------------------
!         isotropic synchrotron emission
             if(abs(theta_local(jv,kv)).gt.pi/2)then
                    i = 0
!                    rnum = drand(0)
                    rnum = fibran()
  125               i  =  i + 1
                    if (eps_tot(i, jph, kph).lt.rnum .and. i.lt.n_vol) &
     &                 go to 125
!
                   if (i.lt.n_vol) then
!                      xnu = E_ph(i)
!     1                    + drand(0)*(E_ph(i+1) - E_ph(i))
                      xnu = E_ph(i)&
     &                    + fibran()*(E_ph(i+1) - E_ph(i))
                   else
                      xnu = E_ph(i) 
                   endif
!
!                    wmu = 2.d0*drand(0) - 1.d0
!                    phi = 4.4d1/7.d0*drand(0)
!
                    wmu = 2.d0*fibran() - 1.d0
                    phi = 2.d0*pi*fibran()
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
!
                    if (jph.eq.1) then          
                        zpre = z(1)*fibran()
                    else
                        zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                    endif
                    
!                    psi = drand(0)
                    psi = fibran()
                    if (kph.eq.1) then
                       rpre = dsqrt((rmin**2)+ psi*(r(kph)**2&
     &                                        - rmin**2))
                    else
                       rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2&
     &                                        - r(kph-1)**2))  
                    endif
!------------------------------------
!         anisotropic synchrotron emission
!         photons are assumed to travel at the same direction as the emitting electrons
             else
!     get the angle between the photon and the B
                    wmu = 2.d0*fibran() - 1.d0
                    phi = 2.d0*pi*fibran()
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                    cosa_eb = wmu*dcos(theta_local(jv,kv))+&
     &       dsqrt(1-wmu**2)*dsin(theta_local(jv,kv))*dcos(phi-pi/2.d0)
                    sina_eb2 = 1.d0-cosa_eb**2
                do while (fibran().gt.sina_eb2)
                    wmu = 2.d0*fibran() - 1.d0
                    phi = 2.d0*pi*fibran()
                    if (wmu.gt.9.9999999d-1) wmu = 9.9999999d-1
                    if (wmu.lt.-9.9999999d-1) wmu = -9.9999999d-1
                    cosa_eb = wmu*dcos(theta_local(jv,kv))+&
     &       dsqrt(1-wmu**2)*dsin(theta_local(jv,kv))*dcos(phi-pi/2.d0)
                    sina_eb2 = 1.d0-cosa_eb**2
                end do
                sina_eb = dsqrt(sina_eb2)

!      get the e-B angular bin of the photon
                 i_a=0
                do while(dsin(pi/num_a*i_a).le.sina_eb.and.i_a.lt.num_a)
                   i_a = i_a+1
                enddo
!      get the energy bin of the photon
                 i = 1
                 rnum = fibran()
                do while (eps_a(i_a,i,jph,kph).lt.rnum.and.i.lt.n_vol)
                   i = i+1
                enddo
!      get the energy of the photon
                if (i.lt.n_vol) then
                   xnu = E_ph(i)&
     &                 + fibran()*(E_ph(i+1) - E_ph(i))
                else
                   xnu = E_ph(i) 
                endif
!      get the position of the photon
                if (jph.eq.1) then          
                    zpre = z(1)*fibran()
                else
                    zpre = z(jph-1) + fibran()*(z(jph) - z(jph-1))
                endif
                 
                psi = fibran()
                if (kph.eq.1) then
                   rpre = dsqrt((rmin**2)+ psi*(r(kph)**2 - rmin**2))
                else
                  rpre = dsqrt((r(kph-1)**2)+ psi*(r(kph)**2&
     &                                    - r(kph-1)**2))  
                endif
             endif
          endif
      
!
!     Determine photon group of new photon energy
!     in photon energy bin structure for spectrum
!                   
                jbot = 1
                jtop = nphtotal + 1
                hubot = 1.000001*hu(1)
                hutop =.999999*hu(jtop)
                if( xnu .ge. hutop ) then
                   jgpsp = nphtotal
                   go to 300
                endif
                if( xnu .le. hubot ) then
                   jgpsp = 0
                   go to 300
                endif
!
  140           continue
                jmid = (jbot+jtop)/2
                if( jmid .eq. jbot ) go to 160
                if(xnu .eq. hu(jmid)) go to 160
                if(xnu .lt. hu(jmid)) then
                   jtop = jmid
                else
                   jbot = jmid
                endif
                go to 140
!
  160           continue
                jgpsp = jmid
!
  300           continue
!
!
!     Determine photon group of new photon energy
!   in photon energy bin structure for light curves
!
                jgplc = 0
                do 400 m = 1, nph_lc
                   if ((xnu.gt.Elcmin(m)).and.&
     &                 (xnu.le.Elcmax(m))) then
                      jgplc = m
                      goto 600
                   endif
  400           continue
!
  600           continue
!
!
!    Determine angular bin of photon
!
                do 700 n = 1, nmu
                   if (wmu.le.mu(n)) then
                      jgpmu = n
                      goto 750
                   endif
 700            continue
                jgpmu = nmu
!
 750            continue
                emitype = 0
!
!
                call imctrk2d(-1)
!
!
 450        continue
            if(kv.eq.1) then
!             write(*,*) 'end jv=',jv,' rand=',fibran()
             endif
!
!
 800  format(' rpre = ',e14.7,', zpre = ',e14.7)
 801  format(' xnu = ',e14.7)
 802  format(' phi = ',e14.7,',  wmu = ',e14.7)
 803  format(' jgpsp = ',i3,',  jgplc = ',i2)
 804  format(' jgpmu = ',i2)
 805  format(' ew = ',e14.7)
 806  format('     F = ',e14.7,', xnu = ',e14.7)
 807  format('     F = ',e14.7,', Delta (F) = ',e14.7)
 808  format('  dcen = ',e14.7)
!
!
!
      return
      end
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 13:30:44 EDT 2006
! version: 2
! Name: J. Finke
! Changed 'write' statments.      
