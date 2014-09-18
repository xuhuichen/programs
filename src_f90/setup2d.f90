!     everything is done in master node. Then some of them are broadcasted to
!     slave nodes 
      include 'nrutil.f90'
      include 'nrtype.f90'
      include 'nrface.f90'

      subroutine setup(setupflag)
      USE nrface, ONLY : bessk
      implicit none
      include 'mpif.h'
      include 'general.pa'
      include 'commonblock.f90'
!
      integer setupflag
      integer i, j, k, n, m, m1, nunit_temp
      integer status(MPI_STATUS_SIZE)
!
      double precision delj, delk
      double precision dmu
      double precision dE, priE, secE
      double precision dist_min
!____________________________________________________________________
      double precision sum_iniE
      double precision inj_sum, inj_E, inject_ne(num_nt), inj_y,&
     &                 gamma0(num_nt), vol_tot, beta !local varible
!_____________________________________________________________________

!
!___________________________________________________________
!
!
!     most of setup is done by the master node only.
      if(myid.eq.master) then
       if(setupflag==1) goto 2001
!
!
!
!     Call imcdate to load Compton cross-section averages
!       and average energy exchange data into memory
!
!
      call imcdate
!
!
!      Calculate Compton reflection matrix
!
 2001 continue
      call Pref_calc
      call Wref_calc
      if(setupflag==1) goto 2002
!

!
!        Initialize time and time step
!
      ncycle = 0
      time = 0.d0
!      dist_min = dmax1((2.d0*(r(nr) - rmin)), z(nz))
      dist_min = dmin1((r(nr)-rmin)/nr,z(nz)/nz) ! Xuhui 5/11/09
      dt(1) = mcdt*dist_min/inj_v
      dt(2) = dt(1)
      write(4,*)
      write(4, 100) dt(1)
 100  format('Initial time step: ',e14.7,' s')
!
!
!      Calculate spatial zoning (linear)
!
      j = 1
      k = 1
      delj = z(nz)/dble(nz)
      delk= (r(nr) - rmin)/dble(nr)
      zmin = 0.d0
      z(1) = delj
      r(1) = rmin + delk
!
      do 120 j = 1, nz-1
 120  z(j+1) = z(j) + delj
!
      do 130 k = 1, nr-1
 130  r(k+1) = r(k) + delk
!
      write(4,*)
      do 133 j = 1, nz
 133  write(4,810) j,z(j)
      do 136 k = 1, nr
 136  write(4,815) k,r(k)
!
      vol(1,1) = pi*(r(1)**2 - rmin**2)*z(1)
      zsurf(1,1) = 2.d0*pi*((r(1) + rmin)*z(1) &
     &                     + (r(1)**2 - rmin**2))
      do 138 k = 2, nr
         vol(1,k) = pi*(r(k)**2 - r(k-1)**2)*z(1)
         zsurf(1,k) = 2.d0*pi*((r(k) + r(k-1))*z(1) &
     &                        + (r(k)**2 - r(k-1)**2))
 138  continue
      do 143 j = 2, nz
         vol(j, 1) = pi*(r(1)**2 - rmin**2)&
     &                        *(z(j) - z(j-1))
         zsurf(j, 1) = 2.d0*pi*((r(1) + rmin)*(z(j) - z(j-1))&
     &                         + (r(1)**2 - rmin**2))
         do 142 k = 2, nr
            vol(j,k) = pi*(r(k)**2 - r(k-1)**2)&
     &                          *(z(j) - z(j-1))
            zsurf(j, k) = 2.d0*pi*((r(k) + r(k-1))&
     &                             *(z(j) - z(j-1))&
     &                            + (r(k)**2 - r(k-1)**2))
 142     continue
 143  continue
!
      Asurfu(1) = pi*(r(1)**2 - rmin**2)
      Asurfl(1) = Asurfu(1)
      Asurfi(1) = 2.d0*pi*rmin*z(1)
      Asurfo(1) = 2.d0*pi*r(nr)*z(1)
!
      do 1100 j = 2, nz
         Asurfi(j) = 2.d0*pi*rmin*(z(j) - z(j-1))
         Asurfo(j) = 2.d0*pi*r(nr)*(z(j) - z(j-1))
 1100 continue
      do 1110 k = 2, nr
         Asurfu(k) = pi*(r(k)**2 - r(k-1)**2)
         Asurfl(k) = Asurfu(k)
 1110 continue
!
!      write(*,*) 'Partial surfaces calculated.'
!
!
!       Calculate parameters of nonthermal distribution:
!
 2002 continue
      sum_iniE = 0.d0
      do 1150 j = 1, nz
         do 1140 k = 1, nr
            f_pair(j,k) = 0.d0
            do 1130 i = 1, num_nt
               dne_pa(j, k, i) = 0.d0
               dnp_pa(j, k, i) = 0.d0
 1130       continue
!            write(*,*) 'j = ',j,'; k = ',k,'; calling gam_min'
            call gam_min(j, k, tea(j,k))
!            write(*,*) 'j = ',j,'; k = ',k,'; calling P_nontherm'
            call P_nontherm(j, k)
!            write(*,*) 'Returned from P_nontherm.'
             do i = 1, num_nt-1
                sum_iniE = sum_iniE + &
     &        f_nt(j,k,i)*(gnt(i+1)-gnt(i))*(gnt(i)+1)*n_e(j,k)*vol(j,k)
             enddo
 1140    continue
 1150 continue
      write(*,*) 'Initial Energy amount:', sum_iniE*8.186d-7
      if(setupflag==1) goto 2003
!
!      write(*,*) 'Nonthermal parameters calculated.'

!     relativistic Mawellian distribution for pick up
        nunit_temp=15
        open(unit=nunit_temp,file='electrons/pick_dis.dat',status='unknown')
        do i = 1, num_nt
          beta = dsqrt(1-1/(gnt(i)+1)**2)
          pick_dis(i) = beta/(1-beta**2)/inj_gg/bessk(2,real(1/inj_gg))*dexp(-(gnt(i)+1)/inj_gg)
          write(nunit_temp,'(e14.7)')dmax1(pick_dis(i),1.d-30)
        enddo
        close(nunit_temp)
!
!
!       Calculating the set-up angle for the problem 
!                     (linear)
!
        dmu = 2.d0/dble(nmu)
        mu(1) = dmu - 1.d0
!        
        do 145 n = 2, nmu
          mu(n) = mu(n-1)+dmu
 145    continue
!
        write(4,*)
        do 150 n = 1, nmu
 150    write(4,805) n, mu(n)
!
!
!  Calculating the Energy set-up grid for the problem 
!              (logarithmic)
!        
        i = 1
        do 160 m = 1, nphreg
           priE  = log(Ephmax(m)/Ephmin(m))
           secE = priE/dble(nphbins(m))
           dE    = exp(secE)
           hu(i) = Ephmin(m)
           do 170  m1=1, nphbins(m)
              i = i+1
              hu(i) = hu(i-1)*dE
  170      continue
  160   continue
!
!         write(4,*)
!         do 180 i = 1, nphtotal+1
!  180    write(4,800) i, hu(i)
!
!
!          Array initializations
!
         do 210 j = 1, nz
            do 200 k = 1, nr
               ecens(j,k) = 0.d0
               npcen(j,k) = 0
  200       continue
  210    continue
!
         do 220 k = 1, nr
            Ed_in(k) = 0.d0
            Ed_abs(k) = 0.d0
            Ed_ref(k) = 0.d0
  220    continue
!
!
!      Initialize arrays for internal photon fields
!        for gamma-gamma absorption calculations
!
         dE = dexp(log(1.d2)/dble(n_gg))
         E_gg(1) = 5.d1
         do 300 i = 1, n_gg
            if (i.gt.1) E_gg(i) = E_gg(i-1)*dE
            do 280 j = 1, nz
               do 260 k = 1, nr
                  n_ph(i, j, k) = 0.d0
                  k_gg(i, j, k) = 0.d0
 260           continue
 280        continue     
 300     continue
!
!
!
!      Initialize arrays for internal photon
!       fields for Compton loss calculations
!
      dE = dexp(log(1.d20)/dble(nphfield))
      E_field(1) = 1.d-10
      do i = 2, nphfield
        E_field(i) = E_field(i-1)*dE
      enddo
!
!
!     Calculate table of IC loss kernel values
!
      call IC_loss
!
!
!     Initialize spectral and light curve output arrays
!
         do 400 n = 1, nmu
            do 380 i = 1, nphtotal
 380          fout(n, i) = 0.d0
 400     continue
!
!
!     Initialize sums for temperature distribution output
!     pick up rate is balancing the escape to keep the electron density
!
 2003 continue

         r_acc_peak=1.d40   ! r_acc_peak is only used by the master in update
         time_sum = 0.d0
         do 450 j = 1, nz
            do 440 k = 1, nr
               T_sum(j,k) = 0.d0
               pick_rate(j,k)=n_e(j,k)/(r_esc*z(nz)/c_light)
               if(r_acc_peak.gt.r_acc(j,k))r_acc_peak=r_acc(j,k)
 440        continue
 450     continue

!   
!
!     end master-only part; master and slaves call setup_bcast.
      endif

      call setup_bcast
!
!   
!
!
!          Output formats for log file
!
  800    format('hu(',i3,') = ',e14.7,' keV')
  805    format('mu(',i3,') = ',e14.7)
  810    format('Vertical zone boundary z(',i2,') = ',e14.7,' cm')
  815    format('Radial zone boundary r(',i2,') = ',e14.7,' cm')
!
      return
      end
!
!
! ==========================================================
        FUNCTION bessi0_s(x)
        USE nrtype; USE nrutil, ONLY : poly
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: bessi0_s
        ! Returns the modified Bessel function I0(x) for any real x.
        REAL(SP) :: ax
        REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
        3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
        0.45813e-2_dp/) ! Accumulate polynomials in double precision.
        REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
        0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
        -0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
        0.392377e-2_dp/)
        ax=abs(x)
        if (ax < 3.75) then !Polynomial fit.
        bessi0_s=poly(real((x/3.75_sp)**2,dp),p)
        else
        bessi0_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
        end if
        END FUNCTION bessi0_s

!-------------
        FUNCTION bessk0_s(x)
        USE nrtype; USE nrutil, ONLY : assert,poly
        USE nrface, ONLY : bessi0
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: bessk0_s
        ! Returns the modified Bessel function K0(x) for positive real x.
        REAL(DP) :: y !Accumulate polynomials in double precision.
        REAL(DP), DIMENSION(7) :: p = (/-0.57721566_dp,0.42278420_dp,&
        0.23069756_dp,0.3488590e-1_dp,0.262698e-2_dp,0.10750e-3_dp,&
        0.74e-5_dp/)
        REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,-0.7832358e-1_dp,&
        0.2189568e-1_dp,-0.1062446e-1_dp,0.587872e-2_dp,&
        -0.251540e-2_dp,0.53208e-3_dp/)
        call assert(x > 0.0, "bessk0_s arg")
        if (x <= 2.0) then ! Polynomial fit.
        y=x*x/4.0_sp
        bessk0_s=(-log(x/2.0_sp)*bessi0(x))+poly(y,p)
        else
        y=(2.0_sp/x)
        bessk0_s=(exp(-x)/sqrt(x))*poly(y,q)
        end if
        END FUNCTION bessk0_s
! __________________________________
        FUNCTION bessi1_s(x)
        USE nrtype; USE nrutil, ONLY : poly
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: bessi1_s
        ! Returns the modified Bessel function I1(x) for any real x.
        REAL(SP) :: ax
        REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
        0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
        0.301532e-2_dp,0.32411e-3_dp/)
        ! Accumulate polynomials in double precision.
        REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
        -0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
        0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
        -0.420059e-2_dp/)
        ax=abs(x)
        if (ax < 3.75) then ! Polynomial fit.
        bessi1_s=ax*poly(real((x/3.75_sp)**2,dp),p)
        else
        bessi1_s=(exp(ax)/sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
        end if
        if (x < 0.0) bessi1_s=-bessi1_s
        END FUNCTION bessi1_s
!----------------------------
        FUNCTION bessk1_s(x)
        USE nrtype; USE nrutil, ONLY : assert,poly
        USE nrface, ONLY : bessi1
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: x
        REAL(SP) :: bessk1_s
        ! Returns the modified Bessel function K1(x) for positive real x.
        REAL(DP) :: y ! Accumulate polynomials in double precision.
        REAL(DP), DIMENSION(7) :: p = (/1.0_dp,0.15443144_dp,&
        -0.67278579_dp,-0.18156897_dp,-0.1919402e-1_dp,&
        -0.110404e-2_dp,-0.4686e-4_dp/)
        REAL(DP), DIMENSION(7) :: q = (/1.25331414_dp,0.23498619_dp,&
        -0.3655620e-1_dp,0.1504268e-1_dp,-0.780353e-2_dp,&
        0.325614e-2_dp,-0.68245e-3_dp/)
        call assert(x > 0.0, "bessk1 arg")
        if (x <= 2.0) then ! Polynomial fit.
        y=x*x/4.0_sp
        bessk1_s=(log(x/2.0_sp)*bessi1(x))+(1.0_sp/x)*poly(y,p)
        else
        y=2.0_sp/x
        bessk1_s=(exp(-x)/sqrt(x))*poly(y,q)
        end if
        END FUNCTION bessk1_s


!        Returns the modified Bessel function Kn(x) for positive x and n ≥ 2.
        FUNCTION bessk_s(n,x)
        USE nrtype; USE nrutil, ONLY : assert
        USE nrface, ONLY : bessk0,bessk1
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN) :: n
        REAL(sp), INTENT(IN) :: x
        REAL(sp) :: bessk_s
        INTEGER(I4B) :: j
        REAL(sp) :: bk,bkm,bkp,tox
        call assert(n >= 2, x > 0.0, "bessk args")
        tox=2.0_sp/x
        bkm=bessk0(x) ! Upward recurrence for all x...
        bk=bessk1(x)
        do j=1,n-1 ! ...and here it is.
        bkp=bkm+j*tox*bk
        bkm=bk
        bk=bkp
        end do
        bessk_s=bk
        END FUNCTION bessk_s
