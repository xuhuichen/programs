!
!
!
!     This subroutine initializes the heating rates and 
!     dispersion coefficients for Coulomb scattering off of 
!     protons and electrons (i.e., Moeller scattering).  These
!     quantities are read in from files and stored in the 
!     arrays lib_dg_ce, lib_disp_ce, lib_dg_cp, and lib_disp_cp.
!     The gamma-1 grid from the files is stored in rate_gm1.
!                        JDF Sept. 2004
      subroutine Coulomb
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
      logical ex
!
      integer i, j, k
      double precision avg
      double precision  kTe, kTp
      integer          num_steps
!
      character *30 name_dge, name_dgp
!
!
      data name_dge/'rates/dge0005.dat'/
      data name_dgp/'rates/dgp0005.dat'/
!
      rate_gm1(1)  = 0.100178E+01
      rate_gm1(2)  = 0.100200E+01
      rate_gm1(3)  = 0.100224E+01
      rate_gm1(4)  = 0.100252E+01
      rate_gm1(5)  = 0.100282E+01
      rate_gm1(6)  = 0.100317E+01
      rate_gm1(7)  = 0.100356E+01
      rate_gm1(8)  = 0.100399E+01
      rate_gm1(9)  = 0.100448E+01
      rate_gm1(10) = 0.100502E+01
      rate_gm1(11) = 0.100564E+01
      rate_gm1(12) = 0.100632E+01
      rate_gm1(13) = 0.100709E+01
      rate_gm1(14) = 0.100796E+01
      rate_gm1(15) = 0.100893E+01
      rate_gm1(16) = 0.101002E+01
      rate_gm1(17) = 0.101124E+01
      rate_gm1(18) = 0.101262E+01
      rate_gm1(19) = 0.101415E+01
      rate_gm1(20) = 0.101588E+01
      rate_gm1(21) = 0.101782E+01
      rate_gm1(22) = 0.101999E+01
      rate_gm1(23) = 0.102243E+01
      rate_gm1(24) = 0.102517E+01
      rate_gm1(25) = 0.102824E+01
      rate_gm1(26) = 0.103169E+01
      rate_gm1(27) = 0.103555E+01
      rate_gm1(28) = 0.103989E+01
      rate_gm1(29) = 0.104475E+01
      rate_gm1(30) = 0.105021E+01
      rate_gm1(31) = 0.105634E+01
      rate_gm1(32) = 0.106321E+01
      rate_gm1(33) = 0.107093E+01
      rate_gm1(34) = 0.107958E+01
      rate_gm1(35) = 0.108929E+01
      rate_gm1(36) = 0.110018E+01
      rate_gm1(37) = 0.111240E+01
      rate_gm1(38) = 0.112612E+01
      rate_gm1(39) = 0.114150E+01
      rate_gm1(40) = 0.115877E+01
      rate_gm1(41) = 0.117814E+01
      rate_gm1(42) = 0.119987E+01
      rate_gm1(43) = 0.122425E+01
      rate_gm1(44) = 0.125161E+01
      rate_gm1(45) = 0.128231E+01
      rate_gm1(46) = 0.131675E+01
      rate_gm1(47) = 0.135539E+01
      rate_gm1(48) = 0.139875E+01
      rate_gm1(49) = 0.144740E+01
      rate_gm1(50) = 0.150198E+01
      rate_gm1(51) = 0.156322E+01
      rate_gm1(52) = 0.163194E+01
      rate_gm1(53) = 0.170903E+01
      rate_gm1(54) = 0.179553E+01
      rate_gm1(55) = 0.189259E+01
      rate_gm1(56) = 0.200148E+01
      rate_gm1(57) = 0.212367E+01
      rate_gm1(58) = 0.226075E+01
      rate_gm1(59) = 0.241456E+01
      rate_gm1(60) = 0.258714E+01
      rate_gm1(61) = 0.278077E+01
      rate_gm1(62) = 0.299803E+01
      rate_gm1(63) = 0.324179E+01
      rate_gm1(64) = 0.351528E+01
      rate_gm1(65) = 0.382215E+01
      rate_gm1(66) = 0.416645E+01
      rate_gm1(67) = 0.455276E+01
      rate_gm1(68) = 0.498620E+01
      rate_gm1(69) = 0.547251E+01
      rate_gm1(70) = 0.601816E+01
      rate_gm1(71) = 0.663037E+01
      rate_gm1(72) = 0.731728E+01
      rate_gm1(73) = 0.808799E+01
      rate_gm1(74) = 0.895272E+01
      rate_gm1(75) = 0.992295E+01
      rate_gm1(76) = 0.110116E+02
      rate_gm1(77) = 0.122330E+02
      rate_gm1(78) = 0.136034E+02
      rate_gm1(79) = 0.151410E+02
      rate_gm1(80) = 0.168662E+02
      rate_gm1(81) = 0.188019E+02
      rate_gm1(82) = 0.209737E+02
      rate_gm1(83) = 0.234105E+02
      rate_gm1(84) = 0.261446E+02
      rate_gm1(85) = 0.292122E+02
      rate_gm1(86) = 0.326541E+02
      rate_gm1(87) = 0.365159E+02
      rate_gm1(88) = 0.408488E+02
      rate_gm1(89) = 0.457104E+02
      rate_gm1(90) = 0.511651E+02
      rate_gm1(91) = 0.572852E+02
      rate_gm1(92) = 0.641520E+02
      rate_gm1(93) = 0.718565E+02
      rate_gm1(94) = 0.805010E+02
      rate_gm1(95) = 0.902002E+02
      rate_gm1(96) = 0.101083E+03
      rate_gm1(97) = 0.113293E+03
      rate_gm1(98) = 0.126992E+03
      rate_gm1(99) = 0.142363E+03
      rate_gm1(100)= 0.159610E+03
!
!     Find the min and max electron and proton temperatures 
!     of all the zones.  Let tea go from 50 below to 50 above
!     the min and max tea.
      tea_min=1.e10
      tea_max=0.
      tna_min=1.e10
      tna_max=0
      do 10 k = 1, nr
         do 11 j = 1, nz
            if( tea_min.gt.tea(j,k) )  tea_min=tea(j,k)
            if( tea_max.lt.tea(j,k) )  tea_max=tea(j,k)
            if( tna_min.gt.tna(j,k) )  tna_min=tna(j,k)
            if( tna_max.lt.tna(j,k) )  tna_max=tna(j,k)
 11      continue
 10   continue
      tea_max=tea_max+50.
      if( tea_max.ge.750.) then
         write(*,*) 'In subroutine Coulomb, '
         write(*,*) 'electron temperature is too high!'
         stop
      endif
      if( tea_min.le.55. )  then 
         tea_min=5.
      else 
         tea_min=tea_min-50.
      endif
!
!      write(*,*) 'tea_min:  ', tea_min
!      write(*,*) 'tea_max:  ', tea_max
!      write(*,*) 'tna_min:  ', tna_min
!      write(*,*) 'tna_max:  ', tna_max
!
!     Loop over all electron temperatures
      num_steps=int(tea_max-tea_min)
      if(num_steps.gt.num_temp_max) num_steps=num_temp_max

!      num_steps=11
      do 20 i=1, num_steps
         kTe=i+tea_min
         call rate_file_name(kTe, name_dge)
!         write(4,*) 'kTe:  ', kTe
!         write(4,*) 'name_dge:  ', name_dge
         inquire(file=name_dge, EXIST=ex)
        if(ex) then
           call read_rate_file(lib_dg_ce, lib_disp_ce, name_dge, i)
        else
!          The marker for not having rates for a certain
!          temperature is the first rate is equal to zero.
           lib_dg_ce(1,i)=0.0
           lib_disp_ce(1,i)=0.0
        endif
!        write(4,*) 'lib_dg_ce(1,i),i:  ', lib_dg_ce(1,i), i
 20   continue
!
!
!     This if round tna_min to the nearest 10 KeV.
      if(int(tna_min/10.d0).eq.int(tna_min/10.d0+.5d0)) then
         tna_min=int(tna_min/10)*10.d0
      else
         tna_min=int(tna_min/10)*10.d0 + 10.d0
      endif
      num_steps=nr*nz
      if(num_steps.gt.num_temp_max) num_steps=num_temp_max
!      num_steps=5
!     loop over the proton temperatures in all zones
      do 22 j=1, nz
         do 21 k=1, nr
            i = k+(j-1)*nr
            if( i.gt.num_steps ) goto 23
!           This rounds tna to the nearest 10 KeV and saves
!           it as kTp
            if(int(tna(j,k)/10.d0).eq.int(tna(j,k)/10.d0+.5d0)) then
               kTp=int(tna(j,k)/10)*10.d0
            else
               kTp=int(tna(j,k)/10)*10.d0 + 10.d0
            endif
            call rate_file_name(kTp, name_dgp)
!            write(4,*) 'kTp:  ', kTp
!            write(4,*) 'tna:  ', tna(j,k)
!            write(4,*) 'name_dgp:  ', name_dgp
!            write(4,*) 'jk code number:  ', k+(j-1)*nr
            inquire(file=name_dgp, EXIST=ex)
            if(ex) then
               call read_rate_file(lib_dg_cp, lib_disp_cp, name_dgp, i)
            else
!              The marker for not having rates for a certain
!              temperature is the rate is equal to zero.
               lib_dg_cp(1,i)=0.0
               lib_disp_cp(1,i)=0.0
            endif
 21      continue
 22   continue
 23   continue
!
!
      return
      end
!
!
!     This subroutine gets the name of a rate file,
!     i.e., a dgeXXXX.dat or dgpXXXX.dat file.
      subroutine rate_file_name(kT, name)
!
      double precision kT
      integer temp_i, i_5, i_10000, i_1000, i_100, i_10, i_1
      character *30 name
!

!
      if(name(9:9).eq.'e') then
!        create electron file name
         temp_i = int(kT)
         i_1000 = (temp_i-1000)/1000
         name(10:10) = char(48+i_1000)
         i_100 = (temp_i - 1000*i_1000)/100
         name(11:11) = char(48+i_100)
         i_10 = (temp_i - 1000*i_1000 - 100*i_100)/10
         name(12:12) = char(48+i_10)
         i_1 = temp_i - 1000*i_1000 - 100*i_100 - 10*i_10
         name(13:13) = char(48+i_1)
      else
!        create proton file name
         temp_i = int(kT/10)
         i_5 = temp_i/10000
         if(i_5.gt.0) then
            name(9:9) = char(48+i_5)
         else
            name(9:9) = 'p'
         endif
         i_10000 = (temp_i - 10000*i_5)/1000
         name(10:10) = char(48+i_10000)
         i_1000 = (temp_i - 10000*i_5 - 1000*i_10000)/100
         name(11:11) = char(48+i_1000)
         i_100 = (temp_i - 10000*i_5 - 1000*i_10000 - 100*i_1000)/10
         name(12:12) = char(48+i_100)
         i_10 = temp_i -10000*i_5 - 1000*i_10000 - 100*i_1000 - 10*i_100
         name(13:13) = char(48+i_10)
      endif
      return
      end
!
!
!
!     This subroutine reads in values for dg and disp
!     (which represent the heating/cooling rate and 
!     dispersion coefficient, respectively)
!     from a rate file (dgeXXXX.dat or dgpXXXX.dat file.
      subroutine read_rate_file(dg, disp, name, i)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!
      integer i, m
!     i is the index for temperature, m is the index
!     for energy bin
      double precision g_read, dg(num_nt, num_temp_max), 
     1                 disp(num_nt, num_temp_max)
      character *30 name
!     format for reading from rate file
 65   format(e14.7,1x,e14.7,1x,e14.7)
!
      open(unit=17, file=name, status='unknown')
      do 250 m = 1, num_nt
         read(17, 65) g_read, dg(m,i), disp(m,i)
 250  continue
      close(17)

      return
      end
!
!
!
