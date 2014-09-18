      subroutine graphics
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
!lc	This routine is called at the end of the problem by routine XEC
!lc	It prepares and prints information that may be read and displayed
!lc	by GNUPLOT (a graphing and plotting program)
!
!
      integer nunit_sp, nunit_ph, nunit_lc, lgp,&
     &        n_format, n, rl, r1, r2, nunit_temp,&
     &        n5_sentinel, n5, nt2, n10, n100, n1000
      integer j, k
!

      double precision delnu, foutsp(nmumax), foutpd(nmumax), &
     &                 Ebar, time0
      double precision dN, sum_N, dNT, sum_NT, T_av
      double precision dN_r(kmax), dNT_r(kmax), sumN_r(kmax)
      double precision sumNT_r(kmax), tav_r(kmax)

!
      character *30 fname
      character *30 form(kmax+1), t_profile
!

!
!     Format strings
!
      do k = 1,kmax+1
        write(form(k), '("(e14.6,"i2"(1x,e14.6))")'),kmax
      enddo
       write(*,*)'test, form 1 and 2', form(1), form(2)
!      data form(1) / '(e14.6,1x,e14.6)' / 
!      data form(2) / '(e14.6,2(1x,e14.6))' /
!      data form(3) / '(e14.6,3(1x,e14.6))' /
!      data form(4) / '(e14.6,4(1x,e14.6))' /
!      data form(5) / '(e14.6,5(1x,e14.6))' /
!      data form(6) / '(e14.6,6(1x,e14.6))' /
!      data form(7) / '(e14.6,7(1x,e14.6))' /
!      data form(8) / '(e14.6,8(1x,e14.6))' /
!      data form(9) / '(e14.6,9(1x,e14.6))' /
!      data form(10) / '(e14.6,10(1x,e14.6))' /
!      data form(11) / '(e14.6,11(1x,e14.6))' /
!      data form(12) / '(e14.6,12(1x,e14.6))' /
!      data form(13) / '(e14.6,13(1x,e14.6))' /
!      data form(14) / '(e14.6,14(1x,e14.6))' /
!      data form(15) / '(e14.6,15(1x,e14.6))' /
!      data form(16) / '(e14.6,16(1x,e14.6))' /
!      data form(17) / '(e14.6,17(1x,e14.6))' /
!      data form(18) / '(e14.6,18(1x,e14.6))' /
!      data form(19) / '(e14.6,19(1x,e14.6))' /
!      data form(20) / '(e14.6,20(1x,e14.6))' /
!      data form(21) / '(e14.6,21(1x,e14.6))' /
!      data form(22) / '(e14.6,22(1x,e14.6))' /
!      data form(23) / '(e14.6,23(1x,e14.6))' /
!      data form(24) / '(e14.6,24(1x,e14.6))' /
!      data form(25) / '(e14.6,25(1x,e14.6))' /
!      data form(26) / '(e14.6,26(1x,e14.6))' /
!      data form(27) / '(e14.6,27(1x,e14.6))' /
!      data form(28) / '(e14.6,28(1x,e14.6))' /
!      data form(29) / '(e14.6,29(1x,e14.6))' /
!      data form(30) / '(e14.6,30(1x,e14.6))' /
!      data form(31) / '(e14.6,31(1x,e14.6))' /
!      data form(32) / '(e14.6,32(1x,e14.6))' /
!      data form(33) / '(e14.6,33(1x,e14.6))' /
!      data form(34) / '(e14.6,34(1x,e14.6))' /
!      data form(35) / '(e14.6,35(1x,e14.6))' /
!      data form(36) / '(e14.6,36(1x,e14.6))' /
!      data form(37) / '(e14.6,37(1x,e14.6))' /
!      data form(38) / '(e14.6,38(1x,e14.6))' /
!      data form(39) / '(e14.6,39(1x,e14.6))' /
!      data form(40) / '(e14.6,40(1x,e14.6))' /
!      data form(41) / '(e14.6,41(1x,e14.6))' /
!      data form(42) / '(e14.6,42(1x,e14.6))' /
!      data form(43) / '(e14.6,43(1x,e14.6))' /
!      data form(44) / '(e14.6,44(1x,e14.6))' /
!      data form(45) / '(e14.6,45(1x,e14.6))' /
!      data form(46) / '(e14.6,46(1x,e14.6))' /
!      data form(47) / '(e14.6,47(1x,e14.6))' /
!      data form(48) / '(e14.6,48(1x,e14.6))' /
!      data form(49) / '(e14.6,49(1x,e14.6))' /
!      data form(50) / '(e14.6,50(1x,e14.6))' /
!      data form(51) / '(e14.6,51(1x,e14.6))' /
!      data form(52) / '(e14.6,52(1x,e14.6))' /
!      data form(53) / '(e14.6,53(1x,e14.6))' /
!      data form(54) / '(e14.6,54(1x,e14.6))' /
!      data form(55) / '(e14.6,55(1x,e14.6))' /
!      data form(56) / '(e14.6,56(1x,e14.6))' /
!      data form(57) / '(e14.6,57(1x,e14.6))' /
!      data form(58) / '(e14.6,58(1x,e14.6))' /
!      data form(59) / '(e14.6,59(1x,e14.6))' /
!      data form(60) / '(e14.6,60(1x,e14.6))' /
!      data form(61) / '(e14.6,61(1x,e14.6))' /
!      data form(62) / '(e14.6,62(1x,e14.6))' /
!      data form(63) / '(e14.6,63(1x,e14.6))' /
!      data form(64) / '(e14.6,64(1x,e14.6))' /
!      data form(65) / '(e14.6,65(1x,e14.6))' /
!      data form(66) / '(e14.6,66(1x,e14.6))' /
!      data form(67) / '(e14.6,67(1x,e14.6))' /
!      data form(68) / '(e14.6,68(1x,e14.6))' /
!      data form(69) / '(e14.6,69(1x,e14.6))' /
!      data form(70) / '(e14.6,70(1x,e14.6))' /
!      data form(71) / '(e14.6,71(1x,e14.6))' /
!      data form(72) / '(e14.6,72(1x,e14.6))' /
!      data form(73) / '(e14.6,73(1x,e14.6))' /
!      data form(74) / '(e14.6,74(1x,e14.6))' /
!      data form(75) / '(e14.6,75(1x,e14.6))' /
!      data form(76) / '(e14.6,76(1x,e14.6))' /
!      data form(77) / '(e14.6,77(1x,e14.6))' /
!      data form(78) / '(e14.6,78(1x,e14.6))' /
!      data form(79) / '(e14.6,79(1x,e14.6))' /
!      data form(80) / '(e14.6,80(1x,e14.6))' /
!      data form(81) / '(e14.6,81(1x,e14.6))' /
!      data form(82) / '(e14.6,82(1x,e14.6))' /
!      data form(83) / '(e14.6,83(1x,e14.6))' /
!      data form(84) / '(e14.6,84(1x,e14.6))' /
!      data form(85) / '(e14.6,85(1x,e14.6))' /
!      data form(86) / '(e14.6,86(1x,e14.6))' /
!      data form(87) / '(e14.6,87(1x,e14.6))' /
!      data form(88) / '(e14.6,88(1x,e14.6))' /
!      data form(89) / '(e14.6,89(1x,e14.6))' /
!      data form(90) / '(e14.6,90(1x,e14.6))' /
!      data form(91) / '(e14.6,91(1x,e14.6))' /
!      data form(92) / '(e14.6,92(1x,e14.6))' /
!      data form(93) / '(e14.6,93(1x,e14.6))' /
!      data form(94) / '(e14.6,94(1x,e14.6))' /
!      data form(95) / '(e14.6,95(1x,e14.6))' /
!      data form(96) / '(e14.6,96(1x,e14.6))' /
!      data form(97) / '(e14.6,97(1x,e14.6))' /
!      data form(98) / '(e14.6,98(1x,e14.6))' /
!      data form(99) / '(e14.6,99(1x,e14.6))' /
!
      data t_profile / 'temp/profile0001.dat' /
!
      nunit_sp = 12
      nunit_ph = 13
      nunit_lc = 14
      nunit_temp = 15
      nt2 = 16
!
      open (nunit_sp,file=spname,STATUS='UNKNOWN')
      open (nunit_ph,file=phname,STATUS='UNKNOWN')
!
      do 10 lgp=1,nphtotal
         delnu   =  hu(lgp+1) - hu(lgp)
         Ebar = 5.d-1*(hu(lgp+1) + hu(lgp))
         do 3 n = 1, nmu
  	     foutsp(n) = fout(n, lgp)/(delnu*(time+dt(1)))
	     foutpd(n) = 6.25d8*foutsp(n)/Ebar
!
	     if (foutsp(n) .le. 1.e-20) foutsp(n) = 1.e-20
	     if (foutpd(n) .le. 1.e-20) foutpd(n) = 1.e-20
 3       continue
!
!         write(*,*) 'hu = ',hu(lgp)
!         write(*, fmt=form(nmu)) hu(lgp), (foutsp(n), n=1,nmu)
!
         write(nunit_sp, fmt=form(nmu)) hu(lgp),(foutsp(n), n=1,nmu)
         write(nunit_sp, fmt=form(nmu)) hu(lgp+1),(foutsp(n), n=1,nmu)
         write(nunit_ph, fmt=form(nmu)) hu(lgp),(foutpd(n), n=1,nmu)
         write(nunit_ph, fmt=form(nmu)) hu(lgp+1),(foutpd(n), n=1,nmu)
!
 10   continue
!
      close(nunit_sp)
      close(nunit_ph)
!
!      write(*,*) 'Spectra written.'
!
!
      time0 = dmax1(1.d-20, time)
      rl = 17*(nph_lc+1)
      r1 = 2*ncycle - 1
      r2 = 2*ncycle
      do 30 n = 1, nmu
         do 20 lgp=1,nph_lc
            delnu   =  Elcmax(lgp) - Elcmin(lgp)
            foutsp(lgp) = edout(n, lgp)/delnu
            if (foutsp(lgp).lt.1.d-20) foutsp(lgp) = 1.d-20
!
!            write(*,*) 'lgp = ',lgp
!            write(*,*) 'lc_flux = ',foutsp(lgp)
!
 20      continue
!
         if (ncycle.eq.1) then
            open(unit=nunit_lc,file=lcname(n))
         else
            open(unit=nunit_lc,file=lcname(n),access='append')

         endif
!         write(nunit_lc, fmt=form(nph_lc), rec=r1) time0, (foutsp(lgp), 
!     1                                 lgp=1,nph_lc)
!         write(nunit_lc, fmt=form(nph_lc), rec=r2) (time + dt(1)), 
!     1                              (foutsp(lgp), lgp=1,nph_lc)
         write(nunit_lc, fmt=form(nph_lc)) time0, (foutsp(lgp), &
     &                                 lgp=1,nph_lc)
         write(nunit_lc, fmt=form(nph_lc)) (time + dt(1)), &
     &                              (foutsp(lgp), lgp=1,nph_lc)
!

         close(nunit_lc)
  30  continue
!
!
      sum_N = 0.d0
      sum_NT = 0.d0
      do 50 k = 1, nr
         sumN_r(k) = 0.d0
         sumNT_r(k) = 0.d0
         do 40 j = 1, nz
            dN = vol(j,k)*n_e(j,k)
            dNT = tea(j,k)*dN
            sum_N = sum_N + dN
            sum_NT = sum_NT + dNT
            sumN_r(k) = sumN_r(k) + dN
            sumNT_r(k) = sumNT_r(k) + dNT
  40     continue
         Tav_r(k) = sumNT_r(k)/sumN_r(k)
  50  continue
      T_av = sum_NT/sum_N
!
!      write(nunit_temp, 60) time0, T_av
!      write(nunit_temp, 60) (time + dt(1)), T_av
  60  format (e14.7,1x,e14.7)
!
       write(nunit_temp, fmt=form(1+nr)) time0, T_av, (tav_r(k), k=1,nr)
       write(nunit_temp, fmt=form(1+nr)) (time + dt(1)), T_av, &
     &                               (tav_r(k), k=1,nr)
!
!
       do 65 j = 1, nz
          do 64 k = 1, nr
             T_sum(j,k) = T_sum(j,k) + dt(2)*tea(j,k)
 64       continue
 65    continue
       time_sum = time_sum + dt(2)
       n5 = ncycle/5
       n5_sentinel = ncycle - 5*n5
!
       if (n5_sentinel.eq.0) then
          n10 = ncycle/10
          n100 = ncycle/100
          n1000 = ncycle/1000
          t_profile(16:16) = char(48 + ncycle - 10*n10)
          t_profile(15:15) = char(48 + n10 - 10*n100)
          t_profile(14:14) = char(48 + n100 - 10*n1000)
          t_profile(13:13) = char(48 + n1000)
          open(nt2, file=t_profile, status='unknown')
!
          do 90 j = nz, 1, -1
             do 70 k = 1, nr
 70          T_sum(j,k) = T_sum(j,k)/time_sum
             if (nr.gt.1) then
                write(nt2, fmt=form(nr-1)) (T_sum(j,k), k=1,nr)
             else
                write(nt2, *) T_sum(j,1)
             endif
             do 80 k = 1, nr
 80          T_sum(j,k) = 0.d0
 90       continue
!
          close(nt2)
      endif
!
      if (n5_sentinel.eq.0) time_sum = 0.d0
!
      return
      end
!
