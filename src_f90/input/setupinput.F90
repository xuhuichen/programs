program main
 implicit none
 character(len=79) buftext(20)
 character(len=81) buffer
 integer mhdtimestepmax,nr,nz
 integer ntime,tbbu(1,9),tbbl(1,9),spec_switch,nphreg,nmu,nph_lc,nphbins(9),nst,rseed,rand_switch,&
    cr_sent,upper_sent,dh_sentinel,pair_switch,T_const,cf_sentinel,&
    inj_switch,g2var_switch,esc_sw,pick_sw,split2,split3
 integer i,j,k, m, t, nin
 character(len=79) filename
 character(len=2) xj,xk
 character(len=3) xi
 character(len=30) u_fname(1,9),l_fname(1,9),spname,phname,lcname(1),eventfile,temp_file
 double precision tea,tna,n_e,B_input,B_field,theta_b,amxwl,gmin,gmax,p_nth,q_turb,turb_lev,r_acc,&
                  acc_thick, acc_rwide, z_acc_dis
 double precision z(99),r(99),rmin,dz,dr,star_switch,Rstar,dist_star,tstop,mcdt,t0(9),t1(9),Ephmin(9),Ephmax(9),&
    Elcmin(9),Elcmax(9),&
    r_flare,z_flare,t_flare,sigma_r,sigma_z,sigma_t,flare_amp,r_esc,inj_g1,inj_g2,inj_p,inj_t,inj_L,dummy,&
    inj_gg,inj_sigma,g_bulk,inj_v,R_blr,fr_blr,R_ir,fr_ir,R_disk,d_jet,split1,spl3_trg
 integer ep_switch
 integer inhom_acc  ! switch to control whether to use inhomogeneous acceleration

 double precision mhdtimestep,currenttimestop
 double precision constanta,constantc,constanttau
 double precision,parameter :: pi2=6.2831853072d0
 double precision,parameter :: c_light=2.9979245620d10

!       Input formats
!
 100  format(a80,i2)
 101  format(a80,i3)
 102  format(a80,i10) 
 110  format(a80,e14.7)
 112  format(a80,a30)
 
 buftext(14)=" 0. Current MHD time stop                                     currenttimestop = "
 buftext(1) =" 1. Electron temperature in zone [keV]                                    tea = "
 buftext(2) =" 2. Proton temperature in zone [keV]                                      tna = "
 buftext(3) =" 3. Particle density in zone [cm^(-3)]                                    n_e = "
 buftext(4) =" 4. Magn. field sent. (0 = spec., 1 = ep. w. el., 2 = ep. w. pr.)   ep_switch = "
 buftext(5) =" 5. Magnetic field [G]                                                B_field = "
 buftext(6) =" 6. The angle between the helical magnetic field and z                theta_b = "
 buftext(7) =" 7. Maxwellian fraction in zone                                         amxwl = "
 buftext(8) =" 8. Low-energy cut-off of nonthermal electron population:                gmin = "
 buftext(9) =" 9. High-energy cut-off of nonthermal electron population:               gmax = "
 buftext(10)=" 10. Non-thermal electron spectral index:                               p_nth = "
 buftext(11)=" 11. Turbulence sp.index                                          q_turb(j,k) = "
 buftext(12)=" 12. Turbulence Level (deltaB/B0)^2                            turb_lev (j,k) = "
 buftext(13)=" 13. Acceleratin time scale (Z/c)                                 r_acc (j,k) = "
 
 constanta=5.5d0/81.d0
 constantc=1.d-2
 constanttau=2.1d0

 tea=1.0d1
 tna=1.0d2
 n_e=3.d01
 ep_switch=0
 B_input=4.d-2
 theta_b=45.0d0
 amxwl=0.0d0
 gmin=8.d1
 gmax=1.d3
 p_nth=4.0d0
 q_turb=1.66666667d0
 turb_lev=1.0d-20
 r_acc=2.34d-4
 
! thinckness of the acceleration region
 acc_thick=2.0/20
 acc_rwide=2.0/15
 z_acc_dis=10.0/20
 inhom_acc = 0
! 0 is homogeneous, 1 is bottom acceleration, 2 is shear accleration at outer layer,
! 3 is acceleration at center zones
! 21 is with outer 1/3 no B. 22 is 21 with shear acceleration at outer layer
      nin=5
      open(unit=nin, file='input.dat', status='unknown')
      read(nin, 101) buffer, nz
      read(nin, 101) buffer, nr
      read(nin, 110) buffer, z(nz)
      read(nin, 110) buffer, rmin
      read(nin, 110) buffer, r(nr)
      dz=z(nz)/nz
      dr=r(nr)/nr
      read(nin, 101) buffer, star_switch
      if(star_switch.eq.1) then
         read(nin, 110) buffer, Rstar
         read(nin, 110) buffer, dist_star
      endif
      read(nin, 110) buffer, tstop
      read(nin, 110) buffer, mcdt
      read(nin, 101) buffer, ntime
      do t = 1, ntime
         read(nin, 110) buffer, t0(t)
         read(nin, 110) buffer, t1(t)
         read(nin, 110) buffer, tbbu(1,t)
         read(nin, 112) buffer, u_fname(1, t)
         read(nin, 110) buffer, tbbl(1,t)
         read(nin, 112) buffer, l_fname(1, t)
      enddo
      read(nin,100) buffer, spec_switch
      read(nin,100) buffer, nphreg
      do m = 1, nphreg
         read(nin,110) buffer, Ephmin(m)
         read(nin,110) buffer, Ephmax(m)
         read(nin,101) buffer, nphbins(m)
      enddo
      read(nin,101) buffer, nmu
      read(nin,101) buffer, nph_lc
      do m = 1, nph_lc
         read(nin,110) buffer, Elcmin(m)
         read(nin,110) buffer, Elcmax(m)
      enddo
      read(nin,112) buffer, spname
      read(nin,112) buffer, phname
      read(nin,112) buffer, lcname(1)
      read(nin,112) buffer, eventfile
      read(nin,112) buffer, temp_file
      read(nin,102) buffer, nst
      read(nin,102) buffer, rseed
      read(nin,102) buffer, rand_switch
      read(nin,100) buffer, cr_sent
      read(nin,100) buffer, upper_sent
      read(nin,100) buffer, dh_sentinel
      read(nin,100) buffer, pair_switch
      read(nin,100) buffer, T_const
      read(nin,100) buffer, cf_sentinel
      read(nin,110) buffer, r_flare
      read(nin,110) buffer, z_flare
      read(nin,110) buffer, t_flare
      read(nin,110) buffer, sigma_r
      read(nin,110) buffer, sigma_z
      read(nin,110) buffer, sigma_t
      read(nin,110) buffer, flare_amp
      read(nin,100)buffer,inj_switch
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
      close(nin)


 
 write(*,*) "mhdtimestepmax"
 read(*,*) mhdtimestepmax
 mhdtimestep=tstop/dble(mhdtimestepmax)

!   example minput file
    open(unit=10,file='minput.dat',access="sequential",status="unknown",action="readwrite")
    write(10,"(A80,e14.7)") buftext(14),currenttimestop
    write(10,"(A80,e14.7)") buftext(1),tea
    write(10,"(A80,e14.7)") buftext(2),tna
    write(10,"(A80,e14.7)") buftext(3),n_e
    write(10,"(A80,i2)")    buftext(4),ep_switch
    write(10,"(A80,e14.7)") buftext(5),B_input
    write(10,"(A80,e14.7)") buftext(6),theta_b
    write(10,"(A80,e14.7)") buftext(7),amxwl
    write(10,"(A80,e14.7)") buftext(8),gmin
    write(10,"(A80,e14.7)") buftext(9),gmax
    write(10,"(A80,e14.7)") buftext(10),p_nth
    write(10,"(A80,e14.7)") buftext(11),q_turb
    write(10,"(A80,e14.7)") buftext(12),turb_lev
    write(10,"(A80,e14.7)") buftext(13),r_acc
    close(10)
!---------------------------------------------
 filename="input_00_00_000.dat"
 do i=1,mhdtimestepmax
  write(xi,"(i3.3)") i
  filename(13:15)=trim(xi)
  currenttimestop=dble(i)*mhdtimestep
  do j=1,nz
   do k=1,nr
    write(xj,"(i2.2)") j
    write(xk,"(i2.2)") k
    filename(7:8)=trim(xj)
    filename(10:11)=trim(xk)
    open(unit=10,file=filename,access="sequential",status="unknown",action="readwrite")
!   B_field=constantc*sqrt(4.d0+4.d0*pi2*(j-0.5d0)**4.d0)*exp(-constanta*(j-0.5d0)**2.d0)
!   theta_b=atan((4.d0*pi2)**0.25d0*(j-0.5d0)/(1.d0-sqrt(pi2)*(j-0.5d0)*(j-0.5d0)))/3.141592657d0*180.d0
!   r_acc=constanttau/(1.d0+exp(-constanta*(j-0.5d0)**2.d0)*sqrt(1.d0+(sqrt(pi2)/2.d0-3.d0*constanta+0.5d0/sqrt(pi2)*constanta**2.d0)* &
!         (j-0.5d0)**2.d0+(sqrt(pi2)*constanta+3.d0*constanta**2.d0)*(j-0.5d0)**4.d0+sqrt(pi2)/2.d0*constanta**2.d0*(j-0.5d0)**6.d0))
     if(inj_switch.eq.23 .and. (currenttimestop-mhdtimestep-inj_t).gt.(j-1)*dz/inj_v.and.&
       (currenttimestop-mhdtimestep-inj_t).lt.(j-1+int(acc_thick*nz+0.5d0))*dz/inj_v)then
       B_field = B_input*dsqrt(flare_amp)
       write(*,*)'Amplified B_field=',B_field,'j=',j,'k=',k,'i=',i
     else
       B_field = B_input
     endif
    write(10,"(e14.7)") currenttimestop
    write(10,"(e14.7)") tea
    write(10,"(e14.7)") tna
    if(inhom_acc/10.eq.2.and. (j.le.nz/6.or.j.gt.nz*5/6.or.k.gt.nr*2/3))then
      write(10,"(e14.7)") n_e/2.d0
    else
      write(10,"(e14.7)") n_e
    endif
    write(10,"(i2)")    ep_switch
    if(inhom_acc/10.eq.2.and. (j.le.nz/6.or.j.gt.nz*5/6.or.k.gt.nr*2/3))then
      write(10,"(e14.7)") 0.d0
    else
      write(10,"(e14.7)") B_field
    endif
    write(10,"(e14.7)") theta_b
    write(10,"(e14.7)") amxwl
    write(10,"(e14.7)") gmin
    write(10,"(e14.7)") gmax
    write(10,"(e14.7)") p_nth
    write(10,"(e14.7)") q_turb
    write(10,"(e14.7)") turb_lev
    if(inhom_acc.eq.1.and.j.gt.int(acc_thick*nz+0.5d0))then
      write(10,"(e14.7)") 1.d40
    elseif(inhom_acc.eq.2.and. k.le.nr-int(acc_thick*nz+0.5d0))then
      write(10,"(e14.7)") 1.d40
    elseif(inhom_acc.eq.3.and. (k.gt.int(acc_rwide*nr+0.5).or.&
           j.lt.int(1+(z_acc_dis-0.5*acc_thick)*nz+0.4).or.&
           j.gt.int(1+(z_acc_dis+0.5*acc_thick)*nz-0.4)))then
      write(10,"(e14.7)") 1.d40
    elseif(inhom_acc.eq.21.and. (j.le.nz/6.or.j.gt.nz/6+int(acc_thick*nz+0.5d0).or.k.gt.nr*2/3))then
      write(10,"(e14.7)") 1.d40
    elseif(inhom_acc.eq.22.and. (j.le.nz/6.or.j.gt.nz*5/6.or.k.ne.nr*2/3))then
      write(10,"(e14.7)") 1.d40
    else
      write(10,"(e14.7)") r_acc
    endif
    close(10)
   end do
  end do
 end do
 
end program
