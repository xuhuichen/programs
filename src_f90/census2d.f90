      subroutine write_cens(entries, nunit_c)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
      integer entries, nunit_c
      integer i, n, lwad, lwai, write_flag
!
!      
!
    
          lwad = 0
          lwai = 0   
          do 100 n = 1, entries
             write(nunit_c, 200) (dbufout(lwad+i), i=1,6)
             write(nunit_c, 210) (ibufout(lwai+i), i=1,6)
             lwad = lwad + 6
             lwai = lwai + 6
 100      continue
!
200   format (6e14.7)
210   format (6i5) 
!
      return
      end
!
!
!
      subroutine read_cens(entries,nunit_c)
      implicit none
      include 'general.pa'
      include 'commonblock.f90'
!
      integer entries, nunit_c
      integer i, n, lwad, lwai
!
!
!
!      write(*,*) 'in read_cens'
          lwad = 0
          lwai = 0
          do 100 n = 1, entries
             read(nunit_c, 200) (dbufout(lwad+i), i=1,6)
             read(nunit_c, 210) (ibufout(lwai+i), i=1,6)
         lwad = lwad + 6
         lwai = lwai + 6
 100      continue
!
200   format (6e14.7)
210   format (6i5)
!
      return
      end
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Tue Jun 13 14:50:09 EDT 2006
! version: 2
! Name: J. Finke
! Removed variable nmax which is not used.  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccc 
! Fri Jun 16 12:24:11 EDT 2006
! version: 3
! Name: J. Finke
! fibran seed is now stored in census files. 
! imcfield is now determinable.     
!cccccccccccccccccccccccccccccccccccccccccccccccccc
! Tue Sep 23 15:47 CDT 2008
! Name: Xuhui Chen
! When read_flag and write_flag is 0, census file
! read and write is turned off. The information
! is just stored in buffer.
