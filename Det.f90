      recursive SUBROUTINE deter(n,na,tab,ndeter,accu,compt)

!cc This subroutine generate, recursively, all the determinants
!cc           of a given number of sites and value of ms
!cc n is the number of sites
!cc na is the number of alpha spins
!cc tab contain the determinants of dimension n (m "1" values and n-m "0" values)
!cc ndeter is the number of determinants
!cc accu and accu1 are accumulator variables
!cc compt is a counter of the number of determinants generated 

      integer, intent(in) :: n,na,ndeter
      integer, intent(out) :: tab(1:ndeter),accu,compt
      integer  :: accu1

	
      if (na.gt.n) then
       write(6,*) 'Impossible'
       stop
      endif
      
      if (na == 0) then
       accu=accu+0
       compt=compt+1
       tab(compt)=accu
      else if (n == na) then
       do i=0,n-1
        accu=accu+2**i
       enddo
       compt=compt+1
       tab(compt)=accu

      else
       accu1=accu
       accu=accu+2**(n-1)
       call deter(n-1,na-1,tab,ndeter,accu,compt)
       
       call deter(n-1,na,tab,ndeter,accu1,compt)       
      endif     

      end SUBROUTINE deter

