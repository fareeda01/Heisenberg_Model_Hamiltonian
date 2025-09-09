FUNCTION ValCoupl(debug, nsites, i, j, det1, det2, coupling) RESULT(val)
    implicit none
    logical*1, intent(in) :: debug
    logical*1, intent(in) :: det1(nsites), det2(nsites)
    LOGICAL*1 :: det3(nsites)
    integer, intent(in) :: i, j, nsites
    real*8, intent(in) :: coupling(nsites, nsites)
    integer :: k, l, diffspin, position_diff(2)
    real*8 :: val

    ! Initialize variables
    val = 0.0
    diffspin = 0

    ! **Case 1: Diagonal term (i == j)**
    if (i == j) then
        if (debug) then
            write(6,*) 'Diagonal term: Calculating energy'
        endif
        do k = 1, nsites - 1
            do l = k+1, nsites
                if (det1(k) .neqv. det2(l)) then
                        val = val + coupling(k, l)  ! **Parallel spins contribute +J**
                endif
            end do
        end do
    else 
        do k=1,nsites
        det3(k)=xor(det1(k),det2(k))
        if (det3(k)) then
            diffspin=diffspin+1
            if (diffspin .le. 2) then
                position_diff(diffspin)=k
            endif
        endif
        enddo

    ! **Case 2: Off-diagonal term (Spin exchange)**
        if(debug) then
            write(6,*)
            write(6,*) 'det1', det1, 'det2', det2, 'xor', det3
           endif
           if(mod(diffspin,2).eq.1)then
            else if (diffspin.eq.0)then
                    write(6,*) 'Error, this is a diagonal term'
            else if (diffspin.lt.2) then
                    write (6,*) 'More than 2 differences, no coupling'
            else 
                endif
              endif
            if (diffspin.eq.2)then
                val=-coupling(position_diff(1),position_diff(2))
             endif
End FUNCTION ValCoupl
