FUNCTION QEqual(diff, x, y) RESULT(qeq)
    implicit none
    logical :: qeq
    real*8, intent(in) :: diff, x, y

    if (abs(x - y) .le. abs(diff)) then
        qeq = .true.
    else
        qeq = .false.
    end if
END FUNCTION QEqual