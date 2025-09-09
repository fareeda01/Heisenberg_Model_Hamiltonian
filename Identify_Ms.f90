SUBROUTINE identify_Ms(debug, nms_spaces, size_ms_subspace, energies)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: debug
    INTEGER, INTENT(IN) :: nms_spaces, size_ms_subspace(nms_spaces)
    REAL*8, INTENT(INOUT) :: energies(nms_spaces+1, 0:maxval(size_ms_subspace))
    REAL*8, ALLOCATABLE :: mstab1(:), mstab2(:), tmpmstab(:)
    REAL*8 :: diff
    INTEGER :: i, j, k, l
    LOGICAL :: found

    ! Allocate array for S values
    ALLOCATE(tmpmstab(size_ms_subspace(nms_spaces)))
    diff = 1.0D-8

    ! Initialize S values
    tmpmstab = energies(1, 0)

    ! Loop over Ms spaces
    DO i = 1, nms_spaces - 1
        ALLOCATE(mstab1(size_ms_subspace(i)))
        ALLOCATE(mstab2(size_ms_subspace(i+1)))

        ! Initialize first space
        IF (i == 1) THEN
            mstab1(1) = energies(i, 0)
        ELSE
            mstab1 = tmpmstab(1:size_ms_subspace(i))
        END IF

        k = 1
        DO j = 1, size_ms_subspace(i)
            found = .FALSE.
            DO WHILE (.NOT. found)
                IF (.NOT. QEqual(diff, energies(i,j), energies(i+1,k))) THEN
                    mstab2(k) = energies(i+1,0)
                    IF (debug) WRITE(6,*) 'Different energies, k', k
                ELSE
                    mstab2(k) = mstab1(j)
                    found = .TRUE.
                    IF (debug) THEN
                        WRITE(6,*) 'Same energies, k', k
                        WRITE(6,*) 'S value ', mstab1(j)
                    END IF
                END IF
                k = k + 1
            END DO
        END DO

        ! Assign remaining S values
        IF (k <= size_ms_subspace(i+1)) THEN
            DO l = k, size_ms_subspace(i+1)
                mstab2(l) = energies(i+1,0)
            END DO
        END IF

        tmpmstab(1:size_ms_subspace(i+1)) = mstab2

        DEALLOCATE(mstab1, mstab2)
    END DO

    ! Assign final S values
    energies(nms_spaces+1, 1:size_ms_subspace(nms_spaces)) = tmpmstab

    DEALLOCATE(tmpmstab)

97    format(a30,F8.3)
98    format(a16,i3,a5,F8.1,a6,i3)
99    format(F5.1,a2)

CONTAINS

    FUNCTION QEqual(tolerance, value1, value2)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: tolerance, value1, value2
        LOGICAL :: QEqual

        QEqual = ABS(value1 - value2) <= tolerance
    END FUNCTION QEqual

END SUBROUTINE identify_Ms

