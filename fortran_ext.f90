! find Z2 rank of integer matrix
SUBROUTINE Z2RANK(RNK, MAT, NR, NC)
! input: MAT(NR,NC) - Z2 matrix
! output: RNK - Z2 rank
! caller must ensure MAT contains only 0 and 1
    INTEGER :: RNK
    INTEGER :: MAT(NR, NC)
    INTEGER :: NR, NC
    ! local variable
    INTEGER :: R, I, J, K
    LOGICAL :: FOUND
    
    R = 1 ! current row index
    DO I = 1, NC ! run through cols
        IF (R > NR) THEN ! row exhausted first
            RNK = NR ! row rank is full
            RETURN ! early return
        END IF
        IF (MAT(R, I) == 0) THEN ! need to find pivot
            FOUND = .FALSE. ! set a flag
            DO K = R + 1, NR ! in the rest of rows
                IF (MAT(K, I) == 1) THEN ! MAT(K,I) nonzero
                    FOUND = .TRUE.
                    EXIT ! exit pivot searching
                END IF
            END DO ! K pivot searching
            IF (FOUND) THEN ! if pivot found in K
                MAT([R, K], :) = MAT([K, R], :) ! raw swap
            ELSE ! if pivot not found
                CYCLE ! down with this col
            END IF
        END IF
        ! pivot has moved to MAT(R, I), perform GE
        DO J = R + 1, NR
            IF (MAT(J, I) == 1) THEN ! if MAT(J,I) nonzero
                MAT(J, I:) = MODULO(MAT(J, I:) + MAT(R, I:), 2)
            END IF
        END DO
        R = R + 1 ! rank inc
    END DO ! I, through cols
    RNK = R - 1 ! return rank
END SUBROUTINE Z2RANK

! ========== for debug use ==========
!SUBROUTINE PRINT_MAT(MAT, NR, NC)
!	INTEGER, INTENT(IN) :: MAT(NR,NC)
!	INTEGER :: NR, NC, I
!	
!	PRINT *, '-----------------------'
!	DO I = 1, NR
!		PRINT *, '[', MAT(I,:), ']'
!	END DO
!END SUBROUTINE PRINT_MAT
!
!PROGRAM MAIN
!	INTEGER :: A(3,3)
!	INTEGER :: R
!	
!	A = RESHAPE([0,1,1,1,0,1,1,1,0],[3,3])
!	CALL Z2RANK(R, A, 3, 3)
!	PRINT *, R
!END PROGRAM MAIN