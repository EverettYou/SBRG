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
! Pauli index processor
MODULE IPU
    ! input - A, B, output - C
    INTEGER, ALLOCATABLE :: A(:), B(:), C(:)
    INTEGER :: NA, NB, NC
CONTAINS
    SUBROUTINE SETAB(A0, B0, NA0, NB0)
        INTEGER :: A0(NA0), B0(NB0)
        INTEGER :: NA0, NB0
        
        NA = NA0
        IF (ALLOCATED(A)) THEN
            IF (SIZE(A) >= NA) THEN
                A(:NA) = A0
            ELSE
                DEALLOCATE(A)
                A = A0
            END IF
        ELSE
            A = A0
        END IF
        NB = NB0
        IF (ALLOCATED(B)) THEN
            IF (SIZE(B) >= NB) THEN
                B(:NB) = B0
            ELSE
                DEALLOCATE(B)
                B = B0
            END IF
        ELSE
            B = B0
        END IF
    END SUBROUTINE SETAB
    ! count num of duplicates between two sorted arrays
    SUBROUTINE OVERLAP(DUPS)
        INTEGER :: DUPS
        ! local variables
        INTEGER :: IA, IB, DIFF
        
        DUPS = 0
        IA = 1
        IB = 1
        DO WHILE (IA <= NA .AND. IB <= NB)
            DIFF = A(IA) - B(IB)
            IF (DIFF == 0) THEN ! A == B
                DUPS = DUPS + 1
                IA = IA + 1
                IB = IB + 1
            ELSE IF (DIFF < 0) THEN ! A < B
                IA = IA + 1
            ELSE ! A > B
                IB = IB + 1
            END IF
        END DO 
        ! one list exhausted, count will not inc
    END SUBROUTINE OVERLAP
    ! merge sorted arrays by Z2 fusion
    SUBROUTINE MERGE()
        ! local variable
        INTEGER :: IA, IB, IC, DIFF
    
        IF (ALLOCATED(C)) THEN
            IF (SIZE(C) < NA + NB) THEN
                DEALLOCATE(C)
                ALLOCATE(C(NA + NB))
            END IF
        ELSE
            ALLOCATE(C(NA + NB))
        END IF
        IA = 1
        IB = 1
        IC = 1
        DO WHILE (IA <= NA .AND. IB <= NB)
            DIFF = A(IA) - B(IB)
            IF (DIFF == 0) THEN ! A == B
                IA = IA + 1
                IB = IB + 1
            ELSE IF (DIFF < 0) THEN ! A < B
                C(IC) = A(IA)
                IA = IA + 1
                IC = IC + 1
            ELSE ! A > B
                C(IC) = B(IB)
                IB = IB + 1
                IC = IC + 1
            END IF
        END DO 
        ! one list exhausted, append the rest of the other
        IF (IA <= NA) THEN
            NC = IC + NA - IA
            C(IC:NC) = A(IA:NA)
        ELSE
            NC = IC + NB - IB
            C(IC:NC) = B(IB:NB)
        END IF
    END SUBROUTINE MERGE
    ! combine sorted arrays
    SUBROUTINE COMBINE()
        ! local variable
        INTEGER :: IA, IB, IC, DIFF
    
        IF (ALLOCATED(C)) THEN
            IF (SIZE(C) < NA + NB) THEN
                DEALLOCATE(C)
                ALLOCATE(C(NA + NB))
            END IF
        ELSE
            ALLOCATE(C(NA + NB))
        END IF
        IA = 1
        IB = 1
        IC = 1
        DO WHILE (IA <= NA .AND. IB <= NB)
            DIFF = A(IA) - B(IB)
            IF (DIFF == 0) THEN ! A == B
                C(IC) = A(IA)
                IA = IA + 1
                IB = IB + 1
                IC = IC + 1
            ELSE IF (DIFF < 0) THEN ! A < B
                C(IC) = A(IA)
                IA = IA + 1
                IC = IC + 1
            ELSE ! A > B
                C(IC) = B(IB)
                IB = IB + 1
                IC = IC + 1
            END IF
        END DO 
        ! one list exhausted, append the rest of the other
        IF (IA <= NA) THEN
            NC = IC + NA - IA
            C(IC:NC) = A(IA:NA)
        ELSE
            NC = IC + NB - IB
            C(IC:NC) = B(IB:NB)
        END IF
    END SUBROUTINE COMBINE
    ! intersect sorted arrays
    SUBROUTINE INTERSECT()
        ! local variable
        INTEGER :: IA, IB, IC, DIFF
        
        NC = MAX(NA, NB)
        IF (ALLOCATED(C)) THEN
            IF (SIZE(C) < NC) THEN
                DEALLOCATE(C)
                ALLOCATE(C(NC))
            END IF
        ELSE
            ALLOCATE(C(NC))
        END IF
        IA = 1
        IB = 1
        IC = 1
        DO WHILE (IA <= NA .AND. IB <= NB)
            DIFF = A(IA) - B(IB)
            IF (DIFF == 0) THEN ! A == B
                C(IC) = A(IA)
                IA = IA + 1
                IB = IB + 1
                IC = IC + 1
            ELSE IF (DIFF < 0) THEN ! A < B
                IA = IA + 1
            ELSE ! A > B
                IB = IB + 1
            END IF
        END DO 
        ! one list exhausted, discard the rest of the other
        NC = IC - 1
    END SUBROUTINE INTERSECT
END MODULE IPU

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