        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 19 10:18:35 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SHASTA__genmod
          INTERFACE 
            SUBROUTINE SHASTA(A,B,L,DX,VEL,DTAU,N,BLEF,BRIG,WLEF,WRIG)
              REAL(KIND=4) :: A(5002)
              REAL(KIND=4) :: B(5002)
              REAL(KIND=4) :: L(5002)
              REAL(KIND=4) :: DX
              REAL(KIND=4) :: VEL(5002)
              REAL(KIND=4) :: DTAU
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: BLEF
              REAL(KIND=4) :: BRIG
              REAL(KIND=4) :: WLEF
              REAL(KIND=4) :: WRIG
            END SUBROUTINE SHASTA
          END INTERFACE 
        END MODULE SHASTA__genmod
