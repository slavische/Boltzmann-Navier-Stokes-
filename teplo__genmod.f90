        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 19 10:18:35 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TEPLO__genmod
          INTERFACE 
            SUBROUTINE TEPLO(KM,L1,D1,USL,USP,DTAU,TL,TP,N,Y,DXX)
              REAL(KIND=4) :: KM(5002)
              REAL(KIND=4) :: L1(5002)
              REAL(KIND=4) :: D1(5002)
              INTEGER(KIND=4) :: USL
              INTEGER(KIND=4) :: USP
              REAL(KIND=4) :: DTAU
              REAL(KIND=4) :: TL
              REAL(KIND=4) :: TP
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: Y(0:5002)
              REAL(KIND=4) :: DXX(5002)
            END SUBROUTINE TEPLO
          END INTERFACE 
        END MODULE TEPLO__genmod
