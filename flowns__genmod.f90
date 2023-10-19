        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 19 10:18:35 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FLOWNS__genmod
          INTERFACE 
            SUBROUTINE FLOWNS(RL_B,TL_B,UL_B,RHOA,TEMP,VEL,CVA,N)
              REAL(KIND=4) :: RL_B
              REAL(KIND=4) :: TL_B
              REAL(KIND=4) :: UL_B
              REAL(KIND=4) :: RHOA(5002)
              REAL(KIND=4) :: TEMP(5002)
              REAL(KIND=4) :: VEL(5002)
              REAL(KIND=4) :: CVA
              INTEGER(KIND=4) :: N
            END SUBROUTINE FLOWNS
          END INTERFACE 
        END MODULE FLOWNS__genmod
