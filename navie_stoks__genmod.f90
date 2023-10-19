        !COMPILER-GENERATED INTERFACE MODULE: Thu Oct 19 10:18:35 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NAVIE_STOKS__genmod
          INTERFACE 
            SUBROUTINE NAVIE_STOKS(RL_B,TL_B,UL_B,RL_NS,TL_NS,UL_NS,DTAU&
     &,RHOA,TEMP,CVA,VEL,XX,N,RR_B,TR_B,UR_B,RR_NS,TR_NS,UR_NS)
              REAL(KIND=4) :: RL_B
              REAL(KIND=4) :: TL_B
              REAL(KIND=4) :: UL_B
              REAL(KIND=4) :: RL_NS
              REAL(KIND=4) :: TL_NS
              REAL(KIND=4) :: UL_NS
              REAL(KIND=4) :: DTAU
              REAL(KIND=4) :: RHOA(5002)
              REAL(KIND=4) :: TEMP(5002)
              REAL(KIND=4) :: CVA
              REAL(KIND=4) :: VEL(5002)
              REAL(KIND=4) :: XX(5002)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: RR_B
              REAL(KIND=4) :: TR_B
              REAL(KIND=4) :: UR_B
              REAL(KIND=4) :: RR_NS
              REAL(KIND=4) :: TR_NS
              REAL(KIND=4) :: UR_NS
            END SUBROUTINE NAVIE_STOKS
          END INTERFACE 
        END MODULE NAVIE_STOKS__genmod
