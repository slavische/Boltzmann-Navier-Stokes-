
!  One-component gas with solid parts. Integro_c 2004.
!  Temperature gradient in solid parts.

PROGRAM BOLTZ_NAVIE
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   parameter (rx=5002)
   real rhoa(rx),temp(rx),vel(rx),xx(rx)

   COMMON/BDISTFUNC/DISTFUNC1(M0,MFS,MX),DISTFUNC2(M0,MFS,MX)
   COMMON/BDIST/DIST(M0,MFS,MX)

   COMMON/BTIME/NTAU,KTI0,TIMAX/BDTAU/DTAU,DTR/BTI/TIME/BNVS/NVS
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX/BIKS/IKS/BSTIM/KTI/BNKC/NK
   COMMON/BDIAMETR/DIAM(M0,MX)/BTSLD/TSLD(MX)/BSORB/SORB(MX)
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)

   common/nav_st/rl_ns,tl_ns,ul_ns,rr_ns,tr_ns,ur_ns

   IKS=1

   MAXKTI=1000
   TIME=0.
   MT=0
   !
   !
   !        OPEN(2,FILE='r1.dat',STATUS='UNKNOWN',
   !     *  ACCESS='SEQUENTIAL',FORM='FORMATTED')
   !
   CALL FLOW
   !
   CALL ZEROIT
   CALL COLLISION


   DIST=DISTFUNC1
   NVS=1

   CALL MOMENT
   CALL FPRINT

   rl_b=dens(1,nk)
   tl_b=t(1,nk)
   ul_b=u0(1,nk)

   DIST=DISTFUNC2
   NVS=2

   CALL MOMENT
   CALL FPRINT

   rr_b=dens(1,2)
   tr_b=t(1,2)
   ur_b=u0(1,2)


   call flowns(rl_b,tl_b,ul_b,rhoa,temp,vel,cva,n)

   rl_ns=rhoa(1)
   tl_ns=temp(1)
   ul_ns=vel(1)

   rr_ns=rhoa(n)
   tr_ns=temp(n)
   ur_ns=vel(n)


   LKTI=1./DTR
   KTI=0
   KTII=LKTI


   !      PAUSE

   DO 1 NSTEP=1,MAXKTI

      CALL SYSTEM_CLOCK(COUNT=M1,COUNT_RATE=M2)
      XM1=M1
      D=(XM1-MT)/M2
      MTACT=M1-MT
      MT=M1
      PRINT *,' '
      WRITE(*,555) D,TIME,KTI,KTII
555   FORMAT(3X,'Real step time [sec] =',F10.5,3X, &
         'Time =',F9.3,2I10////)
      !--------------------------------------------------------------------

      DIST=DISTFUNC1
      NVS=1
      CALL MOMENT

      DO 21 IT=1,NTAU
         CALL FUNCAL
         CALL MOMENT
21    CONTINUE
      CALL REL
      CALL MOMENT
      CALL FPRINT


      rl_b=dens(1,nk)
      tl_b=t(1,nk)
      ul_b=u0(1,nk)


      DISTFUNC1=DIST
      DIST=DISTFUNC2
      CALL MOMENT
      NVS=2

      DO 22 IT=1,NTAU
         CALL FUNCAL
         CALL MOMENT
22    CONTINUE
      CALL REL
      CALL MOMENT
      CALL FPRINT

      DISTFUNC2=DIST
      TIME=TIME+DTR

      rr_b=dens(1,2)
      tr_b=t(1,2)
      ur_b=u0(1,2)

      !--------------------------------------------------------------------
      dtau=dtr
      call navie_stoks(rl_b,tl_b,ul_b,rl_ns,tl_ns,ul_ns, &
         dtau,rhoa,temp,cva,vel,xx,n,rr_b,tr_b,ur_b,rr_ns,tr_ns,ur_ns)

      rl_ns=rhoa(1)
      tl_ns=temp(1)
      ul_ns=vel(1)

      rr_ns=rhoa(n)
      tr_ns=temp(n)
      ur_ns=vel(n)

      call fprint_ns(xx,rhoa,temp,vel,n)

      !        PAUSE

      KTI=KTI+1
      IF(KTI.EQ.KTII) THEN
         OPEN(1,FILE='f1_',STATUS='UNKNOWN', &
            ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
         WRITE(1) DISTFUNC1,DISTFUNC2,TIME,DIAM,TSLD,SORB
         CLOSE(1)

         OPEN(15,FILE='init1_',STATUS='UNKNOWN', &
            ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
         WRITE(15) rhoa,temp,vel
         CLOSE(15)

         KTII=KTII+LKTI
         PRINT *,'Record of the Function'
         IF(TIME.GT.10000) STOP
      END IF
1  CONTINUE

   STOP
END


SUBROUTINE FLOW

   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BINTS/A1,A2,B1,B2,C1,C2/BNUM/NVWM,NVM,NUM,NWM,NMAX
   COMMON/BUVMAX0/UM0,VM0/BUVR0/UR0,VR0/BUV0/U00,V00/BR/R,R1
   COMMON/BDTAU/DTAU,DTR/BDSP/DU,DV,DW
   COMMON/BRNT/RNC,TR/BSNT/RS,TS/BSNT1/RS1,TS1
   COMMON/BJ/JY,JS,JD,JZ,JYZ/BTIME/NTAU,KTI0,TIMAX
   COMMON/BK/NK4,NK5,NK6,DDX,DDY,DDZ,XDEL1,YDEL1,ZDEL1
   COMMON/IR/IRAND,LARGE,MULT/BXX/X/BNKC/NK
   COMMON/BTW1/ZW1,TW1,VW1/BTW2/ZW2,TW2,VW2
   COMMON/KOROB/NA1,NA2,NA3,NA4,NA5,NP,XNP,BT/BMIT/MIT
   COMMON/BJSS/JSS
   COMMON/XDIAM/DX1,DX2
   COMMON/BPRM/PR1,PR2,PR3,PR4,PR5,NR6,PR7,NR8,NR81,PR9,NR10,PR21
   COMMON/BSOLID/TSOLID,FSOLID,VSLD,XSOLID,CL
   COMMON/BGAS/GAS,DGAS
   COMMON/BCORR/CORR

   CORR=1

   TRANSF=1.

   GAS=2.4E+25
   SLD=1.E+24
   DGAS=3.E-10

   ! N gas = 2.4e+25
   ! N solid = 1.e24
   ! V(sph-cyl)/Vsph=1+3*Lcyl/(2*Dsph)=8.5 for Lcyl=500, Dsph=100

   XSOLID=100./3.
   !       CL=500./3.7
   CL=0.
   !       FSOLID=0.0025259
   !       FSOLID=0.0001
   FSOLID=SLD/GAS
   TSOLID=1.05
!  FOR MANY-DIMENTIONAL PROBLEM:
   !       JY=MY-1
   !       JZ=MZ-1
   !       JYZ=JY*JZ
   !       JS=JY
   !       JSS=JZ
   !       JD=7
   !       NK=MX-1
!  FOR ONE-DIMENTIONAL PROBLEM:
   JY=1
   JZ=1
   JYZ=JY*JZ
   JS=JY
   JSS=JZ
   NK=21
   !
   XDEL1=0.5
   YDEL1=XDEL1
   ZDEL1=XDEL1
   DDX=1.00
   DDY=1.00
   DDZ=1.00
   !
   MIT=1
   !
   NA1=1
   !
   !       NA2=11
   !       NA3=4
   !       NA4=5
   !       NA5=3
   !       NP=13
   !       XNP=13.
   !  --------------
   !       NA2=14
   !       NA3=12
   !       NA4=7
   !       NA5=6
   !       NP=23
   !       XNP=23.
   !  --------------
   !       NA2=42
   !       NA3=15
   !       NA4=47
   !       NA5=13
   !       NP=53
   !       XNP=53.
   !  --------------
   !       NA2=24
   !       NA3=71
   !       NA4=88
   !       NA5=92
   !       NP=101
   !       XNP=101.
   ! ---------------
   NA2=42
   NA3=172
   NA4=60
   NA5=132
   NP=199
   XNP=199.
   !  --------------
   !        NA2=277
   !        NA3=286
   !        NA4=16
   !        NA5=134
   !        NP=307
   !        XNP=307.
   !  ---------------
   !       NA2=499
   !       NA3=53
   !       NA4=297
   !       NA5=194
   !       NP=523
   !       XNP=523.
   !-----------------------------
   NUM=12
   U00=0.
   NUP=6
   NUN=6
   DU=0.75
   A=NUM*DU
   A1=-DU*NUN
   A2=DU*NUP
   NVM=NUM
   V00=0.
   DV=DU
   N=NVM/2
   NVN=N
   NVP=N
   NWM=NUM
   NWN=N
   NWP=N
   NVWM=NVM*NWM
   B2=N*DV
   B1=-B2
   DW=DU
   C2=N*DW
   C1=-C2
   !
   R1=N*DU
   R=R1*R1
   BT=SQRT(2.)*(A2-A1)*(B2-B1)*(C2-C1)
   !       BT=(3.14159/SQRT(2.))*(A2-A1)*(B2-B1)*(C2-C1)
   !       BT1=4*XSOLID*SQRT(2.)*3.14159*R1**3/3.
   DTAU=0.8*XDEL1/A2
   !------------------------
!       DTAU=0.01
   !------------------------
   NTAU=1
   DTR=DTAU*NTAU
   TIMAX=100000.
   KTI0=1

   UM0=0.
   VM0=0.
   UR0=0.
   VR0=0.

   ! DATA FOR TIME=0

   RNC=0.8
   TR=1.0

   ! PARAMETERS OF THE INTERFACE

   RS=1.0
   TS=1.0

   RS1=0.8
   TS1=1.0

   ! PARAMETERS OF THE METAL SURFACE

   TW1=1.0
   ZW1=1.0

   TW2=1.0
   ZW2=1.0
   !
   TSOLID=TSOLID*TRANSF

   RND=0.333
   IRAND=137462873
   MULT=65539
   LARGE=2147483647
   X=1.

   RETURN
END


SUBROUTINE FPRINT
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX),&
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX),&
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   COMMON/BTLS/TLS(MX)/BNVS/NVS
   COMMON/BROST/ROSTS(MX)/BSORB/SORB(MX)/BTSLD/TSLD(MX)
   COMMON/BXY/XL(MX),YL(MY),ZL(MZ)/BTI/TIME
   COMMON/XDIAM/DX1,DX2
   COMMON/BK/NK4,NK5,NK6,DDX,DDY,DDZ,XDEL1,YDEL1,ZDEL1
   COMMON/BNKC/NK/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BDTAU/DTAU,DTR
   COMMON/BFFW/FW1(20),FW2(20),FW3(MFS),UFW(MFS), &
      VFW(MFS),WFW(MFS)/BNUM/NVWM,NVM,NUM,NWM,NMAX
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BDSP/DU,DV,DW
   COMMON/BKOEF/XKOEF/BSTIM/KTI
   COMMON/BDIAMETR/DIAM(M0,MX)/BGAS/GAS,DGAS
   COMMON/BTEST/UTEST(1000),VTEST(1000),WTEST(1000)
   CHARACTER H0,H1,H2,H3,H4,H5,H6,H7

   SMAS=0.
   FL=0.

   S2=SQRT(2.)
   B=1.
   LR=(NK-1)*XDEL1
   LN=1
   XKN=1./LR

   DO 3 J=1,JYZ
      DO 3 I=2,NK
         SMAS=SMAS+DENS(J,I)
3  CONTINUE

   DO 31 L=2,NK
      FLUX=XJP(1,L)-XJN(1,L)
      FL=FL+FLUX
31 CONTINUE

   SMAS=SMAS/(NK-1)
   FL=FL/(NK-1)/S2

   H0=''
   H1='X'
   H2='N'
   H3='T'
   H4='U'
   H5='S'
   H7='P'
   H8='  Q'
   H9='  JM'
   H10=' T*      '

   IF(NVS.EQ.1) THEN

      OPEN(2,FILE='d_Ln.dat',STATUS='UNKNOWN', &
         ACCESS='SEQUENTIAL',FORM='FORMATTED')
      DNX=0.

   ELSE

      OPEN(2,FILE='d_Rn.dat',STATUS='UNKNOWN', &
         ACCESS='SEQUENTIAL',FORM='FORMATTED')
      !        DNX=90.
      DNX=490.

   END IF

   WRITE(2,557) TIME,DTR,SMAS
   WRITE(2,560) H1,H2,H3,H4,H7,H8,H9

   !       WRITE(*,557) TIME,DTR,SMAS
   WRITE(*,560) H1,H2,H3,H4,H7,H8,H9

   DO 778 L=2,NK
      FLUX=XJP(1,L)-XJN(1,L)
      WRITE(2,550) (XL(L)+DNX)/LR,DENS(1,L),T(1,L), &
         U0(1,L),P(1,L),QX(1,L),FLUX
778 CONTINUE

   DO 777 L=2,NK,LN
      FLUX=XJP(1,L)-XJN(1,L)
      WRITE(*,550) (XL(L)+DNX)/LR,DENS(1,L),T(1,L), &
         U0(1,L),P(1,L),QX(1,L),FLUX
777 CONTINUE


550 FORMAT(2X,F7.4,6F10.5)
557 FORMAT(/3X,'time: ',F8.3,' dt:',F10.5,' mass:',F8.3)
560 FORMAT(3X,A4,6A10)

   CLOSE(2)

   !            WRITE(*,*) 'DENSITY',TIME
   !            WRITE(*,11) ((DENS(I,J),I=2,JY),J=2,NK)

   RETURN
END


SUBROUTINE ZEROIT
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDISTFUNC/DISTFUNC1(M0,MFS,MX),DISTFUNC2(M0,MFS,MX)
   COMMON/BDIST/DIST(M0,MFS,MX)
   COMMON/BR/R,R1/BINTS/A1,A2,B1,B2,C1,C2
   COMMON/BNKC/NK/BUV0/U00,V00/BTI/TIME
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BK/NK4,NK5,NK6,DDX,DDY,DDZ,XDEL1,YDEL1,ZDEL1
   COMMON/BXDEL/XDEL(MX),YDEL(MY),ZDEL(MZ)/BXY/XL(MX),YL(MY),ZL(MZ)
   COMMON/BZERO/DELT,DENS0,PXX0,PYY0,PZZ0/BDSP/DU,DV,DW
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BTW1/ZW1,TW1,VW1/BTW2/ZW2,TW2,VW2
   COMMON/BSNT/RS,TS/INDEX/IND(20,20,20)
   COMMON/BIKS/IKS
   COMMON/BTEMPL1/TL1(M0),VL1(M0)
   COMMON/BSOLID/TSOLID,FSOLID,VSLD,XSOLID,CL
   COMMON/BDIAMETR/DIAM(M0,MX)
   COMMON/BTLS/TLS(MX)
   COMMON/BROST/ROSTS(MX)/BSORB/SORB(MX)/BTSLD/TSLD(MX)

   NK1=NK+1
   JY1=JY
   SU(1)=A1+DU/2
   DO 1 I=2,NUM
      SU(I)=SU(I-1)+DU
1  CONTINUE
   SV(1)=B1+DV/2
   DO 2 J=2,NVM
      SV(J)=SV(J-1)+DV
2  CONTINUE
   SW(1)=C1+DW/2
   DO 3 K=2,NWM
      SW(K)=SW(K-1)+DW
3  CONTINUE
   !       WRITE(*,105) SU,SV,SW
105 FORMAT(1X,(8F8.3))

   DELT=DU*DV*DW

   DO 101 IU=1,20
      DO 101 IV=1,20
         DO 101 IW=1,20
            IND(IU,IV,IW)=1
101 CONTINUE

   NF=0
   VSLD=0.
   VW1=0.
   VW2=0.
   DO 4 IU=1,NUM
      U=SU(IU)
      DO 5 IV=1,NVM
         V=SV(IV)
         DO 5 IW=1,NWM
            W=SW(IW)
            NF=NF+1
            IND(IU,IV,IW)=NF
            TT=U*U+V*V+W*W
            FMR=FMAXR(U,V,W)
            !		   +FMAX(U,V,W)
            DO 6 L=1,NK1
               DO 6 J=1,JYZ
                  DISTFUNC1(J,NF,L)=FMR
                  DISTFUNC2(J,NF,L)=FMR
6           CONTINUE
            IF(W.GE.0) THEN
               VW1=VW1+W*EXP(-0.5/TW1*TT)/(6.28318*TW1)**(1.5)
               VW2=VW2+W*EXP(-0.5/TW2*TT)/(6.28318*TW2)**(1.5)
               VSLD=VSLD+W*EXP(-0.5/TSOLID*TT)/(6.28318*TSOLID)**(1.5)
            END IF
5     CONTINUE
4  CONTINUE

   NMAX=NF+1

   DO 7 L=1,NK1
      DO 7 J=1,JYZ
         DISTFUNC1(J,NMAX,L)=DISTFUNC1(J,NF,L)
         DISTFUNC2(J,NMAX,L)=DISTFUNC2(J,NF,L)
7  CONTINUE

   !
   DO 81 L=1,NK
      TLS(L)=1.
      ROSTS(L)=0.
      SORB(L)=0.
      TSLD(L)=1.
      DO 81 J=1,JYZ
         DIAM(J,L)=XSOLID
81 CONTINUE

   IF(IKS.EQ.2) THEN
      !******************************************************
      OPEN(1,FILE='f1_',STATUS='OLD', &
         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      READ(1) DISTFUNC1,DISTFUNC2,TIME,DIAM,TSLD,SORB
      CLOSE(1)
      !******************************************************
   END IF

   XDEL(1)=XDEL1
   XL(1)=-XDEL1/2

   DO 9 M=2,NK+1
      XDEL(M)=XDEL1
9  XL(M)=XL(M-1)+XDEL(M)

   !       YDEL(1)=YDEL1
   !       YL(1)=YDEL1/2

   !       DO 10 M=1,JY
   !         YDEL(M+1)=YDEL(M)*DDY
   !  10     YL(M+1)=YL(M)+0.5*(YDEL(M)+YDEL(M+1))

   !       ZDEL(1)=ZDEL1
   !       ZL(1)=0.5*ZDEL1
   !       DO 70 M=1,JZ
   !       ZDEL(M+1)=ZDEL(M)*DDZ
   !   70  ZL(M+1)=ZL(M)+0.5*(ZDEL(M)+ZDEL(M+1))
   !

   DENS0=1.
   PXX0=1.
   PYY0=1.
   PZZ0=1.
   VW1=VW1*DELT
   VW2=VW2*DELT
   VSLD=VSLD*DELT

   RETURN
END


SUBROUTINE MOMENT
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDIST/DIST(M0,MFS,MX)
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX
   COMMON/BUV0/U00,V00
   COMMON/BNKC/NK/BR/R,R1/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BZERO/DELT,DENS0,PXX0,PYY0,PZZ0
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   COMMON/BFFW/FW1(20),FW2(20),FW3(MFS),UFW(MFS), &
      VFW(MFS),WFW(MFS)
   COMMON/BKOEF/XKOEF
   DIMENSION A11(M0),A12(M0),A13(M0),A14(M0),A15(M0), &
      A22(M0),A33(M0),A44(M0),XN(M0),XP(M0),YP(M0),YN(M0),&
      ZP(M0),AQ(M0),AE(M0)

   NK1=NK+1
   NK2=NK/2+1
   MSS=MF0
   DO 1 N=1,20
      FW1(N)=0
1  FW2(N)=0

   DO 2 L=1,NK1
      DO 3 J=1,JYZ
         A11(J)=0.
         A12(J)=0.
         A13(J)=0.
         A14(J)=0.
         A15(J)=0.
         A22(J)=0.
         A33(J)=0.
         A44(J)=0.
         AQ(J)=0.
         AE(J)=0.
         XN(J)=0.
         XP(J)=0.
         ZP(J)=0.
         YP(J)=0.
         YN(J)=0.
3     CONTINUE
      NF=0
      KLF=0
      NUMF=NUM/2+1
      DO 4 IU=1,NUM
         U=SU(IU)
         DO 4 IV=1,NVM
            V=SV(IV)
            DO 4 IW=1,NWM
               W=SW(IW)
               NF=NF+1
               C=U*U+V*V+W*W
               IF(IV.EQ.NUMF) THEN
                  IF(IW.EQ.NUMF) THEN
                     KLF=KLF+1
                     UFW(KLF)=U
                     FW1(KLF)=DIST(1,NF,2)
                     FW2(KLF)=DIST(1,NF,2)
                  END IF
               END IF
               DO 5 J=1,JYZ
                  FIK=DIST(J,NF,L)
                  IF(U.LT.0) XN(J)=XN(J)-U*FIK
                  IF(U.GE.0) XP(J)=XP(J)+U*FIK
                  IF(V.LT.0) YN(J)=YN(J)-V*FIK
                  IF(V.GE.0) YP(J)=YP(J)+V*FIK
                  IF(W.GT.0) ZP(J)=ZP(J)+W*FIK
                  A11(J)=A11(J)+FIK
                  A12(J)=A12(J)+U*FIK
                  A13(J)=A13(J)+V*FIK
                  A14(J)=A14(J)+W*FIK
                  A15(J)=A15(J)+C*FIK
                  A22(J)=A22(J)+U*U*FIK
                  A33(J)=A33(J)+V*V*FIK
                  A44(J)=A44(J)+W*W*FIK
5              CONTINUE
4     CONTINUE
      DO 6 J=1,JYZ
         DE=A11(J)*DELT
         DENS(J,L)=DE
         XJN(J,L)=XN(J)*DELT
         XJP(J,L)=XP(J)*DELT
         YJN(J,L)=YN(J)*DELT
         YJP(J,L)=YP(J)*DELT
         ZJP(J,L)=ZP(J)*DELT
         XJX(J,L)=A12(J)*DELT
         E0(J,L)=A15(J)*DELT/DE
         U0(J,L)=A12(J)*DELT/DE
         V0(J,L)=A13(J)*DELT/DE
         W0(J,L)=A14(J)*DELT/DE
         PXX(J,L)=A22(J)*DELT-U0(J,L)*U0(J,L)*DE
         PYY(J,L)=A33(J)*DELT-V0(J,L)*V0(J,L)*DE
         PZZ(J,L)=A44(J)*DELT-W0(J,L)*W0(J,L)*DE
         P(J,L)=(PXX(J,L)+PYY(J,L)+PZZ(J,L))/3.
         TX(J,L)=PXX(J,L)/DE
         TY(J,L)=PYY(J,L)/DE
         TZ(J,L)=PZZ(J,L)/DE
         T(J,L)=P(J,L)/DE
6     CONTINUE
      NF=0
      KLF=0
      DO 7 IU=1,NUM
         U=SU(IU)
         DO 7 IV=1,NVM
            V=SV(IV)
            DO 7 IW=1,NWM
               W=SW(IW)
               NF=NF+1
               DO 8 J=1,JYZ
                  AQ(J)=AQ(J)+((U-U0(J,L))**2+(V-V0(J,L))**2+&
                     (W-W0(J,L))**2)*(U-U0(J,L))*DIST(J,NF,L)
                  AE(J)=AE(J)+(U*U+V*V+W*W)*U*DIST(J,NF,L)
8              CONTINUE
7     CONTINUE
      DO 9 J=1,JYZ
         QX(J,L)=ABS(AQ(J)*DELT/2.)
         EX(J,L)=AE(J)*DELT/2.
9     CONTINUE
2  CONTINUE

   RETURN
END

SUBROUTINE FUNCAL
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDIST/DIST(M0,MFS,MX)/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BSP/SU(20),SV(20),SW(20)/BUV0/U00,V00
   COMMON/BDTAU/DTAU,DTR/BR/R,R1/BFYZ/FY(MY,MZ),FZ(MY,MZ)
   COMMON/BXDEL/XDEL(MX),YDEL(MY),ZDEL(MZ)
   COMMON/BDSP/DU,DV,DW/BUVW/U,V,W,NF,KF,NFR/BNKC/NK
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX/BHF1/HF1(MY,MZ)/BNVS/NVS
   DIMENSION HF(MY,MZ)
   XD=XDEL(1)
   !      YD=YDEL(1)
   !      ZD=ZDEL(1)
   DO 2 L=2,NK
      LK=NK-L+2
      NF=0
      DO 3 IU=1,NUM
         U=SU(IU)
         UA=ABS(U)
         IF(U.LT.0) THEN
            MI=L
            MI1=L+1
         ELSE
            MI=LK
            MI1=LK-1
         END IF
         DO 31 IV=1,NVM
            V=SV(IV)
            !             VA=ABS(V)
            DO 31 IW=1,NWM
               W=SW(IW)
               !             WA=ABS(W)
               NF=NF+1
               DO 4 I=1,JZ
                  DO 4 K=1,JY
                     J=(I-1)*JY+K
                     HF(K,I)=DIST(J,NF,MI)
                     HF1(K,I)=DIST(J,NF,MI1)
4              CONTINUE
               !   1-DIMENSIONAL

               IF((L.EQ.NK).AND.(U.LT.0)) CALL B
               IF((LK.EQ.2).AND.(U.GE.0)) CALL A

               !   2-DIMENSIONAL
               !               IF(V.LT.0) THEN
               !                 CALL G4(MI)
               !                 DO 5 I=1,JZ
               !                 DO 5 K=1,JY-1
               !  5                FY(K,I)=HF(K+1,I,KF)
               !               ELSE
               !                 NVR=NWM*(NVM*IU-IV)+IW
               !                 KFR=NSN(NVR)-IS
               !                 DO 6 I=1,JZ
               !   6               FY(1,I)=HF(1,I,KFR)
               !                 CALL G3(MI)
               !                 DO 6 I=1,JY
               !                 DO 6 K=2,JY
               !  6                FY(K,I)=HF(K-1,I,KF)
               !               END IF
               !  3-DIMENSIONAL
               !                IF(W.LT.0) THEN
               !                  CALL G9(MI)
               !                  DO 25 I=1,JZ-1
               !                  DO 25 K=1,JY
               !  25              FZ(K,I)=HF(K,I+1,KF)
               !              ELSE
               !                  NWR=NWM*(NVM*(IU-1)+IV)-IW+1
               !                  KFR=NSN(NWR)-IS
               !                  DO 30 K=1,JY
               !                  FZ(K,1)=HF(K,1,KFR)
               !                  DO 30 I=2,JZ
               !  30              FZ(K,I)=HF(K,I-1,KF)
               !              END IF
               DO 7 I=1,JZ
                  DO 7 K=1,JY
                     J=(I-1)*JY+K
                     FIK=HF(K,I)
                     !FOR 1-DIM.
                     DIST(J,NF,MI)=FIK-DTAU*UA*(FIK-HF1(K,I))/XD
                     !     &                +VA*(FIK-FY(K,I))/YD)
7              CONTINUE
31       CONTINUE
3     CONTINUE
2  CONTINUE
   RETURN
END


SUBROUTINE A
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BUVW/U,V,W,NF,KF,NFR/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BHF1/HF1(MY,MZ)/BNVS/NVS
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   COMMON/BTW1/ZW1,TW1,VW1/BSNT/RS,TS
   common/nav_st/rl_ns,tl_ns,ul_ns,rr_ns,tr_ns,ur_ns

   !       DENS_IN=2.
   !       FM=DENS_IN*FMAX(U,V,W)

   IF(NVS.EQ.1) THEN

      !        FM=EXP(-0.5/TW1*(U*U+V*V+W*W))/(6.28318*TW1)**(3./2)
      FM=EXP(-0.5/TS*(U*U+V*V+W*W))/(6.28318*TS)**(1.5)*RS
      DO 1 I=1,JZ
         DO 1 K=1,JY
            !        J=(I-1)*JY+K
            !        RN=XJN(1,2)/VW1
            !        HF1(K,I)=FM*RN
            HF1(K,I)=FM
1     CONTINUE

   ELSE

      FM=rr_ns*EXP(-0.5/tr_ns*((U-ur_ns)**2+V*V+W*W))/(6.28318*tr_ns)**(3./2)
      DO 2 I=1,JZ
         DO 2 K=1,JY
            HF1(K,I)=FM
2     CONTINUE

   END IF

   RETURN
END


SUBROUTINE B
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BHF1/HF1(MY,MZ)/BSNT/RS,TS/BNVS/NVS/BSNT1/RS1,TS1
   COMMON/BNKC/NK/BUVW/U,V,W,NF,KF,NFR/BTW2/ZW2,TW2,VW2
   COMMON/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   common/nav_st/rl_ns,tl_ns,ul_ns,rr_ns,tr_ns,ur_ns

   !       FM=FMAXR(U,V,W)
   !       FM=FMAX(U,V,W)

   IF(NVS.EQ.1) THEN

      FM=rl_ns*EXP(-0.5/tl_ns*((U-ul_ns)**2+V*V+W*W))/(6.28318*tl_ns)**(3./2)
      DO 1 I=1,JZ
         DO 1 K=1,JY
            HF1(K,I)=FM
1     CONTINUE

   ELSE

      !       FM=EXP(-0.5/TW2*(U*U+V*V+W*W))/(6.28318*TW2)**(3./2)
      FM=EXP(-0.5/TS1*(U*U+V*V+W*W))/(6.28318*TS1)**(1.5)*RS1
      DO 2 I=1,JZ
         DO 2 K=1,JY
            J=(I-1)*JY+K
            !       RN=XJP(J,NK)/VW2
            HF1(K,I)=FM
            !       HF1(K,I)=0.00000000001

2     CONTINUE

   END IF

   RETURN
END


FUNCTION FMAX(U,V,W)
   COMMON/BUVMAX0/UM0,VM0
   FMAX=EXP(-0.5*((U-UM0)**2+(V-VM0)**2+W*W))/(6.28318)**(3./2)
   RETURN
END


FUNCTION FMAXR(U,V,W)
   COMMON/BRNT/RNC,TR
   FMAXR=EXP(-0.5/TR*(U*U+V*V+W*W))/(6.28318*TR)**(3./2)*RNC
   RETURN
END


FUNCTION FMSOLID(U,V,W)
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BSOLID/TSOLID,FSOLID,VSOLID,XSOLID,CL
   FMSOLID=EXP(-0.5/TSOLID*(U*U+V*V+W*W))/(6.28318*TSOLID)**(3./2)
   RETURN
END


SUBROUTINE COLLISION
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BXKOR/XKOR(4,1525),IUKOR(1525),IVKOR(1525),IWKOR(1525)
   COMMON/BINTS/A1,A2,B1,B2,C1,C2
   COMMON/BR/R,R1/BUV0/U00,V00/BDSP/DU,DV,DW
   COMMON/KOROB/NA1,NA2,NA3,NA4,NA5,NP,XNP,BT
   DIMENSION IUKO(1525),IVKO(1525),IWKO(1525)

   R2=2*(R1-DU/3)
   R21=2*R1

   OPEN(5,FILE='dk.dat',STATUS='UNKNOWN', &
      ACCESS='SEQUENTIAL',FORM='FORMATTED')

   DO 100 I=1,NP

      XX1=NA1*I/XNP
      XX1=XX1-INT(XX1)
      XKOR(1,I)=(XX1-0.5)*R2+U00
      XX2=NA2*I/XNP
      XX2=XX2-INT(XX2)
      XKOR(2,I)=(XX2-0.5)*R2
      XX3=NA3*I/XNP
      XX3=XX3-INT(XX3)
      XKOR(3,I)=(XX3-0.5)*R2

      U1=XKOR(1,I)
      V1=XKOR(2,I)
      W1=XKOR(3,I)

      IU1=INT((U1-A1)/DU+1)
      IV1=INT((V1-B1)/DV+1)
      IW1=INT((W1-C1)/DW+1)

      XKOR(1,I)=SU(IU1)
      XKOR(2,I)=SV(IV1)
      XKOR(3,I)=SW(IW1)
      IUKOR(I)=IU1
      IVKOR(I)=IV1
      IWKOR(I)=IW1

      WRITE(5,5) I,XKOR(1,I),XKOR(2,I),XKOR(3,I), IUKOR(I),IVKOR(I),IWKOR(I)

100 CONTINUE

   M=NP

   DO 101 I=1,NP

      XX1=NA1*I/XNP
      XX1=XX1-INT(XX1)
      XK1=(XX1-0.5)*R21+U00
      !         XK1=A2-XK1+A1
      XX2=NA2*I/XNP
      XX2=XX2-INT(XX2)
      XK2=(XX2-0.5)*R21
      !         XK2=B2-XK2+B1
      XX3=NA3*I/XNP
      XX3=XX3-INT(XX3)
      XK3=(XX3-0.5)*R21
      !         XK3=C2-XK3+C1

      U1=XK1
      V1=XK2
      W1=XK3

      IU1=INT((U1-A1)/DU+1)
      IV1=INT((V1-B1)/DV+1)
      IW1=INT((W1-C1)/DW+1)

      DO 120 J=1,NP

         M1=IUKOR(J)
         M2=IVKOR(J)
         M3=IWKOR(J)

         IF((IU1.EQ.M1).AND.(IV1.EQ.M2).AND.(IW1.EQ.M3)) GO TO 101

120   CONTINUE

      M=M+1
      XKOR(1,M)=SU(IU1)
      XKOR(2,M)=SV(IV1)
      XKOR(3,M)=SW(IW1)
      IUKOR(M)=IU1
      IVKOR(M)=IV1
      IWKOR(M)=IW1

      WRITE(5,5) M,XKOR(1,M),XKOR(2,M),XKOR(3,M), IUKOR(M),IVKOR(M),IWKOR(M)

101 CONTINUE

   !        NP=M

5  FORMAT(1X,I7,3F10.5,3I7)
7  FORMAT(1X,7I7)

   CLOSE(5)
   RETURN
END


SUBROUTINE TRANSFER
   COMMON/BUVWS/US(24),VS(24),WS(24),US1(24),VS1(24),WS1(24)
   COMMON/BUVW/U,V,W,NF,KF,NFR/BUVW1/U1,V1,W1

   DX=U1-U
   DY=V1-V
   DZ=W1-W
   !----------------------
   US(1)=U
   VS(1)=V
   WS(1)=W
   US1(1)=U1
   VS1(1)=V1
   WS1(1)=W1

   US(2)=U+DX
   VS(2)=V
   WS(2)=W
   US1(2)=U1-DX
   VS1(2)=V1
   WS1(2)=W1

   US(3)=U
   VS(3)=V+DY
   WS(3)=W
   US1(3)=U1
   VS1(3)=V1-DY
   WS1(3)=W1

   US(4)=U
   VS(4)=V
   WS(4)=W+DZ
   US1(4)=U1
   VS1(4)=V1
   WS1(4)=W1-DZ
   !----------------------
   !         GO TO 10

   XC=DX/2+U
   YC=DY/2+V
   ZC=DZ/2+W
   CPZ=XC+YC
   CMZ=-XC+YC
   CPY=XC+ZC
   CMY=-XC+ZC
   CPX=YC+ZC
   CMX=-YC+ZC

   US(5)=CPZ-VS(1)
   VS(5)=CMZ+US(1)
   WS(5)=WS(1)
   US1(5)=CPZ-VS1(1)
   VS1(5)=CMZ+US1(1)
   WS1(5)=WS1(1)

   US(6)=CPZ-VS(2)
   VS(6)=CMZ+US(2)
   WS(6)=WS(2)
   US1(6)=CPZ-VS1(2)
   VS1(6)=CMZ+US1(2)
   WS1(6)=WS1(2)

   US(7)=CPZ-VS(3)
   VS(7)=CMZ+US(3)
   WS(7)=WS(3)
   US1(7)=CPZ-VS1(3)
   VS1(7)=CMZ+US1(3)
   WS1(7)=WS1(3)

   US(8)=CPZ-VS(4)
   VS(8)=CMZ+US(4)
   WS(8)=WS(4)
   US1(8)=CPZ-VS1(4)
   VS1(8)=CMZ+US1(4)
   WS1(8)=WS1(4)

   US(9)=CPY-WS(1)
   VS(9)=VS(1)
   WS(9)=CMY+US(1)
   US1(9)=CPY-WS1(1)
   VS1(9)=VS1(1)
   WS1(9)=CMY+US1(1)

   US(10)=CPY-WS(2)
   VS(10)=VS(2)
   WS(10)=CMY+US(2)
   US1(10)=CPY-WS1(2)
   VS1(10)=VS1(2)
   WS1(10)=CMY+US1(2)

   US(11)=CPY-WS(3)
   VS(11)=VS(3)
   WS(11)=CMY+US(3)
   US1(11)=CPY-WS1(3)
   VS1(11)=VS1(3)
   WS1(11)=CMY+US1(3)

   US(12)=CPY-WS(4)
   VS(12)=VS(4)
   WS(12)=CMY+US(4)
   US1(12)=CPY-WS1(4)
   VS1(12)=VS1(4)
   WS1(12)=CMY+US1(4)

   US(13)=US(1)
   VS(13)=CPX-WS(1)
   WS(13)=CMX+VS(1)
   US1(13)=US1(1)
   VS1(13)=CPX-WS1(1)
   WS1(13)=CMX+VS1(1)

   US(14)=US(2)
   VS(14)=CPX-WS(2)
   WS(14)=CMX+VS(2)
   US1(14)=US1(2)
   VS1(14)=CPX-WS1(2)
   WS1(14)=CMX+VS1(2)

   US(15)=US(3)
   VS(15)=CPX-WS(3)
   WS(15)=CMX+VS(3)
   US1(15)=US1(3)
   VS1(15)=CPX-WS1(3)
   WS1(15)=CMX+VS1(3)

   US(16)=US(4)
   VS(16)=CPX-WS(4)
   WS(16)=CMX+VS(4)
   US1(16)=US1(4)
   VS1(16)=CPX-WS1(4)
   WS1(16)=CMX+VS1(4)

   GO TO 10

   US(17)=US(5)
   VS(17)=CPX-WS(5)
   WS(17)=CMX+VS(5)
   US1(17)=US1(5)
   VS1(17)=CPX-WS1(5)
   WS1(17)=CMX+VS1(5)

   US(18)=US(6)
   VS(18)=CPX-WS(6)
   WS(18)=CMX+VS(6)
   US1(18)=US1(6)
   VS1(18)=CPX-WS1(6)
   WS1(18)=CMX+VS1(6)

   US(19)=US(7)
   VS(19)=CPX-WS(7)
   WS(19)=CMX+VS(7)
   US1(19)=US1(7)
   VS1(19)=CPX-WS1(7)
   WS1(19)=CMX+VS1(7)

   US(20)=US(8)
   VS(20)=CPX-WS(8)
   WS(20)=CMX+VS(8)
   US1(20)=US1(8)
   VS1(20)=CPX-WS1(8)
   WS1(20)=CMX+VS1(8)

   US(21)=US(9)
   VS(21)=CPX-WS(9)
   WS(21)=CMX+VS(9)
   US1(21)=US1(9)
   VS1(21)=CPX-WS1(9)
   WS1(21)=CMX+VS1(9)

   US(22)=US(10)
   VS(22)=CPX-WS(10)
   WS(22)=CMX+VS(10)
   US1(22)=US1(10)
   VS1(22)=CPX-WS1(10)
   WS1(22)=CMX+VS1(10)

   US(23)=US(11)
   VS(23)=CPX-WS(11)
   WS(23)=CMX+VS(11)
   US1(23)=US1(11)
   VS1(23)=CPX-WS1(11)
   WS1(23)=CMX+VS1(11)

   US(24)=US(12)
   VS(24)=CPX-WS(12)
   WS(24)=CMX+VS(12)
   US1(24)=US1(12)
   VS1(24)=CPX-WS1(12)
   WS1(24)=CMX+VS1(12)

10 CONTINUE

   RETURN
END


SUBROUTINE INTEGR_C
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDIST/DIST(M0,MFS,MX)
   COMMON/BENU/EN(M0,MX),ENU(M0,MX)
   COMMON/BDFMAX/DFMAX(M0,500,MX)
   COMMON/BNKC/NK/BIUIVIW/IU,IV,IW
   COMMON/BUVW/U,V,W,NF,KF,NFR/BUUVVWW/UU,VV,WW
   COMMON/BINTS/A1,A2,B1,B2,C1,C2/BDSP/DU,DV,DW
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX
   COMMON/BSP/SU(20),SV(20),SW(20)/INDEX/IND(20,20,20)
   COMMON/BUVW1/U1,V1,W1
   COMMON/BUVWS/US(24),VS(24),WS(24),US1(24),VS1(24),WS1(24)
   COMMON/BR/R,R1/BUV0/U00,V00
   COMMON/KOROB/NA1,NA2,NA3,NA4,NA5,NP,XNP,BT
   COMMON/BXX/X/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BXKOR/XKOR(4,1525),IUKOR(1525),IVKOR(1525),IWKOR(1525)
   COMMON/BTEST/UTEST(1000),VTEST(1000),WTEST(1000)
   DIMENSION TFN(M0,MX),TFNU(M0,MX)
   !       DIMENSION FUNCS(M0,MX),FUNCS1(M0,MX)

   R2=R1-DU/2
   !       R2=R1
   NCOLL=16
   EPS=1E-5

   J=1

   DO 1 L=2,NK
      !       DO 1 J=1,JYZ
      EN(J,L)=0.
      ENU(J,L)=0.
      TFNU(J,L)=0
      TFN(J,L)=0
1  CONTINUE
   !=====================================
   IP=0
   DO 2 I=1,NP

      U1=XKOR(1,I)
      V1=XKOR(2,I)
      W1=XKOR(3,I)
      IU1=IUKOR(I)
      IV1=IVKOR(I)
      IW1=IWKOR(I)
      !----------------------

      CALL TRANSFER

      !----------------------
      NCL=0

      DO 5 N=1,NCOLL

         AUS=ABS(US(N))
         AVS=ABS(VS(N))
         AWS=ABS(WS(N))
         AUS1=ABS(US1(N))
         AVS1=ABS(VS1(N))
         AWS1=ABS(WS1(N))

         !           K=0
         IF((AUS.GT.R2).OR.(AVS.GT.R2).OR.(AWS.GT.R2).OR. &
            (AUS1.GT.R2).OR.(AVS1.GT.R2).OR.(AWS1.GT.R2)) THEN
            !           NFS=NMAX
            !           NFS1=NMAX
            !           NCL=NCL+1
            !           GO TO 55
            GO TO 5
         END IF

         IUS=INT((US(N)-A1)/DU+1)
         IVS=INT((VS(N)-B1)/DV+1)
         IWS=INT((WS(N)-C1)/DW+1)

         X1=IUS-(US(N)+DU/2-A1)/DU
         IF(ABS(X1).LE.EPS) THEN
            X2=IVS-(VS(N)+DV/2-B1)/DV
            IF(ABS(X2).LE.EPS) THEN
               X3=IWS-(WS(N)+DW/2-C1)/DW
               IF(ABS(X3).LE.EPS) THEN

                  IUS1=INT((US1(N)-A1)/DU+1)
                  IVS1=INT((VS1(N)-B1)/DV+1)
                  IWS1=INT((WS1(N)-C1)/DW+1)

                  !	     IF(K.EQ.0) THEN

                  Q=SQRT(ABS((U1-U)*(US(N)-U)+(V1-V)*(VS(N)-V)+(W1-W)*(WS(N)-W)))
                  !               IF(Q.EQ.0.) GO TO 5

                  NCL=NCL+1

                  NF1=NWM*(NVM*(IU1-1)+IV1-1)+IW1
                  NFS=NWM*(NVM*(IUS-1)+IVS-1)+IWS
                  NFS1=NWM*(NVM*(IUS1-1)+IVS1-1)+IWS1

55                CONTINUE

                  DO 71 L=2,NK
                     !           DO 71 J=1,JYZ

                     TFN(J,L)=TFN(J,L)+Q*DIST(J,NFS,L)*DIST(J,NFS1,L)
                     TFNU(J,L)=TFNU(J,L)+Q*DIST(J,NF1,L)

71                CONTINUE

                  !           END IF

               END IF
            END IF
         END IF

5     CONTINUE

      IF(NCL.EQ.0) GO TO 2

      DO 8 L=2,NK
         !           DO 8 J=1,JYZ
         TFNU(J,L)=TFNU(J,L)/NCL
         TFN(J,L)=TFN(J,L)/NCL
8     CONTINUE

      DO 9 L=2,NK
         !           DO 9 J=1,JYZ
         ENU(J,L)=ENU(J,L)+TFNU(J,L)
         EN(J,L)=EN(J,L)+TFN(J,L)
9     CONTINUE

      IP=IP+1

2  CONTINUE

!       CLOSE(7)

   VOLL=BT/IP

   DO 10 L=2,NK
      !       DO 10 J=1,JYZ
      EN(J,L)=EN(J,L)*VOLL
      ENU(J,L)=ENU(J,L)*VOLL
10 CONTINUE

   RETURN
END


SUBROUTINE REL
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDIST/DIST(M0,MFS,MX)
   COMMON/BENU/EN(M0,MX),ENU(M0,MX)
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BDSP/DU,DV,DW/BDTAU/DTAU,DTR
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX/BIUIVIW/IU,IV,IW
   COMMON/BNKC/NK/INDEX/IND(20,20,20)
   COMMON/BJ/JY,JS,JD,JZ,JYZ/BR/R,R1/BUV0/U00,V00
   COMMON/BUVW/U,V,W,NF,KF,NFR/BUUVVWW/UU,VV,WW/BSNNU/SN,SNU
   COMMON/BCORR/CORR/BXX/X
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   DIMENSION DFF(M0,MX),DIF(M0,MX),SDNS(M0,MX),DIST_C(M0,MFS,MX)

   SN=0.
   SNU=0.
   DELT=4.*DU*DV*DW
   !         DELT=2*DU*DV*DW
   !         DELT=DU*DV*DW
   NK2=NK/2
   NV2=NVM/2+1
   NW2=NWM/2+1

   !         CALL VELMAX

   DO 2 L=2,NK
      DO 2 J=1,JYZ
         SDNS(J,L)=0.
2  CONTINUE

   DO 10 L=1,NK
      DO 10 N=1,NMAX
         DO 10 J=1,JYZ
            DIST_C(J,N,L)=DIST(J,N,L)
10 CONTINUE

   NF=0.

   !           RND=RNDM(X)
   !           KR=1
   !           IF(RND.LE.(0.5)) KR=4
   !           X=X+1

   DO 3 IU=1,NUM
      U=SU(IU)
      DO 3 IV=NV2,NVM
         !         DO 3 IV=1,NVM
         V=SV(IV)
         DO 3 IW=NW2,NWM
            !         DO 3 IW=1,NWM
            W=SW(IW)
            NF=IND(IU,IV,IW)

            DO 7 L=2,NK
               DO 7 J=1,JYZ
                  DIF(J,L)=DIST(J,NF,L)
                  DFF(J,L)=0.
7           CONTINUE

            CALL INTEGR_C

            DO 4 L=2,NK
               DO 4 J=1,JYZ

                  S=(DTR*EN(J,L)+DIF(J,L))/(1+DTR*ENU(J,L))
                  !             S=DIF(J,L)+DTR*(EN(J,L)-DIF(J,L)*ENU(J,L))
                  DFF(J,L)=S
                  SDNS(J,L)=SDNS(J,L)+S*DELT
4           CONTINUE
            IVX=NVM-IV+1
            IWX=NWM-IW+1
            NR1=IND(IU,IVX,IW)
            NR2=IND(IU,IV,IWX)
            NR3=IND(IU,IVX,IWX)
            DO 11 L=2,NK
               DO 11 J=1,JYZ
                  FJL=DFF(J,L)
                  ! 1-dim
                  DIST_C(J,NF,L)=FJL
                  DIST_C(J,NR1,L)=FJL
                  DIST_C(J,NR2,L)=FJL
                  DIST_C(J,NR3,L)=FJL
                  ! 2-dim
                  !             DIST_C(J,NF,L)=FJL
                  !             DIST_C(J,NR2,L)=FJL
                  ! 3-dim
                  !             DIST_C(J,NF,L)=FJL

11          CONTINUE

3  CONTINUE

   !       IMAX=NUM*NVM*NWM
   !       SN=SN/IMAX
   !       SNU=SNU/IMAX


   DO 12 L=2,NK
      DO 12 J=1,JYZ
         FCORR=DENS(J,L)/SDNS(J,L)
         IF(CORR.EQ.0) FCORR=1.
         DO 12 N=1,NMAX
            DIST(J,N,L)=DIST_C(J,N,L)*FCORR
            DIST(J,N,L)=DIST_C(J,N,L)
12 CONTINUE

   RETURN
END


FUNCTION RNDM(X)
   COMMON/IR/IRAND,LARGE,MULT
   RMAX=1.
   IRAND=IRAND*MULT
   IF(IRAND.LT.0) IRAND=(IRAND+LARGE)+1
   RNDM=RMAX*(IRAND-1)/(LARGE-1)
   RETURN
END


SUBROUTINE VELMAX
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDFMAX/DFMAX(M0,500,MX)
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BDSP/DU,DV,DW
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX/MASS/DM,SDM
   COMMON/BNKC/NK/BJ/JY,JS,JD,JZ,JYZ
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX),&
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX),&
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)

   NV=500
   DO 2 L=1,NK
      DO 2 N=1,NV
         DO 2 J=1,JYZ
            DFMAX(J,N,L)=1.E-10
2  CONTINUE

   DPU=NUM*DU
   DPV=NVM*DV
   DPW=NWM*DW
   DNU=DPU/2
   DNV=DPV/2
   DNW=DPW/2
   NUM2=NUM/2
   NVM2=NVM/2
   NWM2=NWM/2

   DO 3 IU=1,NUM2
      UP=SU(IU)+DPU
      UN=SU(IU)-DNU
      DO 3 IV=1,NVM
         V=SV(IV)
         DO 3 IW=1,NWM
            W=SW(IW)
            DO 7 L=2,NK
               DO 7 J=1,JYZ
                  DF=DENS(J,L)
                  TF=T(J,L)
                  UVWP=(UP-U0(J,L))**2+(V-V0(J,L))**2+(W-W0(J,L))**2
                  NP=INT(SQRT(UVWP)+0.5)
                  DFMAX(J,NP,L)=DF*EXP(-0.5/TF*UVWP)/(6.28318*TF)**(1.5)
                  UVWN=(UN-U0(J,L))**2+(V-V0(J,L))**2+(W-W0(J,L))**2
                  NN=INT(SQRT(UVWP)+0.5)
                  DFMAX(J,NN,L)=DF*EXP(-0.5/TF*UVWN)/(6.28318*TF)**(1.5)
7           CONTINUE
3  CONTINUE
   DO 4 IU=1,NUM
      U=SU(IU)
      DO 4 IV=1,NVM2
         VP=SV(IV)+DPV
         VN=SV(IV)-DNV
         DO 4 IW=1,NWM
            W=SW(IW)
            DO 8 L=2,NK
               DO 8 J=1,JYZ
                  DF=DENS(J,L)
                  TF=T(J,L)
                  UVWP=(U-U0(J,L))**2+(VP-V0(J,L))**2+(W-W0(J,L))**2
                  NP=INT(SQRT(UVWP)+0.5)
                  DFMAX(J,NP,L)=DF*EXP(-0.5/TF*UVWP)/(6.28318*TF)**(1.5)
                  UVWN=(U-U0(J,L))**2+(VN-V0(J,L))**2+(W-W0(J,L))**2
                  NN=INT(SQRT(UVWN)+0.5)
                  DFMAX(J,NN,L)=DF*EXP(-0.5/TF*UVWN)/(6.28318*TF)**(1.5)
8           CONTINUE
4  CONTINUE
   DO 5 IU=1,NUM
      U=SU(IU)
      DO 5 IV=1,NVM2
         V=SV(IV)
         DO 5 IW=1,NWM
            WP=SW(IW)+DPW
            WN=SW(IW)-DNW
            DO 9 L=2,NK
               DO 9 J=1,JYZ
                  DF=DENS(J,L)
                  TF=T(J,L)
                  UVWP=(U-U0(J,L))**2+(V-V0(J,L))**2+(WP-W0(J,L))**2
                  NP=INT(SQRT(UVWP)+0.5)
                  DFMAX(J,NP,L)=DF*EXP(-0.5/TF*UVWP)/(6.28318*TF)**(1.5)
                  UVWN=(U-U0(J,L))**2+(V-V0(J,L))**2+(WN-W0(J,L))**2
                  NN=INT(SQRT(UVWN)+0.5)
                  DFMAX(J,NN,L)=DF*EXP(-0.5/TF*UVWN)/(6.28318*TF)**(1.5)
9           CONTINUE
5  CONTINUE
   RETURN
END


SUBROUTINE SOLID_t
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BDIST/DIST(M0,MFS,MX)
   COMMON/BDTAU/DTAU,DTR
   COMMON/BSP/SU(20),SV(20),SW(20)
   COMMON/BNUM/NVWM,NVM,NUM,NWM,NMAX
   COMMON/BNKC/NK/BSOLID/TSOLID,FSOLID,VSLD,XSOLID,CL
   COMMON/BJ/JY,JS,JD,JZ,JYZ/BR/R,R1/BUV0/U00,V00
   COMMON/BUVW/U,V,W,NF,KF,NFR/BLOK1/BT1
   COMMON/BZERO/DELT,DENS0,PXX0,PYY0,PZZ0
   COMMON/BDIAMETR/DIAM(M0,MX)/BGAS/GAS,DGAS
   COMMON/BXY/XL(MX),YL(MY),ZL(MZ)
   COMMON/BK/NK4,NK5,NK6,DDX,DDY,DDZ,XDEL1,YDEL1,ZDEL1
   COMMON/BTW1/ZW1,TW1,VW1/BTW2/ZW2,TW2,VW2
   COMMON/BTLS/TLS(MX)
   COMMON/BROST/ROSTS(MX)/BSORB/SORB(MX)/BTSLD/TSLD(MX)
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)
   DIMENSION DIFN(MFS),DIFR(MFS)
   
   DMIN=13.51
   !-----------------------
   GAS=2.4E+25
   SLD=1.E+24
   DGAS=3.E-10
   BETA=1.
   QSORB=1.45E+7
   Q1M=14.5*3.35E-21
   TEPLOEMK=994.
   RO=2.1E+3
   RMS=0.5E-8
   SURF=4.*3.14159*(RMS**2)
   SLDV=4.*3.14159*(RMS**3)/3.
   SLDM=SLDV*RO
   TEPLSLD=TEPLOEMK*SLDM
   RH=8.31/0.002
   FL=1.03E-7
   T00=300.
   VBASE=SQRT(RH*T00)
   DS=0.00001
   DS_=DS*GAS
   SQ2=SQRT(2.)
   DT_=DTR*FL/VBASE
   CRSLD=DS_*SQ2*VBASE*SURF*DT_
   !       T_PREDEL=1.1
   !-----------------------

   !      MAX DS=0.065

   !-----------------------
   VVV=(4/3.)*3.14159*(DGAS**3)/8.
   UC=0.
   !-----------------------
   !       DO 1 L=2,NK
   !       IF(ABS(SORB(L-1)-SORB(L)).LE.(0.05)) GO TO 2
   !   1   CONTINUE
   !   2   NSK=L-1
   !-----------------------
   NSK=NK
   !       PRINT *,"L=",NSK
   !
   XDIFR1=FSOLID*DTR
   XDIFR2=FSOLID*CL*DTR/2
   DIAM1=DIAM(1,2)
   XLX=XL(NK)+XDEL1/2.

   DO 33 J=1,JYZ
      DO 33 L=2,NSK

         SUMROST1=0.
         SUMROST2=0.
         UU=U0(J,L)
         VV=V0(J,L)
         WW=W0(J,L)

         DO 34 N=1,NMAX
            DIFN(N)=0.
            DIFR(N)=0.
34       CONTINUE

         NF=0

         SN1=0.
         SN2=0.

         DM1=DIAM(J,L)
         DM2=DIAM(J,L)**2

         !               TSL=TSOLID
         !               TSL=TW1-XL(L)/XLX*(TW1-TW2)
         !               TLS(L)=TSL

         TSL=TSLD(L)
         CONDENSAT=BETA
         CALL TPREDEL(L,TSL,T_PREDEL,CONDENSAT)

         XFUNC1=XDIFR1*DM2*SQRT(2.*T(J,L))
         !               PRINT *,XFUNC1

         CRSLD1=CRSLD*SQRT(TSL)

         CRS=SORB(L)-CRSLD1

         !       PRINT *,CRS,CRSLD,SORB(L)

         IF(CRS.LE.0.) THEN
            DNS=DS*SORB(L)/CRSLD1
         ELSE
            DNS=DS
         END IF

         DO 3 IU=1,NUM
            U=SU(IU)
            DO 3 IV=1,NVM
               V=SV(IV)
               DO 3 IW=1,NWM
                  W=SW(IW)
                  NF=NF+1
                  FSURF1=EXP(-0.5/TSL*((U-UC)**2+V*V+W*W))/(6.28318*TSL)**(3./2)

                  !              Q1=SQRT((U-UC)**2+V*V+W*W)
                  !              Q2=SQRT(V*V+W*W)
                  !			 DFR=DIST(J,NF,L)*(XDIFR1*DM2*Q1+XDIFR2*DM1*Q2)

                  DFR=DIST(J,NF,L)*XFUNC1
                  DIFN(NF)=DIST(J,NF,L)-DFR
                  DEN1=2*DFR*DELT

                  FUNCS=DNS*EXP(-0.5/TSL*((U-UC)**2+V*V+W*W))/(6.28318*TSL)**(3./2)
                  SUMROST1=SUMROST1+CONDENSAT*DFR
                  SUMROST2=SUMROST2+FUNCS
                  !			 *(XDIFR1*DM2*Q1+XDIFR2*DM1*Q2)
                  SUMROST=SUMROST1-SUMROST2

                  SN1=SN1+DFR


                  NF1=0
                  DEN_NF=0.
                  DO 4 IU1=1,NUM
                     U1=SU(IU1)
                     DO 4 IV1=1,NVM
                        V1=SV(IV1)
                        DO 4 IW1=1,NWM
                           W1=SW(IW1)
                           NF1=NF1+1
                           PLOSK=U*U1+V*V1+W*W1
                           IF(PLOSK.LE.0.) THEN
                              !                   FUNC=DEN1*EXP(-0.5/TSL*((U1-UC)**2+V1*V1+W1*W1))/(6.28318*TSL)**(3./2)
                              !                   DIFR(NF1)=DIFR(NF1)+FUNC
                              DEN_NF=DEN_NF+DIST(J,NF1,L)
                           END IF
4                 CONTINUE

                  DEN_NF=2*DEN_NF*DELT
                  DIFR(NF)=DEN_NF*XFUNC1*FSURF1

3        CONTINUE

         SUMROST=SUMROST*DELT/FSOLID
         ROSTS(L)=SUMROST
         SUMROST2=SUMROST2*DELT

         D=DIAM(J,L)*DGAS
         VAL=VVV*SUMROST
         SORB(L)=SORB(L)+SUMROST

         QSORB=SUMROST*Q1M
         TSLD(L)=TSLD(L)+QSORB/TEPLSLD/T00

         !	    CYL=CL*DGAS
         !	    SURF=3.14159*(D**2)+3.14159*D*CYL
         !          ROST=VAL/SURF
         !	    D=D+2*ROST
         !          DIAM(J,L)=D/DGAS
         !	    IF(DIAM(J,L).LT.DMIN) DIAM(J,L)=DMIN

         NF=0.
         SN2=0.

         DO 5 IU=1,NUM
            U=SU(IU)
            DO 5 IV=1,NVM
               V=SV(IV)
               DO 5 IW=1,NWM
                  W=SW(IW)
                  NF=NF+1
                  SN2=SN2+DIFR(NF)
                  FUNCS=SUMROST2*EXP(-0.5/TSL*((U-UC)**2+V*V+W*W))/(6.28318*TSL)**(3./2)

                  DIST(J,NF,L)=DIFN(NF)+DIFR(NF)*(1-CONDENSAT)+FUNCS
5        CONTINUE

         XCORR=SN1/SN2
         !       PRINT *,SN1,SN2,XCORR

         GOTO 33
         !
         !
         !    CORRECTION
         !
         SN2=0.
         NF=0.
         DO 25 IU=1,NUM
            U=SU(IU)
            DO 25 IV=1,NVM
               V=SV(IV)
               DO 25 IW=1,NWM
                  W=SW(IW)
                  NF=NF+1
                  SN2=SN2+DIFR(NF)*XCORR
                  DIST(J,NF,L)=DIFN(NF)+DIFR(NF)*XCORR
25       CONTINUE

33 CONTINUE

   RETURN
END




SUBROUTINE TPREDEL(L,TSL,T_PREDEL,CONDENSAT)
   PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
   COMMON/BMOM/DENS(M0,MX),U0(M0,MX),V0(M0,MX),W0(M0,MX), &
      E0(M0,MX),PXX(M0,MX),PYY(M0,MX),PZZ(M0,MX),P(M0,MX), &
      TX(M0,MX),TY(M0,MX),TZ(M0,MX),T(M0,MX),XM(M0,MX), &
      XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
      ZJP(M0,MX),YJN(M0,MX),QX(M0,MX),EX(M0,MX)

   !       C=1.005
   !       T_PREDEL=C*P(1,L)
   !       IF(TSL.GE.T_PREDEL) CONDENSAT=0.00001

   BET1=1.
   T1=1.
   BET2=0.
   T2=1.1
   CONDENSAT=(BET1-BET2)*((T2-TSL)/(T2-T1))+BET2
   IF(CONDENSAT.LT.0.) CONDENSAT=0.

   RETURN
END



!       SUBROUTINE G3(L)
!       PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
!       COMMON/BTW1/ZW1,TW1,VW1
!	 COMMON/BMOM1/DENS1(M0,MX),U01(M0,MX),V01(M0,MX),W01(M0,MX),
!     * E01(M0,MX),PXX1(M0,MX),PYY1(M0,MX),PZZ1(M0,MX),P1(M0,MX),
!     * TX1(M0,MX),TY1(M0,MX),TZ1(M0,MX),T1(M0,MX),XM1(M0,MX),
!     * XJN1(M0,MX),XJP1(M0,MX),YJP1(M0,MX),XJX1(M0,MX),
!     * ZJP1(M0,MX),YJN1(M0,MX),QX1(M0,MX),EX1(M0,MX)
!	 COMMON/BMOM2/DENS2(M0,MX),U02(M0,MX),V02(M0,MX),W02(M0,MX),
!     * E02(M0,MX),PXX2(M0,MX),PYY2(M0,MX),PZZ2(M0,MX),P2(M0,MX),
!     * TX2(M0,MX),TY2(M0,MX),TZ2(M0,MX),T2(M0,MX),XM2(M0,MX),
!     * XJN2(M0,MX),XJP2(M0,MX),YJP2(M0,MX),XJX2(M0,MX),
!     * ZJP2(M0,MX),YJN2(M0,MX),QX2(M0,MX),EX2(M0,MX)
!       COMMON/BFYZ/FY(MY,MZ),FZ(MY,MZ)
!       COMMON/BJ/JY,JS,JD,JZ,JYZ/BUVW/U,V,W,NF,KF,NFR
!       COMMON/BTEMPL3/TL3(MX),VL3(MX)
!       FR=FMAXR(U,V,W)
!       DO 1 I=1,JZ
!       J=I*JY
!       RN=ZW1*YJP1(J,L)/VW1
!       TTL=TW1
!        RN=YJN(J,L)/VL3(L)
!        TTL=TL3(L)
!      FY(1,I)=FR
!       FY(1,I)=EXP(-0.5/TTL*(U*U+V*V+W*W))/(6.28318*TTL)**(1.5)*RN
!  1    CONTINUE
!       RETURN
!       END



!       SUBROUTINE G4(L)
!       PARAMETER(MFS=1730,MF0=260,M0=1,MX=42,MY=1,MZ=1)
!       COMMON/BTW1/ZW1,TW1,VW1
!	 COMMON/BMOM1/DENS1(M0,MX),U01(M0,MX),V01(M0,MX),W01(M0,MX),
!     * E01(M0,MX),PXX1(M0,MX),PYY1(M0,MX),PZZ1(M0,MX),P1(M0,MX),
!     * TX1(M0,MX),TY1(M0,MX),TZ1(M0,MX),T1(M0,MX),XM1(M0,MX),
!     * XJN1(M0,MX),XJP1(M0,MX),YJP1(M0,MX),XJX1(M0,MX),
!     * ZJP1(M0,MX),YJN1(M0,MX),QX1(M0,MX),EX1(M0,MX)
!	 COMMON/BMOM2/DENS2(M0,MX),U02(M0,MX),V02(M0,MX),W02(M0,MX),
!     * E02(M0,MX),PXX2(M0,MX),PYY2(M0,MX),PZZ2(M0,MX),P2(M0,MX),
!     * TX2(M0,MX),TY2(M0,MX),TZ2(M0,MX),T2(M0,MX),XM2(M0,MX),
!     * XJN2(M0,MX),XJP2(M0,MX),YJP2(M0,MX),XJX2(M0,MX),
!     * ZJP2(M0,MX),YJN2(M0,MX),QX2(M0,MX),EX2(M0,MX)
!       COMMON/BFYZ/FY(MY,MZ),FZ(MY,MZ)
!       COMMON/BJ/JY,JS,JD,JZ,JYZ/BUVW/U,V,W,NF,KF,NFR
!       COMMON/BTEMPL4/TL4(MX),VL4(MX)
!       FR=FMAXR(U,V,W)
!       DO 1 I=1,JZ
!       J=I*JY

!       RN=ZW1*YJP1(J,L)/VW1
!       TTL=TW1
!       FY(JY,I)=FR
!       FY(JY,I)=EXP(-0.5/TTL*(U*U+V*V+W*W))/(6.28318*TTL)**(1.5)*RN
!  1    CONTINUE
!       RETURN
!       END
!
!
!
!      SUBROUTINE G5
!        PARAMETER(MFS=282,MF0=65,M0=1,MX=42,MY=1,MZ=1)
!      COMMON/BUVW/U,V,W,NF,KF,NFR/BJ/JY,JS,JD,JZ,JYZ/BNKC/NK
!      COMMON/BHF1/HF1(MY,MZ,MF0)/BMOM4/XJN(M0,MX),XJP(M0,MX), &
!       YJP(M0,MX),XJX(M0,MX),ZJP(M0,MX)/BTW2/ZW2,TW2,VW2
!      COMMON/BJSS/JSS

! FOR "PLUME" WITH ORIFACE
!      JS1=JS+1
!      JSS1=JSS+1
!      FMR=FMAXR(U,V,W)
!      DO 5 I=1,JSS
!      DO 5 K=1,JS
!      HF1(K,I,KF)=FMAX(U,V,W)
!   5  CONTINUE
!      DO 10 I=JSS1,JZ
!      DO 10 K=JS1,JY
!     J=(I-1)*JY+K
!     RN=ZW2*XJN(J,2)/VW2
!     HF1(K,I,KF)=EXP(-0.5/TW2*(U*U+V*V+W*W))/(6.28318*TW2)**(3./2)*RN
!      HF1(K,I,KF)=FMR
!  10  CONTINUE
!      DO 15 I=1,JSS
!      DO 15 K=JS1,JY
!     J=(I-1)*JY+K
!     RN=ZW2*XJN(J,2)/VW2
!     HF1(K,I,KF)=EXP(-0.5/TW2*(U*U+V*V+W*W))/(6.28318*TW2)**(3./2)*RN
!      HF1(K,I,KF)=FMR
!  15  CONTINUE
!      DO 20 I=JSS1,JZ
!      DO 20 K=1,JS
!     J=(I-1)*JY+K
!     RN=ZW2*XJN(J,2)/VW2
!     HF1(K,I,KF)=EXP(-0.5/TW2*(U*U+V*V+W*W))/(6.28318*TW2)**(3./2)*RN
!      HF1(K,I,KF)=FMR
!  20  CONTINUE
!      RETURN
!      END


!      SUBROUTINE G6
!        PARAMETER(MFS=282,MF0=65,M0=1,MX=42,MY=1,MZ=1)
!      COMMON/BHF1/HF1(MY,MZ,MF0)/BSNT/RS,TS
!      COMMON/BNKC/NK/BUVW/U,V,W,NF,KF,NFR/BTW1/ZW1,TW1,VW1
!      COMMON/BMOM4/XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX), &
!       ZJP(M0,MX)/BJ/JY,JS,JD,JZ,JYZ
!      DO 10 I=1,JZ
!      DO 10 K=1,JY
!      J=(I-1)*JY+K
!      RN=ZW1*XJP(J,NK)/VW1
!      HF1(K,I,KF)=EXP(-0.5/TW1*(U*U+V*V+W*W))/(6.28318*TW1)**(3./2)*RN
!      HF1(K,I,KF)=EXP(-0.5/TS*(U*U+V*V+W*W))/(6.28318*TS)**(3./2)*RS
!      HF1(K,I,KF)=FMAXR(U,V,W)
!  10  CONTINUE
!      RETURN
!      END


!      SUBROUTINE G9(L)
!      PARAMETER(MFS=282,MF0=65,M0=1,MX=42,MY=1,MZ=1)
!      COMMON/BFYZ/FY(MY,MZ),FZ(MY,MZ)/BJ/JY,JS,JD,JZ,JYZ
!      COMMON/BMOM4/XJN(M0,MX),XJP(M0,MX),YJP(M0,MX),XJX(M0,MX),ZJP(M0,MX)
!      COMMON/BTW1/ZW1,TW1,VW1/BUVW/U,V,W,NF,KF,NFR
!      DO 10 K=1,JY
!      J=(JZ-1)*JY+K
!      RN=ZW1*ZJP(J,L)/VW1
!      FZ(K,JZ)=EXP(-0.5/TW1*(U*U+V*V+W*W))/(6.28318*TW1)**(3./2)*RN
!      FZ(K,JZ)=FMAXR(U,V,W)
!  10  CONTINUE
!      RETURN
!      END


subroutine flowns(rl_b,tl_b,ul_b,rhoa,temp,vel,cva,n)
   parameter (rx=5002)
   integer i,n,ktt
   real gama,cva,umix0,tmix0,rhomix0
   real rhoa(rx),temp(rx),vel(rx)
   real rl_b,tl_b,ul_b
   n=960

   ktt=1
   umix0=0.
   tmix0=1.
   rhomix0=1.


   ! ***** <a-Ar> **** <b-He> ******
   !  mamb=40./4.
   !  dadb=3.67/2.2
   !
   gama=5./3.
   !  gamb=5./3.
   ! ***** <a-Ar> **** <b-He> ******
   ! ***** <a-H20> ***** <b-N2> *****
   !  mamb=2.99/4.65
   !  dadb=4.34/3.7

   !  gama=9./7.
   !  gamb=7./5.
   ! ***** <a-H20> ***** <b-N2> *****
   cva=1./(gama-1.)


   IF(Ktt.EQ.2) THEN
      !******************************************************
      OPEN(15,FILE='init1_',STATUS='OLD', &
         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      READ(15) rhoa,temp,vel
      CLOSE(15)
      !******************************************************
   else
      do i=1, n
         rhoa(i)=rl_b
         vel(i)=ul_b
         temp(i)=tl_b
      enddo
   END IF

   open(1,file='init.dat',status='unknown')

   do i=1,n,1
      write(1,99) rhoa(i),temp(i),vel(i)
99    FORMAT(2F9.5,f11.7)
   end do

   close(2)
   return
end


subroutine fprint_ns(xx,rhoa,temp,vel,n)
   parameter (rx=5002)
   !  parameter (tend=3.47,kk=100)
   integer j
   real rhoa(rx),temp(rx),vel(rx),xx(rx)

   open(2,file='ev_con-n.dat',status='unknown')

   do j=1,n,1
      write(2,99)xx(j)/10., rhoa(j),temp(j),vel(j),vel(j)*rhoa(j)
99    FORMAT(F7.4,2F9.5,2f11.7)
   end do

   close(2)

   return
end


subroutine navie_stoks(rl_b,tl_b,ul_b,rl_ns,tl_ns,ul_ns, &
   dtau,rhoa,temp,cva,vel,xx,n,rr_b,tr_b,ur_b,rr_ns,tr_ns,ur_ns)

   real cva,cpa,ca,dx,time,dtau,umix0,tmix
   real rl_b,tl_b,ul_b,rl_ns,tl_ns,ul_ns
   real rr_b,tr_b,ur_b,rr_ns,tr_ns,ur_ns
   !  parameter (rx=1002,pi=3.1415926)
   parameter (rx=5002)
   !  parameter (tend=3.47,kk=100)
   integer j,n,i,jr,usl,usp
   real kk1,kk2
   real ei,gama,vsound,vmax
   real sumrha,pi,sumrhoe
   real rhoa(rx),rhou(rx),rhoe(rx),temp(rx),vel(rx),xx(rx)
   real rho(rx),preje(rx),rho1(rx)
   real a(rx),b(rx),l(rx)
   real a1(rx),c1(rx),b1(rx),f1(rx),y(0:rx),dxx(rx)
   real l1(rx),km(rx),kpl(rx),d1(rx),diffm(rx)
   real ulef,urig,wlef,wrig,blef,brig
   real rhol,rhor,tlef,trig,cmix(rx),cvmd(rx),cvmix(rx),kmd(rx)
   real rhoal,rhoar,rhoalef,rhoarig
   real tl,tp,ulefa,uriga,caa(rx)
   real flowdiffa(rx),kmdif(rx),kmu(rx),kmp(rx)
   real dudx(rx),flowdiffa1(rx)
   real kmda(rx),d2rhadx(rx)

   !  kk1=2.
   !  kk2=700.
   n=960

   !  rhoalef=1.974
   !  rhoarig=1.026

   !  rhoalef=1.1627
   !  rhoarig=0.995

   tlef=tl_b
   trig=tr_b

   !  open(2,file='he_tepl.dat',status='unknown')
   !  open(4,file='voda-azotflow.dat',status='unknown')

   pi=3.1415926

   jr=1

   !  umix0=0.

   !  dx=0.165
   dx=0.5
   !  dx=0.3

   !  tmix=1.

   !  dtau=0.03
   !  dtau=0.0001

   !  mamb=10./1.
   !  dadb=3.67/2.2

   ! ***** <a-Ar> **** <b-He> ******
   !  mamb=40./4.
   !  dadb=3.67/2.2

   !  gama=5./3.
   !  gamb=5./3.
   ! ***** <a-Ar> **** <b-He> ******
   ! ***** <a-H20> ***** <b-N2> *****
   !  mamb=2.99/4.65
   !  dadb=4.34/3.7

   !  gama=9./7.
   !  gamb=7./5.
   ! ***** <a-H20> ***** <b-N2> *****
   !  cva=1./(gama-1.)
   !  cvb=mamb/(gamb-1.)


   !	vsound=dsqrt(gamb*tmix)
   !	dtau=0.3*dx/dabs(vsound)/kk1

   ulefa=ul_b
   uriga=ur_b

   xx(1)=10.+dx/2.
   do j=2,n
      xx(j)=xx(j-1)+dx
   enddo

   !	dxx(1)=dx/2.

   do j=1,n+1
      dxx(j)=dx
   enddo

   !	dxx(n+1)=dx/2.

   do i=1, n
      !		rhoa(i)=1.0
      !		u(i)=umix0
      rhou(i)=rhoa(i)*vel(i)
      !		temp(i)=tmix
      rhoe(i)=cva*temp(i)*rhoa(i)
      preje(i)=rhoa(i)*temp(i)
   enddo

   !	time=0.

   rhoal=rl_b
   rhoar=rr_b

   do j=1,n
      a(j)=rhoe(j)
      b(j)=vel(j)
      l(j)=rhoa(j)*temp(j)
   enddo

   blef=ul_b
   brig=ur_b

   !		wlef=rhoal*tlef*cva*ulefa*dtau/dx
   wlef=rl_b*tl_b*cva*ul_b*dtau/dx
   wrig=rr_b*tr_b*cva*ur_b*dtau/dx
   !		wrig=rhoar*trig*cva*uriga*dtau/dx

   call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)

   do j=1,n
      rhoe(j)=a(j)
   enddo

   ! ****************************************************
   do j=1,n
      a(j)=rhou(j)
      b(j)=rhoa(j)*temp(j)
      l(j)=1.
   enddo


   blef=rl_b*tl_b
   brig=rr_b*tr_b


   !		wlef=rhoal*ulefa*ulefa*dtau/dx
   wlef=rl_b*ul_b*ul_b*dtau/dx
   wrig=rr_b*ur_b*ur_b*dtau/dx
   !		wrig=rhoar*uriga*uriga*dtau/dx


   call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)

   do j=1,n
      rhou(j)=a(j)
   enddo

   ! ******************************************************

   do j=1,n
      a(j)=rhoa(j)
      b(j)=0.
      l(j)=1.
   enddo

   blef=0.
   brig=0.


   !		wlef=rhoal*ulefa*dtau/dx
   wlef=rl_b*ul_b*dtau/dx
   wrig=rr_b*ur_b*dtau/dx
   !		wrig=rhoar*uriga*dtau/dx

   call shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)

   do j=1,n
      rhoa(j)=a(j)
   enddo

   ! ******************************************************
   do j=1,n
      vel(j)=rhou(j)/rhoa(j)
      temp(j)=rhoe(j)/rhoa(j)/cva
   enddo

   ! *****************momentum equation **************************
   do j=2,n
      kmu(j)=5./16.*sqrt(2.*pi)*(sqrt(temp(j))+sqrt(temp(j-1)))/2.
   enddo

   do j=1,n
      l1(j)=rhoa(j)
      y(j)=vel(j)
      d1(j)=0.
   enddo

   usl=1
   usp=1

   !	tl=2.*u(1)-u(2)
   !	tp=2.*u(n)-u(n-1)

   !	tl=ulef
   !	tp=urig

   tl=ul_b
   tp=ur_b

   dtau=4./3.*dtau

   call teplo(kmu,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

   do j=1,n
      vel(j)=y(j)
   enddo

   !	rl_ns=rhoa(1)
   !	ul_ns=y(0)

   dtau=3./4.*dtau
   !	kmu=km

   ! **********energy equation ************************************
   usl=1
   usp=1

   tl=tl_b
   tp=tr_b

   do j=2,n
      kmp(j)=5./16.*sqrt(2.*pi)*5./2.*cva*(sqrt(temp(j))+sqrt(temp(j-1)))/2.
   enddo

   do j=1,n
      l1(j)=rhoa(j)*cva
      y(j)=temp(j)
      d1(j)=0.
   enddo

   do j=2,n-1
      dudx(j)=(vel(j+1)-vel(j-1))/2./dxx(j)
   enddo

   dudx(1)=(vel(2)-ul_b)/2./dxx(2)
   dudx(n)=(ur_b-vel(n-1))/2./dxx(n-1)

   call teplo(kmp,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)

   do j=1,n
      temp(j)=y(j)
   enddo

   do j=1,n
      temp(j)=dtau*(4./3.*kmu(j)*dudx(j)*dudx(j))+temp(j)
   enddo

   ! **************************************************


   !	do i=1, n
   !		rhou(i)=rhoa(i)*u(i)
   !		rhoe(i)=rhoa(i)*temp(i)*cva
   !	enddo

   rl_ns=rhoa(1)
   tl_ns=temp(1)
   ul_ns=vel(1)

   rr_ns=rhoa(n)
   tr_ns=temp(n)
   ur_ns=vel(n)

   return
end



subroutine shasta(a,b,l,dx,vel,dtau,n,blef,brig,wlef,wrig)
   integer j,n
   parameter (rx=5002)
   real dtau,dtdx,dx
   real qplus(rx),qmin(rx),splus(rx),smin(rx),daudx(rx),dbdx(rx)
   real a1(rx),a2(rx),a3(rx),a4(rx)
   real a(rx),b(rx),l(rx),vel(rx)
   real ac(rx),dp(rx),dpl(rx),fc(rx)
   real koef1,koef2,koef3,koef4,koef5
   real blef,brig,wlef,wrig

   dtdx=dtau/dx

   do j=2,n-1
      qplus(j)=(1./2.-vel(j)*dtdx)/(1.+(vel(j+1)-vel(j))*dtdx)
      qmin(j)=(1./2.+vel(j)*dtdx)/(1.-(vel(j-1)-vel(j))*dtdx)

      splus(j)=l(j)*(b(j+1)-b(j))*dtdx
      smin(j)=l(j)*(b(j)-b(j-1))*dtdx

      a1(j)=qplus(j)*qplus(j)*(a(j+1)-a(j))/2.
      a2(j)=qmin(j)*qmin(j)*(a(j-1)-a(j))/2.
      a3(j)=qplus(j)*(a(j)-splus(j))
      a4(j)=qmin(j)*(a(j)-smin(j))
   enddo

   do j=2,n-1
      ac(j)=a1(j)+a2(j)+a3(j)+a4(j)
      !  	a(j)=a1(j)+a2(j)+a3(j)+a4(j)
   enddo

   qplus(1)=(1./2.-vel(1)*dtdx)/(1.+(vel(2)-vel(1))*dtdx)
   splus(1)=l(1)*(b(2)-blef)*dtdx
   a1(1)=qplus(1)*qplus(1)*(a(2)-a(1))/2.
   a3(1)=qplus(1)*(a(1)-splus(1))

   !	a(1)=a1(1)+a3(1)+a(1)/2.+wlef
   a(1)=a1(1)+a3(1)+a(1)/2.+wlef+(blef-b(1))*l(1)*dtdx

   qmin(n)=(1./2.+vel(n)*dtdx)/(1.-(vel(n-1)-vel(n))*dtdx)
   smin(n)=l(n)*(brig-b(n-1))*dtdx
   a2(n)=qmin(n)*qmin(n)*(a(n-1)-a(n))/2.
   a4(n)=qmin(n)*(a(n)-smin(n))

   !	a(n)=a2(n)+a4(n)+a(n)/2.-wrig
   a(n)=a2(n)+a4(n)+a(n)/2.-wrig+(b(n)-brig)*l(n)*dtdx

   ! ************ anti diffusion step ***********************
   do j=2,n-1
      dp(j)=ac(j+1)-ac(j)
   enddo

   do j=2,n-1
      if(dp(j)<0.)then
         dpl(j)=-1.
      elseif (dp(j)>0.)then
         dpl(j)=1.
      else
         dpl(j)=0.
      endif
   enddo

   do j=2,n-1

      koef1=dpl(j)*dp(j-1)
      koef2=abs(dp(j))
      koef2=koef2/8.
      koef3=dpl(j)*dp(j+1)

      koef4=min1(koef1,koef2,koef3)
      koef5=max1(0.0,koef4)

      fc(j)=dpl(j)*koef5

   enddo

   fc(1)=0.
   fc(n)=0.
   fc(n-1)=0.

   do j=2,n-1
      !  a(j)=ac(j)
      a(j)=ac(j)-(fc(j)-fc(j-1))
   enddo
   !************ end of anti diffusion step ***********************

   return
end subroutine shasta


subroutine teplo(km,l1,d1,usl,usp,dtau,tl,tp,n,y,dxx)
   parameter (rx=5002)
   real a1(rx),b1(rx),c1(rx),l1(rx),alpha(rx),betta(rx)
   real y(0:rx),f1(rx),d1(rx),km(rx),dxx(rx)
   real dtau,tl,tp
   real kappa2,kappa1
   integer j,n,usl,usp

   km(1)=km(2)
   !  km(n)=km(n-1)
   km(n+1)=km(n)
   do j=1,n
      a1(j)=2.*dtau*km(j)/dxx(j)/(dxx(j+1)+dxx(j))
      b1(j)=2.*dtau*km(j+1)/dxx(j+1)/(dxx(j+1)+dxx(j))
      c1(j)=a1(j)+b1(j)+l1(j)
      f1(j)=l1(j)*y(j)+d1(j)
   enddo

   if(usl==2)then
      alpha(1)=1.
      betta(1)=0.
   else
      alpha(1)=0.
      betta(1)=tl
   endif

   do j=1,n
      alpha(j+1)=b1(j)/(c1(j)-a1(j)*alpha(j))
      betta(j+1)=(a1(j)*betta(j)+f1(j))/(c1(j)-a1(j)*alpha(j))
   enddo

   if(usp==2)then
      kappa2=1.
      kappa1=0.
   else
      kappa2=0.
      kappa1=tp
   endif

   y(n+1)=(kappa2*betta(n+1)+kappa1)/(1.-kappa2*alpha(n+1))
   do j=n,0,-1
      y(j)=alpha(j+1)*y(j+1)+betta(j+1)
   enddo

   return
end subroutine teplo
