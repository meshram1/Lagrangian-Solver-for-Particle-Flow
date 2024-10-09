      PROGRAM LAGRANGIAN_2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN= 5000000, NB= 100000)
!  
      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK


      OPEN(UNIT = 25, FILE = "particle1.dat") 
cc    file for particle data
      OPEN(UNIT = 30, FILE = "flow.in")
cc    input file
      OPEN(UNIT = 41, FILE = 'x_particle.dat')
      OPEN(UNIT = 42, FILE = "y_particle.dat")
      OPEN(UNIT = 43, FILE = "z_particle.dat")
      OPEN(UNIT = 46, FILE = "u_particle.dat")
      OPEN(UNIT = 47, FILE = "v_particle.dat")
      OPEN(UNIT = 48, FILE = "w_particle.dat")
cc    the next few files are files with the coordinates of the particle
      OPEN(UNIT = 81, FILE = "XY_particle.dat")
      OPEN(UNIT = 85, FILE = "XY_mm_particle.dat")
      OPEN(UNIT = 84, FILE = "XZ.dat")
      OPEN(UNIT = 87, FILE = "XZ_mm.dat")

      PI = 22.0d0/7.0d0
      READ(30,*)R_P,AMOL,DT,NITER,RO_P,X_P_INI,Y_P_INI,Z_P_INI      
      AM_P = RO_P*(4.0d0/3.0d0)*PI*R_P**3.0d0

      CALL BOUNDARY

C      PRINT *, NODES
CCV INITIALISATION
      X_P_OLD = X_P_INI
      Y_P_OLD = Y_P_INI !!max_y = 20.3 mm at throat
      Z_P_OLD = Z_P_INI

      U_P_OLD = 0.0d0
      V_P_OLD = 0.0d0
      W_P_OLD = 0.0d0

      X_P = X_P_OLD
      Y_P = Y_P_OLD
      Z_P = Z_P_OLD

      U_P = U_P_OLD
      V_P = V_P_OLD
      W_P = W_P_OLD
      CHECK = 0
      TIME = 0.0d0
C INITIALISATION (AGAIN)
      DO ITER = 1,NITER
       TIME = TIME + DT   

!      DT = 1e-3


C RK2 FOR VELOCITY COMPUTATION
        CALL CLOSEST_NODE(X_P_OLD,Y_P_OLD,Z_P_OLD)
        CALL DOMAIN
        CALL INTERPOLATION
        CALL SOURCES


C IF THE PARTICLE REACHES BOUNDARY OF WALL OR OUTLET WE STOP AND RECORD ITS POSITION.

         
       IF(CHECK.EQ.1)THEN
        X_P_OLD = X(MIN_INDEX)        
        Y_P_OLD = Y(MIN_INDEX)
        Z_P_OLD = Z(MIN_INDEX)
        X_P_MM = X_P_OLD*1000.d0        
        Y_P_MM = Y_P_OLD*1000.d0
        Z_P_MM = Z_P_OLD*1000.d0
                       
        WRITE(25,*)TIME,X_P_OLD,Y_P_OLD,Z_P_OLD,U_P_OLD,V_P_OLD,W_P_OLD
        WRITE(41,*)TIME,X_P_OLD
        WRITE(42,*)TIME,Y_P_OLD
        WRITE(43,*)TIME,Z_P_OLD
        WRITE(46,*)TIME,U_P_OLD
        WRITE(47,*)TIME,V_P_OLD
        WRITE(48,*)TIME,W_P_OLD
        WRITE(81,*)X_P_OLD,Y_P_OLD        
        WRITE(84,*)X_P_OLD,Z_P_OLD
        WRITE(85,*)X_P_MM,Y_P_MM
        GOTO 100
        
       ELSEIF(CHECK.EQ.2)THEN
        GOTO 100
C RK2 FIRST STEP FOR VELOCITIES
       ELSEIF(CHECK.EQ.0)THEN
        CALL CLOSEST_NODE(X_P_OLD,Y_P_OLD,Z_P_OLD)
        CALL DOMAIN
        CALL INTERPOLATION
        CALL SOURCES

        AK1_U_P = DT*DRAG_X/AM_P
        U_P = U_P_OLD + 0.5d0*AK1_U_P

        AK1_V_P = DT*DRAG_Y/AM_P
        V_P = V_P_OLD + 0.5d0*AK1_V_P
        
        AK1_W_P = DT*DRAG_Z/AM_P
        W_P = W_P_OLD + 0.5d0*AK1_W_P

C RK2 SECOND STEP VELOCITIES
        CALL CLOSEST_NODE(X_P_OLD,Y_P_OLD,Z_P_OLD)
        CALL DOMAIN
        CALL INTERPOLATION 
        CALL SOURCES

        AK2_U_P = DT*DRAG_X/AM_P
        U_P = U_P_OLD + AK2_U_P

        AK2_V_P = DT*DRAG_Y/AM_P
        V_P = V_P_OLD + AK2_V_P

        AK2_W_P = DT*DRAG_Z/AM_P
        W_P = W_P_OLD + AK2_W_P


C MOVE PARTICLES USING UPDATED VELOCITIES
        X_P = X_P_OLD + U_P*DT
        Y_P = Y_P_OLD + V_P*DT
        Z_P = Z_P_OLD + Z_P*DT

C       UPDATING OLD VALUES WITH THE NEW ONES FOR VELOCITIES

        U_P_OLD = U_P
        V_P_OLD = V_P
        W_P_OLD = W_P
        X_P_OLD = X_P
        Y_P_OLD = Y_P
        Z_P_OLD = Z_P

        X_P_MM = X_P_OLD*1000.d0        
        Y_P_MM = Y_P_OLD*1000.d0
        Z_P_MM = Z_P_OLD*1000.d0
                       
        WRITE(25,*)TIME,X_P_OLD,Y_P_OLD,Z_P_OLD,U_P_OLD,V_P_OLD,W_P_OLD
        WRITE(41,*)TIME,X_P_OLD
        WRITE(42,*)TIME,Y_P_OLD
        WRITE(43,*)TIME,Z_P_OLD
        WRITE(46,*)TIME,U_P_OLD
        WRITE(47,*)TIME,V_P_OLD
        WRITE(48,*)TIME,W_P_OLD
        WRITE(81,*)X_P_OLD,Y_P_OLD        
        WRITE(84,*)X_P_OLD,Z_P_OLD
        WRITE(85,*)X_P_MM,Y_P_MM

c        PRINT *,TIME,U_P_OLD,V_P_OLD,W_P_OLD
       ENDIF
       
       CALL CLOSEST_NODE(X_P_OLD,Y_P_OLD,Z_P_OLD)
c       WRITE(53,*)MIN_INDEX,X_P_OLD,Y_P_OLD,Z_P_OLD
       DO I = 1,NB_T5
            IF(MIN_INDEX.EQ.BN_WALL(I))THEN !!! CLOSE TO A NODE AT THE BOUNDARY
                U_P = 0.0d0
                V_P = 0.0d0
                W_P = 0.0d0
                U_P_OLD = 0.0d0
                V_P_OLD = 0.0d0
                W_P_OLD = 0.0d0
c                DRAG_X = 0.0d0
c                DRAG_Y = 0.0d0
c                DRAG_Z = 0.0d0
                CHECK = 1
            ENDIF
        ENDDO

         DO I = 1,NB_T3
            IF(MIN_INDEX.EQ.BN_OUTLET(I))THEN !!! CLOSE TO A NODE AT THE BOUNDARY
                U_P = 0.0d0
                V_P = 0.0d0
                W_P = 0.0d0
                U_P_OLD = 0.0d0
                V_P_OLD = 0.0d0
                W_P_OLD = 0.0d0
c                DRAG_X = 0.0d0
c                DRAG_Y = 0.0d0
c                DRAG_Z = 0.0d0
                CHECK = 2
            ENDIF
         ENDDO
      ENDDO
100   CONTINUE
         CLOSE(25)
         CLOSE(41)
         CLOSE(42)
         CLOSE(43)
         CLOSE(46)
         CLOSE(47)
         CLOSE(48)
         CLOSE(68)
C      RETURN
      END
    
C ---------------------------------------------------------------------------
      SUBROUTINE BOUNDARY !!! PROVIDES THE BOUNDARY NODES AND NUMBER OF NODES AT WALL, AT INLET AND OUTLET.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NN= 5000000, NB= 100000)
    
      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK

    
      OPEN(unit=10,file='grid.dat')
      OPEN(unit=11,file='bc.in')
      OPEN(UNIT=51, FILE = 'BOUNDARIES_WALL.dat')
      OPEN(UNIT=52, FILE = 'BOUNDARIES_OUTLET.dat')
      OPEN(UNIT = 66, FILE = "z_plane.dat")
c      OPEN(UNIT = 61, FILE = "boundary_coor.dat")

      READ(10,*)NODES,NELES
      DO I = 1, NODES
       READ(10,*)N,X(N),Y(N),Z(N)
       X(N) = X(N)*0.001d0
       Y(N) = Y(N)*0.001d0
       Z(N) = Z(N)*0.001d0
c       PRINT *,X(N) its reading x,
      ENDDO
      DO I = 1, NODES
        P_MAX_Z = Z_P_INI + 0.00005
        P_MIN_Z = Z_P_INI - 0.00005
      IF((Z(I).LE. P_MAX_Z).AND.(Z(I).GE.P_MIN_Z))THEN
        WRITE(66,*)X(I),Y(I)
        PRINT *,X(I),Y(I)
       END IF   
      END DO

      DO I = 1, NELES
       READ(10,*)NEL,NOD(NEL,1),NOD(NEL,2),NOD(NEL,3),NOD(NEL,4),
     *  NC1(NEL),NC2(NEL),NC3(NEL),NC4(NEL)
C     PRINT *,NEL
      ENDDO

      READ(11,*)NGHOSTS
        DO I = 1, NGHOSTS
          READ(11,*)NPARENT(I),NGHOST(I),NTYPE(I),NSIDE(I)
        ENDDO

      CLOSE(10)
      CLOSE(11)
      NB_T5 = 0
      NB_T2 = 0
      NB_T3 = 0

      DO I = 1, NGHOSTS
        NP = NPARENT(I)
        NT = NTYPE(I)
        NS = NSIDE(I)
         IF (NS.EQ.1)THEN
            NODE1 = NOD(NP,2)
            NODE2 = NOD(NP,3)
            NODE3 = NOD(NP,4)
         ELSEIF (NS.EQ.2)THEN
            NODE1 = NOD(NP,1)
            NODE2 = NOD(NP,3)
            NODE3 = NOD(NP,4)
         ELSEIF (NS.EQ.3)THEN
            NODE1 = NOD(NP,1)
            NODE2 = NOD(NP,2)
            NODE3 = NOD(NP,4)
         ELSEIF (NS.EQ.4)THEN
            NODE1 = NOD(NP,1)
            NODE2 = NOD(NP,2)
            NODE3 = NOD(NP,3)
         ENDIF

         IF(NT.EQ.5)THEN    
            BN_WALL(NB_T5+1) = NODE1
            BN_WALL(NB_T5+2) = NODE2
            BN_WALL(NB_T5+3) = NODE3
            NB_T5 = NB_T5 + 3
         ELSEIF(NT.EQ.2)THEN
            BN_INLET(NB_T2+1) = NODE1
            BN_INLET(NB_T2+2) = NODE2
            BN_INLET(NB_T2+3) = NODE3
            NB_T2 = NB_T2 + 3
         ELSEIF(NT.EQ.3)THEN
            BN_OUTLET(NB_T3+1) = NODE1
            BN_OUTLET(NB_T3+2) = NODE2
            BN_OUTLET(NB_T3+3) = NODE3
            NB_T3 = NB_T3 + 3
         ENDIF
C         PRINT *, NODE1,NODE2,NODE3
      ENDDO
c      PRINT *,NB_T5,NB_T3
      DO I = 1,NB_T5
        WRITE(51,*)NB_T5,BN_WALL(I)
c     WRITE(61,*)X(INT(BN_WALL(I))), Y(INT(BN_WALL(I))),
c     * Z(INT(BN_WALL(I)))
      ENDDO
    
      CLOSE(51)
      CLOSE(61)
      CLOSE(62)
      CLOSE(63)
      DO I = 1,NB_T3
        WRITE(52,*)NB_T3,BN_OUTLET(I)
      ENDDO
      CLOSE(52)
      CLOSE(10)
      CLOSE(11)
      CLOSE(66)
      RETURN
      END
!!!  WE HAVE COMPUTED THE BOUNDARY NODES AT THE WALL AND THE INLET AND OUTLET
!!!   PROVIDES: BN_OUTLET, BN_INLET, BN_WALL

C ------------------------------------------------------------------------------

!!! PROGRAM TO FIND THE CLOSEST NODE (GIVES INDEX)
      SUBROUTINE CLOSEST_NODE (XPAR,YPAR,ZPAR)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NN= 5000000, NB= 100000)

      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK


      DMIN = +1.0e12
      DO I = 1, NODES
       DIST(I) = DSQRT((X(I)-XPAR)**2 + (Y(I)-YPAR)**2 + (Z(I)-ZPAR)**2)
        IF (DIST(I).LE.DMIN) THEN
            DMIN = DIST(I)
            MIN_INDEX = I
        ENDIF
      ENDDO

      RETURN
      END

C --------------------------------------------------------------------------------
! TO FIND THE INTERPOLATED VALUES OF FLUID PROPERTIES AT THE PARTICLE LOCATION FROM THE NODES SUBROUTINE INTERPOLATION
    
      SUBROUTINE INTERPOLATION
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NN= 5000000, NB= 100000)

      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK



      OPEN(UNIT=20,FILE='t_plot.dat')


      DO I = 1,NODES
        READ(20,*)FJX,FJY,FJZ,U(I),V(I),W(I),P(I),RO(I),AM(I),T(I)
      END DO

      NOD1 = MIN_INDEX
! SINCE THE CLOSEST NODE WOULD BE THE FIRST OF THE SORTED NOD_SD_ARRAY
      NOD2 = NOD_SD(2)
      NOD3 = NOD_SD(3)
      NOD4 = NOD_SD(4)
      NOD5 = NOD_SD(5)

      X1 = X(NOD1)
      X2 = X(NOD2)
      X3 = X(NOD3)
      X4 = X(NOD4)
      X5 = X(NOD5)

      Y1 = Y(NOD1)
      Y2 = Y(NOD2)
      Y3 = Y(NOD3)
      Y4 = Y(NOD4)
      Y5 = Y(NOD5)

      Z1 = Z(NOD1)
      Z2 = Z(NOD2)
      Z3 = Z(NOD3)
      Z4 = Z(NOD4)
      Z5 = Z(NOD5)

      U1 = U(NOD1)
      U2 = U(NOD2)
      U3 = U(NOD3)
      U4 = U(NOD4)
      U5 = U(NOD5)
    
      V1 = V(NOD1)
      V2 = V(NOD2)
      V3 = V(NOD3)
      V4 = V(NOD4)
      V5 = V(NOD5)

      W1 = W(NOD1)
      W2 = W(NOD2)
      W3 = W(NOD3)
      W4 = W(NOD4)
      W5 = W(NOD5)

      RO1 = RO(NOD1)
      RO2 = RO(NOD2)
      RO3 = RO(NOD3)
      RO4 = RO(NOD4)
      RO5 = RO(NOD5)
    
      T1 = T(NOD1)
      T2 = T(NOD2)
      T3 = T(NOD3)
      T4 = T(NOD4)
      T5 = T(NOD5)

      D1 = DSQRT((X1 - X_P)**2 + (Y1 - Y_P)**2 + (Z1 - Z_P)**2)
      D2 = DSQRT((X2 - X_P)**2 + (Y2 - Y_P)**2 + (Z2 - Z_P)**2)
      D3 = DSQRT((X3 - X_P)**2 + (Y3 - Y_P)**2 + (Z3 - Z_P)**2)
      D4 = DSQRT((X4 - X_P)**2 + (Y4 - Y_P)**2 + (Z4 - Z_P)**2)
      D5 = DSQRT((X5 - X_P)**2 + (Y5 - Y_P)**2 + (Z5 - Z_P)**2)

      WT1 = 1/D1
      WT2 = 1/D2
      WT3 = 1/D3
      WT4 = 1/D4
      WT5 = 1/D5
      
      WTSUM = (WT1 + WT2 + WT3 + WT4 + WT5)
      
      U_F = (WT1*U1 + WT2*U2 + WT3*U3 + WT4*U4 + WT5*U5)/(WTSUM)
      V_F = (WT1*V1 + WT2*V2 + WT3*V3 + WT4*V4 + WT5*V5)/(WTSUM)
      W_F = (WT1*W1 + WT2*W2 + WT3*W3 + WT4*W4 + WT5*W5)/(WTSUM)
      RO_F = (WT1*RO1 + WT2*RO2 + WT3*RO3 + WT4*RO4 + WT5*RO5)/(WTSUM)
      T_F = (WT1*T1 + WT2*T2 + WT3*T3 + WT4*T4 + WT4*T5)/(WTSUM)
      
      
      CLOSE(19)
      CLOSE(20)
      RETURN
      END
C------------------------------------------------------------------------------------------------------------


      SUBROUTINE SOURCES
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NN= 5000000, NB= 100000)

      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK


        IF (CHECK.EQ.1) THEN
          DRAG_X = 0.0D0
          DRAG_Y = 0.0D0
          DRAG_Z = 0.0D0
          U_P = 0.0D0
          V_P = 0.0D0
          W_P = 0.0D0
          U_F = 0.0D0
          V_F = 0.0D0
          W_F = 0.0D0
      ELSEIF(CHECK.EQ.0)THEN
        C43 = 4.0d0/3.0d0
        C13 = 1.0d0/3.0d0

    !!! CALL INTERPOLATION

            
        REL_U = U_F - U_P !
        REL_V = V_F - V_P !
        REL_W = W_F - W_P !
!        REL_T = T_F - T_P !

        CMU = 11.8480e-8*DSQRT(AMOL)
        AMU = CMU*T_F**0.6d0
        
        REYNOLDS= (RO_F*R_P*2.0d0*DSQRT(REL_U**2+REL_V**2+REL_W**2))/AMU !!!NEED TO INPUT PARTICLE RADIUS
        
        IF(REYNOLDS.LE.0.1d0)THEN
            CD = 24.0d0/(REYNOLDS + 1e-3)
        ELSEIF ((REYNOLDS.GE.0.1d0).AND.(REYNOLDS.LE.1000))THEN
            CD = 0.48d0 + 24*((REYNOLDS)**(-0.85d0))
        ELSEIF (REYNOLDS.GE.1000)THEN
            CD = 0.48d0
        ENDIF
        PI = 22.0d0/7.0d0
        A = 0.5d0*CD*RO_F*PI*R_P**2.0d0
        DRAG_X = A*REL_U*DABS(REL_U)
        DRAG_Y = A*REL_V*DABS(REL_V)
        DRAG_Z = A*REL_W*DABS(REL_W)
        ENDIF
C        PRINT*,REL_U,REL_V,REL_W
C        PRINT*,DRAG_X, DRAG_Y, DRAG_Z
C        PRINT*,REYNOLDS
        PRINT*," "
        RETURN
        END

c-----------------------------------------------------------------------------------------------------------------------

      SUBROUTINE DOMAIN
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NN= 5000000, NB= 100000)

      COMMON/ABC/ X(NN), Y(NN), Z(NN), XC(NN), YC(NN), ZC(NN), 
     *  NOD(NN,4), NC1(NN), NC2(NN), NC3(NN), NC4(NN), VOL(NN),
     *  DL(NN), SC1X(NN), SC1Y(NN), SC1Z(NN), SC2X(NN), SC2Y(NN),
     *  SC2Z(NN), SC3X(NN), SC3Y(NN), SC3Z(NN),SC4X(NN), SC4Y(NN),
     *  SC4Z(NN), ITEST(NN),NODES,NELES,DIST(NN),IND_CLOSEST_NODE(NB),
     *  N1(NN),N2(NN),N3(NN),N4(NN),N5(NN),NOD_ELEM(NN),NOD_SD(NB),
     *  IND_NOD_CELL(NB,NB),X_SD(NB),Y_SD(NB),Z_SD(NB), MIN_INDEX

      COMMON/BBC/ NPARENT(NB), NGHOST(NB), NTYPE(NB), NSIDE(NB),
     * NGHOSTS, NB_T5 ,NB_T2, NB_T3,BN_OUTLET(NB), BN_INLET(NB),
     *  BN_WALL(NB),X_P_INI,Y_P_INI,Z_P_INI
    
      COMMON/CBC/ U(NN),V(NN) ,W(NN),P(NN),RO(NN), AM(NN), T(NN),
     * U_F, V_F, W_F, RO_F, T_F, P_F, M_F, X_P, Y_P, Z_P, U_P, V_P, W_P, 
     * DRAG_X, DRAG_Y, DRAG_Z,
     * X_P_OLD, Y_P_OLD,Z_P_OLD, U_P_OLD, V_P_OLD, W_P_OLD 

      COMMON/PROPERTIES/P0, T0, GAMA, AMOL, PRANDTL, CFL, NITER, ISTART,
     *  AMACH, CV_G, CP_G, CP_P, RR, CMU, S0, RO0, RO_P, AM_P, R_P,CHECK

c      INTEGER, ALLOCATABLE :: IND_NOD_CELL(:,:)
           
      OPEN(UNIT = 88, FILE = "coordinates_input.in")
	  OPEN(UNIT = 89, FILE = "subdomain.dat")

      DO I = 1, NELES
        READ(88,*) J1,J2,J3,NOD(I,1),NOD(I,2),NOD(I,3),NOD(I,4)
      ENDDO

      NOD_ELEM = 0  !Gives no. of elements sharing a node
      DO I = 1,NELES
	   DO J = 1,4		
	     NQ = NOD(I,J)
	     NOD_ELEM(NQ) = NOD_ELEM(NQ) + 1	! no. of cells for each node	
		 ENDDO
	    ENDDO

      MAX_SHARE_CELL = MAXVAL(NOD_ELEM)
c      ALLOCATE (IND_NOD_CELL(NODES,MAX_SHARE_CELL))
      
cc Now we have read the nodes to each cell, We find maximum number of cells connected.
      NOD_ELEM = 0
      IND_NOD_CELL = 0
C      DO I = 1,NODES
C        DO J = 1,MAX_SHARE_CELL
C          IND_NOD_CELL(I,J) = 0
C        ENDDO
C      ENDDO
      DO I = 1,NELES
        DO J = 1,4
          NQ = NOD(I,J)
          NOD_ELEM(NQ) = NOD_ELEM(NQ) + 1
          IND_NOD_CELL(NQ,NOD_ELEM(NQ)) = I
        ENDDO
      ENDDO

      ! Applying the closest node algorithm to find the closest node
      M1 = 0
      DO I = 1,NOD_ELEM(MIN_INDEX)
      NEL1 = IND_NOD_CELL(MIN_INDEX,I) 
      ! Defining Domain Nodes as the nodal points of the subdomain
       DO J = 1,4
        NODE1= NOD(NEL1,J)
        M1 = M1+1
        WRITE(89,*) NODE1, X(NODE1),Y(NODE1),Z(NODE1)!NOD_ELEM(MIN_INDEX),MIN_INDEX,J 
        NOD_SD(M1) = NODE1
        X_SD(M1)=X(NODE1)
        Y_SD(M1)=Y(NODE1)
        Z_SD(M1)=Z(NODE1)
       ENDDO
      ENDDO
   
      TOTAL_SUBNODES = NOD_ELEM(MIN_INDEX)*4.0D0
C      PRINT *, TOTAL_SUBNODES,M1
C      DO I = 1, TOTAL_SUBNODES
C       READ(89,*)NOD_SD(I),X_SD(I),Y_SD(I),Z_SD(I)
C       PRINT *, NOD_SD(I),X_SD(I),Y_SD(I),Z_SD(I)
C      END DO

      ! SO THE SUBDOMAIN NODES AND THEIR COORDINATES IS KNOWN
      DMIN = +1.0e+12
      PARA = 4
      !DO I = 1, NODES
!       DIST(I) = DSQRT((X(I)-XPAR)**2 + (Y(I)-YPAR)**2 + (Z(I)-ZPAR)**2)
!        IF (DIST(I).LE.DMIN) THEN
!            DMIN = DIST(I)
!            MIN_INDEX = I
!        ENDIF
!      ENDDO
        DO I = 1, TOTAL_SUBNODES
        DIST(I) = DSQRT((X(MIN_INDEX)-X_SD(I))**2 +
     * (Y(MIN_INDEX)-Y_SD(I))**2 + (Z(MIN_INDEX) - Z_SD(I))**2)
         ENDDO
        TEMP = 0
        ! SORTING NODE_sD BASED ON DISTANCE
        DO I = 1, TOTAL_SUBNODES - 1
         DO J = I + 1, TOTAL_SUBNODES
           IF (DIST(I).GT.DIST(J)) THEN
             TEMP = NOD_SD(I)
             NOD_SD(I) = NOD_SD(J)
             NOD_SD(J) = TEMP
           END IF
         ENDDO         
        ENDDO

      CLOSE(89)
      CLOSE(88)
      END