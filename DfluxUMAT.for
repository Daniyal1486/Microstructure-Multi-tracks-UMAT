      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
      REAL thermcon,alpha,vel,R0,R,pi,dist
     & x,y,z,power,lamda,abs,AI,M,Eff,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,u
      power=300.0
      vel=0.5
C ********************
      Eff=2. ! 11 this added t othis code
      pi=3.141592654
      abs=0.41
      R0=0.00100
      u=0.01
      t1=u/vel
      t2=t1+u/vel
      t3=t2+u/vel
      t4=t3+u/vel
      t5=t4+u/vel
      t6=t5+u/vel
      t7=t6+u/vel
      t8=t7+u/vel                                             ! four passes on each section :
      t9=t8+u/vel
      t10=t9+u/vel
      AI=(Eff*power)/(pi*(R0*R0)) ! 11
      IF (Time(1) .GT. 0 .AND. TIME(1) .LE. t1) THEN          ! first pass
      dist=vel*Time(1)
      x=COORDS(1)-dist
      y=COORDS(3)
      z=COORDS(2)-0.0006
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t1 .AND. TIME(1) .LE. t2) THEN    ! second pass
      dist=vel*t1+vel*(Time(1)-t1)
      x=COORDS(1)+0.01-dist          ! 0.01 = u
      y=COORDS(3)
      z=COORDS(2)-0.0012
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t2 .AND. TIME(1) .LE. t3) THEN   ! third pass
      dist=vel*t1+vel*(t2-t1)+vel*(Time(1)-t2)
      x=COORDS(1)+0.02-dist
      y=COORDS(3)
      z=COORDS(2)-0.0018
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t3 .AND. TIME(1) .LE. t4) THEN   ! fourth pass          
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(Time(1)-t3)
      x=COORDS(1)+0.03-dist
      y=COORDS(3)
      z=COORDS(2)-0.0024
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t4 .AND. TIME(1) .LE. t5) THEN   ! fifth pass           New beginning point
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(Time(1)-t4)
      x=COORDS(1)+0.03-dist
      y=COORDS(3)
      z=COORDS(2)-0.0024
      R=sqrt(x*x+z*z)
      power=450.0                                                  ! new process parameters
      vel=0.5
      Eff=2. ! 11 this added t othis code
      abs=0.81
      R0=0.00100
      AI=(Eff*power)/(pi*(R0*R0)) ! 11
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t5 .AND. TIME(1) .LE. t6) THEN   ! sixth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(Time(1)-t5)
      x=COORDS(1)+0.04-dist
      y=COORDS(3)
      z=COORDS(2)-0.0006
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t6 .AND. TIME(1) .LE. t7) THEN   ! seventh pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(Time(1)-t6)
      x=COORDS(1)+0.05-dist
      y=COORDS(3)
      z=COORDS(2)-0.0012
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t7 .AND. TIME(1) .LE. t8) THEN   ! eighth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(Time(1)-t7)
      x=COORDS(1)+0.06-dist
      y=COORDS(3)
      z=COORDS(2)-0.0018
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))

      FLUX(2)=0.
      JLTYP=0. 
      END IF
      IF (Time(1) .GT. t8) THEN
      FLUX(1)=0
      FLUX(2)=0.
      JLTYP=0.
      END IF
      RETURN
      END
c----------------------------------------------   
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),D(3,3),
     4 Et(6),DEt(6),Etr(6),DSTATEV(1)

      

      double precision E1, E2, v           ! E1 = EAustenite            ! E2 = EMartensite
      double precision  Denominator, slopeA, slopeM
      double precision pi, Tinit, Af, As, Mf, Ms, aaust, amart

c  initialzing :      User Material  (Mechanical Constants in property Module)
      E1=props(1)              ! E1 = EAustenite
      E2=props(2)              ! E2 = EMartensite
      v=props(3)
      
      
      pi=3.14
      Tinit=298.00
      Af=548.15
      As=496.15
      Ms=457.15
      Mf=348.15
            
      AA=(E1*(1-v))/(1-v-2*(v**2))
      AM=(E2*(1-v))/(1-v-2*(v**2))
      BA=(E1*v)/(1-v-2*(v**2))
      BM=(E2*v)/(1-v-2*(v**2))
      CA=E1/(2*(1+v))
      CM=E2/(2*(1+v))
      
c  calculation of recoverable transformation strain :                           1
      Etr(1)=0.03
      Etr(2)=0.03
      Etr(3)=0.03
      Etr(4)=0.
      Etr(5)=0.
      Etr(6)=0.
      
      
      DSTATEV(1)=0.005
      
       DO K2=1, NTENS
        Etr(K2)=Etr(K2)
       END DO
      
       
c  calculation of recoverable transformation strain :      
      DO K1=1, NTENS
      DO K2=1, NTENS
       DEt(K2)=DSTATEV(1)*Etr(K2)
      END DO
      END DO
      
      DO K1=1, NTENS
      DO K2=1, NTENS
       Et(K2)=Et(K2)+DEt(K2)
      END DO
      END DO
      
      
      
c  calculation of DDSDDE : 
c  most of the elements of DDSDDE are zero . so first we generate a zero matrix and then edit it .   
      do k1=1, NTENS
      do k2=1, NTENS
         DDSDDE(k1,k2)=0
      end do
      end do
c  the matrix relating the stresses and strains is symmetric
      
      
      
c     editing nonzero elements of DDSDDE :
      DDSDDE(1,1)=AA+statev(1)*(AM-AA)
      DDSDDE(1,2)=BA+statev(1)*(BM-BA)
      DDSDDE(1,3)=BA+statev(1)*(BM-BA)
      DDSDDE(2,1)=BA+statev(1)*(BM-BA)
      DDSDDE(2,2)=AA+statev(1)*(AM-AA)
      DDSDDE(2,3)=BA+statev(1)*(BM-BA)
      DDSDDE(3,1)=BA+statev(1)*(BM-BA)
      DDSDDE(3,2)=BA+statev(1)*(BM-BA)
      DDSDDE(3,3)=AA+statev(1)*(AM-AA)
      DDSDDE(4,4)=CA+statev(1)*(CM-CA)
      DDSDDE(5,5)=CA+statev(1)*(CM-CA)
      DDSDDE(6,6)=CA+statev(1)*(CM-CA)
      

c  calculation of stress :
      DO K1=1, NTENS
      DO K2=1, NTENS
    
       STRESS(K1)=STRESS(K1)+DDSDDE(K1, K2)*DSTRAN(K2)        
      END DO
      END DO
      
      
      slopeA=8000000.0
      slopeM=8000000.0
      aaust=pi*(Af-As)
      amart=pi*(Ms-Mf)
      xx=(aaust)*(TEMP-As-(STRESS(K1)/slopeA))
      yy=(amart)*(TEMP-Mf-(STRESS(K1)/slopeM))
      MVF0=1.0
      MVF1=0.0
      
c  using state variables :
c  for the use of SDV's , the Martensite Volume Fraction (statev(1)) is coded .
c  we don't edit the calculated stresses and DDSDDE .

c every statev can be seen in visualization module as contour plots :      


      IF (TEMP .GT. 298 .AND. TEMP .LT. 496) THEN
        statev(1)=(0.50)*(COS(yy))+0.50
      ELSE IF (TEMP .GE. 496) THEN
        statev(1)=(0.50)*(COS(xx)+1)
   
      END IF
     
      

      RETURN
      END