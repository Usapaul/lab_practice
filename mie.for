c*********** Light Scattering by Spherical Particles *************
c                                                                *
c   Calculations of extinction, scattering, absorption, etc      *
c  efficiency factors for homogeneous spheres (the Mie theory).  *
c................................................................*
c  Input data:          filename:   mie.dat                      *
c   ri = n-k*i: complex index of refraction                      *
c           Nx: number of sizes                                  *
c           x1: size 1                                           *
c           x2: size 2                                           *
c          ...: ...                                              *
c................................................................*
c  Output data:         filename:   mie.out                      *
c   Qext : extinction factors                                    *
c   Qsca : scattering factors                                    *
c   Qabs : absorption factors                                    *
c   Qbk  : backscattering factors                                *
c   Qpr  : radiation pressure factors                            *
c  albedo: particle's albedo                                     *
c   g    : asymmetry factor                                      *
c................................................................*
c NB! In order to treat very large particles,                    *
c     one needs to enlarge the parameter NTERMS.                 *
c................................................................*
c created by N.V. Voshchinnikov                                  *
c with a support of the Volkswagen Foundation (Germany)          *
c (c) 1989/98 Astronomical Institute, St.Petersburg University   *
c*****************************************************************
c
      parameter(nterms=300)
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMPLEX*16 ri
      DIMENSION xx(50)
   20 FORMAT(2D14.10)                                                 
   21 FORMAT(3X,'SPHERES: homogeneous')                                 
   22 FORMAT(3X,'THEORY:  exact')                                        
  223 FORMAT(I4)                                          
   24 FORMAT(4X,'x',7X,'Qext',4X,'Qsca',4X,'Qabs',4X,'Qbk',5x,'Qpr',
     *       4X,'Albedo',3x,'g')
   25 FORMAT(13X,F8.4,9X,I3,6X,2D15.4,F10.2,5X,D12.3,F11.3)             
   29 FORMAT(I4/(F8.4))                                                 
   30 FORMAT(3x,'m = '2F14.10,'*i')
   31 FORMAT(F8.4)                                                      
   32 FORMAT(1H  )                                                      
   33 FORMAT(1X,8(F8.4))
  233 FORMAT(4X,i4,2F10.3,f11.1)                                             
   34 FORMAT( 8X,'x = ',F7.3,5x,'Qext = ',f8.4)          
   36 FORMAT(1X,64 ('='))
   37 FORMAT(1X,60 ('*')//)                                             
   38 FORMAT(1X,10I4)                                                   
   41 FORMAT(1X,64('.'))
   42 FORMAT(1X,64('-'))
   44 FORMAT(1X,A8/(1X,2F8.4))                                          
   48 FORMAT( 5X,F8.4,2X,2(F10.3,3X),3X,F8.2)                           
c
c*                  INPUT
c
        print *,'start Mie'
      open (unit=05,file='mie.dat',status='old',access='sequential')
      open (unit=07,file='mie.out',status='unknown',access='append')
*      PI = 4d0 * datan(1d0)
      READ (5,20) ri
      READ (5,29) nx,(xx(I),I=1,Nx)
c
c-------------------------------------------------------------------
c***              * Spheres *
c-------------------------------------------------------------------
c*                  Efficiency factors
c
      qext=0d0
      qsca=0d0
      qabs=0d0
      qbk =0d0
      qpr =0d0
      alb =0d0
      g   =0d0
c
c*                  Exact theory
c
      write (7,22)
      write (7,30) ri
      write (7,41)
      write (7,24)
      write (7,42)
      write (*,42)

        do 276 i = 1, nx
        x=xx(i)
        if(x.le.1d-6) go to 444
        call shexq(ri,x,qext,qsca,qabs,qbk,qpr,alb,g)
  444   continue
        write (7,33) x,qext,qsca,qabs,qbk,qpr,alb,g
        write (*,*) x, qext
  276   continue
c
 1000 continue
      write (7,42)
      write (*,42)
      STOP
      END
c--------------------------------------------------------------------
c **********   shexq - Sphers: homogeneous
c                      Theory: exact
c                      Results: efficiency factors
c--------------------------------------------------------------------
      SUBROUTINE shexq(RI,X,QEXT,QSCA,qabs,qbk,qpr,alb,g)
      parameter(nterms=900000)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)                      
      DIMENSION RU(2*nterms), BESJ(nterms), BESY(nterms),
     *          RA(nterms), RB(nterms)
      AX=1.0D0/X                                                        
      NUM=NM(X)                                                         
      NUM2=NM(cDSQRT(RI*DCONJG(RI))*X)
        if(num.gt.nterms) then
        write(*,*) 'nterms, num=', nterms, num
        pause
        stop
        end if
        if(num2.gt.2*nterms) then
        write(*,*) '2*nterms, num2=', 2*nterms, num2
        pause
        stop
        end if


      CALL AA(AX,RI,NUM2,RU)                                            
      CALL BESSEL(AX,NUM,BESJ,BESY)                                     
      CALL AB(AX,RI,NUM,NUM1,RU,BESJ,BESY,RA,RB)
      CALL QQ1(AX,NUM1,QEXT,QSCA,qbk,qpr,RA,RB)                             
      qabs=qext-qsca
      alb=qsca/qext
      g=(qext-qpr)/qsca
      RETURN                                                            
      END                                                               
c--------------------------------------------------------------------
c NM-auxiliary function for AA & BESSEL
c    (number NM is calculated using X)
c see: Trudy Astronom. Observ. LGU V.28,P.14,1971
c    for X>1 value of NM was raised
c August 1989, AO LGU
c--------------------------------------------------------------------
      FUNCTION NM(X)
      real*8 x
      IF(X.LT.1) GO TO 11
      IF(X.GT.100.and.X.lt.50000) GO TO 12
      IF(X.GT.50000) GO TO 13
      NM=1.25*X+15.5
      RETURN
   11 NM=7.5*X+9.0 
      RETURN
   12 NM=1.0625*X+28.5
      RETURN
   13 NM=1.05*X+50.5
      RETURN
      END
c--------------------------------------------------------------------
c AA-subroutine for calculations of the ratio of the derivative
c    to the function for Bessel functions of half order with
c    the complex argument: J'(N)/J(N).
c    The calculations are given by the recursive expression
c    ``from top to bottom'' beginning from N=NUM.
c    RU-array of results.
c    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
c    RI - complex refractive index.
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE AA(A,RI,NUM,RU)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION RU(NUM)
      S=A/RI
      RU(NUM)=(NUM+1.0D0)*S
      NUM1=NUM-1
      DO 13 J=1,NUM1
      I=NUM-J
      I1=I+1
      S1=I1*S
   13 RU(I)=S1-1.0D0/(RU(I1)+S1)
      RETURN
      END
c--------------------------------------------------------------------
c BESSEL-subroutine for calculations of the Bessel functions
c    of half order and first (J(N)) second (Y(N)) kinds with
c    the real argument X (A=1/X).
c    The calculations are started from N=NUM.
c    Desription of method see:
c    V.M.Loskutov, Trudy Astronom. Observ. LGU V.28,P.14,1971
c    BESJ-array of functions J(N), BESY-array of functions Y(N)
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE BESSEL(A,NUM,BESJ,BESY)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1) 
      BESJ(NUM+1)=0.0D0
      BESJ(NUM)=0.1D-11
      N=2*NUM+1
      NUM1=NUM-1
      DO 11 I=1,NUM1
      N=N-2
      I1=NUM-I
   11 BESJ(I1)=N*A*BESJ(I1+1)-BESJ(I1+2)
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*DSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
   14 BESY(I)=(2.0D0*I2+1.0D0)*A*BESY(I-1)-BESY(I2)
      RETURN
      END
c--------------------------------------------------------------------
c AB-subroutine for calculations of the complex coefficients
c    A(N), B(N) for spherical particles.
c    A=1/X (X=2*PI*A(particle radius)/LAMBDA - size parameter).
c    RI - complex refractive index.
c    The coefficients are calculated up to the number NUM1.LE.NUM,
c    for which |A(N)**2+B(N)**2|.LE.10**(-60)
c    RA-array of coefficients A(N), RB-array of coefficients B(N)
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE AB(A,RI,NUM,NUM1,RU,BESJ,BESY,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION RU(NUM),BESJ(NUM),BESY(NUM),RA(NUM),RB(NUM)
      S3=(0.0D0,1.0D0)
      DO 11 I=1,NUM-1
      N=I+1
      S=RU(I)/RI+I*A
      S1=S*BESJ(N)-BESJ(I)
      S2=S*BESY(N)-BESY(I)
      RA(I)=S1/(S1-S3*S2)
      S=RU(I)*RI+I*A
      S1=S*BESJ(N)-BESJ(I)
      S2=S*BESY(N)-BESY(I)
      RB(I)=S1/(S1-S3*S2)
      P=RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I))
      IF (P.LE.1D-60) GO TO 12
   11 CONTINUE
   12 NUM1=I
      RETURN
      END
c--------------------------------------------------------------------
c QQ1-subroutine for calculations of the efficiency factors for
c     extinction (QEXT), scattering (QSCA), backscattering (QBK)
c     and radiation pressure (QPR) for spherical particles.
c August 1989, AO LGU
c--------------------------------------------------------------------
      SUBROUTINE QQ1(A,NUM,QEXT,QSCA,qbk,qpr,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION RA(NUM),RB(NUM)
      B=2.0D0*A*A
      C=0.0D0
      D=0.0D0
      s=(0d0,0d0)
      r=(0d0,0d0)
      N=1
      DO 11 I=1,NUM-1
      N=N+2      
      r=r+(i+0.5d0)*(-1)**i*(ra(i)-rb(i))      
      s=s+i*(i+2d0)/(i+1d0)*(ra(i)*dconjg(ra(i+1))
     *  +rb(i)*dconjg(rb(i+1)))+n/i/(i+1d0)*(ra(i)*dconjg(rb(i)))    
      C=C+N*(RA(I)+RB(I))
   11 D=D+N*(RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I)))
      QEXT=B*C
      QSCA=B*D                          
      qbk=2d0*b*r*dconjg(r)
      qpr=qext-2d0*b*s
      RETURN
      END
c=== eof ===
