program mie

implicit none

integer, parameter :: N = 837, Nsizes = 30

real(8), dimension(N) :: n_list, k_list, l_list
real(8), dimension(Nsizes,N) :: Q_list, x_list
real(8), dimension(Nsizes) :: sizes
complex(8) :: ri

real(8) :: qext, QSCA, qabs, qbk, alb, g
real(8) :: Qpr

real(8), parameter :: PI = 4d0 * datan(1d0)

integer :: i, j


     

qext=0d0
qsca=0d0
qabs=0d0
qbk =0d0
qpr =0d0
alb =0d0
g   =0d0





forall (i=1:Nsizes) sizes(i) = 1.5**(12-i) ! mkm

open(100,file='silic.dat',status='old')
do i = 1, N
    read(100,*) l_list(i), n_list(i), k_list(i)
end do
close(100)

k_list = -k_list

! l_list = l_list * 1e4 ! потому что в файле сантиметры, а не мкм

forall (i=1:N) x_list(:,i) = 2.d0 * 4 * datan(1.d0) * sizes / l_list(i)

! do i = 1, N
    ! do j = 1, Nsizes
        ! write(*,*) x_list(j,i), sizes(j), l_list(i)
    ! end do
! end do


do i = 1, N
    write(*,*) 'ITER:   ', i
    ri = cmplx(n_list(i),k_list(i))
    do j = 1, Nsizes
        if (x_list(j,i) > 1d-6) then
            if (x_list(j,i) < 1d4) then
                call shexq(ri,x_list(j,i),QEXT,QSCA,qabs,qbk,Qpr,alb,g)
            else
                Qpr = -99999.d0
            end if
        else 
            Qpr = 0.d0
        end if 
        Q_list(j,i) = Qpr
    end do
end do


open(555,file='my_out.dat')
do i = 1, N
    write(555,*) l_list(i), Q_list(:,i)
end do
close(555)

do i = 1, Nsizes
    write(*,*) sizes(i)
end do


end program mie


      SUBROUTINE shexq(RI,X,QEXT,QSCA,qabs,qbk,qpr,alb,g)
      parameter(nterms=900000)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)                      
      DIMENSION RU(2*nterms), BESJ(nterms), BESY(nterms),RA(nterms), RB(nterms)
      AX=1.0D0/X                                                        
      NUM=NM(X)                                                         
      NUM2=NM(REAL(cDSQRT(RI*DCONJG(RI)))*X)
        if(num.gt.nterms) then
        write(*,*) 'nterms, num=', nterms, num
        qpr = -99999.0d0
        RETURN
        end if
        if(num2.gt.2*nterms) then
        write(*,*) '2*nterms, num2=', 2*nterms, num2
        qpr = -99999.0d0
        RETURN
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





      FUNCTION NM(X)
      real*8 x
      IF(X.LT.1) GO TO 11
      IF(X.GT.100) GO TO 12
      IF(X.GT.100d0.and.X.lt.50000d0) GO TO 12
      IF(X.Ge.50000d0) GO TO 13
      NM=1.25*X+15.5
      RETURN
   11 NM=7.5*X+9.0 
      RETURN
   12 NM=1.0625*X+28.5
      RETURN
   13 NM=1.05*X+50.5
      RETURN
      END









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
      s=s+i*(i+2d0)/(i+1d0)*(ra(i)*dconjg(ra(i+1))+rb(i)*dconjg(rb(i+1)))+n/i/(i+1d0)*(ra(i)*dconjg(rb(i)))    
      C=C+N*(RA(I)+RB(I))
   11 D=D+N*(RA(I)*DCONJG(RA(I))+RB(I)*DCONJG(RB(I)))
      QEXT=B*C
      QSCA=B*D                          
      qbk=2d0*b*r*dconjg(r)
      qpr=qext-2d0*b*s
      RETURN
      END


