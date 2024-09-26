PROGRAM Trapezoidal
  IMPLICIT NONE
  REAL*8, ALLOCATABLE :: f(:)
  REAL*8, ALLOCATABLE :: x(:)
  REAL*8, ALLOCATABLE :: df(:)
  REAL*8 :: error, c, h
  REAL*8, PARAMETER :: a1=0.5d0, a2=0.5d0

  INTEGER :: i, n, ierr, j, m
  REAL*8  :: intf,dl,lambda

!  WRITE(*,*)"Trapezoidal Rule (for equal grid size)"
  OPEN(1,FILE="int_dhdl.dat",STATUS="old",FORM="formatted")
  OPEN(10,FILE="deltaf.dat",STATUS="replace")
!  OPEN(20,FILE="error.dat",STATUS="replace")
  i=0
  DO 
    READ(1,*,IOSTAT=ierr) 
    IF(ierr/=0)EXIT
    i=i+1
  END DO
  n=i  !Get number of points

  !PRINT *, " Number of points = ", n


  ALLOCATE(x(n))
  ALLOCATE(f(n))
  ALLOCATE(df(n))  
  REWIND(1)

 ! read values of the function
  DO i=1,n
    READ(1,*)x(i), f(i)!, df(i)
  END DO
!
! Loop for trapezoidal integration
 
  lambda=0.d0
  m=1
  DO j=1,n-1
     m=m+1   
     intf=0.d0
     DO i=1,m-1
       dl=x(i+1)-x(i)
       intf=intf+a1*(f(i)+f(i+1))*dl
     END DO
     lambda=lambda+dl
     WRITE(10,*) lambda, " ", intf
  END DO
 !Loop for estimating the overall error
  m=0
  DO j=1,n
     m=m+1
     error=0.d0
     DO i=1,m
       IF(i.eq.1)THEN
         c=x(i+1)-x(i)
       ELSE IF(i.eq.m) THEN
         c=x(i)-x(i-1)
       ELSE
         c=x(i+1)-x(i-1)
       END IF
       error=error + (a1**2)*(c**2)*((df(i)**2))
     END DO
     error=sqrt(error)
     !WRITE(20,*), error
  END DO
  !PRINT*, "Statistical error =", error   

  DEALLOCATE(f)
  DEALLOCATE(x)
  DEALLOCATE(df)
  PRINT *, "XXXX", intf

END PROGRAM Trapezoidal
