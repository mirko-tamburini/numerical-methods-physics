PROGRAM RungeKutta

 IMPLICIT NONE
 INTEGER var,N,i
 REAL xi,xf
 REAL, DIMENSION(:), ALLOCATABLE :: y0
 REAL, DIMENSION(:,:), ALLOCATABLE :: z
 REAL, DIMENSION(:), ALLOCATABLE :: t
 EXTERNAL derivs

 OPEN (Unit=2,file='Results',status='unknown')
 var=2                ! number of variables
 N=10                 ! number of step
 xi=0.d0
 xf=2.d0
 
 ALLOCATE (t(N+1))
 ALLOCATE (z(var,N+1))
 ALLOCATE (y0(var))
 y0(1)=1.d0
 y0(2)=1.d0

 CALL rkdumb(y0,var,xi,xf,N,derivs,t,z)
 
 DO i=1,(N+1)
    WRITE(2,*)t(i),z(1,i),z(2,i)
 ENDDO

 CLOSE(unit=2)

 STOP
END PROGRAM RungeKutta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs,xx,y)
INTEGER nstep,nvar,NMAX,NSTPMX
PARAMETER (NMAX=2,NSTPMX=200) !Maximum number of functions and 
!maximum number of values to be stored. (usually set NMAX=50 and NSTPMX=200)

REAL x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
EXTERNAL derivs
!COMMON /path/ xx,y          ! Storage of results.
! USES rk4
!Starting from initial values vstart(1:nvar) known at x1 use fourth-order Runge-Kutta to
!advance nstep equal increments to x2. The user-supplied subroutine derivs(x,v,dvdx)
!evaluates derivatives. Results are stored in the common block path. Be sure to dimension
!the common block appropriately.
INTEGER i,k
REAL h,x,dv(NMAX),v(NMAX)
DO i=1,nvar              ! Load starting values.
    v(i)=vstart(i)
    y(i,1)=v(i)
ENDDO
xx(1)=x1
x=x1
h=(x2-x1)/nstep
DO k=1,nstep              !Take nstep steps.
    CALL derivs(x,v,dv)
    CALL rk4(v,dv,nvar,x,h,v,derivs)
    IF((x+h)==x) STOP 'stepsize not significant in rkdumb'
    x=x+h
    xx(k+1)=x               ! Store intermediate steps.
    DO i=1,nvar
        y(i,k+1)=v(i)  
    ENDDO
!    PRINT*,y(1,k)
ENDDO
RETURN
END SUBROUTINE rkdumb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
INTEGER n,NMAX
REAL h,x,dydx(n),y(n),yout(n)
EXTERNAL derivs
PARAMETER (NMAX=2) !Set to the maximum number of functions.
!Given values for the variables y(1:n) and their derivatives dydx(1:n) known at x, use
!the fourth-order Runge-Kutta method to advance the solution over an interval h and return
!the incremented variables as yout(1:n), which need not be a distinct array from y. The
!user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
INTEGER i
REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
hh=h*0.5
h6=h/6.
xh=x+hh
DO i=1,n                     !First step.
    yt(i)=y(i)+hh*dydx(i)
ENDDO
CALL derivs(xh,yt,dyt)          !Second step.
DO i=1,n
    yt(i)=y(i)+hh*dyt(i)
ENDDO
CALL derivs(xh,yt,dym)          !Third step.
DO i=1,n
    yt(i)=y(i)+h*dym(i)
    dym(i)=dyt(i)+dym(i)
ENDDO
CALL derivs(x+h,yt,dyt)         !Fourth step.
DO i=1,n                     !Accumulate increments with proper weights.
    yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
ENDDO
RETURN
END SUBROUTINE rk4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION derivs(x,y,dydx)
REAL x,y(2),dydx(2)
dydx(1)=y(1)
dydx(2)=1.d0
RETURN
END FUNCTION derivs