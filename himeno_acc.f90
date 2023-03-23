! Simple f90 port of original f77 code
! This code add OpenACC offload
! Douglas Shanks, HPE , 2023
!
! This benchmark test program is measuring a cpu performance
! of floating point operation by a Poisson equation solver.
!
! If you have any question, please ask me via email.
! written by Ryutaro HIMENO, November 26, 2001.
! Version 3.0
! ----------------------------------------------
! Ryutaro Himeno, Dr. of Eng.
! Head of Computer Information Division,
! RIKEN (The Institute of Pysical and Chemical Research)
! Email : himeno@postman.riken.go.jp
! -----------------------------------------------------------
! You can adjust the size of this benchmark code to fit your target
! computer. In that case, please chose following sets of
! (mimax,mjmax,mkmax):
! small : 65,33,33
! small : 129,65,65
! midium: 257,129,129
! large : 513,257,257
! ext.large: 1025,513,513
! This program is to measure a computer performance in MFLOPS
! by using a kernel which appears in a linear solver of pressure
! Poisson eq. which appears in an incompressible Navier-Stokes solver.
! A point-Jacobi method is employed in this solver as this method can
! be easyly vectrized and be parallelized.
! ------------------
! Finite-difference method, curvilinear coodinate system
! Vectorizable and parallelizable on each grid point
! No. of grid points : imax x jmax x kmax including boundaries
! ------------------
! A,B,C:coefficient matrix, wrk1: source term of Poisson equation
! wrk2 : working area, OMEGA : relaxation parameter
! BND:control variable for boundaries and objects ( = 0 or 1)
! P: pressure
! -------------------

PROGRAM himenobmtxp

IMPLICIT REAL*4(a-h,o-z)

INCLUDE 'mpif.h'
INCLUDE 'param.h'

!     ttarget specifys the measuring period in sec
PARAMETER (ttarget=60.0)

REAL*8  cpu,cpu0,cpu1,xmflops2,flop

omega=0.8
mx= mx0-1
my= my0-1
mz= mz0-1

!C Initializing communicator
CALL initcomm

!C Initializaing computational index
CALL initmax(mx,my,mz,it)

!C Initializing matrixes
CALL initmt(mz,it)
IF(id == 0) THEN
  WRITE(*,*) 'Sequential version array size'
  WRITE(*,*) ' mimax=',mx0,' mjmax=',my0,' mkmax=',mz0
  WRITE(*,*) 'Parallel version  array size'
  WRITE(*,*) ' mimax=',mimax,' mjmax=',mjmax,' mkmax=',mkmax
  WRITE(*,*) ' imax=',imax,' jmax=',jmax,' kmax=',kmax
  WRITE(*,*) ' I-decomp= ',ndx,' J-decomp= ',ndy, ' K-decomp= ',ndz
  WRITE(*,*)
END IF

!C Start measuring

nn=10
IF(id == 0) THEN
  WRITE(*,*) ' Start rehearsal measurement process.'
  WRITE(*,*) ' Measure the performance in 3 times.'
END IF

gosa= 0.0
cpu= 0.0
CALL mpi_barrier(mpi_comm_world,ierr)
cpu0= mpi_wtime()
! Jacobi iteration
CALL jacobi(nn,gosa)
cpu1= mpi_wtime() - cpu0

CALL mpi_allreduce(cpu1, cpu,  &
    1, mpi_real8,  &
    mpi_max, mpi_comm_world,  &
    ierr)

flop=REAL(mx-2)*REAL(my-2)*REAL(mz-2)*34.0
IF(cpu /= 0.0) xmflops2=flop/cpu*1.0D-6*REAL(nn)
IF(id == 0) THEN
  WRITE(*,*) '  MFLOPS:',xmflops2,'  time(s):',cpu,gosa
END IF
!      nn= int(ttarget/(cpu/3.0))
!     Number of iterations forced to 50 for the XL case on a 2x2x2 decomposition
nn=50

!     end the test loop
IF(id == 0) THEN
  WRITE(*,*) 'Now, start the actual measurement process.'
  WRITE(*,*) 'The loop will be excuted in',nn,' times.'
  WRITE(*,*) 'This will take about one minute.'
  WRITE(*,*) 'Wait for a while.'
END IF

gosa= 0.0
cpu= 0.0
CALL mpi_barrier(mpi_comm_world,ierr)
cpu0= mpi_wtime()
! Jacobi iteration
CALL jacobi(nn,gosa)
cpu1= mpi_wtime() - cpu0

CALL mpi_reduce(cpu1, cpu,  &
    1, mpi_real8,  &
    mpi_max, 0,  &
    mpi_comm_world, ierr)

IF(id == 0) THEN
  IF(cpu /= 0.0)  xmflops2=flop*1.0D-6/cpu*REAL(nn)
  
  WRITE(*,*) ' Loop executed for ',nn,' times'
  WRITE(*,*) ' Gosa :',gosa
  WRITE(*,*) ' MFLOPS:',xmflops2, '  time(s):',cpu
  score=xmflops2/82.84
  WRITE(*,*) ' Score based on Pentium III 600MHz :',score
END IF
CALL mpi_finalize(ierr)

STOP
END PROGRAM himenobmtxp


!**************************************************************

SUBROUTINE initmt(mz,it)
!**************************************************************

IMPLICIT REAL*4(a-h,o-z)

INTEGER, INTENT(IN OUT)                  :: mz
INTEGER, INTENT(IN OUT)                  :: it

INCLUDE 'param.h'

DO k=1,mkmax
  DO j=1,mjmax
    DO i=1,mimax
      a(i,j,k,1)=0.0
      a(i,j,k,2)=0.0
      a(i,j,k,3)=0.0
      a(i,j,k,4)=0.0
      b(i,j,k,1)=0.0
      b(i,j,k,2)=0.0
      b(i,j,k,3)=0.0
      c(i,j,k,1)=0.0
      c(i,j,k,2)=0.0
      c(i,j,k,3)=0.0
      p(i,j,k)=0.0
      wrk1(i,j,k)=0.0
      wrk2(i,j,k)=0.0
      bnd(i,j,k)=0.0
    END DO
  END DO
END DO

DO k=1,kmax
  DO j=1,jmax
    DO i=1,imax
      a(i,j,k,1)=1.0
      a(i,j,k,2)=1.0
      a(i,j,k,3)=1.0
      a(i,j,k,4)=1.0/6.0
      b(i,j,k,1)=0.0
      b(i,j,k,2)=0.0
      b(i,j,k,3)=0.0
      c(i,j,k,1)=1.0
      c(i,j,k,2)=1.0
      c(i,j,k,3)=1.0
      p(i,j,k)=FLOAT((k-1+it)*(k-1+it)) /FLOAT((mz-1)*(mz-1))
      wrk1(i,j,k)=0.0
      wrk2(i,j,k)=0.0
      bnd(i,j,k)=1.0
    END DO
  END DO
END DO

RETURN
END SUBROUTINE initmt

!*************************************************************

SUBROUTINE jacobi(nn,gosa)
!*************************************************************

IMPLICIT REAL*4(a-h,o-z)

INTEGER, INTENT(IN)                      :: nn
REAL, INTENT(OUT)                        :: gosa

INTEGER*4 istat

INCLUDE 'mpif.h'
INCLUDE 'param.h'

!$acc data copyin(a,b,c,p,wrk2,wrk1,bnd)
DO loop=1,nn
  gosa=0.0
  wgosa=0.0

!$acc kernels
!$acc parallel loop collapse(3) reduction(+: wgosa) private(s0,ss) 
DO k=2,kmax-1
  DO j=2,jmax-1
    DO i=2,imax-1
      s0=a(i,j,k,1)*p(i+1,j,k)+a(i,j,k,2)*p(i,j+1,k) +a(i,j,k,3)*p(i,j,k+1)  &
          +b(i,j,k,1)*(p(i+1,j+1,k)-p(i+1,j-1,k) -p(i-1,j+1,k)+p(i-1,j-1,k))  &
          +b(i,j,k,2)*(p(i,j+1,k+1)-p(i,j-1,k+1) -p(i,j+1,k-1)+p(i,j-1,k-1))  &
          +b(i,j,k,3)*(p(i+1,j,k+1)-p(i-1,j,k+1) -p(i+1,j,k-1)+p(i-1,j,k-1))  &
          +c(i,j,k,1)*p(i-1,j,k)+c(i,j,k,2)*p(i,j-1,k)  &
          +c(i,j,k,3)*p(i,j,k-1)+wrk1(i,j,k)
      ss=(s0*a(i,j,k,4)-p(i,j,k))*bnd(i,j,k)
      wgosa=wgosa+ss*ss
      wrk2(i,j,k)=p(i,j,k)+omega *ss
    END DO
  END DO
END DO
!$acc end parallel loop

!$acc parallel loop collapse(3) 
DO k=2,kmax-1
  DO j=2,jmax-1
    DO i=2,imax-1
      p(i,j,k)=wrk2(i,j,k)
    END DO
  END DO
END DO
!$acc end parallel loop
!$acc end kernels

CALL sendp(ndx,ndy,ndz)

CALL mpi_allreduce(wgosa, gosa,  &
    1, mpi_real4,  &
    mpi_sum, mpi_comm_world,  &
    ierr)

END DO
!$acc end data
!C End of iteration
RETURN
END SUBROUTINE jacobi



SUBROUTINE initcomm

IMPLICIT REAL*4(a-h,o-z)

INCLUDE 'mpif.h'
INCLUDE 'param.h'

LOGICAL :: ipd(3),ir
DIMENSION  idm(3)

CALL mpi_init(ierr)
CALL mpi_comm_size(mpi_comm_world,npe,ierr)
CALL mpi_comm_rank(mpi_comm_world,id,ierr)

IF(ndx*ndy*ndz /= npe) THEN
  IF(id == 0) THEN
    WRITE(*,*) 'Invalid number of PE'
    WRITE(*,*) 'Please check partitioning pattern'
    WRITE(*,*) '                 or number of  PE'
  END IF
  CALL mpi_finalize(ierr)
  STOP
END IF

icomm= mpi_comm_world

idm(1)= ndx
idm(2)= ndy
idm(3)= ndz

ipd(1)= .false.
ipd(2)= .false.
ipd(3)= .false.
ir= .false.

CALL mpi_cart_create(icomm, ndims,  &
    idm, ipd,  &
    ir, mpi_comm_cart,  &
    ierr)
CALL mpi_cart_get(mpi_comm_cart, ndims,  &
    idm, ipd,  &
    iop, ierr)


IF(ndz > 1) THEN
  CALL mpi_cart_shift(mpi_comm_cart, 2,  &
      1, npz(1),  &
      npz(2), ierr)
END IF

IF(ndy > 1) THEN
  CALL mpi_cart_shift(mpi_comm_cart, 1,  &
      1, npy(1),  &
      npy(2), ierr)
END IF

IF(ndx > 1) THEN
  CALL mpi_cart_shift(mpi_comm_cart, 0,  &
      1, npx(1),  &
      npx(2), ierr)
END IF

RETURN
END SUBROUTINE initcomm



SUBROUTINE initmax(mx,my,mz,ks)

IMPLICIT REAL*4(a-h,o-z)

INTEGER, INTENT(IN)                      :: mx
INTEGER, INTENT(IN)                      :: my
INTEGER, INTENT(IN)                      :: mz
INTEGER, INTENT(OUT)                     :: ks

INCLUDE 'param.h'
INCLUDE 'mpif.h'

INTEGER :: itmp
INTEGER :: mx1(0:ndx),my1(0:ndy),mz1(0:ndz)
INTEGER :: mx2(0:ndx),my2(0:ndy),mz2(0:ndz)

!C    define imax, communication direction
itmp= mx/ndx
mx1(0)= 0
DO  i=1,ndx
  IF(i <= MOD(mx,ndx)) THEN
    mx1(i)= mx1(i-1) + itmp + 1
  ELSE
    mx1(i)= mx1(i-1) + itmp
  END IF
END DO
DO i=0,ndx-1
  mx2(i)= mx1(i+1) - mx1(i)
  IF(i /= 0)     mx2(i)= mx2(i) + 1
  IF(i /= ndx-1) mx2(i)= mx2(i) + 1
END DO

itmp= my/ndy
my1(0)= 0
DO  i=1,ndy
  IF(i <= MOD(my,ndy)) THEN
    my1(i)= my1(i-1) + itmp + 1
  ELSE
    my1(i)= my1(i-1) + itmp
  END IF
END DO
DO i=0,ndy-1
  my2(i)= my1(i+1) - my1(i)
  IF(i /= 0)      my2(i)= my2(i) + 1
  IF(i /= ndy-1)  my2(i)= my2(i) + 1
END DO

itmp= mz/ndz
mz1(0)= 0
DO  i=1,ndz
  IF(i <= MOD(mz,ndz)) THEN
    mz1(i)= mz1(i-1) + itmp + 1
  ELSE
    mz1(i)= mz1(i-1) + itmp
  END IF
END DO
DO i=0,ndz-1
  mz2(i)= mz1(i+1) - mz1(i)
  IF(i /= 0)      mz2(i)= mz2(i) + 1
  IF(i /= ndz-1)  mz2(i)= mz2(i) + 1
END DO

imax= mx2(iop(1))
jmax= my2(iop(2))
kmax= mz2(iop(3))

IF(iop(3) == 0) THEN
  ks= mz1(iop(3))
ELSE
  ks= mz1(iop(3)) - 1
END IF

!     j-k plane  divied by i-direction
IF(ndx > 1) THEN
  CALL mpi_type_vector(jmax*kmax, 1,  &
      mimax, mpi_real4,  &
      jkvec, ierr)
  CALL mpi_type_commit(jkvec, ierr)
END IF

!     i-k plane  divied by j-direction
IF(ndy > 1) THEN
  CALL mpi_type_vector(kmax, imax,  &
      mimax*mjmax, mpi_real4,  &
      ikvec, ierr)
  CALL mpi_type_commit(ikvec, ierr)
END IF

!     new vector k-direction
IF(ndz > 1) THEN
  CALL mpi_type_vector(jmax, imax,  &
      mimax, mpi_real4,  &
      ijvec, ierr)
  CALL mpi_type_commit(ijvec, ierr)
END IF

RETURN
END SUBROUTINE initmax



SUBROUTINE sendp(ndx,ndy,ndz)

IMPLICIT REAL*4(a-h,o-z)

INTEGER, INTENT(IN)                  :: ndx
INTEGER, INTENT(IN)                  :: ndy
INTEGER, INTENT(IN)                  :: ndz

IF(ndz > 1) THEN
  CALL sendp3()
END IF

IF(ndy > 1) THEN
  CALL sendp2()
END IF

IF(ndx > 1) THEN
  CALL sendp1()
END IF

RETURN
END SUBROUTINE sendp



SUBROUTINE sendp3()

IMPLICIT REAL*4(a-h,o-z)

INCLUDE 'mpif.h'
INCLUDE 'param.h'

DIMENSION ist(mpi_status_size,0:3),ireq(0:3)
DATA ireq /4*mpi_request_null/

CALL mpi_irecv(p(1,1,kmax), 1,  &
    ijvec, npz(2),  &
    1, mpi_comm_cart,  &
    ireq(3), ierr)

CALL mpi_irecv(p(1,1,1), 1,  &
    ijvec, npz(1),  &
    2, mpi_comm_cart,  &
    ireq(2), ierr)

CALL mpi_isend(p(1,1,2), 1,  &
    ijvec, npz(1),  &
    1, mpi_comm_cart,  &
    ireq(0), ierr)

CALL mpi_isend(p(1,1,kmax-1), 1,  &
    ijvec, npz(2),  &
    2, mpi_comm_cart,  &
    ireq(1), ierr)

CALL mpi_waitall(4, ireq,  &
    ist, ierr)

RETURN
END SUBROUTINE sendp3



SUBROUTINE sendp2()

IMPLICIT REAL*4(a-h,o-z)

INCLUDE 'mpif.h'
INCLUDE 'param.h'

DIMENSION ist(mpi_status_size,0:3),ireq(0:3)
DATA ireq /4*mpi_request_null/

CALL mpi_irecv(p(1,1,1), 1,  &
    ikvec, npy(1),  &
    2, mpi_comm_cart,  &
    ireq(3), ierr)

CALL mpi_irecv(p(1,jmax,1), 1,  &
    ikvec, npy(2),  &
    1, mpi_comm_cart,  &
    ireq(2), ierr)

CALL mpi_isend(p(1,2,1), 1,  &
    ikvec, npy(1),  &
    1, mpi_comm_cart,  &
    ireq(0), ierr)

CALL mpi_isend(p(1,jmax-1,1), 1,  &
    ikvec, npy(2),  &
    2, mpi_comm_cart,  &
    ireq(1), ierr)

CALL mpi_waitall(4, ireq,  &
    ist, ierr)

RETURN
END SUBROUTINE sendp2



SUBROUTINE sendp1()

IMPLICIT REAL*4(a-h,o-z)

INCLUDE 'mpif.h'
INCLUDE 'param.h'

DIMENSION ist(mpi_status_size,0:3),ireq(0:3)
DATA ireq /4*mpi_request_null/

CALL mpi_irecv(p(1,1,1), 1,  &
    jkvec, npx(1),  &
    2, mpi_comm_cart,  &
    ireq(3), ierr)

CALL mpi_irecv(p(imax,1,1), 1,  &
    jkvec, npx(2),  &
    1, mpi_comm_cart,  &
    ireq(2), ierr)

CALL mpi_isend(p(2,1,1), 1,  &
    jkvec, npx(1),  &
    1, mpi_comm_cart,  &
    ireq(0), ierr)

CALL mpi_isend(p(imax-1,1,1), 1,  &
    jkvec, npx(2),  &
    2, mpi_comm_cart,  &
    ireq(1), ierr)

CALL mpi_waitall(4, ireq,  &
    ist, ierr)

RETURN
END SUBROUTINE sendp1

