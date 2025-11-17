MODULE SystemParameters
  LOGICAL  :: PBC
  INTEGER  :: ell, NumTot
  INTEGER,        DIMENSION(:,:), ALLOCATABLE :: iSpin
  DOUBLE PRECISION :: Lambda, gfield, LambdaAmul
END MODULE SystemParameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM Ising
  USE SystemParameters

  IMPLICIT NONE
  LOGICAL      :: ctrlDav
  CHARACTER*80 :: CodeFile
  INTEGER      :: ii, jj, itemp, cnt, Kmax
  DOUBLE PRECISION :: CpuTime, Time
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: Energy, MagnetX, MagnetZ
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: GS1
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: Ham

  OPEN (Unit=1,file='chain.in',status='old')
     READ (1,*) ell       ! Chain length
     READ (1,*) gfield    ! transverse magnetic field
     READ (1,*) Lambda    ! longitudinal field
     READ (1,*) PBC       ! Type of boundary condirions (.true. -> PBC,   .false. -> OBC)
     READ (1,*) ctrlDav   ! Type of diagonalization     (.true. -> Davidson,   .false. -> Lapack full diag)
   CLOSE (unit=1)

!  WRITE (codefile,101) ell, gfield, Lambda
!101 FORMAT('_L',I2.2,'_g',ES7.1,'_lam',ES7.1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Initialization of the relevant spin registers  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! In this part we codify in iSpin each basis element ii = 1, 2^ell in binary notation
  ! For example, suppose that ell = 3. Then NumTot = 8. The various basis elements are:
  ! 1  ->  |000>
  ! 2  ->  |100>
  ! 3  ->  |010>
  ! 4  ->  |110>
  ! 5  ->  |001>
  ! 6  ->  |101>
  ! 7  ->  |011>
  ! 8  ->  |111>
  
  NumTot = 2**ell
  ALLOCATE (iSpin(NumTot,ell))
  iSpin = 0
  DO ii = 1,NumTot
     itemp = ii-1
     DO jj = 1,ell
        iSpin(ii,jj) = Mod(itemp,2)
        itemp = itemp/2
     ENDDO
  ENDDO

!  PRINT *,'Hilbert space dimension', NumTot
!  DO ii = 1,NumTot
!     WRITE (6,102) ii, iSpin(ii,:)
!  ENDDO
!102 FORMAT (I6, 4x, 30(I2))
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  Diagonalization of the INITIAL Hamiltonian(s)  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE (GS1(NumTot), MagnetX(ell), MagnetZ(ell))
  
  IF (.not.(ctrlDav)) THEN      ! Finds the initial g.s. by Exact Diagonalization
     Kmax = NumTot
     ALLOCATE (Energy(NumTot), Ham(NumTot,NumTot))
     
! This subroutine builds up the full Hamiltonian of the model (2^L x 2^L)
     CALL BuildHam(Lambda, gfield, Ham)
     
! This subroutine calls the required LAPACK routines to diagonalize a complex Hermitian matrix
! see:  http://www.netlib.org/lapack/explore-3.1.1-html/zheev.f.html
     CALL Diagonalization(Ham, NumTot, Energy)

! This stores the ground-state wave-function in the complex variable GS1     
     GS1 = Ham(:,1)
     
  ELSEIF (CtrlDav) THEN         ! Finds the initial g.s. by Davidson
     Kmax = 3
     ALLOCATE (Energy(Kmax))
     LambdaAmul = Lambda

! This routine calls the DAVIDSON package to find the ground state of a large sparse Hermitian matrix
! without constructing the full Hamiltonian matrix, but giving only the instructions
! to calculate:  |PsiOut> = Ham * |PsiIn>    [see subroutine AMUL]
! see: http://www.staff.science.uu.nl/~sleij101/JD_software/jd.html
     CALL DAVIDSON(NumTot, Kmax, Energy, GS1)
  ENDIF

  CALL Magnet(GS1, MagnetX, MagnetZ)
  
! Prints out the lowest three energy levels
  PRINT *,'First three energy levels:   ', Energy(1:3)

! Prints out the magnetizations along X and Z
  IF (PBC) THEN
     PRINT *,'GS magnetization along X:   ', MagnetX(1)
     PRINT *,'GS magnetization along Z:   ', MagnetZ(1), MagnetZ(1)/Lambda
  ELSE
     WRITE (6,111)  MagnetX
     WRITE (6,112)  MagnetZ     
111  FORMAT ('GS magnetization profile along X:   ', 100(ES15.8,3x))
112  FORMAT ('GS magnetization profile along Z:   ', 100(ES15.8,3x))
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL cpu_time(cputime)
  PRINT *,'CPU time for the full process:',CPUtime


  STOP
END PROGRAM ISING

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BuildHam(hh, gg, HamOut)
  USE SystemParameters

  IMPLICIT NONE
  INTEGER          :: iHam, Exc, ii
  DOUBLE PRECISION :: hh, gg
  DOUBLE COMPLEX, DIMENSION(NumTot,NumTot) :: HamOut

  HamOut = (0.d0,0.d0)

!  Sigma_z Sigma_z [coupling]
  DO iHam=1,ell-1
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == iSpin(ii,iHam+1))    HamOut(ii,ii) = HamOut(ii,ii) - 1.d0
        IF (iSpin(ii,iHam) /= iSpin(ii,iHam+1))    HamOut(ii,ii) = HamOut(ii,ii) + 1.d0
     ENDDO
  ENDDO
  IF (PBC) THEN   ! This implements periodic boundary conditions
     DO ii=1,NumTot
        IF (iSpin(ii,ell) == iSpin(ii,1))    HamOut(ii,ii) = HamOut(ii,ii) - 1.d0
        IF (iSpin(ii,ell) /= iSpin(ii,1))    HamOut(ii,ii) = HamOut(ii,ii) + 1.d0
     ENDDO
  ENDIF

!  Sigma_z  [longitudinal field]
  DO iHam=1,ell
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == 1)   HamOut(ii,ii) = HamOut(ii,ii) - hh
        IF (iSpin(ii,iHam) == 0)   HamOut(ii,ii) = HamOut(ii,ii) + hh
     ENDDO
  ENDDO

!  Sigma_x  [transverse field]
  DO iHam = 1,ell
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == 1)  Exc = ii - 2**(iHam-1)
        IF (iSpin(ii,iHam) == 0)  Exc = ii + 2**(iHam-1)
        HamOut(Exc,ii) = HamOut(Exc,ii) + gg
     ENDDO
  ENDDO

END SUBROUTINE BuildHam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Diagonalization(Matr, Dim, Ener)
  USE SystemParameters

  IMPLICIT NONE
  INTEGER :: Dim, LWork, LRWork, Info
  DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: Work
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWork
  DOUBLE PRECISION, DIMENSION(Dim)     :: Ener
  DOUBLE COMPLEX,   DIMENSION(Dim,Dim) :: Matr

  LWork  = 22*(2*Dim-1)
  LRWork = 3*Dim-2
  ALLOCATE (Work(Lwork), RWork(LRWork))
  Work   = 0.d0
  RWork  = 0

  CALL ZHEEV( 'V', 'U', Dim, Matr, Dim, Ener, Work, LWork, RWork, Info )
  IF (Info /= 0)  STOP 'Problems in diagonalization'
!  PRINT *,'Workspace dimension:',work(1),lwork
  DEALLOCATE (Work, RWork)

END SUBROUTINE Diagonalization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Magnet(Psi, MagX, MagZ)
  USE SystemParameters

  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(ell)    :: MagX, MagZ
  DOUBLE COMPLEX,   DIMENSION(NumTot) :: Psi
  DOUBLE COMPLEX   :: Mx_sum
  DOUBLE PRECISION :: MagnetZ
  INTEGER          :: iSite, ii, Exc, Mag_ii
  
  MagX = 0.d0
  MagZ = 0.d0
  DO iSite = 1,ell
     Mx_sum = 0.d0
     DO ii = 1,NumTot
        IF (iSpin(ii,iSite) == 1)   MagZ(iSite) = MagZ(iSite) + Abs(Psi(ii))**2
        IF (iSpin(ii,iSite) == 0)   MagZ(iSite) = MagZ(iSite) - Abs(Psi(ii))**2

        IF (iSpin(ii,iSite) == 1)   Exc = ii -2**(iSite-1)
        IF (iSpin(ii,iSite) == 0)   Exc = ii +2**(iSite-1)
        Mx_sum = Mx_sum + dConjg(Psi(ii)) * Psi(Exc)
     ENDDO
     IF (Abs(aImag(Mx_sum)) > 1.d-10)  STOP 'Non real magnetization'
     MagX(iSite) = dReal(Mx_sum)
  ENDDO

  DO ii = 1,NumTot
     Mag_ii = 0
     DO iSite = 1,ell
        IF (iSpin(ii,iSite) == 1)   Mag_ii = Mag_ii + 1
        IF (iSpin(ii,iSite) == 0)   Mag_ii = Mag_ii - 1
     ENDDO
     MagnetZ = MagnetZ + Abs(1.d0*Mag_ii) * Abs(Psi(ii))**2
  ENDDO
  PRINT *,'Symmetry-broken magnetization along Z:   ', MagnetZ/ell
  
END SUBROUTINE Magnet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Davidson(n, Kmax, Eigenvalue, EigenVector)
  USE SystemParameters

  IMPLICIT NONE
  LOGICAL           UseGuess,Wanted
  INTEGER           kmax,jmax,jmin,maxstep,method,m,l,maxnmv,order,testspace,j,lwork,istate,ii,n,kmaxuser
  DOUBLE PRECISION  tol,lock,targetEn,Norm,emin,etemp
  DOUBLE PRECISION, DIMENSION(Kmax)             :: EigenValue
  DOUBLE COMPLEX,   DIMENSION(n)                :: EigenVector
  DOUBLE COMPLEX,   DIMENSION(:),   ALLOCATABLE :: alpha,beta,tmp,residu
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: eivec,zwork


  !!  INIZIALIZATION OF PARAMETERS  !!
  Useguess = .false.
  KMaxUser = KMax
  targetEn = -5.d0*ell
  tol = 1.d-9    ! Tolerance of the eigensolutions: $\Vert \beta H_{SB} x - \alpha x \vert$
  maxnmv = 100    ! Maximum number of matvecs in cgstab or gmres (very problem dependent; typically choose in [5-100])
  wanted = .true. ! If wanted=.true. then computes the converged eigenvectors
  order = -1      ! Selection criterion for Ritz values:  0 (nearest to target);  -1 (smallest real part)
  IF(order == 0)  testspace = 3 ! put 3 if a reasonable value for target is known, else take 2
  IF(order /= 0)  testspace = 2

  IF (3*KmaxUser <= 20) jmax=20          ! Maximum size of the search space:
  IF (3*KmaxUser >  20) jmax=3*KmaxUser
  jmin=2*KmaxUser                        ! Minimum size of the search space
  maxstep = 1000                         ! Maximum number of Jacobi-Davidson iterations
  lock = 1.d-12                          ! Tracking parameter
  method = 2                             ! Method for the linear equation solver  1: gmres(m)  2: cgstab(l)
  m = 30                                 ! Maximum dimension of searchspace for gmres(m):
  l= 2                                   ! Degree of gmres-polynomial in bi-cgstab(l):
  IF (method == 1) lwork =  4 +  m  + 5*jmax + 3*KmaxUser  ! Size of workspace
  IF (method == 2) lwork = 10 + 6*l + 5*jmax + 3*KmaxUser  !KmaxUser is used since Kmax = 1 gives problems ...!
  !!  END OF INIZIALIZATION  !!

  ALLOCATE (alpha(jmax), beta(jmax), eivec(n,Kmax))
  Alpha=0.d0
  Beta=0.d0
  EiVec=0.d0
  ALLOCATE (tmp(n), residu(n), zwork(n,lwork))
  tmp=0.d0
  residu=0.d0
  zwork=0.d0

  CALL JDQZ(ALPHA, BETA, EIVEC, wanted, n, targetEn, tol, Kmax, jmax, jmin, method, m, l, maxnmv, maxstep, &
            lock, order, testspace, zwork, lwork, UseGuess )

  !     Computes the norms of the residuals:
  DO j = 1, Kmax
     CALL AMUL  ( n, eivec(1,j), residu )
     CALL ZSCAL ( n, beta(j), residu, 1 )
     CALL BMUL  ( n, eivec(1,j), tmp )
     CALL ZAXPY ( n, -alpha(j), tmp, 1, residu, 1 )
  ENDDO
  DEALLOCATE (zwork,tmp,residu)
  Eigenvalue(1:Kmax) = dReal(alpha(1:Kmax)/beta(1:Kmax))

  !     Calculates the smallest eigenvalue (ground state)
  emin=eigenvalue(1)
  istate = 1
  DO ii=2,Kmax
     IF (eigenvalue(ii) < emin) THEN
        emin=eigenvalue(ii)
        istate = ii
     ENDIF
  ENDDO
  IF (istate /= 1) THEN
     etemp=eigenvalue(1)
     eigenvalue(1)=eigenvalue(istate)
     eigenvalue(istate)=etemp
  ENDIF
  DEALLOCATE (alpha,beta)

!  print *,'istate',istate
!  Chooses the eigenvector corresponding to the selected eigenvalue
  EigenVector = eivec(:,istate)
  Norm = Sum(dConjg(EigenVector)*(EigenVector))
  EigenVector = EigenVector/(Norm**0.5d0)
  DEALLOCATE (eivec)

END SUBROUTINE Davidson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AMUL(n, PsiIn, PsiOut)
  USE SystemParameters

  IMPLICIT NONE
  INTEGER  :: n, iHam, Exc, ii
  DOUBLE COMPLEX, DIMENSION(NumTot), INTENT(IN)  :: PsiIn
  DOUBLE COMPLEX, DIMENSION(NumTot), INTENT(OUT) :: PsiOut

  PsiOut = (0.d0,0.d0)

!  Sigma_z Sigma_z [coupling]
  DO iHam=1,ell-1
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == iSpin(ii,iHam+1))    PsiOut(ii) = PsiOut(ii) - PsiIn(ii)
        IF (iSpin(ii,iHam) /= iSpin(ii,iHam+1))    PsiOut(ii) = PsiOut(ii) + PsiIn(ii)
     ENDDO
  ENDDO
  
  IF (PBC) THEN   ! This implements periodic boundary conditions
     DO ii=1,NumTot
        IF (iSpin(ii,ell) == iSpin(ii,1))    PsiOut(ii) = PsiOut(ii) - PsiIn(ii)
        IF (iSpin(ii,ell) /= iSpin(ii,1))    PsiOut(ii) = PsiOut(ii) + PsiIn(ii)
     ENDDO
  ENDIF
  
!  Sigma_z  [longitudinal field]
  DO iHam=1,ell
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == 1)   PsiOut(ii) = PsiOut(ii) - LambdaAmul*PsiIn(ii)
        IF (iSpin(ii,iHam) == 0)   PsiOut(ii) = PsiOut(ii) + LambdaAmul*PsiIn(ii)
     ENDDO
  ENDDO

!  Sigma_x  [transverse field]
  Do iHam = 1,ell
     DO ii=1,NumTot
        IF (iSpin(ii,iHam) == 1)   Exc = ii -2**(iHam-1)
        IF (iSpin(ii,iHam) == 0)   Exc = ii +2**(iHam-1)
        PsiOut(Exc) = PsiOut(Exc) + gfield*PsiIn(ii)
     ENDDO
  ENDDO

END SUBROUTINE AMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE BMUL( neq, q, r ) !  Auxiliary subroutine needed by Davidson
  IMPLICIT NONE
  INTEGER :: neq
  DOUBLE COMPLEX :: q(neq),r(neq)
    
  r=q

END SUBROUTINE BMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Guess(N,X)
  IMPLICIT NONE
  INTEGER :: n
  DOUBLE COMPLEX :: X( * )

  X(1:n) = X(1:n)!PsiGuess(1:n) ! NOT IMPLEMENTED

END SUBROUTINE Guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PRECON (neq,psi)  ! NOT IMPLEMENTED
  !     This subroutine computes $\psi = K \psi$, where $K$ is a matrix which mimics the action of 
  !     $(H - \tau \mathbb{1})^{-1}$. Here H is the Hamiltonian to be diagonalized, $\tau$ is the target 
  !     of the Davidson method, namely the value near which the eigenvalues are sought.
  !     A zeroth order approximation is typically used: $K_{i,j} = \delta_{i,j} (H_{i,i}-\tau)^{-1}$
  IMPLICIT NONE
  INTEGER :: neq
  DOUBLE COMPLEX :: psi(neq)

  psi=psi

END SUBROUTINE PRECON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

