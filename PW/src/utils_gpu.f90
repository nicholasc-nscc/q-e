!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! General-purpose routines for scalar products, printing,
! linear-algebra operators for exact-exchange and localization
!
!----------------------------------------------------------------------
SUBROUTINE matcalc_gpu( label, DoE, PrtMat, ninner, n, m, U, V, mat, ee )
  !------------------------------------------------------------------
  !! Compute the (n,n) matrix representation \(\langle U|V\rangle\)
  !! and its weighted trace (energy) from \(V(m,n)\) and \(U(m,n)\).
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : current_k, wg
  USE gvect,                ONLY : gstart
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: label
  !! it specifies the meaning of the output
  LOGICAL, INTENT(IN) :: DoE
  !! if TRUE calculate the trace
  INTEGER, INTENT(IN) :: PrtMat
  !! printing index
  INTEGER, INTENT(IN) :: ninner
  !! inner dimension in the matrix product
  INTEGER, INTENT(IN) :: n
  !! second dimension of U
  INTEGER, INTENT(IN) :: m
  !! second dimension of V
  COMPLEX(DP), INTENT(IN) :: U(ninner,n)
  !! input - U matrix
  COMPLEX(DP), INTENT(IN) :: V(ninner,m)
  !! input - V matrix
  REAL(DP), INTENT(OUT) :: mat(n,m)
  !! output matrix \(\langle U|V\rangle\)
  REAL(DP), INTENT(OUT) :: ee
  !! the weighted trace (energy) of the product
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: U, V, mat
#endif
  !
  ! ... local variables
  !
  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock_gpu('matcalc')

  string = 'M-'
  mat = 0.0_dp
  CALL MYDGEMM( 'C', 'N', n, m, 2*ninner, 2.0_DP, U, 2*ninner, V, 2*ninner, 0.0_DP, mat, n )
  IF ( gstart == 2 ) CALL MYDGER( n, m, -1.0_DP, U, 2*ninner, V, 2*ninner, mat, n )
  CALL mp_sum( mat( :, 1:m ), intra_bgrp_comm )

  IF( PrtMat > 1 ) CALL errore('matcalc_gpu', 'cannot print matrix', 1)

  IF(DoE) THEN
     IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
     string = 'E-'
     ee = 0.0_dp
     !$acc parallel loop reduction(+:ee) copyin(wg)
     DO i = 1,n
        ee = ee + wg(i,current_k)*mat(i,i)
     ENDDO
     IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock_gpu('matcalc')

END SUBROUTINE matcalc_gpu
!
!--------------------------------------------------------------------------
SUBROUTINE matcalc_k_gpu (label, DoE, PrtMat, ik, ninner, n, m, U, V, mat, ee)
  !------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE io_global,ONLY : stdout
  USE wvfct,                ONLY : wg, npwx
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ik, ninner, n, m
  COMPLEX(dp), INTENT(IN) :: U(ninner,n), V(ninner,m)
  COMPLEX(dp), INTENT(OUT):: mat(n,m)
  REAL(DP), INTENT(OUT) :: ee
  CHARACTER(len=*), INTENT(IN) :: label
#if defined(__CUDA)
  attributes(DEVICE) :: U, V, mat
#endif
  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock_gpu('matcalc')

  string = 'M-'
  mat = (0.0_dp, 0.0_dp)
  CALL MYZGEMM( 'C', 'N', n, m, ninner, (1.0_DP,0.0_DP), U, ninner, V, ninner, (0.0_DP,0.0_DP), mat, n )
  CALL mp_sum( mat( :, 1:m ), intra_bgrp_comm )

  IF( PrtMat > 1 ) CALL errore('matcalc_k_gpu', 'cannot print matrix', 1)

  IF(DoE) THEN
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    string = 'E-'
    ee = 0.0_dp
    !$acc parallel loop reduction(+:ee) copyin(wg)
    DO i = 1,n
      ee = ee + wg(i,ik)*DBLE(mat(i,i))
    ENDDO
    IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock_gpu('matcalc')

END SUBROUTINE matcalc_k_gpu
! NSCC
!-----------------------------------------------------------------------------
SUBROUTINE MatCholInv_gpu( MShape, n, A )
  !--------------------------------------------------------------------------
  !! Given a (real, positive definite) matrix A, returns the Cholesky factor
  !! in A (only the Lower Triangular part of the input matrix is considered).
  !!
  !! Given a real square matrix A, returns its inverse in the same shape
  !! as the input matrix, as indicated by MShape.

  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n
  !! matrix dimension
  REAL(DP), INTENT(INOUT):: A(n,n)
  !! the input-output matrix
  CHARACTER(LEN=1) :: MShape
  !! L: A is Lower Triangular (allocated in square shape);  
  !! U: A is Upper Triangular (allocated in square shape);  
  !! G: A is a general matrix.
  !
  ! ... local variables
  !
  INTEGER :: INFO, LWORK
  INTEGER, ALLOCATABLE :: IPIV(:)
  REAL(DP), ALLOCATABLE :: WORK(:)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: A
#endif


  ! MatChol
  IF(MShape.eq.'L'.or.MShape.eq.'U') then 
    INFO = -1
    CALL MYDPOTRF( Mshape, n, A, n, INFO )
    CALL errinfo( 'MYDPOTRF', 'Cholesky failed in MatCholInv.', INFO )
  ELSEIF(MShape.eq.'G') then
    CALL errinfo( 'MYDPOTRF', 'Mshape not implemented in MatCholInv.', INFO )
  ELSE
    call errore('MatCholInv', 'Wrong MShape.', 1) 
  END IF

  ! MatInv
  IF(MShape.eq.'L'.or.MShape.eq.'U') then 
    INFO = -1
    CALL MYDTRTRI( MShape, DIAG, n, A, n, INFO )
    CALL errinfo('MYDTRTRI','inversion failed in MatCholInv.',INFO)
  ELSEIF(MShape.eq.'G') then 
    CALL errinfo( 'MYDGETRI', 'Mshape not implemented in MatCholInv.', INFO )
    ! LWORK = 3*n
    ! ALLOCATE( IPIV(n), WORK(LWORK) ) 
    ! INFO = -1
    ! CALL MYDGETRF( n, n, A, n, IPIV, INFO )
    ! CALL errinfo('MYDGETRF','LU decomposition failed in MatCholInv.',INFO)
    ! INFO = -1
    ! CALL MYDGETRI( n, A, n, IPIV, WORK, LWORK, INFO )
    ! CALL errinfo('MYDGETRI','inversion failed in MatCholInv.',INFO)
    ! DEALLOCATE( IPIV, WORK ) 
  ELSE
    call errore('MatCholInv', 'Wrong MShape.', 1) 
  END IF 

END SUBROUTINE MatCholInv_gpu
!---------------------------------------------------------------------------
SUBROUTINE invchol_k_gpu( n, A )
  !---------------------------------------------------------------------------
  !! Given a matrix A, returns the inverse of the Cholesky decomposition of A
  !! for Cholesky matrices.
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n
  !! the matrix dimension
  COMPLEX(dp), INTENT(IN):: A(n,n)
  !! the input matrix
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: A
#endif
  !
  INTEGER :: INFO

  INFO = -1
  CALL MYZPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('ZPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL MYZTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('ZTRTRI','inversion failed in invchol.',INFO)
  Call MatSymm_k_gpu('L','L',A, n)

END SUBROUTINE invchol_k_gpu
!
!----------------------------------------------------------------
SUBROUTINE MatSymm_k_gpu( MShape, How, Mat, n )
  !---------------------------------------------------------------
  ! Symmetrize the (square) matrix Mat - complex output.
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=1) :: How
  !! U: copying the upper block into the lower block;  
  !! L: copying the lower block into the upper block;  
  !! S: averaging
  CHARACTER(LEN=1) :: MShape
  !! U: return the Upper Triangular (Zeros in Lower);  
  !! L: return the Lower Triangular (Zeros in Upper);  
  !! S: return the Square symmetric matrix
  INTEGER :: n
  !! the matrix dimension
  COMPLEX(DP) :: Mat(n,n)
  !! input/output matrix
  !
  ! ... local variables
  !
  INTEGER :: i, j 
  COMPLEX(DP), ALLOCATABLE :: MatT(:,:)
  REAL(DP), PARAMETER :: Zero=0.0d0, Two=2.0d0
#if defined(__CUDA)
  ATTRIBUTES(DEVICE) :: Mat, MatT
#endif
  ALLOCATE( MatT(n,n) )
 !TODO: I am here
! Properly fill the lower triangular of MatT
  MatT = (Zero,Zero) 
  IF(How.eq.'L') then ! use lower
    !$cuf kernel do(2)
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n 
        MatT(j,i) = Mat(j,i)
      end do        
    end do        
  ELSE IF( How.eq.'U' ) then ! use upper
    !$cuf kernel do(2)
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n
        MatT(j,i) = Mat(i,j)
      end do        
    end do        
  ELSE IF( How.eq.'S' ) then ! use average 
    !$cuf kernel do(2)
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n
        MatT(j,i) = (Mat(i,j) + Mat(j,i))  / Two
      end do        
    end do        
  ELSE
    Call errore('MatSymm_k_gpu','Wrong How in MatSymm_k_gpu.',1)
  END IF 

! Properly copy the results in Mat
  Mat = (Zero,Zero) 
  IF(MShape.eq.'L') then ! return lower 
    !$cuf kernel do(2)
    do i = 1, n
      do j = 1, n
        Mat(i,j) = MatT(i,j)
      enddo
    end do
  ELSE IF(MShape.eq.'U') then ! return upper 
    !$cuf kernel do(2)
    do i = 1, n
      Mat(i,i) = MatT(i,i)
      do j = i+1, n
        Mat(i,j) = MatT(j,i)   
      end do        
    end do        
  ELSE IF(MShape.eq.'S') then ! return square
    !$cuf kernel do(2)
    do i = 1, n
      do j = 1, n
        Mat(i,j) = MatT(i,j)
      enddo
    end do
    !$cuf kernel do(2)  
    do i = 1, n
      do j = i+1, n
        Mat(i,j) = MatT(j,i)   
      end do        
    end do        
  ELSE
    Call errore('MatSymm_k_gpu','Wrong MShape in MatSymm_k_gpu.',1)
  END IF 

  DEALLOCATE( MatT )

END SUBROUTINE MatSymm_k_gpu