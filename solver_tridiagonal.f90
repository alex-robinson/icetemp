
module solver_tridiagonal

   use defs, only : prec

   implicit none 

contains

!***************************************************************
!
! tridiag : routine de resolution d'un systeme tri-diagonal
!  M * U = R   ou M est tridiag.  (ssdiag :A, diag:B, surdiag : C)
!  U est la solution recherchee
!
!*****************************************************************
subroutine tridiag(A,B,C,R,U,n,ifail)

IMPLICIT NONE

integer, intent(in) :: n
real(prec),dimension(n), intent(in)    :: a       ! sous-diagonale l'indice 1 ne sert pas
real(prec),dimension(n), intent(inout) :: b       ! diagonale
real(prec),dimension(n), intent(in)    :: c       ! sur-diagonale, l'indice n ne sert pas
real(prec),dimension(n), intent(in)    :: r       ! vecteur membre de droite
real(prec),dimension(n), intent(out)   :: u       ! solution
integer,                 intent(out)   :: ifail ! permet de detecter une erreur

! Local variables 
real(prec),dimension(n) :: gam       ! tableau de travail.
real(prec) :: BET
integer :: JJ

ifail=0 

if (abs(B(1)).lt.1.e-20) then 
   B(1)=1.0 
   ifail=1 
endif

BET=1.0/B(1) 
U(1)=R(1)*BET

do jj=2,n
   GAM(jj)=C(jj-1)*BET
   BET=1./(B(jj)-A(jj)*GAM(jj))

   if (abs(BET)>1.0e20) then
      BET=1.0 
      ifail=1 
   endif
 
  U(jj) = (R(JJ)-A(JJ)*U(JJ-1))*BET  
end do

do JJ=N-1,1,-1
   U(JJ)=U(JJ)-GAM(JJ+1)*U(JJ+1) 
end do

if (ifail.eq.1) then 
   write(*,*) "tridiag:: Error: An element of B is null."
endif

end subroutine tridiag 

    subroutine solve_tridiag(a,b,c,d,x,n)
        ! Solve a tridiagonal system of equations
        !   a - sub-diagonal (means it is the diagonal below the main diagonal)
        !   b - the main diagonal
        !   c - sup-diagonal (means it is the diagonal above the main diagonal)
        !   d - right part
        !   x - the answer
        !   n - number of equations
        ! Note that the index i here is one based, in other words i=1,2,...,n 
        ! where n is the number of unknowns.
        ! Adapted from: https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm#Fortran_90
        ! Implements the "Thomas algorithm" for Gaussian Elimination:
        ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

        implicit none
        
        integer,intent(in) :: n
        real(prec),dimension(n),intent(in) :: a,b,c,d
        real(prec),dimension(n),intent(out) :: x
        real(prec),dimension(n) :: cp,dp
        real(prec) :: m
        integer i

        ! Initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
        
        ! Solve for vectors c-prime and d-prime
        do i = 2,n
            m = b(i)-cp(i-1)*a(i)
            cp(i) = c(i)/m
            dp(i) = (d(i)-dp(i-1)*a(i))/m
        end do
        
        ! Initialize x
        x(n) = dp(n)
        
        ! Solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
        end do

        return 

    end subroutine solve_tridiag

!   FUNCTION tridiagonal_solve(low_diagonal, central_diagonal, upper_diagonal, rhs, string_error_message) RESULT(solution)
!     ! Lapack tridiagonal solver (in double precision):
!     ! Matrix system solver for tridiagonal matrices. 
!     ! low_diagonal     = lower diagonal elements (j,j-1) of the matrix
!     ! central_diagonal = diagonal elements (j,j) of the matrix
!     ! upper_diagonal   = upper diagonal elements (j,j+1) of the matrix
!     ! rhs              = right hand side of the matrix equation
!     USE configuration_module, ONLY: dp, C
!     IMPLICIT NONE

!     ! Input variables:
!     REAL(dp), DIMENSION(2:C%NZ  ), INTENT(IN) :: low_diagonal(:)
!     REAL(dp), DIMENSION(1:C%NZ  ), INTENT(IN) :: central_diagonal(:)
!     REAL(dp), DIMENSION(1:C%NZ-1), INTENT(IN) :: upper_diagonal(:)
!     REAL(dp), DIMENSION(1:C%NZ  ), INTENT(IN) :: rhs(:)
!     CHARACTER(LEN=*),              INTENT(IN) :: string_error_message

!     ! Result variables:
!     REAL(dp), DIMENSION(1:C%NZ  )             :: solution(size(central_diagonal,1))
    
!     ! Local variables:     
!     REAL(dp), DIMENSION(2:C%NZ  )             :: low_diagonal_copy(:)
!     REAL(dp), DIMENSION(1:C%NZ  )             :: central_diagonal_copy(:)
!     REAL(dp), DIMENSION(1:C%NZ-1)             :: upper_diagonal_copy(:)
!     INTEGER                                   :: info

!     ! External subroutines:      
!     EXTERNAL DGTSV ! Lapack routine that solves tridiagonal systems (in double precision).

!     ! The LAPACK solver will overwrite the rhs with the solution. Therefore we 
!     ! first copy the rhs in the solution vector:
!     solution = rhs

!     ! The LAPACK solver will change the elements in the matrix, therefore we copy 
!     ! these INTENT(IN) in local variables:
!     low_diagonal_copy     = low_diagonal
!     central_diagonal_copy = central_diagonal
!     upper_diagonal_copy   = upper_diagonal

!     CALL DGTSV(C%NZ, 1, low_diagonal_copy, central_diagonal_copy, upper_diagonal_copy, solution, C%NZ, info)
!     ! Check if solver was successful:
!     IF(info /= 0) THEN
!      WRITE(UNIT=*, FMT='(3A, I5)') ' In the module ', string_error_message, ' in the function tridiagonal_solve: info=', info
!      STOP ' DGTSV problem with tridiagonal system, --STOPPED'
!     END IF

!     RETURN
!   END FUNCTION tridiagonal_solve

end module solver_tridiagonal
