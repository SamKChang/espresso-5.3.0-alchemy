!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!
!----------------------------------------------------------------------------
SUBROUTINE alchemy( npw, npwx, nvec, nvecx, npol, evc, ethr, &
                    uspp, e, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  USE kinds,         ONLY : DP
  USE mp_bands,      ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp, nbgrp
  USE mp,            ONLY : mp_sum, mp_bcast
  USE control_flags, ONLY : alchemy_pred
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  !
  REAL(DP), EXTERNAL :: ddot

  !
  ! EXTERNAL  h_psi,    s_psi,    g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  CALL start_clock( 'alchemy' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'alchemy', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ALLOCATE(  psi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' alchemy ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate sc ', ABS(ierr) )
  ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate hc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' alchemy ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = 0
  nbase  = nvec
  conv   = .TRUE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,:,1:nvec) = evc(:,:,1:nvec)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi( npwx, npw, nvec, psi, hpsi )
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi( npwx, npw, nvec, psi, spsi )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  !
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi, kdmx, hpsi, kdmx, ZERO, hc, nvecx )
  !
  CALL mp_sum( hc( :, 1:nbase ), intra_bgrp_comm )
  !
!  IF ( uspp ) THEN
!     !
!     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
!                 psi, kdmx, spsi, kdmx, ZERO, sc, nvecx )
!     !     
!  ELSE
!     !
!     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
!                 psi, kdmx, psi, kdmx, ZERO, sc, nvecx )
!     !
!  END IF
!  !
!  CALL mp_sum( sc( :, 1:nbase ), intra_bgrp_comm )
!  !
  IF ( lrot ) THEN
     DO n = 1, nbase
        !
        e(n) = REAL( hc(n,n) ) ! <<< Alchemical derivatives!!!
        !
!        vc(n,n) = ONE
        !
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
!     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     !
!     e(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  CALL stop_clock( 'alchemy' )
  !
  RETURN
  !
END SUBROUTINE alchemy
