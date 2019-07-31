  !/****h* utils
  !*  NAME
  !*    utils - Collection of fortran subroutines
  !*
  !*  SYNOPSIS
  !*    This group of subroutines is a collection of routines
  !*    that is used by the parafem solvers.
  !*
  !*    These routines are separated from the solver
  !*    fortran files to improve the structure of the software.
  !*
  !*  FUNCTION
  !*    This is a group of subroutines is required for the
  !*    the parafem-adapter. 
  !*
  !*    The routines provide adittional functionality to parafem.
  !*
  !*    Subroutine           Purpose
  !*
  !*    finddiagprecon       Generate K,M and the diagonal preconditioner
  !*    gravity_loads        Applies gravity loads
  !*
  !*    write_largestrain    Write out the data for largestrain problems
  !*    write_smallstrain    Write out the data for smallstrain problems
  !*
  !*    Functions            Purpose
  !*
  !*  AUTHOR
  !* 	  S.Hewitt
  !*  COPYRIGHT
  !*    (c) University of Manchester 1996-2017
  !******
  !*/

  SUBROUTINE find_diag_precon(store_km_pp,store_mm_pp,g_coord_pp,&
                              numVar,sProp,diag_precon_pp,time_step)
  !/****f* utils/find_diag_precon
  !*  NAME
  !*    SUBROUTINE: find_diag_precon
  !*
  !*  SYNOPSIS
  !*    Usage: find_diag_precon(store_km_pp_OF_,store_mm_pp_OF_,    &
  !*                            numSchemes_,diag_precon_pp_OF_);
  !*
  !*  FUNCTION
  !*    Calculates the diagonal preconditioner vector used
  !*    in the PCG solution.
  !*
  !*  INPUTS
  !*    store_km_pp (ntot,ntot,nels_pp)  - Stiffness Matrix
  !*    store_mm_pp (ntot,ntot,nels_pp)  - Mass Matrix
  !*
  !*    numVar      (a1 b1 theta dTim)   - Numerical Variables
  !*    time_step                        - time step
  !*    sProp       (e,v,rho)            - Rheology Properties
  !*
  !*  OUTPUT
  !*   	diag_precon_pp (neq_pp)          - Diagonal Preconditioner
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !******
  !*/

  USE mpi_wrapper;    USE precision; USE global_variables;
  USE mp_interface;   USE input;     USE output;
  USE loading;        USE timing;    USE maths;
  USE gather_scatter; USE steering;  USE new_library;

  IMPLICIT NONE

!------------------------------------------------------------------------------
! 1. Declare variables
!------------------------------------------------------------------------------

  INTEGER,PARAMETER         :: nodof=3,ndim=3,nst=6,nod=8
  REAL(iwp),PARAMETER       :: zero=0.0_iwp

  INTEGER                   :: i,k,iel,ndof,nip

  REAL(iwp),INTENT(INOUT)   :: store_km_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: store_mm_pp(ndim*nod,ndim*nod,nels_pp)
  REAL(iwp),INTENT(INOUT)   :: g_coord_pp(nod,ndim,nels_pp)
  REAL(iwp),INTENT(IN)      :: numVar(4),sProp(3),time_step

  REAL(iwp),INTENT(INOUT)   :: diag_precon_pp(neq_pp)

  REAL(iwp)                 :: alpha1,beta1,theta
  REAL(iwp)                 :: dtim,c3,c4
  REAL(iwp)                 :: det,e,rho,v,volume

  CHARACTER(LEN=15)         :: element

  LOGICAL                   :: consistent=.TRUE.

!----------------------------------------------------------------------
! 2. Declare dynamic arrays
!----------------------------------------------------------------------

  REAL(iwp),ALLOCATABLE :: diag_precon_tmp(:,:)
  REAL(iwp),ALLOCATABLE :: points(:,:),dee(:,:),weights(:)
  REAL(iwp),ALLOCATABLE :: jac(:,:),der(:,:),deriv(:,:),bee(:,:)
  REAL(iwp),ALLOCATABLE :: fun(:),emm(:,:),ecm(:,:)


  IF(numpe .EQ. 1)PRINT*,"Finding Diagonal Preconditioner: "

  ! Set Parameters
  alpha1  =  numVar(1)
  beta1   =  numVar(2)
  theta   =  numVar(3)
  dtim    =  time_step !numVar(4)
  c3      =  alpha1+1._iwp/(theta*dtim)
  c4      =  beta1+theta*dtim
  ndof    =  ntot

!----------------------------------------------------------------------
! 6. Element Stiffness and Mass Integration
!----------------------------------------------------------------------

  ! Variables required for writing and partioning
  nip=8; element="hexahedron"

  ALLOCATE(points(nip,ndim),dee(nst,nst),jac(ndim,ndim))
  ALLOCATE(der(ndim,nod),deriv(ndim,nod),bee(nst,ntot))
  ALLOCATE(weights(nip),ecm(ntot,ntot),emm(ntot,ntot),fun(nod))

  ! Youngs Modulus, Poissons Ratio and Density
  e=sProp(1);v=sProp(2);rho=sProp(3);

  ! [D] : Stress-Strain Matrix
  dee  =  zero
  CALL deemat(dee,e,v)
  CALL sample(element,points,weights)

  ! Clean [K] and [M]
  store_km_pp  =  zero
  store_mm_pp  =  zero

  elements_1: DO iel=1,nels_pp
    volume   =  zero
    emm      =  zero
    ecm      =  zero

    gauss_points_1: DO i=1,nip
      ! [N] : Shape Functions
      CALL shape_der(der,points,i)

      ! [J] : Jacobian, global to local derivative
      jac=MATMUL(der,g_coord_pp(:,:,iel))

      ! det|J|
      det=determinant(jac)
      CALL invert(jac)
      deriv=matmul(jac,der)

      ! [B] : Strain - Displacement Matrix
      CALL beemat(bee,deriv)

      ! Integration for [K]
      store_km_pp(:,:,iel)=store_km_pp(:,:,iel) +                     &
          MATMUL(MATMUL(TRANSPOSE(bee),dee),bee) * det*weights(i)
     volume=volume+det*weights(i)
     CALL shape_fun(fun,points,i)

     ! Integration for [M]
     IF(consistent)THEN
       CALL ecmat(ecm,fun,ntot,nodof)
       ecm  =  ecm*det*weights(i)*rho
       emm  =  emm+ecm
     END IF
    END DO gauss_points_1
    IF(.NOT.consistent)THEN
      DO i=1,ntot; emm(i,i)=volume*rho/13._iwp; END DO
      DO i=1,19,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=2,20,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=3,21,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=37,55,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=38,56,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
      DO i=39,57,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
    END IF
    store_mm_pp(:,:,iel)  =  emm
  END DO elements_1

  DEALLOCATE(points,dee,jac,der,deriv,bee)
  DEALLOCATE(weights,ecm,emm,fun)

!----------------------------------------------------------------------
! 2. Calculate Diagonal preconditioner
!----------------------------------------------------------------------
  ALLOCATE(diag_precon_tmp(ntot,nels_pp))

  diag_precon_pp  =  zero
  diag_precon_tmp =  zero

  elements_2: DO iel=1,nels_pp
    DO k=1,ndof
      diag_precon_tmp(k,iel)=diag_precon_tmp(k,iel)    +             &
                  store_mm_pp(k,k,iel)*c3+store_km_pp(k,k,iel)*c4
    END DO
  END DO elements_2

  CALL scatter(diag_precon_pp,diag_precon_tmp);

  DEALLOCATE(diag_precon_tmp)

  diag_precon_pp=1.0_iwp/diag_precon_pp

  END SUBROUTINE

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------


  SUBROUTINE gravity_loads(gravlo_pp,specWeight,nodof,nod,ndim,g_coord_pp)

  !/****f* parafemutils/gloads
  !*  NAME
  !*    SUBROUTINE: gloads
  !*
  !*  SYNOPSIS
  !*    Usage: gravity_loads(gravlo_pp,specWeight,nodof,nod,ndim,g_coord_pp);
  !*
  !*  FUNCTION
  !*    Supplies a gravitational force vector
  !*
  !*  INPUTS
  !*    specWeight   - Specific Weight (density  * gravitaional loading)
  !*    nodof        - Number of degrees of freedom per node
  !*    nod          - Number of Nodes per Element
  !*    ndim         - Number of dimensions
  !*    g_coord_pp   - Coordinates per processor
  !*
  !*  OUTPUT
  !*    gravlo_pp   (neq_pp)     - Vector of loads
  !*
  !*  WARNINGS
  !*    * Default direction is negative Y
  !*    * REST must be rearranged before function call
  !*
  !*  TODO
  !*    Only implemented for 8noded hexahedral elements  
  !*
  !*  AUTHOR
  !*    S. Hewitt
  !*
  !******
  !*
  !*/

  USE mpi_wrapper  
  USE precision;   USE global_variables; USE mp_interface;
  USE input;       USE output;           USE loading;
  USE timing;      USE maths;            USE gather_scatter;
  USE steering;    USE new_library;

  IMPLICIT NONE

  INTEGER                 :: i,iel,ndof,nodof,nip

  INTEGER,INTENT(IN)      :: nod,ndim

  REAL(iwp),PARAMETER     :: zero=0.0_iwp

  REAL(iwp)               :: det
  REAL(iwp)               :: pmul_pp(ntot,nels_pp)

  REAL(iwp),INTENT(IN)    :: specWeight(3),g_coord_pp(nod,ndim,nels_pp)

  REAL(iwp),INTENT(INOUT) :: gravlo_pp(neq_pp)

  CHARACTER(LEN=15)       :: element

  REAL(iwp),ALLOCATABLE   :: points(:,:),weights(:)
  REAL(iwp),ALLOCATABLE   :: fun(:),der(:,:),jac(:,:)

  !----------  Apply loading to integration Points ---------
  element="hexahedron" ; ndof = nodof*nod; nip = 8

  IF(numpe==1)PRINT*,"Calculating gravity loads:"

  ALLOCATE(points(nip,ndim))
  ALLOCATE(weights(nip),fun(nod))
  ALLOCATE(der(ndim,nod),jac(ndim,ndim))

  pmul_pp=zero

  CALL sample(element,points,weights)

  elements_1: DO iel=1,nels_pp
    gauss_points_1: DO i=1,nip
      CALL shape_der(der,points,i)
      jac=MATMUL(der,g_coord_pp(:,:,iel))
      det=determinant(jac)
      CALL shape_fun(fun,points,i)
      pmul_pp(1:ndof-2:3,iel)=pmul_pp(1:ndof-2:3,iel)+fun(:)*det*weights(i)
      pmul_pp(2:ndof-1:3,iel)=pmul_pp(2:ndof-1:3,iel)+fun(:)*det*weights(i)
      pmul_pp(3:ndof-0:3,iel)=pmul_pp(3:ndof-0:3,iel)+fun(:)*det*weights(i)
    END DO gauss_points_1

    pmul_pp(1:ndof-2:3,iel)=pmul_pp(1:ndof-2:3,iel)*specWeight(1)
    pmul_pp(2:ndof-1:3,iel)=pmul_pp(2:ndof-1:3,iel)*specWeight(2)
    pmul_pp(3:ndof-0:3,iel)=pmul_pp(3:ndof-0:3,iel)*specWeight(3)
  END DO elements_1

  gravlo_pp = zero
  CALL scatter(gravlo_pp,pmul_pp)


END SUBROUTINE gravity_loads

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  SUBROUTINE write_largestrain(job_name,nn,nr,loaded_nodes,timest,  &
                               nr_timest,inewton,nr_iters)
  !/****f* utils/write_largestrain
  !*  NAME
  !*    SUBROUTINE: write_largestrain
  !*
  !*  SYNOPSIS
  !*    Usage:      CALL write_largestrain(timest)
  !*
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and
  !*    some performance data
  !*
  !*  INPUTS
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    21.03.2018
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE mpi_wrapper;    USE precision;  USE global_variables;
  USE mp_interface;   USE input;      USE output;
  USE loading;        USE timing;     USE maths;
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: nn,nr,inewton
  INTEGER, INTENT(IN)            :: loaded_nodes,nr_iters(inewton,1)
  INTEGER                        :: i  
  REAL(iwp), INTENT(IN)          :: timest(17)
  REAL(iwp), INTENT(IN)          :: nr_timest(10,20)

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  CHARACTER(LEN=50)              :: fname

  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".run"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "

    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes
    END IF
!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                    ",&
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Load the Structure                          ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Set Starting conditions                     ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(17)-timest(1)))*100

    WRITE(11,'(/A,I1,A)')   "NEWTON RAPHSON ITERATIONS (",inewton,")"
    WRITE(11,'(/A)')   "PCG ITERATIONS "
    DO i=1,inewton
        WRITE(11,'(I4)') nr_iters(i,1)
    ENDDO

    WRITE(11,'(A,F12.6,F8.2)') "Gather Displacement Increment               ",&
                           SUM(nr_timest(:,1)),                               &
                           ((SUM(nr_timest(:,1)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Build Stiffness and Mass Matricies          ",&
                           SUM(nr_timest(:,2)),                               &
                           ((SUM(nr_timest(:,2)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Set internal Forces                         ",&
                           SUM(nr_timest(:,3)),                               &
                           ((SUM(nr_timest(:,3)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Set Newmark and Residual                    ",&
                           SUM(nr_timest(:,4)),                               &
                           ((SUM(nr_timest(:,4)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           SUM(nr_timest(:,5)),                               &
                           ((SUM(nr_timest(:,5)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           SUM(nr_timest(:,6)),                               &
                          ((SUM(nr_timest(:,6)))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Check convergence                           ",&
                           SUM(nr_timest(:,7)),                             &
                          ((SUM(nr_timest(:,7)))/(timest(17)-timest(1)))*100
    WRITE(11,'(/A,F12.6,F8.2)') "Total Time in N-R loop                      ",&
                           timest(15)-timest(5),                               &
                          ((timest(15)-timest(5))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Update Velocity and Acceleration            ",&
                            timest(16)-timest(15),                            &
                          ((timest(16)-timest(15))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Gather Data to pass out                     ",&
                            timest(17)-timest(16),                            &
                          ((timest(17)-timest(16))/(timest(17)-timest(1)))*100
    WRITE(11,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(17)-timest(1),"  100.00"
    CLOSE(11)

  END IF

  RETURN

  END SUBROUTINE write_largestrain 

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------

  SUBROUTINE write_smallstrain(job_name,timest,iters)
  !/****f* utils/write_smallstrain
  !*  NAME
  !*    SUBROUTINE: write_smallstrain
  !*
  !*  SYNOPSIS
  !*    Usage:      CALL write_smallstrain(timest)
  !*
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and
  !*    some performance data
  !*
  !*  INPUTS
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*
  !*  CREATION DATE
  !*    21.03.2018
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE mpi_wrapper;    USE precision;  USE global_variables;
  USE mp_interface;   USE input;      USE output;
  USE loading;        USE timing;     USE maths;
  USE gather_scatter; USE steering;   USE new_library;
  USE large_strain;

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  REAL(iwp), INTENT(IN)          :: timest(8)
  INTEGER,INTENT(IN)             :: iters

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  CHARACTER(LEN=50)              :: fname

  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".run"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')

!------------------------------------------------------------------------------
! 2. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                    ",&
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Load the Structure                          ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Set Starting conditions                     ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Set Displacement & Velocity                 ",&
                            timest(5)-timest(4),                            &
                          ((timest(5)-timest(4))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Build RHS                                   ",&
                            timest(6)-timest(5),                            &
                          ((timest(6)-timest(5))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)')"Solve Equations                             ",&
                            timest(7)-timest(6),                            &
                          ((timest(7)-timest(6))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,F8.2)') "Gather data to pass out                     ",&
                            timest(8)-timest(7),                            &
                          ((timest(8)-timest(7))/(timest(8)-timest(1)))*100
    WRITE(11,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(8)-timest(1),"  100.00"

    WRITE(11,'(A,12I3)')       "Number of Iterations                        ",&
                            iters
    CLOSE(11)

  END IF

  RETURN

  END SUBROUTINE write_smallstrain 

  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
