MODULE INPUT_PRECICE

 USE mpi_wrapper
 USE precision
 USE mp_interface

 CONTAINS

  SUBROUTINE READ_XX24(job_name,numpe,alpha1,beta1,e,element,beta,delta,       &
                       limit,loaded_nodes,mesh,nels,nip,nn,nod,dtim,           &
                       npri,nr,nres,partitioner,rho,tol_pcg,tol_nr,v)

  !/****f* parafemutils/read_xx24
  !*  NAME
  !*    SUBROUTINE: read_xx24
  !*  SYNOPSIS
  !*    Usage:      CALL read_xx24()
  !*  FUNCTION
  !*    Master processor reads the general data for the problem and broadcasts
  !*    it to the slave processors.
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : File name that contains the data to be read
  !*
  !*    The following scalar integer argument has the INTENT(IN) attribute:
  !*
  !*    numpe                  : Processor ID of calling processor
  !*
  !*    The following scalar real arguments have the INTENT(INOUT) attribute:
  !*
  !*    beta                   : Newmark parameters (0.25)
  !*    delta                  : Newmark parameters (0.5)
  !*    alpha1                 : Rayleigh damping parameter
  !*    beta1                  : Rayleigh damping parameter
  !*    e                      : Young's modulus
  !*    omega                  : Intermediate value
  !*    rho                    : Density
  !*    theta                  : Parameter in "theta" integrator
  !*    tol_pcg                : Convergence tolerance for PCG
  !*    tol_nr                 : Convergence of Newton-Raphson
  !*    v                      : Poisson's ratio
  !*
  !*    The following scalar integer arguments have an INTENT(INOUT) attribute:
  !*
  !*    limit                  : Maximum number of PCG iterations allowed
  !*    loaded_nodes           : Number of nodes with applied forces
  !*    mesh                   : Mesh numbering scheme
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    nels                   : Total number of elements
  !*    nip                    : Number of integration points
  !*    nn                     : Number of nodes in the mesh
  !*    nod                    : Number of nodes per element
  !*    npri                   : Print interval
  !*    nr                     : Number of restrained nodes
  !*    partitioner            : Type of partitioning
  !*                           : 1 = Smith and Griffiths internal partitioning
  !*                           : 2 = External partitioning with .psize file
  !*
  !*    The following scalar character argument has an INTENT(INOUT) attribute:
  !*
  !*    element                : Element type
  !*                           : Values: 'hexahedron' or 'tetrahedron'
  !*
  !*  AUTHOR
  !*    Sam Hewitt
  !*  CREATION DATE
  !*
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010-13
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Need to add some error traps
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)    :: job_name
  CHARACTER(LEN=15), INTENT(INOUT) :: element
  CHARACTER(LEN=50)                :: fname
  INTEGER, INTENT(IN)              :: numpe
  INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nod,nip,loaded_nodes
  INTEGER, INTENT(INOUT)           :: npri,limit,mesh,partitioner
  REAL(iwp), INTENT(INOUT)         :: rho,e,v,beta,delta,alpha1,beta1
  REAL(iwp), INTENT(INOUT)         :: tol_pcg,tol_nr,dtim

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

  INTEGER                          :: bufsize,ier,integer_store(11)
  REAL(iwp)                        :: real_store(10)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) element,mesh,partitioner,nels,nn,nr,nip,nod,loaded_nodes,nres, &
               rho,e,v,beta,delta,alpha1,beta1,dtim,npri,tol_pcg,tol_nr,limit
    CLOSE(10)

    integer_store      = 0

    integer_store(1)   = mesh
    integer_store(2)   = nels
    integer_store(3)   = nn
    integer_store(4)   = nr
    integer_store(5)   = nip
    integer_store(6)   = nod
    integer_store(7)   = loaded_nodes
    integer_store(8)   = npri
    integer_store(9)   = limit
    integer_store(10)  = partitioner
    integer_store(11)  = nres

    real_store         = 0.0_iwp

    real_store(1)      = rho
    real_store(2)      = e
    real_store(3)      = v
    real_store(4)      = beta
    real_store(5)      = delta
    real_store(6)      = alpha1
    real_store(7)      = beta1
    real_store(8)      = tol_pcg
    real_store(9)      = tol_nr
    real_store(10)     = dtim

  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------

  bufsize = 11
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  bufsize = 10
  CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

  bufsize = 15
  CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

  IF (numpe/=1) THEN

    mesh         = integer_store(1)
    nels         = integer_store(2)
    nn           = integer_store(3)
    nr           = integer_store(4)
    nip          = integer_store(5)
    nod          = integer_store(6)
    loaded_nodes = integer_store(7)
    npri         = integer_store(8)
    limit        = integer_store(9)
    partitioner  = integer_store(10)
    nres         = integer_store(11)

    rho          = real_store(1)
    e            = real_store(2)
    v            = real_store(3)
    beta         = real_store(4)
    delta        = real_store(5)
    alpha1       = real_store(6)
    beta1        = real_store(7)
    tol_pcg      = real_store(8)
    tol_nr       = real_store(9)
    dtim         = real_store(10)
  END IF

  RETURN
END SUBROUTINE READ_XX24

END MODULE INPUT_PRECICE
