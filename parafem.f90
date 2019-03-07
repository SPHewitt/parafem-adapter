PROGRAM main

    USE mpi_wrapper;    USE precision;  USE global_variables;
    USE mp_interface;   USE input;      USE output;
    USE loading;        USE timing;     USE maths;
    USE gather_scatter; USE steering;   USE new_library;
    USE large_strain;

  IMPLICIT NONE

  !----------------------------------------------------------------------
  ! PRECICE VARIABLES
  !----------------------------------------------------------------------
  CHARACTER*50                    :: config, participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, vertexID, bool, vertexSize
  REAL                            :: dtlimit
  REAL, DIMENSION(:), ALLOCATABLE :: vertex

  !----------------------------------------------------------------------
  ! PARAFEM VARIABLES
  !----------------------------------------------------------------------
  REAL(iwp),PARAMETER      :: zero=0.0_iwp

  REAL(iwp)                :: e,v,alpha1,beta1,omega,rho,theta,tol

  INTEGER                  :: iel,nip,npes_pp,partitioner,ndof
  INTEGER                  :: i,j,k,limit,loaded_nodes,meshgen
  INTEGER                  :: ndim,nels,nn,nod,nodof,npri,nr,nres
  INTEGER                  :: nstep,nlen

  INTEGER,ALLOCATABLE      :: g_g_pp(:,:), rest(:,:),g_num_pp(:,:)
  REAL(iwp),ALLOCATABLE    :: g_coord_pp(:,:,:)

  CHARACTER(LEN=50)        :: argv,fname
  CHARACTER(LEN=15)        :: element


  ! Constants in f90 have to be prefilled with blanks to be compatible with preCICE
  writeInitialData(1:50)='                                                  '
  readItCheckp(1:50)='                                                  '
  writeItCheckp(1:50)='                                                  '

  CALL precicef_action_write_initial_data(writeInitialData)
  CALL precicef_action_read_iter_checkp(readItCheckp)
  CALL precicef_action_write_iter_checkp(writeItCheckp)

  WRITE (*,*) 'ParaFEM: Starting solver...'
  CALL getarg(2, config)
  CALL getarg(3, participantName)
  CALL getarg(4, meshName)


  print*, config,participantName, meshName
  rank = 0
  commsize = 1
  CALL precicef_create(participantName, config, rank, commsize)

  !----------------------------------------------------------------------
  ! !. PARAFEM Input and Initialisation
  !----------------------------------------------------------------------

  CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)

  CALL read_p129(argv,numpe,alpha1,beta1,e,element,limit,loaded_nodes,&
                 meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,omega,    &
                 partitioner,rho,theta,tol,v)


  ! Calculate number of Elements per Core
  CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

  ! Degrees of Freedon per Element
  ndof  =  nod*nodof
  ntot  =  ndof

  ALLOCATE(g_num_pp(nod,nels_pp),g_coord_pp(nod,ndim,nels_pp),        &
          rest(nr,nodof+1),g_g_pp(ntot,nels_pp))
  g_num_pp = zero; g_coord_pp = zero; rest = zero;

  CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
  IF(meshgen==2) CALL abaqus2sg(element,g_num_pp)
  CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(argv,numpe,rest)

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  zero
  neq     =  zero

  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
    !CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Stable but slow
    CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Unstable but Fast
  END DO elements_0

  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp

  ! Most failures occur in this routine
  CALL calc_npes_pp(npes,npes_pp)

  CALL make_ggl(npes_pp,npes,g_g_pp)

  !----------------------------------------------------------------------
  ! !. PARAFEM Input and Initialisation END
  !----------------------------------------------------------------------


  ! Allocate parafem mesh boundary
  CALL precicef_get_dims(dimensions)

  !TODO: This is primitive does not allow, internal loads
  !      to be applied on top of wetted surface loads
  !      Fix this by writing out a second file with wetted surface
  !      during mesh gen phase

  ! Read in wetted surface, (xx23.lds)
  fname = argv(1:INDEX(argv," ")-1) // ".lds"
  vertexSize = loaded_nodes
  ALLOCATE(vertex(dimensions))
  CALL precicef_get_mesh_id(meshName, meshID)
  OPEN(1,file=fname)
  DO i=1,vertexSize
    READ(1,*) vertexID, vertex(1), vertex(2), vertex(3)
    CALL precicef_set_vertex(meshID, vertex, vertexID)
  ENDDO
  CLOSE(1)
  DEALLOCATE(vertex)

  CALL precicef_initialize(dtlimit)

  CALL precicef_action_required(writeInitialData, bool)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'ParaFEM: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)

    CALL precicef_action_required(writeItCheckp, bool)
    IF (bool.EQ.1) THEN
      CALL dummy()
      WRITE (*,*) 'ParaFEM: Writing iteration checkpoint'
      CALL precicef_fulfilled_action(writeItCheckp)
    ENDIF

    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

    CALL precicef_action_required(readItCheckp, bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'ParaFEM: Reading iteration checkpoint'
      CALL precicef_fulfilled_action(readItCheckp)
    ELSE
      WRITE (*,*) 'ParaFEM: Advancing in time'
    ENDIF

  ENDDO

  CALL precicef_finalize()
  WRITE (*,*) 'ParaFEM: Finalize...'

END PROGRAM
