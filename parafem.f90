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
  CHARACTER*50                    :: config, participantName, meshName
  CHARACTER*50                    :: writeInitialData, readItCheckp
  CHARACTER*50                    :: writeItCheckp, dataName
  INTEGER                         :: rank, commsize, ongoing, dimensions
  INTEGER                         :: meshID, vertexID, bool, vertexSize
  INTEGER                         :: displID,forceID,hasData
  REAL(iwp)                       :: dtlimit,dt
  REAL(iwp), DIMENSION(:), ALLOCATABLE :: vertex,forces,displacements,nodes

  !----------------------------------------------------------------------
  ! PARAFEM VARIABLES
  !----------------------------------------------------------------------
  REAL(iwp),PARAMETER      :: zero=0.0_iwp

  REAL(iwp)                :: e,v,alpha1,beta1,omega,rho,theta,tol

  INTEGER                  :: iel,nip,npes_pp,partitioner,ndof
  INTEGER                  :: i,j,k,limit,loaded_nodes,meshgen
  INTEGER                  :: ndim,nels,nn,nod,nodof,npri,nr,nres
  INTEGER                  :: nstep,nlen,flag

  INTEGER,ALLOCATABLE      :: g_g_pp(:,:), rest(:,:),g_num_pp(:,:)

  REAL(iwp),ALLOCATABLE    :: g_coord_pp(:,:,:)

  REAL(iwp)                :: mat_prop(3),num_var(5)

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


  CALL find_pe_procs(numpe,npes); CALL getname(argv,nlen)

  rank = numpe-1
  commsize = npes

  CALL precicef_create(participantName, config, rank, commsize)

  !----------------------------------------------------------------------
  ! !. PARAFEM Input and Initialisation
  !----------------------------------------------------------------------

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

  ALLOCATE(vertex(dimensions),nodes(vertexSize))

  CALL precicef_get_mesh_id(meshName, meshID)

  OPEN(1,file=fname)
  DO i=1,vertexSize
    READ(1,*) vertexID, vertex(1), vertex(2), vertex(3)
    nodes(i)=vertexID
    CALL precicef_set_vertex(meshID, vertex, vertexID)
  ENDDO
  CLOSE(1)
  DEALLOCATE(vertex)
  ALLOCATE(forces(vertexSize*dimensions),displacements(vertexSize*dimensions))

  displID=0
  forceID=0

  !dataName='Displacements0'
  !print*,dataName
  
  !CALL precicef_get_data_id(dataName,meshID,displID,50)

  dataName='Forces'
  print*,dataName
 CALL precicef_has_data(dataName, hasData, 50, meshID)
 print*,"HELL"
 CALL precicef_get_data_id(dataName,meshID,forceID,7)

  forces=0.0
  displacements=0.0

  CALL precicef_initialize(dtlimit)

  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  !! HARDCODED properties for the sake of quick and easy test Case
  num_var(1)=alpha1
  num_var(2)=beta1
  num_var(3)=theta
  num_var(4)=tol
  num_var(5)=limit

  mat_prop(1)=e
  mat_prop(2)=v
  mat_prop(3)=rho

  dt = 0.001

    print*,"JHELLLo"
  !----------------------------------------------------------------------
  ! Time Loop
  !----------------------------------------------------------------------
  DO WHILE (ongoing.NE.0)
    CALL precicef_read_bvdata(forceID,vertexSize,nodes,forces)

    CALL runl(nodes,displacements,num_var,mat_prop,nr,loaded_nodes,dt, &
               g_g_pp,g_num_pp,g_coord_pp,flag)

    CALL precicef_write_bvdata(displID,vertexSize,nodes,displacements)

    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)


  ENDDO

  CALL precicef_finalize()
  WRITE (*,*) 'ParaFEM: Finalize...'

END PROGRAM
