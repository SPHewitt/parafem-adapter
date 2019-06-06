PROGRAM main

    USE mpi_wrapper;    USE precision;  USE global_variables;
    USE mp_interface;   USE input;      USE output;
    USE loading;        USE timing;     USE maths;
    USE gather_scatter; USE steering;   USE new_library;
    USE large_strain;   USE input_precice

  IMPLICIT NONE

  !----------------------------------------------------------------------
  ! PRECICE VARIABLES
  !----------------------------------------------------------------------
  CHARACTER*50             :: config, participantName, meshName
  CHARACTER*50             :: writeInitialData, readItCheckp
  CHARACTER*50             :: writeItCheckp, dataName

  INTEGER                  :: rank, commsize, ongoing, dimensions
  INTEGER                  :: meshID, vertexID, bool, vertexSize
  INTEGER                  :: displID,forceID,localVertexSize

  INTEGER,ALLOCATABLE      :: nodes(:),vertexIDs(:)
  INTEGER,ALLOCATABLE      :: local_nodes(:)

  REAL(iwp)                :: dtlimit,dt
  REAL(iwp),ALLOCATABLE    :: vertex(:),forces(:),displacements(:)
  REAL(iwp),ALLOCATABLE    :: verticies(:),local_verticies(:)
  !----------------------------------------------------------------------
  ! PARAFEM VARIABLES
  !----------------------------------------------------------------------
  CHARACTER(LEN=50)        :: argv,fname
  CHARACTER(LEN=15)        :: element

  INTEGER                  :: iel,nip,npes_pp,partitioner,ndof
  INTEGER                  :: i,j,limit,loaded_nodes,meshgen
  INTEGER                  :: ndim,nels,nn,nod,nodof,npri,nr,nres
  INTEGER                  :: nlen,flag,count,int_var(2)

  INTEGER,ALLOCATABLE      :: g_g_pp(:,:), rest(:,:),g_num_pp(:,:)

  REAL(iwp),PARAMETER      :: zero=0.0_iwp

  REAL(iwp)                :: e,v,alpha1,beta1,rho
  REAL(iwp)                :: beta,delta,tol_pcg,tol_nr
  REAL(iwp)                :: mat_prop(3),real_var(6)

  REAL(iwp),ALLOCATABLE    :: g_coord_pp(:,:,:)

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

  CALL READ_XX24(argv,numpe,alpha1,beta1,e,element,beta,delta,       &
                       limit,loaded_nodes,meshgen,nels,nip,nn,nod,dt,     &
                       npri,nr,nres,partitioner,rho,tol_pcg,tol_nr,v)
   real_var(1)=beta
   real_var(2)=delta
   real_var(3)=alpha1
   real_var(4)=beta1
   real_var(5)=tol_pcg
   real_var(6)=tol_nr
   int_var(1)=limit
   int_var(2)=npri

   mat_prop(1)=e
   mat_prop(2)=v
   mat_prop(3)=rho

   ! Calculate number of Elements per Core
  CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

  ! Degrees of Freedon per Element
  nodof = 3
  ndof  = nod*nodof
  ntot  = ndof
  ndim  = 3

  ALLOCATE(g_num_pp(nod,nels_pp))
  ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
  ALLOCATE(rest(nr,nodof+1))
  ALLOCATE(g_g_pp(ntot,nels_pp))

  g_num_pp = 0; g_coord_pp = zero; rest = 0; g_g_pp=0

  CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)
  IF(meshgen==2) CALL abaqus2sg(element,g_num_pp)
  CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
  CALL read_rest(argv,numpe,rest)

  ! Rearrange the rest Array
  CALL rearrange(rest)

  ! Clean arrays
  g_g_pp  =  0 
  neq     =  0 

  ! Find the global Steering Matrix
  elements_0: DO iel=1,nels_pp
    CALL find_g(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Stable but slow
   ! CALL find_g3(g_num_pp(:,iel),g_g_pp(:,iel),rest) ! Unstable but Fast
  END DO elements_0

  ! Build GGL Array
  neq  =  MAXVAL(g_g_pp)
  neq  =  max_p(neq)
  CALL calc_neq_pp
  CALL calc_npes_pp(npes,npes_pp)
  CALL make_ggl(npes_pp,npes,g_g_pp)

  !----------------------------------------------------------------------
  ! !. PARAFEM Input and Initialisation END
  !----------------------------------------------------------------------

  ! Allocate parafem mesh boundary
  CALL precicef_get_dims(dimensions)

  ! TODO: The current method reads the wetted surface from the .lds file.
  !       The mg2d program was bodged so that instead of writing loads
  !       it writes coordinates, this needs to be changed eventually.

  ! TODO: The current method is awfully primitive, it reads the file to get
  !       a list of nodes and count the loacal nodes it then uses another loop
  !       to populate the local nodes and verticies. It needs optimising however
  !       It doesnt currently pose a bottle neck as it occurs only once but
  !       memory issues could become a factor on bigger problems.

  ! Read in wetted surface, (xx23.lds)
  fname = argv(1:INDEX(argv," ")-1) // ".lds"
  vertexSize = loaded_nodes

  ALLOCATE(vertex(dimensions))
  ALLOCATE(nodes(vertexSize))
  ALLOCATE(verticies(dimensions*vertexSize))

  nodes=0; verticies=0; localVertexSize=0
  CALL precicef_get_mesh_id(meshName, meshID)

  OPEN(1,file=fname)
  DO i=1,vertexSize
    READ(1,*) vertexID, vertex(1), vertex(2), vertex(3)
    nodes(i)=vertexID
    verticies(i*dimensions-2)=vertex(1)
    verticies(i*dimensions-1)=vertex(2)
    verticies(i*dimensions-0)=vertex(3)
    flag=0
     DO iel=1,nels_pp
       IF(flag==1)EXIT
       DO j=1,nod
         IF(vertexID==g_num_pp(j,iel))THEN
           localVertexSize=localVertexSize+1
           flag=1
           EXIT
         ENDIF
       ENDDO
     ENDDO
  ENDDO
  CLOSE(1)

  ! Distribute the loaded nodes to local arrays

 ALLOCATE(vertexIDs(localVertexSize))
 ALLOCATE(local_nodes(localVertexSize))
 ALLOCATE(local_verticies(localVertexSize*dimensions))
 vertexIDs=0;
 count=1
 DO i =1, vertexSize
   flag=0
   vertexID=nodes(i)
   DO iel=1,nels_pp
     IF(flag==1)EXIT
     DO j=1,nod
       IF(vertexID==g_num_pp(j,iel))THEN
         local_nodes(count)=vertexID
         local_verticies(count*dimensions-2)=verticies(i*dimensions-2)
         local_verticies(count*dimensions-1)=verticies(i*dimensions-1)
         local_verticies(count*dimensions-0)=verticies(i*dimensions-0)
         count=count+1
         flag=1
         EXIT
       ENDIF
     ENDDO
   ENDDO
 ENDDO

 print*,numpe,count
 !print*,numpe,local_nodes(localVertexSize-5:localVertexSize)

  CALL precicef_set_vertices(meshID,localVertexSize,local_verticies,vertexIDs)
  DEALLOCATE(vertex)

  ALLOCATE(forces(localVertexSize*dimensions))
  ALLOCATE(displacements(localVertexSize*dimensions))

  forces=0;  displacements=0;  displID=0;  forceID=0

  dataName='Displacements0'
  CALL precicef_get_data_id(dataName,meshID,displID)

  dataName='Forces0'
  CALL precicef_get_data_id(dataName,meshID,forceID)

  forces=zero
  displacements=zero

  CALL precicef_initialize(dtlimit)
  CALL precicef_action_required(writeInitialData, bool)
  CALL precicef_initialize_data()
  CALL precicef_ongoing(ongoing)

  !----------------------------------------------------------------------
  ! Time Loop
  !----------------------------------------------------------------------
  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)

    CALL precicef_action_required(writeItCheckp, bool)
      IF (bool.EQ.1) THEN
        IF(numpe==1)WRITE(*,*)"Writing iteration checkpoint"
        CALL precicef_fulfilled_action(writeItCheckp)
      ENDIF

    CALL precicef_read_bvdata(forceID,localVertexSize,vertexIDs,forces)

    ! if bool == 0 Advance in time otherwise

    !CALL runl(nodes,forces,real_var,int_var,mat_prop,nr,loaded_nodes,dt, &
    !           g_g_pp,g_num_pp,g_coord_pp,flag,displacements)

    !CALL runnl(local_nodes,forces,real_var,int_var,mat_prop,nr,localVertexSize,dt, &
    !           g_g_pp,g_num_pp,g_coord_pp,bool,displacements,nn)
    
    CALL runnl_bd(local_nodes,forces,real_var,int_var,mat_prop,nr,localVertexSize,dt, &
               g_g_pp,g_num_pp,g_coord_pp,bool,displacements,nn)

    CALL precicef_write_bvdata(displID,localVertexSize,vertexIDs,displacements)

    CALL precicef_advance(dtlimit)

    CALL precicef_action_required(readItCheckp, bool)
    IF (bool.EQ.1) THEN
      IF(numpe==1)WRITE(*,*)"Reading iteration checkpoint"
      CALL precicef_fulfilled_action(readItCheckp)
    ELSE
      IF(numpe==1)WRITE(*,*) 'Advancing in time'

    ENDIF

    CALL precicef_ongoing(ongoing)

  ENDDO

  CALL precicef_finalize()
  WRITE (*,*) 'ParaFEM: Finalize...'

END PROGRAM
