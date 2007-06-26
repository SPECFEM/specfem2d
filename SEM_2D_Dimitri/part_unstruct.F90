module part_unstruct

  implicit none
  
  include './constants_unstruct.h'

  
contains


  !-----------------------------------------------
  ! Read the mesh and storing it in array 'elmnts' (which is allocated here).
  ! 'num_start' is used to have the numbering of the nodes starting at '0'.
  ! 'nelmnts' is the number of elements, 'nnodes' is the number of nodes in the mesh.
  !-----------------------------------------------
  subroutine read_mesh(filename, nelmnts, elmnts, nnodes, num_start)

    character(len=256), intent(in)  :: filename
    integer, intent(out)  :: nelmnts
    integer, intent(out)  :: nnodes
    integer, dimension(:), pointer  :: elmnts
    integer, intent(out)  :: num_start
    
    integer  :: i

    print *, trim(filename)

    open(unit=990, file=trim(filename), form='formatted' , status='old', action='read')
    read(990,*) nelmnts
    allocate(elmnts(0:ESIZE*nelmnts-1))
    do i = 0, nelmnts-1
       read(990,*) elmnts(i*ESIZE), elmnts(i*ESIZE+1), elmnts(i*ESIZE+2), elmnts(i*ESIZE+3) 

    end do
    close(990)

    num_start = minval(elmnts)
    elmnts(:) = elmnts(:) - num_start
    nnodes = maxval(elmnts) + 1


  end subroutine read_mesh



  subroutine read_nodes_coords(filename, nnodes, nodes_coords)
    
    character(len=256), intent(in)  :: filename
    integer, intent(out)  :: nnodes
    double precision, dimension(:,:), pointer  :: nodes_coords

    integer  :: i

    print *, trim(filename)

    open(unit=991, file=trim(filename), form='formatted' , status='old', action='read')
    read(991,*) nnodes
    allocate(nodes_coords(2,nnodes))
    do i = 1, nnodes
       read(991,*) nodes_coords(1,i), nodes_coords(2,i)

    end do
    close(991)

  end subroutine read_nodes_coords



  subroutine read_mat(filename, nelmnts, num_material)
    
    character(len=256), intent(in)  :: filename
    integer, intent(in)  :: nelmnts
    integer, dimension(1:nelmnts), intent(out)  :: num_material

    integer  :: i

    print *, trim(filename)

    open(unit=992, file=trim(filename), form='formatted' , status='old', action='read')
    do i = 1, nelmnts
       read(992,*) num_material(i)

    end do
    close(992)

  end subroutine read_mat



  subroutine read_acoustic_surface(filename, nelem_acoustic_surface, acoustic_surface, &
       nelmnts, num_material, ANISOTROPIC_MATERIAL, nb_materials, icodemat, cs, num_start)
    
    include './constants.h'

    character(len=256), intent(in)  :: filename
    integer, intent(out)  :: nelem_acoustic_surface
    integer, dimension(:,:), pointer  :: acoustic_surface 
    integer, intent(in)  :: nelmnts
    integer, dimension(0:nelmnts-1)  :: num_material
    integer, intent(in)  :: ANISOTROPIC_MATERIAL
    integer, intent(in)  :: nb_materials
    integer, dimension(1:nb_materials), intent(in)  :: icodemat
    double precision, dimension(1:nb_materials), intent(in)  :: cs
    integer, intent(in)  :: num_start

    
    integer, dimension(:,:), allocatable  :: acoustic_surface_tmp 
    integer  :: nelmnts_surface
    integer  :: i
    integer  :: imaterial_number
    
    
    open(unit=993, file=trim(filename), form='formatted' , status='old', action='read')
    read(993,*) nelmnts_surface
    print *, 'POY', nelmnts_surface
    
    allocate(acoustic_surface_tmp(4,nelmnts_surface))
    
    do i = 1, nelmnts_surface
       read(993,*) acoustic_surface_tmp(1,i), acoustic_surface_tmp(2,i), acoustic_surface_tmp(3,i), acoustic_surface_tmp(4,i)
    
    end do    

    close(993)
    acoustic_surface_tmp(1,:) = acoustic_surface_tmp(1,:) - num_start
    acoustic_surface_tmp(3,:) = acoustic_surface_tmp(3,:) - num_start
    acoustic_surface_tmp(4,:) = acoustic_surface_tmp(4,:) - num_start
    
    nelem_acoustic_surface = 0
    do i = 1, nelmnts_surface
       imaterial_number = num_material(acoustic_surface_tmp(1,i))
       if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. cs(imaterial_number) < TINYVAL ) then
          nelem_acoustic_surface = nelem_acoustic_surface + 1
          
       end if
    end do
    
    allocate(acoustic_surface(4,nelem_acoustic_surface))
    
    nelem_acoustic_surface = 0
    do i = 1, nelmnts_surface
       imaterial_number = num_material(acoustic_surface_tmp(1,i))
       if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. cs(imaterial_number) < TINYVAL ) then
          nelem_acoustic_surface = nelem_acoustic_surface + 1
          acoustic_surface(:,nelem_acoustic_surface) = acoustic_surface_tmp(:,i)
       end if
    end do
       
    
  end subroutine read_acoustic_surface
  
  
  
 subroutine read_abs_surface(filename, nelemabs, abs_surface, num_start)
    
    include './constants.h'

    character(len=256), intent(in)  :: filename
    integer, intent(out)  :: nelemabs
    integer, dimension(:,:), pointer  :: abs_surface 
    integer, intent(in)  :: num_start


    integer  :: i
    
    
    open(unit=994, file=trim(filename), form='formatted' , status='old', action='read')
    read(994,*) nelemabs
    
    allocate(abs_surface(4,nelemabs))
    
    do i = 1, nelemabs
       read(994,*) abs_surface(1,i), abs_surface(2,i), abs_surface(3,i), abs_surface(4,i)
    
    end do
    
    close(994)
    
    abs_surface(1,:) = abs_surface(1,:) - num_start
    abs_surface(3,:) = abs_surface(3,:) - num_start
    abs_surface(4,:) = abs_surface(4,:) - num_start


  end subroutine read_abs_surface



  !*************************************************************************
  ! creation du graphe dual (adjacence des elements definie par au moins un noeud commun )
  !**************************************************************************
  ! Version qui construit les adjacences en commencant par dresser la liste des elements pour chaque noeud (vertices).
  ! On passe en plus en argument le tableau alloue qui contiendra la liste des elements par noeud, accompagne de la liste du
  ! nombre d'elements qui contiennent chaque noeud (c'est pour pouvoir faire l'allocation avant sans connaitre la taille exacte).
  subroutine mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts, ncommonnodes)

    integer, intent(in)  :: nelmnts
    integer, intent(in)  :: nnodes
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts 
    integer, dimension(:), pointer  :: xadj
    integer, dimension(:), pointer  :: adjncy
    integer, dimension(:), pointer  :: nnodes_elmnts
    integer, dimension(:), pointer  :: nodes_elmnts
    integer, intent(in)  :: ncommonnodes

    integer  :: i, j, k, l, m, nb_edges
    logical  ::  is_neighbour    
    integer  :: num_node, n
    integer  :: elem_base, elem_target
    integer  :: connectivity
 

    allocate(xadj(0:nelmnts))
    xadj(:) = 0
    allocate(adjncy(0:max_neighbour*nelmnts-1))
    adjncy(:) = 0
    allocate(nnodes_elmnts(0:nnodes-1))
    nnodes_elmnts(:) = 0
    allocate(nodes_elmnts(0:nsize*nnodes-1))
    nodes_elmnts(:) = 0

    nb_edges = 0

    ! remplissage de la liste des elements par noeuds
    do i = 0, esize*nelmnts-1
       nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
       nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1

    end do

    print *, 'nnodes_elmnts'
    !print *, nnodes_elmnts

    ! pour chaque noeud, remplissage du tableau des adjacences, 
    ! avec verification qu'un element ne soit pas compte plusieurs fois comme voisins d'un autre.
    do j = 0, nnodes-1
       do k = 0, nnodes_elmnts(j)-1
          do l = k+1, nnodes_elmnts(j)-1

             connectivity = 0
             elem_base = nodes_elmnts(k+j*nsize)
             elem_target = nodes_elmnts(l+j*nsize)
             do n = 1, esize
                num_node = elmnts(esize*elem_base+n-1)
                do m = 0, nnodes_elmnts(num_node)-1
                   if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                      connectivity = connectivity + 1
                   end if
                end do
             end do
             
             if ( connectivity >=  ncommonnodes) then
                
                is_neighbour = .false.
                
                do m = 0, xadj(nodes_elmnts(k+j*nsize))
                   if ( .not.is_neighbour ) then
                      if ( adjncy(nodes_elmnts(k+j*nsize)*max_neighbour+m) == nodes_elmnts(l+j*nsize) ) then
                         is_neighbour = .true.
                         
                      end if
                   end if
                end do
                if ( .not.is_neighbour ) then
                   adjncy(nodes_elmnts(k+j*nsize)*max_neighbour+xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)
                   xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                   adjncy(nodes_elmnts(l+j*nsize)*max_neighbour+xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)
                   xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                end if
             end if
          end do
       end do
    end do

    ! maintenant que la liste d'adjacence est contruite, on va refondre les tableaux dadjncy et dxadj pour qu'ils 
    ! soient au bon format  
    do i = 0, nelmnts-1
       k = xadj(i)
       xadj(i) = nb_edges
       do j = 0, k-1
          adjncy(nb_edges) = adjncy(i*max_neighbour+j)
          nb_edges = nb_edges + 1
       end do
    end do
    
    xadj(nelmnts) = nb_edges


  end subroutine mesh2dual_ncommonnodes



  subroutine read_weights(nelmnts, vwgt, nb_edges, adjwgt)
    
    integer, intent(in)  :: nelmnts, nb_edges 
    integer, dimension(:), pointer  :: vwgt, adjwgt
    
    allocate(vwgt(0:nelmnts-1))
    allocate(adjwgt(0:nb_edges-1))
    
    vwgt(:) = 1
    adjwgt(:) = 1
    
    
  end subroutine read_weights
  


  !*************************************************************************
  ! construction du tableau de correspondance entre numerotation globale et locale (pour les elements)
  !*************************************************************************/
  subroutine Construct_glob2loc_elmnts(nelmnts, part, nparts, glob2loc_elmnts)

    integer, intent(in)  :: nelmnts, nparts
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(:), pointer  :: glob2loc_elmnts

    integer  :: num_glob, num_part
    integer, dimension(0:nparts-1)  :: num_loc


    allocate(glob2loc_elmnts(0:nelmnts-1))

    do num_part = 0, nparts-1
       num_loc(num_part) = 0

    end do

    do num_glob = 0, nelmnts-1
       num_part = part(num_glob)
       glob2loc_elmnts(num_glob) = num_loc(num_part)
       num_loc(num_part) = num_loc(num_part) + 1

    end do


  end subroutine Construct_glob2loc_elmnts



 !*************************************************************************
  ! construction du tableau de correspondance entre numerotation globale et locale (pour les noeuds)
  !*************************************************************************/
  subroutine Construct_glob2loc_nodes(nelmnts, nnodes, nnodes_elmnts, nodes_elmnts, part, nparts, &
       glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes)

    integer, intent(in)  :: nelmnts, nnodes, nparts
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(0:nnodes-1), intent(in)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1), intent(in)  :: nodes_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes

    integer  :: num_node
    integer  :: el
    integer  ::  num_part
    integer  ::  size_glob2loc_nodes
    integer, dimension(0:nparts-1)  :: parts_node
    integer, dimension(0:nparts-1)  :: num_parts


    allocate(glob2loc_nodes_nparts(0:nnodes))

    size_glob2loc_nodes = 0

!!$    do num_part = 0, nparts-1
!!$       parts_node(num_part) = 0
!!$    end do
    parts_node(:) = 0


    do num_node = 0, nnodes-1
       glob2loc_nodes_nparts(num_node) = size_glob2loc_nodes
       do el = 0, nnodes_elmnts(num_node)-1
          parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1

       end do

       do num_part = 0, nparts-1
          if ( parts_node(num_part) == 1 ) then
             size_glob2loc_nodes = size_glob2loc_nodes + 1
             parts_node(num_part) = 0

          end if
       end do

    end do


    glob2loc_nodes_nparts(nnodes) = size_glob2loc_nodes 
    !print *, glob2loc_nodes_nparts
    !print *, glob2loc_nodes_nparts(nnodes)


    allocate(glob2loc_nodes_parts(0:glob2loc_nodes_nparts(nnodes)-1))
    allocate(glob2loc_nodes(0:glob2loc_nodes_nparts(nnodes)-1))

    glob2loc_nodes(0) = 0

!!$    do num_part = 0, nparts-1
!!$       parts_node(num_part) = 0
!!$       num_parts(num_part) = 0
!!$
!!$    end do
    parts_node(:) = 0
    num_parts(:) = 0
    size_glob2loc_nodes = 0


    do num_node = 0, nnodes-1
       !print *, 'num_node', num_node 
       do el = 0, nnodes_elmnts(num_node)-1
          !print *, 'el_', el
          parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1

       end do
       !print *, 'PPP' 
       do num_part = 0, nparts-1
          !print *, 'num_part', num_part 
          !print *, 'parts_nodes(num_part)', parts_node(num_part) 

          if ( parts_node(num_part) == 1 ) then
             glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
             glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
             size_glob2loc_nodes = size_glob2loc_nodes + 1
             num_parts(num_part) = num_parts(num_part) + 1
             parts_node(num_part) = 0
          end if
          !print *, 'num_part__', num_part 

       end do
    end do


  end subroutine Construct_glob2loc_nodes



  subroutine Construct_interfaces(nelmnts, nparts, part, elmnts, xadj, adjncy, tab_interfaces, &
       tab_size_interfaces, ninterfaces, nb_materials, cs_material, num_material)

    include 'constants.h'

    integer, intent(in)  :: nelmnts, nparts
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts
    integer, dimension(0:nelmnts), intent(in)  :: xadj
    integer, dimension(0:max_neighbour*nelmnts-1), intent(in)  :: adjncy
    integer, dimension(:),pointer  :: tab_size_interfaces, tab_interfaces
    integer, intent(out)  :: ninterfaces
    integer, dimension(1:nelmnts), intent(in)  :: num_material
    double precision, dimension(1:nb_materials), intent(in)  :: cs_material
    integer, intent(in)  :: nb_materials
    

    integer  :: num_part, num_part_bis, el, el_adj, num_interface, num_edge, ncommon_nodes, &
         num_node, num_node_bis
    integer  :: i, j
    logical  :: is_acoustic_el, is_acoustic_el_adj

    print *, '00'
    ninterfaces = 0
    do  i = 0, nparts-1
       do j = i+1, nparts-1
          ninterfaces = ninterfaces + 1
       end do
    end do

    allocate(tab_size_interfaces(0:ninterfaces))
    tab_size_interfaces(:) = 0

    print *, 'num_interface', ninterfaces

    num_interface = 0
    num_edge = 0
    print *, '01'

    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                if ( cs_material(num_material(el+1)) < TINYVAL) then
                   is_acoustic_el = .true.
                else
                   is_acoustic_el = .false.
                end if
                do el_adj = xadj(el), xadj(el+1)-1
                   if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                      is_acoustic_el_adj = .true.
                   else
                      is_acoustic_el_adj = .false.
                   end if
                   if ( (part(adjncy(el_adj)) == num_part_bis) .and. (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                   !if ( (part(adjncy(el_adj)) == num_part_bis) ) then
                      num_edge = num_edge + 1

                   end if
                end do
             end if
          end do
          print *, 'num_interface', num_interface
          tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge 
          num_edge = 0
          num_interface = num_interface + 1

       end do
    end do
    print *, 'POUM',tab_size_interfaces(num_interface)

    num_interface = 0
    num_edge = 0

    print *, '02'

    ! le 5 vient du fait que l'on represente l'interface comme un tableau avec 5 valeurs: element1-element2-valence-noeud1-noeud2
    allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*5-1)))
    tab_interfaces(:) = 0
    print *, '03'

    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                if ( cs_material(num_material(el+1)) < TINYVAL) then
                   is_acoustic_el = .true.
                else
                   is_acoustic_el = .false.
                end if
                do el_adj = xadj(el), xadj(el+1)-1
                   if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                      is_acoustic_el_adj = .true.
                   else
                      is_acoustic_el_adj = .false.
                   end if
                   if ( (part(adjncy(el_adj)) == num_part_bis) .and. (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                   !if ( (part(adjncy(el_adj)) == num_part_bis) ) then
                      tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+0) = el
                      tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+1) = adjncy(el_adj)
                      ncommon_nodes = 0
                      do num_node = 0, 4-1
                         do num_node_bis = 0, 4-1
                            if ( elmnts(el*esize+num_node) == elmnts(adjncy(el_adj)*esize+num_node_bis) ) then
                               tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+3+ncommon_nodes) &
                                    = elmnts(el*esize+num_node)
                               ncommon_nodes = ncommon_nodes + 1
                            end if
                         end do
                      end do
                      if ( ncommon_nodes > 0 ) then
                         tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+2) = ncommon_nodes
                      else 
                         print *, "Erreur lors de la construction des interfaces!!!", ncommon_nodes
                      end if
                      num_edge = num_edge + 1
                   end if
                end do
             end if

          end do
          num_edge = 0
          num_interface = num_interface + 1
       end do
    end do


  end subroutine Construct_interfaces



  subroutine write_glob2loc_nodes_database(IIN_database, iproc, npgeo, nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, & 
       glob2loc_nodes, nnodes, num_phase)
    
    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: nnodes, iproc, num_phase
    integer, intent(inout)  :: npgeo
    
    double precision, dimension(2,nnodes)  :: nodes_coords
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    
    integer  :: i, j

    if ( num_phase == 1 ) then 
       npgeo = 0
       
       do i = 0, nnodes-1
          do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
             if ( glob2loc_nodes_parts(j) == iproc ) then
                npgeo = npgeo + 1
                
             end if
             !write(992,*) i, glob2loc_nodes_parts(j), glob2loc_nodes(j)
             
          end do
       end do
    else
       do i = 0, nnodes-1
          do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
             if ( glob2loc_nodes_parts(j) == iproc ) then
                write(IIN_database,*) glob2loc_nodes(j)+1, nodes_coords(1,i+1), nodes_coords(2,i+1)
             end if
          end do
       end do
    end if
        
  end subroutine Write_glob2loc_nodes_database



  subroutine write_partition_database(IIN_database, iproc, nspec, nelmnts, elmnts, glob2loc_elmnts, glob2loc_nodes_nparts, &
     glob2loc_nodes_parts, glob2loc_nodes, part, num_modele, ngnod, num_phase)

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: nelmnts, num_phase, iproc
    integer, intent(inout)  :: nspec
    integer, dimension(:), pointer  :: part, elmnts, glob2loc_elmnts
    integer, dimension(:)  :: num_modele
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    integer, intent(in)  :: ngnod

    integer  :: i,j,k
    integer, dimension(0:ngnod-1)  :: loc_nodes
      

    if ( num_phase == 1 ) then 
       nspec = 0

       do i = 0, nelmnts-1
          if ( part(i) == iproc ) then
             nspec = nspec + 1
             
          end if
       end do
       
    else 
       do i = 0, nelmnts-1
          if ( part(i) == iproc ) then
             
             do j = 0, ngnod-1
                do k = glob2loc_nodes_nparts(elmnts(i*ngnod+j)), glob2loc_nodes_nparts(elmnts(i*ngnod+j)+1)-1
                  
                   if ( glob2loc_nodes_parts(k) == iproc ) then
                      loc_nodes(j) = glob2loc_nodes(k)
                      
                   end if
                end do
                
             end do
             write(IIN_database,*) glob2loc_elmnts(i)+1, num_modele(i+1), (loc_nodes(k)+1, k=0,ngnod-1)
          end if
       end do
    end if


  end subroutine write_partition_database



  subroutine Write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, nparts, iproc, ninterfaces, &
       my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
       glob2loc_nodes, num_phase)
    
    
    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: iproc
    integer, intent(in)  :: nparts
    integer, intent(in)  :: ninterfaces
    integer, intent(inout)  :: my_ninterface
    integer, dimension(:), pointer  :: tab_size_interfaces
    integer, dimension(:), pointer  :: tab_interfaces
    integer, dimension(0:ninterfaces-1), intent(inout)  :: my_interfaces
    integer, dimension(0:ninterfaces-1), intent(inout)  :: my_nb_interfaces
    integer, dimension(:), pointer  :: glob2loc_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    

    integer, dimension(2)  :: local_nodes
    integer  :: local_elmnt
    integer  :: num_phase
    
    integer  :: i, j, k, l
    integer  :: num_interface

    
    num_interface = 0
    
    if ( num_phase == 1 ) then 
       
       my_interfaces(:) = 0
       my_nb_interfaces(:) = 0
       
       do i = 0, nparts-1
          do j = i+1, nparts-1
             if ( (tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
                  (i == iproc .or. j == iproc) ) then
                my_interfaces(num_interface) = 1
                my_nb_interfaces(num_interface) = tab_size_interfaces(num_interface+1) - tab_size_interfaces(num_interface)
             end if
             num_interface = num_interface + 1
          end do
       end do
       my_ninterface = sum(my_interfaces(:))
       
    else
       
      do i = 0, nparts-1
         do j = i+1, nparts-1
            if ( my_interfaces(num_interface) == 1 ) then
               if ( i == iproc ) then
                  write(IIN_database,*) j, my_nb_interfaces(num_interface)
               else 
                  write(IIN_database,*) i, my_nb_interfaces(num_interface)
               end if
               
               do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
                  if ( i == iproc ) then
                     local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+0))+1
                  else 
                     local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+1))+1
                  end if

                  if ( tab_interfaces(k*5+2) == 1 ) then 
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(1) = glob2loc_nodes(l)+1
                        end if
                     end do
                     
                     write(IIN_database,*) local_elmnt, tab_interfaces(k*5+2), local_nodes(1), -1
                  else
                     if ( tab_interfaces(k*5+2) == 2 ) then
                        do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                             glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                           if ( glob2loc_nodes_parts(l) == iproc ) then
                              local_nodes(1) = glob2loc_nodes(l)+1
                           end if
                        end do
                        do l = glob2loc_nodes_nparts(tab_interfaces(k*5+4)), &
                           glob2loc_nodes_nparts(tab_interfaces(k*5+4)+1)-1
                           if ( glob2loc_nodes_parts(l) == iproc ) then
                              local_nodes(2) = glob2loc_nodes(l)+1
                           end if
                        end do
                        write(IIN_database,*) local_elmnt, tab_interfaces(k*5+2), local_nodes(1), local_nodes(2)
                     else 
                        write(IIN_database,*) "erreur_write_interface_", tab_interfaces(k*5+2)
                     end if
                  end if
               end do
                  
            end if
                   
            num_interface = num_interface + 1
         end do
      end do
      
   end if
    
   
 end subroutine Write_interfaces_database



 subroutine Write_surface_database(IIN_database, nsurface, surface, &
      nsurface_loc, nparts, iproc, glob2loc_elmnts, &
      glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, part, num_phase)
    
    
   integer, intent(in)  :: IIN_database
   integer, intent(in)  :: iproc
   integer, intent(in)  :: nparts
   integer  :: nsurface
   integer  :: nsurface_loc
   integer, dimension(:,:), pointer  :: surface
   
   integer, dimension(:), pointer  :: glob2loc_elmnts
   integer, dimension(:), pointer  :: glob2loc_nodes_nparts
   integer, dimension(:), pointer  :: glob2loc_nodes_parts
   integer, dimension(:), pointer  :: glob2loc_nodes
   integer, dimension(:), pointer  :: part
    
   
   integer, dimension(2)  :: local_nodes
   integer  :: local_elmnt
   integer  :: num_phase
   
   integer  :: i, j, k, l
   integer  :: num_surface
   
   
   
   if ( num_phase == 1 ) then 
       
      nsurface_loc = 0
      
      do i = 1, nsurface
         if ( part(surface(1,i)) == iproc ) then
            nsurface_loc = nsurface_loc + 1
            
         end if
      end do
      
   else
      
      nsurface_loc = 0
      
       do i = 1, nsurface
          if ( part(surface(1,i)) == iproc ) then
             nsurface_loc = nsurface_loc + 1
             
             local_elmnt = glob2loc_elmnts(surface(1,i)) + 1
             
                if ( surface(2,i) == 1 ) then
                   do l = glob2loc_nodes_nparts(surface(3,i)), &
                        glob2loc_nodes_nparts(surface(3,i)+1)-1
                      if ( glob2loc_nodes_parts(l) == iproc ) then
                         local_nodes(1) = glob2loc_nodes(l)+1
                      end if
                   end do
                   
                   write(IIN_database,*) local_elmnt, surface(2,i), local_nodes(1), -1
                end if
                if ( surface(2,i) == 2 ) then
                   do l = glob2loc_nodes_nparts(surface(3,i)), &
                        glob2loc_nodes_nparts(surface(3,i)+1)-1
                      if ( glob2loc_nodes_parts(l) == iproc ) then
                         local_nodes(1) = glob2loc_nodes(l)+1
                      end if
                   end do
                   do l = glob2loc_nodes_nparts(surface(4,i)), &
                        glob2loc_nodes_nparts(surface(4,i)+1)-1
                      if ( glob2loc_nodes_parts(l) == iproc ) then
                         local_nodes(2) = glob2loc_nodes(l)+1
                      end if
                   end do
                   
                   write(IIN_database,*) local_elmnt, surface(2,i), local_nodes(1), local_nodes(2)
                end if
                
             end if
             
          end do
          
       end if
       
       
     end subroutine Write_surface_database



     subroutine merge_abs_boundaries(nelemabs, nelemabs_merge, abs_surface, abs_surface_char, abs_surface_merge, &
          ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
          jbegin_left,jend_left,jbegin_right,jend_right, &
          nodes_coords, &
          nedges_coupled, edges_coupled, nb_materials, cs_material, num_material, &
          nelmnts, &
          elmnts, ngnod)

       implicit none
       include 'constants.h'


       integer, intent(inout)  :: nelemabs
       integer, intent(out)  :: nelemabs_merge
       integer, dimension(:,:), pointer  :: abs_surface
       logical, dimension(:,:), pointer  :: abs_surface_char
       integer, dimension(:), pointer  :: abs_surface_merge
       integer, dimension(:), pointer  :: elmnts
       integer, intent(in)  :: ngnod
       integer, dimension(:), pointer  :: ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
            jbegin_left,jend_left,jbegin_right,jend_right
       double precision, dimension(:,:), pointer  :: nodes_coords
       integer  :: nedges_coupled
       integer, dimension(:,:), pointer  :: edges_coupled
       integer  :: nb_materials
       double precision, dimension(nb_materials), intent(in)  :: cs_material
       integer, dimension(1:nelmnts), intent(in)  :: num_material
       integer  :: nelmnts


       logical, dimension(nb_materials)  :: is_acoustic
       integer  :: num_edge, nedge_bound
       integer  :: num_edge_bis
       integer  :: match
       integer  :: nb_elmnts_abs
       integer  :: i
       integer  :: temp
       integer  :: common_node, other_node
       integer  :: modified_edge_elmnt
       integer  :: modified_edge
       double precision  :: vect_product
       integer  :: iedge, inode1, inode2
       
       
    
       allocate(abs_surface_char(4,nelemabs))
       allocate(abs_surface_merge(nelemabs))
       abs_surface_char(:,:) = .false.
       abs_surface_merge(:) = -1
       
       nedge_bound = nelemabs
       nb_elmnts_abs = 0

       do num_edge = 1, nedge_bound
          
          match = 0
          do i = 1, nb_elmnts_abs
             if ( abs_surface(1,num_edge) == abs_surface_merge(i) ) then
                match = i
                exit
             end if
          end do
          
          if ( match == 0 ) then
             nb_elmnts_abs = nb_elmnts_abs + 1
             match = nb_elmnts_abs
          end if
          !print *, 'match', num_edge,match

          abs_surface_merge(match) = abs_surface(1,num_edge)
          
                   
          if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
               abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1)) ) then
             abs_surface_char(1,match) = .true.             
             
          end if
          
          if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
               abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1)) ) then
             temp = abs_surface(4,num_edge)
             abs_surface(4,num_edge) = abs_surface(3,num_edge)
             abs_surface(3,num_edge) = temp
             abs_surface_char(1,match) = .true.
           
          end if
          
          if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
               abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
             abs_surface_char(4,match) = .true.
             
          end if
          
          if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
               abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
             temp = abs_surface(4,num_edge)
             abs_surface(4,num_edge) = abs_surface(3,num_edge)
             abs_surface(3,num_edge) = temp
             abs_surface_char(4,match) = .true.
           
          end if
          
          if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
               abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2)) ) then
             abs_surface_char(2,match) = .true.

          end if
          
          if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
               abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2)) ) then
             temp = abs_surface(4,num_edge)
             abs_surface(4,num_edge) = abs_surface(3,num_edge)
             abs_surface(3,num_edge) = temp
             abs_surface_char(2,match) = .true.
           
          end if
          
          if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
               abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
             temp = abs_surface(4,num_edge)
             abs_surface(4,num_edge) = abs_surface(3,num_edge)
             abs_surface(3,num_edge) = temp
             abs_surface_char(3,match) = .true.
             
          end if
          
          if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
               abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
             abs_surface_char(3,match) = .true.
             
          end if
         
       end do
       
       nelemabs_merge = nb_elmnts_abs
       
       
       allocate(ibegin_bottom(nelemabs_merge))
       allocate(iend_bottom(nelemabs_merge))       
       allocate(jbegin_right(nelemabs_merge))
       allocate(jend_right(nelemabs_merge))
       allocate(ibegin_top(nelemabs_merge))
       allocate(iend_top(nelemabs_merge))
       allocate(jbegin_left(nelemabs_merge))
       allocate(jend_left(nelemabs_merge))

       ibegin_bottom(:) = 1
       jbegin_right(:) = 1
       ibegin_top(:) = 1
       jbegin_left(:) = 1
       iend_bottom(:) = NGLLX
       jend_right(:) = NGLLZ
       iend_top(:) = NGLLX
       jend_left(:) = NGLLZ


!!$       do num_edge = 1, nedge_bound
!!$          do num_edge_bis = 1, nedge_bound
!!$             
!!$              modified_edge = 0
!!$              
!!$              if ( abs_surface(1,num_edge) /= abs_surface(1,num_edge_bis) ) then 
!!$
!!$                 vect_product = ( (nodes_coords(1,abs_surface(4,num_edge)+1) - nodes_coords(1,abs_surface(3,num_edge)+1) ) &
!!$                      * (nodes_coords(2,abs_surface(4,num_edge_bis)+1) - nodes_coords(2,abs_surface(3,num_edge_bis)+1) ) &
!!$                      - (nodes_coords(2,abs_surface(4,num_edge)+1) - nodes_coords(2,abs_surface(3,num_edge)+1) ) &
!!$                      * (nodes_coords(1,abs_surface(4,num_edge_bis)+1) - nodes_coords(1,abs_surface(3,num_edge_bis)+1) ) )
!!$                 if ( abs(vect_product) > TINYVAL*1000000 ) then
!!$                 !if ( abs(vect_product) > 10 ) then
!!$                    
!!$                    if ( abs_surface(3,num_edge) == abs_surface(3,num_edge_bis) ) then
!!$                       common_node = abs_surface(3,num_edge)
!!$                       if ( abs_surface(1,num_edge) < abs_surface(1,num_edge_bis) ) then
!!$                          modified_edge_elmnt = abs_surface(1,num_edge_bis)
!!$                          other_node = abs_surface(4,num_edge_bis)
!!$                          modified_edge = num_edge_bis
!!$                       else
!!$                          modified_edge_elmnt = abs_surface(1,num_edge)
!!$                          other_node = abs_surface(4,num_edge)
!!$                          modified_edge = num_edge
!!$                       end if
!!$                    end if
!!$                    
!!$                    if ( abs_surface(3,num_edge) == abs_surface(4,num_edge_bis) ) then
!!$                       common_node = abs_surface(3,num_edge)
!!$                       if ( abs_surface(1,num_edge) < abs_surface(1,num_edge_bis) ) then
!!$                          modified_edge_elmnt = abs_surface(1,num_edge_bis)
!!$                          other_node = abs_surface(3,num_edge_bis)
!!$                          modified_edge = num_edge_bis
!!$                       else
!!$                          modified_edge_elmnt = abs_surface(1,num_edge)
!!$                          other_node = abs_surface(4,num_edge)
!!$                          modified_edge = num_edge
!!$                       end if
!!$                    end if
!!$                    
!!$                    if ( abs_surface(4,num_edge) == abs_surface(3,num_edge_bis) ) then
!!$                       common_node = abs_surface(4,num_edge)
!!$                       if ( abs_surface(1,num_edge) < abs_surface(1,num_edge_bis) ) then
!!$                          modified_edge_elmnt = abs_surface(1,num_edge_bis)
!!$                          other_node = abs_surface(4,num_edge_bis)
!!$                          modified_edge = num_edge_bis
!!$                       else
!!$                          modified_edge_elmnt = abs_surface(1,num_edge)
!!$                          other_node = abs_surface(3,num_edge)
!!$                          modified_edge = num_edge
!!$                       end if
!!$                    end if
!!$                    
!!$                    if ( abs_surface(4,num_edge) == abs_surface(4,num_edge_bis) ) then
!!$                       common_node = abs_surface(4,num_edge)
!!$                       if ( abs_surface(1,num_edge) < abs_surface(1,num_edge_bis) ) then
!!$                          modified_edge_elmnt = abs_surface(1,num_edge_bis)
!!$                          other_node = abs_surface(3,num_edge_bis)
!!$                          modified_edge = num_edge_bis
!!$                       else
!!$                          modified_edge_elmnt = abs_surface(1,num_edge)
!!$                          other_node = abs_surface(3,num_edge)
!!$                          modified_edge = num_edge
!!$                       end if
!!$                       
!!$                    end if
!!$                                       
!!$
!!$                    if ( modified_edge /= 0 ) then
!!$                       print *, 'SSSSSSS',  nodes_coords(1,common_node+1), nodes_coords(2,common_node+1), &
!!$                            common_node, modified_edge, modified_edge_elmnt, other_node
!!$                       match = 0
!!$                       do i = 1, nelemabs_merge
!!$                          if ( abs_surface(1,modified_edge) == abs_surface_merge(i) ) then
!!$                             match = i
!!$                             exit
!!$                          end if
!!$                       end do
!!$                       
!!$                       if ( abs_surface(3,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
!!$                            abs_surface(4,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+1) )  then
!!$                          if ( common_node == abs_surface(3,modified_edge) ) then
!!$                             ibegin_bottom(match) = 2
!!$                          else
!!$                             iend_bottom(match) = NGLLX - 1 
!!$                          end if
!!$                       end if
!!$                       
!!$                       if ( abs_surface(3,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
!!$                            abs_surface(4,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+2) )  then
!!$                          if ( common_node == abs_surface(3,modified_edge) ) then
!!$                             jbegin_right(match) = 2
!!$                          else
!!$                             jend_right(match) = NGLLZ - 1 
!!$                          end if
!!$                       end if
!!$                       
!!$                       if ( abs_surface(3,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+3) .and. &
!!$                            abs_surface(4,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+2) )  then
!!$                          if ( common_node == abs_surface(3,modified_edge) ) then
!!$                             ibegin_top(match) = 2
!!$                          else
!!$                             iend_top(match) = NGLLX - 1 
!!$                          end if
!!$                       end if
!!$                       
!!$                       if ( abs_surface(3,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
!!$                            abs_surface(4,modified_edge) == elmnts(ngnod*abs_surface_merge(match)+3) )  then
!!$                          if ( common_node == abs_surface(3,modified_edge) ) then
!!$                             jbegin_left(match) = 2
!!$                          else
!!$                             jend_left(match) = NGLLZ - 1 
!!$                          end if
!!$                       end if
!!$                       
!!$                    end if
!!$
!!$                 end if
!!$
!!$              end if
!!$              
!!$           end do
!!$           
!!$        end do

        
        is_acoustic(:) = .false.
        do i = 1, nb_materials
           if (cs_material(i) < TINYVAL) then
              is_acoustic(i) = .true.
           end if
           
        end do


        do num_edge = 1, nedge_bound

           match = 0
           do i = 1, nelemabs_merge
              if ( abs_surface(1,num_edge) == abs_surface_merge(i) ) then
                 match = i
                 exit
              end if
           end do

           if ( is_acoustic(num_material(abs_surface(1,num_edge)+1)) ) then
              
              do iedge = 1, nedges_coupled
                 
                 do inode1 = 0, 3
                    if ( abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1) ) then
                       do inode2 = 0, 3
                          if ( abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2) ) then
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) )  then 
                                ibegin_bottom(match) = 2
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then 
                                jbegin_right(match) = 2
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then 
                                ibegin_top(match) = 2
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) )  then 
                                jbegin_left(match) = 2
                                
                             end if
                             
                          end if
                       end do
                       
                    end if
                    
                    if ( abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1) ) then
                       do inode2 = 0, 3
                          if ( abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2) ) then
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) )  then 
                                iend_bottom(match) = NGLLX - 1
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then 
                                jend_right(match) = NGLLZ - 1
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then 
                                iend_top(match) = NGLLX - 1
                                
                             end if
                             if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                                  abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) )  then 
                                jend_left(match) = NGLLZ - 1
                                
                             end if
                          end if
                       end do
                       
                    end if
                    
                 end do
                 
                 
              end do
              
           end if

        end do

!!$print *, 'AAAAAAAAa'
!!$        do iedge = 1, nedges_coupled
!!$           print *, nodes_coords(1,elmnts(ngnod*edges_coupled(1,iedge))+1), &
!!$                nodes_coords(2,elmnts(ngnod*edges_coupled(1,iedge))+1), 0
!!$           print *, nodes_coords(1,elmnts(ngnod*edges_coupled(2,iedge))+1), &
!!$                nodes_coords(2,elmnts(ngnod*edges_coupled(2,iedge))+1), 0
!!$           
!!$        end do
!!$print *, 'AAAAAAAAa'       

     end subroutine merge_abs_boundaries
     
     
     
     subroutine write_abs_merge_database(IIN_database, nelemabs_merge, nelemabs_loc, &
          abs_surface_char, abs_surface_merge, &
          ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
          jbegin_left,jend_left,jbegin_right,jend_right, &
          glob2loc_elmnts, part, iproc, num_phase)
       
       implicit none
       
       integer, intent(in)  :: IIN_database
       integer, intent(in)  :: nelemabs_merge
       integer, intent(inout)  :: nelemabs_loc
       logical, dimension(:,:), pointer  :: abs_surface_char
       integer, dimension(:), pointer  :: abs_surface_merge
       integer, dimension(:), pointer  :: glob2loc_elmnts
       integer, dimension(:), pointer  :: part
       integer, intent(in)  :: iproc
       integer, intent(in)  :: num_phase
       integer, dimension(:), pointer  :: ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
            jbegin_left,jend_left,jbegin_right,jend_right
       
       integer  :: i
       
       
       if ( num_phase == 1 ) then
          nelemabs_loc = 0
          do i = 1, nelemabs_merge
             if ( part(abs_surface_merge(i)) == iproc ) then
                nelemabs_loc = nelemabs_loc + 1
             end if
          end do
       else
          do i = 1, nelemabs_merge
             if ( part(abs_surface_merge(i)) == iproc ) then
                
                write(IIN_database,*) glob2loc_elmnts(abs_surface_merge(i))+1, abs_surface_char(1,i), &
                     abs_surface_char(2,i), abs_surface_char(3,i), abs_surface_char(4,i), &
                     ibegin_bottom(i), iend_bottom(i), & 
                     jbegin_right(i), jend_right(i), & 
                     ibegin_top(i), iend_top(i), & 
                     jbegin_left(i), jend_left(i)
                
             end if
             
          end do
       end if
       
       
     end subroutine write_abs_merge_database
     
    
     
#ifdef USE_METIS
     subroutine Part_metis(nelmnts, xadj, adjncy, vwgt, adjwgt, nparts, nb_edges, edgecut, part, metis_options)

    integer, intent(in)  :: nelmnts, nparts, nb_edges 
    integer, intent(inout)  :: edgecut
    integer, dimension(0:nelmnts), intent(in)  :: xadj
    integer, dimension(0:max_neighbour*nelmnts-1), intent(in)  :: adjncy
    integer, dimension(0:nelmnts-1), intent(in)  :: vwgt
    integer, dimension(0:nb_edges-1), intent(in)  :: adjwgt
    integer, dimension(:), pointer  :: part
    integer, dimension(0:4)  :: metis_options

    integer  :: wgtflag
    integer  :: num_start

    num_start = 0
    wgtflag = 0

    print *, 'avant', edgecut
    call METIS_PartGraphRecursive(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, num_start, nparts, &
         metis_options, edgecut, part(0));
    !call METIS_PartGraphVKway(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, num_start, nparts, &
    !     options, edgecut, part(0));
    print *, 'apres', edgecut


  end subroutine Part_metis
#endif



#ifdef USE_SCOTCH
  subroutine Part_scotch(nelmnts, xadj, adjncy, vwgt, adjwgt, nparts, nedges, edgecut, part, scotch_strategy)
    
    include 'scotchf.h'

    integer, intent(in)  :: nelmnts, nparts, nedges 
    integer, intent(inout)  :: edgecut
    integer, dimension(0:nelmnts), intent(in)  :: xadj
    integer, dimension(0:max_neighbour*nelmnts-1), intent(in)  :: adjncy
    integer, dimension(0:nelmnts-1), intent(in)  :: vwgt
    integer, dimension(:), pointer  :: adjwgt
    integer, dimension(:), pointer  :: part

    double precision, dimension(SCOTCH_GRAPHDIM)  :: SCOTCHGRAPH
    double precision, dimension(SCOTCH_STRATDIM)  :: SCOTCHSTRAT
    character(len=256), intent(in)  :: scotch_strategy
    integer  :: IERR

    integer  :: wgtflag
    integer, dimension(0:4) :: options
    integer  :: num_start

    
    call scotchfstratinit (SCOTCHSTRAT(1), IERR)
     IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot initialize strat'
       STOP
    ENDIF
    
    !scotch_strategy="b{strat=x}"
    !scotch_strategy='b{strat=h{pass=10000}}'
    !scotch_strategy='b{strat=m{asc=b{strat=g{pass=100}},low=g{pass=100},type=h,vert=10}}'
    !scotch_strategy='b{strat=f{bal=0.0001,move=1000000,pass=-1}}'
    !scotch_strategy='b{strat=g}'
    !scotch_strategy='g'
    !scotch_strategy='b{strat=hf{move=1000}}'
    !scotch_strategy='b{strat=m{vert=4,asc=h|g,low=f}}'
    !scotch_strategy='b{strat=m{asc=g,low=g,rat=1.0,type=h,vert=4}}'
    !scotch_strategy='b{strat=m{asc=g,low=f{bal=0.001,move=10000,pass=500},rat=0.5,type=h,vert=4}}'
    call scotchfstratgraphmap (SCOTCHSTRAT(1), trim(scotch_strategy), IERR)
    !call scotchfstratgraphbipart (SCOTCHSTRAT(1), trim(scotch_strategy), IERR)
     IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot build strat'
       STOP
    ENDIF
    
    CALL SCOTCHFGRAPHINIT (SCOTCHGRAPH (1), IERR)
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot initialize graph'
       STOP
    ENDIF

    CALL SCOTCHFGRAPHBUILD (SCOTCHGRAPH (1), 0, nelmnts, &
         xadj (0), xadj (0), &
         xadj (0), xadj (0), &
         nedges, &
       !  adjncy (0), adjncy (0), IERR) 
         adjncy (0), adjwgt (0), IERR) 
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot build graph'
       STOP
    ENDIF
    
    CALL SCOTCHFGRAPHCHECK (SCOTCHGRAPH (1), IERR)
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Invalid check'
       STOP
    ENDIF
    
    call scotchfgraphpart (SCOTCHGRAPH (1), nparts, SCOTCHSTRAT(1), part(0), IERR)
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot part graph'
       STOP
    ENDIF

    CALL SCOTCHFGRAPHEXIT (SCOTCHGRAPH (1), IERR)
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot destroy graph'
       STOP
    ENDIF

    call scotchfstratexit (SCOTCHSTRAT(1), IERR)
    IF (IERR .NE. 0) THEN
       PRINT *, 'ERROR : MAIN : Cannot destroy strat'
       STOP
    ENDIF

    
  end subroutine Part_scotch
#endif



subroutine acoustic_elastic_repartitioning (nelmnts, nnodes, elmnts, nb_materials, cs_material, num_material, &
     nproc, part, nedges_coupled, edges_coupled)
       
  implicit none
  include 'constants.h'

  integer, intent(in)  :: nelmnts, nnodes, nproc, nb_materials
  double precision, dimension(nb_materials), intent(in)  :: cs_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material
  integer, dimension(:), pointer  :: elmnts
  integer, dimension(:), pointer :: part
  integer, intent(out)  :: nedges_coupled
  integer, dimension(:,:), pointer  :: edges_coupled
  
  
  logical, dimension(nb_materials)  :: is_acoustic
  integer, dimension(:), pointer  :: xadj
  integer, dimension(:), pointer  :: adjncy
  integer, dimension(:), pointer  :: nodes_elmnts
  integer, dimension(:), pointer  :: nnodes_elmnts
  
  integer  :: i, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned
  

  integer  :: aaa


  is_acoustic(:) = .false.
  do i = 1, nb_materials
     if (cs_material(i) < TINYVAL) then
        is_acoustic(i) = .true.
     end if
     
  end do


  call mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts,2)

!!$  aaa = 1600
!!$  print *, 'POUM', xadj(aaa), xadj(aaa+1) - 1
!!$  do i =  xadj(aaa),  xadj(aaa+1) -1
!!$     print *, i, adjncy(i)
!!$  end do
!!$  stop

  nedges_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj(el), xadj(el+1) - 1
           if ( .not. is_acoustic(num_material(adjncy(el_adj)+1)) ) then
              nedges_coupled = nedges_coupled + 1
           end if
        
        end do
     end if
  end do
  
  print *, 'nedges_coupled', nedges_coupled

  allocate(edges_coupled(2,nedges_coupled))
  
  nedges_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj(el), xadj(el+1) - 1
           if ( .not. is_acoustic(num_material(adjncy(el_adj)+1)) ) then
              nedges_coupled = nedges_coupled + 1
              edges_coupled(1,nedges_coupled) = el
              edges_coupled(2,nedges_coupled) = adjncy(el_adj)
           end if
           
        end do
     end if
  end do
  
  do i = 1, nedges_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_coupled
        if ( part(edges_coupled(1,num_edge)) /= part(edges_coupled(2,num_edge)) ) then
           if ( part(edges_coupled(1,num_edge)) < part(edges_coupled(2,num_edge)) ) then
              part(edges_coupled(2,num_edge)) = part(edges_coupled(1,num_edge))
           else
              part(edges_coupled(1,num_edge)) = part(edges_coupled(2,num_edge))
           end if
           is_repartitioned = .true.
        end if
        
     end do
     if ( .not. is_repartitioned ) then
        exit
     end if
  end do


end subroutine acoustic_elastic_repartitioning
  


subroutine write_acoustic_elastic_edges_database(IIN_database, nedges_coupled, nedges_coupled_loc, &
     edges_coupled, glob2loc_elmnts, nelmnts,part, iproc, num_phase)
       
  implicit none
  
  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: nedges_coupled
  integer, intent(inout)  :: nedges_coupled_loc
  integer, dimension(:,:), pointer  :: edges_coupled
  integer, dimension(:), pointer  :: glob2loc_elmnts
  integer, intent(in)  :: nelmnts
  integer, dimension(:), pointer  :: part
  integer, intent(in)  :: iproc
  integer, intent(in)  :: num_phase
  
  
  integer  :: i
  
  
  if ( num_phase == 1 ) then
     nedges_coupled_loc = 0
     do i = 1, nedges_coupled
        if ( part(edges_coupled(1,i)) == iproc ) then
           nedges_coupled_loc = nedges_coupled_loc + 1
        end if
     end do
  else
     do i = 1, nedges_coupled
        if ( part(edges_coupled(1,i)) == iproc ) then
           write(IIN_database,*) glob2loc_elmnts(edges_coupled(1,i))+1, glob2loc_elmnts(edges_coupled(2,i))+1 
           
        end if
        
     end do
  end if
  
  
end subroutine write_acoustic_elastic_edges_database



end module part_unstruct
