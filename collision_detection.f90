module collision_detection
  use vector_operations
  implicit none

  type :: AABB
    real, dimension(3) :: min_point
    real, dimension(3) :: max_point
  end type AABB

  type :: BVHNode
    type(AABB) :: bounding_box
    integer :: left_child
    integer :: right_child
    integer :: body_index
  end type BVHNode

  type :: CollisionPair
    integer :: body1_index
    integer :: body2_index
  end type CollisionPair

contains

  function create_aabb(position, half_size) result(aabb)
    real, dimension(3), intent(in) :: position, half_size
    type(AABB) :: aabb

    aabb%min_point = position - half_size
    aabb%max_point = position + half_size
  end function create_aabb

  function aabb_overlap(aabb1, aabb2) result(overlap)
    type(AABB), intent(in) :: aabb1, aabb2
    logical :: overlap
    integer :: i

    overlap = .true.
    do i = 1, 3
      if (aabb1%max_point(i) < aabb2%min_point(i) .or. aabb1%min_point(i) > aabb2%max_point(i)) then
        overlap = .false.
        return
      end if
    end do
  end function aabb_overlap

  subroutine build_bvh(bodies, nodes, num_bodies)
    type(RigidBody), dimension(:), intent(in) :: bodies
    type(BVHNode), dimension(:), allocatable, intent(out) :: nodes
    integer, intent(in) :: num_bodies
    integer :: i, node_count

    allocate(nodes(2*num_bodies-1))
    node_count = 1

    do i = 1, num_bodies
      nodes(i)%bounding_box = create_aabb(bodies(i)%position, [1.0, 1.0, 1.0])  ! Assuming unit cube for simplicity
      nodes(i)%body_index = i
      nodes(i)%left_child = -1
      nodes(i)%right_child = -1
    end do

    call recursive_build_bvh(nodes, 1, num_bodies, node_count)
  end subroutine build_bvh

  recursive subroutine recursive_build_bvh(nodes, start, end, node_index)
    type(BVHNode), dimension(:), intent(inout) :: nodes
    integer, intent(in) :: start, end
    integer, intent(inout) :: node_index
    integer :: mid, i
    
    if (start == end) return

    mid = (start + end) / 2
    nodes(node_index)%left_child = start
    nodes(node_index)%right_child = mid + 1

    if (start < mid) then
      node_index = node_index + 1
      call recursive_build_bvh(nodes, start, mid, node_index)
    end if

    if (mid + 1 < end) then
      node_index = node_index + 1
      call recursive_build_bvh(nodes, mid + 1, end, node_index)
    end if

    nodes(node_index)%bounding_box = nodes(nodes(node_index)%left_child)%bounding_box
    do i = nodes(node_index)%left_child + 1, nodes(node_index)%right_child
      nodes(node_index)%bounding_box = merge_aabb(nodes(node_index)%bounding_box, nodes(i)%bounding_box)
    end do
  end subroutine recursive_build_bvh

  function merge_aabb(aabb1, aabb2) result(merged)
    type(AABB), intent(in) :: aabb1, aabb2
    type(AABB) :: merged
    integer :: i

    do i = 1, 3
      merged%min_point(i) = min(aabb1%min_point(i), aabb2%min_point(i))
      merged%max_point(i) = max(aabb1%max_point(i), aabb2%max_point(i))
    end do
  end function merge_aabb

  subroutine broad_phase_collision(nodes, collision_pairs, num_pairs)
    type(BVHNode), dimension(:), intent(in) :: nodes
    type(CollisionPair), dimension(:), allocatable, intent(out) :: collision_pairs
    integer, intent(out) :: num_pairs
    integer :: max_pairs, i, j

    max_pairs = size(nodes) * (size(nodes) - 1) / 2
    allocate(collision_pairs(max_pairs))
    num_pairs = 0

    call recursive_check_collision(nodes, 1, 1, collision_pairs, num_pairs)
  end subroutine broad_phase_collision

  recursive subroutine recursive_check_collision(nodes, node1, node2, collision_pairs, num_pairs)
    type(BVHNode), dimension(:), intent(in) :: nodes
    integer, intent(in) :: node1, node2
    type(CollisionPair), dimension(:), intent(inout) :: collision_pairs
    integer, intent(inout) :: num_pairs

    if (.not. aabb_overlap(nodes(node1)%bounding_box, nodes(node2)%bounding_box)) return

    if (nodes(node1)%body_index /= -1 .and. nodes(node2)%body_index /= -1) then
      if (nodes(node1)%body_index < nodes(node2)%body_index) then
        num_pairs = num_pairs + 1
        collision_pairs(num_pairs)%body1_index = nodes(node1)%body_index
        collision_pairs(num_pairs)%body2_index = nodes(node2)%body_index
      end if
      return
    end if

    if (nodes(node1)%body_index == -1) then
      call recursive_check_collision(nodes, nodes(node1)%left_child, node2, collision_pairs, num_pairs)
      call recursive_check_collision(nodes, nodes(node1)%right_child, node2, collision_pairs, num_pairs)
    else
      call recursive_check_collision(nodes, node1, nodes(node2)%left_child, collision_pairs, num_pairs)
      call recursive_check_collision(nodes, node1, nodes(node2)%right_child, collision_pairs, num_pairs)
    end if
  end subroutine recursive_check_collision

  function gjk_collision(body1, body2) result(colliding)
    type(RigidBody), intent(in) :: body1, body2
    logical :: colliding
    real, dimension(3) :: direction, support
    real, dimension(3, 4) :: simplex
    integer :: num_vertices

    direction = body2%position - body1%position
    support = get_support(body1, body2, direction)
    simplex(:, 1) = support
    num_vertices = 1
    direction = -support

    do while (.true.)
      support = get_support(body1, body2, direction)
      if (dot_product(support, direction) <= 0) then
        colliding = .false.
        return
      end if

      num_vertices = num_vertices + 1
      simplex(:, num_vertices) = support

      if (do_simplex(simplex, num_vertices, direction)) then
        colliding = .true.
        return
      end if
    end do
  end function gjk_collision

  function get_support(body1, body2, direction) result(support)
    type(RigidBody), intent(in) :: body1, body2
    real, dimension(3), intent(in) :: direction
    real, dimension(3) :: support
    
    support = furthest_point(body2, direction) - furthest_point(body1, -direction)
  end function get_support

  function furthest_point(body, direction) result(point)
    type(RigidBody), intent(in) :: body
    real, dimension(3), intent(in) :: direction
    real, dimension(3) :: point
    
    ! Assuming the body is a cube with side length 2 centered at its position
    point = body%position + sign(1.0, direction)
  end function furthest_point

  function do_simplex(simplex, num_vertices, direction) result(contains_origin)
    real, dimension(3, 4), intent(inout) :: simplex
    integer, intent(inout) :: num_vertices
    real, dimension(3), intent(inout) :: direction
    logical :: contains_origin

    select case (num_vertices)
      case (2)
        contains_origin = do_line(simplex, num_vertices, direction)
      case (3)
        contains_origin = do_triangle(simplex, num_vertices, direction)
      case (4)
        contains_origin = do_tetrahedron(simplex, num_vertices, direction)
    end select
  end function do_simplex

  function do_line(simplex, num_vertices, direction) result(contains_origin)
    real, dimension(3, 4), intent(inout) :: simplex
    integer, intent(inout) :: num_vertices
    real, dimension(3), intent(inout) :: direction
    logical :: contains_origin
    real, dimension(3) :: ab, ao

    ab = simplex(:, 1) - simplex(:, 2)
    ao = -simplex(:, 2)

    if (dot_product(ab, ao) > 0) then
      direction = cross_product(cross_product(ab, ao), ab)
    else
      simplex(:, 1) = simplex(:, 2)
      num_vertices = 1
      direction = ao
    end if

    contains_origin = .false.
  end function do_line

  function do_triangle(simplex, num_vertices, direction) result(contains_origin)
    real, dimension(3, 4), intent(inout) :: simplex
    integer, intent(inout) :: num_vertices
    real, dimension(3), intent(inout) :: direction
    logical :: contains_origin
    real, dimension(3) :: a, b, c, ab, ac, ao, abc

    a = simplex(:, 3)
    b = simplex(:, 2)
    c = simplex(:, 1)

    ab = b - a
    ac = c - a
    ao = -a

    abc = cross_product(ab, ac)

    if (dot_product(cross_product(abc, ac), ao) > 0) then
      if (dot_product(ac, ao) > 0) then
        simplex(:, 1) = simplex(:, 3)
        simplex(:, 2) = simplex(:, 1)
        num_vertices = 2
        direction = cross_product(cross_product(ac, ao), ac)
      else
        simplex(:, 2) = simplex(:, 3)
        num_vertices = 2
        direction = cross_product(cross_product(ab, ao), ab)
      end if
    else if (dot_product(cross_product(ab, abc), ao) > 0) then
      simplex(:, 2) = simplex(:, 3)
      num_vertices = 2
      direction = cross_product(cross_product(ab, ao), ab)
    else
      contains_origin = .true.
      return
    end if

    contains_origin = .false.
  end function do_triangle

  function do_tetrahedron(simplex, num_vertices, direction) result(contains_origin)
    real, dimension(3, 4), intent(inout) :: simplex
    integer, intent(inout) :: num_vertices
    real, dimension(3), intent(inout) :: direction
    logical :: contains_origin
    real, dimension(3) :: a, b, c, d, ab, ac, ad, ao, abc, acd, adb

    a = simplex(:, 4)
    b = simplex(:, 3)
    c = simplex(:, 2)
    d = simplex(:, 1)

    ab = b - a
    ac = c - a
    ad = d - a
    ao = -a

    abc = cross_product(ab, ac)
    acd = cross_product(ac, ad)
    adb = cross_product(ad, ab)

    if (dot_product(abc, ao) > 0) then
      simplex(:, 1) = simplex(:, 4)
      simplex(:, 2) = simplex(:, 3)
      simplex(:, 3) = simplex(:, 2)
      num_vertices = 3
      direction = abc
      contains_origin = .false.
    else if (dot_product(acd, ao) > 0) then
      simplex(:, 2) = simplex(:, 4)
      simplex(:, 3) = simplex(:, 2)
      num_vertices = 3
      direction = acd
      contains_origin = .false.
    else if (dot_product(adb, ao) > 0) then
      simplex(:, 3) = simplex(:, 4)
      simplex(:, 2) = simplex(:, 3)
      num_vertices = 3
      direction = adb
      contains_origin = .false.
    else
      contains_origin = .true.
    end if
  end function do_tetraheddon

end module collision_detection
