module energy_grid

  use constants,        only: MAX_LINE_LEN
  use global
  use list_header,      only: ListReal
  use output,           only: write_message

contains

!===============================================================================
! UNIONIZED_GRID creates a single unionized energy grid combined from each
! nuclide of each material. Right now, the grid for each nuclide is added into a
! linked list one at a time with an effective insertion sort. Could be done with
! a hash for all energy points and then a quicksort at the end (what hash
! function to use?)
!===============================================================================

  subroutine unionized_grid()

    integer :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc => null()

    message = "Creating unionized energy grid..."
    call write_message(5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      call add_grid_points(list, nuc % energy)
    end do

    ! Set size of unionized energy grid 
    n_grid = list % size() 

    ! create allocated array from linked list
    allocate(e_grid(n_grid))
    do i = 1, n_grid
      e_grid(i) = list % get_item(i)
    end do

    ! delete linked list and dictionary
    call list % clear()
    deallocate(list)

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

  end subroutine unionized_grid

!===============================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a linked list
! of points already stored from previous arrays.
!===============================================================================

  subroutine add_grid_points(list, energy)

    type(ListReal), pointer :: list
    real(8), intent(in) :: energy(:)

    integer :: i       ! index in energy array
    integer :: n       ! size of energy array
    integer :: current ! current index 
    real(8) :: E       ! actual energy value

    i = 1
    n = size(energy)

    ! If the original list is empty, we need to allocate the first element and
    ! store first energy point
    if (.not. associated(list)) then
      allocate(list)
      do i = 1, n
        call list % append(energy(i))
      end do
      return
    end if

    ! Set current index to beginning of the list 
    current = 1

    do while (i <= n)
      E = energy(i)

      ! If we've reached the end of the grid energy list, add the remaining
      ! energy points to the end
      if (current > list % size()) then
        ! Finish remaining energies
        do while (i <= n)
          call list % append(energy(i))
          i = i + 1
        end do
        exit
      end if

      if (E < list % get_item(current)) then

        ! Insert new energy in this position
        call list % insert(current, E)

        ! Advance index in linked list and in new energy grid
        i = i + 1
        current = current + 1

      elseif (E == list % get_item(current)) then
        ! Found the exact same energy, no need to store duplicates so just
        ! skip and move to next index
        i = i + 1
        current = current + 1
      else
        current = current + 1
      end if

    end do

  end subroutine add_grid_points

!===============================================================================
! GRID_POINTERS creates an array of pointers (ints) for each nuclide to link
! each point on the nuclide energy grid to one on the unionized energy grid
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for nuclide energy grid
    integer :: index_e      ! index on union energy grid
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    type(Nuclide), pointer :: nuc => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      allocate(nuc % grid_index(n_grid))

      index_e = 1
      energy = nuc % energy(index_e)

      do j = 1, n_grid
        union_energy = e_grid(j)
        if (union_energy >= energy .and. index_e < nuc % n_grid) then
          index_e = index_e + 1
          energy = nuc % energy(index_e)
        end if
        nuc % grid_index(j) = index_e - 1
      end do
    end do

  end subroutine grid_pointers

!===============================================================================
! UNIVERSAL_GRID creates a single universal energy grid combined from each
! nuclide of each material.
!===============================================================================

  subroutine universal_grid()

    integer :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc => null()

    message = "Creating universal energy grid..."
    call write_message(5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      call add_grid_points(list, nuc % energy)
    end do

    ! Set size of unionized energy grid 
    n_grid = list % size()

    ! Create allocated array from linked list
    allocate(e_grid(n_grid))
    do i = 1, n_grid
      e_grid(i) = list % get_item(i)
    end do

    ! Delete linked list and dictionary
    call list % clear()
    deallocate(list)

    !!!!! Add grid pointers until I think of a better way to find the cross
    ! section value for each reaction at the current particle energy
    call grid_pointers()

    ! Allocate universal energy grid. Ignoring heating for now.
    allocate(xs_grid(5, n_nuclides_total, n_grid))
    call add_nuclide_xs()

  end subroutine universal_grid

!===============================================================================
! ADD_NUCLIDE_XS adds the cross section values for each nuclide to the universal
! energy grid, interpolating where the nuclide grid does not have the
! corresponding energy in the universal grid
!===============================================================================

  subroutine add_nuclide_xs()

    integer :: i ! index in nuclides array
    integer :: j ! index in universal energy grid
    integer :: k ! index in nuclide energy grid
    integer :: m ! index in reactions array
    real(8) :: f ! interpolation factor
    type(Nuclide), pointer :: nuc => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)

      ! Set nuclide energy index and universal energy index to start
      k = 1
      j = 1

      do while (j <= n_grid)

        ! If we have reached the end of the nuclide energy grid ?? 
        if (k == nuc % n_grid) then

          xs_grid(XS_TOTAL, i, j:n_grid)        = xs_grid(XS_TOTAL, i, j - 2)
          xs_grid(XS_ELASTIC, i, j:n_grid)      = xs_grid(XS_ELASTIC, i, j - 2)
          xs_grid(XS_ABSORPTION, i, j:n_grid)   = xs_grid(XS_ABSORPTION, i, j - 2)  

          if (nuc % fissionable) then
            xs_grid(XS_NU_FISSION, i, j:n_grid) = xs_grid(XS_NU_FISSION, i, j - 2)
            xs_grid(XS_FISSION, i, j:n_grid)    = xs_grid(XS_FISSION, i, j - 2)
          end if

          exit

        ! If nuclide energy grid contains this point on the universal grid, add
        ! the corresponding cross section values from each reaction type
        elseif (nuc % energy(k) == e_grid(j)) then

          xs_grid(XS_TOTAL, i, j)        = nuc % total(k)
          xs_grid(XS_ELASTIC, i, j)      = nuc % elastic(k)
          xs_grid(XS_ABSORPTION, i, j)   = nuc % absorption(k)

          if (nuc % fissionable) then
            xs_grid(XS_NU_FISSION, i, j) = nuc % nu_fission(k)
            xs_grid(XS_FISSION, i, j)    = nuc % fission(k)
          end if

          k = k + 1
          j = j + 1

        ! While energy on nuclide grid is greater than energy on universal grid,
        ! interpolate the microscopic cross sections in between
        else
          do while (e_grid(j) < nuc % energy(k))

            ! Calculate interpolation factor
            f = (e_grid(j) - nuc % energy(k - 1)) / (nuc % energy(k) - nuc % energy(k - 1))

            ! Calculate the microscopic nuclide total cross sections
            xs_grid(XS_TOTAL, i, j) = (ONE - f) * nuc % total(k - 1) &
                 + f * nuc % total(k)

            ! Calculate the microscopic nuclide elastic cross section
            xs_grid(XS_ELASTIC, i, j) = (ONE - f) * nuc % elastic(k - 1) &
                 + f * nuc % elastic(k)

            ! Calculate the microscopic nuclide absorption cross section
            xs_grid(XS_ABSORPTION, i, j) = (ONE - f) * nuc % absorption(k - 1) &
                 + f * nuc % absorption(k)

            if (nuc % fissionable) then

              ! Calculate the microscopic nuclide fission cross section
              xs_grid(XS_FISSION, i, j) = (ONE - f) * nuc % fission(k - 1) &
                   + f * nuc % fission(k)

              ! Calculate the microscopic nuclide nu-fission cross section
              xs_grid(XS_NU_FISSION, i, j) = (ONE - f) * nuc % nu_fission(k - 1) &
                   + f * nuc % nu_fission(k)
            end if

            j = j + 1

          end do
        end if
      end do

    ! Deallocate the nuclide energy grid and cross section arrays
    call nuc % xs_clear()

    end do

  end subroutine add_nuclide_xs

end module energy_grid
