module grid_size

  use ace_header,        only: Nuclide, Reaction
  use constants
  use error,             only: fatal_error
  use global
  use output,            only: write_message
  use search,            only: binary_search
  use string,            only: to_str

  implicit none

contains

!===============================================================================
! GET_XS_SIZES determines and prints the current sizes of xs and energy grids.
! Unallocated arrays will print zero size.
!===============================================================================

  subroutine get_grid_sizes()
    integer(8) :: n_xs_tot   ! sum over nuclides of total cross section sizes
    integer(8) :: n_xs_main  ! sum over nuclides of tot, abs, el, (fis) sizes
    integer(8) :: n_xs_rxn   ! sum over nuclides of all reactions' sigmas
    integer(8) :: n_xs_urr   ! size of URR energy grid and probability tables
    integer(8) :: n_xs_2nd   ! size of secondary distributions for reactions
    integer(8) :: n_xs_other ! size of other xs uncounted previously
    integer(8) :: n_grid     ! sum of all grid index sizes
    integer(8) :: n_e        ! sum over nuclides of the local energy grid sizes
    integer(8) :: n_ueg      ! union energy grid size
    integer(8) :: n_aeg      ! fractional cascading augmented energy grid size
    integer(8) :: n_ptrs     ! fractional cascading pointers size
    integer(8) :: n_tot      ! total size
    integer(8) :: n_rel      ! total relevant size
    integer :: sz_real       ! size of real(8) in B
    integer :: sz_int        ! size of int in B
    integer :: i_nuclide     ! iteration index over nuclides
    integer :: i_rxn         ! iteration index over reactions
    integer :: i_2nd         ! iteration index over secondary energy distribution
    integer :: KB            ! conversion from B to KB
    integer :: MB            ! conversion from B to MB
    character(len=10) :: num ! temporary storage
    type(Nuclide), pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()

    ! Initialize to zero
    n_rel = ZERO
    n_xs_tot = ZERO
    n_xs_main = ZERO
    n_xs_rxn = ZERO
    n_xs_urr = ZERO
    n_xs_2nd = ZERO
    n_xs_other = ZERO
    n_grid = ZERO
    n_e = ZERO
    n_aeg = ZERO
    n_ptrs = ZERO
    n_ueg = ZERO

    ! Conversion sizes
    sz_real = 8
    sz_int = 4
    KB = 1024
    MB = KB * KB

    ! Loop over nuclides and accumulate sizes
    do i_nuclide = 1, n_nuclides_total
      nuc => nuclides(i_nuclide)
      ! Total
      n_xs_tot = n_xs_tot + sz_real * size(nuc % total)

      ! Main
      n_xs_main = n_xs_main + sz_real * &
        (size(nuc % total) + size(nuc % elastic))
      if (allocated(nuc % absorption)) then
        n_xs_main = n_xs_main + sz_real * size(nuc % absorption)
      end if
      if (nuc % fissionable) then
        n_xs_main = n_xs_main + sz_real * &
          (size(nuc % fission) + size(nuc % nu_fission))
      end if

      ! Reactions
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        if (allocated(rxn % sigma)) then
          n_xs_rxn = n_xs_rxn + sz_real * size(rxn % sigma)
        end if
      end do

      ! Secondary
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        if (allocated(rxn % adist % energy)) then
          n_xs_2nd = n_xs_2nd + sz_real * size(rxn % adist % energy)
        end if
        if (allocated(rxn % adist % type)) then
          n_xs_2nd = n_xs_2nd + sz_int * size(rxn % adist % type)
        end if
        if (allocated(rxn % adist % location)) then
          n_xs_2nd = n_xs_2nd + sz_int * size(rxn % adist % location)
        end if
        if (allocated(rxn % adist % data)) then
          n_xs_2nd = n_xs_2nd + sz_real * size(rxn % adist % data)
        end if
        if (associated(rxn % edist)) then
          if (allocated(rxn % edist % data)) then
            n_xs_2nd = n_xs_2nd + sz_real * size(rxn % edist % data)
          end if
        end if
      end do

      ! URR
      if (nuc % urr_present) then
        n_xs_urr = n_xs_urr + sz_real * ( &
          size(nuc % urr_data % energy) + size(nuc % urr_data % prob, 1) * &
          size(nuc % urr_data % prob, 2) * size(nuc % urr_data % prob, 3))
      end if

      ! Other
      if (nuc % fissionable .and. allocated(nuc % nu_t_data)) then
        n_xs_other = n_xs_other + sz_real * size(nuc % nu_t_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_p_data)) then
        n_xs_other = n_xs_other + sz_real * size(nuc % nu_p_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_d_data)) then
        n_xs_other = n_xs_other + sz_real * size(nuc % nu_d_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_d_precursor_data)) then
        n_xs_other = n_xs_other + sz_real * size(nuc % nu_d_precursor_data)
      end if
      if (nuc % fissionable .and. associated(nuc % nu_d_edist)) then
        do i_2nd = 1, size(nuc % nu_d_edist)
          if (allocated(nuc % nu_d_edist(i_2nd) % data)) then
            n_xs_other = n_xs_other + &
              sz_real * size(nuc % nu_d_edist(i_2nd) % data)
          end if
        end do
      end if

      ! Grid index
      if (allocated(nuc % grid_index)) then
        n_grid = n_grid + sz_int * size(nuc % grid_index)
      end if

      ! Nuclide energy grid
      n_e = n_e + sz_real * size(nuc % energy)

      ! Augmented energy grid
      if (allocated(nuc % aug_energy)) then
        n_aeg = n_aeg + sz_real * size(nuc % aug_energy)
      end if

      if (allocated(nuc % aug_pointers)) then
        n_ptrs = n_ptrs + sz_int * size(nuc % aug_pointers)
      end if

    end do

    ! Union energy grid
    if (allocated(e_grid)) then
      n_ueg = sz_real * size(e_grid)
    end if

    ! Cumulative (overflow issues)
    n_tot = int(n_xs_main, 8) + int(n_xs_rxn, 8) + int(n_xs_urr, 8) + &
         int(n_xs_2nd, 8) + int(n_xs_other, 8) + int(n_grid, 8) + int(n_e, 8) + &
         int(n_ueg, 8) + int(n_aeg, 8) + int(n_ptrs, 8)

    n_rel = int(n_xs_main, 8) + int(n_grid, 8) + int(n_ueg, 8) + &
         int(n_aeg, 8) + int(n_ptrs, 8) + int(n_e, 8)

    ! Convert from B to KB
    n_xs_tot = n_xs_tot / KB
    n_xs_main = n_xs_main / KB
    n_xs_rxn = n_xs_rxn  / KB
    n_xs_urr = n_xs_urr / KB
    n_xs_2nd = n_xs_2nd / KB
    n_xs_other = n_xs_other / KB
    n_grid = n_grid / KB
    n_e = n_e / KB
    n_ueg = n_ueg / KB
    n_aeg = n_aeg / KB
    n_ptrs = n_ptrs / KB
    n_tot = n_tot / MB
    n_rel = n_rel / MB

    ! Print sizes
    write(num, '(i10)') n_xs_tot
    message = 'Size of all total xs is           ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_xs_main
    message = 'Size of all main xs is            ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_xs_rxn
    message = "Size of all reactions' sigmas is  " // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_xs_urr
    message = 'Size of all urr data is           ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_xs_2nd
    message = "Size of all rxns' 2ndary xs is    " // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_xs_other
    message = 'Size of other fission data is     ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_grid
    message = "Size of nuclides' index grids is  " // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_e
    message = "Size of nuclides' energy grids is " // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_ueg
    message = 'Size of unionized energy grid is  ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_aeg
    message = 'Size of augmented energy grid is  ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_ptrs
    message = 'Size of grid pointers is          ' // num // ' KB'
    call write_message(5)

    write(num, '(i10)') n_tot
    message = 'Total size is                     ' // num // ' MB'
    call write_message(5)

    write(num, '(i10)') n_rel
    message = 'Total energy and xs grid size is  ' // num // ' MB'
    call write_message(5)

  end subroutine get_grid_sizes

end module
