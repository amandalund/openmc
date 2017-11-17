module tracking

  use constants
  use cross_section,      only: calculate_xs
  use error,              only: fatal_error, warning
  use geometry_header,    only: cells
  use geometry,           only: find_cell, distance_to_boundary, cross_surface, &
                                cross_lattice, check_cell_overlap
  use output,             only: write_message
  use material_header,    only: materials
  use message_passing
  use mgxs_header
  use nuclide_header
  use particle_header,    only: LocalCoord, Particle
  use physics,            only: collision
  use physics_mg,         only: collision_mg
  use random_lcg,         only: prn
  use settings
  use simulation_header
  use string,             only: to_str
  use tally_header
  use tally,              only: score_analog_tally, score_tracklength_tally, &
                                score_collision_tally, score_surface_current, &
                                score_track_derivative, score_surface_tally, &
                                score_collision_derivative, zero_flux_derivs
  use track_output,       only: initialize_particle_track, write_particle_track, &
                                add_particle_track, finalize_particle_track

  implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), intent(inout) :: p

    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! sampled distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: tau_hat                ! optical depth used in newton's method
    real(8) :: optical_depth          ! optical depth integral
    real(8) :: pnc                    ! probilitiy of no collision
    real(8) :: xyz_orig(3)            ! original xyz coordinate of the particle

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      call write_message("Simulating Particle " // trim(to_str(p % id)))
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
!$omp atomic
    total_weight = total_weight + p % wgt

    ! Force calculation of cross-sections by setting last energy to zero
    if (run_CE) then
      micro_xs % last_E = ZERO
    end if

    ! Prepare to write out particle track.
    if (p % write_track) then
      call initialize_particle_track()
    endif

    ! Every particle starts with no accumulated flux derivative.
    if (active_tallies % size() > 0) call zero_flux_derivs()

    EVENT_LOOP: do
      ! Store pre-collision particle properties
      p % last_wgt = p % wgt
      p % last_E   = p % E
      p % last_uvw = p % coord(1) % uvw
      p % last_xyz = p % coord(1) % xyz

      ! If the cell hasn't been determined based on the particle's location,
      ! initiate a search for the current cell. This generally happens at the
      ! beginning of the history and again for any secondary particles
      if (p % coord(p % n_coord) % cell == NONE) then
        call find_cell(p, found_cell)
        if (.not. found_cell) then
          call fatal_error("Could not locate particle " // trim(to_str(p % id)))
        end if

        ! set birth cell attribute
        if (p % cell_born == NONE) p % cell_born = p % coord(p % n_coord) % cell
      end if

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

      if (check_overlaps) call check_cell_overlap(p)

      ! Calculate microscopic and macroscopic cross sections
      if (run_CE) then
        ! If the material is the same as the last material and the temperature
        ! hasn't changed, we don't need to lookup cross sections again.
        if (p % material /= p % last_material .or. &
             p % sqrtkT /= p % last_sqrtkT) call calculate_xs(p)
      else
        ! Since the MGXS can be angle dependent, this needs to be done
        ! After every collision for the MGXS mode
        if (p % material /= MATERIAL_VOID) then
          ! Update the temperature index
          call macro_xs(p % material) % obj % find_temperature(p % sqrtkT)
          ! Get the data
          call macro_xs(p % material) % obj % calculate_xs(p % g, &
               p % coord(p % n_coord) % uvw, material_xs)
        else
          material_xs % total      = ZERO
          material_xs % absorption = ZERO
          material_xs % nu_fission = ZERO
        end if

        ! Finally, update the particle group while we have already checked for
        ! if multi-group
        p % last_g = p % g
      end if

      ! Find the distance to the nearest boundary
      call distance_to_boundary(p, d_boundary, surface_crossed, &
           lattice_translation, next_level)

      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY

      ! Continuous material tracking
      else if (materials(p % material) % continuous_num_density) then
        ! Integrate along complete neutron flight path
        call simpsons_path_integration(p, optical_depth, xs_t, d_boundary)

        ! Sample the collision for analytic density
        pnc = exp(-optical_depth)

        ! No Collision
        if (prn() <= pnc) then
          d_collision = INFINITY

        ! Collision
        else
          ! Sample optical depth
          tau_hat = -log(ONE - (ONE - pnc) * prn())

          ! Get the flight distance for sampled optical depth
          call estimate_flight_distance(xs_t, d_boundary, tau_hat, d_collision)

          ! Move particle to the point of collision so we can make sure that
          ! cross sections are updated for later tallying of keff and user
          ! tallies
          xyz_orig = p % coord(p % n_coord) % xyz
          call move_particle_coord(p % coord(p % n_coord), d_collision)
          call calculate_xs(p)

          ! Move particle back
          p % coord(p % n_coord) % xyz = xyz_orig
        end if

      else
        d_collision = -log(prn()) / material_xs % total
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)

      ! Advance particle
      do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) then
        call score_tracklength_tally(p, distance)
      end if

      ! Score track-length estimate of k-eff
      if (run_mode == MODE_EIGENVALUE) then
        global_tally_tracklength = global_tally_tracklength + p % wgt * &
             distance * material_xs % nu_fission
      end if

      ! Score flux derivative accumulators for differential tallies.
      if (active_tallies % size() > 0) call score_track_derivative(p, distance)

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        if (next_level > 0) p % n_coord = next_level

        ! Saving previous cell data
        do j = 1, p % n_coord
          p % last_cell(j) = p % coord(j) % cell
        end do
        p % last_n_coord = p % n_coord

        p % coord(p % n_coord) % cell = NONE
        if (any(lattice_translation /= 0)) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(p, lattice_translation)
          p % event = EVENT_LATTICE
        else
          ! Particle crosses surface
          p % surface = surface_crossed

          call cross_surface(p)
          p % event = EVENT_SURFACE
        end if
        ! Score cell to cell partial currents
        if(active_surface_tallies % size() > 0) call score_surface_tally(p)
      else
        ! ====================================================================
        ! PARTICLE HAS COLLISION

        ! Score collision estimate of keff
        if (run_mode == MODE_EIGENVALUE) then
          global_tally_collision = global_tally_collision + p % wgt * &
               material_xs % nu_fission / material_xs % total
        end if

        ! score surface current tallies -- this has to be done before the collision
        ! since the direction of the particle will change and we need to use the
        ! pre-collision direction to figure out what mesh surfaces were crossed

        if (active_current_tallies % size() > 0) call score_surface_current(p)

        ! Clear surface component
        p % surface = NONE

        if (run_CE) then
          call collision(p)
        else
          call collision_mg(p)
        end if

        ! Score collision estimator tallies -- this is done after a collision
        ! has occurred rather than before because we need information on the
        ! outgoing energy for any tallies with an outgoing energy filter
        if (active_collision_tallies % size() > 0) call score_collision_tally(p)
        if (active_analog_tallies % size() > 0) call score_analog_tally(p)

        ! Reset banked weight during collision
        p % n_bank   = 0
        p % wgt_bank = ZERO
        p % n_delayed_bank(:) = 0

        ! Reset fission logical
        p % fission = .false.

        ! Save coordinates for tallying purposes
        p % last_xyz_current = p % coord(1) % xyz

        ! Set last material to none since cross sections will need to be
        ! re-evaluated
        p % last_material = NONE

        ! Set all uvws to base level -- right now, after a collision, only the
        ! base level uvws are changed
        do j = 1, p % n_coord - 1
          if (p % coord(j + 1) % rotated) then
            ! If next level is rotated, apply rotation matrix
            p % coord(j + 1) % uvw = matmul(cells(p % coord(j) % cell) % &
                 rotation_matrix, p % coord(j) % uvw)
          else
            ! Otherwise, copy this level's direction
            p % coord(j + 1) % uvw = p % coord(j) % uvw
          end if
        end do

        ! Score flux derivative accumulators for differential tallies.
        if (active_tallies % size() > 0) call score_collision_derivative(p)
      end if

      ! If particle has too many events, display warning and kill it
      n_event = n_event + 1
      if (n_event == MAX_EVENTS) then
        if (master) call warning("Particle " // trim(to_str(p%id)) &
             &// " underwent maximum number of events.")
        p % alive = .false.
      end if

      ! Check for secondary particles if this particle is dead
      if (.not. p % alive) then
        if (p % n_secondary > 0) then
          call p % initialize_from_source(p % secondary_bank(p % n_secondary), &
                                          run_CE, energy_bin_avg)
          p % n_secondary = p % n_secondary - 1
          n_event = 0

          ! Enter new particle in particle track file
          if (p % write_track) call add_particle_track()
        else
          exit EVENT_LOOP
        end if
      end if
    end do EVENT_LOOP

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

!===============================================================================
! MOVE_PARTICLE_COORD moves particle along integration path
!===============================================================================

  subroutine move_particle_coord(coord, ds)

    type(LocalCoord), intent(inout) :: coord
    real(8),          intent(in)    :: ds    ! path length to move particle

    ! Move particle along path
    coord % xyz(1) = coord % xyz(1) + ds*coord % uvw(1);
    coord % xyz(2) = coord % xyz(2) + ds*coord % uvw(2);
    coord % xyz(3) = coord % xyz(3) + ds*coord % uvw(3);

  end subroutine move_particle_coord

!===============================================================================
! SIMPSONS_PATH_INTEGRATION uses Simpson's rule to integrate the total cross
! section along the particle path
!===============================================================================

  subroutine simpsons_path_integration(p, optical_depth, xs_t, distance)

    type(Particle), intent(inout) :: p
    real(8),        intent(inout) :: optical_depth
    real(8),        intent(inout) :: xs_t(:)
    real(8),        intent(in)    :: distance

    integer           :: i
    real(8)           :: ds       ! differential path length
    real(8)           :: new_temp ! temperature at current position
    character(len=90) :: format   ! format of the dbg output
    character(len=90) :: fname    ! filename to print total cross
    type(LocalCoord)  :: coord    ! coordinate level of continuous transport

    ! We need to get the lowest coordinate level. This is an assumption of the
    ! method.
    coord = p % coord(p % n_coord)

    ! Calculate differential path length
    ds = distance / real(num_intervals,8)

    do i = 1, num_intervals + 1
      ! Recalculate the cross section
      call calculate_xs(p)

      ! Save the total cross section
      xs_t(i) = material_xs % total

      ! Move particle along path
      call move_particle_coord(p % coord(p % n_coord), ds)
    end do

    p % coord(p % n_coord) = coord

    optical_depth = ZERO

    do i = 1, num_intervals - 1, 2
      ! Accumulate integral
      optical_depth = optical_depth + &
           TWO * ds / 6.0_8 * (xs_t(i) + FOUR * xs_t(i+1) + xs_t(i+2))
    end do

  end subroutine simpsons_path_integration

!===============================================================================
! ESTIMATE_FLIGHT_DISTANCE
!===============================================================================

  subroutine estimate_flight_distance(xs_t, distance, tau_hat, s)

    real(8), intent(inout) :: s     ! estimated flight distance
    real(8), intent(in) :: distance ! total distance to integrate
    real(8), intent(in) :: xs_t(:)  ! total cross sections at each point on path
    real(8), intent(in) :: tau_hat  ! sampled optical depth

    integer :: i
    logical :: tau_overrun   ! search for the correct bin failed?
    real(8) :: ds            ! spacing between points, assumed constant
    real(8) :: delta_tau_hat ! delta tau hat that we are solving for in the bin
    real(8) :: optical_depth
    real(8) :: a, b, c, d    ! polynomial expansion values
    real(8) :: m             ! slope of a linear guess for quadratic polynomial
    real(8) :: dds           ! differential length in the particular pin

    ! Calculate differential path length
    ds = distance / real(num_intervals,8)

    ! Set variables for loop
    i = 1
    optical_depth = ZERO
    tau_overrun = .false.

    ! Find the index of the cell that goes over the sampled path integration
    do while(optical_depth <= tau_hat)
      if (i > (num_intervals - 1)) then
        tau_overrun = .true.
        optical_depth = tau_hat
      else
        optical_depth = optical_depth + TWO * ds / 6.0_8 * (xs_t(i) + FOUR * &
             xs_t(i+1) + xs_t(i+2))
        i = i + 2
      end if
    end do

    ! Subtract off the last delta optical depth that we added
    i = i - 2
    optical_depth = optical_depth - TWO * ds / 6.0_8 * (xs_t(i) + FOUR * &
         xs_t(i+1) + xs_t(i+2))

    ! Calculate the delta in the optical depth
    delta_tau_hat = tau_hat - optical_depth

    if (tau_overrun) then
      call warning("Search for the optical depth bin failed.")
      s = distance
    else
      ! Determine the coefficients for a second order expansion in that bin
      ! sigma_t = ax^2 + bx + c
      ! Note that we are shifting the portion of the curve so that the first
      ! point lies at zero. This makes integration much easier.
      c = xs_t(i)
      b = (TWO * xs_t(i+1) - HALF * xs_t(i+2) - 1.5_8 * c) / ds
      a = (xs_t(i+2) - c - TWO * ds * b) / (FOUR * ds * ds)

      ! Define quantities for solution of the cubic equation
      b = b / TWO
      a = a / THREE
      d = -delta_tau_hat

      ! If the polynomial is not truly cubic, the soluiton will blow up.
      ! Choose between integrated polynomial of cubic, quadratic, and linear
      if (abs(a) > 1.e-10) then
        call get_cubic_root(a, b, c, d, ZERO, TWO*ds, dds)
        s = ds * (i - ONE) + dds

      else if (abs(b) > 1.e-10) then
        ! Get the best guest for root based on a linear fit
        m = (xs_t(i+2) - xs_t(i)) / (TWO * ds)
        call get_quadratic_root(HALF*m, c, d, ZERO, TWO*ds, dds)
        s = ds * (i - ONE) + dds

      else
        ! This is the case where the cross section is constant.
        s = ds * (i - ONE) - d / c

        ! Put some checks in for now
        if (s < (ds * (i - ONE)) .or. s > (ds * (i + ONE)) ) then
          call fatal_error("Estimated flight distance for the constant cross &
               &section case lies out of bounds.")
        end if
      end if
    end if

    if (s > distance .or. s < ZERO) then
      call fatal_error("The estimated flight distance is greater than the &
           &distance to boundary or less than zero.")
    endif

  end subroutine estimate_flight_distance

!===============================================================================
! GET_QUADTRATIC_ROOT returns the first quadratic root in the interval
! [lower_b, upper_b]
!===============================================================================

  subroutine get_quadratic_root(a, b, c, lower_b, upper_b, root)

    real(8), intent(inout) :: root    ! root of the equation
    real(8), intent(in)    :: a, b, c ! coefficients in ax^2 + bx + c
    real(8), intent(in)    :: lower_b ! lower bound of search
    real(8), intent(in)    :: upper_b ! upper bound of search

    real(8) :: q              ! for numeric solution
    real(8) :: c_1, c_2       ! the two roots
    real(8) :: inner_sqrt     ! temporary variable to check if real
    real(8) :: e1, e2, e3, e4 ! errors when root goes over bound
    real(8) :: min_error      ! min error used to determine which root should
                              ! be assigned

    inner_sqrt = b * b - FOUR * a * c

    if (inner_sqrt < ZERO) then
      call fatal_error("Non-real roots for quadratic equation.")
    else
      q = -ONE / TWO * (b + sign(ONE,b) * sqrt(inner_sqrt))
      c_1 = q / a
      c_2 = c / q

      if (c_1 >= lower_b .and. c_1 <= upper_b) then
        root = c_1
      else if (c_2 >= lower_b .and. c_2 <= upper_b) then
        root = c_2
      else
        e1 = abs(c_1 - lower_b)
        e2 = abs(c_1 - upper_b)
        e3 = abs(c_2 - lower_b)
        e4 = abs(c_2 - upper_b)
        min_error = min(e1, e2, e3, e4)
        if (abs(e1) == min_error .or. abs(e3) == min_error) then
          root = lower_b
        else if (abs(e2) == min_error .or. abs(e4) == min_error) then
          root = upper_b
        else
          call fatal_error("Invalid case for handling quadratic bounds overload.")
        end if

        call warning("Quadratic roots are not within 1.e-5 of bounds")
      end if
    end if

  end subroutine get_quadratic_root

!===============================================================================
! GET_CUBIC_ROOT
!===============================================================================

  subroutine get_cubic_root(a, b, c, d, lower_b, upper_b, root)

    real(8), intent(inout) :: root    ! root of the equation
    real(8), intent(in) :: a, b, c, d ! coefficients in ax^3 + bx^2 + cx + d = 0
    real(8), intent(in) :: lower_b    ! lower bound of search
    real(8), intent(in) :: upper_b    ! upper bound of search

    real(8) :: aa, bb, cc          ! coefficients in x^3 + aa x^2 + bb x + cc = 0
    real(8) :: q, r, dd, ee, theta ! temporary variables needed in root finding

    ! Calculate coefficients as they appear in the numerical recipes book
    aa = b / a
    bb = c / a
    cc = d / a

    q = (aa**2 - THREE * bb) / 9.0_8
    r = (TWO * aa**3 - 9.0_8*aa*bb + 27.0_8*cc) / 54.0_8

    if (r**2 < q**3) then
       theta = acos(r / sqrt(q**3))
       root = -TWO * sqrt(q) * cos(theta / THREE) - aa / THREE
       if (root < lower_b .or. root > upper_b) then
          root = -TWO * sqrt(q) * cos((theta + TWO * PI) / THREE) - aa / THREE
          if (root < lower_b .or. root > upper_b) then
             root = -TWO * sqrt(q) * cos((theta - TWO * PI) / THREE) - aa / THREE
             if (root < lower_b .or. root > upper_b) then
                if (root < lower_b) then
                   root = lower_b
                else
                   root = upper_b
                end if
                call warning("Acceptable cubic root not found.")
             end if
          end if
       end if
    else
       if (r >= ZERO) then
          dd = -(abs(r) + sqrt(abs(r**2 - q**3)))**(ONE / THREE)
       else
          dd = (abs(r) + sqrt(abs(r**2 - q**3)))**(ONE / THREE)
       end if

       if (dd == ZERO) then
          ee = ZERO
       else
          ee = q / dd
       end if

       root = (dd + ee) - aa / THREE
       if (root < lower_b .or. root > upper_b) then
          if (root < lower_b) then
             root = lower_b
          else
             root = upper_b
          endif
          call warning("Acceptable cubic root not found.")
       end if
    end if

  end subroutine get_cubic_root

end module tracking
