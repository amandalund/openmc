module dict_header

!===============================================================================
! DICT_HEADER module
!
! This module provides an implementation of a dictionary that has (key,value)
! pairs. This data structure is used to provide lookup features, e.g. cells and
! surfaces by name.
!
! The original version was roughly based on capabilities in the 'flibs' open
! source package. However, it was rewritten from scratch so that it could be
! used stand-alone without relying on the implementation of lists. As with
! lists, it was considered writing a single dictionary used unlimited
! polymorphism, but again compiler support is spotty and doesn't always prevent
! duplication of code.
!===============================================================================

  implicit none

  integer, parameter          :: EMPTY           = -huge(0)
  integer, parameter          :: KEY_CHAR_LENGTH = 255
  integer, parameter, private :: MIN_SIZE        = 8
  integer, parameter, private :: HASH_MULTIPLIER = 31
  real(8), parameter, private :: MAX_LOAD_FACTOR = 1.

!===============================================================================
! DICTENTRY* contains (key,value) pairs and a pointer to the next (key,value)
! pair
!===============================================================================

  type DictEntryCI
    type(DictEntryCI), pointer :: next => null()
    character(len=KEY_CHAR_LENGTH) :: key
    integer :: value
  end type DictEntryCI

  type DictEntryII
    type(DictEntryII), pointer :: next => null()
    integer :: key
    integer :: value
  end type DictEntryII

!===============================================================================
! HASHLIST* types contain a single pointer to a linked list of (key,value)
! pairs. This type is necesssary so that the Dict types can be dynamically
! allocated.
!===============================================================================

  type, private :: HashListCI
    type(DictEntryCI), pointer :: list => null()
  end type HashListCI

  type, private :: HashListII
    type(DictEntryII), pointer :: list => null()
  end type HashListII

!===============================================================================
! DICT* is a dictionary of (key,value) pairs with convenience methods as
! type-bound procedures. DictCharInt has character(*) keys and integer values,
! and DictIntInt has integer keys and values.
!===============================================================================

  type, public :: DictCharInt
    integer, private :: entries = 0
    type(HashListCI), allocatable, private :: table(:)
  contains
    procedure :: add => add_ci
    procedure :: get => get_ci
    procedure :: has => has_ci
    procedure :: keys => keys_ci
    procedure :: size => size_ci
    procedure :: clear => clear_ci
    procedure, private :: get_elem => get_elem_ci
    procedure, private :: resize => resize_ci
    procedure, private :: hash => hash_ci
  end type DictCharInt

  type, public :: DictIntInt
    integer, private :: entries = 0
    type(HashListII), allocatable, private :: table(:)
  contains
    procedure :: add => add_ii
    procedure :: get => get_ii
    procedure :: has => has_ii
    procedure :: keys => keys_ii
    procedure :: size => size_ii
    procedure :: clear => clear_ii
    procedure, private :: get_elem => get_elem_ii
    procedure, private :: resize => resize_ii
    procedure, private :: hash => hash_ii
  end type DictIntInt

contains

!===============================================================================
! ADD adds a (key,value) entry to a dictionary. If the key is already in the
! dictionary, the value is replaced by the new specified value.
!===============================================================================

  subroutine add_ci(this, key, value)

    class(DictCharInt) :: this
    character(*), intent(in) :: key
    integer,      intent(in) :: value

    integer :: hash
    type(DictEntryCI), pointer :: elem => null()
    type(DictEntryCI), pointer :: new_elem => null()

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
    else if (real(this % entries+1, 8)/size(this % table) > MAX_LOAD_FACTOR) then
      call this % resize()
    end if

    hash = this % hash(key)

    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) then
        elem % value = value
        return
      end if
      elem => elem % next
    end do

    ! Create new element
    allocate(new_elem)
    new_elem % key = key
    new_elem % value = value

    ! Add element to front of list
    new_elem % next => this % table(hash) % list
    this % table(hash) % list => new_elem
    this % entries = this % entries + 1

  end subroutine add_ci

  subroutine add_ii(this, key, value)

    class(DictIntInt) :: this
    integer, intent(in) :: key
    integer, intent(in) :: value

    integer :: hash
    type(DictEntryII), pointer :: elem => null()
    type(DictEntryII), pointer :: new_elem => null()

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
    else if (real(this % entries+1, 8)/size(this % table) > MAX_LOAD_FACTOR) then
      call this % resize()
    end if

    hash = this % hash(key)

    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) then
        elem % value = value
        return
      end if
      elem => elem % next
    end do

    ! Create new element
    allocate(new_elem)
    new_elem % key = key
    new_elem % value = value

    ! Add element to front of list
    new_elem % next => this % table(hash) % list
    this % table(hash) % list => new_elem
    this % entries = this % entries + 1

  end subroutine add_ii

!===============================================================================
! GET_ELEM returns a pointer to the (key,value) pair for a given key. This
! method is private.
!===============================================================================

  function get_elem_ci(this, key) result(elem)

    class(DictCharInt)         :: this
    character(*), intent(in)   :: key
    type(DictEntryCI), pointer :: elem

    integer :: hash

    ! Check for dictionary not being allocated
    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
    end if

    hash = this % hash(key)
    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) exit
      elem => elem % next
    end do

  end function get_elem_ci

  function get_elem_ii(this, key) result(elem)

    class(DictIntInt)          :: this
    integer, intent(in)        :: key
    type(DictEntryII), pointer :: elem

    integer :: hash

    ! Check for dictionary not being allocated
    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
    end if

    hash = this % hash(key)
    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) exit
      elem => elem % next
    end do

  end function get_elem_ii

!===============================================================================
! RESIZE allocates a new hash table to accomodate the number of entries and
! reinserts all of the entries into the new table. This method is private.
!===============================================================================

  subroutine resize_ci(this)

    class(DictCharInt) :: this

    integer :: i
    integer :: hash
    integer :: new_size
    type(HashListCI), allocatable :: table(:)
    type(DictEntryCI), pointer :: elem
    type(DictEntryCI), pointer :: next

    new_size = 2 * size(this % table)

    call move_alloc(this % table, table)
    allocate(this % table(new_size))

    ! Rehash each entry into the new table
    do i = 1, size(table)
      elem => table(i) % list
      do while (associated(elem))
        hash = this % hash(elem % key)
        next => elem % next
        elem % next => this % table(hash) % list
        this % table(hash) % list => elem
        elem => next
      end do
    end do
    deallocate(table)

  end subroutine resize_ci

  subroutine resize_ii(this)

    class(DictIntInt) :: this

    integer :: i
    integer :: hash
    integer :: new_size
    type(HashListII), allocatable :: table(:)
    type(DictEntryII), pointer :: elem
    type(DictEntryII), pointer :: next

    new_size = 2 * size(this % table)

    call move_alloc(this % table, table)
    allocate(this % table(new_size))

    ! Rehash each entry into the new table
    do i = 1, size(table)
      elem => table(i) % list
      do while (associated(elem))
        hash = this % hash(elem % key)
        next => elem % next
        elem % next => this % table(hash) % list
        this % table(hash) % list => elem
        elem => next
      end do
    end do
    deallocate(table)

  end subroutine resize_ii

!===============================================================================
! GET returns the value matching a given key. If the dictionary does not
! contain the key, the value EMPTY is returned.
!===============================================================================

  function get_ci(this, key) result(value)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer                  :: value

    type(DictEntryCI), pointer :: elem

    elem => this % get_elem(key)

    if (associated(elem)) then
      value = elem % value
    else
      value = EMPTY
    end if

  end function get_ci

  function get_ii(this, key) result(value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer             :: value

    type(DictEntryII), pointer :: elem

    elem => this % get_elem(key)

    if (associated(elem)) then
      value = elem % value
    else
      value = EMPTY
    end if

  end function get_ii

!===============================================================================
! HAS determines whether a dictionary has a (key,value) pair with a
! given key.
!===============================================================================

  function has_ci(this, key) result(has)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    logical                  :: has

    type(DictEntryCI), pointer :: elem

    elem => this % get_elem(key)
    has = associated(elem)

  end function has_ci

  function has_ii(this, key) result(has)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    logical             :: has

    type(DictEntryII), pointer  :: elem

    elem => this % get_elem(key)
    has = associated(elem)

  end function has_ii

!===============================================================================
! HASH returns the hash value for a given key
!===============================================================================

  function hash_ci(this, key) result(hash)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer                  :: hash

    integer :: i

    hash = 0

    do i = 1, len_trim(key)
      hash = HASH_MULTIPLIER * hash + ichar(key(i:i))
    end do

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    hash = 1 + mod(abs(hash-1), size(this % table))

  end function hash_ci

  function hash_ii(this, key) result(hash)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer             :: hash

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    hash = 1 + mod(abs(key-1), size(this % table))

  end function hash_ii

!===============================================================================
! KEYS returns a pointer to a linked list of all (key,value) pairs
!===============================================================================

  function keys_ci(this) result(keys)
    class(DictCharInt)         :: this
    type(DictEntryCI), pointer :: keys

    integer :: i
    type(DictEntryCI), pointer :: current => null()
    type(DictEntryCI), pointer :: elem => null()

    keys => null()

    do i = 1, size(this % table)
      ! Get pointer to start of bucket i
      elem => this % table(i) % list

      do while (associated(elem))
        ! Allocate (key,value) pair
        if (.not. associated(keys)) then
          allocate(keys)
          current => keys
        else
          allocate(current % next)
          current => current % next
        end if

        ! Copy (key,value) pair
        current % key   = elem % key
        current % value = elem % value

        ! Move to next element in bucket i
        elem => elem % next
      end do
    end do

  end function keys_ci

  function keys_ii(this) result(keys)
    class(DictIntInt)          :: this
    type(DictEntryII), pointer :: keys

    integer :: i
    type(DictEntryII), pointer :: current => null()
    type(DictEntryII), pointer :: elem => null()

    keys => null()

    do i = 1, size(this % table)
      ! Get pointer to start of bucket i
      elem => this % table(i) % list

      do while (associated(elem))
        ! Allocate (key,value) pair
        if (.not. associated(keys)) then
          allocate(keys)
          current => keys
        else
          allocate(current % next)
          current => current % next
        end if

        ! Copy (key,value) pair
        current % key   = elem % key
        current % value = elem % value

        ! Move to next element in bucket i
        elem => elem % next
      end do
    end do

  end function keys_ii

!===============================================================================
! CLEAR Deletes and deallocates the dictionary item
!===============================================================================

  subroutine clear_ci(this)

    class(DictCharInt) :: this

    integer :: i
    type(DictEntryCI), pointer :: current
    type(DictEntryCI), pointer :: next

    if (allocated(this % table)) then
      do i = 1, size(this % table)
        current => this % table(i) % list
        do while (associated(current))
          if (associated(current % next)) then
            next => current % next
          else
            nullify(next)
          end if
          deallocate(current)
          current => next
        end do
        if (associated(this % table(i) % list)) &
             nullify(this % table(i) % list)
      end do
      deallocate(this % table)
    end if

  end subroutine clear_ci

  subroutine clear_ii(this)

    class(DictIntInt) :: this

    integer :: i
    type(DictEntryII), pointer :: current
    type(DictEntryII), pointer :: next

    if (allocated(this % table)) then
      do i = 1, size(this % table)
        current => this % table(i) % list
        do while (associated(current))
          if (associated(current % next)) then
            next => current % next
          else
            nullify(next)
          end if
          deallocate(current)
          current => next
        end do
        if (associated(this % table(i) % list)) &
             nullify(this % table(i) % list)
      end do
      deallocate(this % table)
    end if

  end subroutine clear_ii

!===============================================================================
! SIZE returns the number of entries in the dictionary
!===============================================================================

  pure function size_ci(this) result(size)

    class(DictCharInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ci

  pure function size_ii(this) result(size)

    class(DictIntInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ii

end module dict_header
