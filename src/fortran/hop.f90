program xyz_hop_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    integer, parameter :: LINESIZE = 512

    type frame_data
        integer :: n_atoms
        double precision :: H(3,3), Hinv(3,3)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=10), allocatable :: sym(:)
    end type frame_data

    ! === Input Parameters ===
    character(len=512) :: XYZ_FNAME
    character(len=8)   :: SPECIES
    integer            :: HOP_WINDOW_SIZE
    double precision   :: TIME_STEP_PS

    ! === Data Storage ===
    type(frame_data), allocatable :: trajectory(:)
    integer :: num_frames, num_valid_times, n_species
    double precision, allocatable :: hop_values_matrix(:,:)
    character(len=512) :: out_fname
    integer :: i

    ! === Read Input File ===
    open(unit=99, file='hop_xyz.inp', status='old', action='read')
    read(99, '(A)') XYZ_FNAME
    read(99, '(A)') SPECIES
    read(99, *) HOP_WINDOW_SIZE
    read(99, *) TIME_STEP_PS
    close(99)

    ! === 1. Read Trajectory Data ===
    call read_trajectory_from_xyz(XYZ_FNAME, trajectory, num_frames)

    ! === 2. Calculate Hop Function using the C code's logic ===
    call calculate_hop(trajectory, num_frames, SPECIES, HOP_WINDOW_SIZE, &
                                hop_values_matrix, num_valid_times, n_species)

    if (num_valid_times <= 0) then
        print *, 'WARNING: No valid hop data calculated. Stopping.'
        stop
    end if

    ! === 3. Write Output File ===
    write(out_fname, '(A,A,A)') './hop_per_atom_', trim(SPECIES), '.dat'
    open(unit=50, file=trim(out_fname), status='replace')
    write(50, '(A,A)') '# Per-atom hop function h_i(t) for species: ', trim(SPECIES)
    write(50, '(A)') '# CALCULATION METHOD: C-code logic (moving reference point)'
    write(50, '(A,I0,A,F10.4,A)') '# Window size = ', HOP_WINDOW_SIZE, ' frames, Timestep = ', TIME_STEP_PS, ' ps'
    block
        character(len=:), allocatable :: header_str, fmt_str
        header_str = '# Time (ps)      '
        do i = 1, n_species
            write(fmt_str, '(A,I0,A)') 'h_atom_', i, ' (arb. units)'
            header_str = trim(header_str) // ' ' // adjustl(fmt_str)
        end do
        write(50, '(A)') trim(header_str)
    end block
    do i = 1, num_valid_times
        block
            double precision :: time_val
            time_val = (dble(i) + dble(HOP_WINDOW_SIZE) - 1.0d0) * TIME_STEP_PS
            write(50, '(F17.8, tr1)', advance='no') time_val
            write(50, '(1000(F20.12))') hop_values_matrix(:, i)
        end block
    end do
    close(50)
    print *, '--- Analysis complete. Output saved to: ', trim(out_fname), ' ---'

    ! Deallocate memory
    call deallocate_traj(trajectory, num_frames)
    if (allocated(hop_values_matrix)) deallocate(hop_values_matrix)

contains
    subroutine calculate_hop(trajectory, num_frames, species, hop_window_size, &
                                     hop_out_matrix, num_valid_times, n_species_out)
        implicit none
        type(frame_data), intent(in) :: trajectory(:)
        integer, intent(in)          :: num_frames
        character(len=*), intent(in) :: species
        integer, intent(in)          :: hop_window_size
        double precision, allocatable, intent(out) :: hop_out_matrix(:,:)
        integer, intent(out)         :: num_valid_times, n_species_out

        integer :: i, j, t, dt, p1, t_center, output_idx
        integer :: n_species, window_half
        integer, allocatable :: species_indices(:)
        double precision :: sum_dA, sum_dB, dA, dB
        double precision, allocatable :: avg_pos_A_series(:,:,:), avg_pos_B_series(:,:,:)
        double precision :: current_pos(3), ref_pos(3), dr(3)

        ! --- Initialization ---
        window_half = hop_window_size / 2 
        if (num_frames < 2 * hop_window_size) then
            write(error_unit,*) 'WARNING: Not enough frames for C-logic calculation.'
            num_valid_times = 0
            n_species_out = 0
            return
        end if

        n_species = count(trajectory(1)%sym == trim(species))
        if (n_species == 0) then
            write(error_unit,*) 'WARNING: Species not found: ', trim(species)
            num_valid_times = 0; n_species_out = 0; return
        end if
        n_species_out = n_species
        allocate(species_indices(n_species))
        species_indices = pack([(i, i=1,trajectory(1)%n_atoms)], trajectory(1)%sym == trim(species))

          print *, "Phase 1: Pre-calculating moving average position series..."
        allocate(avg_pos_A_series(3, n_species, num_frames))
        allocate(avg_pos_B_series(3, n_species, num_frames))
        avg_pos_A_series = 0.0d0
        avg_pos_B_series = 0.0d0

        do j = 1, n_species
            p1 = species_indices(j)
            ! Loop over all possible center times for the moving averages
            do t = window_half, num_frames - window_half + 1
                ! Avg A: average of [t - window_half + 1, t]
                call get_avg_pbc(trajectory, p1, t - window_half + 1, t, avg_pos_A_series(:, j, t))
                ! Avg B: average of [t, t + window_half - 1]
                call get_avg_pbc(trajectory, p1, t, t + window_half - 1, avg_pos_B_series(:, j, t))
            end do
        end do
        
        print *, "Phase 2: Calculating hop values with moving reference..."
        num_valid_times = num_frames - (2 * hop_window_size) + 2
        allocate(hop_out_matrix(n_species, num_valid_times))
        hop_out_matrix = 0.0d0
        output_idx = 0

        do t_center = hop_window_size, num_frames - hop_window_size + 1
            output_idx = output_idx + 1
            do j = 1, n_species
                p1 = species_indices(j)
                sum_dA = 0.0d0
                sum_dB = 0.0d0

                do dt = 0, window_half - 1
                    ! --- dA calculation part ---
                    ! Compare pos at (t_center-dt) with avg_B at (t_center-dt)
                    current_pos = [trajectory(t_center-dt)%x(p1), trajectory(t_center-dt)%y(p1), trajectory(t_center-dt)%z(p1)]
                    ref_pos = avg_pos_B_series(:, j, t_center-dt)
                    dr = current_pos - ref_pos
                    call min_image_dr(trajectory(t_center-dt)%H, trajectory(t_center-dt)%Hinv, dr(1), dr(2), dr(3))
                    sum_dA = sum_dA + dot_product(dr, dr)

                    ! --- dB calculation part ---
                    ! Compare pos at (t_center+dt) with avg_A at (t_center+dt)i
                    current_pos = [trajectory(t_center+dt)%x(p1), trajectory(t_center+dt)%y(p1), trajectory(t_center+dt)%z(p1)]
                    ref_pos = avg_pos_A_series(:, j, t_center+dt)
                    dr = current_pos - ref_pos
                    call min_image_dr(trajectory(t_center+dt)%H, trajectory(t_center+dt)%Hinv, dr(1), dr(2), dr(3))
                    sum_dB = sum_dB + dot_product(dr, dr)
                end do

                dA = sum_dA / dble(window_half)
                dB = sum_dB / dble(window_half)

                if (dA < 0.0d0) dA = 0.0d0
                if (dB < 0.0d0) dB = 0.0d0
                hop_out_matrix(j, output_idx) = sqrt(dA * dB)
            end do
        end do

        deallocate(species_indices, avg_pos_A_series, avg_pos_B_series)

    end subroutine calculate_hop

    ! ==============================================================================
    ! Helper and File Reading Subroutines
    ! ==============================================================================
    
    subroutine read_trajectory_from_xyz(xyz_fname, trajectory, num_frames)
        implicit none
        character(len=*), intent(in)              :: xyz_fname
        type(frame_data), allocatable, intent(out) :: trajectory(:)
        integer, intent(out)                       :: num_frames

        integer :: unit_num, i, j, stat, current_natoms
        character(len=LINESIZE) :: header_line, lattice_str
        integer :: p1, p2

        open(newunit=unit_num, file=trim(xyz_fname), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write(error_unit, '(A,A)') 'FATAL: Could not open file: ', trim(xyz_fname)
            stop
        end if

        num_frames = 0
        do
            read(unit_num, *, iostat=stat) current_natoms
            if (stat /= 0) exit
            read(unit_num, '(A)') header_line
            do i = 1, current_natoms
                read(unit_num, '(A)') header_line
            end do
            num_frames = num_frames + 1
        end do
        rewind(unit_num)

        if (num_frames == 0) then
            write(error_unit,*) 'FATAL: No frames found in file ', trim(xyz_fname)
            close(unit_num); stop
        end if

        allocate(trajectory(num_frames), stat=stat)
        if (stat /= 0) stop 'Allocation error for trajectory array'

        do i = 1, num_frames
            read(unit_num, *) trajectory(i)%n_atoms
            read(unit_num, '(A)') header_line
            p1 = index(header_line, 'Lattice="')
            if (p1 > 0) then
                p2 = index(header_line(p1+9:), '"')
                if (p2 > 0) then
                    lattice_str = header_line(p1+9 : (p1+8)+p2-1)
                    call read_lattice_9(lattice_str, trajectory(i)%H)
                    call invert3x3(trajectory(i)%H, trajectory(i)%Hinv)
                else
                    write(error_unit,*) 'FATAL: Malformed Lattice string in frame ', i; stop
                end if
            else
                write(error_unit,*) 'FATAL: Lattice information not found in frame ', i; stop
            end if
            allocate(trajectory(i)%x(trajectory(i)%n_atoms), &
                     trajectory(i)%y(trajectory(i)%n_atoms), &
                     trajectory(i)%z(trajectory(i)%n_atoms), &
                     trajectory(i)%sym(trajectory(i)%n_atoms), stat=stat)
            if (stat /= 0) stop 'Allocation error for atom arrays in frame'
            do j = 1, trajectory(i)%n_atoms
                read(unit_num, *) trajectory(i)%sym(j), trajectory(i)%x(j), trajectory(i)%y(j), trajectory(i)%z(j)
            end do
        end do
        close(unit_num)
        print *, 'Successfully read ', num_frames, ' frames from ', trim(xyz_fname)
    end subroutine read_trajectory_from_xyz

    subroutine get_avg_pbc(trajectory, atom_idx, t_start, t_end, avg_pos)
        implicit none
        type(frame_data), intent(in) :: trajectory(:)
        integer, intent(in)           :: atom_idx, t_start, t_end
        double precision, intent(out) :: avg_pos(3)
        double precision :: sum_pos(3), current_pos(3), prev_pos(3), dr(3)
        integer :: t, n_points
        
        n_points = t_end - t_start + 1
        if (n_points <= 0) then; avg_pos = 0.0d0; return; end if
        
        prev_pos = [trajectory(t_start)%x(atom_idx), trajectory(t_start)%y(atom_idx), trajectory(t_start)%z(atom_idx)]
        sum_pos = prev_pos
        
        do t = t_start + 1, t_end
            current_pos = [trajectory(t)%x(atom_idx), trajectory(t)%y(atom_idx), trajectory(t)%z(atom_idx)]
            dr = current_pos - prev_pos
            call min_image_dr(trajectory(t)%H, trajectory(t)%Hinv, dr(1), dr(2), dr(3))
            current_pos = prev_pos + dr
            sum_pos = sum_pos + current_pos
            prev_pos = current_pos
        end do
        avg_pos = sum_pos / dble(n_points)
    end subroutine get_avg_pbc

    subroutine deallocate_traj(trajectory, num_frames)
        implicit none
        type(frame_data), allocatable, intent(inout) :: trajectory(:)
        integer, intent(in) :: num_frames
        integer :: i
        if (.not. allocated(trajectory)) return
        do i = 1, size(trajectory)
            if (allocated(trajectory(i)%x)) deallocate(trajectory(i)%x, trajectory(i)%y, trajectory(i)%z, trajectory(i)%sym)
        end do
        deallocate(trajectory)
    end subroutine deallocate_traj

    subroutine read_lattice_9(str, H)
        character(len=*), intent(in)  :: str
        double precision, intent(out) :: H(3,3)
        double precision :: a(9)
        integer :: ios
        read(str, *, iostat=ios) a
        if (ios /= 0) stop 'FATAL: cannot parse 9 lattice numbers.'
        H = reshape(a, [3,3], order=[2,1])
    end subroutine read_lattice_9

    subroutine invert3x3(A, Ainv)
        double precision, intent(in)  :: A(3,3)
        double precision, intent(out) :: Ainv(3,3)
        double precision :: detA
        detA = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        if (abs(detA) < 1.d-12) stop 'FATAL: singular cell matrix.'
        Ainv(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2)) / detA
        Ainv(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2)) / detA
        Ainv(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2)) / detA
        Ainv(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1)) / detA
        Ainv(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1)) / detA
        Ainv(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1)) / detA
        Ainv(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1)) / detA
        Ainv(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1)) / detA
        Ainv(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / detA
    end subroutine invert3x3

    subroutine min_image_dr(H, Hinv, dx, dy, dz)
        double precision, intent(in)     :: H(3,3), Hinv(3,3)
        double precision, intent(inout)  :: dx, dy, dz
        double precision :: s(3), dr(3)
        dr = [dx, dy, dz]
        s = matmul(Hinv, dr)
        s = s - anint(s)
        dr = matmul(H, s)
        dx = dr(1); dy = dr(2); dz = dr(3)
    end subroutine min_image_dr

end program xyz_hop_analyzer