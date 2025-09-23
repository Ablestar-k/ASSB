program ensemble_hop_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    ! Derived type for storing data of a single trajectory frame
    type frame_data
        integer :: n_atoms
        double precision :: H(3,3), Hinv(3,3)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=10), allocatable :: sym(:)
    end type frame_data

    character(len=256) :: SIMULATION_NAME
    integer :: NUM_ENSEMBLES
    character(len=8) :: SPECIES
    integer :: HOP_WINDOW_SIZE  ! Total number of frames for the two windows (A and B), must be even
    double precision :: TIME_STEP_PS

    integer :: i, j, t, stat, num_valid_times
    character(len=512) :: xyz_fname, dir_part, file_part, out_fname
    
    ! === Data storage ===
    double precision, allocatable :: all_hop_values(:,:)
    double precision, allocatable :: current_hop_t(:)
    double precision, allocatable :: mean_hop_arr(:), std_hop_arr(:)

    ! === Read Input File ===
    open(unit=99, file='hop_all.inp', status='old', action='read')
    read(99, *) SIMULATION_NAME, NUM_ENSEMBLES
    read(99, *) SPECIES
    read(99, *) HOP_WINDOW_SIZE
    read(99, *) TIME_STEP_PS
    close(99)

    if (mod(HOP_WINDOW_SIZE, 2) /= 0) then
        write(error_unit,*) 'FATAL: HOP_WINDOW_SIZE must be an even integer.'
        stop
    end if

    print *, '--- Starting Ensemble Hop Function h(t) Analysis ---'
    print *, 'Simulation Series Name : ', trim(SIMULATION_NAME)
    print *, 'Number of Ensembles    : ', NUM_ENSEMBLES
    print *, 'Species to Analyze     : ', trim(SPECIES)
    print *, 'Hop Function Window    : ', HOP_WINDOW_SIZE, ' frames'
    print *, 'Time Step              : ', TIME_STEP_PS, ' ps'
    print *, '----------------------------------------------------'


    ! =====================================
    !  Main loop over each ensembles
    ! =====================================
    do i = 1, NUM_ENSEMBLES
        write(dir_part, '(A,A,A,I0)') '../dump/dump_', trim(SIMULATION_NAME), '_', i
        write(file_part, '(I0,A)') i, '_product.xyz'
        xyz_fname = trim(dir_part) // '/' // trim(file_part)
        print *, 'Processing Ensemble ', i, ': ', trim(xyz_fname)

        if (i == 1) then
            call calculate_hop_single(xyz_fname, SPECIES, HOP_WINDOW_SIZE, current_hop_t, num_valid_times)
            if (num_valid_times > 0) then
                allocate(all_hop_values(num_valid_times, NUM_ENSEMBLES), stat=stat)
                if (stat /= 0) stop 'Allocation error for all_hop_values'
                all_hop_values = 0.0d0
                all_hop_values(:, i) = current_hop_t(:)
                deallocate(current_hop_t)

            else
                print *, 'WARNING: No valid hop data from first ensemble. Stopping.'
                stop
            end if

        else
            call calculate_hop_single(xyz_fname, SPECIES, HOP_WINDOW_SIZE, current_hop_t, num_valid_times)
            if (size(current_hop_t) == size(all_hop_values, 1)) then
                 all_hop_values(:, i) = current_hop_t(:)
            else
                 print *, 'WARNING: Trajectory length mismatch in ensemble', i, '. Skipping.'
                 all_hop_values(:, i) = -1.0d0 
            end if
            if (allocated(current_hop_t)) deallocate(current_hop_t)
        end if
    end do

    print *, ''
    print *, '--- Aggregating results and calculating statistics ---'

    allocate(mean_hop_arr(num_valid_times), std_hop_arr(num_valid_times))

    ! === Statistics calculation ===
    do t = 1, num_valid_times
        block_stats: block
            double precision :: mean_val, sum_sq, std_val
            double precision, allocatable :: valid_data(:)
            integer :: valid_count
            
            valid_count = count(all_hop_values(t, :) >= 0.0d0)
            if (valid_count > 0) then
                allocate(valid_data(valid_count))
                valid_data = pack(all_hop_values(t,:), all_hop_values(t,:) >= 0.0d0)
                mean_val = sum(valid_data) / dble(valid_count)
                
                if (valid_count > 1) then
                    sum_sq = sum((valid_data - mean_val)**2)
                    std_val = sqrt(sum_sq / dble(valid_count - 1))
                else
                    std_val = 0.0d0
                end if
                deallocate(valid_data)

            else
                mean_val = 0.0d0
                std_val = 0.0d0
            end if

            mean_hop_arr(t) = mean_val
            std_hop_arr(t) = std_val
        end block block_stats
    end do

    ! === Write Output File ===
    write(out_fname, '(A,A,A)') '../result/hop_avg_', trim(SPECIES), '.dat'
    open(unit=50, file=trim(out_fname), status='replace')
    
    write(50, '(A,A)') '# Ensemble-averaged hop function h(t) for species: ', trim(SPECIES)
    write(50, '(A,I0,A,F10.4,A)') '# Window size = ', HOP_WINDOW_SIZE, ' frames, Timestep = ', TIME_STEP_PS, ' ps'
    write(50, '(A23, A25, A25)') '# Time (ps)', 'h(t)_mean (Angstrom^2)', 'h(t)_std (Angstrom^2)'
    
    do t = 1, num_valid_times
        block_write: block
            double precision :: time_val
            time_val = (dble(t) + dble(HOP_WINDOW_SIZE)/2.0d0) * TIME_STEP_PS
            write(50, '(F24.8, F25.12, F25.12)') time_val, mean_hop_arr(t), std_hop_arr(t)
        end block block_write
    end do
    close(50)

    print *, '--- Analysis complete. Output saved to: ', trim(out_fname), ' ---'

    deallocate(all_hop_values, mean_hop_arr, std_hop_arr)

contains

! ==============================================================================
!  Calculates the species-averaged hop function h(t) for a single trajectory
! ==============================================================================
subroutine calculate_hop_single(xyz_fname, species, hop_window_size, hop_out, num_valid_times)
    implicit none
    character(len=*), intent(in)  :: xyz_fname, species
    integer, intent(in)           :: hop_window_size
    double precision, allocatable, intent(out) :: hop_out(:)
    integer, intent(out)          :: num_valid_times

    integer :: ios, i, j, k, p1, p2, num_frames, t_idx, n_species
    character(len=1024) :: header_line, aline, lattice_str
    type(frame_data), allocatable :: trajectory(:)
    integer, allocatable :: species_indices(:)
    double precision, allocatable :: atom_hop_values(:,:)
    integer :: window_half

    ! === Initialization ===
    ! CORRECTED: Removed the problematic initial allocation `hop_out = [0.0d0]`
    num_valid_times = 0

    open(unit=10, file=trim(xyz_fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'WARNING: Could not open file, skipping: ', trim(xyz_fname)
        return
    end if

    ! --- First pass: count frames ---
    num_frames = 0
    do
        read(10,*, iostat=ios) p1
        if (ios /= 0) exit
        read(10,'(A)') header_line
        do i = 1, p1
            read(10,'(A)') aline
        end do
        num_frames = num_frames + 1
    end do
    rewind(10)
    
    if (num_frames < hop_window_size) then
        write(error_unit,*) 'WARNING: Not enough frames for hop calculation in ', trim(xyz_fname)
        close(10); return
    end if

    ! === Second pass: Load all trajectory data into memory ===
    allocate(trajectory(num_frames))
    do i = 1, num_frames
        read(10,*) trajectory(i)%n_atoms
        read(10,'(A)') header_line
        allocate(trajectory(i)%x(trajectory(i)%n_atoms), &
                 trajectory(i)%y(trajectory(i)%n_atoms), &
                 trajectory(i)%z(trajectory(i)%n_atoms), &
                 trajectory(i)%sym(trajectory(i)%n_atoms))

        do j = 1, trajectory(i)%n_atoms
            read(10,'(A)') aline
            read(aline,*) trajectory(i)%sym(j), trajectory(i)%x(j), trajectory(i)%y(j), trajectory(i)%z(j)
        end do

        p1 = index(header_line, 'Lattice="')
        if (p1 > 0) then
            p2 = index(header_line(p1+9:), '"')
            lattice_str = header_line(p1+9 : (p1+8)+p2-1)
            call read_lattice_9(lattice_str, trajectory(i)%H)
            call invert3x3(trajectory(i)%H, trajectory(i)%Hinv)
        else
            write(error_unit,*) 'FATAL: Lattice information not found in frame ', i
            stop
        end if
    end do
    close(10)
    
    ! === Find indices of the target species ===
    n_species = count(trajectory(1)%sym == trim(species))
    if (n_species == 0) then
        write(error_unit,*) 'WARNING: Species not found in file: ', trim(xyz_fname)
        call deallocate_traj(trajectory, num_frames)
        return
    end if

    allocate(species_indices(n_species))
    species_indices = pack([(i, i=1,trajectory(1)%n_atoms)], trajectory(1)%sym == trim(species))

    ! === Main Hop Function Calculation ===
    window_half = hop_window_size / 2
    num_valid_times = num_frames - hop_window_size + 1
    allocate(atom_hop_values(n_species, num_valid_times))
    atom_hop_values = 0.0d0

    ! Loop over each atom of the target species
    do j = 1, n_species
        p1 = species_indices(j)
        ! Loop over each time t that can be the center of a full window
        do t_idx = 1, num_valid_times
            block_hop_calc: block
                integer :: t_start, t_center
                double precision :: avg_pos_A(3), avg_pos_B(3)
                double precision :: term_A, term_B, h_val
                
                t_start = t_idx
                t_center = t_start + window_half

                ! Get average positions in window A and B
                call get_avg_pbc(trajectory, p1, t_start, t_center - 1, avg_pos_A)
                call get_avg_pbc(trajectory, p1, t_center, t_start + hop_window_size - 1, avg_pos_B)
                
                ! Calculate term A = <(r - avg_pos_B)^2>_A
                call get_msd_from_point_pbc(trajectory, p1, t_start, t_center - 1, avg_pos_B, term_A)

                ! Calculate term B = <(r - avg_pos_A)^2>_B
                call get_msd_from_point_pbc(trajectory, p1, t_center, t_start + hop_window_size - 1, avg_pos_A, term_B)
                
                h_val = sqrt(term_A * term_B)
                atom_hop_values(j, t_idx) = h_val
            end block block_hop_calc
        end do
    end do

    ! Average over all atoms for each time step
    allocate(hop_out(num_valid_times))
    do t_idx = 1, num_valid_times
        hop_out(t_idx) = sum(atom_hop_values(:, t_idx)) / dble(n_species)
    end do

    deallocate(species_indices, atom_hop_values)
    call deallocate_traj(trajectory, num_frames)

end subroutine calculate_hop_single

! ==============================================================================
! Helper Subroutines
! ==============================================================================

! --- Calculates average position in a trajectory slice with PBC ---
subroutine get_avg_pbc(trajectory, atom_idx, t_start, t_end, avg_pos)
    implicit none
    type(frame_data), intent(in) :: trajectory(:) 
    integer, intent(in)          :: atom_idx, t_start, t_end
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

! --- Calculates MSD from a fixed point over a trajectory slice with PBC ---
subroutine get_msd_from_point_pbc(trajectory, atom_idx, t_start, t_end, ref_pos, msd)
    implicit none
    type(frame_data), intent(in) :: trajectory(:)
    integer, intent(in)          :: atom_idx, t_start, t_end
    double precision, intent(in)  :: ref_pos(3)
    double precision, intent(out) :: msd
    double precision :: sum_sq_dist, dr(3)
    integer :: t, n_points
    
    n_points = t_end - t_start + 1
    if (n_points <= 0) then; msd = 0.0d0; return; end if

    sum_sq_dist = 0.0d0
    do t = t_start, t_end
        dr = [trajectory(t)%x(atom_idx), trajectory(t)%y(atom_idx), trajectory(t)%z(atom_idx)] - ref_pos
        call min_image_dr(trajectory(t)%H, trajectory(t)%Hinv, dr(1), dr(2), dr(3))
        sum_sq_dist = sum_sq_dist + dot_product(dr, dr)
    end do
    msd = sum_sq_dist / dble(n_points)
end subroutine get_msd_from_point_pbc

! --- Deallocates trajectory data ---
subroutine deallocate_traj(trajectory, num_frames)
    implicit none
    type(frame_data), allocatable, intent(inout) :: trajectory(:)
    integer, intent(in) :: num_frames
    integer :: i
    do i = 1, num_frames
        if (allocated(trajectory(i)%x)) deallocate(trajectory(i)%x, trajectory(i)%y, trajectory(i)%z, trajectory(i)%sym)
    end do
    if (allocated(trajectory)) deallocate(trajectory)
end subroutine deallocate_traj

! --- Parses 9 lattice vector components from a string ---
subroutine read_lattice_9(str, H)
    character(len=*), intent(in)  :: str
    double precision, intent(out) :: H(3,3)
    double precision :: a(9)
    integer :: ios
    read(str, *, iostat=ios) a
    if (ios /= 0) stop 'FATAL: cannot parse 9 lattice numbers.'
    H = reshape(a, [3,3], order=[2,1]) ! More concise way to assign matrix
end subroutine read_lattice_9

! --- Inverts a 3x3 matrix ---
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
    ! Convert to fractional coordinates
    s = matmul(Hinv, dr)
    ! Apply nearest image convention
    s = s - anint(s)
    ! Convert back to Cartesian coordinates
    dr = matmul(H, s)
    dx = dr(1); dy = dr(2); dz = dr(3)
end subroutine min_image_dr

end program ensemble_hop_analyzer
