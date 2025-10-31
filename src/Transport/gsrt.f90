program ensemble_vanhove_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    character(len=256) :: SIMULATION_NAME
    integer :: NUM_ENSEMBLES
    character(len=8) :: SPECIES
    double precision :: BINSIZE, MAX_R
    integer :: T_INTERVAL, N_T_DELTAS
    integer, allocatable :: T_DELTAS(:)
    integer :: i, j, k, n_bins, stat
    character(len=512) :: xyz_fname, dir_part, file_part, out_fname
    double precision, allocatable :: all_gs_values(:,:,:)
    double precision, allocatable :: current_gs(:,:)
    double precision, allocatable :: mean_gs_arr(:,:), std_gs_arr(:,:)

    ! === Read Input File ===
    open(unit=99, file='gsrt.inp', status='old', action='read')
    read(99, *) SIMULATION_NAME, NUM_ENSEMBLES
    read(99, *) SPECIES
    read(99, *) BINSIZE, MAX_R
    read(99, *) T_INTERVAL
    read(99, *) N_T_DELTAS
    allocate(T_DELTAS(N_T_DELTAS))
    read(99, *) T_DELTAS(:)
    close(99)

    n_bins = ceiling(MAX_R / BINSIZE)

    ! === Allocate Memory ===
    allocate(all_gs_values(n_bins, N_T_DELTAS, NUM_ENSEMBLES), stat=stat)
    if (stat /= 0) stop 'Allocation error for all_gs_values'
    allocate(current_gs(n_bins, N_T_DELTAS), stat=stat)
    if (stat /= 0) stop 'Allocation error for current_gs'
    all_gs_values = 0.0d0

    print *, '--- Starting Ensemble Van Hove Gs(r,t) Analysis ---'
    print *, 'Species: ', trim(SPECIES)
    print *, 'Time Deltas (frames): ', T_DELTAS

    ! =====================================
    !  Main loop over each ensemble member
    ! =====================================
    do i = 1, NUM_ENSEMBLES
        write(dir_part, '(A,A,A,I0)') '../../dump/dump_', trim(SIMULATION_NAME), '_', i
        write(file_part, '(I0,A)') i, '_product_unwrapped.xyz'
        xyz_fname = trim(dir_part) // '/' // trim(file_part)
        print *, 'Processing: ', trim(xyz_fname)

        call calculate_vanhove_single(xyz_fname, SPECIES, n_bins, MAX_R, T_DELTAS, current_gs)
        all_gs_values(:, :, i) = current_gs(:, :)
    end do

    print *, ''
    print *, '--- Aggregating results and calculating statistics ---'

    allocate(mean_gs_arr(n_bins, N_T_DELTAS), std_gs_arr(n_bins, N_T_DELTAS))

    ! === Statistics calculation ===
    do k = 1, N_T_DELTAS
        do i = 1, n_bins
            block_stats: block
                double precision :: mean_val, sum_sq, std_val
                mean_val = sum(all_gs_values(i, k, :)) / dble(NUM_ENSEMBLES)
                sum_sq = sum((all_gs_values(i, k, :) - mean_val)**2)
                
                if (NUM_ENSEMBLES > 1) then
                    std_val = sqrt(sum_sq / dble(NUM_ENSEMBLES - 1))
                else
                    std_val = 0.0d0
                end if
                mean_gs_arr(i, k) = mean_val
                std_gs_arr(i, k) = std_val
            end block block_stats
        end do
    end do

    ! === Write Output File ===
    write(out_fname, '(A,A,A)') '../../result/gsrt_ensemble_average_', trim(SPECIES), '.dat'
    open(unit=50, file=trim(out_fname), status='replace')
    
    do k = 1, N_T_DELTAS
        write(50, '(A, I0, A)') '# t_delta = ', T_DELTAS(k) * T_INTERVAL, '(ps)'
        write(50, '(A24, A20, A20)') '# r (Angstrom)', 'Gs(r,t)_mean', 'Gs(r,t)_std'
        do i = 1, n_bins
            block_write: block
                double precision :: r_center
                r_center = (dble(i) - 0.5d0) * BINSIZE
                write(50, '(F25.8, F20.8, F20.8)') r_center, mean_gs_arr(i, k), std_gs_arr(i, k)
            end block block_write
        end do
        write(50,*) ''
        write(50,*) ''
    end do
    close(50)

    print *, '--- Analysis complete. Output saved to: ', trim(out_fname), ' ---'

    deallocate(T_DELTAS, all_gs_values, current_gs, mean_gs_arr, std_gs_arr)

contains

! ==========================================================
!  Calculates Gs(r,t) for a single trajectory file
! ==========================================================
subroutine calculate_vanhove_single(xyz_fname, species, n_bins, max_r, t_deltas, gs_out)
    implicit none
    character(len=*), intent(in)  :: xyz_fname, species
    integer, intent(in)           :: n_bins
    double precision, intent(in)  :: max_r
    integer, intent(in)           :: t_deltas(:)
    double precision, intent(out) :: gs_out(n_bins, size(t_deltas))

    type frame_data
        integer :: n_atoms
        double precision :: H(3,3), Hinv(3,3)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=10), allocatable :: sym(:)
    end type frame_data
    
    integer :: ios, i, j, k, p1, p2, num_frames, t0_idx, tf_idx
    integer :: n_t, n_species
    integer :: bin_index
    double precision :: binsize, dist, pi, r, shell_vol
    double precision :: dx, dy, dz
    character(len=1024) :: header_line, aline, lattice_str
    
    type(frame_data), allocatable :: trajectory(:)
    integer(8), allocatable :: hist_sum(:,:)
    integer, allocatable :: n_origins(:)
    integer, allocatable :: species_indices(:)

    ! === Initialization ===
    pi = 4.0d0 * datan(1.0d0)
    binsize = max_r / dble(n_bins)
    n_t = size(t_deltas)
    gs_out = 0.0d0

    open(unit=10, file=trim(xyz_fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,A)') 'WARNING: Could not open file, skipping: ', trim(xyz_fname)
        return
    end if
    num_frames = 0
    rewind(10)
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
    
    if (num_frames == 0) then
        write(error_unit,*) 'WARNING: No frames found in ', trim(xyz_fname)
        close(10); return
    end if

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
            if (.not. safe_parse_atom_line(aline, trajectory(i)%sym(j), &
                trajectory(i)%x(j), trajectory(i)%y(j), trajectory(i)%z(j))) then
                write(error_unit,*) 'FATAL: Bad atom line in pre-read step.'
                stop
            end if
        end do
        p1 = index(header_line, 'Lattice="')
        p2 = index(header_line(p1+9:), '"')
        lattice_str = header_line(p1+9 : (p1+8)+p2-1)
        call read_lattice_9(lattice_str, trajectory(i)%H)
        block_inv: block
            double precision :: detH
            call invert3x3(trajectory(i)%H, trajectory(i)%Hinv, detH)
        end block block_inv
    end do
    close(10)
    
    n_species = 0
    do i = 1, trajectory(1)%n_atoms
        if (trim(trajectory(1)%sym(i)) == trim(species)) then
            n_species = n_species + 1
        end if
    end do

    if (n_species == 0) then
        write(error_unit,*) 'WARNING: Species not found in file: ', trim(xyz_fname)
        deallocate(trajectory); return
    end if
    allocate(species_indices(n_species))
    k = 0
    do i = 1, trajectory(1)%n_atoms
        if (trim(trajectory(1)%sym(i)) == trim(species)) then
            k = k + 1
            species_indices(k) = i
        end if
    end do

    allocate(hist_sum(n_bins, n_t), n_origins(n_t))
    hist_sum = 0
    n_origins = 0
    do t0_idx = 1, num_frames
        do k = 1, n_t
            tf_idx = t0_idx + t_deltas(k)
            if (tf_idx > num_frames) cycle
            n_origins(k) = n_origins(k) + 1
            do i = 1, n_species
                p1 = species_indices(i)
                dx = trajectory(tf_idx)%x(p1) - trajectory(t0_idx)%x(p1)
                dy = trajectory(tf_idx)%y(p1) - trajectory(t0_idx)%y(p1)
                dz = trajectory(tf_idx)%z(p1) - trajectory(t0_idx)%z(p1)

                ! **********************************************************************
                ! *** CRITICAL FIX 1: Removed call to min_image_dr.                   ***
                ! *** This code NOW ASSUMES the input trajectory file is UNWRAPPED,   ***
                ! *** which is required for Gs(r,t).                                  ***
                ! **********************************************************************
                ! call min_image_dr(trajectory(t0_idx)%H, trajectory(t0_idx)%Hinv, dx, dy, dz)
                
                dist = sqrt(dx*dx + dy*dy + dz*dz)
                if (dist < max_r) then
                    bin_index = ceiling(dist / binsize)
                    if (bin_index > 0 .and. bin_index <= n_bins) then
                         hist_sum(bin_index, k) = hist_sum(bin_index, k) + 1
                    end if
                end if
            end do
        end do
    end do

    do k = 1, n_t
        if (n_origins(k) == 0) cycle
        do i = 1, n_bins
            r = (dble(i) - 0.5d0) * binsize
            if (r > 1.d-8) then
                shell_vol = 4.0d0 * pi * r**2 * binsize
                gs_out(i, k) = (dble(hist_sum(i,k)) / dble(n_origins(k))) &
                             / dble(n_species) / shell_vol
            else
                gs_out(i, k) = 0.0d0
            end if
        end do
    end do

    deallocate(hist_sum, n_origins, species_indices)
    do i = 1, num_frames
        deallocate(trajectory(i)%x, trajectory(i)%y, trajectory(i)%z, trajectory(i)%sym)
    end do
    deallocate(trajectory)
end subroutine calculate_vanhove_single

! ===================================================================
! Subroutines
! ===================================================================
subroutine read_lattice_9(str, H)
    character(len=*), intent(in)  :: str
    double precision, intent(out) :: H(3,3)
    double precision :: a(9)
    integer :: ios
    a = 0.d0
    read(str, *, iostat=ios) a
    if (ios /= 0) stop 'FATAL: cannot parse 9 lattice numbers.'
    H(1,1)=a(1); H(1,2)=a(2); H(1,3)=a(3)
    H(2,1)=a(4); H(2,2)=a(5); H(2,3)=a(6)
    H(3,1)=a(7); H(3,2)=a(8); H(3,3)=a(9)
end subroutine read_lattice_9

subroutine invert3x3(A, Ainv, detA)
    double precision, intent(in)  :: A(3,3)
    double precision, intent(out) :: Ainv(3,3), detA
    double precision :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    c11 =  A(2,2)*A(3,3) - A(2,3)*A(3,2)
    c12 =-(A(2,1)*A(3,3) - A(2,3)*A(3,1))
    c13 =  A(2,1)*A(3,2) - A(2,2)*A(3,1)
    c21 =-(A(1,2)*A(3,3) - A(1,3)*A(3,2))
    c22 =  A(1,1)*A(3,3) - A(1,3)*A(3,1)
    c23 =-(A(1,1)*A(3,2) - A(1,2)*A(3,1))
    c31 =  A(1,2)*A(2,3) - A(1,3)*A(2,2)
    c32 =-(A(1,1)*A(2,3) - A(1,3)*A(2,1))
    c33 =  A(1,1)*A(2,2) - A(1,2)*A(2,1)
    detA = A(1,1)*c11 + A(1,2)*c12 + A(1,3)*c13
    if (abs(detA) < 1.d-12) stop 'FATAL: singular cell matrix.'
    Ainv = transpose(reshape([c11,c21,c31, c12,c22,c32, c13,c23,c33], [3,3])) / detA
end subroutine invert3x3

! **************************************************
! *** CRITICAL FIX 2: Corrected PBC logic below (transposed) ***
! *** (Note: This function is no longer used by calculate_vanhove_single, 
! *** but is corrected here for future use.)
! **************************************************
subroutine min_image_dr(H, Hinv, dx, dy, dz)
    double precision, intent(in)     :: H(3,3), Hinv(3,3)
    double precision, intent(inout)  :: dx, dy, dz
    double precision :: sx, sy, sz

    ! s = (Hinv)^T * dr (Corrected: Hinv(j,i))
    sx = Hinv(1,1)*dx + Hinv(2,1)*dy + Hinv(3,1)*dz
    sy = Hinv(1,2)*dx + Hinv(2,2)*dy + Hinv(3,2)*dz
    sz = Hinv(1,3)*dx + Hinv(2,3)*dy + Hinv(3,3)*dz

    sx = sx - dnint(sx);  sy = sy - dnint(sy);  sz = sz - dnint(sz)

    ! dr = H^T * s_mic (Corrected: H(j,i))
    dx = H(1,1)*sx + H(2,1)*sy + H(3,1)*sz
    dy = H(1,2)*sx + H(2,2)*sy + H(3,2)*sz
    dz = H(1,3)*sx + H(2,3)*sy + H(3,3)*sz
end subroutine min_image_dr

logical function safe_parse_atom_line(aline, sym, x, y, z)
    character(len=*), intent(in)  :: aline
    character(len=*), intent(out) :: sym
    double precision, intent(out) :: x, y, z
    character(len=1024) :: buf
    integer :: p, ios
    safe_parse_atom_line = .false.
    buf = adjustl(aline)
    p = index(buf, ' ')
    if (p <= 1) return
    sym = adjustl(buf(1:p-1))
    read(buf(p+1:), *, iostat=ios) x, y, z
    if (ios /= 0) return
    safe_parse_atom_line = .true.
end function safe_parse_atom_line

end program ensemble_vanhove_analyzer