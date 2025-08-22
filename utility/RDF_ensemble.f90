program ensemble_rdf_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    character(len=20), parameter :: ENSEMBLE_DIR_PATTERN = 'dump_NTOC_{}'
    character(len=20), parameter :: FILENAME_PATTERN     = '{}_product.xyz'
    integer, parameter          :: NUM_ENSEMBLES         = 5
    double precision, parameter :: BINSIZE               = 0.05d0
    double precision, parameter :: MAX_R                 = 10.0d0

    character(len=8) :: SPECIES1, SPECIES2

    integer :: i, j, n_bins, stat
    double precision :: r_center, mean_gr, std_gr, sum_sq
    character(len=512) :: xyz_fname, dir_part, file_part, out_fname

    double precision, allocatable :: gr_values(:,:)     
    double precision, allocatable :: current_gr(:)       
    double precision, allocatable :: mean_gr_arr(:), std_gr_arr(:)

    open(unit=99, file='g_r.inp', status='old', action='read')
    read(99, *) SPECIES1, SPECIES2
    close(99)

    n_bins = ceiling(MAX_R / BINSIZE)

    allocate(gr_values(n_bins, NUM_ENSEMBLES), stat=stat)
    if (stat /= 0) stop 'Allocation error for gr_values'
    allocate(current_gr(n_bins), stat=stat)
    if (stat /= 0) stop 'Allocation error for current_gr'

    gr_values = 0.0d0

    print *, '--- Starting Ensemble RDF Analysis ---'

    ! =============================
    ! Main loop over ensemble members
    ! =============================
    do i = 1, NUM_ENSEMBLES
        write(dir_part, '(A,I0,A)') 'dump_NTOC_', i
        write(file_part, '(I0,A)')  i, '_product.xyz'
        xyz_fname = trim(dir_part)//'/'//trim(file_part)
        print *, 'Processing: ', trim(xyz_fname)

        call calculate_rdf_single(xyz_fname, SPECIES1, SPECIES2, n_bins, MAX_R, current_gr)
        gr_values(:, i) = current_gr(:)
    end do

    print *, ''
    print *, '--- Aggregating results and calculating statistics ---'

    allocate(mean_gr_arr(n_bins), std_gr_arr(n_bins))

    do i = 1, n_bins
        mean_gr = sum(gr_values(i, :)) / dble(NUM_ENSEMBLES)
        sum_sq  = 0.0d0
        do j = 1, NUM_ENSEMBLES
            sum_sq = sum_sq + (gr_values(i, j) - mean_gr)**2
        end do
        if (NUM_ENSEMBLES > 1) then
            std_gr = sqrt(sum_sq / dble(NUM_ENSEMBLES - 1))
        else
            std_gr = 0.0d0
        end if
        mean_gr_arr(i) = mean_gr
        std_gr_arr(i)  = std_gr
    end do

    write(out_fname, '(A,A,A,A,A)') 'gr_fortran_ensemble_average_', trim(SPECIES1), '_', trim(SPECIES2), '.dat'
    open(unit=50, file=trim(out_fname), status='replace')
    write(50, '(A24, A20, A20)') '# r (Angstrom)', 'g(r)_mean', 'g(r)_std'
    do i = 1, n_bins
        r_center = (dble(i) - 0.5d0) * BINSIZE
        write(50, '(F25.8, F20.8, F20.8)') r_center, mean_gr_arr(i), std_gr_arr(i)
    end do
    close(50)

    print *, '--- Analysis complete. Output saved to: ', trim(out_fname), ' ---'

    deallocate(gr_values, current_gr, mean_gr_arr, std_gr_arr)

contains

! ===============================
! Robust single-file RDF calculator
!  - Safe parsing of atom lines (symbol + x y z only)
!  - Supports non-orthogonal cells via 3x3 H matrix and its inverse
!  - Normalization consistent for AA (i<j) and AB pairs
! ===============================
subroutine calculate_rdf_single(xyz_fname, s1, s2, n_bins, max_r, gr_out)
    implicit none
    character(len=*), intent(in)  :: xyz_fname, s1, s2
    integer,          intent(in)  :: n_bins
    double precision, intent(in)  :: max_r
    double precision, intent(out) :: gr_out(n_bins)

    integer :: ios, i, j, k, num_mol, total_frames
    integer :: bin_index, p1, p2
    integer :: stored_n_species1, stored_n_species2
    double precision :: binsize, dist, pi
    double precision :: H(3,3), Hinv(3,3), detH, stored_V
    double precision :: dx, dy, dz
    logical :: same_species, first_frame, frame_bad
    character(len=1024) :: header_line, aline, lattice_str
    double precision :: r, shell_vol, norm_factor

    double precision, allocatable :: x_crd(:), y_crd(:), z_crd(:)
    character(len=10), allocatable :: atom_symbols(:)
    integer, allocatable :: index_s1(:), index_s2(:)
    integer(8), allocatable :: hist_sum(:)

    pi = 4.0d0 * datan(1.0d0)
    total_frames = 0
    first_frame  = .true.
    binsize      = max_r / dble(n_bins)
    stored_V     = 0.0d0
    stored_n_species1 = 0
    stored_n_species2 = 0

    allocate(hist_sum(n_bins), stat=ios)
    if (ios /= 0) stop 'Allocation error for hist_sum in subroutine'
    hist_sum = 0

    open(unit=10, file=trim(xyz_fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit, '(A,1x,A)') 'WARNING: Could not open file, skipping:', trim(xyz_fname)
        gr_out = 0.0d0; deallocate(hist_sum); return
    end if

    frame_loop: do
        read(10, *, iostat=ios) num_mol
        if (ios /= 0) exit frame_loop
        read(10, '(A)', iostat=ios) header_line
        if (ios /= 0) exit frame_loop

        frame_bad = .false.
        allocate(x_crd(num_mol), y_crd(num_mol), z_crd(num_mol), atom_symbols(num_mol), stat=ios)
        if (ios /= 0) then
            write(error_unit,*) 'Allocation error in frame; skipping.'
            frame_bad = .true.
        else
            do i = 1, num_mol
                read(10,'(A)', iostat=ios) aline
                if (ios /= 0) then
                    frame_bad = .true.; exit
                end if
                if (.not. safe_parse_atom_line(aline, atom_symbols(i), x_crd(i), y_crd(i), z_crd(i))) then
                    write(error_unit,'(A)') 'WARNING: bad atom line skipped (dropping frame): '//trim(aline)
                    frame_bad = .true.; exit
                end if
            end do
        end if

        if (frame_bad) then
            if (allocated(x_crd)) deallocate(x_crd, y_crd, z_crd, atom_symbols)
            cycle frame_loop
        end if

        ! --- Parse Lattice 3x3 matrix H from header ---
        p1 = index(header_line, 'Lattice="')
        if (p1 > 0) then
            p2 = index(header_line(p1+9:), '"')
            if (p2 > 0) then
                lattice_str = header_line(p1+9 : (p1+8)+p2-1)
                call read_lattice_9(lattice_str, H)
                call invert3x3(H, Hinv, detH)
                if (detH <= 0.d0) then
                    write(error_unit,*) 'WARNING: non-positive cell volume; dropping frame.'
                    deallocate(x_crd, y_crd, z_crd, atom_symbols)
                    cycle frame_loop
                end if
            else
                write(error_unit,*) 'WARNING: malformed Lattice string; dropping frame.'
                deallocate(x_crd, y_crd, z_crd, atom_symbols)
                cycle frame_loop
            end if
        else
            write(error_unit,*) 'WARNING: Lattice not found; dropping frame.'
            deallocate(x_crd, y_crd, z_crd, atom_symbols)
            cycle frame_loop
        end if

        if (first_frame) then
            stored_V = detH
            do i = 1, num_mol
                if (trim(atom_symbols(i)) == trim(s1)) stored_n_species1 = stored_n_species1 + 1
                if (trim(atom_symbols(i)) == trim(s2)) stored_n_species2 = stored_n_species2 + 1
            end do
            first_frame = .false.
        end if

        if (stored_n_species1 == 0 .or. stored_n_species2 == 0) then
            deallocate(x_crd, y_crd, z_crd, atom_symbols)
            cycle frame_loop
        end if

        allocate(index_s1(stored_n_species1), index_s2(stored_n_species2))
        k = 0
        do i = 1, num_mol
            if (trim(atom_symbols(i)) == trim(s1)) then
                k = k + 1; index_s1(k) = i
            end if
        end do
        k = 0
        do i = 1, num_mol
            if (trim(atom_symbols(i)) == trim(s2)) then
                k = k + 1; index_s2(k) = i
            end if
        end do

        same_species = (trim(s1) == trim(s2))

        if (same_species) then
            do i = 1, stored_n_species1 - 1
                do j = i + 1, stored_n_species1
                    dx = x_crd(index_s1(i)) - x_crd(index_s1(j))
                    dy = y_crd(index_s1(i)) - y_crd(index_s1(j))
                    dz = z_crd(index_s1(i)) - z_crd(index_s1(j))
                    call min_image_dr(H, Hinv, dx, dy, dz)
                    dist = sqrt(dx*dx + dy*dy + dz*dz)
                    if (dist < max_r) then
                        bin_index = ceiling(dist / binsize)
                        if (bin_index > 0 .and. bin_index <= n_bins) hist_sum(bin_index) = hist_sum(bin_index) + 1
                    end if
                end do
            end do
        else
            do i = 1, stored_n_species1
                do j = 1, stored_n_species2
                    dx = x_crd(index_s1(i)) - x_crd(index_s2(j))
                    dy = y_crd(index_s1(i)) - y_crd(index_s2(j))
                    dz = z_crd(index_s1(i)) - z_crd(index_s2(j))
                    call min_image_dr(H, Hinv, dx, dy, dz)
                    dist = sqrt(dx*dx + dy*dy + dz*dz)
                    if (dist < max_r) then
                        bin_index = ceiling(dist / binsize)
                        if (bin_index > 0 .and. bin_index <= n_bins) hist_sum(bin_index) = hist_sum(bin_index) + 1
                    end if
                end do
            end do
        end if

        total_frames = total_frames + 1
        deallocate(x_crd, y_crd, z_crd, atom_symbols, index_s1, index_s2)
    end do frame_loop

    close(10)

    if (total_frames == 0 .or. stored_V <= 0.0d0) then
        gr_out = 0.0d0; deallocate(hist_sum); return
    end if

    same_species = (trim(s1) == trim(s2))
    if (same_species) then
        if (stored_n_species1 <= 1) then
            norm_factor = 0.0d0
        else
            norm_factor = (dble(stored_n_species1) * (dble(stored_n_species1) - 1.0d0) / 2.0d0) / stored_V
        end if
    else
        if (stored_n_species1 == 0 .or. stored_n_species2 == 0) then
            norm_factor = 0.0d0
        else
            norm_factor = (dble(stored_n_species1) * dble(stored_n_species2)) / stored_V
        end if
    end if

    if (norm_factor == 0.0d0) then
        gr_out = 0.0d0; deallocate(hist_sum); return
    end if

    do i = 1, n_bins
        r = (dble(i) - 0.5d0) * binsize
        if (r == 0.0d0) then
            gr_out(i) = 0.0d0
        else
            shell_vol = 4.0d0 * pi * r**2 * binsize
            gr_out(i) = (dble(hist_sum(i)) / dble(total_frames)) / (shell_vol * norm_factor)
        end if
    end do

    deallocate(hist_sum)
end subroutine calculate_rdf_single

! ===============================
! Helpers
! ===============================

pure subroutine read_lattice_9(str, H)
    character(len=*), intent(in)  :: str
    double precision,  intent(out) :: H(3,3)
    double precision :: a(9)
    integer :: ios
    a = 0.d0
    read(str, *, iostat=ios) a
    if (ios /= 0) stop 'FATAL: cannot parse 9 lattice numbers.'
    H(1,1)=a(1); H(1,2)=a(2); H(1,3)=a(3)
    H(2,1)=a(4); H(2,2)=a(5); H(2,3)=a(6)
    H(3,1)=a(7); H(3,2)=a(8); H(3,3)=a(9)
end subroutine read_lattice_9

pure subroutine invert3x3(A, Ainv, detA)
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

pure subroutine min_image_dr(H, Hinv, dx, dy, dz)
    double precision, intent(in)    :: H(3,3), Hinv(3,3)
    double precision, intent(inout) :: dx, dy, dz
    double precision :: sx, sy, sz
    sx = Hinv(1,1)*dx + Hinv(1,2)*dy + Hinv(1,3)*dz
    sy = Hinv(2,1)*dx + Hinv(2,2)*dy + Hinv(2,3)*dz
    sz = Hinv(3,1)*dx + Hinv(3,2)*dy + Hinv(3,3)*dz
    sx = sx - dnint(sx);  sy = sy - dnint(sy);  sz = sz - dnint(sz)
    dx = H(1,1)*sx + H(1,2)*sy + H(1,3)*sz
    dy = H(2,1)*sx + H(2,2)*sy + H(2,3)*sz
    dz = H(3,1)*sx + H(3,2)*sy + H(3,3)*sz
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

end program ensemble_rdf_analyzer
