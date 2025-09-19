program ensemble_hop_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    ! Input parameters
    character(len=256) :: SIMULATION_NAME
    integer :: NUM_ENSEMBLES
    character(len=8) :: SPECIES
    double precision :: TSTEP_FS
    integer :: DELTA_T, MAX_T_DELTA

    ! Internal
    integer :: i, j, dt, stat
    character(len=512) :: xyz_fname, dir_part, file_part, out_fname
    double precision, allocatable :: all_hop(:,:), mean_hop(:), std_hop(:)

    ! --- Read input file ---
    open(unit=99, file='hop_all.inp', status='old', action='read')
    read(99, *) SIMULATION_NAME, NUM_ENSEMBLES
    read(99, *) SPECIES
    read(99, *) TSTEP_FS, DELTA_T, MAX_T_DELTA
    close(99)

    allocate(all_hop(MAX_T_DELTA, NUM_ENSEMBLES), stat=stat)
    if (stat /= 0) stop 'Allocation error for all_hop'
    all_hop = 0.0d0

    print *, '--- Starting Ensemble Hop Function Analysis ---'
    print *, 'Species: ', trim(SPECIES)

    ! Loop over ensembles
    do i = 1, NUM_ENSEMBLES
        write(dir_part, '(A,A,A,I0)') '../dump/dump_', trim(SIMULATION_NAME), '_', i
        write(file_part, '(I0,A,A)') i, '_product', '.xyz'
        xyz_fname = trim(dir_part) // '/' // trim(file_part)

        print *, 'Processing: ', trim(xyz_fname)

        call calculate_hop_single(xyz_fname, SPECIES, DELTA_T, MAX_T_DELTA, all_hop(:,i))
    end do

    ! Ensemble average
    allocate(mean_hop(MAX_T_DELTA), std_hop(MAX_T_DELTA))
    do dt = 1, MAX_T_DELTA
        mean_hop(dt) = sum(all_hop(dt,:)) / dble(NUM_ENSEMBLES)
        if (NUM_ENSEMBLES > 1) then
            std_hop(dt) = sqrt(sum((all_hop(dt,:) - mean_hop(dt))**2) / dble(NUM_ENSEMBLES-1))
        else
            std_hop(dt) = 0.0d0
        end if
    end do

    ! Write output
    write(out_fname, '(A,A,A)') '../result/hop_all_', trim(SPECIES), '.dat'
    open(unit=50, file=trim(out_fname), status='replace')
    write(50, '(A)') '# Δt(fs)   hop_mean   hop_std'
    do dt = 1, MAX_T_DELTA
        write(50,'(I10, 2F20.8)') dt*DELTA_T*int(TSTEP_FS), mean_hop(dt), std_hop(dt)
    end do
    close(50)

    print *, '--- Hop analysis complete. Output: ', trim(out_fname)

    deallocate(all_hop, mean_hop, std_hop)

contains

! ==========================================================
! Subroutine: Hop calculation for single trajectory
! ==========================================================
subroutine calculate_hop_single(xyz_fname, species, DELTA_T, MAX_T_DELTA, hop_out)
    implicit none
    character(len=*), intent(in) :: xyz_fname, species
    integer, intent(in) :: DELTA_T, MAX_T_DELTA
    double precision, intent(out) :: hop_out(MAX_T_DELTA)

    type frame_data
        integer :: n_atoms
        double precision :: H(3,3), Hinv(3,3)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=10), allocatable :: sym(:)
    end type frame_data

    type(frame_data), allocatable :: trajectory(:)
    integer :: num_frames, ios, i, j, t0, dt, p
    integer :: n_species
    integer, allocatable :: species_idx(:)
    double precision :: hval, sum_h, meanA(3), meanB(3)

    ! === Load trajectory ===
    call load_xyz(xyz_fname, trajectory, num_frames)

    ! Identify species
    n_species = 0
    do i = 1, trajectory(1)%n_atoms
        if (trim(trajectory(1)%sym(i)) == trim(species)) n_species = n_species + 1
    end do
    if (n_species == 0) return
    allocate(species_idx(n_species))
    j = 0
    do i = 1, trajectory(1)%n_atoms
        if (trim(trajectory(1)%sym(i)) == trim(species)) then
            j=j+1; species_idx(j)=i
        end if
    end do

    hop_out = 0.0d0

    ! Main loop
    do dt = 1, MAX_T_DELTA
        sum_h = 0.0d0
        do t0 = DELTA_T, num_frames-DELTA_T
            do i = 1, n_species
                p = species_idx(i)
                call compute_means(trajectory, p, t0, DELTA_T, meanA, meanB)
                call compute_hval(trajectory, p, t0, DELTA_T, meanA, meanB, hval)
                sum_h = sum_h + hval
            end do
        end do
        hop_out(dt) = sum_h / dble(n_species * (num_frames-2*DELTA_T))
    end do

    deallocate(species_idx)
    do i=1,num_frames
        deallocate(trajectory(i)%x, trajectory(i)%y, trajectory(i)%z, trajectory(i)%sym)
    end do
    deallocate(trajectory)
end subroutine calculate_hop_single

! ==========================================================
! Subroutines: Hop function helper
! ==========================================================
subroutine compute_means(trajectory, p, t0, DELTA_T, meanA, meanB)
    type(frame_data), intent(in) :: trajectory(:)
    integer, intent(in) :: p, t0, DELTA_T
    double precision, intent(out) :: meanA(3), meanB(3)
    integer :: i

    meanA = 0.0d0; meanB = 0.0d0
    do i = t0-DELTA_T/2, t0
        meanA(1)=meanA(1)+trajectory(i)%x(p)
        meanA(2)=meanA(2)+trajectory(i)%y(p)
        meanA(3)=meanA(3)+trajectory(i)%z(p)
    end do
    meanA = meanA / dble(DELTA_T/2+1)

    do i = t0, t0+DELTA_T/2
        meanB(1)=meanB(1)+trajectory(i)%x(p)
        meanB(2)=meanB(2)+trajectory(i)%y(p)
        meanB(3)=meanB(3)+trajectory(i)%z(p)
    end do
    meanB = meanB / dble(DELTA_T/2+1)
end subroutine compute_means

subroutine compute_hval(trajectory, p, t0, DELTA_T, meanA, meanB, hval)
    type(frame_data), intent(in) :: trajectory(:)
    integer, intent(in) :: p, t0, DELTA_T
    double precision, intent(in) :: meanA(3), meanB(3)
    double precision, intent(out) :: hval
    integer :: i
    double precision :: dA, dB, sumA, sumB

    sumA=0.0d0; sumB=0.0d0
    do i = t0-DELTA_T/2, t0
        dA = (trajectory(i)%x(p)-meanB(1))**2 + (trajectory(i)%y(p)-meanB(2))**2 + (trajectory(i)%z(p)-meanB(3))**2
        sumA = sumA + dA
    end do
    sumA = sumA / dble(DELTA_T/2+1)

    do i = t0, t0+DELTA_T/2
        dB = (trajectory(i)%x(p)-meanA(1))**2 + (trajectory(i)%y(p)-meanA(2))**2 + (trajectory(i)%z(p)-meanA(3))**2
        sumB = sumB + dB
    end do
    sumB = sumB / dble(DELTA_T/2+1)

    hval = sqrt(sumA * sumB)
end subroutine compute_hval

! ==========================================================
! Subroutines: XYZ Reader
! ==========================================================
subroutine load_xyz(xyz_fname, trajectory, num_frames)
    implicit none
    character(len=*), intent(in) :: xyz_fname
    type(frame_data), allocatable, intent(out) :: trajectory(:)
    integer, intent(out) :: num_frames

    type frame_data
        integer :: n_atoms
        double precision :: H(3,3), Hinv(3,3)
        double precision, allocatable :: x(:), y(:), z(:)
        character(len=10), allocatable :: sym(:)
    end type frame_data

    character(len=1024) :: header_line, aline, lattice_str
    integer :: ios, p1, p2, i, j
    double precision :: detH

    ! Count frames
    open(unit=10, file=trim(xyz_fname), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(error_unit,'(A,A)') 'WARNING: Cannot open ', trim(xyz_fname)
        num_frames=0; return
    end if
    num_frames=0
    do
        read(10,*,iostat=ios) p1
        if (ios /= 0) exit
        read(10,'(A)') header_line
        do i=1,p1
            read(10,'(A)') aline
        end do
        num_frames=num_frames+1
    end do
    rewind(10)

    allocate(trajectory(num_frames))
    do i=1,num_frames
        read(10,*) trajectory(i)%n_atoms
        read(10,'(A)') header_line
        allocate(trajectory(i)%x(trajectory(i)%n_atoms), &
                 trajectory(i)%y(trajectory(i)%n_atoms), &
                 trajectory(i)%z(trajectory(i)%n_atoms), &
                 trajectory(i)%sym(trajectory(i)%n_atoms))
        do j=1,trajectory(i)%n_atoms
            read(10,'(A)') aline
            call parse_atom_line(aline, trajectory(i)%sym(j), &
                                 trajectory(i)%x(j), trajectory(i)%y(j), trajectory(i)%z(j))
        end do
        p1 = index(header_line, 'Lattice="')
        p2 = index(header_line(p1+9:), '"')
        lattice_str = header_line(p1+9 : (p1+8)+p2-1)
        call read_lattice_9(lattice_str, trajectory(i)%H)
        call invert3x3(trajectory(i)%H, trajectory(i)%Hinv, detH)
    end do
    close(10)
end subroutine load_xyz

subroutine parse_atom_line(aline, sym, x, y, z)
    character(len=*), intent(in) :: aline
    character(len=*), intent(out) :: sym
    double precision, intent(out) :: x, y, z
    character(len=1024) :: buf
    integer :: p, ios
    buf=adjustl(aline)
    p=index(buf,' ')
    sym=adjustl(buf(1:p-1))
    read(buf(p+1:),*,iostat=ios) x,y,z
    if (ios /= 0) stop 'Bad atom line.'
end subroutine parse_atom_line

subroutine read_lattice_9(str, H)
    character(len=*), intent(in) :: str
    double precision, intent(out) :: H(3,3)
    double precision :: a(9)
    integer :: ios
    read(str,*,iostat=ios) a
    if (ios /= 0) stop 'Bad lattice numbers.'
    H = reshape(a, [3,3])
end subroutine read_lattice_9

subroutine invert3x3(A, Ainv, detA)
    double precision, intent(in) :: A(3,3)
    double precision, intent(out) :: Ainv(3,3), detA
    double precision :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    c11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
    c12=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    c13=A(2,1)*A(3,2)-A(2,2)*A(3,1)
    c21=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    c22=A(1,1)*A(3,3)-A(1,3)*A(3,1)
    c23=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    c31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
    c32=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    c33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    detA=A(1,1)*c11+A(1,2)*c12+A(1,3)*c13
    if (abs(detA)<1.d-12) stop 'Singular cell matrix.'
    Ainv=transpose(reshape([c11,c21,c31,c12,c22,c32,c13,c23,c33],[3,3]))/detA
end subroutine invert3x3

end program ensemble_hop_analyzer

