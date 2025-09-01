program ensemble_msd_analyzer
    use, intrinsic :: iso_fortran_env, only: real64, error_unit
    implicit none

    integer :: NUM_ENSEMBLES
    integer :: TOTAL_TIME_MAX
    integer :: NUM_TOTAL_PARTICLES
    double precision :: TIME_STEP_FS
    character(len=100) :: SIMULATION_NAME
    character(len=5) :: TARGET_SPECIES_NAME

    integer :: e, dt, stat, uout, ios, num_target_particles
    double precision :: time_ps, mean_msd_val, std_msd_val, sum_sq
    
    double precision, allocatable :: disp_sum_total(:), msd_per_ensemble(:,:)
    double precision, allocatable :: mean_msd(:), std_msd(:)
    double precision, allocatable :: current_disp_sum(:)
    
    integer, allocatable :: actual_times(:), total_origins(:)
    integer, allocatable :: target_indices(:)

    character(len=512) :: xyz_fname, dir_part, file_part

    print *, "--- Starting Ensemble MSD Analysis ---"

    open(unit=10, file='ensemble_msd.inp', status='old', action='read', iostat=ios)
    read(10,*) SIMULATION_NAME, TARGET_SPECIES_NAME, TIME_STEP_FS, TOTAL_TIME_MAX, NUM_ENSEMBLES, NUM_TOTAL_PARTICLES
    close(10)

    write(dir_part, '(A,I0,A)') 'dump_', SIMULATION_NAME, '_', 1
    write(file_part, '(I0,A)') 1, '_product.xyz'
    xyz_fname = '../dump/' // trim(dir_part) // '/' // trim(file_part)
    call prescan_trajectory(xyz_fname, NUM_TOTAL_PARTICLES, TARGET_SPECIES_NAME, &
                              num_target_particles, target_indices)
    print *, "Found", num_target_particles, " target particles ('"//trim(TARGET_SPECIES_NAME)//"') from pre-scan."


    allocate(disp_sum_total(TOTAL_TIME_MAX), stat=stat); if (stat /= 0) stop 'alloc disp_sum_total failed'
    allocate(msd_per_ensemble(TOTAL_TIME_MAX, NUM_ENSEMBLES), stat=stat); if (stat /= 0) stop 'alloc msd_per_ensemble failed'
    allocate(actual_times(NUM_ENSEMBLES), stat=stat); if (stat /= 0) stop 'alloc actual_times failed'
    
    disp_sum_total = 0.0d0
    msd_per_ensemble = 0.0d0

    ! Loop over all ensembles
    do e = 1, NUM_ENSEMBLES
        allocate(current_disp_sum(TOTAL_TIME_MAX), stat=stat)
        
        write(dir_part, '(A,I0,A)') 'dump_', SIMULATION_NAME, '_', e
        write(file_part, '(I0,A)') e, '_product.xyz'
        xyz_fname = 'dump/' // trim(dir_part) // '/' // trim(file_part)
        
        print *, "Processing Ensemble #", e, ": ", trim(xyz_fname)

        call calculate_disp_sum_single(xyz_fname, TOTAL_TIME_MAX, NUM_TOTAL_PARTICLES, &
                                       target_indices, num_target_particles, &
                                       current_disp_sum, actual_times(e))
        
        disp_sum_total(1:actual_times(e)) = disp_sum_total(1:actual_times(e)) + current_disp_sum(1:actual_times(e))
        
        do dt = 1, actual_times(e) - 1
            msd_per_ensemble(dt, e) = current_disp_sum(dt) / (dble(num_target_particles) * dble(actual_times(e) - dt))
        end do
        deallocate(current_disp_sum)
    end do

    print *, ""
    print *, "--- Aggregating results and calculating statistics ---"

    allocate(mean_msd(TOTAL_TIME_MAX), std_msd(TOTAL_TIME_MAX), total_origins(TOTAL_TIME_MAX))
    total_origins = 0

    ! Calculate total number of origins for each dt
    do dt = 1, TOTAL_TIME_MAX - 1
        do e = 1, NUM_ENSEMBLES
            if (actual_times(e) > dt) then
                total_origins(dt) = total_origins(dt) + (actual_times(e) - dt)
            end if
        end do
    end do

    ! Calculate mean MSD and standard deviation
    do dt = 1, TOTAL_TIME_MAX - 1
        if (total_origins(dt) > 0 .and. num_target_particles > 0) then
            mean_msd(dt) = disp_sum_total(dt) / (dble(num_target_particles) * dble(total_origins(dt)))
        else
            mean_msd(dt) = 0.0d0
        end if
        
        sum_sq = 0.0d0
        do e = 1, NUM_ENSEMBLES
            sum_sq = sum_sq + (msd_per_ensemble(dt, e) - mean_msd(dt))**2
        end do
        
        if (NUM_ENSEMBLES > 1) then
            std_msd(dt) = sqrt(sum_sq / dble(NUM_ENSEMBLES - 1))
        else
            std_msd(dt) = 0.0d0
        end if
    end do

    ! Write final output file
    open(newunit=uout, file='MSD_ensemble_average.dat', status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'Error creating output file.'

    write(uout, '(A18, A24, A24)') '# Time (ps)', 'MSD_mean (Angstrom^2)', 'MSD_std (Angstrom^2)'
    do dt = 1, TOTAL_TIME_MAX - 1
        time_ps = dble(dt) * TIME_STEP_FS / 1000.0d0
        write(uout, '(F18.6, F24.8, F24.8)') time_ps, mean_msd(dt), std_msd(dt)
    end do
    close(uout)

    print *, "--- Analysis complete. Output saved to: MSD_ensemble_average.dat ---"

    deallocate(disp_sum_total, msd_per_ensemble, actual_times, mean_msd, std_msd, total_origins, target_indices)

contains

    subroutine prescan_trajectory(trajectory_file_name, num_total_particles, target_species_name, &
                                  num_target_particles, target_indices)
        implicit none
        character(len=*), intent(in) :: trajectory_file_name, target_species_name
        integer, intent(in) :: num_total_particles
        integer, intent(out) :: num_target_particles
        integer, allocatable, intent(out) :: target_indices(:)

        integer :: i, ios, uin, natoms_in_frame
        character(len=10) :: species_name_from_file
        character(len=2048) :: header_line, line_buffer

        open(newunit=uin, file=trim(trajectory_file_name), status='old', action='read', iostat=ios)
        if (ios /= 0) stop 'Pre-scan: Could not open file.'

        read(uin,*,iostat=ios) natoms_in_frame; if (ios /= 0) stop 'pre-scan fail: atom count'
        read(uin,'(A)',iostat=ios) header_line; if (ios /= 0) stop 'pre-scan fail: header'
        
        num_target_particles = 0
        do i = 1, num_total_particles
            read(uin, '(A)') line_buffer
            read(line_buffer, *) species_name_from_file
            if (trim(species_name_from_file) == trim(target_species_name)) then
                num_target_particles = num_target_particles + 1
            end if
        end do
        if (num_target_particles == 0) stop 'Error: Target species not found.'
        
        allocate(target_indices(num_target_particles), stat=ios); if (ios/=0) stop 'Alloc target_indices failed'
        
        rewind(uin)
        read(uin,*,iostat=ios) natoms_in_frame
        read(uin,'(A)',iostat=ios) header_line
        
        num_target_particles = 0
        do i = 1, num_total_particles
            read(uin, '(A)') line_buffer
            read(line_buffer, *) species_name_from_file
            if (trim(species_name_from_file) == trim(target_species_name)) then
                num_target_particles = num_target_particles + 1
                target_indices(num_target_particles) = i
            end if
        end do
        close(uin)
    end subroutine prescan_trajectory

    subroutine calculate_disp_sum_single(trajectory_file_name, total_time_max, num_total_particles, &
                                         target_indices, num_target_particles, disp_sum_out, t_actual)
        implicit none
        character(len=*), intent(in) :: trajectory_file_name
        integer, intent(in) :: total_time_max, num_total_particles, num_target_particles
        integer, intent(in) :: target_indices(num_target_particles)
        double precision, intent(out) :: disp_sum_out(total_time_max)
        integer, intent(out) :: t_actual

        integer :: t, i, dt, ios, uin, natoms_in_frame
        double precision, allocatable :: r(:,:,:), r_unwrapped(:,:,:), cell(:,:,:)
        double precision :: inv_lattice_matrix(3,3), displacement(3), frac_displacement(3)
        double precision :: ax, ay, az, bx, by, bz, cx, cy, cz
        character(len=10) :: species_name_from_file
        character(len=2048) :: header_line, line_buffer
        character(len=512) :: number_string
        integer :: istart, iend

        allocate(r(num_total_particles, total_time_max, 3), stat=ios); if (ios/=0) stop 'Alloc r failed'
        allocate(r_unwrapped(num_total_particles, total_time_max, 3), stat=ios); if (ios/=0) stop 'Alloc r_unwrapped failed'
        allocate(cell(3,3,total_time_max), stat=ios); if (ios/=0) stop 'Alloc cell failed'
        disp_sum_out = 0.0d0

        open(newunit=uin, file=trim(trajectory_file_name), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(error_unit, *) "WARNING: Could not open file, skipping: ", trim(trajectory_file_name)
            t_actual = 0
            return
        end if

        t_actual = 0
        do t = 1, total_time_max
            read(uin,*,iostat=ios) natoms_in_frame
            if (ios /= 0) exit
            read(uin,'(A)',iostat=ios) header_line
            if (ios /= 0) exit
            istart = index(header_line, 'Lattice="') + 9
            iend = index(header_line(istart:), '"')
            number_string = header_line(istart : istart + iend - 2)
            read(number_string, *, iostat=ios) ax,ay,az,bx,by,bz,cx,cy,cz
            if (ios /= 0) stop 'Failed to parse lattice matrix'
            cell(1,1,t)=ax; cell(2,1,t)=ay; cell(3,1,t)=az
            cell(1,2,t)=bx; cell(2,2,t)=by; cell(3,2,t)=bz
            cell(1,3,t)=cx; cell(2,3,t)=cy; cell(3,3,t)=cz
            do i = 1, num_total_particles
                read(uin, '(A)') line_buffer
                read(line_buffer, *, iostat=ios) species_name_from_file, r(i, t, 1), r(i, t, 2), r(i, t, 3)
                if (ios /= 0) stop 'Error reading coordinates.'
            end do
            t_actual = t
        end do
        close(uin)

        if (t_actual == 0) return
        print *, "    -> Read", t_actual, "frames from trajectory."
        
        r_unwrapped(:, 1, :) = r(:, 1, :)
        do t = 2, t_actual
            call invert_matrix_3x3(cell(:,:,t), inv_lattice_matrix)
            do i = 1, num_total_particles
                displacement(:) = r(i, t, :) - r(i, t-1, :)
                frac_displacement(:) = matmul(inv_lattice_matrix, displacement)
                frac_displacement(:) = frac_displacement(:) - anint(frac_displacement(:))
                displacement(:) = matmul(cell(:,:,t), frac_displacement)
                r_unwrapped(i, t, :) = r_unwrapped(i, t-1, :) + displacement(:)
            end do
        end do

        do dt = 1, t_actual - 1
            do t = 1, t_actual - dt
                disp_sum_out(dt) = disp_sum_out(dt) + &
                    sum( ( r_unwrapped(target_indices, t+dt, :) - r_unwrapped(target_indices, t, :) )**2 )
            end do
        end do

        deallocate(r, r_unwrapped, cell)
    end subroutine calculate_disp_sum_single

    subroutine invert_matrix_3x3(A, A_inv)
        double precision, intent(in)  :: A(3,3)
        double precision, intent(out) :: A_inv(3,3)
        double precision :: det_A
        det_A = A(1,1) * (A(2,2)*A(3,3) - A(3,2)*A(2,3)) - &
                A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
                A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        if (abs(det_A) < 1.0d-12) stop 'Error: Matrix is singular.'
        A_inv(1,1) = (A(2,2)*A(3,3) - A(3,2)*A(2,3)) / det_A
        A_inv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / det_A
        A_inv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / det_A
        A_inv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / det_A
        A_inv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / det_A
        A_inv(2,3) = (A(2,1)*A(1,3) - A(1,1)*A(2,3)) / det_A
        A_inv(3,1) = (A(2,1)*A(3,2) - A(3,1)*A(2,2)) / det_A
        A_inv(3,2) = (A(3,1)*A(1,2) - A(1,1)*A(3,2)) / det_A
        A_inv(3,3) = (A(1,1)*A(2,2) - A(2,1)*A(1,2)) / det_A
    end subroutine invert_matrix_3x3

end program ensemble_msd_analyzer
