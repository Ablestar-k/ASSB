program MSD_ensemble
    implicit none

    integer :: i, j, k
    integer :: n_ensemble, n_mol, n_bin
    integer :: ios, unit_in, unit_out
    integer :: frame_count

    double precision :: time_diff
    double precision :: bin_size

    double precision :: atom_x, atom_y, atom_z
    double precision :: fx, fy, fz
    double precision :: mx, my, mz

    ! r(Ensemble index, Current time, coordinate)
    ! MSD(Ensemble index, dt)
    double precision, allocatable :: r(:,:,:)
    double precision, allocatable :: msd(:,:)

    character(len=100) :: header, dummy_line
    character(len=50) :: simulation_name
    character(len=10) :: species_name
    character(len=100) :: atom_index
    character(len=100) :: xyz_fname
    character(len=100) :: input_file, output_file

    open(unit=10, file='MSD.inp', status='old', action='read', iostat=ios)
    read(10,*) simulation_name, time_diff, species_name, n_ensemble, n_mol, n_bin, bin_size 

    unit_in = 20
    write(xyz_fname, '(A,A,I0,A)') '../dump/', simulation_name, 1, '_product_nve.xyz'
    PRINT *, 'Initializing with file: ', TRIM(xyz_fname)

    OPEN(UNIT=unit_in, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error opening first file. IOSTAT =', ios
        STOP 1
    END IF

    ! Count frames to determine array sizes
    frame_count = 0
    DO
        READ(unit_in, '(A)') header ! Read and discard header
        do i = 1, n_mol
            READ(unit_in, '(A)', IOSTAT=ios) dummy_line
            IF (ios /= 0) EXIT
        END DO

        frame_count = frame_count + 1

        IF (ios /= 0) EXIT
    END DO

    rewind(unit_in)
    
    PRINT *, 'Number of frames found: ', frame_count

    ! Allocate arrays
    ALLOCATE(r(n_ensemble, frame_count, 3), &
             msd(n_ensemble, n_bin), stat=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error: Memory allocation failed.'
        STOP 2
    END IF  
    PRINT *, 'Allocated arrays for MSD calculation.'

    ! Read coordinates from the ensemble files
    DO i = 1, n_ensemble
        write(xyz_fname, '(A,A,I0,A)') '../dump/', simulation_name, i, '_product_nve.xyz'
        PRINT *, 'Processing file: ', TRIM(xyz_fname)

        OPEN(UNIT=unit_in, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, 'Error opening file, skipping: ', TRIM(xyz_fname)
            CYCLE
        END IF

        DO j = 1, frame_count
            READ(unit_in, '(A)') header ! Read header
            do k = 1, n_mol
                IF (ios /= 0) THEN
                    PRINT *, 'Error reading atom line at frame', j, 'from file:', TRIM(xyz_fname)
                    EXIT
                END IF
                READ(unit_in, *, IOSTAT=ios) atom_index, atom_x, atom_y, atom_z, 
                
                r(i,j,1), r(i,j,2), r(i,j,3)

                IF (ios /= 0) THEN
                    PRINT *, 'Error reading data at frame', j, 'from file:', TRIM(xyz_fname)
                    EXIT
                END IF
            END DO
        END DO
        CLOSE(unit_in)
    END DO

    PRINT *, 'Finished reading coordinates from all ensemble files.'

    ! Calculate MSD for each ensemble
    DO i = 1, n_ensemble
        PRINT *, 'Calculating MSD for ensemble:', i
        DO j = 1, n_bin
            msd(i,j) = 0.0
        END DO

        DO j = 1, frame_count - 1
            DO k = j + 1, frame_count
                time_diff = (k - j) * time_diff
                ! Calculate MSD for each bin
                DO n_bin = 1, n_bin
                    msd(i,n_bin) = msd(i,n_bin) + &
                        ((r(i,k,1) - r(i,j,1))**2 + (r(i,k,2) - r(i,j,2))**2 + (r(i,k,3) - r(i,j,3))**2)
                END DO
            END DO
        END DO

        ! Normalize by the number of pairs and bin size
        DO n_bin = 1, n_bin
            msd(i,n_bin) = msd(i,n_bin) / ((frame_count - 1) * bin_size)
        END DO
    END DO

    PRINT *, 'Finished calculating MSD for all ensembles.'

    ! Write MSD results to output file
    unit_out = 22
    output_file = '../dump/msd_results.txt'
    OPEN(UNIT=unit_out, FILE=output_file, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error opening output file: ', output_file
        STOP 3
    END IF  
    WRITE(unit_out, '(A)') 'Ensemble MSD Results'
    WRITE(unit_out, '(A,I0,A)') 'Number of Ensembles: ', n_ensemble
    WRITE(unit_out, '(A,I0,A)') 'Number of Frames: ', frame_count
    WRITE(unit_out, '(A,I0,A)') 'Number of Bins: ', n_bin
    WRITE(unit_out, '(A,F10.5,A)') 'Bin Size: ', bin_size, ' Angstroms' 
    WRITE(unit_out, '(A)') 'MSD Results:'
    DO i = 1, n_ensemble
        WRITE(unit_out, '(A,I0,A)') 'Ensemble ', i, ':'
        DO j = 1, n_bin
            WRITE(unit_out, '(I0,F10.5)') j, msd(i,j)
        END DO
    END DO  
    CLOSE(unit_out)

    PRINT *, 'MSD results written to: ', output_file
    PRINT *, 'Program completed successfully.'

END PROGRAM MSD_ensemble
