program MSD_ensemble
    implicit none

    integer, parameter :: NUM_ENSEMBLES = 5

    integer :: i, j, k
    integer :: n_ensemble, n_bin
    integer :: ios, unit_in, unit_out
    integer :: frame_count

    double precision :: time_diff,
    double precision :: bin_size, 

    character(len=10) :: species_name
    character(len=100) :: xyz_fname
    character(len=100) :: input_file, output_file

    unit_in = 11
    write(xyz_fname, '(A,A,I0,A)') '../dump/', 'NTOC_ver1_', 1, '_product_nve.thermo'
    PRINT *, 'Initializing with file: ', TRIM(xyz_fname)

    OPEN(UNIT=unit_in, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error opening first file. IOSTAT =', ios
        STOP 1
    END IF

    ! Count frames to determine array sizes
    frame_count = 0
    READ(unit_in, '(A)') header ! Read and discard header
    DO
        READ(unit_in, '(A)', IOSTAT=ios) dummy_line
        IF (ios /= 0) EXIT
        frame_count = frame_count + 1
    END DO
    PRINT *, 'Number of frames found: ', frame_count