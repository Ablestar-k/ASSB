PROGRAM ensemble_avg
    IMPLICIT NONE

    INTEGER, PARAMETER :: NUM_ENSEMBLES = 5
    INTEGER            :: i, j, ios, frame_count, unit_in, unit_out
    CHARACTER(LEN=200) :: xyz_fname, OUTPUT_FILE
    CHARACTER(LEN=500) :: dummy_line, header

    ! Allocate arrays for ensemble data
    ! (Ensemble index, Data)
    DOUBLE PRECISION, ALLOCATABLE :: step(:), time_ps(:)
    DOUBLE PRECISION, ALLOCATABLE :: temp(:,:), pressure_atm(:,:), pe(:,:), ke(:,:), Etot(:,:), vol(:,:), density(:,:)

    ! Variables for calculation
    DOUBLE PRECISION :: temp_sum, pressure_atm_sum, pe_sum, ke_sum, Etot_sum, vol_sum, density_sum
    DOUBLE PRECISION :: temp_sq_sum, pressure_atm_sq_sum, pe_sq_sum, ke_sq_sum, Etot_sq_sum, vol_sq_sum, density_sq_sum
    DOUBLE PRECISION :: temp_avg, pressure_atm_avg, pe_avg, ke_avg, Etot_avg, vol_avg, density_avg
    DOUBLE PRECISION :: temp_std, pressure_atm_std, pe_std, ke_std, Etot_std, vol_std, density_std

    ! == 1. Initialization and Memory Allocation based on the first file ==
    unit_in = 11
    write(xyz_fname, '(A,A,I0,A)') '../thermo/', 'NTOC_ver1_', 1, '_product_nve.thermo'
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

    ! Allocate arrays
    ALLOCATE(step(frame_count), time_ps(frame_count))
    ALLOCATE(temp(NUM_ENSEMBLES, frame_count), &
             pressure_atm(NUM_ENSEMBLES, frame_count), &
             pe(NUM_ENSEMBLES, frame_count), &
             ke(NUM_ENSEMBLES, frame_count), &
             Etot(NUM_ENSEMBLES, frame_count), &
             vol(NUM_ENSEMBLES, frame_count), &
             density(NUM_ENSEMBLES, frame_count), stat=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error: Memory allocation failed.'
        STOP 2
    END IF

    ! Read step and time from the first file
    REWIND(unit_in)
    READ(unit_in, '(A)') header ! Re-read header
    DO j = 1, frame_count
        READ(unit_in, *) step(j), time_ps(j)
    END DO
    CLOSE(unit_in)

    ! == 2. Read data from all ensemble files ==
    DO i = 1, NUM_ENSEMBLES
        write(xyz_fname, '(A,A,I0,A)') '../thermo/', 'NTOC_ver1_', i, '_product_nve.thermo'
        PRINT *, 'Processing file: ', TRIM(xyz_fname)

        OPEN(UNIT=unit_in, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, 'Error opening file, skipping: ', TRIM(xyz_fname)
            CYCLE
        END IF

        READ(unit_in, '(A)') header ! Read header

        DO j = 1, frame_count
            READ(unit_in, *, IOSTAT=ios) dummy_line, dummy_line, temp(i,j), pressure_atm(i,j), pe(i,j), ke(i,j), &
            Etot(i,j), vol(i,j), density(i,j)
            IF (ios /= 0) THEN
                PRINT *, 'Error reading data at frame', j, 'from file:', TRIM(xyz_fname)
                EXIT
            END IF
        END DO
        CLOSE(unit_in)
    END DO

    ! == 3. Calculate averages and standard deviations and write to file ==
    unit_out = 22
    OUTPUT_FILE = '../thermo/ensemble_avg.txt'
    OPEN(UNIT=unit_out, FILE=OUTPUT_FILE, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
        PRINT *, 'Error opening output file: ', OUTPUT_FILE
        STOP 3
    END IF

    WRITE(unit_out, '(A)') '# Ensemble Average Results:'
    WRITE(unit_out, '(A,A,A,A,A,A,A,A)') '#      Step      Time(ps) ', &
        '     <Temp>      Temp_std ', '        <P>         P_std ', '       <PE>        PE_std ', &
        '       <KE>        KE_std ', '      <Etot>      Etot_std ', '      <Vol>       Vol_std ', &
        '    <Density>   Density_std'

    DO j = 1, frame_count
        ! Initialize sum variables for each frame
        temp_sum = 0.d0; pressure_atm_sum = 0.d0; pe_sum = 0.d0; ke_sum = 0.d0
        Etot_sum = 0.d0; vol_sum = 0.d0; density_sum = 0.d0
        temp_sq_sum = 0.d0; pressure_atm_sq_sum = 0.d0; pe_sq_sum = 0.d0; ke_sq_sum = 0.d0
        Etot_sq_sum = 0.d0; vol_sq_sum = 0.d0; density_sq_sum = 0.d0

        DO i = 1, NUM_ENSEMBLES
            temp_sum = temp_sum + temp(i,j)
            pressure_atm_sum = pressure_atm_sum + pressure_atm(i,j)
            pe_sum = pe_sum + pe(i,j)
            ke_sum = ke_sum + ke(i,j)
            Etot_sum = Etot_sum + Etot(i,j)
            vol_sum = vol_sum + vol(i,j)
            density_sum = density_sum + density(i,j)

            temp_sq_sum = temp_sq_sum + temp(i,j)**2
            pressure_atm_sq_sum = pressure_atm_sq_sum + pressure_atm(i,j)**2
            pe_sq_sum = pe_sq_sum + pe(i,j)**2
            ke_sq_sum = ke_sq_sum + ke(i,j)**2
            Etot_sq_sum = Etot_sq_sum + Etot(i,j)**2
            vol_sq_sum = vol_sq_sum + vol(i,j)**2
            density_sq_sum = density_sq_sum + density(i,j)**2
        END DO

        ! Calculate average
        temp_avg = temp_sum / NUM_ENSEMBLES
        pressure_atm_avg = pressure_atm_sum / NUM_ENSEMBLES
        pe_avg = pe_sum / NUM_ENSEMBLES
        ke_avg = ke_sum / NUM_ENSEMBLES
        Etot_avg = Etot_sum / NUM_ENSEMBLES
        vol_avg = vol_sum / NUM_ENSEMBLES
        density_avg = density_sum / NUM_ENSEMBLES

        ! Calculate standard deviation: sqrt(<x^2> - <x>^2)
        temp_std = SQRT(MAX(0.d0, temp_sq_sum / NUM_ENSEMBLES - temp_avg**2))
        pressure_atm_std = SQRT(MAX(0.d0, pressure_atm_sq_sum / NUM_ENSEMBLES - pressure_atm_avg**2))
        pe_std = SQRT(MAX(0.d0, pe_sq_sum / NUM_ENSEMBLES - pe_avg**2))
        ke_std = SQRT(MAX(0.d0, ke_sq_sum / NUM_ENSEMBLES - ke_avg**2))
        Etot_std = SQRT(MAX(0.d0, Etot_sq_sum / NUM_ENSEMBLES - Etot_avg**2))
        vol_std = SQRT(MAX(0.d0, vol_sq_sum / NUM_ENSEMBLES - vol_avg**2))
        density_std = SQRT(MAX(0.d0, density_sq_sum / NUM_ENSEMBLES - density_avg**2))

        ! Write results for the current frame
        WRITE(unit_out, '(I11, F13.4, 14(1X, E14.6))') &
            step(j), time_ps(j), temp_avg, temp_std, pressure_atm_avg, pressure_atm_std, &
            pe_avg, pe_std, ke_avg, ke_std, Etot_avg, Etot_std, &
            vol_avg, vol_std, density_avg, density_std
    END DO

    CLOSE(unit_out)
    PRINT *, 'Ensemble average calculations completed successfully.'
    PRINT *, 'Results are saved in: ', TRIM(OUTPUT_FILE)

END PROGRAM ensemble_avg