PROGRAM ensemble_avg
    IMPLICIT NONE
    
    INTEGER :: NUM_ENSEMBLES, i, j, ios, frame_count, unit_in, unit_out
    CHARACTER(LEN=200) :: SIMULATION_NAME, BASE_FILE_NAME, xyz_fname, OUTPUT_FILE
    CHARACTER(LEN=500) :: dummy_line, header

    DOUBLE PRECISION, ALLOCATABLE :: step(:), time_ps(:)
    DOUBLE PRECISION, ALLOCATABLE :: temp(:,:), pressure_atm(:,:), pe(:,:), ke(:,:), Etot(:,:), vol(:,:), density(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: temp_time_avg(:), pressure_atm_time_avg(:), pe_time_avg(:), ke_time_avg(:), Etot_time_avg(:)
    DOUBLE PRECISION, ALLOCATABLE :: vol_time_avg(:), density_time_avg(:)

    DOUBLE PRECISION :: temp_sum, pressure_atm_sum, pe_sum, ke_sum, Etot_sum, vol_sum, density_sum
    DOUBLE PRECISION :: temp_sq_sum, pressure_atm_sq_sum, pe_sq_sum, ke_sq_sum, Etot_sq_sum, vol_sq_sum, density_sq_sum
    DOUBLE PRECISION :: temp_ensemble_sum, pressure_atm_ensemble_sum, pe_ensemble_sum, ke_ensemble_sum, Etot_ensemble_sum
    DOUBLE PRECISION :: vol_ensemble_sum, density_ensemble_sum
    DOUBLE PRECISION :: temp_ensemble_avg, pressure_atm_avg, pe_avg, ke_avg, Etot_avg, vol_avg, density_avg
    DOUBLE PRECISION :: temp_std, pressure_atm_std, pe_std, ke_std, Etot_std, vol_std, density_std

    ! === 1. Read input settings ===
    OPEN(UNIT=10, FILE='ensemble.inp', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    READ(10, *) SIMULATION_NAME
    READ(10, *) NUM_ENSEMBLES
    READ(10, *) BASE_FILE_NAME
    CLOSE(10)

    PRINT *, 'Simulation Name: ', TRIM(SIMULATION_NAME)
    PRINT *, 'Number of Ensembles: ', NUM_ENSEMBLES
    PRINT *, 'Base Filename: ', TRIM(BASE_FILE_NAME)

    ! === 2. Determine frame count from first file ===
    WRITE(xyz_fname, '(A,"thermo_",A,I0,"/",A,I0,A)') '../thermo/', TRIM(SIMULATION_NAME), 1, &
    TRIM(SIMULATION_NAME), 1, TRIM(BASE_FILE_NAME)

    OPEN(UNIT=11, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) STOP 'Cannot open first file.'
    
    READ(11, '(A)') header
    frame_count = 0
    DO
        READ(11, '(A)', IOSTAT=ios) dummy_line
        IF (ios /= 0) EXIT
        frame_count = frame_count + 1
    END DO
    CLOSE(11)
    PRINT *, 'Frames detected:', frame_count

    ALLOCATE(step(frame_count), time_ps(frame_count))
    ALLOCATE(temp(NUM_ENSEMBLES, frame_count), pressure_atm(NUM_ENSEMBLES, frame_count), &
             pe(NUM_ENSEMBLES, frame_count), ke(NUM_ENSEMBLES, frame_count), Etot(NUM_ENSEMBLES, frame_count), &
             vol(NUM_ENSEMBLES, frame_count), density(NUM_ENSEMBLES, frame_count))
    ALLOCATE(temp_time_avg(NUM_ENSEMBLES), pressure_atm_time_avg(NUM_ENSEMBLES), pe_time_avg(NUM_ENSEMBLES), &
             ke_time_avg(NUM_ENSEMBLES), Etot_time_avg(NUM_ENSEMBLES), vol_time_avg(NUM_ENSEMBLES), density_time_avg(NUM_ENSEMBLES))

    ! === 3. Read each ensemble file ===
    DO i = 1, NUM_ENSEMBLES
        WRITE(xyz_fname, '(A,"thermo_",A,I0,"/",A,I0,A)') '../thermo/', TRIM(SIMULATION_NAME), i, &
        TRIM(SIMULATION_NAME), i, TRIM(BASE_FILE_NAME)
        OPEN(UNIT=11, FILE=TRIM(xyz_fname), STATUS='OLD', ACTION='READ', IOSTAT=ios)
        IF (ios /= 0) THEN
            PRINT *, 'Skipping missing file:', TRIM(xyz_fname)
            CYCLE
        END IF
        READ(11, '(A)') header
        DO j = 1, frame_count
            READ(11, *, IOSTAT=ios) step(j), time_ps(j), temp(i,j), pressure_atm(i,j), pe(i,j), ke(i,j), &
            Etot(i,j), vol(i,j), density(i,j)

            IF (ios /= 0) EXIT
        END DO
        CLOSE(11)
    END DO

    ! === 4. Calculate time averages per ensemble ===
    DO i = 1, NUM_ENSEMBLES
        temp_sum = 0.d0; pressure_atm_sum = 0.d0; pe_sum = 0.d0; ke_sum = 0.d0
        Etot_sum = 0.d0; vol_sum = 0.d0; density_sum = 0.d0
        DO j = 1, frame_count
            temp_sum = temp_sum + temp(i,j)
            pressure_atm_sum = pressure_atm_sum + pressure_atm(i,j)
            pe_sum = pe_sum + pe(i,j)
            ke_sum = ke_sum + ke(i,j)
            Etot_sum = Etot_sum + Etot(i,j)
            vol_sum = vol_sum + vol(i,j)
            density_sum = density_sum + density(i,j)
        END DO
        temp_time_avg(i) = temp_sum / frame_count
        pressure_atm_time_avg(i) = pressure_atm_sum / frame_count
        pe_time_avg(i) = pe_sum / frame_count
        ke_time_avg(i) = ke_sum / frame_count
        Etot_time_avg(i) = Etot_sum / frame_count
        vol_time_avg(i) = vol_sum / frame_count
        density_time_avg(i) = density_sum / frame_count
    END DO

    ! === 5. Ensemble average & standard deviation ===
    temp_ensemble_sum = 0.d0; pressure_atm_ensemble_sum = 0.d0; pe_ensemble_sum = 0.d0
    ke_ensemble_sum = 0.d0; Etot_ensemble_sum = 0.d0; vol_ensemble_sum = 0.d0; density_ensemble_sum = 0.d0
    temp_sq_sum = 0.d0; pressure_atm_sq_sum = 0.d0; pe_sq_sum = 0.d0; ke_sq_sum = 0.d0
    Etot_sq_sum = 0.d0; vol_sq_sum = 0.d0; density_sq_sum = 0.d0

    DO i = 1, NUM_ENSEMBLES
        temp_ensemble_sum = temp_ensemble_sum + temp_time_avg(i)
        pressure_atm_ensemble_sum = pressure_atm_ensemble_sum + pressure_atm_time_avg(i)
        pe_ensemble_sum = pe_ensemble_sum + pe_time_avg(i)
        ke_ensemble_sum = ke_ensemble_sum + ke_time_avg(i)
        Etot_ensemble_sum = Etot_ensemble_sum + Etot_time_avg(i)
        vol_ensemble_sum = vol_ensemble_sum + vol_time_avg(i)
        density_ensemble_sum = density_ensemble_sum + density_time_avg(i)

        temp_sq_sum = temp_sq_sum + temp_time_avg(i)**2
        pressure_atm_sq_sum = pressure_atm_sq_sum + pressure_atm_time_avg(i)**2
        pe_sq_sum = pe_sq_sum + pe_time_avg(i)**2
        ke_sq_sum = ke_sq_sum + ke_time_avg(i)**2
        Etot_sq_sum = Etot_sq_sum + Etot_time_avg(i)**2
        vol_sq_sum = vol_sq_sum + vol_time_avg(i)**2
        density_sq_sum = density_sq_sum + density_time_avg(i)**2
    END DO

    temp_ensemble_avg = temp_ensemble_sum / NUM_ENSEMBLES
    pressure_atm_avg  = pressure_atm_ensemble_sum / NUM_ENSEMBLES
    pe_avg            = pe_ensemble_sum / NUM_ENSEMBLES
    ke_avg            = ke_ensemble_sum / NUM_ENSEMBLES
    Etot_avg          = Etot_ensemble_sum / NUM_ENSEMBLES
    vol_avg           = vol_ensemble_sum / NUM_ENSEMBLES
    density_avg       = density_ensemble_sum / NUM_ENSEMBLES

    temp_std = SQRT(MAX(0.d0, temp_sq_sum / NUM_ENSEMBLES - temp_ensemble_avg**2))
    pressure_atm_std = SQRT(MAX(0.d0, pressure_atm_sq_sum / NUM_ENSEMBLES - pressure_atm_avg**2))
    pe_std = SQRT(MAX(0.d0, pe_sq_sum / NUM_ENSEMBLES - pe_avg**2))
    ke_std = SQRT(MAX(0.d0, ke_sq_sum / NUM_ENSEMBLES - ke_avg**2))
    Etot_std = SQRT(MAX(0.d0, Etot_sq_sum / NUM_ENSEMBLES - Etot_avg**2))
    vol_std = SQRT(MAX(0.d0, vol_sq_sum / NUM_ENSEMBLES - vol_avg**2))
    density_std = SQRT(MAX(0.d0, density_sq_sum / NUM_ENSEMBLES - density_avg**2))

    ! === 6. Write output ===
    OUTPUT_FILE = '../result/Thermo_ensemble.txt'
    OPEN(UNIT=22, FILE=OUTPUT_FILE, STATUS='REPLACE', ACTION='WRITE')
    WRITE(22, '(A)') '# Ensemble average and std results'
    WRITE(22, '(A)') '# i  Temp_avg  Temp_std  Press_avg  Press_std  PE_avg  KE_avg  Etot_avg  Vol_avg  Density_avg'
    DO i = 1, NUM_ENSEMBLES
        WRITE(22, '(I3, 9(1X,F14.6))') i, temp_time_avg(i), temp_std, pressure_atm_time_avg(i), pressure_atm_std, &
            pe_time_avg(i), ke_time_avg(i), Etot_time_avg(i), vol_time_avg(i), density_time_avg(i)
    END DO
    CLOSE(22)
    PRINT *, 'Output saved in: ', TRIM(OUTPUT_FILE)

END PROGRAM ensemble_avg
