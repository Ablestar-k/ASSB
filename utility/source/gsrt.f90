PROGRAM gsrt
    
    implicit none

    open(unit=99, file='gsrt.inp', status='old', action='read')
    read(99, *) SIMULATION_NAME, NUM_ENSEMBLES, BINSIZE, MAX_R 
    read(99, *) SPECIES
    close(99)

    n_bins = ceiling(MAX_R / BINSIZE)
    allocate(gsrt_values(n_bins, NUM_ENSEMBLES), stat=stat)
    if (stat /= 0) stop 'Allocation error for gsrt_values'
    allocate(current_gsrt(n_bins), stat=stat)
    if (stat /= 0) stop 'Allocation error for current_gsrt'

    