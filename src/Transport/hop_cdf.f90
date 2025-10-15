PROGRAM HOP_CDF
    implicit none

   ! Input parameters
    character(len=256) :: SIMULATION_NAME
    integer :: NUM_ENSEMBLES
    character(len=8) :: SPECIES
    double precision :: BINSIZE
    integer :: H_DELTA, MAX_H_DELTA

    character(len=512) :: data_fname, dir_part, file_part, out_fname, BASE_FILE_NAME

    open(unit=99, file='hop_pdf.inp', status='old', action='read')
    read(99, '(A)') SIMULATION_NAME
    read(99, *) NUM_ENSEMBLES
    read(99, *) NUM_SPECIES
    read(99, *) BINSIZE
    read(99, *) H_DELTA
    read(99, *) MAX_H_DELTA
    read(99, '(A)') BASE_FILE_NAME
    close(99)

    n_bins = int(MAX_H_DELTA / BINSIZE) + 1

    allocate(hop_hist(n_bins), stat = stat)
    hop_hist = 0
    if (stat /= 0) then
        print *, 'Error allocating memory for hop_hist'
        stop
    end if  

    print *, 'Simulation Name: ', SIMULATION_NAME
    print *, 'Number of Ensembles: ', NUM_ENSEMBLES
    print *, 'Species: ', SPECIES
    print *, 'Bin Size: ', BINSIZE
    print *, 'H Delta: ', H_DELTA
    print *, 'Max H Delta: ', MAX_H_DELTA
    

    do i = 1, NUM_SPECIES
        write(data_fname, '(A,I0,A)') '../result/hop/', i , '_hop.dat'
        print *, 'Processing: ', trim(data_fname)

        call calculate_hop_cdf(data_fname, H_DELTA, MAX_H_DELTA, BINSIZE, hop_hist, n_bins)
        all__values(:, :, i) = current_gs(:, :)
    end do


SUBROUTINE calculate_hop_cdf(data_fname, H_DELTA, MAX_H_DELTA, BINSIZE, hop_hist, n_bins)
    implicit none
    character(len=*), intent(in) :: data_fname
    integer, intent(in) :: H_DELTA, MAX_H_DELTA, n_bins
    double precision, intent(in) :: BINSIZE
    integer, dimension(n_bins), intent(inout) :: hop_hist

    integer :: iunit, ios, num_lines
    double precision :: hop_value
    integer :: bin_index

    ! Open the data file
    open(unit=10, file=data_fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(data_fname)
        return
    end if

    ! Read through the file and populate the histogram
    num_lines = 0
    do
        read(10, *, iostat=ios) hop_value
        if (ios /= 0) exit
        num_lines = num_lines + 1

        if (hop_value >= 0.0d0 .and. hop_value <= dble(MAX_H_DELTA)) then
            bin_index = int(hop_value / BINSIZE) + 1
            if (bin_index >= 1 .and. bin_index <= n_bins) then
                hop_hist(bin_index) = hop_hist(bin_index) + 1
            end if
        end if
    end do

    close(10)

    print *, 'Processed ', num_lines, ' lines from ', trim(data_fname)
END SUBROUTINE calculate_hop_cdf

END PROGRAM HOP_CDF