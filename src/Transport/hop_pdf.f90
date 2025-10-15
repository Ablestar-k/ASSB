PROGRAM HOP_PDF
    implicit none

   ! Input parameters
    character(len=256) :: SIMULATION_NAME
    integer :: NUM_ENSEMBLES
    character(len=8) :: SPECIES
    double precision :: BINSIZE, MAX_R
    integer :: T_DELTA, MAX_T_DELTA 

    character(len=512) :: xyz_fname, dir_part, file_part, out_fname

    open(unit=99, file='hop_pdf.inp', status='old', action='read')