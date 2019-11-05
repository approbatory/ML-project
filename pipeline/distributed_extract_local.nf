#!/usr/bin/env nextflow

params.movie_file = "$HOME/20Hz_data/c14m4/c14m4d15/_data/c14m4d15_gfix_rm_cr_mc_cr_norm40_dff_cechunk.hdf5"
params.dataset_name = "/Data/Images"
params.partition = 4

movie = Channel.from(file(params.movie_file))

process extract_part {
    executor 'local' 
    //queue 'normal'
    memory '4 GB'
    cpus 1
    time '3h'
    //module 'matlab'
    maxForks 2

    input:
    file M from movie
    each part_index from Channel.from(1..(params.partition*params.partition))

    output:
    set part_index, "extract_part_*.mat" into extract_part

    """
    matlab -nodisplay -r "addpath('~/ML-project/pipeline'); launch_extract_part('$M', '$params.dataset_name', 'testname', $part_index, $params.partition); exit"
    """
}
