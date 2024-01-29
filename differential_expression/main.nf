if (!params.processes){
params.processes = 2
}

if (params.help) {

log.info"""
==============================================
Challenge 1 - Quantify the size and variation of knock out effect
==============================================
Usage:
Run the pipeline with default parameters:
nextflow run main.nf -profile docker

Run with user parameters:

nextflow run main.nf -profile docker --input {input.dir}

For more informations on the parameters, please see the README.md file located in the dockers repository 'example_data' folder.
Mandatory arguments:
--input                 Input directory where the RDS files are
--output 				Output directory where the analysis will be stored (Defaults to {input}/Analysis)

Optional arguments:
-- processes            Number of processes to run in parallel. If the execution gives an error exit of 137, that means
                        you ran out of memory. Try reducing the number of processes in parallel. Default 2
Flags:
--help                  Display this message
"""

exit 1
} else {

log.info """\
==============================================
CHALLENGE 1 - Quantify the size and variation of knock out effect
==============================================
input directory: ${params.input}
output directory: ${params.output}
Number of parallel processes: ${params.processes}
"""

}

""" Find within the input directory all the RDS files available """
input_resolved = file(params.input).resolve("*.rds")
input_files = Channel.fromPath(input_resolved)

output_directory = file(params.output)


process rename_sgrna {

    maxForks 2

    input:
    path input_file

    output:
    path "${task.workDir}/tmp/${input_file}"

    """
    Rscript /app/0_rename_sgrna.R -f $input_file
    """
}

process select_top_genotypes {

    maxForks ${params.processes}

    input:
    path input_file
    path output_directory

    output:
    path "${output_directory}/*.rds"

    """
    Rscript /app/1_select_top_genotypes.R -f $input_file -o $output_directory
    """
}


process cell_type_proportions_D18 {


    input:
    path input_file
    path output_directory

    output:
    path "${output_directory}/D18.ctproption.pdf"

    when:
    input_file.name =~ /^Sample_L.*/

    script:
    """
    Rscript /app/1b_d18_cellt_proportions.R -f $input_file -o $output_directory
    """
}


workflow {
    tmp_files = rename_sgrna(input_files)
    tmp_files.view()
    """
    counts = select_top_genotypes(input_files, output_directory)
    ct_d18_pdf = cell_type_proportions_D18(input_files, output_directory)
    """
}