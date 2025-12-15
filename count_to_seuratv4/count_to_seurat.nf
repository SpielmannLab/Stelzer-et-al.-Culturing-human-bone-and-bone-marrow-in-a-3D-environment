// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scVelocity
// Submit example: nextflow run count_to_seurat.nf -params-file count_to_seurat_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------DEFINE GLOBAL PARAMETERS----------------*/

// scripts
code_count_to_seurat="$SCRATCH/read_10x.R"
code_add_velocyto_to_seurat="$SCRATCH/add_velocyto_to_seurat.R"

/*--------------MAIN WORKFLOW----------------*/

workflow {
    // write the parameters to file
    write_params()
    
    // create a channel for parallelising the processing of all samples
    samples_channel=Channel.fromList(params.samples)
        .view()

    // convert 10x count matrix to seurat object
    read_10x(code_count_to_seurat,
        samples_channel.map{it[0..5]})
    
    // based on user input, do velocyto or not
    if (params.perform_velocyto) {
        
        // run velocyto to get spliced/unspliced reads (as a loom file) from bam file
        run_velocyto(samples_channel, 
            params.file_repeatMasker, 
            params.file_genesGTF)
        
        // create a channel with seurat object and loom file
        seurat_velocyto_channel=read_10x.out
            .join(run_velocyto.out)
            .view()

        // and spliced/unspliced "assays" into seurat object
        add_velocyto_to_seurat(code_add_velocyto_to_seurat,
            seurat_velocyto_channel)
    }
}

/*--------------PROCESSES----------------*/

process write_params{ // write input parameters to file
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        file "*.txt"
    """
    echo \$(date) > parameters_${params.id}.txt
    echo "Parameters used to create the Seurat object from the count matrix:" >> parameters_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_${params.id}.txt
    """
}

process read_10x{ // create seurat object from count matrix output by cell ranger pipeline
    publishDir params.outfolder+'1_seurat_10x', mode: 'copy', overwrite: true
    input:
        path "code.R"
        tuple val(condition), val(treatment), val(replicate), val(batch), val(nickname), path(path_count_matrix)
    output:
        tuple val("${condition}-${treatment}-${replicate}-${batch}"), file("*.rds")
    """
    Rscript code.R --infolder=${path_count_matrix} --condition=${condition} --treatment=${treatment} --replicate=${replicate} --batch=${batch} --nickname=${nickname}
    """
}

process run_velocyto{ // run velocyto to create spliced/unspliced counts using bam files from cellranger pipeline
    stageInMode 'copy'
    publishDir params.outfolder+'2_velocyto', mode: 'copy', overwrite: true
    input:
        tuple val(condition), val(treatment), val(replicate), val(batch), val(nickname), path(path_count_matrix), path(path_bam)
        path file_repeatMask
        path file_genesGTF
    output:
        tuple val("${condition}-${treatment}-${replicate}-${batch}"), file("*.loom")
    """
    sampleid="${condition}-${treatment}-${replicate}-${batch}"
    velocyto run -m ${file_repeatMask} -o "./" -b "${path_count_matrix}/barcodes.tsv.gz" -e \$sampleid $path_bam $file_genesGTF
    """
}

process add_velocyto_to_seurat{ // add spliced/unspliced counts into the Seurat object
    publishDir params.outfolder+'3_seurat_10x_n_velocyto', mode: 'copy', overwrite: true
    input:
        path "code.R"
        tuple val(project), path(file_sc_obj), path(file_velocyto)
    output:
        file "*_velocyto.rds"
    """
    Rscript code.R --file_sc_obj=$file_sc_obj --file_velocyto=$file_velocyto
    """
}