// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scVelocity
// Submit example: nextflow run sc_multi_sample.nf -params-file sc_multi_sample.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------DEFINE GLOBAL PARAMETERS----------------*/
// scripts
code_merge_integrate = "$SCRATCH/src/merge_or_integrate.R"
code_dim_reduction = "$SCRATCH/src/dimensionality_reduction.R"
code_cluster = "$SCRATCH/src/cluster.R"
code_de_genes = "$SCRATCH/src/de_genes.R"

/*--------------MAIN WORKFLOW----------------*/

workflow {
    write_params()
    
    inputfiles = Channel
        .of(params.in_seurat_rds)

    merge_integrate(code_merge_integrate,
        inputfiles)

    dim_reduction(code_dim_reduction,
        merge_integrate.out)

    cluster(code_cluster,
        dim_reduction.out.scobj)

    de_genes(code_de_genes,
        cluster.out.scobj_cluster)
}

/*--------------PROCESSES----------------*/

process write_params{ // write input parameters to file
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        path "*.txt"
    """
    echo \$(date) > parameters_sc_multi_sample_${params.id}.txt
    echo "Parameters used to merge/integrate:" >> parameters_sc_multi_sample_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_sc_multi_sample_${params.id}.txt
    """
}

process merge_integrate{ // merge or integrate Seurat datasets
    publishDir params.outfolder+'merge_integrate', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path sc_file
    output:
        file  "*.rds"
    """
    Rscript code.R --task=${params.task} --method=${params.method} --project=${params.project} --nhvg=${params.nhvg} --ncores=${params.ncores} --npcs=${params.npcs} --integrate_seurat_algorithm=${params.integrate_seurat_algorithm} --reference=${params.reference} --integrateBy=${params.integrateBy} --sample.tree=${params.sampleTree} --nclust=${params.nclust} --covars=${params.covars} ${sc_file}
    """
}

process dim_reduction{
    publishDir params.outfolder+'dim_reduc', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit: scobj
        path "*.pdf"
    """
    Rscript code.R --scobject=${scobj} --npcs='${params.npcs}' --task=${params.task} --ncores='${params.ncores}'
    """
}

process cluster{
    publishDir params.outfolder+'cluster', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit:scobj_cluster
        path "*.pdf"
    """
    Rscript code.R --scobject=${scobj} --res='${params.res}' --task=${params.task} --ncores='${params.ncores}'
    """
}

process de_genes{
    publishDir params.outfolder+'cluster', mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        path "code.R"
        path scobj
    output:
        path '*.tsv'
        path '*.pdf' 
    """
    Rscript code.R --sc_file=${scobj} --test.use='${params.test_use}' --min.cells.group='${params.min_cell_group}' --min.pct='${params.min_pct}' --logfc.threshold='${params.logfc_threshold}' --resolution='${params.res_de_genes}' --features='${params.features}' --no.of.pages='${params.no_of_pages}' --ncores='${params.ncores}'
    """
}
