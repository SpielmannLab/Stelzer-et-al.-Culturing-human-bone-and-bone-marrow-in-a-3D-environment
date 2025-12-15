// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scVelocity
// Submit example: nextflow run sc_p_sample.nf -params-file sc_p_sample.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------DEFINE GLOBAL PARAMETERS----------------*/

// scripts
code_qc = "$SCRATCH/src/qc_sc.R"
code_filter_qc = "$SCRATCH/src/filter_qc.R"
code_run_scrublet = "$SCRATCH/src/run_scrublet.py"
code_filter_doublets = "$SCRATCH/src/filter_doublets.R"
code_normalize = "$SCRATCH/src/normalize.R"
code_dim_reduction = "$SCRATCH/src/dimensionality_reduction.R"
code_cluster = "$SCRATCH/src/cluster.R"
code_de_genes = "$SCRATCH/src/de_genes.R"

/*--------------MAIN WORKFLOW----------------*/

workflow {
    write_params()
    
    scobj_channel=Channel.fromList(params.in_seurat_rds)

    qc_unfiltered(code_qc,
        scobj_channel)

    filter_qc(code_filter_qc,
        qc_unfiltered.out.scobj)

    scrublet_filter_doublet(code_run_scrublet,
        code_filter_doublets,
        filter_qc.out.sc_mtx_obj)

    qc_filtered(code_qc,
        scrublet_filter_doublet.out.scobj)

    normalize(code_normalize, 
        qc_filtered.out.scobj)

    dim_reduction(code_dim_reduction,
        normalize.out.scobj)

    cluster(code_cluster,
        dim_reduction.out.scobj)

    de_genes(code_de_genes,
        cluster.out.scobj)
}

/*--------------PROCESSES----------------*/

process write_params{ // write input parameters to file
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        path "*.txt"
    """
    echo \$(date) > parameters_sc_p_sample_${params.id}.txt
    echo "Parameters used to analyze the seurat object for each sample:" >> parameters_sc_p_sample_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_sc_p_sample_${params.id}.txt
    """
}

process qc_unfiltered { // create QC plots for unfiltered seurat object
    publishDir params.outfolder+'4_seurat_unfiltered_qc', mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit: scobj
        path "*.pdf"
    """
    Rscript code.R --scobject=${scobj} --dataset='RAW'
    """
}

process filter_qc{ // filter the seurat object based on user parameters
    publishDir params.outfolder+'5_seurat_filtered',pattern:'*.pdf', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.pdf"
        tuple path("*.mtx"), path("*.rds"), emit: sc_mtx_obj
    """
    Rscript code.R --sc_file=${scobj} --mincount_p_gene='${params.mincount_p_gene}' --maxcount_p_gene='${params.maxcount_p_gene}' --mincell_p_gene='${params.mincell_p_gene}' --maxcell_p_gene='${params.maxcell_p_gene}' --mincount_p_cell='${params.mincount_p_cell}' --maxcount_p_cell='${params.maxcount_p_cell}' --mingene_p_cell='${params.mingene_p_cell}' --maxgene_p_cell='${params.maxgene_p_cell}' --maxpct_mt='${params.maxpct_mt}' --maxpct_rb='${params.maxpct_rb}' --rm_mt='${params.rm_mt}' --rm_rb='${params.rm_rb}'
    """
}

process scrublet_filter_doublet{
    publishDir params.outfolder+'6_seurat_filtered', pattern:"*.pdf", mode: 'copy', overwrite: true
    input:
        path "code.py"
        path "code.R"
        tuple path(sc_mtx), path(scobj)
    output:
        path "*.csv"
        path "*.rds", emit: scobj
        path "*.pdf"
    """
    python code.py ${sc_mtx} '${params.npcs}' '${params.exp_db_rate}'
    Rscript code.R --sc_file=${scobj} --threshold='${params.threshold}'
    """
}

process qc_filtered{
    publishDir params.outfolder+'7_seurat_filtered_qc', mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        path "code.R"
        path scobj
    output:
        path "*.pdf"
        path "*.rds", emit: scobj
    """
    Rscript code.R --scobject=${scobj} --dataset='Filtered'
    """
}

process normalize{
    errorStrategy 'ignore'
    publishDir params.outfolder+'8_seurat_normalized', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit: scobj
    """
    Rscript code.R --scobject=${scobj} --method='${params.method}' --nhvg='${params.nhvg}' --covars='${params.covars}' --genes_regress='${params.genes_regress}' --ncores='${params.ncores}'
    """
}

process dim_reduction{
    publishDir params.outfolder+'9_seurat_dim_red', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit: scobj
        path "*.pdf"
    """
    Rscript code.R --scobject=${scobj} --npcs='${params.npcs}' --ncores='${params.ncores}'
    """
}

process cluster{
    publishDir params.outfolder+'10_seurat_clustered', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path scobj
    output:
        path "*.rds", emit: scobj
        path "*.pdf"
    """
    Rscript code.R --scobject=${scobj} --res='${params.res}' --ncores='${params.ncores}'
    """
}

process de_genes{
    publishDir params.outfolder+'11_de', mode: 'copy', overwrite: true
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
