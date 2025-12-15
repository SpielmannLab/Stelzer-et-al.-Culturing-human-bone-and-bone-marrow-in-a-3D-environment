// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment ScRNA
// Submit example: nextflow run make_plots.nf -params-file make_plots_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------MAIN WORKFLOW----------------*/
workflow {
    write_params()

    sc_obj_channel = Channel.fromList(params.in_seurat_rds)
        .view()

    if (params.plot_cell_composition) {
        cell_composition(sc_obj_channel)
    }

    if (params.plot_cell_embedding) {
        cell_embedding(sc_obj_channel)
    }

    if (params.plot_genes_in_umap) {
        genes_in_umap(sc_obj_channel)
    }

    if (params.plot_genes_in_ridgeplot) {
        genes_in_ridgeplot(sc_obj_channel)
    }

    if (params.plot_genes_in_violinplot) {
        genes_in_violinplot(sc_obj_channel)
    }

    if (params.plot_genes_in_dotplot) {
        genes_in_dotplot(sc_obj_channel)
    }

    if (params.plot_gvg) {
        gene_vs_gene(sc_obj_channel, params.gvg_gene_pairs)
    }
}

/*--------------PROCESSES----------------*/

process write_params{ // write input parameters to file
    output:
        path "*.txt"
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_sc_p_sample_${suffix}.txt
    echo "!{params}" | tr , '\n' >> parameters_sc_p_sample_${suffix}.txt
    '''
}

process cell_composition {
    label "medium"
    input:
        path file_sc_obj
    output:
        path '*.tsv'
        path '*.pdf' optional true
    """
    cell_composition.R --file_sc_obj=${file_sc_obj}  --grouping_var="${params.grouping_var}" --xvar="${params.xvar}" --correct_grouping_var_imbalance="${params.correct_grouping_var_imbalance}" --figure_height="${params.figure_height}"
    """
}

process cell_embedding {
    label "medium"
    input:
        path file_sc_obj
    output:
        path '*.pdf'
    """
    cells_in_embedding.R --file_sc_obj=${file_sc_obj} --group_by=${params.cell_embedding_group_by} --color_values=${params.cell_embedding_color_values} --split_by=${params.cell_embedding_split_by}  --pt_size=${params.cell_embedding_pt_size} --shuffle_pts=${params.shuffle_pts} --width=${params.cell_embedding_width} --height=${params.cell_embedding_height} --reduction_to_use=${params.cell_embedding_reduction_to_use} --aspect_ratio=${params.cell_embedding_aspect_ratio}
    """
}

process genes_in_umap {
    label "medium"
    input:
        path file_sc_obj
    output:
        path '*.pdf'
    """
    genes_in_umap.R --file_sc_obj=${file_sc_obj} --genes=${params.genes} --split_by=${params.split_by} --assay=${params.assay} --genes_per_file=${params.umap_genes_per_file} --pt_size=${params.umap_pt_size} --cell_order=${params.umap_cell_order} --width=${params.umap_width} --height=${params.umap_height}
    """
}

process genes_in_violinplot {
    label "large"
    input:
        path file_sc_obj
    output:
        path '*.pdf'
    """
    genes_in_violinplot.R --file_sc_obj=${file_sc_obj} --genes=${params.genes} --group_by=${params.group_by} --assay=${params.assay} --color_values=${params.vln_color_values} --genes_per_file=${params.vln_genes_per_file} --pt_size=${params.vln_pt_size} --width=${params.vln_width} --height=${params.vln_height} --do_statistics=${params.vln_do_statistics}
    """
}

process genes_in_ridgeplot {
    label "medium"
    debug true
    input:
        path file_sc_obj
    output:
        path '*.pdf'
    """
    genes_in_ridgeplot.R \
        --file_sc_obj=${file_sc_obj} \
        --genes=${params.genes} \
        --group_by=${params.group_by} \
        --split_by=${params.split_by} \
        --sort_splits_by=${params.sort_splits_by} \
        --assay=${params.assay} \
        --color_values=${params.ridge_color_values} \
        --bandwidth=${params.ridge_bandwidth} \
        --genes_per_file=${params.ridge_genes_per_file} \
        --width=${params.ridge_width} \
        --height=${params.ridge_height}
    """
}

process genes_in_dotplot {
    label "medium"
    input:
        path file_sc_obj
    output:
        path '*_DotPlot.pdf'
    """
    genes_in_dotplot.R --file_sc_obj=${file_sc_obj} --genes=${params.genes} --group_by=${params.group_by} --split_by=${params.split_by} --assay=${params.assay} --width=${params.dot_width} --height=${params.dot_height} --angle_x_text=${params.angle_x_text}
    """
}

process gene_vs_gene {
    label "medium"
    input:
        path file_sc_obj
        each gene_pair
    output:
        path '*_gvg_*.pdf'
    """
    gene_vs_gene.R --file_sc_obj=${file_sc_obj} --assay=${params.gvg_assay} --gene_pair=${gene_pair} --group_by=${params.gvg_group_by} --pt_size=${params.gvg_pt_size} --width=${params.gvg_width} --height=${params.gvg_height}
    """
}
