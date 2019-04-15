### add barcodes from R1 to R2 #########

rule fastq_barcode:
        input: ## remember that we swapped reads[] in internals.snakefile in this workflow!!!
            R2 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R1 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R2_barcoded = "FASTQ_barcoded/{sample}"+reads[0]+".fastq.gz"
        params:
            UMI_length = UMI_length,
            UMI_offset = UMI_offset,
            CELLI_length = CELLI_length,
            CELLI_offset = CELLI_offset
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell:"""
            paste <(paste - - - - < <(zcat {input.R1}))   <(paste - - - - < <(zcat {input.R2})) | \
            tr '\t' '\n' | \
            awk -v CBAR_LEN={params.CELLI_length} -v UMI_LEN={params.UMI_length} \
                -v CBAR_OFFSET={params.CELLI_offset} -v UMI_OFFSET={params.UMI_offset} ' \
           	BEGIN{{
				for(n=0;n<256;n++) 
					phred33[sprintf("%c",n)]=n-33
 			}}
			{{
			if (NR%8==2) 
		 		{{CB=substr($0,CBAR_OFFSET,CBAR_LEN);
	 			UMI=substr($0,UMI_OFFSET,UMI_LEN);
	 			if (CBAR_OFFSET+CBAR_LEN>UMI_OFFSET+UMI_LEN)
	 			TAIL=substr($0,CBAR_OFFSET+CBAR_LEN); else
	 			TAIL=substr($0,UMI_OFFSET+UMI_LEN);
	 			NUMT=gsub(/T/,"#",TAIL);
			}} 
			if (NR%8==4) 
				{{split(substr($0,UMI_LEN+1,CBAR_LEN),CBQA,"");
	 			split(substr($0,1,UMI_LEN),UMIQA,""); 
	 			QUAL_UMI=0; QUAL_CB=0;
 	 			for (i=1;i<=length(UMIQA);i++)
					{{QUAL_UMI+=phred33[UMIQA[i]]}}; 
	 			for (i=1;i<=length(CBQA);i++)
					{{QUAL_CB+=phred33[CBQA[i]]}};
				}}; 
			OFS=" ";
			if (NR%8==5) 
				{{$1=$1":SC:"CB":"sprintf("%.0f",QUAL_CB/CBAR_LEN)":UMI:"UMI":"sprintf("%.0f",QUAL_UMI/UMI_LEN)":"NUMT":"length(TAIL);print $0;
				}}; 
			if (NR%8==0 || NR%8>5) print $0}}' | pigz -c -p 8 > {output.R2_barcoded}
            """
    
    
rule sc_bam_featureCounts_genomic:
    input:
        bam = mapping_prg+"/{sample}.bam",
        gtf = "Annotation/genes.filtered.gtf"
    output:
        counts = "Counts/{sample}.raw_counts.txt",
        counts_summary = "Counts/{sample}.featureCounts_summary.txt"
    params:
        count_script = workflow.basedir+"/scRNAseq_bam_featureCounts.sh",
        bc_file = barcode_file,
        lib_type = library_type
    threads: 
        5
    conda: CONDA_RNASEQ_ENV
    shell:
        """
        {params.count_script} {input.bam} {input.gtf} {params.bc_file} {wildcards.sample} {params.lib_type} ${{TMPDIR}} {threads} 1>{output.counts} 2>{output.counts_summary};       
        """


rule extract_scale_counts:
    input:
        counts = "Counts/{sample}.raw_counts.txt"
    output:
        corrected = "Counts/{sample}.corrected.txt",
        umis = "Counts/{sample}.umis.txt",
        reads = "Counts/{sample}.reads.txt"
    params:
        count_script = os.path.join(maindir, "shared", "tools", "correct_sc_counts.py"),
        UMI_length = UMI_length
    log:
        out = "Counts/logs/extract_counts.{sample}.out",
        err = "Counts/logs/extract_counts.{sample}.err"
    conda: CONDA_RNASEQ_ENV
    shell: """
        {params.count_script} --umiLength {params.UMI_length} \
            {input.counts} {output.reads} {output.umis} {output.corrected} > {log.out} 2> {log.err}
        """


rule combine_sample_counts:
    input:
        expand("Counts/{sample}.corrected.txt",sample = samples),
    output:
        merged_matrix = "Results/all_samples.gencode_genomic.corrected_merged.csv",
        used_cell_names_file = "Results/all_samples.used_cells.tsv"
    params:
        merge_script = workflow.basedir+"/scRNAseq_merge_coutt_files2.R",
        split = split_lib,
        sample_cell_names = str(cell_names or '')
    conda: CONDA_scRNASEQ_ENV
    shell:
        "Rscript {params.merge_script} Counts/ {output.merged_matrix} {output.used_cell_names_file} {params.split} {params.sample_cell_names} """


rule filter_cells_monocle:
        input: 
            merged_matrix = "Results/all_samples.gencode_genomic.corrected_merged.csv"
        output:
            metrics_tab = "Filtered_cells_monocle/metrics.tab.RData"
        params:
            wdir=os.path.join(outdir,"Filtered_cells_monocle"),
            fINpath=lambda wildcards,input: os.path.join(outdir,input.merged_matrix)
        log:
            err='Filtered_cells_monocle/logs/filt_cells.err',
            out='Filtered_cells_monocle/logs/filt_cells.out'
        threads: 1
        conda: CONDA_scRNASEQ_ENV
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'scRNAseq_cell_filter_monocle.R ') + "{params.wdir} {params.fINpath} 1>{log.out} 2>{log.err}"


rule cluster_cells_monocle:
        input: 
            metrics_tab = "Filtered_cells_monocle/metrics.tab.RData"
        output:
            dummy="Filtered_cells_monocle/sessionInfo.txt"
        params:
            wdir=os.path.join(outdir,"Filtered_cells_monocle"),
            fINpath=lambda wildcards,input: os.path.join(outdir,input.metrics_tab),
            err='Filtered_cells_monocle/logs/cluster_cells.err',
            out='Filtered_cells_monocle/logs/cluster_cells.out',
            metric=cell_filter_metric
        threads: 1
        conda: CONDA_scRNASEQ_ENV
        shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'scRNAseq_select_threshold_cluster_monocle.R ') + "{params.wdir} {params.fINpath} {params.metric} 1>{params.out} 2>{params.err}"

rule monocle_report:
    input: 
        dummy="Filtered_cells_monocle/sessionInfo.txt"
    output:
        html='Filtered_cells_monocle/Stats_report.html'
    params:
        statdir=os.path.join(outdir,"Filtered_cells_monocle"),
        rmd_in=os.path.join(workflow_rscripts,"scRNAseq_monocle_stats_report.Rmd"),
        rmd_out=os.path.join(outdir, "scRNAseq_monocle_stats_report.Rmd"),
        outFull=lambda wildcards,output: os.path.join(outdir,output.html),
        metric=cell_filter_metric
    log:
        err='Filtered_cells_monocle/logs/stats_report.err',
        out='Filtered_cells_monocle/logs/stats_report.out'
    conda: CONDA_RMD_ENV
    threads: 1
    shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\",metric=\"{params.metric}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"


if not skipRaceID:
    rule filter_cells_raceid:
            input: 
                merged_matrix = "Results/all_samples.gencode_genomic.corrected_merged.csv"
            output:
                metrics_tab = "Filtered_cells_RaceID/metrics.tab.RData"
            params:
                wdir=os.path.join(outdir,"Filtered_cells_RaceID"),
                fINpath=lambda wildcards,input: os.path.join(outdir,input.merged_matrix)
            log:
                err='Filtered_cells_RaceID/logs/filt_cells.err',
                out='Filtered_cells_RaceID/logs/filt_cells.out'
            threads: 1
            conda: CONDA_scRNASEQ_ENV
            shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'scRNAseq_cell_filter_raceid.R ') + "{params.wdir} {params.fINpath} 1>{log.out} 2>{log.err}"


    rule cluster_cells_raceid:
            input: 
                metrics_tab = "Filtered_cells_RaceID/metrics.tab.RData"
            output:
                dummy="Filtered_cells_RaceID/sessionInfo.txt"
            params:
                wdir=os.path.join(outdir,"Filtered_cells_RaceID"),
                fINpath=lambda wildcards,input: os.path.join(outdir,input.metrics_tab),
                err='Filtered_cells_RaceID/logs/cluster_cells.err',
                out='Filtered_cells_RaceID/logs/cluster_cells.out',
                metric=cell_filter_metric
            threads: 1
            conda: CONDA_scRNASEQ_ENV
            shell: "Rscript --no-save --no-restore " + os.path.join(workflow_rscripts,'scRNAseq_select_threshold_cluster_raceid.R ') + "{params.wdir} {params.fINpath} {params.metric} 1>{params.out} 2>{params.err}"


    rule raceid_report:
        input: 
            dummy="Filtered_cells_RaceID/sessionInfo.txt"
        output:
            html='Filtered_cells_RaceID/Stats_report.html'
        params:
            statdir=os.path.join(outdir,"Filtered_cells_RaceID"),
            rmd_in=os.path.join(workflow_rscripts,"scRNAseq_raceid_stats_report.Rmd"),
            rmd_out=os.path.join(outdir, "scRNAseq_raceid_stats_report.Rmd"),
            outFull=lambda wildcards,output: os.path.join(outdir,output.html),
            metric=cell_filter_metric
        log:
            err='Filtered_cells_RaceID/logs/stats_report.err',
            out='Filtered_cells_RaceID/logs/stats_report.out'
        conda: CONDA_RMD_ENV
        threads: 1
        shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\",metric=\"{params.metric}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"


rule sc_QC_metrics:
    input:
        expand("Counts/{sample}.featureCounts_summary.txt",sample=samples),
        cell_names_merged = "Results/all_samples.used_cells.tsv"
    output:
        summary = "QC_report/QC_report.all_samples.libstats_reads.tsv"
    params: 
        in_dir = outdir+"/Counts/",
        cellsum_dir = "QC_report/data/",
        out_dir = outdir+"/QC_report/",
        plot_script = workflow.basedir+"/scRNAseq_QC_metrics2.R",
        out_prefix = "QC_report/QC_report.all_samples",
        plot_format = plot_format,
        split = split_lib
    conda: CONDA_RNASEQ_ENV
    shell:
        ""+workflow.basedir+"/scRNAseq_QC_metrics.sh {params.in_dir} {params.out_dir} >{output.summary};"
        " Rscript {params.plot_script} {params.cellsum_dir} {params.out_prefix} {params.split} {input.cell_names_merged} {params.plot_format}"
