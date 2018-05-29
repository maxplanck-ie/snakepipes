
####---- Common deeptools commands used for multiple workflows ------#####


# bamcompare
bamcompare_log2_cmd = """
    bamCompare -b1 {input.chip_bam} \
               -b2 {input.control_bam} \
               -o {output} \
               --operation log2 \
               --scaleFactorsMethod readCount \
               {params.ignoreForNorm} \
               --binSize {params.bw_binsize} \
               -p {threads} \
               {params.read_extension} \
               {params.blacklist} > {log.out} 2> {log.err}
    """

# bamcompare subtract
bamcompare_subtract_cmd = """
    bamCompare -b1 {input.chip_bam} \
               -b2 {input.control_bam} \
               -o {output} \
               --operation subtract \
               --scaleFactorsMethod readCount \
               --normalizeUsing RPGC \
               --effectiveGenomeSize {params.genome_size} \
               {params.ignoreForNorm} \
               --binSize {params.bw_binsize} \
               -p {threads} \
               {params.read_extension} \
               {params.blacklist}  > {log.out} 2> {log.err}
    """

# bamCoverage RAW
bamcov_raw_cmd = """
    bamCoverage -b {input.bam} \
                -o {output} \
                --binSize {params.bw_binsize} \
                -p {threads}  > {log.out} 2> {log.err}
    """

# bamCoverage RPKM
bamcov_RPKM_cmd = """
    bamCoverage -b {input.bam} -o {output} --binSize {params.bw_binsize} -p {threads} --normalizeUsing RPKM  > {log.out} 2> {log.err}
    """

# bamCoverage CHIP
bamcov_cmd = """
    bamCoverage -b {input.bam} \
                -o {output} \
                --binSize {params.bw_binsize} \
                -p {threads} \
                --normalizeUsing RPGC \
                --effectiveGenomeSize {params.genome_size} \
                {params.ignoreForNorm} \
                {params.blacklist} \
                {params.read_extension}  > {log.out} 2> {log.err}
    """

## computeGC bias (DNA), requires params.median_fragment_length
gcbias_cmd = """
    computeGCBias -b {input.bam} \
                --biasPlot {output.png} \
                --GCbiasFrequenciesFile {output.tsv} \
                --effectiveGenomeSize {params.genome_size} \
                --genome {params.genome_2bit} \
                {params.median_fragment_length}  \
                --sampleSize 10000000 \
                {params.blacklist} \
                -p {threads}  > {log.out} 2> {log.err}
    """

# plot Enrichment (RNAseq)
plotEnrich_cmd = """
    plotEnrichment -p {threads} \
                   -b {input.bam} \
                   --BED {input.gtf} {input.gtf2} \
                   --plotFile {output.png} \
                   --labels {params.labels} \
                   --plotTitle 'Fraction of reads in regions' \
                   --outRawCounts {output.tsv} \
                   --variableScales  > {log.out} 2> {log.err}
    """

# plot Enrichment (ChIPSeq)
plotEnrich_chip_cmd = """
    plotEnrichment \
        -b {input.bams} \
        --BED {params.genes_gtf} \
        --plotFile {output.png} \
        --labels {params.labels} \
        --plotTitle 'Sigal enrichment (fraction of reads) without duplicates' \
        --outRawCounts {output.tsv} \
        --variableScales \
        {params.blacklist} \
        -p {threads} \
        {params.read_extension} \
        --ignoreDuplicates  > {log.out} 2> {log.err}
    """

#plot fingerprint (ChIP-seq)
plotFingerprint_cmd = """
    plotFingerprint \
            -b {input.bams} \
            --labels {params.labels} \
            --plotTitle 'Cumulative read counts per bin without duplicates' \
            --ignoreDuplicates \
            --outQualityMetrics {output.metrics} \
            -p {threads} \
            {params.blacklist} \
            {params.png} \
            {params.read_extension} \
            {params.jsd}  > {log.out} 2> {log.err}
    """


# multiBAMsum ChIP
multiBamSummary_cmd = """
    multiBamSummary bins \
                    -b {input.bams} \
                    -o {output} \
                    --labels {params.labels} \
                    --binSize 1000 \
                    {params.blacklist} \
                    -p {threads} \
                    {params.read_extension}  > {log.out} 2> {log.err}
    """

# multiBWsum RNA
multiBWsum_bed_cmd = """
    multiBigwigSummary BED-file \
                --BED {input.bed} \
                -b {input.bw} \
                -o {output} \
                --labels {params.labels} \
                --binSize 1000 \
                -p {threads}  > {log.out} 2> {log.err}
    """

## plot Corr (both), requires params.label
plotCorr_cmd = """
    plotCorrelation \
                -in {input} \
                {params.plotcmd} \
                --corMethod pearson \
                --whatToPlot heatmap \
                --skipZeros \
                --plotTitle 'Pearson correlation of {params.title} coverage' \
                --outFileCorMatrix {output} \
                --colorMap PuBuGn \
                --plotNumbers  > {log.out} 2> {log.err}
    """

## plot Corr Spearman (both), requires params.label
plotCorrSP_cmd = """
    plotCorrelation \
        -in {input} \
        {params.plotcmd} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotTitle 'Spearman correlation of {params.title} coverage' \
        --outFileCorMatrix {output} \
        --colorMap PuBuGn \
        --plotNumbers  > {log.out} 2> {log.err}
    """

# plot PCA (both), requires params.label
plotPCA_cmd = """
    plotPCA -in {input} \
            {params.plotcmd} \
            --transpose \
            --outFileNameData {output} \
            -T 'PCA of {params.title} coverage'  > {log.out} 2> {log.err}
    """

# plot Coverage
plotCoverage_cmd = """
    plotCoverage -b {input.bams} \
                 {params.plotcmd} \
                 --outRawCounts {output} \
                 --labels {params.labels} \
                 --plotTitle 'Genome fragment coverage without duplicates' \
                 -p {threads} \
                 {params.read_extension} \
                 --ignoreDuplicates  > {log.out} 2> {log.err}
    """

#EstimateReadFiltering
estimateReadFiltering_cmd = """
    estimateReadFiltering -b {input.bam} -o {output}  > {log.out} 2> {log.err}
    """

#bamPEFragmentSize
bamPEFragmentSize_cmd = """
    bamPEFragmentSize --bamfiles {input.bams} \
    {params.plotcmd} \
    --table {output} -p {threads} > {log.out} 2> {log.err}
    """
