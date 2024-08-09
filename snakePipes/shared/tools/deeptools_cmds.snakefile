
####---- Common deeptools commands used for multiple workflows ------#####


# bamcompare
bamcompare_log2_cmd = """
    bamCompare -b1 {input.chip_bam} \
               -b2 {input.control_bam} \
               -o {output} \
               --operation log2 \
               {params.scaleFactors} \
               {params.ignoreForNorm} \
               --binSize {params.bwBinSize} \
               -p {threads} \
               {params.read_extension} \
               {params.blacklist}
    """

# bamcompare subtract
bamcompare_subtract_cmd = """
    bamCompare -b1 {input.chip_bam} \
               -b2 {input.control_bam} \
               -o {output} \
               --operation subtract \
               {params.scaleFactors} \
               {params.ignoreForNorm} \
               --binSize {params.bwBinSize} \
               -p {threads} \
               {params.read_extension} \
               {params.blacklist}
    """

# bamCoverage RAW
bamcov_raw_cmd = """
    bamCoverage -b {input.bam} \
                -o {output} \
                --binSize {params.bwBinSize} \
                -p {threads}
    """

# bamCoverage RPKM
bamcov_RPKM_cmd = """
    bamCoverage -b {input.bam} \
        -o {output} --binSize {params.bwBinSize} \
        -p {threads} --normalizeUsing RPKM {params.ignoreForNorm} \
        {params.blacklist}
    """

# bamCoverage three-prime-seq non-unique, strand-specific mappings
# samflagExclude 128 == exclude second read in pair
# uses params.direction == forward or reverse
bamcov_3pseq_cmd = """
    bamCoverage -b {input.bam} \
        -o {output.bw_fwd} --binSize {params.bwBinSize} \
        --Offset 1 --samFlagExclude 128 --filterRNAstrand {params.direction} \
        -p {threads} --skipNAs
    """


# bamCoverage RNA-seq unique mappings
bamcov_unique_cmd = """
    bamCoverage -b {input.bam} \
        -o {output.bw_fwd} --binSize {params.bwBinSize} \
        --minMappingQuality 10 --samFlagExclude 2304 --filterRNAstrand forward \
        -p {threads}
    bamCoverage -b {input.bam} \
        -o {output.bw_rev} --binSize {params.bwBinSize} \
        --minMappingQuality 10 --samFlagExclude 2304 --filterRNAstrand reverse \
        -p {threads}
    """


# bamCoverage CHIP
bamcov_cmd = """
    bamCoverage -b {input.bam} \
                -o {output} \
                --binSize {params.bwBinSize} \
                -p {threads} \
                --normalizeUsing RPGC \
                --effectiveGenomeSize {params.genome_size} \
                {params.ignoreForNorm} \
                {params.blacklist} \
                {params.read_extension}
"""

bamcov_spikein_cmd = """
    bamCoverage -b {input.bam} \
                -o {output} \
                --binSize {params.bwBinSize} \
                -p {threads} \
                --normalizeUsing None \
                --effectiveGenomeSize {params.genome_size} \
                {params.ignoreForNorm} \
                {params.blacklist} \
                {params.scaling_factors} \
                {params.read_extension}
    """

## computeGC bias (DNA), requires params.median_fragment_length
gcbias_cmd = """
    computeGCBias -b {input.bam} \
                --biasPlot {output.png} \
                --GCbiasFrequenciesFile {output.tsv} \
                --effectiveGenomeSize {params.genome_size} \
                --genome {params.genome_2bit} \
                {params.median_fragment_length}  \
                --sampleSize {params.sampleSize} \
                {params.blacklist} \
                -p {threads}
    """

# plot Enrichment (RNAseq)
plotEnrich_cmd = """
    plotEnrichment -p {threads} \
                   -b {input.bam} \
                   {params.plotcmd} \
                   --BED {input.gtf} \
                   --labels {params.labels} \
                   --plotTitle 'Fraction of reads in regions' \
                   --outRawCounts {output} \
                   --variableScales
    """

# plot Enrichment (ChIPSeq)
plotEnrich_chip_cmd = """
    plotEnrichment \
        -b {input.bams} \
        --BED {params.genes_gtf} \
        --plotFile {output.png} \
        {params.labels} \
        --plotTitle 'Signal enrichment (fraction of reads) without duplicates' \
        --outRawCounts {output.tsv} \
        --variableScales \
        {params.blacklist} \
        -p {threads} \
        {params.read_extension} \
        --ignoreDuplicates
    """

#plot fingerprint (ChIPseq)
plotFingerprint_cmd = """
    plotFingerprint \
            -b {input.bams} \
            {params.labels} \
            --plotTitle 'Cumulative read counts per bin without duplicates' \
            --ignoreDuplicates \
            --outQualityMetrics {output.metrics} \
            -p {threads} \
            {params.blacklist} \
            {params.png} \
            {params.read_extension} \
            {params.jsd}
    """


# multiBAMsum ChIP
multiBamSummary_cmd = """
    multiBamSummary bins \
                    -b {input.bams} \
                    -o {output.npz} \
                    --labels {params.labels} \
                    {params.blacklist} \
                    {params.scaling_factors} \
                    {params.binSize} \
                    {params.spikein_region} \
                    -p {threads} \
                    {params.read_extension}
    """

# multiBAMsum ChIP with spikein
multiBamSummary_spikein_cmd = """
    multiBamSummary BED-file \
                    --BED {input.bed} \
                    -b {input.bams} \
                    -o {output.npz} \
                    --labels {params.labels} \
                    {params.blacklist} \
                    {params.scaling_factors} \
                    {params.binSize} \
                    -p {threads} \
                    {params.read_extension}
    """

# multiBWsum RNA
multiBWsum_bed_cmd = """
    multiBigwigSummary BED-file \
                --BED {input.bed} \
                -b {input.bw} \
                -o {output} \
                --labels {params.labels} \
                --binSize 1000 \
                -p {threads}
    """

# multiBamSum RNA
multiBamSum_bed_cmd = """
    multiBamSummary BED-file \
                --BED {input.bed} \
                -b {input.bam} \
                -o {output.npz} \
                --labels {params.labels} \
                --binSize 100 \
                --scalingFactors {output.scalingFactors} \
                {params.blacklist} \
                -p {threads}
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
                --plotNumbers
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
        --plotNumbers
    """

# plot PCA (both), requires params.label
plotPCA_cmd = """
    plotPCA -in {input} \
            {params.plotcmd} \
            --transpose \
            --outFileNameData {output} \
            --plotWidth 10 \
            -T 'PCA of {params.title} coverage'
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
                 --ignoreDuplicates
    """

#EstimateReadFiltering
estimateReadFiltering_cmd = """
    estimateReadFiltering -b {input.bam} \
        -o {output}
    """

#bamPEFragmentSize
bamPEFragmentSize_cmd = """
    bamPEFragmentSize \
        --bamfiles {input.bams} \
        --binSize 1000000 \
        {params.plotcmd} \
        --table {output} -p {threads}
    """
