.libPaths(R.home("library"))

logfile = file(snakemake@log[["err"]], open="w")
sink(logfile, type="message")
snakemake@source(snakemake@params[["importfunc"]])

# Read in the sample sheet
options(stringsAsFactors=FALSE,na.rm=TRUE)
sampleSheet = read.table(snakemake@params[["sampleSheet"]], sep="\t", as.is=TRUE, header=TRUE)
if (!"PlottingID" %in% colnames(sampleSheet)) {
    sampleSheet$PlottingID<-sampleSheet$name
}

# Read in the expression files iff they're listed in sampleSheet
expressionFiles = snakemake@input[["expressionFiles"]]
expressionDFs = lapply(sampleSheet$name, function(sampleName) {
    idx = which(endsWith(expressionFiles, sprintf("/%s.CpG.filt2.bed", sampleName)))[1]
    fname = expressionFiles[idx]
    tab = read.delim(fname, header=TRUE)
    tab = tab[,c("Beta", "ms")]
    colnames(tab)[1] = sampleName
    return(tab)
})

# Merge files
fullMatrix = Reduce(function(...) merge(..., all=T, by="ms", sort=FALSE), expressionDFs)
rm(expressionDFs)

if (nrow(fullMatrix) == 0) {
    print_sessionInfo("None of the single CpG sites passed the filtering.")
} else {
    # reorder colums so group members are next to each other
    # Treatment or WT should go first, otherwise lexographic ordering
    cnames = colnames(fullMatrix)
    if ("Mut" %in% sampleSheet$condition) {
        WT = which(!sampleSheet$condition == "Mut")
        Mut = which(sampleSheet$condition == "Mut")
        fullMatrix = fullMatrix[,c("ms", cnames[WT + 1], cnames[Mut + 1])]
        sampleSheet = sampleSheet[c(WT, Mut),]
    } else if ("Treatment" %in% sampleSheet$condition) {
        CTRL = which(!sampleSheet$condition == "Treatment")
        Treatment = which(sampleSheet$condition == "Treatment")
        fullMatrix = fullMatrix[,c("ms", cnames[CTRL + 1], cnames[Treatment + 1])]
        sampleSheet = sampleSheet[c(CTRL, Treatment),]
    } else {
        baseLevel = levels(factor(sampleSheet$condition))[1]
        CTRL = which(!sampleSheet$condition == baseLevel)
        Treatment = which(sampleSheet$condition == baseLevel)
        fullMatrix = fullMatrix[,c("ms", cnames[CTRL + 1], cnames[Treatment + 1])]
        sampleSheet = sampleSheet[c(CTRL, Treatment),]
    }

    fullMatrix$chr = gsub("_.+", "", fullMatrix$ms)
    fullMatrix$pos = gsub(".+_", "", fullMatrix$ms)
    fullMatrix = fullMatrix[,c(ncol(fullMatrix) - 1, ncol(fullMatrix), 2:(ncol(fullMatrix)-2))]

    # Change the column names to sample_group
    snames = paste(sampleSheet$name, sampleSheet$condition, sep="_")
    colnames(fullMatrix) = c("chr", "pos", snames)

    write.table(fullMatrix, file=snakemake@output[["MetileneIN"]], sep="\t", row.names=FALSE, quote=FALSE)
}

print_sessionInfo("Analysis completed succesfully.")
sink(type="message")
close(logfile)
