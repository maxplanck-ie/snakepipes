
## at this point the script simply links the workflows in the main dir, but this is
## just a placeholder to make the point that the installation is coming in future

for workflow in ChIP-seq RNA-seq HiC DNA-mapping ATAC-seq WGBS ; do
ln -s workflows/${workflow}/${workflow} .
done
# seperate link for scRNAseq
workflow=scRNAseq
ln -s workflows/${workflow}/${workflow}-mapcount .
