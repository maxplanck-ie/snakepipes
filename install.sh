
## at this point the script simply links the workflows in the main dir, but this is
## just a placeholder to make the point that the installation is coming in future

for workflow in ChIP-seq RNA-seq HiC DNA-mapping ATAC-seq WGBS ; do
ln -s workflows/${workflow}/${workflow} . &&
ln -s shared/shared_environment.yaml workflows/${workflow}/shared_environment.yaml
done
# seperate link for scRNAseq
workflow=scRNAseq
ln -s workflows/${workflow}/${workflow}-mapcount .
ln -s shared/shared_environment.yaml workflows/scRNAseq/shared_environment.yaml
