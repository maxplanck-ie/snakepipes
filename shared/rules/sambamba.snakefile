##sambamba is used for marking up duplications

rule sambamba_markdup:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           mapping_prg+"/{sample}.bam"
       threads: 10       
       shell:
           sambamba_path+"sambamba_v0.6.6 "
           "markdup -t {threads} --sort-buffer-size=6000 "
           "{input} "
           "{output}"

rule sambamba_flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.dup.txt"
       shell:
           sambamba_path+"sambamba_v0.6.6 flagstat -p"
           " {input}"
           " | tee {output}"
rule convert_flagstat_output:
      input: 
         "Sambamba/{sample}.dup.txt"
      output:
         "Sambamba/{sample}.dup.converted.tsv"
      run:
         shell("sed -n '1p;4p;5p' {input} > Sambamba/{wildcards.sample}.temp"), ##TODO Change it to a 'while' loop
         shell("paste -s -d'\t' <(awk 'FNR==1' Sambamba/{wildcards.sample}.temp | cut -d' ' -f1) <(awk 'FNR==1' Sambamba/{wildcards.sample}.temp | cut -d' ' -f3) <(awk 'FNR==2' Sambamba/{wildcards.sample}.temp | cut -d' ' -f1) <(awk 'FNR==2' Sambamba/{wildcards.sample}.temp | cut -d' ' -f3)  <(awk 'FNR==3' Sambamba/{wildcards.sample}.temp | cut -d' ' -f1)  <(awk 'FNR==3' Sambamba/{wildcards.sample}.temp | cut -d' ' -f3) | tr '\n' '\t' | sed 's/^/{wildcards.sample}\t/' > Sambamba/{wildcards.sample}.temp1"),
         shell(" sed -e '$a\\' Sambamba/{wildcards.sample}.temp1 > {output}")
         shell("rm Sambamba/{wildcards.sample}.temp"),
         shell("rm Sambamba/{wildcards.sample}.temp1"),

####sambamba flagstat is run here######
rule report_flagstat_all_data:
      input:
         expand("Sambamba/{sample}.dup.converted.tsv",sample=samples)
      output:
         "Sambamba/flagstat_report_all.tsv"
      run:
         header = (  "Sample\t"
                     "Total QC-passed reads(P)\t"
                     " Total QC-failed reads(F)\t"
                     " Duplicated P\t"
                     " Duplicated F\t"
                     " Mapped P\t"
                     " Mapped F"
                   )
         shell('cat {input} >temp'),
         shell('cat temp | cat <(echo "'+header+'") - | tee {output}')

