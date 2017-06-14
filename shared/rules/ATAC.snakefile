rule callOpenChromatin:
    input:
        bam = None
    output:
        pass
    params:
        pass
    threads: 6
    log:
        pass
    shell: # or run:
        ## macs2
        macs2_path+'macs2 callpeak'

rule sortByName:
    input:
        bam
    output:
        bam
    params:
        byQuery='-n'
    threads: 6
    log:
        pass
    shell:
        samtools_path+"samtools sort -n -@ {threads} {input} -o {output}"

rule reads2cutfragments:
    input:
        bam
    output:
        bedpe
    params:
        pass
    threads: 1
    log:
        pass
    shell:
        pass

rule filterFragments:
    input:
        bedpe
    output:
        bedpe
    params:
        cutoff = 147
    threads:1
    log:
        pass
    shell:
        pass

rule fragmentSizeDistribution:
    input: bedpe
    output: fragdistr
    params: pass
    threads: 1
    log: pass
    shell: pass
