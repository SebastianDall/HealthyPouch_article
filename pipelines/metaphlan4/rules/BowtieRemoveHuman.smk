rule Bowtie2:
    input:
        R1     = OUT_DIR+"{barcode}/fastp/{sample}_R1_fastp.fastq",
        R2     = OUT_DIR+"{barcode}/fastp/{sample}_R2_fastp.fastq"
    resources:
        walltime = "9800",
        nodetype = "thinnode",
        mem = "20"
    threads:
        8
    params:
        bowtie_index = "databases/GRCh38.d1.vd1/bowtie2build/230302_bowtie2build_2.5.0_GRCh38.d1.vd1",
        tempdir = TMP_DIR+"tmp/{barcode}/{sample}/bowtie_tmp/",
        very_sensitive = "--very-sensitive",
    output:
        bam     = temp(TMP_DIR+"tmp/{barcode}/{sample}_NonHuman.bam"),
    log: 
        bowtie2  = OUT_DIR+"{barcode}/log/{sample}_bowtie2_"+getCurrentTime()+".log",
        view     = OUT_DIR+"{barcode}/log/{sample}_samtoolsView_"+getCurrentTime()+".log",
        sort     = OUT_DIR+"{barcode}/log/{sample}_samtoolsSort_"+getCurrentTime()+".log",
    envmodules:
        "tools",
        "samtools/1.16",
        "bowtie2/2.5.0"
    shell:
        '''
        mkdir -p {params.tempdir} && \
        bowtie2 \
        -x {params.bowtie_index} \
        {params.very_sensitive} \
        --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log.bowtie2} | \
        samtools view -hb -@ {threads} -S - | \
        samtools view -hb -@ {threads} -f 12 -F 256 2> {log.view} | \
        samtools sort -n -o {output.bam} -T {params.tempdir} -o {output}  2> {log.sort} && \
        rm -r {params.tempdir}
        '''

rule SamtoolsFastq:
    input:
        bam     = TMP_DIR+"tmp/{barcode}/{sample}_NonHuman.bam"
    threads:
        8
    resources:
        walltime = "9800",
        nodetype = "thinnode",
        mem = "20"
    output:
        R1      = temp(TMP_DIR+"tmp/{barcode}/{sample}_R1_NonHuman.fastq"),
        R2      = temp(TMP_DIR+"tmp/{barcode}/{sample}_R2_NonHuman.fastq")
    log: 
        samtoolsFastq  = OUT_DIR+"{barcode}/log/{sample}_samtoolsFastq_"+getCurrentTime()+".log"
    envmodules:
        "tools",
        "samtools/1.16"
    shell:
        '''
        samtools fastq \
        -1 {output.R1} \
        -2 {output.R2} \
        -0 /dev/null -s /dev/null \
        -n -@ {threads} {input.bam} &> {log.samtoolsFastq}
        '''

rule CompressFastq:
    input:
        R1      = TMP_DIR+"tmp/{barcode}/{sample}_R1_NonHuman.fastq",
        R2      = TMP_DIR+"tmp/{barcode}/{sample}_R2_NonHuman.fastq"
    threads:
        2
    resources:
        walltime = "9800",
        nodetype = "thinnode",
        mem = "20"
    output:
        R1      = protected(OUT_DIR+"{barcode}/{sample}_R1_NonHuman.fastq.gz"),
        R2      = protected(OUT_DIR+"{barcode}/{sample}_R2_NonHuman.fastq.gz")
    shell:
        '''
        gzip -cvf {input.R1} > {output.R1} && \
        gzip -cvf {input.R2} > {output.R2}
        '''


