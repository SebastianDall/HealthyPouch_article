rule CountReadsInFastq:
    input:
        fastq      = OUT_DIR+"{barcode}/{sample}_R1_NonHuman.fastq.gz"
    threads:
        2
    resources:
        walltime = "1200",
        nodetype = "thinnode",
        mem = "4G"
    output:
        txt      = OUT_DIR+"{barcode}/{sample}_NonHuman_ReadCount.txt"


    shell:
        '''
        if  [ -s {input.fastq} ]
        then
            Count=$(echo $(zcat {input.fastq} | wc -l)/4|bc)
        else
            Count=0
        fi
        echo -e "{wildcards.barcode}\t{wildcards.sample}\t${{Count}}" > {output.txt}
        '''



rule CollectReadCounts:
    input:
        files = GetReadCountFiles#expand(OUT_DIR+"{barcode}/{sample}_NonHuman_ReadCount.txt", zip, barcode=BARCODES, sample=SAMPLES)
    params:
        files = GetReadCountFiles#expand(OUT_DIR+"{barcode}/{sample}_NonHuman_ReadCount.txt", zip, barcode=BARCODES, sample=SAMPLES)
    output:
        txt = OUT_DIR+"fastqs_NonHuman_ReadCount.txt"
    resources:
        walltime = "7200",
        nodetype = "thinnode",
        mem = "4G",
    threads: 1
    run:
        with open(output.txt, "w") as out:
            for file in params.files:
                with open(file, 'r') as f:
                    for line in f:
                        out.write(line.strip()+"\n")
                f.close()




rule ConcatenateFastqFiles:
    input:
        R1      = GetSampleFastqFilesR1,
        R2      = GetSampleFastqFilesR2
    output:
        R1      = OUT_DIR+"{barcode}/{barcode}_R1_NonHuman_Combined.fastq.gz",
        R2      = OUT_DIR+"{barcode}/{barcode}_R2_NonHuman_Combined.fastq.gz"
    resources:
        walltime = "7200",
        nodetype = "thinnode",
        mem = "4G",
    threads: 1
    shell:
        '''
        cat {input.R1} > {output.R1} && \
        cat {input.R2} > {output.R2}
        '''

