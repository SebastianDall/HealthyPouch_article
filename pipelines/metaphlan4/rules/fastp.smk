rule FastpFilter:
    input:
        R1 = get_fastq_R1,
        R2 = get_fastq_R2
    threads:
        8
    resources:
        walltime = "7200",
        mem = "40G",
        nodetype = "thinnode"
    params:
        correction = "--correction",
        detect_adapter_for_pe = "--detect_adapter_for_pe",
        dedup = "--dedup",
        cut_right = "--cut_right",
        overrepresentation_analysis = "--overrepresentation_analysis",
        cut_right_window_size = 4,
        cut_right_mean_quality = 20,
        average_qual = 25,
        length_required = 100,
        dup_calc_accuracy = 6
    output:
        R1     = OUT_DIR+"{barcode}/fastp/{sample}_R1_fastp.fastq",
        R2     = OUT_DIR+"{barcode}/fastp/{sample}_R2_fastp.fastq",
        json   = OUT_DIR+"{barcode}/fastp/qc/{sample}_fastp.json",
        html   = OUT_DIR+"{barcode}/fastp/qc/{sample}_fastp.html"
    log: 
        fastp  = OUT_DIR+"{barcode}/log/{sample}_fastp_"+getCurrentTime()+".log",
    envmodules:
        "tools",
        "fastp/0.23.2"
    shell:
        '''
        fastp \
        --in1 {input.R1} \
        --in2 {input.R2} \
        --out1 {output.R1} \
        --out2 {output.R2}  \
        {params.correction} \
        {params.detect_adapter_for_pe} \
        {params.dedup} \
        {params.cut_right} \
        {params.overrepresentation_analysis} \
        --cut_right_window_size {params.cut_right_window_size} \
        --cut_right_mean_quality {params.cut_right_mean_quality} \
        --average_qual {params.average_qual} \
        --length_required {params.length_required} \
        --dup_calc_accuracy {params.dup_calc_accuracy} \
        --thread {threads} \
        --json {output.json} \
        --html {output.html} &> {log.fastp}
        '''
