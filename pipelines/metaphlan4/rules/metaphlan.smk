# DB used: mpa_vJan21_CHOCOPhlAn_202103
rule MetaPhlAn:
    conda:
        "../_envs/metaphlan.yaml"
    input:
        R1=OUT_DIR + "{barcode}/{barcode}_R1_NonHuman_Combined_Subsampled_" + str(CUTOFF) + ".fastq.gz",
        R2=OUT_DIR + "{barcode}/{barcode}_R2_NonHuman_Combined_Subsampled_" + str(CUTOFF) + ".fastq.gz",
    params:
        index=METAPHLAN_INDEX,
        db=METAPHLAN_DB,
    threads: 8
    resources:
        walltime = "86000",
        mem = "34GB",
        nodetype = "thinnode",
    output:
        bowtieindex=OUT_DIR + "{barcode}/{barcode}_NonHuman_Combined_Subsampled_" + str(CUTOFF) + "_bowtie2.bz2",
        txt=OUT_DIR + "{barcode}/{barcode}_NonHuman_Combined_Subsampled_" + str(CUTOFF) + "_profile.txt",
    log:
        MetaPhlAn=OUT_DIR + "{barcode}/log/{barcode}_bowtie2_" + getCurrentTime()+ ".log",
    shell:
        """
        metaphlan {input.R1},{input.R2} \
        --input_type fastq \
        --nproc {threads} \
        --bowtie2db {params.db} \
        --index {params.index} \
        --bowtie2out {output.bowtieindex} \
        --unclassified_estimation \
        -o {output.txt} #&> {log.MetaPhlAn}
        """

# envmodules:
#         "tools",
#         "bowtie2/2.5.0"
#         # 'metaphlan/4.0.3',


rule MergeMetaPhlAnProfiles:
    conda:
        "../_envs/metaphlan.yaml"
    input:
        files=expand(
            OUT_DIR + "{barcode}/{barcode}_NonHuman_Combined_Subsampled_" + str(CUTOFF) + "_profile.txt",
            barcode=SUBSAMPLED_BARCODES,
        ),
    threads: 2,
    resources:
        walltime = "86000",
        mem = "20GB",
        nodetype = "thinnode",
    output:
        txt=OUT_DIR + "MetaPhlAn_4.1.0_Combined_NonHuman_Subsampled_" + str(CUTOFF) + "_profile.txt",
    # envmodules:
    #     "tools",
    #     "metaphlan/4.0.3",
    shell:
        """
        merge_metaphlan_tables.py \
        {input.files} > {output.txt}
        """
