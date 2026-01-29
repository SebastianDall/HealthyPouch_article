import glob
import os
import pandas as pd


SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config.yaml")

path = config["bins_path"]
mags = os.listdir(path)
base_mags = [m.rsplit(".fa")[0] for m in mags]

rule all:
  input:
    os.path.join("bakta", "annotation_summary.tsv"),
    # os.path.join("bakta", "mefinder.tsv")

# rule mefinder:
#   conda: "mobile_element_finder_env"
#   input:
#     os.path.join(path, "{mag}.fa")
#   output:
#     os.path.join(config["outdir"], "{mag}", "mobile_element_finder", "{mag}_mef.csv.csv.csv")
#   threads: 4
#   resources:
#     walltime = "1:00:00",
#     mem = "15G",
#   params:
#     db = config["mgedb"],
#     o = os.path.join(config["outdir"], "{mag}", "mobile_element_finder", "{mag}_mef.csv.csv")
#   shell:
#     """
#     mefinder find -c {input} --db-path {params.db} --threads {threads} {params.o}
#     """

# rule concat_mefinder:
#   input:
#     files = expand(os.path.join(config["outdir"], "{mag}", "mobile_element_finder", "{mag}_mef.csv.csv.csv"), mag = base_mags)
#   output:
#     o = os.path.join(config["outdir"], "mefinder.tsv")
#   threads: 1
#   resources:
#     walltime = "1:00:00",
#     mem = "15G",
#   run:
#     dfs = []

#     for f in input.files:
#       b = pd.read_csv(f, skiprows=4)
#       mag_name = os.path.basename(f).replace("_mef.csv.csv")
#       b["bin"] = mag_name
#       dfs.append(b)

#     df = pd.concat(dfs)
#     df.to_csv(output.o, index = False)



rule bakta:
  conda: "bakta"
  input:
    os.path.join(path, "{mag}.fa")
  output:
    os.path.join(config["outdir"], "{mag}", "{mag}.tsv"),
  threads: 4
  resources:
    walltime = "1:00:00",
    mem = "15G",
  params:
    db = config["baktadb"],
    out = os.path.join(config["outdir"], "{mag}"),
  shell:
    """
      bakta --db {params.db} --threads {threads} -o {params.out} --force \
      --skip-crispr \
      --skip-pseudo \
      --skip-gap \
      --skip-ori \
      --skip-plot \
      --force \
      --keep-contig-headers \
      --meta
      {input}
    """

rule barrnap:
  conda: "barnap"
  input:
    os.path.join(path, "{mag}.fa")
  output:
    os.path.join(config["outdir"], "{mag}", "{mag}_barrnap.tsv"),
  threads: 4
  resources:
    walltime = "1:00:00",
    mem = "15G",
  shell:
    """
    BIN={wildcards.mag}
    echo "bin,barrnap_16s,barrnap_23s,barrnap_5s,barrnap_16s_len,barrnap_23s_len,barrnap_5s_len" > {output}

    barrnap --threads {threads} --kingdom bac --quiet {input} | \
    awk -F "\t" -v bin=$BIN '
    BEGIN {{
        a["16S"]=0; a["23S"]=0; a["5S"]=0
        len["16S"]=0; len["23S"]=0; len["5S"]=0
      }}
      /^[^#]/ {{
        split($9, attrs, ";")
        for (i in attrs) {{
          if (attrs[i] ~ /^Name=/) {{
            sub(/Name=/, "", attrs[i])
            sub(/_rRNA/, "", attrs[i])
            gene = attrs[i]
            break
          }}
        }}
        gene_len = $5 - $4 + 1
        a[gene]++
        if (len[gene] < gene_len) len[gene] = gene_len
      }}
      END {{
        OFS=","
        print bin, a["16S"], a["23S"], a["5S"], len["16S"], len["23S"], len["5S"]
      }}
      ' >> {output}
    """

rule trnascan:
  envmodules:
    "tools",
    "trnascan-se/2.0.12"
  input:
    os.path.join(path, "{mag}.fa")
  output:
    out = os.path.join(config["outdir"], "{mag}", "{mag}_trna.tsv"),
    temp_out = temp(os.path.join(config["outdir"], "{mag}", "trna_{mag}.txt")),
    temp_stats = temp(os.path.join(config["outdir"], "{mag}", "stats_{mag}.txt")),
  threads: 4
  resources:
    walltime = "1:00:00",
    mem = "15G",
  params:
    kingdom = "A",
  shell:
    """
     tRNAscan-SE -{params.kingdom} \
      -o {output.temp_out} \
      -m {output.temp_stats} \
      -d {input} \
      --thread {threads}
    
    # Remove first 3 header lines and add MAG name as last column
    sed -e '1,3d' -e "s/$/\t{wildcards.mag}/g" {output.temp_out} > {output.out}
    """


rule barrnap_sum:
  input:
    files = expand(os.path.join(config["outdir"], "{mag}", "{mag}_barrnap.tsv"), mag = base_mags)
  output:
    o = os.path.join(config["outdir"], "barrnap_stats.csv")
  resources:
    walltime = "1:00:00",
    mem = "5G",
  run:
    dfs = []

    for f in input.files:
      b = pd.read_csv(f)
      dfs.append(b)

    df = pd.concat(dfs)
    df.to_csv(output.o, index = False)
    

rule trnascan_sum:
  input:
    expand(os.path.join(config["outdir"], "{mag}", "{mag}_trna.tsv"), mag = base_mags)
  output:
    os.path.join(config["outdir"], "trna_stats.csv")
  resources:
    walltime = "1:00:00",
    mem = "5G",
  run:
    import pandas as pd
    
    # Read all tRNA files and concatenate
    dfs = []
    for file in input:
        df = pd.read_csv(file, sep='\t', header=None, comment='#')
        dfs.append(df)
    
    all_trna = pd.concat(dfs, ignore_index=True)
    
    # Extract tRNA type (column 4, 0-indexed) and bin name (column 10)
    trna_data = all_trna[[4, 10]].copy()
    trna_data.columns = ['trna_type', 'bin']
    
    # Remove undetermined and suppressed tRNAs
    trna_data = trna_data[~trna_data['trna_type'].str.contains('Undet|Sup', na=False)]
    
    # Count unique tRNA types per bin
    unique_counts = trna_data.drop_duplicates().groupby('bin').size().reset_index()
    unique_counts.columns = ['bin', 'custom_trna_uniq']
    
    # Save to CSV
    unique_counts.to_csv(output[0], index=False)

rule bakta_sum:
  input:
    expand(os.path.join(config["outdir"], "{mag}", "{mag}.tsv"), mag = base_mags)
  output:
    os.path.join(config["outdir"], "bakta_stats.csv")
  resources:
    walltime = "1:00:00",
    mem = "5G",
  shell:
    """
    echo "bin,bakta_trna_all,bakta_trna_uniq,bakta_16s,bakta_23s,bakta_5s,bakta_16s_len,bakta_23s_len,bakta_5s_len" > {output}
    
    for tsv_file in {input}; do
        dir=$(dirname $tsv_file)
        name=$(basename $dir)

        tRNA_all=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $tsv_file | grep -c "trn" -) || true
        tRNA_uniq=$(awk -F "\t" '{{ if ($2 == "tRNA") {{print $7}} }}' $tsv_file | sort -u - | grep -c "trn" -) || true
        rRNA_16S=$(awk -F "\t" '{{ if ($7 == "rrs") {{print $2}} }}' $tsv_file | grep -c "rRNA" -) || true
        rRNA_23S=$(awk -F "\t" '{{ if ($7 == "rrl") {{print $2}} }}' $tsv_file | grep -c "rRNA" -) || true
        rRNA_5S=$(awk -F "\t" '{{ if ($7 == "rrf") {{print $2}} }}' $tsv_file | grep -c "rRNA" -) || true

        rRNA_16S_len=$(awk -F "\t" '{{ if ($7 == "rrs" && $2 == "rRNA") {{ len = $4 - $3 + 1; if (len > max) max = len }} }} END {{ print max+0 }}' $tsv_file)
        rRNA_23S_len=$(awk -F "\t" '{{ if ($7 == "rrl" && $2 == "rRNA") {{ len = $4 - $3 + 1; if (len > max) max = len }} }} END {{ print max+0 }}' $tsv_file)
        rRNA_5S_len=$(awk -F "\t" '{{ if ($7 == "rrf" && $2 == "rRNA") {{ len = $4 - $3 + 1; if (len > max) max = len }} }} END {{ print max+0 }}' $tsv_file)
        echo "$name,$tRNA_all,$tRNA_uniq,$rRNA_16S,$rRNA_23S,$rRNA_5S,$rRNA_16S_len,$rRNA_23S_len,$rRNA_5S_len" >> {output}
    done
    """

rule annotation_aggregate:
  input:
    bakta=os.path.join(config["outdir"], "bakta_stats.csv"),
    barrnap=os.path.join(config["outdir"], "barrnap_stats.csv"),
    trnascan=os.path.join(config["outdir"], "trna_stats.csv")
  output:
    o = os.path.join(config["outdir"], "annotation_summary.tsv")
  resources:
    walltime = "1:00:00",
    mem = "5G",
  run:
    bakta = pd.read_csv(input.bakta)
    barrnap = pd.read_csv(input.barrnap)
    trnascan = pd.read_csv(input.trnascan)

    # Merge all three dataframes on 'bin' column
    merged = bakta.merge(barrnap, on='bin', how='outer')
    merged = merged.merge(trnascan, on='bin', how='outer')

    # Save to TSV
    merged.to_csv(output.o, sep='\t', index=False)    
