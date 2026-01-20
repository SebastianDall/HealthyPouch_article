import os
import glob
import datetime

SNAKEDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(SNAKEDIR, "..", "config", "config.yaml")

rule all:
  input:
    expand(os.path.join(config["outdir"], "assembly_stats", "{sample}_stats.csv"), sample = config["samples"]),
    os.path.join(config["outdir"], "assembly_stats", "stats_combined.csv")


rule assembly_stats:
  conda: "python311",
  input:
    i = lambda wc: config["samples"][wc.sample]["assembly"],
    b = lambda wc: config["samples"][wc.sample]["contig_bin"],
  output:
    os.path.join(config["outdir"], "assembly_stats", "{sample}_stats.csv")
  threads: 1,
  resources:
      walltime = "1:00:00:00",
      mem = "10G",
      nodetype = "thinnode",
  params:
    script = os.path.join(SNAKEDIR, "scripts", "assembly_stats.py")
  shell:
    """
    python {params.script} -i {input.i} -b {input.b} -o {output}
    """

rule gather_stats:
  input:
    stats = expand(os.path.join(config["outdir"], "assembly_stats", "{sample}_stats.csv"), sample = config["samples"]),
  output:
    os.path.join(config["outdir"], "assembly_stats", "stats_combined.csv"),
  threads: 1,
  resources:
      walltime = "1:00:00:00",
      mem = "10G",
      nodetype = "thinnode",
  run:
    import pandas as pd
    df_list = [pd.read_csv(path) for path in input.stats]
    # concatenate into a single DataFrame
    summary_df = pd.concat(df_list, ignore_index=True)
    # write out the merged summary
    summary_df.to_csv(output[0], index=False)
