# MiMAG pipeline

A pipeline that runs `barrnap`, `bakta`, and `tRNAscan-ce` for a directory of bins.

`bakta` requires a database, which should be downloaded first and the path to the database should be provided in the [config]("./config/config.example.yaml") along with the directory with bins.

The summary file produced will contain the following:
The `custom_trna_uniq` is based on the tRNAs found by either bakta or tRNAscan-ce

```csv
bin	bakta_trna_all	bakta_trna_uniq	bakta_16s	bakta_23s	bakta_5s	bakta_16s_len	bakta_23s_len	bakta_5s_len	barrnap_16s	barrnap_23s	barrnap_5s	barrnap_16s_len	barrnap_23s_len	barrnap_5s_len	custom_trna_uniq
bin.2.21	30	20	1	1	1	1523	3012	116	1	1	1	1520	3006	93	19
bin.1.108	47	21	4	4	4	1533	2893	118	4	4	4	1529	2889	96	22
bin.1.24	37	18	0	1	2	0	1303	116	0	1	2	0	1298	110	18
```


