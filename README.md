

![](./scripts/celltag.final.gif)

# CellTag Workflow

This repository contains CellTag workflow, written in Python.

You may clone this repository and run the jupyter notebook file (`celltagging_workflow.ipynb`) to go through CellTag workflow.

**Paper**: The article can be found [here](https://www.nature.com/articles/s41586-018-0744-4).

**GEO dataset**: Here is the link to the GEO DataSet ([GSE99915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99915)) which contains all of the sequencing data we generated.

**CellTag Libraries**: The CellTag libraries are available from AddGene [here](https://www.addgene.org/pooled-library/morris-lab-celltag/). This link also contains whitelists for each CellTag library as well as sequence information for each library. The CellTag libraries available at AddGene are labeled CellTag-V1, CellTag-V2, and CellTag-V3. These labels correspond to CellTag<sup>MEF</sup>, CellTag<sup>D3</sup>, and CellTag<sup>D13</sup> respectively.

From our GEO DataSet we will select one timepoint to use as an example of our CellTag processing and analysis.
We will use data collected from Timecourse 1 at Day 15. The link to the SRA for this sample is [here](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR7347033).
We are using the data from Day 15 because we should be able to call clones using each of the CellTag versions (CellTag<sup>MEF</sup>, CellTag<sup>Day3</sup>, and CellTag<sup>Day13</sup>)

### Do the following steps in bash to get CellTag Matrix Expression
1. Download BAM file for samples at day 15
2. Extract celltag reads from BAM file
3. Parse reads to extract required information
I've also added the extracted file in the repository.


```python
# #bash
# wget https://sra-pub-src-1.s3.amazonaws.com/SRR7347033/hf1.d15.possorted_genome_bam.bam.1

# samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'GGT[ACTG]{8}GAATTC' > v1.celltag.reads.out
# samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'GTGATG[ACTG]{8}GAATTC' > v2.celltag.reads.out
# samtools view hf1.d15.possorted_genome_bam.bam | grep -P 'TGTACG[ACTG]{8}GAATTC' > v3.celltag.reads.out

# ./scripts/celltag.parse.reads.10x.sh -v tagregex="CCGGT([ACTG]{8})GAATTC" v1.celltag.reads.out > v1.celltag.parsed.tsv
# ./scripts/celltag.parse.reads.10x.sh -v tagregex="GTGATG([ACTG]{8})GAATTC" v2.celltag.reads.out > v2.celltag.parsed.tsv
# ./scripts/celltag.parse.reads.10x.sh -v tagregex="TGTACG([ACTG]{8})GAATTC" v3.celltag.reads.out > v3.celltag.parsed.tsv

# Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v1.celltag.parsed.tsv hf1.d15.v1
# Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v2.celltag.parsed.tsv hf1.d15.v2
# Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v3.celltag.parsed.tsv hf1.d15.v3
```


```python
import pyreadr #!pip install pyradr
import numpy as np
import pandas as pd
import celltagging_utils as ct
```

### Read data from celltag matrix


```python
mef = pyreadr.read_r("./celltag_matrix/hf1.d15.v1.celltag.matrix.Rds")[None]
d3 = pyreadr.read_r("./celltag_matrix/hf1.d15.v2.celltag.matrix.Rds")[None]
d13 = pyreadr.read_r("./celltag_matrix/hf1.d15.v3.celltag.matrix.Rds")[None]

mef.set_index('Cell.BC',inplace=True)
d3.set_index('Cell.BC',inplace=True)
d13.set_index('Cell.BC',inplace=True)

mef.shape, d3.shape, d13.shape
```




    ((3812, 6319), (3812, 8246), (3812, 4630))




```python
mef.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>AAAAAAGA</th>
      <th>AAAAAAGC</th>
      <th>AAAAAATA</th>
      <th>AAAAACTC</th>
      <th>AAAAACTG</th>
      <th>AAAAAGAC</th>
      <th>AAAAAGCC</th>
      <th>AAAAAGCG</th>
      <th>AAAAAGGG</th>
      <th>AAAACTAA</th>
      <th>...</th>
      <th>TTTTTCGG</th>
      <th>TTTTTCTA</th>
      <th>TTTTTGAT</th>
      <th>TTTTTGCA</th>
      <th>TTTTTGTT</th>
      <th>TTTTTTAT</th>
      <th>TTTTTTCC</th>
      <th>TTTTTTCT</th>
      <th>TTTTTTGG</th>
      <th>TTTTTTTT</th>
    </tr>
    <tr>
      <th>Cell.BC</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>AAACCTGAGTATGACA</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AAACCTGCAGCCTATA</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AAACCTGGTAAGTAGT</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AAACCTGTCACAACGT</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AAACCTGTCCGCGCAA</th>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 6319 columns</p>
</div>



### Get stats


```python
ct.get_stats(mef)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Cell_UMI_Counts</th>
      <th>CellTags_per_Cell</th>
      <th>CellTag_UMI_Counts</th>
      <th>Cells_per_CellTag</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>count</th>
      <td>3812.000000</td>
      <td>3812.000000</td>
      <td>6319.000000</td>
      <td>6319.000000</td>
    </tr>
    <tr>
      <th>mean</th>
      <td>54.328437</td>
      <td>5.917629</td>
      <td>32.774173</td>
      <td>3.569869</td>
    </tr>
    <tr>
      <th>std</th>
      <td>96.041066</td>
      <td>6.887981</td>
      <td>379.364072</td>
      <td>16.430824</td>
    </tr>
    <tr>
      <th>min</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>25%</th>
      <td>2.000000</td>
      <td>2.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>50%</th>
      <td>20.000000</td>
      <td>4.000000</td>
      <td>1.000000</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>75%</th>
      <td>64.000000</td>
      <td>8.000000</td>
      <td>4.000000</td>
      <td>2.000000</td>
    </tr>
    <tr>
      <th>max</th>
      <td>1549.000000</td>
      <td>84.000000</td>
      <td>20226.000000</td>
      <td>569.000000</td>
    </tr>
  </tbody>
</table>
</div>



### Binarize data
Convert the CellTag UMI count matrices into binary matrices and any Cell Barcode/CellTag pair with a UMI count less than a cutoff will be disregarded


```python
mef_bin = ct.single_cell_data_binarization(mef, 2)
d3_bin = ct.single_cell_data_binarization(d3, 2)
d13_bin = ct.single_cell_data_binarization(d13, 2)
```

### Filter Tags (White list)
‘Whitelisting’ is performed to remove PCR and sequencing artifacts that are not corrected in the previous step. This whitelisting consists of filtering out CellTags that are not detected from sequencing of the original complex CellTag library.


```python
mef_filt = ct.single_cell_data_whitelist(mef_bin, "./whitelist/V1.CellTag.Whitelist.csv")
d3_filt = ct.single_cell_data_whitelist(d3_bin, "./whitelist/V2.CellTag.Whitelist.csv")
d13_filt = ct.single_cell_data_whitelist(d13_bin, "./whitelist/V3.CellTag.Whitelist.csv")
print(mef_filt.shape, d3_filt.shape, d13_filt.shape)
```

    (3812, 3256) (3812, 2537) (3812, 1981)


### Filter Cells
Cells with >20 CellTags (likely to correspond to cell multiplets) and less than two unique CellTags per cell are filtered out.


```python
mef_filt = ct.metric_based_filtering(mef_filt, 20, "less")
d3_filt = ct.metric_based_filtering(d3_filt, 20, "less")
d13_filt = ct.metric_based_filtering(d13_filt, 20, "less")

mef_filt = ct.metric_based_filtering(mef_filt, 2, "greater")
d3_filt = ct.metric_based_filtering(d3_filt, 2, "greater")
d13_filt = ct.metric_based_filtering(d13_filt, 2, "greater")
print(mef_filt.shape, d3_filt.shape, d13_filt.shape)
```

    (1852, 3256) (1733, 2537) (717, 1981)


### Jaccard Analysis (Cell-Cell)
Clone calling is performed where Jaccard coefficient scores were calculated to assess the similarity of CellTag expression signatures in all cells in a pairwise manner, thereby identifying clonally related cells.


```python
mef_sim = ct.jaccard_analysis(mef_filt, id='mef')
d3_sim = ct.jaccard_analysis(d3_filt, id='d3')
d13_sim = ct.jaccard_analysis(d13_filt, id='d13')
```

    Calculating Jaccard Similarities:   0%|          | 0/1852 [00:00<?, ?it/s]
    Calculating Jaccard Similarities:   0%|          | 0/1733 [00:00<?, ?it/s]
    Calculating Jaccard Similarities:   0%|          | 0/717 [00:00<?, ?it/s]


### Clone Calling
Cells with a given threshold for similarity (here 0.7) will be get the same clone index. Celones with one cell will be disregarded.


```python
mef_clones, mef_clone_size = ct.clone_calling(mef_sim, "./hf1.d15.v1.clones.csv", 0.7)
d3_clones, d3_clone_size = ct.clone_calling(d3_sim, "./hf1.d15.v2.clones.csv", 0.7)
d13_clones, d13_clone_size = ct.clone_calling(d13_sim, "./hf1.d15.v13.clones.csv", 0.7)
print(mef_clones.head())
print(mef_clone_size.head())
```

       clone_id      cell_barcode
    0         1  CTAGCCTAGTGTCCAT
    1         1  AGCTTGAAGTACACCT
    2         1  CTGAAGTAGAGTTGGC
    3         1  CAGATCACAACTTGAC
    4         1  ATCCACCTCGAATGCT
       Clone_ID  Frequency
    0         2        228
    1         4        202
    2        16         65
    3        13         56
    4         8         53


### Lineage and Visualization


```python
mef_clones.rename(columns={mef_clones.columns[0]: "CellTagV1"}, inplace=True)
d3_clones.rename(columns={d3_clones.columns[0]: "CellTagV2"}, inplace=True)
d13_clones.rename(columns={d13_clones.columns[0]: "CellTagV3"}, inplace=True)
```


```python
clone_cells = pd.concat([mef_clones.cell_barcode, d3_clones.cell_barcode, d13_clones.cell_barcode]).unique()
celltag_data = pd.DataFrame(index=clone_cells, columns=["CellTagV1", "CellTagV2", "CellTagV3"])
```


```python
celltag_data.loc[mef_clones['cell_barcode'], "CellTagV1"] = mef_clones["CellTagV1"].values
celltag_data.loc[d3_clones['cell_barcode'], "CellTagV2"] = d3_clones["CellTagV2"].values
celltag_data.loc[d13_clones['cell_barcode'], "CellTagV3"] = d13_clones["CellTagV3"].values
celltag_data.index = celltag_data.index.map(lambda x: x + "-1")
celltag_data
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CellTagV1</th>
      <th>CellTagV2</th>
      <th>CellTagV3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CTAGCCTAGTGTCCAT-1</th>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>AGCTTGAAGTACACCT-1</th>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>CTGAAGTAGAGTTGGC-1</th>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>CAGATCACAACTTGAC-1</th>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>ATCCACCTCGAATGCT-1</th>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>GGAGCAAGTTACCGAT-1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>45</td>
    </tr>
    <tr>
      <th>TCTCATAGTCCGTTAA-1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>48</td>
    </tr>
    <tr>
      <th>TGACAACCAAAGGAAG-1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>48</td>
    </tr>
    <tr>
      <th>TGCGCAGAGATATGCA-1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>50</td>
    </tr>
    <tr>
      <th>TTAGGACGTAAGTTCC-1</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>50</td>
    </tr>
  </tbody>
</table>
<p>1957 rows × 3 columns</p>
</div>



### Create Linklist of Relation between Cells and CellTags


```python
link_list = ct.convert_cell_tag_matrix_to_link_list(celltag_data)
nodes = ct.get_nodes_from_link_list(link_list)
```

    Preprocessing data..
    Cells that have CellTagV1: 1402
    Cells that have CellTagV2: 1174
    Cells that have CellTagV3: 100
    find connection between [celltag -> cells]...
    find hidden links [CellTagV2 -> CellTagV3], or [CellTagV1 -> CellTagV3]...
    find hidden links [CellTagV1 -> CellTagV2]...
    finished



```python
link_list
```



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source</th>
      <th>target</th>
      <th>tag</th>
      <th>target_unmodified</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CellTagV3_30</td>
      <td>GCGGGTTTCATGCTCC-1_V3</td>
      <td>CellTagV3</td>
      <td>GCGGGTTTCATGCTCC-1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CellTagV3_46</td>
      <td>GGTGAAGCAGTCACTA-1_V3</td>
      <td>CellTagV3</td>
      <td>GGTGAAGCAGTCACTA-1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>CellTagV3_12</td>
      <td>ATCATCTCAATGGAGC-1_V3</td>
      <td>CellTagV3</td>
      <td>ATCATCTCAATGGAGC-1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>CellTagV3_3</td>
      <td>AGGCCACAGCGGCTTC-1_V3</td>
      <td>CellTagV3</td>
      <td>AGGCCACAGCGGCTTC-1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>CellTagV3_30</td>
      <td>CGAGCACAGCACCGTC-1_V3</td>
      <td>CellTagV3</td>
      <td>CGAGCACAGCACCGTC-1</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2802</th>
      <td>CellTagV2_176</td>
      <td>CellTagV3_39</td>
      <td>CellTagV2</td>
      <td>CellTagV3_39</td>
    </tr>
    <tr>
      <th>2803</th>
      <td>CellTagV2_17</td>
      <td>CellTagV3_2</td>
      <td>CellTagV2</td>
      <td>CellTagV3_2</td>
    </tr>
    <tr>
      <th>2804</th>
      <td>CellTagV2_100</td>
      <td>CellTagV3_14</td>
      <td>CellTagV2</td>
      <td>CellTagV3_14</td>
    </tr>
    <tr>
      <th>2805</th>
      <td>CellTagV2_114</td>
      <td>CellTagV3_17</td>
      <td>CellTagV2</td>
      <td>CellTagV3_17</td>
    </tr>
    <tr>
      <th>2806</th>
      <td>CellTagV2_148</td>
      <td>CellTagV3_28</td>
      <td>CellTagV2</td>
      <td>CellTagV3_28</td>
    </tr>
  </tbody>
</table>
<p>2807 rows × 4 columns</p>
</div>



