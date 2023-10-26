

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
Cells with a given threshold for similarity (here 0.7) will be get the same clone index. Clones with one cell will be disregarded.


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
- CellTagV1 Clones: 166 (1402 Cells)
- CellTagV2 Clones: 197 (1174 Cells)
- CellTagV3 Clones: 50 (100 Cells)
- Cells in 1 Clone: 1268
- Cells in 2 Clones: 659 (1318 Nodes)
- Cells in 3 Clones: 30 (90 Nodes)
- Unique Cells: 1957
- **Nodes Count: 3089** -> 413 Clones + 2676 Cells -> (166+197+50) + (1268+1318+90)

(1268 + 1318 + 90) = (1402 + 1174 + 100) = 2676


```python
link_list = ct.convert_cell_tag_matrix_to_link_list(celltag_data)
all_nodes = ct.get_nodes_from_link_list(link_list)
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
ref_nodes = ["CellTagV1_95"]
sub_links, sub_nodes = ct.get_subnet(ref_nodes, link_list, all_nodes)
print("Links:", sub_links.shape[0], "\nNodes:",sub_nodes.shape[0])
```

    Links: 25 
    Nodes: 15



```python
G = nx.Graph()
G.add_nodes_from(sub_nodes['nodes'])
edges = [(row['source'], row['target']) for index, row in sub_links.iterrows()]
G.add_edges_from(edges)

# Define the layout of the graph
layout = go.Layout(
    title="Network",
    showlegend=False,
    hovermode='closest',
    width=800,
    height=600
)

plotly_graph = nx.spring_layout(G)  

# Create node positions
node_x = []
node_y = []
for node in plotly_graph:
    x, y = plotly_graph[node]
    node_x.append(x)
    node_y.append(y)

color_mapping = {
    'cell': 'red',
    'clone': 'blue',
}
node_colors = [color_mapping.get(category, 'gray') for category in sub_nodes['type']]

node_trace = go.Scatter(
    x=node_x,
    y=node_y,
    mode='markers+text',  # Add 'text' mode to display labels
    #"hoverinfo='text',
    text=sub_nodes['nodes'],  # Set the labels to node names
    textposition='bottom center',  # Adjust the label position
    marker=dict(
        showscale=False,
        color=node_colors,
        size=10
    )
)

edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = plotly_graph[edge[0]]
    x1, y1 = plotly_graph[edge[1]]
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x,
    y=edge_y,
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines'
)

fig = go.Figure(data=[edge_trace, node_trace], layout=layout)

iplot(fig)
```


<div>                            <div id="c838f629-596e-40c7-8edd-24698cece92c" class="plotly-graph-div" style="height:600px; width:800px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("c838f629-596e-40c7-8edd-24698cece92c")) {                    Plotly.newPlot(                        "c838f629-596e-40c7-8edd-24698cece92c",                        [{"hoverinfo":"none","line":{"color":"#888","width":0.5},"mode":"lines","x":[0.09321163696583713,0.0926191676389283,null,0.2584437430596077,0.0926191676389283,null,-0.09318372261472853,0.0926191676389283,null,-0.1787929718453414,0.0926191676389283,null,0.45462745516849656,0.0926191676389283,null,0.3741666434966264,0.0926191676389283,null,0.423929772012322,0.0926191676389283,null,0.29429248396269336,0.0926191676389283,null,-0.007015646491360902,0.0926191676389283,null,0.1964132867843469,0.0926191676389283,null,0.01662451923388329,0.0926191676389283,null,-0.17287047362115618,0.0926191676389283,null,-0.646728773502388,-1.0,null,-0.646728773502388,-0.860780907271751,null,-0.646728773502388,-0.9611186291800916,null,-0.646728773502388,-0.7476951303143334,null,-0.646728773502388,-0.8207710975427938,null,-0.646728773502388,0.0926191676389283,null,-0.646728773502388,-0.9982347672912837,null,0.5644949379757486,0.480020988225612,null,0.5644949379757486,0.896570315384758,null,0.5644949379757486,0.6515249508887048,null,0.5644949379757486,0.8960244266707574,null,0.5644949379757486,0.7942277922069045,null,0.5644949379757486,0.0926191676389283,null],"y":[0.15065775649522983,-0.15672078388231322,null,-0.3495075967625502,-0.15672078388231322,null,0.08119076913388411,-0.15672078388231322,null,-0.1499880078133015,-0.15672078388231322,null,-0.2874049980761336,-0.15672078388231322,null,-0.4662213872212358,-0.15672078388231322,null,-0.11482433670348534,-0.15672078388231322,null,0.040965968446403926,-0.15672078388231322,null,-0.3934538385088683,-0.15672078388231322,null,-0.5434733105993319,-0.15672078388231322,null,-0.5537612702699515,-0.15672078388231322,null,-0.4038367396026915,-0.15672078388231322,null,-0.05702350750214658,0.09851604842706728,null,-0.05702350750214658,0.1543149627109535,null,-0.05702350750214658,-0.23250230389938603,null,-0.05702350750214658,0.2656552374383664,null,-0.05702350750214658,-0.3352148703337683,null,-0.05702350750214658,-0.15672078388231322,null,-0.05702350750214658,-0.0738506545698931,null,0.38601137101231014,0.6988853942271266,null,0.38601137101231014,0.5180025639893666,null,0.38601137101231014,0.7340167677203039,null,0.38601137101231014,0.33887246674759214,null,0.38601137101231014,0.650694299396453,null,0.38601137101231014,-0.15672078388231322,null],"type":"scatter"},{"marker":{"color":["red","red","red","red","red","red","red","red","red","red","red","red","blue","blue","blue"],"showscale":false,"size":10},"mode":"markers+text","text":["CCGTTCAGTTTGTTTC-1_V1","CTCCTAGTCATATCGG-1_V1","CGGACACCACCGAATT-1_V1","CTACACCGTAGGGACT-1_V1","CGTAGGCAGCGATAGC-1_V1","GAAACTCTCAATACCG-1_V1","GCATGATGTACCCAAT-1_V1","TTGCCGTAGTACACCT-1_V1","ATCGAGTTCCCTAATT-1_V1","CGACCTTTCTACCTGC-1_V1","TTGAACGCATCAGTAC-1_V1","TAGACCACAATGACCT-1_V1","CellTagV2_147","CellTagV2_109","CellTagV1_95"],"textposition":"bottom center","x":[0.09321163696583713,0.2584437430596077,-0.09318372261472853,-0.1787929718453414,0.45462745516849656,0.3741666434966264,0.423929772012322,0.29429248396269336,-0.007015646491360902,0.1964132867843469,0.01662451923388329,-0.17287047362115618,-0.646728773502388,0.5644949379757486,0.0926191676389283,-1.0,0.480020988225612,0.896570315384758,-0.860780907271751,0.6515249508887048,0.8960244266707574,0.7942277922069045,-0.9611186291800916,-0.7476951303143334,-0.8207710975427938,-0.9982347672912837],"y":[0.15065775649522983,-0.3495075967625502,0.08119076913388411,-0.1499880078133015,-0.2874049980761336,-0.4662213872212358,-0.11482433670348534,0.040965968446403926,-0.3934538385088683,-0.5434733105993319,-0.5537612702699515,-0.4038367396026915,-0.05702350750214658,0.38601137101231014,-0.15672078388231322,0.09851604842706728,0.6988853942271266,0.5180025639893666,0.1543149627109535,0.7340167677203039,0.33887246674759214,0.650694299396453,-0.23250230389938603,0.2656552374383664,-0.3352148703337683,-0.0738506545698931],"type":"scatter"}],                        {"height":600,"hovermode":"closest","showlegend":false,"template":{"data":{"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"#E5ECF6","showlakes":true,"showland":true,"subunitcolor":"white"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2}}},"title":{"text":"Network"},"width":800},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('c838f629-596e-40c7-8edd-24698cece92c');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>



```python
celltag_data[celltag_data['CellTagV1']==95]
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
      <th>CCGTTCAGTTTGTTTC-1</th>
      <td>95</td>
      <td>147</td>
      <td>26</td>
    </tr>
    <tr>
      <th>CTCCTAGTCATATCGG-1</th>
      <td>95</td>
      <td>e</td>
      <td>e</td>
    </tr>
    <tr>
      <th>CGGACACCACCGAATT-1</th>
      <td>95</td>
      <td>109</td>
      <td>e</td>
    </tr>
    <tr>
      <th>CTACACCGTAGGGACT-1</th>
      <td>95</td>
      <td>109</td>
      <td>e</td>
    </tr>
    <tr>
      <th>CGTAGGCAGCGATAGC-1</th>
      <td>95</td>
      <td>147</td>
      <td>e</td>
    </tr>
    <tr>
      <th>GAAACTCTCAATACCG-1</th>
      <td>95</td>
      <td>109</td>
      <td>e</td>
    </tr>
    <tr>
      <th>GCATGATGTACCCAAT-1</th>
      <td>95</td>
      <td>e</td>
      <td>e</td>
    </tr>
    <tr>
      <th>TTGCCGTAGTACACCT-1</th>
      <td>95</td>
      <td>109</td>
      <td>e</td>
    </tr>
    <tr>
      <th>ATCGAGTTCCCTAATT-1</th>
      <td>95</td>
      <td>109</td>
      <td>e</td>
    </tr>
    <tr>
      <th>CGACCTTTCTACCTGC-1</th>
      <td>95</td>
      <td>147</td>
      <td>26</td>
    </tr>
    <tr>
      <th>TTGAACGCATCAGTAC-1</th>
      <td>95</td>
      <td>147</td>
      <td>e</td>
    </tr>
    <tr>
      <th>TAGACCACAATGACCT-1</th>
      <td>95</td>
      <td>147</td>
      <td>e</td>
    </tr>
  </tbody>
</table>
</div>


