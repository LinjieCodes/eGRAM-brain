# Human-specific lncRNAs contributed critically to human evolution by distinctly regulating gene expression

The eGRAM program identifies gene modules comprising co-expressed genes and their regulatory lncRNAs based on lncRNA/DNA bindings and gene expression correlations.

There are two kinds of input files - (a) user data and (b) system data. The former has two files - (a) expression profile, (b) lncRNA target prediction. LncRNAs' targets are predicted using LongTarget (or the more rapid version Fasim). System data include KEGG pathway annotation files downloaded from the KEGG websites.

# Requirements
1. **Python**: >=3.7.0

2. **numpy**: >=1.21.6

3. **pandas**: >=1.3.5

4. **scipy**: >=1.7.3

5. **OS**: the eGRAMv2R1 code has been tested on Linux system.

# Data
1. **GTEx_Frontal_cortex_BA9**  --  the gene expression matrix of human frontal cortex.

2. **GTEx_Anterior_cingulate_cortex_BA24**  --  the gene expression matrix of human anterior cingulate cortex.

3. **macaque_BA9**  --  the gene expression matrix of macaque frontal cortex.

4. **macaque_BA24**  --  the gene expression matrix of macaque anterior cingulate cortex.

5. **TableS4-HS_lncRNA_DBS_matrix**  --  the DNA binding matrix (i.e. target genes) of lncRNAs.

6. **PathwayAnnotation**  --  the KEGG pathway annotation.


# Usage
Here is a command line to run the eGRAM program:

```
'Example: python eGRAM.py --t1 100 --t2 60 --m 20 --c 0.5 --f1 TableS4-HS_lncRNA_DBS_matrix --f2 data/GTEx_Frontal_cortex_BA9 --f3 Kegg_NervousSystem --s human --o BA9_module'
```

# Help information
Here is a brief explanation of the command line arguments:

```
Options   Parameters      Functions
t1   Binding affinity threshold  An integer specifying the threshold for determining a strong binding.
t2   Binding affinity threshold  An integer specifying the threshold for determining a qualified binding.
m    Module size                 An integer specifying the minimum size of modules.
c    Correlation threshold       A floating point, indicating the minimum correlation coefficient to determine co-expressed genes.
f1   Binding affinity matrix     A string indicating the binding affinity matrix file.
f2   Gene expression matrix      A string indicating the gene expression matrix file.
o    Output                      A string specifying the output file name.
```

# Bug reports
Please send comments and bug reports to JL.linjie@outlook.com.

# Related website
To obtain details about lncRNA/DNA binding prediction using LongTarget, please go to our website http://www.gaemons.net/LongTarget.
