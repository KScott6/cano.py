# cano.py walkthrough

Inspired by [BUSCO_phylogenomics](https://github.com/jamiemcg/BUSCO_phylogenomics), by Jamie McGowan. Cano.py also creates a species tree from BUSCO output, but instead uses a [partitioned analysis](https://iqtree.github.io/doc/Advanced-Tutorial), where IQ-TREE is used to identify the best substition model for every BUSCO gene used. 

I also recommend checking out [MycoTools](https://github.com/xonq/mycotools) - it is great for downloading and organizing your genome data prior to using cano.py. Cano.py works very well with the "ome"-based naming system that MycoTools uses.

Canopy is currently coded to work seamlessly on SCINet for users in the "arsef" project. If you are on a different system or do not have access to the arsef project, it will take a bit more work to use cano.py (still possible with miminal effort).

What cano.py does:  takes BUSCO output for a set of genomes, finds all Single Copy Orthologs (SCOs), creates individual gene trees in IQ-TREE with the modelfinder option, then runs a partitioned analysis in IQ-TREE - the final result being a species tree with bootstrap values.

**cano.py input:**  BUSCO output of for a set of genomes

**cano.py output:**  IQ-TREE species tree (partitioned analysis), individual IQ-TREE gene tree, various useful metadata files

<br>

---

<br>

# Assumptions

This script assumes you have a collection of genomes and have already completed a BUSCO analysis for each genome (using the same BUSCO database for each genome). Here is an example BUSCO command:

```bash
# SCINet users can use the shared environment "funannotate"
module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate
busco -i acitak1.fna -o "acitak1" -m genome -l "/project/arsef/databases/funannotate_databases/fungi_odb10"
```

Cano.py also assumes you have stored this BUSCO output in a folder named the same as the "ome" code ("acitak1" in the given example). All BUSCO output folders should be stored in the same folder.

SCINet users:  I have run BUSCO ("fungi" database) on the entirety of the MycoTools database. The BUSCO output folder can be found here: /project/arsef/databases/mycotools/database_stats/busco/fungi/by_ome/

The last assumption is that you have access to the arsef shared conda environment "canopy". If you don't, just make sure that you have iqtree (v3.0.1), trimal (v1.5.rev0), and muscle (v5.1) in your PATH.

<br>

# Running cano.py

Provide a list of the "omes" you want to include in the analysis, the path to your BUSCO outputs, the number of threads you want to use, the number of bootstraps desired for the gene tree and species tree creation, and a path to your desired output. I recommend running as a job, as this may take several days (test case:  took 1.5 days for 130 genome with 64 threads).

```bash
module load miniconda
source activate /project/arsef/environments/canopy

python /project/arsef/projects/hypocreales_tree/cano.py \
--omes /project/arsef/projects/hypocreales_tree/testset1/testset1_ome_list.txt \
--busco /project/arsef/databases/mycotools/database_stats/busco/fungi/by_ome/ \
--threads 64 \
--bootstraps 1000 \
--sco-threshold 1.0 \
--output /project/arsef/projects/hypocreales_tree/testset1/canopy/
```

I want to point out the option `sco-threshold`, which is the minimum proportion of genomes with a particular BUSCO to retain it as single-copy (default: 1.0). 

Sometimes, if you have poor quality genomes, a very diverse set of genomes, or a very large set of genomes, you will be unable to find any true SCOs in your dataset. Lower the sco-threshold until you have enough BUSCO genes that make the cut. 

So if you have 100 genomes and set the sco-threshold to 0.9 (90%), only 90 genomes will need to have a particular BUSCO gene complete and single-copy for the BUSCO gene to be treated as an SCO.

<br>

# Explanation of output

`busco_genes_presence_absence_matrix_FULL.tsv` : a matrix indicating the presence/absence of BUSCO genes (single copy and complete) for each genome, for every single BUSCO gene in a given database (e.g. if you provide the fungi odb10 BUSCO output, this will always include 758 fugnal BUSCO genes)

`busco_genes_presence_absence_matrix.tsv` : a matrix indicating the presence/absence of BUSCO genes (single copy and complete) for each genome. 

`single_copy_sequences` : folder containing aa fastas of BUSCO genes, one fasta for each determined SCO/pseudo-SCO

`alignments` : MUSCLE alignments of the fastas from single_copy_sequences

`trimmed_alignments` : alignments trimmed with trimAl

`single_gene_trees` : output of IQ-TREE for each individual aligned+trimmed BUSCO SCO

`concatenated_alignment` : a concatenated alignment fasta prepared for the IQ-TREE partitioned analysis. Also contains the partition file (.nexus) and partitioned coordinates file (.txt)

`final_tree` : output of the IQ-TREE partitioned analysis. The species tree is the .treefile

`final_report.txt` : a useful text file recording parameters used, number of BUSCOs identified as SCOs, runtime, etc. 

`logs` :  some logfiles that are useful in case you have an error (missing BUSCO output, SCO-threshold set too high)

<br>

---

<br>

# Citations

If you used cano.py, please cite all the software wrapped into this script. If you want, you can also cite this webpage for cano.py itself. I have listed the accesspory software versions used in the shared conda environment on SCINet. 

**IQ-TREE (v3.0.1)** : Wong, T.K., Ly-Trong, N., Ren, H., Baños, H., Roger, A.J., Susko, E., Bielow, C., De Maio, N., Goldman, N., Hahn, M.W. and Huttley, G., 2025. IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models.

**trimal (v1.5.rev0 build[2024-05-27])** : Capella-Gutiérrez, S., Silla-Martínez, J.M. and Gabaldón, T., 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), pp.1972-1973.

**Muscle (5.1.linux64)** :  Edgar, R.C., 2022. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. Nature communications, 13(1), p.6968.

---