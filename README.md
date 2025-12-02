<p align="center">
  <img src="extras/canopy_logo.png" alt="cano.py logo" width="400">
</p>

<br>

# cano.py walkthrough

Takes OrthoFinder output **or** BUSCO output and automatically creates phylogenomic trees. 

The BUSCO method was inspired by [BUSCO_phylogenomics](https://github.com/jamiemcg/BUSCO_phylogenomics), by Jamie McGowan. Cano.py also creates a species tree from BUSCO output, but instead uses a [partitioned analysis](https://iqtree.github.io/doc/Advanced-Tutorial), where IQ-TREE is used to identify the best substition model for every BUSCO gene used. 

I also recommend checking out [MycoTools](https://github.com/xonq/mycotools) - it is great for downloading and organizing your genome data prior to using cano.py. Cano.py works very well with the "ome"-based naming system that MycoTools uses.

Canopy is currently coded to work seamlessly on SCINet for users in the "arsef" project. If you are on a different system or do not have access to the arsef project, it will take a bit more work to use cano.py (still possible with miminal effort).

What cano.py does:  with given input, finds all Single Copy Orthologs (SCOs), creates individual gene trees in IQ-TREE with the modelfinder option, then runs a partitioned analysis in IQ-TREE - the final result being a species tree with bootstrap values.

**cano.py input:**  BUSCO output or OrthoFinder output for a set of genomes

**cano.py output:**  IQ-TREE species tree (partitioned analysis), individual IQ-TREE gene tree, various useful metadata files

<br>

---

<br>

# Assumptions

This script assumes you have a collection of genomes and have already completed either a BUSCO analysis for each genome (using the same BUSCO database for each genome), or have completed a standard or scaled OrthoFinder run. The BUSCO input and the OrthoFinder input are mutually exclusive.

This script also assumes that you have access to the arsef shared conda environment "canopy". If you don't, that's fine. Just make sure that you have iqtree (v3.0.1), trimal (v1.5.rev0), and muscle (v5.1) in your PATH.

## BUSCO input 

Cano.py will only take BUSCO results from BUSCO versions >3, due to how the output is formatted.

Cano.py also assumes you have stored this BUSCO output in a folder named the same as the "ome" code ("acitak1" in the given example). For instance, the BUSCO output for acitak1 is stored here: ./busco/fungi/by_ome/acitak1 . All BUSCO output folders should be stored in this same parent folder in a similar manner.

SCINet users:  I have run BUSCO ("fungi" database) on the entirety of the MycoTools database. The BUSCO output folder can be found here: /project/arsef/databases/mycotools/database_stats/busco/fungi/by_ome/

Here is an example BUSCO command, used to generate the BUSCO results:

```bash
# SCINet users can use the shared environment "funannotate"
module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate
busco -i acitak1.fna -o "acitak1" -m genome -l "/project/arsef/databases/funannotate_databases/fungi_odb10"
```

## OrthoFinder input 

Cano.py requires the absolute path to the OrthoFinder output folder called "Single_Copy_Orthologue_Sequences". 


<br>

# Running cano.py

Provide a list of the "omes" you want to include in the analysis, the path to your BUSCO outputs/OrthoFinder output, the number of threads you want to use, the number of bootstraps desired for the gene tree and species tree creation, and a path to your desired output. I recommend running as a job, as this may take several days (test case:  took 1.5 days for 130 Hypocreales genomes with 64 threads).

```bash
module load miniconda
source activate /project/arsef/environments/canopy

python /project/arsef/scripts/cano.py \
--omes testset1_ome_list.txt \
--busco /project/arsef/databases/mycotools/database_stats/busco/fungi/by_ome/ \ # BUSCO input example
# --orthofinder_sco_dir /OrthoFinder/Results_Nov20/Single_Copy_Orthologue_Sequences \ # OrthoFinder input example
--threads 64 \
--bootstraps 1000 \
--seq_min_length 300 \
--sco-threshold 1.0 \
--output /project/arsef/projects/hypocreales_tree/testset1/canopy/
```

Other options:

`sco-threshold` : the minimum proportion of genomes with a particular gene to retain it as single-copy (default: 1.0). 
      Sometimes, if you have poor quality genomes, a very diverse set of genomes, or a very large set of genomes, you will be unable to find any true SCOs in your dataset. Lower the sco-threshold until you have enough genes that make the cut. So if you have 100 genomes and set the sco-threshold to 0.9 (90%), only 90 genomes will need to have a particular gene complete and single-copy for the gene to be treated as an SCO.

`seq_min_length` : the minimum length threshold (aa) that an aligned and trimmed SCO must surpass in order to move forward in the pipeline.

<br>

# Explanation of output

`busco_genes_presence_absence_matrix_FULL.tsv` : a matrix indicating the presence/absence of BUSCO genes (single copy and complete) for each genome, for every single BUSCO gene in a given database (e.g. if you provide the fungi odb10 BUSCO output, this will always include 758 fugnal BUSCO genes). This file will only be generated if BUSCO input is provided. 

`busco_genes_presence_absence_matrix.tsv` : a matrix indicating the presence/absence of BUSCO genes (single copy and complete) for each genome, but only for the BUSCOs that were identified as passing your SCO threshold. This file will only be generated if BUSCO input is provided. 

`single_copy_sequences` : folder containing aa fastas of genes, one fasta for each determined SCO/pseudo-SCO

`alignments` : MUSCLE alignments of the fastas from single_copy_sequences

`trimmed_alignments` : alignments trimmed with trimAl

`single_gene_trees` : output of IQ-TREE for each individual aligned+trimmed SCO

`concatenated_alignment` : a concatenated alignment fasta prepared for the IQ-TREE partitioned analysis. Also contains the partition file (.nexus) and partitioned coordinates file (.txt)

`final_tree` : output of the IQ-TREE partitioned analysis. The .contree file has the ML tree structure and the UF bootstrap values at each node.

`final_report.txt` : a useful text file recording parameters used, number of genes identified as SCOs, runtime, etc. 

`logs` :  some logfiles that are useful in case you have an error.

  `alignment_lengths.tsv` : a list of all the SCOs identified, their legnths, and if they were removed (filtered) or retained (kept) in the pipeline, based on the given `seq_min_length` value. 
  
  `genomes_with_no_buscos.txt` : a list of samples where no BUSCO SCOs were identified. If populated, this means that you should probably remove the offending samples (poor quality, bad matchup, etc). This file will only be generated if BUSCO input is provided. 

  `missing_busco_dirs.txt` : a list of samples that were in your input ome list, but couldn't find their BUSCO output. If populated, if means something is wrong with your folder structure. This file will only be generated if BUSCO input is provided. 

  



<br>

---

<br>

# Citations

If you used cano.py, please cite all the software wrapped into this script. If you want, you can also cite this webpage for cano.py itself. I have listed the accesspory software versions used in the shared conda environment on SCINet. 

**IQ-TREE (v3.0.1)** : Wong, T.K., Ly-Trong, N., Ren, H., Baños, H., Roger, A.J., Susko, E., Bielow, C., De Maio, N., Goldman, N., Hahn, M.W. and Huttley, G., 2025. IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models.

**trimal (v1.5.rev0 build[2024-05-27])** : Capella-Gutiérrez, S., Silla-Martínez, J.M. and Gabaldón, T., 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics, 25(15), pp.1972-1973.

**Muscle (5.1.linux64)** :  Edgar, R.C., 2022. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. Nature communications, 13(1), p.6968.

---