# Defined by Shape: Elucidating the Molecular Recognition of Dynamic Loops with Covalent Ligands
Protein loops harness conformational heterogeneity to perform an array of functions, ranging from catalyzing enzymatic reactions to communicating allosteric signals. Although attractive targets for small molecule modulation, these functional hubs are often considered unligandable due to their lack of well-defined binding pockets and highly dynamic structure. Recent studies, however, have demonstrated the power of covalent chemistry to selectively capture cryptic pockets formed by protein loops. Herein, we leverage machine learning to elucidate the molecular basis of covalent ligand:loop recognition in the transcriptional coactivator Med25. Key to our success was classification by ligand shape prior to model training, which led to descriptive and predictive models. The models were experimentally validated through the synthesis and in vitro testing of novel top-ranked ligands, revealing canonical structure-affinity relationships, including an activity cliff. Further feature analyses identified traditional topological and spatial parameters predictive of binding, and molecular modeling uncovered a potential binding pocket with at least two distinct conformations with high shape complementarity. Collectively, these findings reveal the hidden potential of dynamic loops as specific sites for covalent small molecule modulation, challenging the notion that protein loops are unligandable and demonstrating their capacity for exquisite, shape-based molecular recognition.

## Requirements
Python == 3.7.3 
* RDKit
* pandas
* math

R Statistical Software == 4.0
* stringi
* caret
* ggplot2
* withr
* Formula
* plotmo
* plotrix
* TeachingDemos
* recipes
* dplyr
* crayon
* Matrix
* lattice

## File Descriptions
* Descriptor-BoltzmannAverage.ipynb\
**Description**: Script calculates a Boltzmann average for each 3D descriptor, weighted by the conformational energy of the fragment.

* DataProcessing-ModelBuildingLoop.R\
**Description**: Script removes invariant and highly correlated descriptors, standardizes descriptors, builds machine learning models, and performs statistical analysis. 

* CShape_Search.ipynb\
**Description**: Script searches the Tethering library for C-shaped fragments. The C-shape molecular pattern is defined by a central ring bearing two substituents ortho to one another: a disulfide linker and a second ring system. 

* LinearShape_Search.ipynb\
**Description**: Script searches the Tethering library for linear-shaped fragments. The linear-shape molecular pattern is defined by the presence of a central carbonyl group that is connected to i) aromatic ring(s) and ii) a piperidine ring, which is also bonded to the disulfide linker.

* StratifiedSplitting.ipynb\
**Description**: Script executes a stratified splitting method that partitions the fragment library into four classes and then randomly assigns the fragments to training (75%) and test (25%) sets.

* YScrambling.R\
**Description**: Script executes y-scrambling for the machine learning algorithms and feature selection methods.

* ZINC-CShape_Search.ipynb\
**Description**: Script searches the ZINC Database for C-shaped fragments. The C-shape molecular pattern is defined by a central ring bearing two substituents ortho to one another: a disulfide linker and a second ring system. Additional restraints were added to search the ZINC database: (1) to add the disulfide linker, the fragments were required to have a primary amine or carboxylic acid no more than two atoms from the central ring; (2) functional groups that if present, the fragment was removed: acrylamides, azoles, disulfides, hydrazines, hydrazones, and imidines; (3) if the amine was a part of a double bond or had double bond character the fragments was removed, where the latter was defined as a nitrogen involved in a resonance structure and the double bond was not in a ring (e.g. amide).

* ZINC-CShape_StructureMod.ipynb\
**Description**: Script modifies ZINC fragments to match the modeled structure of the Tethering fragments by adding a linker.

* ZINC-ConsensusPredictions.R\
**Description**: Script predicts Tethering percentages for modeled ZINC fragments.

* ZINC-ApplicabilityDomain.R\
**Description**: Script calculates leverage values and applies y-cutoffs to QSAR consensus predictions.

## Acknowledgements
This work was supported by grants from the Burroughs Wellcome Fund Award to B.S.M. and National Institutes of Health to J.Z.V. (R35GM146888). The authors gratefully acknowledge fellowship support to P.A.S. from Mary E. Galvin Scholars Program at the University of Notre Dame, to J.M.R. from the National Institutes of Health (T32 GM145773), and to C.A.M from Rosen Fellows Program at the University of Michigan Life Sciences Institute and National Science Foundation Graduate Research Fellowship (NSF GRPF). Authors thank Drs. Bill Boggess and Mijoon Lee in the Mass Spectrometry & Proteomics Facility at the University of Notre Dame for assistance with data acquisition and technical support and Ana-Teresa Mendoza and Julia Wendhuda for their assistance with early-stage synthesis efforts. Authors acknowledge the Magnetic Resonance Research Center at the University of Notre Dame for providing access to Bruker 400 MHz, 500 MHz and Varian 600 MHz NMR instruments.

## Questions?
Email us: bmorgan3@nd.edu
