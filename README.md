# Surface MOTIF detection
The overall project is part of the dissertation's final stage. Due to lack of time, this project is on hold.

# Discription
The underlying hypothesis only considers part of the surface residues (srs-MOTIF) to classify the proteins' functionality on the surface of the protein molecules.
To evaluate the hypothesis, we utilized 3D structures (PDBs) of Cancer-related functions that promote (ONGO) or suppress (TSG) cell growth.

This repository contains the projects(approaches) related to my dissertation work. Thus, the codes and reports presented in the warehouse are copyright reserved to ULL.

# Approaches

 ## Excaustive_search ## 

We search all possible MOTIFs on the surface  (srs-MOTIF) of the 3D (PDB) structure.

During the exhaustive search, only the presence of srs-MOTIF in each functionality is evaluated. 
Unfortunately, the lack of structure details in the search procedure makes it difficult to support the hypothesis.

To obtain more details, please check the report "Excaustive_search_Report_Part_from_my_desseration_work" attached in the folder.
 
 ## Statistics_of_Patterns_find_Possible-SRS-MOTIFs ## 

Instead of exhaustively searching, we can use the statistics of the central residue and the residues in the SRS-MOTIF.
The initial statistics are promising and can guide us in this direction.
 
 To obtain more details, please check the report "Stats_SRS_MOTIF_report_Part_from_my_desseration_work" attached in the folder.
 
 Currently, the project is on hold.


 ## Deeplearning_based_approach ## 
  Proposed a CNN-deep learning model to evaluate the hypothesis. 
  First, suspected small structures (srs-MOTIFs) were obtained and represented.
  I suspect that most small structures with no functionality assigned with functionality may produce noise.

  To obtain more details, please check the report "MOTIF_deep_learning_GitHub" attached in the folder.

  ## Future direction
  Due to the continuous effort of the researchers and the rapid increase in deposits of the PDB structures, 
  this hypothesis may be supported based on the idea proposed in this report. 
    First, obtain the functionally relevant structures based on the statistics, then use the deep learning model
  to find the undetected srs-MOTIFs.

  
