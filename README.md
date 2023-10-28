# Surface MOTIF detection
The overall project is part of the dissertation's final stage. Due to lack of time, this project is on hold.

# Discription
The underlying hypothesis only considers part of the surface residues (srs-MOTIF) to classify the proteins' functionality on the surface of the protein molecules.
To evaluate the hypothesis, we utilized 3D structures (PDBs) of Cancer-related functions that promote (ONGO) or suppress (TSG) cell growth.

This repository contains the projects(approaches) related to my dissertation work. Thus, the codes and reports presented in the warehouse are copyright reserved to ULL.

# Approaches

 ## Excaustive search ## 

We search all possible MOTIFs on the surface  (srs-MOTIF) of the 3D (PDB) structure.

During the exhaustive search, only the presence of srs-MOTIF in each functionality is evaluated. 
Unfortunately, the lack of structure details in the search procedure makes it difficult to support the hypothesis.

To obtain more details, please check the report "Excaustive_search_Report_Part_from_my_desseration_work" attached in the folder.
 
 ## Statistics of Patterns to Find Possible SRS-MOTIFs ## 

Instead of exhaustively searching, we can use the statistics of the central residue and the residues in the SRS-MOTIF.
The initial statistics are promising and can guide us in this direction.
 
 To obtain more details, please check the report "Stats_SRS_MOTIF_report_Part_from_my_desseration_work" attached in the folder.
 
 Currently, the project is on hold.


 ## Deep learning based approach ## 
  Proposed a CNN-deep learning model to evaluate the hypothesis. 
  First, suspected small structures (srs-MOTIFs) were obtained and represented.
  I suspect that most small structures with no functionality assigned with functionality may produce noise.

  To obtain more details, please check the report "MOTIF_deep_learning_GitHub" attached in the folder.

  ## Future direction
  Due to the continuous effort of the researchers and the rapid increase in deposits of the PDB structures, 
  this hypothesis may be supported based on the idea proposed in this report. 
    First, obtain the functionally relevant structures based on the statistics, then use the deep learning model
  to find the undetected srs-MOTIFs.
  
## Contact
For more details, feel free to reach me via LinkedIn
www.linkedin.com/in/nishanth-anand

If you find the approaches and the results presented are helpful, please cite.

[1] Anandanadarajah, N., Chu, C. H., & Loganantharaj, R. (2021). An integrated deep learning and dynamic programming method for predicting tumor suppressor genes, oncogenes, and fusion from PDB structures. In Computers in Biology and Medicine (Vol. 133, p. 104323). Elsevier BV. https://doi.org/10.1016/j.compbiomed.2021.104323


