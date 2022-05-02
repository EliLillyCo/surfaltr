# surfaltr
Surface proteins are hydrophobic and remain difficult to study thereby necessitating the use of TM topology prediction methods such as TMHMM (1) and Phobius (2). However, there exists a need for bioinformatic approaches to streamline batch processing of isoforms for comparing and visualizing topologies. To address this gap, we have developed an R package, SurfaltR. It pairs inputted isoforms, either known alternatively spliced or novel, with their APPRIS (3) annotated principal counterparts, predicts their TM topologies using TMHMM or Phobius, and generates a customizable graphical output. Further, SurfaltR facilitates the prioritization of biologically diverse isoform pairs through the incorporation of three different ranking metrics and through protein alignment functions.

# Note: 
If you use surfaltr in published research, please cite this page and possibly a subsequent publication (will be updated later).
 
# Installation: 
As surfaltr is hosted on Bioconductor, please follow installation instructions outlined here: https://bioconductor.org/packages/release/bioc/html/surfaltr.html 
  
# Please refer to the vignette for detailed descriptions of installation, workflow, functions and troubleshooting.


# References
1. 	Sonnhammer EL, von Heijne G, Krogh A. A hidden Markov model for predicting transmembrane helices in protein sequences. Proc Int Conf Intell Syst Mol Biol. 1998;6:175–82. 
2. 	Käll L, Krogh A, Sonnhammer ELL. A combined transmembrane topology and signal peptide prediction method. J Mol Biol. 2004 May 14;338(5):1027–36. 
3. 	Rodriguez JM, Rodriguez-Rivas J, Di Domenico T, Vázquez J, Valencia A, Tress ML. APPRIS 2017: principal isoforms for multiple gene sets. Nucleic Acids Res. 2018 Jan 4;46(D1):D213–7. 
