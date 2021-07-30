# SurfaltR
Surface proteins are hydrophobic and remain difficult to study thereby necessitating the use of TM topology prediction methods such as TMHMM (1) and Phobius (2). However, there exists a need for bioinformatic approaches to streamline batch processing of isoforms for comparing and visualizing topologies. To address this gap, we have developed an R package, SurfaltR. It pairs inputted isoforms, either known alternatively spliced or novel, with their APPRIS (3) annotated principal counterparts, predicts their TM topologies using TMHMM or Phobius, and generates a customizable graphical output. Further, SurfaltR facilitates the prioritization of biologically diverse isoform pairs through the incorporation of three different ranking metrics and through protein alignment functions.

Note: If you use surfaltR in published research, please cite:
<publication> 
 
# Installation 
<add details on installation from github/CRAN>
  
# TMHMM standalone software Installation
In order to be able to use TMHMM R package4 within surfaltR to predict membrane topology, it is important to first ensure that you have TMHMM 2.0 standalone software installed on your computer. To do this, simply navigate to https://services.healthtech.dtu.dk/service.php?TMHMM-2.0, and follow directions for installation of standalone software. In order to install TMHMM 2.0 in your R environment, you will also need the package “tmhmm”. The package should automatically install when you download the surfaltR package. In the event that this does not happen, you can use the following installation code:
install.packages(“tmhmm”)
Once you have obtained your link to install TMHMM 2.0 and successfully loaded the “tmhmm” package, you will need to use the following code to make TMHMM operable within your R development environment:
library(“tmhmm”)
install_tmhmm("https://services.healthtech.dtu.dk/download/28c408dc-ef5e-47ad-a284-66754bcd27f7")
In the code above, be sure to replace the URL shown in the quotation marks with the URL emailed to you after requesting the TMHMM 2.0 download. 

# Phobius Installation
As run_phobius() relies on the Phobius API, a copy of the software does not need to be downloaded on the user’s device. However, in order to ensure that sequences can be adequately processed in the R development environment, the “ragp” package needs to be installed. To install this package, the following code can be used:
devtools::install_github("missuse/ragp")
  
# Please refer to the vignette for detailed descriptions of wrokflow and functions.


# References
1. 	Sonnhammer EL, von Heijne G, Krogh A. A hidden Markov model for predicting transmembrane helices in protein sequences. Proc Int Conf Intell Syst Mol Biol. 1998;6:175–82. 
2. 	Käll L, Krogh A, Sonnhammer ELL. A combined transmembrane topology and signal peptide prediction method. J Mol Biol. 2004 May 14;338(5):1027–36. 
3. 	Rodriguez JM, Rodriguez-Rivas J, Di Domenico T, Vázquez J, Valencia A, Tress ML. APPRIS 2017: principal isoforms for multiple gene sets. Nucleic Acids Res. 2018 Jan 4;46(D1):D213–7. 
