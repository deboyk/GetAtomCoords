---------------
GetAtomCoords
---------------

A MATLAB function to retrieve the coordinates for all atoms in a Materials Studio file (or files), and subsequently calculate the molecular descriptor R3m for each molecule in the file. 

This function allows the user to select data exported from Materials Studio as 3D atomistic files (.xsd file format) and calculate R3m values for each molecule in each frame. A distribution of R3m values is then plotted as histogram. The goal is to determine the distribution of R3m values that result from conformational differences in molecular structures as determined from Molecular Dynamics Simulations performed using Materials Studio. 

This package also contains the following functions: R3mCalculate_auto, AtomConnection, AtomicWeighting, EuclidDistance,InfluenceDistanceMat, MolecInfluenceMatrix, MolecMatrix

These functions are adapted from the R3mCalculate function ([Mathworks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/69891-r3mcalculate?s_tid=prof_contriblnk) and [GitHub](https://github.com/deboyk/R3mCalculate) links). The primary difference is that the R3mCalculate_auto function has been edited so that the user does not need to manually select each file. 


**NOTE:** The 'GetAtomCoords.m' function is the primary function. The other functions and scripts are 
      referenced by this function. The user will only need to run the R3mCalculate function:
      e.g. [r3m_out] = GetAtomCoords();


## References:  
[1] Todeschini, R.; Consonni, V., Molecular Descriptors for Chemoinformatics, Volume 41 (2 Volume Set). 2nd ed.; John Wiley & Sons: 2009.  
[2] Consonni V, Todeschini R, Pavan M, Gramatica P 2002. Structure/response correlations and similarity/diversity analysis by GETAWAY descriptors. 2. Application of 	the novel 3D molecular descriptors to QSAR/QSPR studies. Journal of chemical information and computer sciences  42(3):693-705.  
[3] Tetko, I. V.; Gasteiger, J.; Todeschini, R.; Mauri, A.; Livingstone, D.; Ertl, P.; Palyulin, V. A.; Radchenko, E. V.; Zefirov, N. S.; Makarenko, A. S.  Virtual 	computational chemistry laboratory–design and description. J. Comput. Aided Mol. Des. 2005, 19, (6), 453-463.  
[4] DeBoyace, K.; Buckner, I. S.; Gong, Y.; Ju, T.-c. R.; Wildfong, P. L. D.  Modeling and Prediction of Drug Dispersability in Polyvinylpyrrolidone-Vinyl Acetate 	Copolymer Using a Molecular Descriptor. J. Pharm. Sci. 2018, 107, (1), 334-343.  

R3m can also be calculated online with appoximately 1600 other descriptors at http://www.vcclab.org/lab/edragon/

 
### Author:	
 	Kevin DeBoyace  
 	Wildfong Lab  
 	Duquesne University  
  
 Last Updated: January 2019  
 Uploaded to GitHub: March 2022
