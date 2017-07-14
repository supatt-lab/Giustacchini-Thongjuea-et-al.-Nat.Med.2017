%%##########################################################
%%Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
%%#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
%%#Maintainer: Supat Thongjuea and Alice Giustacchini
%%#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
%%#Journal : Nature Medicine
%%#Year : 2017
%%##########################################################
%%#This is the Matlab script run on MATLAB_R2014b.
%%#The tSNE software was downloaded from https://lvdmaaten.github.io/tsne/
%%################################################################
%%####1. import the data from the text file using MATLAB####
DAT = table2cell(NORMALHSCsandK562);
CML_matrix=cell2mat(DAT);
B=transpose(CML_matrix);

%%##2. Run tSNE with the perplexity =20#####################
tSNE = tsne(B,[],2,20,20)

%%##3. Scatter plot when using the Dim1 and Dim2 as a reduced dimension.
scatter(tSNE(:,1),tSNE(:,2),[])
