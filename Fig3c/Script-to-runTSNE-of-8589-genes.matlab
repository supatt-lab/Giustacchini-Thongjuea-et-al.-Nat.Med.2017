%%##########################################################
%%#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
%%#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
%%#Maintainer: Supat Thongjuea and Alice Giustacchini
%%#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
%%#Journal : Nature Medicine
%%#Year : 2017
%%##########################################################
%%#This is the Matlab script run on MATLAB_R2014b.
%%#The tSNE software was downloaded from https://lvdmaaten.github.io/tsne/
%%##########################################################
%%####1. import the data from the text file using MATLAB####
DAT = table2cell(NormalHSCsandDiagnosisCells);
CML_matrix=cell2mat(DAT);
B=transpose(CML_matrix);

%%####2. Run tSNE with the perplexity =30#####################
%%####The perplexity = 30 was used due to the high number of cells (n=1,086 cells).
%%####The similar pattern of clusters would also be observed using the lower perplexity (e.g. 20 and 25) 
%%############################################################
tSNE = tsne(B,[],2,30,30)

%%####3. Scatter plot when using the Dim1 and Dim2 as a reduced dimension.
scatter(tSNE(:,1),tSNE(:,2),[],bcrablstatusNormalHSCsandDiagnosisCells)
