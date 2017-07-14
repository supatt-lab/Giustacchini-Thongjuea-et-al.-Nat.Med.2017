#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017

setClass(
		Class="SingleCellFluidigm",
		representation(
				biomark_files="vector",
				chip_names="vector",
				data_dir="vector",
				LOD="numeric",
				org.data="data.frame",
				outlier.data="data.frame"
		),
		prototype(
				biomark_files=character(1),
				chip_names=character(1),
				data_dir=character(1),
				LOD=numeric(1),
				org.data=data.frame(),
				outlier.data=data.frame()
		),
		validity = function(object) {
			if(length(object@biomark_files)==1 & object@biomark_files[1]=="")
				return( "Please input your BiomarkFiles!!" )
			##############################################
			if(length(object@chip_names)==1 & object@chip_names[1]=="")
				return( "Please input your unique chip names" )
			##############################################
			if(object@LOD==0)
				return( "Please input the limit of detection (LOD) CT value" )
			if(length(unique(object@biomark_files)) != length(unique(object@chip_names)))
				return( "Your chip names must be unique for a single Biomark file." )
		}
)
