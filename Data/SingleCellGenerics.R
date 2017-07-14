#############################################################
#Author: Supat Thongjuea, MRC Molecular Haematology Unit, Weatherall Institute of Molecular Medicine, University of Oxford, UK 
#Contact email : supat.thongjuea@ndcls.ox.ac.uk or supat.thongjuea@gmail.com
#Maintainer: Supat Thongjuea and Alice Giustacchini
#Title: Single-cell Transcriptomics Uncovers Distinct and Clinically Predictive Molecular Signatures of Stem Cells in Chronic Myeloid Leukemia
#Journal : Nature Medicine
#Year : 2017
####################################
#########  AllGenerics.R
#########
############Common methods##########
setGeneric(
		name="biomarkFiles",
		def=function(object){
			standardGeneric("biomarkFiles")
		})
setMethod("biomarkFiles",
		signature(object = "SingleCellFluidigm"),
		function (object){
			object@biomark_files
		})
########
setGeneric(
		name="SingleCellChipNames",
		def=function(object){
			standardGeneric("SingleCellChipNames")
		})
setMethod("SingleCellChipNames",
		signature(object = "SingleCellFluidigm"),
		function (object){
			object@chip_names
		}
)

#########
setGeneric(
		name="SingleCellGroupNames",
		def=function(object){
			standardGeneric("SingleCellGroupNames")
		}
)
setGeneric(
		name="SingleCellGroupNames<-",
		def=function(object,value){
			standardGeneric("SingleCellGroupNames<-")
		}
)
setMethod("SingleCellGroupNames",
		signature(object = "SingleCellFluidigm"),
		function (object){
			object@group_names
		}
)
setReplaceMethod(
		f="SingleCellGroupNames",
		signature="SingleCellFluidigm",
		definition=function(object,value){
			initialize(object,SingleCellGroupName=value)
		})
#########
setGeneric(
		name="SingleCellDir",
		def=function(object){
			standardGeneric("SingleCellDir")
		})
setMethod("SingleCellDir",
		signature(object = "SingleCellFluidigm"),
		function (object){
			object@data_dir
		}
)
#########
setGeneric(
		name="LOD_PCR",
		def=function(object){
			standardGeneric("LOD_PCR")
		})
setMethod("LOD_PCR",
		signature(object = "SingleCellFluidigm"),
		function (object){
			object@LOD
		}
)
#########
setGeneric(
		name="CtValues",
		def=function(object){
			standardGeneric("CtValues")
		}
)
setGeneric(
		name="CtValues<-",
		def=function(object,value){
			standardGeneric("CtValues<-")
		}
)
setMethod("CtValues",
		signature(object="SingleCellFluidigm"),
		function(object){
			object@org.data
		}
)
setReplaceMethod(
		f="CtValues",
		signature="SingleCellFluidigm",
		definition=function(object,value){
			initialize(object,org.data=value)
		})
###########
###########
setGeneric(
		name="OutlierCtValues",
		def=function(object){
			standardGeneric("OutlierCtValues")
		}
)
setGeneric(
		name="OutlierCtValues<-",
		def=function(object,value){
			standardGeneric("OutlierCtValues<-")
		}
)
setMethod("OutlierCtValues",
		signature(object="SingleCellFluidigm"),
		function(object){
			object@outlier.data
		}
)
setReplaceMethod(
		f="OutlierCtValues",
		signature="SingleCellFluidigm",
		definition=function(object,value){
			initialize(object,outlier.data=value)
		})
###########