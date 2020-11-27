#!/usr/bin/env Rscript

# Copyright by Jan-Philipp Praetorius
# Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# https://www.leibniz-hki.de/en/applied-systems-biology.html
# HKI-Center for Systems Biology of Infection
# Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
# Adolf-Reichwein-Straße 23, 07745 Jena, Germany
# 
# This code is licensed under BSD 2-Clause
# See the LICENSE file provided with this code for the full license.

lp_svm <- function(data, ref){
  
  refIDs <- ref$ID
  
  # prepare reference
  tmp <- data[match(refIDs,data$ID),]
  tmp2 <- cbind(tmp, ref[2])
  tmp <- cbind(tmp, ref[2])
  tmp <- tmp[, -match(c("ID", "cx", "cy", "cx_gray", "cx_green", "cy_gray", "cy_green", "sf"), names(tmp))] #remove ID spalte

  tmp <- subset(tmp, tmp$class > 0)
  tmp2 <- subset(tmp2, tmp2$class > 0)

 tmp <- subset(tmp, tmp$class <= 6)
 tmp2 <- subset(tmp2, tmp2$class <=6)

  tmp <- tmp[-1]
  length <- length(tmp$class)

  # classificiation with SVM
  library(kernlab)
  
  # prepare lern data set
  list <- seq(1, length, by=3)
  list2 <- seq(2, length, by = 3);
  list3 <-  seq(3, length, by = 3);
  learnfungi <- tmp[list,]
  valds <- tmp[list2,]
  valds <- rbind(valds, tmp[list3,])
  
 fungisvm <- ksvm(class ~ ., data=learnfungi,  kernel = "rbfdot", type = "C-bsvc", prob.model=TRUE, C=10, cross = 5)
 fungisvm;

}


classify <- function(datafile, outfile, type){
	#print(datafile)

	data <- read.table(datafile, header=T)
	data130819_red <- read.table("./data130819_red.csv", header=T)
	ref7 <- read.table("./ref7.csv", header=T)
	
	classifier <- lp_svm(data130819_red, ref7);
	
	ids <- data$ID;
	ids <- na.exclude(ids);
	
	
	if(type == "probabilities"){	
		pred <- predict(classifier, data, type="probabilities");
		#pred <- predict(classifier, data, type="votes");
	 	pred <- round(pred, 3);
		result <- cbind(ids, pred);
		write.table(x = result, file=outfile, append = F, sep = "\t", row.names = F, quote = F)
	}
	
	else if(type == "response"){
		pred <- predict(classifier, data, type="response");
	
		result <- cbind(ids, pred);
		write.table(x = result, file=outfile, append = F, sep = "\t", row.names = F, col.names=c("ID", "class"), quote = F)
	}
}

print("Breakpoint 0")
args <- commandArgs(TRUE)
print("Breakpoint 1")
print(args[1])
print(args[2])

print(getwd())

setwd("./src/")

print(getwd())
classify(args[1], args[2], type="response")
