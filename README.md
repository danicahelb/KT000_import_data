# KT000_import_data
Import Apac data from the KT array
##########################################################################################
##########################################################################################

rm(list=ls()) 

setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/") 


##########################################################################################
#Load array data:

csv_files <- list.files("./Data/Raw Data/")

for(i in csv_files){

  temp <- read.csv(paste("./Data/Raw Data/", i, sep=""), skip=3, as.is = c(1,4:13))
  
  data_name <- gsub("_$","",gsub(" ","", gsub("[S/s]lide","", gsub("^new_Apac","",gsub("[.]csv$", "", i)))))
  
  temp$slideN <- gsub("^.+_", "", data_name)
  
  temp$study <- paste("Apac", gsub("_.+$", "", data_name), sep="_")
  
  
  temp$subarrayN <- 1
  temp$subarrayN[temp$Block %in% c("Block5", "Block6", "Block7", "Block8")] <- 2
  temp$subarrayN[temp$Block %in% c("Block9", "Block10", "Block11", "Block12")] <- 3
  temp$subarrayN[temp$Block %in% c("Block13", "Block14", "Block15", "Block16")] <- 4
  temp$subarrayN[temp$Block %in% c("Block17", "Block18", "Block19", "Block20")] <- 5
  
  names(temp)[grep("Raw.Mean", names(temp))] <- "Raw.Mean"
  names(temp)[grep("Background.Mean", names(temp))] <- "Bkgrd.Mean" 
  names(temp)[grep("Foreground.Mean", names(temp))] <- "Frgrd.Mean"
  names(temp)[grep("Raw.Median", names(temp))] <- "Raw.Median"
  names(temp)[grep("Background.Median", names(temp))] <- "Bkgrd.Median" 
  names(temp)[grep("Foreground.Median", names(temp))] <- "Frgrd.Median"
  
  do.call("<-",list(data_name, temp))   
}


##########################################################################################
#merge data in separate slides together

for(i in gsub("^X", "", ls()[grep("^X", ls())])){
  temp <- get(paste("X", i, sep=""))
  print(paste("X", i, sep=""))
  if(exists("d1")){
    d1 <- rbind(d1, temp)
  } else {
    d1 <- temp
  }
}

names(d1)[4] <- "antigen"
names(d1)[5] <- "sampleID"
d1$slideN <- as.numeric(d1$slideN)
for(i in c(6:11)){
  d1[,i] <- gsub(",", "", d1[,i])
  d1[,i] <- as.numeric(d1[,i])
}


##########################################################################################
#determine position

d1$MC <- 1
d1$MC[d1$Block %in% c("Block2", "Block6", "Block10", "Block14", "Block18")] <- 2
d1$MC[d1$Block %in% c("Block3", "Block7", "Block11", "Block15", "Block19")] <- 3
d1$MC[d1$Block %in% c("Block4", "Block8", "Block12", "Block16", "Block20")] <- 4
table(d1$MC)

names(d1)[2] <- "SR"
names(d1)[3] <- "SC"
d1$MR <- d1$subarrayN

d1$position <- paste(d1$MC, d1$SR, sep="_")
d1$position <- paste(d1$position, d1$SC, sep=".")


##########################################################################################
#determine which antigens are in which positions on the array & create unique identifiers for duplicate spots

d1 <- d1[order(d1$study, d1$slideN,d1$MR, d1$SR, d1$MC, d1$SC),]
d1$INDEX <- c(1:256)
d1$antigen <- gsub("_", "-", d1$antigen) 
d1$antigen[grep("lank", d1$antigen)] <- "blank"         # "Blank" and "blank" now are both "blank"
d1$spot <- paste(gsub(" ", "", d1$antigen), d1$INDEX, sep="_")



##########################################################################################
#clean up array data 

d1 <- d1[d1$antigen!="REF",]
d1 <- d1[order(d1$study, d1$slideN, d1$sampleID, d1$antigen, d1$INDEX, d1$spot),
         c("study", "sampleID", "antigen", "spot", "INDEX", "slideN", "subarrayN", "Raw.Mean", "Bkgrd.Mean", "Frgrd.Mean", "Raw.Median",
           "Bkgrd.Median", "Frgrd.Median", "position", "Block", "MR", "MC", "SR", "SC", "Flag", "Annotation")]
rownames(d1) <- 1:nrow(d1)


##########################################################################################
#save array data

write.csv(d1, file=paste("./Data/Processed Data/", gsub("-", "", Sys.Date()), "_", "import_data.csv", sep=""), row.names=F)
save(d1, file=paste("./Data/Processed Data/", gsub("-", "", Sys.Date()), "_", "import_data.RData", sep="")) 
