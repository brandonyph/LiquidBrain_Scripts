#Remove unused memory in the environment 
gc()

#Create Sample csv 
element <- 100000
sample_csv <- matrix(runif(element),ncol=element/100)

#Write csv to HDD/SSD
write.csv(sample_csv,file="sample.csv",row.names = FALSE)

#Create empty matrix for sum
sum <- matrix(0,nrow=100)

#Calculate sum if object is in RAM all the time
system.time(
for(i in 1:100){
  sum[i] <- sum(sample_csv[i,])
}
)
#Remove and reset object 
rm(sum)
sum <- matrix(0,nrow=100)

#Calculate sum if object is in HDD/SDD and is loaded in bit by bit
system.time(
  for(i in 1:100){
    sample_csv_small <- read.csv(file="sample.csv", nrows=i)
    sum[i] <- sum(sample_csv_small)
  }
)