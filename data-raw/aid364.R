# read in response

data <- read.csv("AID_364_response.csv")

# read in descriptor sets compute number of descriptors

desc_lengths <- c()
for(desc in c("BurdenNumbers.csv","Pharmacophores.csv")){
  d <- read.csv(desc)[-1]
  data <- cbind(data, d)
  desc_lengths <- c(desc_lengths, ncol(d))
}

dim(data)

# make the dataset smaller so that it will run faster
aid364 <- data[1:500, ]