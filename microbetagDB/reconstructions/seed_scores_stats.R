
setwd("/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/")

file_list <- list.files("/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/scores_subset/")
file_list <- list.files("/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/scores/")


# Initialize an empty dataframe
competition_df <- data.frame()
cooperation_df <- data.frame()

# Create a list to store the third columns from each file
competition_list <- list()
cooperation_list <- list()

# Loop through the list of files
for (file_name in file_list) {

  # Read the file into a dataframe
  # data <- read.table(paste("scores_subset/", file_name, sep=""), header = FALSE)
  data <- read.table(paste("scores/", file_name, sep=""), header = FALSE)
  
  # Store the third column (column 3) in the list
  competition_list[[file_name]] <- data[, 3]
  cooperation_list[[file_name]] <- data[, 4]
}


# Combine the columns from the list into the pairwise dataframe
competition_df <- do.call(cbind, competition_list)
cooperation_df <- do.call(cbind, cooperation_list)

# Set column names if needed
colnames(competition_df) <- file_list

character_value <- as.character(data[,2])
character_value <- sub("\\.0*$", "", data[,2])

rownames(competition_df) <- character_value
rownames(cooperation_df) <- character_value



# Compute the mean values of each column
comp_means <- colMeans(competition_df, na.rm = TRUE)
coop_means <- colMeans(cooperation_df, na.rm = TRUE)


# Create histograms using all values of competition and cooperation scores
hist(competition_df, main = "Histogram of Competition values", xlab = "Competition mean values", ylab = "Frequency", col = "lightgreen", border = "black")
hist(cooperation_df, main = "Histogram of Co-operation values", xlab = "Co-operation mean values", ylab = "Frequency", col = "lightgreen", border = "black")

# Create a histogram of the mean values
hist(comp_means, main = "Histogram of Competition mean values", xlab = "Competition mean values", ylab = "Frequency", col = "lightblue", border = "black")
hist(coop_means, main = "Histogram of Co-operation mean values", xlab = "Co-operation mean values", ylab = "Frequency", col = "lightblue", border = "black")






# Hierachical clustering -- this will take ages so think of moving to a cluster or something
hc_comp <- hclust(dist(competition_df))
hc_coop <- hclust(dist(cooperation_df))







