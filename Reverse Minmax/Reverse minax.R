Dataframe <- c()

Dataframe$letters <- letters
Dataframe <- data.frame(Dataframe)

# Factorization
Dataframe$letters_factor <- as.factor(Dataframe$letters)

print(class(Dataframe$letters))
print(class(Dataframe$letters_factor))

#Extract levels labels 
levels <- levels(Dataframe$letters_factor)

#Convert levels to numeric 
Dataframe$letters_factor_numeric <- as.numeric(Dataframe$letters_factor)

normalize <- function(x)
{
  return((x - min(x)) / (max(x) - min(x)))
}

# MinMax Normalization
min <-  min(Dataframe$letters_factor_numeric)
max <-  max(Dataframe$letters_factor_numeric)

Dataframe$letters_factor_minmax <-
  normalize(Dataframe$letters_factor_numeric)

# Reverse Minmax
Dataframe$letters_factor_reverse <-
  Dataframe$letters_factor_minmax * (max - min) + min

# Reverse Labeling
temp <- c()
for (i in 1:length(levels)) {
  temp[i] <- gsub(i, levels[i], Dataframe$letters_factor_reverse[i])
}

Dataframe$letters_reverse  <- temp
rm(temp)
