library(readr)

df <- read_csv("1-Input/national-history.csv")

library(dplyr)
df2 <- df %>% select(date, deathIncrease, positiveIncrease)

library(ggplot2)
p1 <- ggplot(df2, aes(x = date, y = deathIncrease)) + geom_line()
p1

ggsave(
  "Plot1.png",
  units = "px",
  dpi = 300,
  width = 1600,
  height = 900,
  path = "./3-Output/"
)

p2 <- ggplot(df2, aes(x = date, y = positiveIncrease)) + geom_line()
p2

ggsave(
  "Plot2.png",
  units = "px",
  dpi = 300,
  width = 1600,
  height = 900,
  path = "./3-Output/"
)
