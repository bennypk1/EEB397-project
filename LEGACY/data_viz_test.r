library(ggplot2)
library(dplyr)

iris %>% ggplot(aes(x = Sepal.Length, y = Sepal.Width)) + geom_point()