#essential:
print("installing essential packages:")
install.packages(c("survival","mgcv","survPen","dplyr"))
#recommended:
print("installing recommended packages:")
install.packages(c("ggplot2","cowplot","scico","gridExtra","ggrepel"))
#optionals:
try(library(psych, character.only = TRUE), silent = TRUE)
try(library(segmented, character.only = TRUE), silent = TRUE) #used for breakpoints
try(library(strucchangeRcpp, character.only = TRUE), silent = TRUE) #used for breakpoints

