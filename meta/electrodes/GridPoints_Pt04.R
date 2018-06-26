library('png')
library('spdep')
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

# --- Pt04

# Check
img <- png::readPNG("C:/Users/mbmhscc4/GitHub/ECoG_Data_Prep/meta/electrodes/Map_Pt03.png")
plot(0, type='n', xlim=c(0,350), ylim=c(0,350), main="Not the best use, but this gives the idea",xlab='x',ylab='y')
rasterImage(img, 0, 0, 350, 350)
points(dro$x,dro$y,cex=2.2)
text(dro$x,dro$y,labels = dro$labels, cex=0.5)
