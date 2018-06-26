GridPoints_Pt08 <- function(plotgrid = FALSE) {
    # --- Pt08
    Pt <- 8

    # Grid E
    grid <- list(
        x = data.frame(
            x = c(174,164,161,162,164,170),
            y = 350-c(277,266,251,236,220,206),
            grid = rep('E', 6),
            index = 1:6,
            labels = paste('E',1:6,sep='')),
        angle = NA
        )
    GE <- grid$x

    # Grid F
    F1 <- c(89, 261)
    F5 <- c(98, 195)
    F16 <- c(142, 271)
    F20 <- c(150, 207)
    grid <- defineGrid(dims = c(4,5), F1, F5, F16, label='F')
    GF <- grid$x

    AllGrids <- rbind(GE,GF)

    AllGrids$subject <- Pt

    # Check
    if (plotgrid) {
        require('png')
        img <- png::readPNG(sprintf("C:/Users/mbmhscc4/GitHub/ECoG_Data_Prep/meta/electrodes/Map_Pt%02d.png",Pt))
        plot(0, type='n', xlim=c(0,350), ylim=c(0,350), main=sprintf("Pt%02d electrode grids",Pt),xlab='x',ylab='y')
        rasterImage(img, 0, 0, 350, 350)
        points(AllGrids$x,AllGrids$y,cex=2.2)
        text(AllGrids$x,AllGrids$y,labels = AllGrids$labels, cex=0.5)
    }

    return(AllGrids)
}
