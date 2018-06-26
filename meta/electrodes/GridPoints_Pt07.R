GridPoints_Pt07 <- function(plotgrid = FALSE) {
    # --- Pt7
    Pt <- 7

    # Grid I
    I1 <- c(158, 276)
    I5 <- c(137, 210)
    I16 <- c(205, 262)
    I20 <- c(184, 195)
    grid <- defineGrid(dims = c(4,5), I1, I5, I16, label='I')
    GI <- grid$x

    # Grid J
    J1 <- c(233, 232)
    J3 <- c(205, 212)
    ((J1[1]-J3[1]) / 3) * 4
    ((J1[2]-J3[2]) / 3) * 4
    J4 <- c(195, 203)
    grid <- defineGrid(dims = c(1,4), J1, J4, label='J')
    GJ <- grid$x

    # Grid K
    K1 <- c(280, 246)
    K6 <- c(213, 199)
    grid <- defineGrid(dims = c(1,6), K1, K6, label='K')
    GK <- grid$x

    AllGrids <- rbind(GI,GJ,GK)

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
