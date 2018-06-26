GridPoints_Pt09 <- function(plotgrid = FALSE) {
    # --- Pt9
    Pt <- 9

    # Grid D
    D1 <- c(165, 282)
    D5 <- c(163, 219)
    D16 <- c(214, 280)
    D20 <- c(211, 218)
    grid <- defineGrid(dims = c(4,5), D1, D5, D16, label='D')
    GD <- grid$x

    # Grid E
    E1 <- c(239, 265)
    E4 <- c(227, 214)
    grid <- defineGrid(dims = c(1,4), E1, E4, label='E')
    GE <- grid$x

    # Grid F
    F1 <- c(255,262)
    F4 <- c(240,212)
    grid <- defineGrid(dims = c(1,4), F1, F4, label='F')
    GF <- grid$x

    # Grid G
    G1 <- c(266,236)
    G4 <- c(247,188)
    grid <- defineGrid(dims = c(1,4), G1, G4, label='G')
    GG <- grid$x
    AllGrids <- rbind(GD,GE,GF,GG)

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
