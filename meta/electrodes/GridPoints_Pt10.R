GridPoints_Pt10 <- function(plotgrid = FALSE) {
    # --- Pt10
    Pt <- 10
    # Grid F
    F1 <- c(156,313)
    F4 <- c(147,257)
    grid <- defineGrid(dims = c(1,4), F1, F4, label='F')
    GF <- grid$x

    # Grid G
    G1 <- c(185,297)
    G4 <- c(163,245)
    grid <- defineGrid(dims = c(1,4), G1, G4, label='G')
    GG <- grid$x

    # Grid H
    H1 <- c(252, 289)
    H6 <- c(186, 234)
    grid <- defineGrid(dims = c(1,6), H1, H6, label='H')
    GH <- grid$x
    AllGrids <- rbind(GF,GG,GH)

    AllGrids$subject <- Pt

    # Check
    if (plotgrid) {
        require('png')
        img <- png::readPNG(sprintf("C:/Users/mbmhscc4/GitHub/ECoG_Data_Prep/meta/electrodes/Map_Pt%d.png",Pt))
        plot(0, type='n', xlim=c(0,350), ylim=c(0,350), main=sprintf("Pt%d electrode grids",Pt),xlab='x',ylab='y')
        rasterImage(img, 0, 0, 350, 350)
        points(AllGrids$x,AllGrids$y,cex=2.2)
        text(AllGrids$x,AllGrids$y,labels = AllGrids$labels, cex=0.5)
    }

    return(AllGrids)
}
