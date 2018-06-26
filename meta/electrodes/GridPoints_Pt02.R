GridPoints_Pt02 <- function(plotgrid = FALSE) {
    # --- Pt2
    Pt <- 2

    # Grid A
    A1 <- c(120, 286)
    A5 <- c(122, 214)
    A16 <- c(171, 287)
    A20 <- c(173, 216)
    grid <- defineGrid(dims = c(4,5), A1, A5, A16, label='A')
    GA <- grid$x

    AllGrids <- rbind(GA)

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
