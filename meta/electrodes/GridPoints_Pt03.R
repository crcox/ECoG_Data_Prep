GridPoints_Pt03 <- function(plotgrid = FALSE) {
    # --- Pt3
    Pt <- 3

    # Grid A
    A1 <- c(182, 258)
    A4 <- c(153, 212)
    grid <- defineGrid(dims = c(1,4), A1, A4, label='A')
    GA <- grid$x

    # Grid B
    B1 <- c(198, 254)
    B5 <- c(161, 191)
    B16 <- c(243, 228)
    B20 <- c(206, 166)
    grid <- defineGrid(dims = c(4,5), B1, B5, B16, label='B')
    GB <- grid$x

    AllGrids <- rbind(GA,GB)

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
