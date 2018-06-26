defineGrid <- function(dims,lowerLeft,upperLeft,lowerRight=NULL,label=NULL,imgSize=c(350,350)) {
    # Note: rotation is around the lowerLeft point.
    require('spdep')
    deg2rad <- function(deg) {(deg * pi) / (180)}
    if (is.null(lowerRight)) {
        lowerRight <- lowerLeft
    }
    lowerLeft[2] <- imgSize[2] - lowerLeft[2]
    upperLeft[2] <- imgSize[2] - upperLeft[2]
    lowerRight[2] <- imgSize[2] - lowerRight[2]

    d <- expand.grid(y=0:(dims[2]-1),x=0:(dims[1]-1))
    # Scale
    h <- norm(lowerLeft-upperLeft, type = '2')
    w <- norm(lowerRight-lowerLeft, type = '2')
    if (max(abs(d$x)) > 0) {
        d$x <- d$x * w/max(abs(d$x))
    }
    if (max(abs(d$y)) > 0) {
        d$y <- d$y * h/max(abs(d$y))
    }

    # Rotate
    x <- c(upperLeft-lowerLeft)
    t <- atan2(x[2],x[1])
    if (t > 0) {
        rad <- t - (pi/2)
    } else {
        rad <- t + (pi/2)
    }
    dr <- as.data.frame(spdep::Rotation(d[,c('x','y')], rad))
    names(dr) <- c('x','y')

    # Offset
    dro <- dr
    xo <- lowerLeft[1]# - dr$x[1] (this will always be zero)
    dro$x <- dr$x + xo
    yo <- lowerLeft[2]# - dr$y[1] (this will always be zero)
    dro$y <- dr$y + yo

    if (is.null(label)) {
        dro$index <- 1:prod(dims)
    } else {
        dro$grid <- label
        dro$index <- 1:prod(dims)
        dro$labels <- paste(label, 1:prod(dims), sep='')
    }
    grid <- list(
        x = dro,
        angle = rad
    )
    return(grid)
}
