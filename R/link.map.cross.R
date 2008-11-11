`link.map.cross` <- 
function (object, chr, max.dist, marker.names = TRUE, horizontal = FALSE, 
    ...) 
{
    dots <- list(...)
    old.xpd <- par("xpd")
    par(xpd = TRUE)
    on.exit(par(xpd = old.xpd))
    map <- pull.map(object)
    if (!missing(chr)) {
        if (any(is.na(pmatch(chr, names(map))))) 
            stop("Some names of chromosome(s) subset do not match names of map.")
        map <- map[chr]
    }
    n.chr <- length(map)
    mt <- list()
    if (!missing(max.dist)) 
        map <- lapply(map, function(el, max.dist) el[el < max.dist], 
            max.dist)
    maxlen <- max(unlist(lapply(map, max)))
    if (!marker.names) {
        chrpos <- 1:n.chr
        thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
        if (!is.na(pmatch("cex", names(dots)))) 
            cex <- dots$cex
        else cex <- par("cex")
        chrpos <- seq(1, n.chr * 2, by = 2)
        thelim <- range(chrpos) + c(-0.35, 1.35)
        for (i in 1:n.chr) {
            mt[[i]] <- map[[i]]
            conv <- par("pin")[2]/maxlen
            for (j in 1:(length(mt[[i]]) - 1)) {
                ch <- mt[[i]][j + 1] * conv - (mt[[i]][j] * conv + 
                  10 * par("csi") * cex/9)
                if (ch < 0) {
                  temp <- mt[[i]][j + 1] * conv + abs(ch)
                  mt[[i]][j + 1] <- temp/conv
                }
            }
        }
        maxlen <- max(unlist(lapply(mt, max)))
    }
    if (horizontal) {
        plot(0, 0, type = "n", xlim = c(0, maxlen), ylim = rev(thelim), 
            yaxs = "i", xlab = "Location (cM)", ylab = "Chromosome", 
            yaxt = "n", ...)
        a <- par("usr")
        for (i in 1:n.chr) {
            if (marker.names) {
                text(mt[[i]], chrpos[i] + 0.5, names(map[[i]]), 
                  srt = 90, adj = c(1, 0.5), ...)
                segments(map[[i]], chrpos[i] + 0.25, map[[i]], 
                  chrpos[i] + 0.3)
                segments(map[[i]], chrpos[i] + 0.3, mt[[i]], 
                  chrpos[i] + 0.4)
                segments(mt[[i]], chrpos[i] + 0.4, mt[[i]], chrpos[i] + 
                  0.45)
            }
            segments(min(map[[i]]), chrpos[i] - 0.02, max(map[[i]]), 
                chrpos[i] - 0.02)
            segments(min(map[[i]]), chrpos[i] + 0.02, max(map[[i]]), 
                chrpos[i] + 0.02)
            segments(map[[i]], chrpos[i] - 0.2, map[[i]], chrpos[i] + 
                0.2)
        }
        for (i in seq(along = chrpos)) axis(side = 2, at = chrpos[i], 
            labels = names(map)[i])
    }
    else {
        plot(0, 0, type = "n", ylim = c(maxlen, 0), xlim = thelim, 
            xaxs = "i", ylab = "Location (cM)", xlab = "Chromosome", 
            axes = FALSE, ...)
        axis(side = 2, ylim = c(maxlen, 0))
        for (i in 1:n.chr) {
            if (marker.names) {
                text(chrpos[i] + 0.5, mt[[i]], names(map[[i]]), 
                  adj = c(0, 0.5), ...)
                segments(chrpos[i] + 0.25, map[[i]], chrpos[i] + 
                  0.3, map[[i]])
                segments(chrpos[i] + 0.3, map[[i]], chrpos[i] + 
                  0.4, mt[[i]])
                segments(chrpos[i] + 0.4, mt[[i]], chrpos[i] + 
                  0.45, mt[[i]])
            }
            segments(chrpos[i] - 0.02, min(map[[i]]), chrpos[i] - 
                0.02, max(map[[i]]), lwd = 1)
            segments(chrpos[i] + 0.02, min(map[[i]]), chrpos[i] + 
                0.02, max(map[[i]]), lwd = 1)
            segments(chrpos[i] - 0.2, map[[i]], chrpos[i] + 0.2, 
                map[[i]])
        }
        axis(side = 1, at = chrpos, labels = names(map))
    }
    if (is.na(pmatch("main", names(dots))) & !as.logical(sys.parent())) 
        title("Genetic Map")
    names(mt) <- names(map)
    invisible(list(mt = mt, map = map, chrpos = chrpos))
}

