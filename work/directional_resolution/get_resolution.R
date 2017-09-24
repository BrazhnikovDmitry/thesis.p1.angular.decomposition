source("ddec_solver.R")
library(oce)

load("dummy_mesh.rda")
mesh <- dummy_mesh

L <- 200e3
wavenum <- 2*pi/L

wavenum_tmp <- array(wavenum, c(dummy_mesh$nx, dummy_mesh$ny))
#gr <- generate_circ_antenna(mesh, 250, 250, wavenum_tmp)
lambda <- L
#gr <- generate_circ_antenna(mesh, 250, 250, wavenum_tmp, num_p = c(12, 24, 12, 6), radius = c(3*lambda/2, lambda / 2, lambda / 4, lambda / 8), shift_rad = c(-pi/5, 0, pi / 5, 2 * pi / 5))
gr <- generate_circ_antenna(mesh, 250, 250, wavenum_tmp, num_p = c(24, 12, 0), radius = c(1.2*lambda/2, 0.5*lambda / 2, lambda/2), shift_rad = c(-pi/5, 0, pi / 5))
#gr$x <- c(gr$x[2:length(gr$x)])
#gr$y <- c(gr$y[2:length(gr$y)])

# triangle
#Lgr <- L/mesh$dx[1]
#gr$x <- as.integer(250 + c(Lgr*cos(pi/6), -Lgr*cos(pi/6), 0))
#gr$y <- as.integer(250 + c(-Lgr*sin(pi/6), -Lgr*sin(pi/6), Lgr))
## cross
#gr$x <- as.integer(250 + c(Lgr*cos(pi/6), -Lgr*cos(pi/6), 0, 0))
#gr$y <- as.integer(250 + c(-Lgr*sin(pi/6), -Lgr*sin(pi/6), Lgr, 0))

#gr$x <- c(100, 200, 300, 250)
#gr$y <- c(100, 300, 200, 250)

np <- length(gr$x)
dist_m <- array(0, c(sum(1:np), 2))
ii <- 1

for (i in seq(np)) {
	for (j in i:length(np)) {
		dist_m[ii, 1] <- mesh$xx[gr$x[i]] - mesh$xx[gr$x[j]]
		dist_m[ii, 2] <- mesh$yy[gr$y[i]] - mesh$yy[gr$y[j]]
		ii <- ii + 1
	}
}

NWN <- 100
kseq <- seq(-2*wavenum, 2*wavenum, len = NWN)
lseq <- seq(-2*wavenum, 2*wavenum, len = NWN)

G <- array(0, dim = c(NWN, NWN))
B <- array(0, dim = c(NWN, NWN))
for (i in seq(NWN)) {
	for (j in seq(NWN)) {
		G[i, j] <- 1 + sum(2*cos(kseq[i]*dist_m[, 1] + lseq[j]*dist_m[, 2]))
		B[i, j] <- 1/np^2*sum( exp(1i*(kseq[i]*(mesh$xx[gr$x] - mesh$xx[250]) + lseq[j]*(mesh$yy[gr$y] - mesh$yy[250]))) )
	}
}

#source("/home/dmitry/Work/ittsunami/ittsunami.R")
m.png("../../figures/dir_resolution.png", w=14, h=5)
par(mfrow = c(1, 3))
imagep(mesh$xx, mesh$yy, mesh$H)
points(mesh$xx[gr$x], mesh$yy[gr$y])
lines(mesh$xx[250] + L/2*cos(seq(0, 2*pi, len = 100)), mesh$yy[250] + L/2*sin(seq(0, 2*pi, len = 100)), lwd = 2)
imagep(kseq, lseq, G)
lines(wavenum*cos(seq(0, 2*pi, len = 100)), wavenum*sin(seq(0, 2*pi, len = 100)), lwd = 2)
imagep(kseq, lseq, abs(B))
lines(wavenum*cos(seq(0, 2*pi, len = 100)), wavenum*sin(seq(0, 2*pi, len = 100)), lwd = 2)
dev.off()
