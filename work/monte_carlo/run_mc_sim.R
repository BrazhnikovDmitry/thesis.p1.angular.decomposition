### The first synthetic experiment is representative of scattering problem
# Four waves: incident (0), reflected (pi), two side lobes (\pm pi/2)
# no phase variation
# with and without noise
###
source("/home/dmitry/Work/ittsunami/ittsunami.R")
source("ddec_solver.R")
source("generate_data.r")

load("dummy_mesh.rda")
mesh <- dummy_mesh
# Physics
wavelen <- 200e3
lambda <- wavelen
om <- 2*pi/(20*60)
f0 <- 0

NUM_trials <- 100
NUM_waves <- 6

# alphas <- seq(0, pi, pi/2, -pi/2, len = NUM_waves)
alphas <- c(0, pi/4, 5*pi/3)
# amp_c <- c(1, 0.5, 0.75)
amp_c <- c(1, 0.75, 0.5, 1.2)
nzs <- 0.0#.05

### Setting up the method
num_dir <- 100
dirs_spectra <- seq(pi/num_dir, pi/2 - pi/num_dir, len=25)
dirs_spectra <- c(dirs_spectra, pi/2 + dirs_spectra, pi + dirs_spectra, 3*pi/2 + dirs_spectra)
# alphas <- dirs_spectra#seq(0, pi, len = NUM_waves)

IP <- as.integer(mesh$nx/2); JP <- as.integer(mesh$ny/2)

Sf_out1 <- array(list(NULL), c(NUM_trials, 1))
Sf_out2 <- array(list(NULL), c(NUM_trials, 1))
err1 <- array(list(NULL), c(NUM_trials, 1))
err2 <- array(list(NULL), c(NUM_trials, 1))

all_amps <- array(0, dim=c(NUM_trials, NUM_waves))
all_phases <- array(0, dim=c(NUM_trials, NUM_waves))

dang <- 10*pi/180

for (i in 1:NUM_trials) {
     ### Generate amps and phases
     # angs <- c(alphas, runif(NUM_waves - length(alphas), min = 0, max = 1)*2*pi)#rep(1, NUM_waves)#
     # amps <- c(amp_c, runif(NUM_waves - length(alphas), min = 0, max = nzs))#rep(1, NUM_waves)#
     # phases <- c(0, 0, 0, runif(NUM_waves - length(alphas), min = 0, max = 2*pi))#rep(0, NUM_waves)#
     angs <- c(alphas, seq(0, 2*pi, len = NUM_trials)[i])
     amps <- amp_c
     phases <- c(0, 0, 0, 0)
     all_phases[i, ] <- phases
     all_amps[i, ] <- amps
     print(amps)
     print(phases)
     tmp_flds <-
          gen_data_rot_wave(mesh,
                   angs,
                   rep(wavelen, NUM_waves),
                   amps,
                   phases,
                   rep(nzs, NUM_waves),
                   om = om,
                   f0 = f0,
                    i0 = IP, j0 = JP)
     p <- tmp_flds$p
     u <- tmp_flds$u
     v <- tmp_flds$v
     
     
     wavenum <- array(2 * pi / wavelen, c(mesh$nx, mesh$nx))
     
     # Run decomposition
     antenna1 <- generate_circ_antenna(mesh, IP, JP, wavenum)
     antenna2 <-
          generate_circ_antenna(
               mesh,
               250,
               250,
               wavenum,
               num_p = c(24, 24, 12, 0),
               radius = c(3*lambda/2, 1.2 * lambda / 2, 0.5 * lambda / 2, lambda / 2),
               shift_rad = c(0, -pi / 5, 0, pi / 5)
          )

     # imagep(Re(p))
     # points(antenna1)
     # points(antenna2, col = "green")

     Sf_out1[[i]]$Sf <- run_solver(antenna1)
     Sf_out2[[i]]$Sf <- run_solver(antenna2)
     # calc errors
     ### Two errors: 1) energy contained in a bin \pm dang + dir scattered wave
     # mean direction in such bin
     err1[[i]] <- calc_errors_1(alphas, amps, dang, dirs_spectra, Sf_out1[[i]]$Sf)
     err2[[i]] <- calc_errors_1(alphas, amps, dang, dirs_spectra, Sf_out2[[i]]$Sf)

     m.png(paste0("dump_pix_mc/", i, ".png"), w = 15, h = 8)
     par(mfrow = c(1, 1))
     plot(dirs_spectra, Sf_out1[[i]]$Sf, ylim = c(0, max(c(Sf_out1[[i]]$Sf, Sf_out2[[i]]$Sf, amps^2))), type = "b")
     lines(dirs_spectra, Sf_out2[[i]]$Sf, col = "blue", pch = 21, bg = "blue", cex = 0.75, type = "b")
     points(angs, amps^2, col = "red")
     abline(v = alphas)
     dev.off()
     # plot(A%*%Sf_out[[i]]$Sf)
     # lines(dt*norm, col = "red")
     # lines(A%*%c(rep(1, 100)), col = "blue")
}
#save.image(file="three_waves.dat")

