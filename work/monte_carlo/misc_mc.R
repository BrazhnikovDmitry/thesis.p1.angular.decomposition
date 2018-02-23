run_solver <- function(antenna, rm = TRUE, poten_constr = T) {
     ddec_model <-
          crt_kernel_matrix_cross_flux_all_ener(
               mesh,
               p,
               u,
               v,
               antenna$wavenum_t,
               dirs_spectra,
               antenna$x,
               antenna$y,
               f0 = f0,#mesh_f_pl$f[mesh_f_pl$i0, mesh_f_pl$j0],
               om = om,#mesh$om
               rm_zeros = rm
          )
     
     norm <- 1#max(abs(ddec_model$dt))
     A <- ddec_model$A; dt <- ddec_model$dt/norm
     # P <- A[ddec_model$pendx, ]
     # Pdt <- dt[ddec_model$pendx]
     G <- array(0, dim=c(dim(A)[2], dim(A)[2]))
     diag(G) <- 1
     M <- list(
          X = A,
          p = rep(0.1, num_dir),
          off = 0,
          S = list(diag(1, num_dir, num_dir)),
          Ain = diag(1, num_dir, num_dir),#matrix(0, num_dir, num_dir),
          bin = rep(0, num_dir),
          C = matrix(0, 0, 0),
          sp = 1,
          y = dt,
          w = dt[] * 0 + 1
     )
     tryCatch({
          # Sf <- lsei(A, dt, NULL, NULL, G, rep(0, dim(A)[2]), type = 2)$X * norm
          Sf <- pcls(M)
          # if (!poten_constr) lsei(A, dt, NULL, NULL, G, rep(0, dim(A)[2]), type = 2)$X * norm
          #      else lsei(A, dt, P, Pdt, G, rep(0, dim(A)[2]), type = 1)$X * norm
     }, error = function(e) {
          sink()
          write(
               paste0("ERROR : point ", 1, 1, " message: ", conditionMessage(e), "\n"),
               file = "errors.txt",
               append = T
          )
     })
     
     return (Sf)
}

### Two errors: 1) energy contained in a bin \pm dang + dir scattered wave
# mean direction in such bin
calc_errors_1 <- function(angs, amps, dang, dirs_spectra, Sf) {
     nang <- length(angs)

     sf_ener <- double(nang)
     sf_dir <- double(nang)
     err_amp <- double(nang)
     err_dir <- double(nang)
     total_sf <- sum(Sf)
     err_tot <- (total_sf - sum(amps^2))/sum(amps^2)

     for (i in seq(nang)) {
          ang_syn <- angs[i]
          ang_syn <- if (ang_syn < 0) ang_syn + 2*pi else ang_syn
          lim_ang <- ang_syn + c(-1, 1)*dang
          idx_sf <- get_idx_spectra(dirs_spectra, lim_ang)
          sf_ener[i] <- sum(Sf[idx_sf])
          sf_dir[i] <- weighted.mean(dirs_spectra[idx_sf], Sf[idx_sf])
          err_amp[i] <- abs((sum(Sf[idx_sf]) - amps[i]^2)/amps[i]^2)
          err_dir[i] <- abs((sf_dir[i] - ang_syn))
     }
     
     return (list(sf_ener = sf_ener, sf_dir = sf_dir, err_amp = err_amp, err_dir = err_dir, err_tot = err_tot))
}

get_fld <- function(dt, fld) {
     out <- array(0, c(length(dt), length(dt[[1]][[fld]])))
     
     for (i in seq(length(dt))) out[i, ] <- dt[[i]][[fld]]
     
     out
}