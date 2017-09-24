library("oce")

om <- 2 * pi / 12.4206
f0 <- -1.0312e-4
ii <- complex(real = 0, imaginary = 1)


gen_data_rot_wave <- function(mesh,
              alphas,
              wavelens,
              amps,
              phis,
              n2s = 0,
              i0 = 1,
              j0 = 1,
              om = 2 * pi / 12.4206 / 3600,
              f0 = -1.0312e-4)
{
     nx <- mesh$nx
     ny <- mesh$ny
     n_alp <- length(alphas)
     
     gg <- expand.grid(mesh$xx, mesh$yy)
     gg[, 1] <- gg[, 1] - mesh$xx[i0]
     gg[, 2] <- gg[, 2] - mesh$yy[j0]
     
     if (is.null(dim(wavelens))) {
          tmp <- array(0, c(nx, ny, n_alp))
          for (i in 1:n_alp)
               tmp[, , i] <- array(wavelens[i], c(nx, ny))
          wavelens <- tmp
     }
     kks <- 2 * pi / wavelens
     
     p <- array(0, dim = c(nx, ny))
     u <- array(0, dim = c(nx, ny))
     v <- array(0, dim = c(nx, ny))
     
     for (i in 1:n_alp) {
          if (n2s[1] != 0) {
               noiz_re <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
               noiz_im <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
          } else {
               noiz_re <- double(nx * ny)
               noiz_im <- noiz_re
          }
          
          re <-
               amps[i] * cos(kks[, , i] * cos(alphas[i]) * gg[, 1] + kks[, , i] * sin(alphas[i]) *
                                  gg[, 2] + phis[i]) + noiz_re
          im <-
               amps[i] * sin(kks[, , i] * cos(alphas[i]) * gg[, 1] + kks[, , i] * sin(alphas[i]) *
                                  gg[, 2] + phis[i]) + noiz_im
          p <- p + re + ii * im
          
          if (any(!is.finite(p))) {
               a <- 1
               print("gggg")
          }
          
          if (n2s[1] != 0) {
               noiz_re <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
               noiz_im <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
          } else {
               noiz_re <- double(nx * ny)
               noiz_im <- noiz_re
          }
          
          ang <- alphas[i]
          amp_cs <-
               kks[, , i] / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(ang) ^ 2
                                                                  + f0 ^ 2 * sin(ang) ^ 2)
          pha_u <- atan2(f0 * sin(ang), om * cos(ang))
          if (n2s[1] != 0) {
               noiz_re <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
               noiz_im <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
          } else {
               noiz_re <- double(nx * ny)
               noiz_im <- noiz_re
          }
          re <-
               amp_cs * (amps[i] * cos(kks[, , i] * cos(ang) * gg[, 1] + kks[, , i] *
                                            sin(ang) * gg[, 2] + phis[i] - pha_u) + noiz_re)
          im <-
               amp_cs * (amps[i] * sin(kks[, , i] * cos(ang) * gg[, 1] + kks[, , i] *
                                            sin(ang) * gg[, 2] + phis[i] - pha_u) + noiz_im)
          u <- u + re + ii * im
          
          amp_cs <-
               kks[, , i] / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(ang) ^ 2
                                                                 + f0 ^ 2 * cos(ang) ^ 2)
          pha_v <- atan2(-f0 * cos(ang), om * sin(ang))
          if (n2s[1] != 0) {
               noiz_re <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
               noiz_im <- sample(-1e6:1e6, nx * ny) / 1e6 * amps[i] * n2s[i]
          } else {
               noiz_re <- double(nx * ny)
               noiz_im <- noiz_re
          }
          re <-
               amp_cs * (amps[i] * cos(kks[, , i] * cos(ang) * gg[, 1] + kks[, , i] *
                                            sin(ang) * gg[, 2] + phis[i] - pha_v) + noiz_re)
          im <-
               amp_cs * (amps[i] * sin(kks[, , i] * cos(ang) * gg[, 1] + kks[, , i] *
                                            sin(ang) * gg[, 2] + phis[i] - pha_v) + noiz_im)
          v <- v + re + ii * im
     }
     
     fx <- 0.5 * Conj(u) * p
     fy <- 0.5 * Conj(v) * p
     
     return (list(
          p = p,
          u = u,
          v = v,
          fx = fx,
          fy = fy
     ))
}