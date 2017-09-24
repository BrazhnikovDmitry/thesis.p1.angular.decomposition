require(limSolve)
require(oce)

crt_kernel_matrix_cross_flux_all_ener <- function(
               mesh,
              dt_pres,
              dt_u,
              dt_v,
              wavenum,
              alpha_s,
              seq_x,
              seq_y,
              plot = F,
              f0 = -1.0312e-4,
              om = 2 * pi / 12.4206 / 3600)
{
     la <- length(alpha_s)

     ns <- length(seq_x)
     ldt <- (ns) ^ 2
     if (ldt < 1) { print('FAIL! BAD INPUT DATA!'); return (NULL) }

     dtc_p <- rep(0, ns)
     dtc_u <- rep(0, ns)
     dtc_v <- rep(0, ns)
     xxc <- rep(0, ns)
     yyc <- rep(0, ns)
     kk <- rep(0, ns)
     l <- 0

     for (i in 1:ns) {
          l <- l + 1
          dtc_p[l] <- dt_pres[seq_x[i], seq_y[i]]
          dtc_u[l] <-  dt_u[seq_x[i], seq_y[i]]
          dtc_v[l] <-  dt_v[seq_x[i], seq_y[i]]
          xxc[l] <- mesh$xx[seq_x[i]]
          yyc[l] <- mesh$yy[seq_y[i]]
          kk[l] <- wavenum[seq_x[i], seq_y[i]]
     }

     ldt <- sum(1:ns)
     b <- rep(0, 10 * ldt)
     kk_x <- rep(0, ldt)
     kk_y <- rep(0, ldt)
     kk_vel <- rep(0, ldt)
     l <- 0

     coef_u_vel <- max(abs(dtc_p)) / max(abs(c(dtc_u, dtc_v)))
     coef_v_vel <- coef_u_vel

     dtc_u <- dtc_u * coef_u_vel
     dtc_v <- dtc_v * coef_v_vel

     for (i in 1:(ns)) for (j in i:(ns)) {
          l <- l + 1
          
          kk_vel[l] <- kk[j]
          kk_x[l] <- kk[i] * xxc[i] - kk[j] * xxc[j]
          kk_y[l] <- kk[i] * yyc[i] - kk[j] * yyc[j]
          
          b[l] <- Re(dtc_p[i] * Conj(dtc_u[j]))
          b[ldt + l] <- Im(dtc_p[i] * Conj(dtc_u[j]))
          b[2 * ldt + l] <- Re(dtc_p[i] * Conj(dtc_v[j]))
          b[3 * ldt + l] <- Im(dtc_p[i] * Conj(dtc_v[j]))
          
          b[4 * ldt + l] <- Re(dtc_p[i] * Conj(dtc_p[j]))
          b[5 * ldt + l] <- Im(dtc_p[i] * Conj(dtc_p[j]))
          b[6 * ldt + l] <- Re(dtc_u[i] * Conj(dtc_u[j]))
          b[7 * ldt + l] <- Im(dtc_u[i] * Conj(dtc_u[j]))
          b[8 * ldt + l] <- Re(dtc_v[i] * Conj(dtc_v[j]))
          b[9 * ldt + l] <- Im(dtc_v[i] * Conj(dtc_v[j]))
     }

     ### Generate linear model aka formulae
     A <- array(0, dim = c(10 * ldt, la))

     for (j in 1:la) {
          ang <- alpha_s[j]
          
          amp_cs <-
               kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(ang) ^ 2 + 
                                                                 f0 ^ 2 * sin(ang) ^ 2) * coef_u_vel
          pha_u <- atan2(f0 * sin(ang), om * cos(ang))
          A[1:ldt, j] <-
               amp_cs * cos(kk_x * cos(ang) + kk_y * sin(ang) - pha_u)
          A[(ldt + 1):(2 * ldt), j] <-
               amp_cs * sin(kk_x * cos(ang) + kk_y * sin(ang) - pha_u)
          
          amp_cs <-
               kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(ang) ^ 2 + 
                                                                 f0 ^ 2 * cos(ang) ^ 2) * coef_v_vel
          pha_v <- atan2(-f0 * cos(ang), om * sin(ang))
          A[(2 * ldt + 1):(3 * ldt), j] <-
               amp_cs * cos(kk_x * cos(ang) + kk_y * sin(ang) - pha_v)
          A[(3 * ldt + 1):(4 * ldt), j] <-
               amp_cs * sin(kk_x * cos(ang) + kk_y * sin(ang) - pha_v)
          
          ### Energy part
          ### "Potential" energy
          A[(4 * ldt + 1):(5 * ldt), j] <-
               cos(kk_x * cos(ang) + kk_y * sin(ang))
          A[(5 * ldt + 1):(6 * ldt), j] <-
               sin(kk_x * cos(ang) + kk_y * sin(ang))
          
          ### "kinetic", u part
          amp_cs <-
               kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(ang) ^ 2 + 
                                                                 f0 ^ 2 * sin(ang) ^ 2) * coef_u_vel
          A[(6 * ldt + 1):(7 * ldt), j] <-
               amp_cs ^ 2 * cos(kk_x * cos(ang) + kk_y * sin(ang))
          A[(7 * ldt + 1):(8 * ldt), j] <-
               amp_cs ^ 2 * sin(kk_x * cos(ang) + kk_y * sin(ang))
          
          ### "kinetic", v part
          amp_cs <-
               kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(ang) ^ 2 + 
                                                                 f0 ^ 2 * cos(ang) ^ 2) * coef_v_vel
          A[(8 * ldt + 1):(9 * ldt), j] <-
               amp_cs ^ 2 * cos(kk_x * cos(ang) + kk_y * sin(ang))
          A[(9 * ldt + 1):(10 * ldt), j] <-
               amp_cs ^ 2 * sin(kk_x * cos(ang) + kk_y * sin(ang))
     }

     return (list(
          A = A,
          dt = b,
          coef_u = coef_u_vel,
          coef_v = coef_v_vel
     ))
}

generate_circ_antenna <- function(
          mesh,
              i,
              j,
              wavenum,
              num_p = c(24, 12, 6),
              lambda =2 * pi / wavenum[i, j],
              radius = c(lambda / 2, lambda / 4, lambda / 8),
              shift_rad = c(0, pi / 5, 2 * pi / 5),
              WKB = T,
              CUT_OFF_DEPTH = 50)
{
     wavenum_t <- wavenum; wavenum_t[, ] <- wavenum[i, j]   # make matrix with uniform distribution of wavenumber
     
     num_circs <- length(radius)
     
     x0 <- i
     y0 <- j
     xx <- x0
     yy <- y0
     
     for (idx in 1:num_circs) {
          delta <- shift_rad[idx]
          theta <- delta + seq(0, 2 * pi, len = num_p[idx])
          
          xx <-
               c(xx, ceiling(x0 + radius[idx] / (mesh$dx[1]) * cos(theta)))
          yy <-
               c(yy, ceiling(y0 + radius[idx] / (mesh$dy[1]) * sin(theta)))
     }
     
     idx <- which (xx >= 1)
     idx <- idx[which(xx[idx] <= mesh$nx)]
     xx <- xx[idx]
     yy <- yy[idx]
     idx <- which (yy >= 1)
     idx <- idx[which(yy[idx] <= mesh$ny)]
     xx <- xx[idx]
     yy <- yy[idx]

     idx <- integer(0)
     for (k in 1:length(xx)) 
          if (mesh$H[xx[k], yy[k]] >= CUT_OFF_DEPTH) idx <- c(idx, k)
     xx <- xx[idx]
     yy <- yy[idx]
     
     if (WKB) {
          npa <- length(xx)
          dk_dx <- double(npa)
          for (idx in 2:npa) {
               tmp_dist <- sqrt((mesh$xx[xx[idx]] - mesh$xx[xx[1]]) ^ 2 + (mesh$yy[yy[idx]] - mesh$yy[yy[1]]) ^ 2)
               dk_dx[idx] <- abs(wavenum[xx[idx], yy[idx]] - wavenum[xx[1], yy[1]]) / tmp_dist
          }
          WKB_crit <- dk_dx * 2 * pi / (wavenum[xx[1], yy[1]]) ^ 2
          
          idx_wkb <- which (WKB_crit <= 1)
          xx <- xx[idx_wkb]
          yy <- yy[idx_wkb]
     }
     
     tmp <- unique(cbind(xx, yy))
     xx <- tmp[, 1]
     yy <- tmp[, 2]
     
     return (list(
          x = xx,
          y = yy,
          wavenum_t = wavenum_t
     ))
}

get_idx_spectra <- function(dirs_spectra, decomp_sect) {
     if (decomp_sect[2] <= 2 * pi) {
          idx <- which (decomp_sect[1] < dirs_spectra)
          idx <- idx[which (decomp_sect[2] >= dirs_spectra[idx])]
     } else {
          idx1 <- which (decomp_sect[1] < dirs_spectra)
          idx2 <- which (decomp_sect[2] - 2 * pi >= dirs_spectra)
          idx <- c(idx1, idx2)
     }
     
     return (idx)
}

decomp_reconst <- function(
              norm,
              Sf,
              ddec_model,
              dirs_spectra,
              decomp_sect,
              num_parms = 10,
              all = T,
              ldt = NULL)
{
     if (all == TRUE) {
          Fx <- double(length(decomp_sect))
          Fy <- double(length(decomp_sect))
     } else {
          Fx <- double(length(decomp_sect) - 1)
          Fy <- double(length(decomp_sect) - 1)
     }
     
     if (is.null(ldt))
          ldt <- dim(ddec_model$A)[1] / num_parms
     
     for (k in 1:(length(decomp_sect) - 1)) {
          idx_dir <- get_idx_spectra(dirs_spectra, decomp_sect[k:(k + 1)])

          ffx <- ddec_model$A[, idx_dir] %*% Sf[idx_dir]

          Fx[k] <- 0.5 * ffx[1] * norm / 1e3 / ddec_model$coef_u
          Fy[k] <- 0.5 * ffx[2 * ldt + 1] * norm / 1e3 / ddec_model$coef_v
          
          if (all == T) {
               k <- length(decomp_sect)
               
               ffx <- ddec_model$A %*% Sf
               Fx[k] <- 0.5 * ffx[1] * norm / 1e3 / ddec_model$coef_u
               Fy[k] <- 0.5 * ffx[2 * ldt + 1] * norm / 1e3 / ddec_model$coef_v
          }
     }
     
     return (list(Fx = Fx, Fy = Fy))
}

get_p_phase <- function(
          mesh,
              dt_pres,
              dt_u,
              dt_v,
              wavenum,
              F_p,
              F_n,
              norms,
              seq_x,
              seq_y,
              plot = F,
              f0 = -1.0312e-4,
              om = 2 * pi / 12.4206 / 3600)
{
     ns <- length(seq_x)
     ldt <- ns
     
     dtc_p <- rep(0, ns)
     dtc_u <- rep(0, ns)
     dtc_v <- rep(0, ns)
     xxc <- rep(0, ns)
     yyc <- rep(0, ns)
     kk_vel <- rep(0, ns)
     l <- 0

     for (i in 1:ns) {
          l <- l + 1
          dtc_p[l] <- dt_pres[seq_x[i], seq_y[i]]
          dtc_u[l] <-  dt_u[seq_x[i], seq_y[i]]
          dtc_v[l] <-  dt_v[seq_x[i], seq_y[i]]
          xxc[l] <- mesh$xx[seq_x[i]]
          yyc[l] <- mesh$yy[seq_y[i]]
          kk_vel[l] <- wavenum[seq_x[i], seq_y[i]]
     }

     coef_u_vel <- max(abs(dtc_p)) / max(abs(c(dtc_u, dtc_v)))
     coef_v_vel <- coef_u_vel
     
     dtc_u <- dtc_u * coef_u_vel
     dtc_v <- dtc_v * coef_v_vel
     
     dt <- double(6 * ns)

     l <- 0
     
     for (i in 1:(ns)) {
          l <- l + 1
          
          dt[l] <- Re(dtc_p[i])
          dt[ldt + l] <- Re(dtc_u[i])
          dt[2 * ldt + l] <- Re(dtc_v[i])
          
          dt[3 * ldt + l] <- Im(dtc_p[i])
          dt[4 * ldt + l] <- Im(dtc_u[i])
          dt[5 * ldt + l] <- Im(dtc_v[i])
     }
     
     angp <- atan2(F_p$Fy, F_p$Fx)
     kp <- cbind(kk_vel * cos(angp), kk_vel * sin(angp))
     
     angn <- atan2(F_n$Fy, F_n$Fx)
     kn <- cbind(kk_vel * cos(angn), kk_vel * sin(angn))

     amp_cs_up <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(angp) ^ 2 + 
                                                            f0 ^ 2 * sin(angp) ^ 2)
     pha_up <- atan2(f0 * sin(angp), om * cos(angp))
     amp_cs_un <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(angn) ^ 2 + 
                                                            f0 ^ 2 * sin(angn) ^ 2)
     pha_un <- atan2(f0 * sin(angn), om * cos(angn))
     
     amp_cs_vp <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(angp) ^ 2 + 
                                                            f0 ^ 2 * cos(angp) ^ 2)
     pha_vp <- atan2(-f0 * cos(angp), om * sin(angp))
     amp_cs_vn <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(angn) ^ 2 + 
                                                            f0 ^ 2 * cos(angn) ^ 2)
     pha_vn <- atan2(-f0 * cos(angn), om * sin(angn))
     
     Ap_u <- sqrt(abs(2 * F_p$Fx * 1e3 / norms / amp_cs_up / cos(pha_up)))
     An_u <- sqrt(abs(2 * F_n$Fx * 1e3 / norms / amp_cs_un / cos(pha_un)))
     
     Ap_v <- sqrt(abs(2 * F_p$Fy * 1e3 / norms / amp_cs_vp / cos(pha_vp)))
     An_v <- sqrt(abs(2 * F_n$Fy * 1e3 / norms / amp_cs_vn / cos(pha_vn)))
     
     A <- array(0, dim = c(6 * ldt, 4))

     for (i in 1:ldt) {
          A[1:ldt, 1] <- Ap_u * cos(kp[, 1] * xxc + kp[, 2] * yyc)
          A[1:ldt, 2] <- -Ap_u * sin(kp[, 1] * xxc + kp[, 2] * yyc)
          A[1:ldt, 3] <- An_u * cos(kn[, 1] * xxc + kn[, 2] * yyc)
          A[1:ldt, 4] <- -An_u * sin(kn[, 1] * xxc + kn[, 2] * yyc)
          
          A[(ldt + 1):(2 * ldt), 1] <-
               amp_cs_up * Ap_u * cos(kp[, 1] * xxc + kp[, 2] * yyc - pha_up) * coef_u_vel
          A[(ldt + 1):(2 * ldt), 2] <-
               -amp_cs_up * Ap_u * sin(kp[, 1] * xxc + kp[, 2] * yyc - pha_up) * coef_u_vel
          A[(ldt + 1):(2 * ldt), 3] <-
               amp_cs_un * An_u * cos(kn[, 1] * xxc + kn[, 2] * yyc - pha_un) * coef_u_vel
          A[(ldt + 1):(2 * ldt), 4] <-
               -amp_cs_un * An_u * sin(kn[, 1] * xxc + kn[, 2] * yyc - pha_un) * coef_u_vel
          
          A[(2 * ldt + 1):(3 * ldt), 1] <-
               amp_cs_vp * Ap_v * cos(kp[, 1] * xxc + kp[, 2] * yyc - pha_vp) * coef_v_vel
          A[(2 * ldt + 1):(3 * ldt), 2] <-
               -amp_cs_vp * Ap_v * sin(kp[, 1] * xxc + kp[, 2] * yyc - pha_vp) * coef_v_vel
          A[(2 * ldt + 1):(3 * ldt), 3] <-
               amp_cs_vn * An_v * cos(kn[, 1] * xxc + kn[, 2] * yyc - pha_vn) * coef_v_vel
          A[(2 * ldt + 1):(3 * ldt), 4] <-
               -amp_cs_vn * An_v * sin(kn[, 1] * xxc + kn[, 2] * yyc - pha_vn) * coef_v_vel
          
          A[(3 * ldt + 1):(4 * ldt), 1] <-
               Ap_u * sin(kp[, 1] * xxc + kp[, 2] * yyc)
          A[(3 * ldt + 1):(4 * ldt), 2] <-
               Ap_u * cos(kp[, 1] * xxc + kp[, 2] * yyc)
          A[(3 * ldt + 1):(4 * ldt), 3] <-
               An_u * sin(kn[, 1] * xxc + kn[, 2] * yyc)
          A[(3 * ldt + 1):(4 * ldt), 4] <-
               An_u * cos(kn[, 1] * xxc + kn[, 2] * yyc)
          
          A[(4 * ldt + 1):(5 * ldt), 1] <-
               amp_cs_up * Ap_u * sin(kp[, 1] * xxc + kp[, 2] * yyc - pha_up) * coef_u_vel
          A[(4 * ldt + 1):(5 * ldt), 2] <-
               amp_cs_up * Ap_u * cos(kp[, 1] * xxc + kp[, 2] * yyc - pha_up) * coef_u_vel
          A[(4 * ldt + 1):(5 * ldt), 3] <-
               amp_cs_un * An_u * sin(kn[, 1] * xxc + kn[, 2] * yyc - pha_un) * coef_u_vel
          A[(4 * ldt + 1):(5 * ldt), 4] <-
               amp_cs_un * An_u * cos(kn[, 1] * xxc + kn[, 2] * yyc - pha_un) * coef_u_vel
          
          A[(5 * ldt + 1):(6 * ldt), 1] <-
               amp_cs_vp * Ap_v * sin(kp[, 1] * xxc + kp[, 2] * yyc - pha_vp) * coef_v_vel
          A[(5 * ldt + 1):(6 * ldt), 2] <-
               amp_cs_vp * Ap_v * cos(kp[, 1] * xxc + kp[, 2] * yyc - pha_vp) * coef_v_vel
          A[(5 * ldt + 1):(6 * ldt), 3] <-
               amp_cs_vn * An_v * sin(kn[, 1] * xxc + kn[, 2] * yyc - pha_vn) * coef_v_vel
          A[(5 * ldt + 1):(6 * ldt), 4] <-
               amp_cs_vn * An_v * cos(kn[, 1] * xxc + kn[, 2] * yyc - pha_vn) * coef_v_vel
     }
     
     return (list(
          dt = dt,
          A = A,
          coef_u = coef_u_vel,
          coef_v = coef_v_vel
     ))
}

get_only_p_phase <- function(
                         mesh,
                             dt_pres,
                             dt_u,
                             dt_v,
                             wavenum,
                             F_p,
                             F_n,
                             norms,
                             seq_x,
                             seq_y,
                             plot = F,
                             f0 = -1.0312e-4,
                             om = 2 * pi / 12.4206 / 3600)
{
     ns <- length(seq_x)
     ldt <- ns
     
     dtc_p <- rep(0, ns)
     dtc_u <- rep(0, ns)
     dtc_v <- rep(0, ns)
     xxc <- rep(0, ns)
     yyc <- rep(0, ns)
     kk <- rep(0, ns)
     l <- 0
     
     for (i in 1:ns) {
          l <- l + 1
          dtc_p[l] <- dt_pres[seq_x[i], seq_y[i]]
          dtc_u[l] <-  dt_u[seq_x[i], seq_y[i]]
          dtc_v[l] <-  dt_v[seq_x[i], seq_y[i]]
          xxc[l] <- mesh$xx[seq_x[i]]
          yyc[l] <- mesh$yy[seq_y[i]]
          kk[l] <- wavenum[seq_x[i], seq_y[i]]
     }
     
     coef_u_vel <- max(abs(dtc_p)) / max(abs(c(dtc_u, dtc_v)))
     coef_v_vel <- coef_u_vel
     
     dtc_u <- dtc_u * coef_u_vel
     dtc_v <- dtc_v * coef_v_vel
     
     dt <- double(2 * ns)
     kk_vel <- kk
     
     l <- 0
     
     for (i in 1:(ns)) {
          l <- l + 1
          
          dt[l] <- Re(dtc_p[i])
          dt[ldt + l] <- Im(dtc_p[i])
     }
     
     angp <- atan2(F_p$Fy, F_p$Fx)
     kp <- cbind(kk_vel * cos(angp), kk_vel * sin(angp))
     
     angn <- atan2(F_n$Fy, F_n$Fx)
     kn <- cbind(kk_vel * cos(angn), kk_vel * sin(angn))

     amp_cs_up <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(angp) ^ 2 + 
                                                            f0 ^ 2 * sin(angp) ^ 2)
     pha_up <- atan2(f0 * sin(angp), om * cos(angp))
     amp_cs_un <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * cos(angn) ^ 2 + 
                                                            f0 ^ 2 * sin(angn) ^ 2)
     pha_un <- atan2(f0 * sin(angn), om * cos(angn))
     
     amp_cs_vp <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(angp) ^ 2 + 
                                                            f0 ^ 2 * cos(angp) ^ 2)
     pha_vp <- atan2(-f0 * cos(angp), om * sin(angp))
     amp_cs_vn <-
          kk_vel / (om ^ 2 - f0 ^ 2) / mesh$rhoConst * sqrt(om ^ 2 * sin(angn) ^ 2 + 
                                                            f0 ^ 2 * cos(angn) ^ 2)
     pha_vn <- atan2(-f0 * cos(angn), om * sin(angn))
     
     Ap_u <- sqrt(abs(2 * F_p$Fx * 1e3 / norms / amp_cs_up / cos(pha_up)))
     An_u <- sqrt(abs(2 * F_n$Fx * 1e3 / norms / amp_cs_un / cos(pha_un)))
     
     Ap_v <- sqrt(abs(2 * F_p$Fy * 1e3 / norms / amp_cs_vp / cos(pha_vp)))
     An_v <- sqrt(abs(2 * F_n$Fy * 1e3 / norms / amp_cs_vn / cos(pha_vn)))
     
     Ap <- mean(c(Ap_u, Ap_v))
     An <- mean(c(An_u, An_v))
     
     A <- array(0, dim = c(2 * ldt, 4))

     for (i in 1:ldt) {
          A[1:ldt, 1] <- Ap * cos(kp[, 1] * xxc + kp[, 2] * yyc)
          A[1:ldt, 2] <- -Ap * sin(kp[, 1] * xxc + kp[, 2] * yyc)
          A[1:ldt, 3] <- An * cos(kn[, 1] * xxc + kn[, 2] * yyc)
          A[1:ldt, 4] <- -An * sin(kn[, 1] * xxc + kn[, 2] * yyc)
          
          A[(ldt + 1):(2 * ldt), 1] <-
               Ap * sin(kp[, 1] * xxc + kp[, 2] * yyc)
          A[(ldt + 1):(2 * ldt), 2] <-
               Ap * cos(kp[, 1] * xxc + kp[, 2] * yyc)
          A[(ldt + 1):(2 * ldt), 3] <-
               An * sin(kn[, 1] * xxc + kn[, 2] * yyc)
          A[(ldt + 1):(2 * ldt), 4] <-
               An * cos(kn[, 1] * xxc + kn[, 2] * yyc)
     }
     
     return (list(
          dt = dt,
          A = A,
          coef_u = coef_u_vel,
          coef_v = coef_v_vel
     ))
}

crt_f_plane <- function(mesh,
                        i0 = NULL,
                        j0 = NULL,
                        ijseq = NULL,
                        ijlims = NULL,
                        xydist = NULL,
                        average = TRUE)
{
     if (!is.null(ijseq)) {
          mesh_f_pl <- crt_f_plane_mesh_ijseq(mesh, i0, j0, ijseq, average = average)
     }

     if (!is.null(ijlims)) {
          i0 <- as.integer(mean(ijlims$i))
          j0 <- as.integer(mean(ijlims$j))
          di <- max(abs(ijlims$i - i0))
          dj <- max(abs(ijlims$j - j0))
          ijseq <- list(di = c(-1, 1)*di, dj = c(-1, 1)*dj)
          mesh_f_pl <- crt_f_plane_mesh_ijseq(mesh, i0, j0, ijseq, average = average)
     }
     
     if (!is.null(xydist)) {
          # SO FAR NOTHING
     }
     return (mesh_f_pl)
}

crt_f_plane_mesh_ijseq <- function(mesh,
                                   i0,
                                   j0,
                                   ijseq = NULL,
                                   xydist = NULL,
                                   average = TRUE)
{
     mesh_f_pl <- list(NULL)
     ### Only ijseq
     lat0 <- mesh$lat[i0, j0]
     lon0 <- mesh$lon[i0, j0]
     
     mesh_f_pl$lat0 <- lat0; mesh_f_pl$lon0 <- lon0
     
     iseq <- max(c(1, i0 + ijseq$di[1])):min(c(mesh$nx, i0 + ijseq$di[2]))
     jseq <- max(c(1, j0 + ijseq$dj[1])):min(c(mesh$ny, j0 + ijseq$dj[2]))
     
     mesh_f_pl$H <- mesh$H[iseq, jseq]
     
     mesh_f_pl$iseq <- iseq
     mesh_f_pl$jseq <- jseq
     mesh_f_pl$nx <- length(iseq)
     mesh_f_pl$ny <- length(jseq)
     
     mesh_f_pl$lats <- mesh$lat[iseq, jseq]
     mesh_f_pl$lons <- mesh$lon[iseq, jseq]

     mesh_f_pl$f <- coriolis(mesh_f_pl$lats)
     mesh_f_pl$rhoConst <- mesh$rhoConst
     
     xy_gr <- geodXy(mesh_f_pl$lons, mesh_f_pl$lats, mesh_f_pl$lon0, mesh_f_pl$lat0)
     mesh_f_pl$xx <- array(xy_gr[,1], dim=c(mesh_f_pl$nx, mesh_f_pl$ny))
     mesh_f_pl$yy <- array(xy_gr[,2], dim=c(mesh_f_pl$nx, mesh_f_pl$ny))
     
     mesh_f_pl$dx <- apply(mesh_f_pl$xx, 2, diff)
     mesh_f_pl$dy <- apply(mesh_f_pl$yy, 1, diff)
     
     mesh_f_pl$aver <- average
     if (average) {
          mesh_f_pl$aver <- TRUE
          mesh_f_pl$xx <- apply(mesh_f_pl$xx, 1, mean)
          mesh_f_pl$yy <- apply(mesh_f_pl$yy, 2, mean)
          
          mesh_f_pl$dx <- diff(mesh_f_pl$xx) #apply(mesh_f_pl$dx, 1, mean)
          mesh_f_pl$dy <- diff(mesh_f_pl$yy) #apply(mesh_f_pl$dy, 1, mean) # note the dimension!
     }
     
     dd_tmp <- sqrt( (mesh_f_pl$lats - mesh_f_pl$lat0)^2 + (mesh_f_pl$lons - mesh_f_pl$lon0)^2 )
     mesh_f_pl$i0 <- which (dd_tmp == min(dd_tmp), arr.ind=T)[1]
     mesh_f_pl$j0 <- which (dd_tmp == min(dd_tmp), arr.ind=T)[2]
     
     return (mesh_f_pl)
}
