### Processing first synthetic experiment
# Error in amplitude
err_amp <- array(0, c(NUM_trials, 4, 2))
err_amp[, , 1] <- get_fld(err1, "err_amp")
err_amp[, , 2] <- get_fld(err2, "err_amp")

par(mfrow = c(1, 2))
imagep(err_amp[,, 1])
imagep(err_amp[,, 2])

which (max(err_amp[,, 1]) == err_amp[,, 1], arr.ind = T)
which (max(err_amp[,, 2]) == err_amp[,, 2], arr.ind = T)

err_dir <- array(0, c(NUM_trials, 4, 2))
err_dir[, , 1] <- get_fld(err1, "err_dir")
err_dir[, , 2] <- get_fld(err2, "err_dir")

par(mfrow = c(1, 2))
imagep(err_dir[,, 1])
imagep(err_dir[,, 2])
