# Keep accidental default graphics output out of the repository during tests.
#
# Some report code opens explicit png/pdf devices, but a plotting call without
# an active explicit device can make R create Rplots.pdf in the working
# directory. A null default PDF device avoids local artifacts while preserving
# explicit output-file tests.
options(device = function(...) grDevices::pdf(file = NULL, ...))
