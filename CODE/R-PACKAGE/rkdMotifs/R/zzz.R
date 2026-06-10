# R/zzz.R
# Explicitly loads libgeco2.so from inst/libs/ before the Rcpp wrapper .so
# so all geco2 symbols are resolved at package load time.
#
# We try both the plain name (libgeco2.so) and the versioned soname
# (libgeco2.so.1) because the embedded soname depends on whether CMake's
# VERSION/SOVERSION was set during the build.

.onLoad <- function(libname, pkgname) {
    ext <- .Platform$dynlib.ext

    candidates <- c(
        system.file("libs", paste0("libgeco2", ext),   package = pkgname),
        system.file("libs", paste0("libgeco2", ext, ".1"), package = pkgname)
    )

    for (so in candidates) {
        if (nzchar(so) && file.exists(so)) {
            dyn.load(so, local = FALSE, now = TRUE)
            return(invisible(NULL))
        }
    }

    warning("rkdMotifs2: libgeco2.so not found in inst/libs/ — ",
            "package may fail to load")
}
