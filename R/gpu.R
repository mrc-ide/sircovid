##' Create a model for the GPU. This requires a working nvcc toolchain
##' and GPU device to work properly. Note that this will cache the
##' compilation within a session, so if you want to change GPU options
##' you will need to restart R. However, if you change the *model*
##' (e.g., by changing the result of options like `rewrite_dims`) the
##' model will recompile. Compilation speed is slow enough that the
##' cachine will be obvious.
##'
##' @title Create GPU model
##'
##' @param model The name of a sircovid model, either `carehomes` or
##'   `basic`
##'
##' @param ... Additional arguments passed to [odin.dust::odin_dust_]
##'   and from there into either odin's options or the gpu options.
##'
##' @param real_t The type to use for floating point numbers. The
##'   default is "float" which differs from the CPU model but is
##'   typically what you want on the GPU.
##'
##' @param gpu The argument passed to [odin.dust::odin_dust_] as
##'   `gpu`. The default here is `TRUE`, but you may want to pass the
##'   results of [dust::dust_cuda_options] in order to control
##'   compilation.
##'
##' @return A model generator, like [sircovid::carehomes] that can run
##'   on the GPU available to this computer.
##'
##' @export
compile_gpu <- function(model = "carehomes", ..., real_t = "float",
                        gpu = TRUE) {
  path <- sircovid_file(sprintf("odin/%s.R", model))
  odin.dust::odin_dust_(path, real_t = real_t, gpu = gpu, ...)
}
