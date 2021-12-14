pkgname <- "GWASr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GWASr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GLambda")
### * GLambda

flush(stderr()); flush(stdout())

### Name: GLambda
### Title: GLambda
### Aliases: GLambda

### ** Examples

GLambda(ADHDMetaP$P)



cleanEx()
nameEx("plot_manhattan")
### * plot_manhattan

flush(stderr()); flush(stdout())

### Name: plot_manhattan
### Title: Manhattan Plot
### Aliases: plot_manhattan

### ** Examples

plot_manhattan(ADHDMeta)



cleanEx()
nameEx("plot_qq")
### * plot_qq

flush(stderr()); flush(stdout())

### Name: plot_qq
### Title: QQ Plot
### Aliases: plot_qq

### ** Examples

plot_qq(ADHDMetaP)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
