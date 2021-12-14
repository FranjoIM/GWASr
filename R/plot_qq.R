##' \emph{plot_qq} takes GWAS summary statistics p-values and creates a ggplot based QQ plot.
##'
##' Basic functionality: \emph{plot_qq} first takes the p-values from specified data frame and then based on sorted p-values creates expectations. Subsequently the data are pruned by observed p-values in a bin-wise manner before getting plotted. This pruning is done to reduce the likelihood of crashing as usual GWAS meta-analyses have large number of datapoints (millions).
##'
##' @title QQ Plot
##' @param data Data frame which contains GWAS summary statistics (REQUIRED)
##' @param showN Shows total number of SNPs on the plot, default is TRUE
##' @param showGL Shows genomic constant lambda on the plot, default is TRUE
##' @param P A column name which contains P values, default is "P", column should be numeric.
##' @param xlab Optional x-axis label, default is NULL.
##' @param ylab Optional y-axis label, default is "-log10(P)"
##' @param main Optional plot title, default is "Manhattan plot"
##' @param pcol String indicating point colors, default is set by the package.
##' @param lcol String indicating reference line color, default is set by the package.
##' @return ggplot object
##' @author Franjo, Mingjing, and Xiaoxiao
##' @export
##' @examples
##' plot_qq(ADHDMetaP)
##' @importFrom dplyr %>%

plot_qq <- function(data=data, showN=TRUE, showGL=TRUE, P="P",
                    xlab=NULL, ylab=NULL, main=NULL, pcol=NULL, lcol=NULL){

  # If dataframe is not provided stop function and show  error
  if(base::is.null(data)){
    base::stop("Data frame is not specified.")
  }

  # Prepreprocessing step
  data <- data %>%
    # Rename variables
    dplyr::rename(Pval={{P}}) %>%
    # Select P only
    dplyr::select(Pval)

  # If p values aren't numeric
  if(base::class(data$Pval) != "numeric"){
    base::stop("P values are not class numeric.")
  }

  # If not provided, define y label
  if(base::is.null(ylab)){
    ylab=base::bquote("Observed"~"-log"[10]~"(P)")}

  # If not provided, define x label
  if(base::is.null(xlab)){
    xlab=base::bquote("Expected"~"-log"[10]~"(P)")}

  # If not provided, define graph title
  if(base::is.null(main)){
    main="QQ Plot"}

  # If not provided, define point color
  if(base::is.null(pcol)){
    pcol="#4092C9"}

  # If not provided, define line color
  if(base::is.null(lcol)){
    lcol="#FF8785"}

  # Calculate genomic constant lambda and number of observations
  if(showGL == TRUE) {
    GLam <- GWASr::GLambda(data$Pval)
    subGL <- base::bquote("\u03BB ="~.(base::round(GLam, digits=3))~"     ")}
  else {subGL <- NULL}

  if(showN == TRUE) {
    GNum <- base::nrow(data)
    subGN <- base::bquote(N[SNP]~"="~.(base::round(GNum, digits=3)))
  }
  else {subGN <- NULL}

  subT <- base::bquote(.(subGL)~.(subGN))

  # Processing step

  Pvals <- base::data.frame(Pobs = -base::log10(base::sort(data$Pval, decreasing=FALSE)))
  Pvals$Pexp <- -base::log10(stats::ppoints(base::length(data$Pval)))

  # Calculate upper limit to plot
  Pmin <- base::ceiling(base::max(-base::log10(data$Pval)))+0.1

  # Prune data to reduce the number of points to plot
  PvalsL <- Pvals %>% dplyr::filter(Pobs >= 3)
  Pvals3 <- Pvals %>% dplyr::filter(Pobs < 3 & Pobs >= 2)
  Pvals2 <- Pvals %>% dplyr::filter(Pobs < 2 & Pobs >= 1)
  Pvals1 <- Pvals %>% dplyr::filter(Pobs < 1)
  Pvals3 <- Pvals3[sample(1:nrow(Pvals3), nrow(PvalsL)),]
  Pvals2 <- Pvals2[sample(1:nrow(Pvals2), nrow(PvalsL)),]
  Pvals1 <- Pvals1[sample(1:nrow(Pvals1), nrow(PvalsL)),]
  Pvals <- rbind(PvalsL, Pvals3, Pvals2, Pvals1)

  # Plotting step
  plotdata <- ggplot2::ggplot(Pvals, ggplot2::aes(y=Pobs, x=Pexp)) +
    # Show all points
    ggplot2::geom_point(color=pcol, alpha=0.8) +
    # Show reference line
    ggplot2::geom_abline(slope=1, color=lcol) +
    # Adjust the axes breaks
    ggplot2::scale_x_continuous(breaks=base::seq(from=0, to=Pmin)) +
    ggplot2::scale_y_continuous(breaks=base::seq(from=0, to=Pmin)) +
    # Add labels
    ggplot2::labs(x=xlab, y=ylab, title=main, subtitle=subT) +
    # Apply bw theme
    ggplot2::theme_bw() +
    # Adjust the theme
    ggplot2::theme(
      # Remove legend
      legend.position="none",
      # Remove panel borders
      panel.border=ggplot2::element_blank(),
      # Remove and minor panel grids
      panel.grid.minor=ggplot2::element_blank(),
      # Adjust major panel grid y
      panel.grid.major=ggplot2::element_line(color="#212529", linetype="dotted"),
      # Set plot and panel backgrounds to white
      plot.background=ggplot2::element_rect(fill="white", color=NA),
      panel.background=ggplot2::element_rect(fill="white", color=NA),
      # Remove x and y ticks
      axis.ticks=ggplot2::element_blank(),
      # Adjust Axis tic labels
      axis.text.x=ggplot2::element_text(face="bold", color="#212529",
                                        margin=grid::unit(base::c(0,0,10,0), "pt")),
      axis.text.y=ggplot2::element_text(face="bold", color="#212529",
                                        margin=grid::unit(base::c(0,0,0,10), "pt")),
      # Adjust axis and graph label
      axis.title=ggplot2::element_text(face="bold", color="#212529"),
      plot.title=ggplot2::element_text(face="bold.italic", color="#212529")
    )

  # Return the plot
  base::print(plotdata)
}
