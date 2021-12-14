##' \emph{plot_manhattan} takes chromosome, base pair, and p-value data and creates a ggplot based Manhattan plot.
##'
##' Basic functionality: \emph{plot_manhattan} takes chromosome and basepair address information, as well as p values from the GWAS summary statistics. Subsequently, data are additionally processed to generate Manhattan plot where chromosome and basepair information are used to generate x-axis, and the \eqn{-log_10 (P)} values are plotted on the y axis.
##'
##' @title Manhattan Plot
##' @param data Data frame which contains GWAS summary statistics (REQUIRED)
##' @param P A column name which contains P values, default is "P", column should be numeric.
##' @param CHR A column name which contains chromosome ID, default is "CHR", column should be numeric or character.
##' @param BP A column name which contains base pair location, default is "BP", column should be numeric.
##' @param maxP A maximum P-value cutoff (sometimes, it's not a good use of memory to plot millions of p values above 0.1), default is 1.
##' @param xlab Optional x-axis label, default is NULL.
##' @param ylab Optional y-axis label, default is "-log10(P)"
##' @param main Optional plot title, default is "Manhattan plot"
##' @param colrs Optional character vector of colors, default is set by package.
##' @param sigA Optional adjustment to stringent p-value reference, default is 3e-8
##' @param sigAcol Option string to adjust the color of sigA line
##' @param sigB Optional adjustment to stringent p-value reference, default is 1e-5
##' @param sigBcol Option string to adjust the color of sigB line
##' @param ylims Optional adjustments to y limits, default adjusts it to data
##' @param ybreaks Optional adjsuutemts to y breaks, default adds whole numbers in increments of 1
##' @return ggplot object
##' @author Franjo, Mingjing, and Xiaoxiao
##' @export
##' @examples
##' plot_manhattan(ADHDMeta)
##' @importFrom dplyr %>%

plot_manhattan <- function(data=data, P="P", CHR="CHR", BP="BP",
                           maxP=NULL, xlab=NULL, ylab=NULL, main=NULL,
                           colrs=NULL, sigA=NULL, sigAcol=NULL, sigB=NULL,
                           sigBcol=NULL, ylims=NULL, ybreaks=NULL){

  # If dataframe is not provided stop function and show  error
  if(base::is.null(data)){
    base::stop("Data frame is not specified.")
  }

  # Preprepreprocessing
  dondata <- data %>%
    # Rename variables
    dplyr::rename(CHR={{CHR}}, P={{P}}, BP={{BP}})

  # If p values aren't numeric
  if(base::class(dondata$P) != "numeric"){
    base::stop("P values are not class numeric.")
  }

  # If bp values aren't numeric
  if(base::class(dondata$BP) != "numeric"){
    base::stop("BP values are not class numeric.")
  }

  # If bp values aren't numeric
  if(!base::class(dondata$CHR) %in% base::c("numeric", "character")){
    base::warning("CHR values are not class numeric or character. This might
                  interfere with data processing")
  }

  # If not provided, define y-limits
  if(base::is.null(ylims)){
    ylims=base::c(base::floor(-base::log10(base::max(base::unlist(data[P])))),
                  -base::log10(base::min(base::unlist(data[P]))))+1}

  # If not provided, define y-breaks
  if(base::is.null(ybreaks)){
    ybreaks=base::seq(from=base::ceiling(-base::log10(base::max(base::unlist(data[P])))),
                      to=base::ceiling(-base::log10(base::min(base::unlist(data[P])))),
                      by=1)}

  # If not provided, define maxP value
  if(base::is.null(maxP)){
    maxP=1}

  # If not provided, define y label
  if(base::is.null(ylab)){
    ylab=base::bquote("-log"[10]~"(P)")}

  # If not provided, define graph title
  if(base::is.null(main)){
    main="Manhattan Plot"}

  # If not provided, define point colors
  if(base::is.null(colrs)){
    colrs=base::rep(base::c("#646E78", "#1E1014"),
                    base::ceiling(base::length(base::unlist(base::unique(data[CHR])))))}
  else{
    colrs=base::rep(colrs,
                    base::ceiling(base::length(base::unlist(base::unique(data[CHR])))))}

  # If not provided, define stringent significance threshold
  if(base::is.null(sigA)){
    sigA=5e-8}

  # If not provided, define lax significance threshold
  if(base::is.null(sigB)){
    sigB=1e-5}

  # If not provided, define stringent significance threshold color
  if(base::is.null(sigAcol)){
    sigAcol="#A61C3C"}

  # If not provided, define lax significance threshold color
  if(base::is.null(sigBcol)){
    sigBcol="#0B0033"}

  # Prepreprocessing step
  dondata <- dondata %>%
    # Facorize and level chromosomes
    dplyr::mutate(CHR=base::factor(CHR,levels=base::sort(base::as.numeric(
      base::unlist(base::unique(CHR))))))

  # Preprocessing step
  dondata <- dondata %>%
    # Compute chromosome size
    dplyr::group_by(CHR) %>%
    dplyr::summarise(chr_len=max(BP)) %>%
    # Calculate cumulative position of each chromosome
    dplyr::mutate(tot=base::cumsum(base::as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    dplyr::left_join(dondata, ., by=base::c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    dplyr::arrange(CHR, BP) %>%
    dplyr::mutate(BPcum=BP+tot) %>%
    # Filter in only P values less than specified
    dplyr::filter(P < maxP)

  # Axis adjsutment
  axisdata <- dondata %>%
    # Find the middle positions for each chromosome
    dplyr::group_by(CHR) %>%
    dplyr::summarise(center=(base::max(BPcum)+base::min(BPcum))/2)

  # Plotting step
  plotdata <- ggplot2::ggplot(dondata, ggplot2::aes(x=BPcum, y=-base::log10(P))) +
    # Show all GWAS p-values
    ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)),alpha=0.8, size=1.3) +
    # Show significane lines
    ggplot2::geom_hline(yintercept=-base::log10(sigA), color=sigAcol) +
    ggplot2::geom_hline(yintercept=-base::log10(sigB), color=sigBcol) +
    # Change colors of points
    ggplot2::scale_color_manual(values=colrs) +
    # Adjust the x-axis
    ggplot2::scale_x_continuous(label=axisdata$CHR,
                                breaks=axisdata$center,
                                expand=base::c(0.03, 0.03)) +
    # Limit the y axis
    ggplot2::coord_cartesian(ylim=ylims) +
    # Adjust the y-axis
    ggplot2::scale_y_continuous(breaks=ybreaks) +
    # Add labels
    ggplot2::labs(x=xlab, y=ylab, title=main) +
    # Apply bw theme
    ggplot2::theme_bw() +
    # Adjust the theme
    ggplot2::theme(
      legend.position="none",
      panel.border=ggplot2::element_blank(),
      panel.grid.major.x=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major.y=ggplot2::element_line(color="#212529", linetype="dotted"),
      plot.background=ggplot2::element_rect(fill="white", color=NA),
      panel.background=ggplot2::element_rect(fill="white", color=NA),
      axis.ticks=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(face="bold", color="#212529"),
      axis.text.y=ggplot2::element_text(face="bold", color="#212529",
                                        margin=grid::unit(base::c(0,0,0,10), "pt")),
      axis.title=ggplot2::element_text(face="bold", color="#212529"),
      plot.title=ggplot2::element_text(face="bold.italic", color="#212529")
    )

  # Return the plot
  base::print(plotdata)
}
