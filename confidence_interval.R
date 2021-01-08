# Confidence Interval script to process nucleotide diversity and tajima's D values per bin
packages = c("ggplot2", "ggrepel", "vcfR")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

### plotting confidence intervals in pi_tajima.pdf file.

bins<-read.table("input.pi.D", header=F)
names(bins)[names(bins) == "V2"] <- "nucleotide_diversity"
names(bins)[names(bins) == "V3"] <- "D"
upper_limit=quantile(bins$D, .025)
upper_limit
lower_limit=quantile(bins$D, .975)
lower_limit
pdf("pi_tajima.pdf")
gg <- ggplot(bins, aes(x=nucleotide_diversity, y=D))+
geom_point()+
  geom_hline(yintercept=lower_limit, linetype="dashed", color = "red")+  # 97.5% interval of D column
  geom_hline(yintercept=upper_limit, linetype="dashed", color = "red")  # 2.5% interval of D column
### geom_label_repel
gg + 
  geom_label_repel(aes(label = V1),
                  data          = subset(bins, D < upper_limit),
                  size          = 2.5,
                  box.padding   = 0.4,
                  point.padding = 0.4,
                  segment.size  = 0.4,
                  segment.color = "grey50",
                  direction     = "x",
                  max.overlaps  = Inf) +
  geom_label_repel(aes(label = V1),
                  data         = subset(bins, D > lower_limit),
                  size          = 2.5,
                  box.padding   = 0.4,
                  point.padding = 0.4,
                  segment.size  = 0.4,
                  segment.color = "grey50",
                  direction     = "x",
                  max.overlaps  = Inf) +
  theme_classic(base_size = 16)
plot(gg)
dev.off()

### getting 2.5 and 97.5% intervals from data using subset function
myvars <- c("V1", "nucleotide_diversity", "D")
upper_bins_subset <- subset(bins, bins$D < upper_limit)
newdata_upper <- upper_bins_subset[myvars]
lower_bins_subset <- subset(bins, bins$D > lower_limit)
newdata_lower <- lower_bins_subset[myvars]
write.table(newdata_upper, file="bins_2.5%_confidence.tab", sep="\t", quote = FALSE, row.names = FALSE)
write.table(newdata_lower, file="bins_97.5%_confidence.tab", sep="\t", quote = FALSE, row.names = FALSE)
proc.time()
sessionInfo()
quit("no")
