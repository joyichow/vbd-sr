library(metafor)
library(dplyr)
library(readxl)

transf_pe <- function(x) (1 - exp(-x)) * 100
pe_to_neg <- function(pe) -log(1 - pe/100)

pe_ticks <- c(-100, -50, 0, 25, 50, 75, 90, 99)
at_vals  <- pe_to_neg(pe_ticks)
at_vals  <- at_vals[is.finite(at_vals)]

# read + prepare data
dat <- read_excel("09042026_vbd_S1 Table.xlsx", sheet = 2) %>%
  rename(wolbachia_strain = wolbachia_type, metric = original_metric) %>%
  mutate(author = trimws(author))

# study IDs 
dat$study_id <- paste(dat$author, dat$year)
dat$study_id[dat$study_id == "dos Santos 2022"] <- paste0("dos Santos 2022 (", dat$virus[dat$study_id == "dos Santos 2022"], ")")
dat$study_id[dat$study_id == "Pinto 2021"] <- paste0("Pinto 2021 (", dat$virus[dat$study_id == "Pinto 2021"], ")")
toledo <- which(dat$author == "Toledo" & dat$year == 2017) # Toledo 2017: similar outcomes, distinguish by arm
dat$study_id[toledo] <- ifelse(grepl("IRS", dat$intervention_specific[toledo]),
                               "Toledo 2017 (IRS)", "Toledo 2017 (ITN)")

# reclassify intervention_class
dat <- dat %>% mutate(intervention_class = case_when(
  intervention_class == "Mosquito modification" & !is.na(wolbachia_strain) ~ "Wolbachia",
  intervention_class == "Mosquito modification" ~ "Sterile insect technique",
  grepl("IRS", intervention_specific) ~ "Indoor residual spraying",
  .default = intervention_class))

# effect sizes
dat$yi  <- log(dat$ratio)
dat$sei <- (log(dat$ratio_uci) - log(as.numeric(dat$ratio_lci))) / (2*1.96)

# Morrison reported one-sided CI
mi <- which(dat$study_id == "Morrison 2022")
dat$sei[mi] <- (log(dat$ratio_uci[mi]) - log(dat$ratio[mi])) / 1.645
dat$ratio_lci[mi] <- exp(dat$yi[mi] - 1.96*dat$sei[mi])

dat$yi_neg <- -dat$yi
dat$exclude_overlap[is.na(dat$exclude_overlap)] <- FALSE
dat$in_primary <- !dat$exclude_overlap
dat$drop_rob   <- dat$rob %in% c("Critical", "High")
dat$design_grp <- ifelse(dat$study_design %in% c("cRCT","RCT"), "cRCT", "NRS")

# main forest plot

make_forest <- function(plot_data, filename, main_title, sg_order,
                        page_width = 13, show_rob = TRUE, base_cex = 1.15) {
  
  # subgroup assignment
  plot_data$sg <- NA
  w <- plot_data$intervention_class == "Wolbachia" & !is.na(plot_data$wolbachia_strain)
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel" & plot_data$virus == "DENV"]  <- "Wolbachia wMel -- DENV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wAlbB"]  <- "Wolbachia wAlbB replacement -- DENV"
  plot_data$sg[w & plot_data$wolbachia_strain == "IIT-SIT"] <- "Wolbachia IIT-SIT -- DENV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel" & plot_data$virus == "CHIKV"] <- "Wolbachia wMel -- CHIKV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel" & plot_data$virus == "ZIKV"]  <- "Wolbachia wMel -- ZIKV"
  plot_data$sg[!w] <- plot_data$intervention_class[!w]
  
  dg <- ifelse(plot_data$study_design %in% c("cRCT","RCT"), "cRCT", "NRS")
  plot_data$subgroup <- paste0(plot_data$sg, " (", dg, ")")
  sg_order <- sg_order[sg_order %in% unique(plot_data$subgroup)]
  
  cur <- 1; ri <- list(); di <- list()
  v <- list(yi=c(), se=c(), lab=c(), row=c(), met=c(), rob=c(), sg=c())
  
  for (sg in rev(sg_order)) {
    d <- plot_data[plot_data$subgroup == sg, ]
    d <- d[order(-d$yi_neg), ]
    k <- nrow(d); if (k == 0) next
    
    drow <- NA
    if (k >= 2) { drow <- cur; cur <- cur+2 }
    srows <- cur:(cur+k-1); cur <- cur+k
    hrow <- cur + 0.8; cur <- cur + 2.5
    ri[[sg]] <- list(h=hrow, s=srows, d=drow)
    
    labs <- d$study_id
    for (i in seq_along(labs))
      if (d$virus[i] != "DENV" && !grepl(d$virus[i], labs[i]))
        labs[i] <- paste0(labs[i], " [", d$virus[i], "]")
    
    v$yi  <- c(v$yi, d$yi_neg);  v$se <- c(v$se, d$sei)
    v$lab <- c(v$lab, labs);     v$row <- c(v$row, srows)
    v$met <- c(v$met, d$metric); v$rob <- c(v$rob, d$rob)
    v$sg  <- c(v$sg, rep(sg, k))
    
    if (k >= 2) {
      fit <- tryCatch(rma(yi=d$yi, sei=d$sei, method="REML", test="knha"),
                      error = function(e) tryCatch(
                        rma(yi=d$yi, sei=d$sei, method="DL", test="knha"),
                        error = function(e2) NULL))
      di[[sg]] <- list(res=fit, row=drow, k=k)
    }
  }
  
  # within-subgroup weights
  wt <- rep(NA_real_, length(v$yi))
  for (sg in names(ri)) {
    idx <- which(v$sg == sg)
    if (sg %in% names(di) && !is.null(di[[sg]]$res))
      wt[idx] <- round(weights(di[[sg]]$res), 1)
    else wt[idx] <- 100
  }
  
  if (show_rob) {
    ilab <- cbind(v$met, v$rob, sprintf("%.1f", wt))
    ix <- c(-7.5, -5.5, -3.5); hd <- c("Metric","RoB","Weight (%)")
    xl <- -13.2
  } else {
    ilab <- cbind(v$met, sprintf("%.1f", wt))
    ix <- c(-9.0, -5.5); hd <- c("Metric","Weight (%)")
    xl <- -13.5
  }
  maxr <- cur
  ph   <- max(7, length(v$yi)*0.65 + 5)
  
  pdf(filename, width=page_width, height=ph, pointsize=14)
  par(mar=c(5, 0.5, 4, 2))
  forest(x=v$yi, sei=v$se, slab=v$lab,
         atransf=transf_pe, at=at_vals,
         alim=c(pe_to_neg(-100), pe_to_neg(99)),
         xlim=c(xl, 7.5), ylim=c(-1, maxr+4),
         rows=v$row, ilab=ilab, ilab.xpos=ix,
         cex=base_cex, header=TRUE, main=main_title,
         xlab="Protective efficacy, PE (%)",
         refline=0, psize=1.2, efac=0.7, digits=0L)
  text(ix, maxr+2.8, hd, cex=base_cex*0.85, font=2)
  
  for (sg in names(ri)) {
    r <- ri[[sg]]
    text(xl, r$h, sg, pos=4, cex=base_cex*0.92, font=2, col="grey30")
    if (!(sg %in% names(di)) || is.null(di[[sg]]$res)) next
    
    res <- di[[sg]]$res
    dcol <- ifelse(grepl("Wolbachia|Sterile", sg), "darkblue", "darkred")
    addpoly(x=-res$beta[1], sei=res$se,
            ci.lb=-res$ci.ub, ci.ub=-res$ci.lb,
            atransf=transf_pe, row=r$d,
            mlab=sprintf("Heterogeneity: I² = %.1f%%, Tau² = %.3f, p %s",
                         res$I2, res$tau2,
                         ifelse(res$QEp < 0.001, "< 0.001", sprintf("= %.3f", res$QEp))),
            cex=base_cex*0.95, col=dcol, border=dcol, digits=0L)
  }
  
  mtext("Random-effects model (REML + HKSJ). Weights are within-subgroup RE weights. PE = (1 - ratio) x 100%.",
        side=1, line=3.8, cex=base_cex*0.8, adj=0)
  dev.off()
  cat("Saved:", filename, "\n")
  invisible(di)
}

# subgroup orderings 
mosq_base <- c("Wolbachia wMel -- DENV", "Wolbachia wAlbB replacement -- DENV",
               "Wolbachia IIT-SIT -- DENV", "Wolbachia wMel -- CHIKV",
               "Wolbachia wMel -- ZIKV", "Sterile insect technique")
other_base <- c("Insecticide-treated materials", "Indoor residual spraying",
                "Autodissemination", "Environmental management",
                "Spatial repellents", "Larvivorous predators",
                "Mosquito repellents", "Integrated programs")
sg_mosq  <- c(paste0(mosq_base, " (cRCT)"),  paste0(mosq_base, " (NRS)"))
sg_other <- c(paste0(other_base, " (cRCT)"), paste0(other_base, " (NRS)"))
is_mosq  <- dat$intervention_class %in% c("Wolbachia", "Sterile insect technique")

# plots A + B 
di_A1 <- make_forest(dat[dat$in_primary & !dat$drop_rob & is_mosq, ],
                     "forest_A_mosqmod_rob.pdf", "Mosquito modification interventions", sg_mosq)

di_B1 <- make_forest(dat[dat$in_primary & !dat$drop_rob & !is_mosq, ],
                     "forest_B_other_rob.pdf", "Other Aedes-borne disease control interventions", sg_other)

di_A2 <- make_forest(dat[dat$in_primary & is_mosq, ],
                     "forest_A_mosqmod_design.pdf", "Mosquito modification interventions", sg_mosq)

di_B2 <- make_forest(dat[dat$in_primary & !is_mosq, ],
                     "forest_B_other.pdf", "Other Aedes-borne disease interventions", sg_other)

# funnel and Egger's 
wall <- dat[dat$in_primary & dat$intervention_class == "Wolbachia" & dat$virus == "DENV", ]
res_wall <- rma(yi = yi, sei = sei, data = wall, method = "REML", test = "knha")
egt <- regtest(res_wall, model = "lm")

pdf("funnel_final.pdf", width = 8, height = 6)
funnel(res_wall, atransf = exp, main = "Funnel Plot: All Wolbachia (DENV outcomes)",
       xlab = "Ratio measure", back = "white", shade = "gray90")
text(res_wall$beta[1], 0, sprintf("Pooled = %.2f", exp(res_wall$beta[1])),
     pos = 4, cex = 0.8, offset = 0.5)
mtext(sprintf("Egger's test: z = %.2f, p = %s",
              egt$zval, ifelse(egt$pval < 0.001, "< 0.001", sprintf("%.3f", egt$pval))),
      side = 1, line = 3.5, cex = 0.75, adj = 0)
dev.off()
