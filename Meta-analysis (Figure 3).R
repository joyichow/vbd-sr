# Nov 01, 2026
library(metafor)
library(dplyr)
library(readxl)

transf_pe <- function(x) (1 - exp(-x)) * 100
pe_to_neg <- function(pe) -log(1 - pe/100)

pe_ticks <- c(-100, -50, 0, 25, 50, 75, 90, 100)
at_vals  <- pe_to_neg(pe_ticks)
at_vals  <- at_vals[is.finite(at_vals)]

dat <- read_excel("28092026_vbd_S1 Table.xlsx", sheet = 2, na = c("", "NA")) %>%
  rename(wolbachia_strain = wolbachia_type, metric = original_metric) %>%
  mutate(author = trimws(author))

dat$study_id <- paste(dat$author, dat$year)
toledo <- which(dat$author == "Toledo" & dat$year == 2017)
dat$study_id[toledo] <- ifelse(grepl("IRS", dat$intervention_specific[toledo]),
                               "Toledo 2017 (IRS)", "Toledo 2017 (ITN)")

dat <- dat %>% mutate(intervention_class = case_when(
  intervention_class == "Mosquito modification" & !is.na(wolbachia_strain) ~ "Wolbachia",
  intervention_class == "Mosquito modification" ~ "Sterile insect technique",
  grepl("IRS", intervention_specific) ~ "Indoor residual spraying",
  .default = intervention_class))

dat$yi  <- log(dat$ratio)
dat$sei <- (log(dat$ratio_uci) - log(as.numeric(dat$ratio_lci))) / (2*1.96)

# Morrison reported one-sided CI
mi <- which(dat$study_id == "Morrison 2022")
dat$virus[mi] <- "ABV" # outcome was any Aedes-borne virus, not dengue alone
dat$sei[mi] <- (log(dat$ratio_uci[mi]) - log(dat$ratio[mi])) / 1.645
dat$ratio_lci[mi] <- exp(dat$yi[mi] - 1.96*dat$sei[mi])

dat$yi_neg <- -dat$yi

# reported CI limits, for plotting only (pooling uses sei)
dat$ci_lo <- -log(dat$ratio_uci)
dat$ci_hi <- -log(as.numeric(dat$ratio_lci))
dat$ci_lo[mi] <- dat$yi_neg[mi] - 1.96*dat$sei[mi] # Morrison: keep the two-sided approximation
dat$exclude_overlap[is.na(dat$exclude_overlap)] <- FALSE
dat$in_primary <- !dat$exclude_overlap
dat$drop_rob   <- dat$rob %in% c("Critical", "High")

# one section = one coloured banner, with its subgroups split into cRCT / NRS
section <- function(label, col, base, display = NULL) {
  groups <- c()
  for (i in seq_along(base)) {
    for (d in c("cRCT", "NRS")) {
      if (is.null(display)) {
        lab <- ifelse(d == "cRCT", "Cluster-randomised trials", "Non-randomised studies")
      } else {
        lab <- paste0(display[i], " (", d, ")")
      }
      groups[paste0(base[i], " (", d, ")")] <- lab
    }
  }
  list(label = label, col = col, groups = groups)
}

make_forest <- function(plot_data, filename, main_title, sections,
                        page_width = 13, show_rob = TRUE, pool = TRUE, base_cex = 1.15,
                        panel = NULL, measure = FALSE) {
  
  plot_data$sg <- NA
  w <- plot_data$intervention_class == "Wolbachia" & !is.na(plot_data$wolbachia_strain)
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel"  & plot_data$virus == "DENV"]  <- "Wolbachia wMel replacement – DENV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wAlbB"]                               <- "Wolbachia wAlbB replacement – DENV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel"  & plot_data$virus == "CHIKV"] <- "Wolbachia wMel replacement – CHIKV"
  plot_data$sg[w & plot_data$wolbachia_strain == "wMel"  & plot_data$virus == "ZIKV"]  <- "Wolbachia wMel replacement – ZIKV"
  plot_data$sg[w & plot_data$wolbachia_strain == "IIT-SIT"]                             <- "Wolbachia wAlbB IIT-SIT suppression – DENV"
  
  plot_data$sg[!w] <- ifelse(
    plot_data$intervention_class[!w] == "Sterile insect technique",
    "Sterile insect technique suppression",
    paste0(plot_data$intervention_class[!w], " -- ", plot_data$virus[!w])
  )
  
  dg <- ifelse(plot_data$study_design %in% c("cRCT","RCT"), "cRCT", "NRS")
  plot_data$subgroup <- paste0(plot_data$sg, " (", dg, ")")
  
  all_keys <- unlist(lapply(sections, function(s) names(s$groups)))
  missed <- setdiff(unique(plot_data$subgroup), all_keys)
  if (length(missed) > 0)
    warning(filename, ": subgroups not plotted: ", paste(missed, collapse = "; "), call. = FALSE)
  
  # rows are built from the bottom of the page upwards
  cur <- 1
  ri <- list()
  banners <- list()
  v <- list(yi=c(), se=c(), lo=c(), hi=c(), lab=c(), row=c(), met=c(), rob=c(), wt=c())
  
  for (s in rev(sections)) {
    keys <- names(s$groups)
    keys <- keys[keys %in% plot_data$subgroup]
    if (length(keys) == 0) next
    
    for (sg in rev(keys)) {
      d <- plot_data[plot_data$subgroup == sg, ]
      d <- d[order(-d$yi_neg), ]
      k <- nrow(d)
      
      fit <- NULL
      if (pool & k >= 2) {
        fit <- try(rma(yi=d$yi, sei=d$sei, method="REML", test="knha"), silent = TRUE)
        if (inherits(fit, "try-error"))
          fit <- try(rma(yi=d$yi, sei=d$sei, method="DL", test="knha"), silent = TRUE)
        if (inherits(fit, "try-error")) fit <- NULL
      }
      
      hetrow <- NA
      drow <- NA
      if (!is.null(fit)) {
        hetrow <- cur
        drow <- cur + 1.1
        cur <- cur + 2.4
      }
      srows <- cur + 0:(k-1)
      cur <- cur + k
      ri[[sg]] <- list(h = cur + 0.2, d = drow, het = hetrow, k = k, fit = fit,
                       lab = s$groups[[sg]], col = s$col)
      cur <- cur + 1.9
      
      if (is.null(fit)) wts <- rep(NA, k) else wts <- weights(fit)
      
      v$yi  <- c(v$yi, d$yi_neg)
      v$se  <- c(v$se, d$sei)
      v$lo  <- c(v$lo, d$ci_lo)
      v$hi  <- c(v$hi, d$ci_hi)
      v$lab <- c(v$lab, d$study_id)
      v$row <- c(v$row, srows)
      v$met <- c(v$met, d$metric)
      v$rob <- c(v$rob, d$rob)
      v$wt  <- c(v$wt, wts)
    }
    
    banners[[length(banners) + 1]] <- list(row = cur - 0.3, label = s$label, col = s$col)
    cur <- cur + 1.7
  }
  
  ilab <- cbind(v$met)
  hd <- "Metric"
  if (show_rob) {
    ilab <- cbind(ilab, v$rob)
    hd <- c(hd, "RoB")
  }
  if (pool) {
    ilab <- cbind(ilab, ifelse(is.na(v$wt), "", sprintf("%.1f", v$wt)))
    hd <- c(hd, "Weight (%)")
  }
  if (length(hd) == 3) ix <- c(-4.7, -3.3, -2.0)
  if (length(hd) == 2) ix <- c(-4.7, -2.8)
  if (length(hd) == 1) ix <- -3.5
  xl <- -8
  
  maxr <- banners[[length(banners)]]$row
  ph   <- max(7, (maxr + 5) * 0.31 + 2)
  
  if (measure) return(ph)
  
  # filename = NULL draws on whatever device is already open (used for the combined figure)
  if (!is.null(filename)) pdf(filename, width=page_width, height=ph, pointsize=14)
  par(mar=c(5, 0.5, 4, 2))
  forest(x=v$yi, ci.lb=v$lo, ci.ub=v$hi, slab=NA,
         atransf=transf_pe, at=at_vals,
         alim=c(pe_to_neg(-100), pe_to_neg(99)),
         xlim=c(xl, 6.6), ylim=c(-1, maxr+3.2),
         rows=v$row, ilab=ilab, ilab.xpos=ix,
         cex=base_cex, header=FALSE, main=main_title,
         xlab="Protective effectiveness/efficacy, PE (%)",
         refline=0, psize=1.2, efac=0.7, digits=0L)
  text(xl, maxr+2.2, "Study", pos=4, cex=base_cex, font=2)
  text(ix, maxr+2.2, hd, cex=base_cex, font=2)
  text(6.6, maxr+2.2, "Estimate [95% CI]", pos=2, cex=base_cex, font=2)
  text(xl, v$row, v$lab, pos=4, cex=base_cex, col="grey40")
  
  usr <- par("usr")
  
  # widths of the study annotations, so the pooled numbers line up with them
  w_est <- max(nchar(round(transf_pe(v$yi))))
  w_lb  <- max(nchar(round(transf_pe(v$lo))))
  w_ub  <- max(nchar(round(transf_pe(v$hi))))
  for (b in banners) {
    segments(usr[1], b$row - 0.65, usr[2], b$row - 0.65, col = b$col)
    text(xl, b$row, b$label, pos=4, cex=base_cex*1.18, font=2, col=b$col)
  }
  
  for (sg in names(ri)) {
    r <- ri[[sg]]
    text(xl, r$h, r$lab, pos=4, cex=base_cex*1.08, font=2)
    if (is.null(r$fit)) next
    res <- r$fit
    
    addpoly(x=-res$beta[1], sei=res$se,
            ci.lb=-res$ci.ub, ci.ub=-res$ci.lb,
            atransf=transf_pe, rows=r$d, mlab="", annotate=FALSE,
            col=r$col, border=r$col)
    pooled <- transf_pe(-c(res$beta[1], res$ci.ub, res$ci.lb))
    pooled <- paste0(formatC(pooled[1], format="f", digits=0, width=w_est), " [",
                     formatC(pooled[2], format="f", digits=0, width=w_lb), ", ",
                     formatC(pooled[3], format="f", digits=0, width=w_ub), "]")
    text(6.6, r$d, pooled, pos=2, cex=base_cex, font=2, col=r$col)
    des <- ifelse(grepl("(cRCT)", sg, fixed = TRUE), "cRCT", "NRS")
    text(xl, r$d, paste0("Pooled estimate (", des, ", k = ", r$k, ")"),
         pos=4, cex=base_cex, font=2, col=r$col)
    
    i2 <- sprintf("%.1f%%", res$I2)
    t2 <- sprintf("%.3f", res$tau2)
    pq <- ifelse(res$QEp < 0.001, "< 0.001", sprintf("= %.3f", res$QEp))
    het <- bquote(paste("Heterogeneity: ", I^2, " = ", .(i2), ", ",
                        tau^2, " = ", .(t2), ", p ", .(pq)))
    text(xl, r$het, het, pos=4, cex=base_cex*0.85)
  }
  
  note <- "PE = (1 - ratio) x 100%. Study CIs as reported."
  if (pool) note <- paste("Random-effects model (REML + HKSJ). cRCTs and NRS are pooled separately; single studies are not pooled. Weights are within-subgroup RE weights.")
  mtext(note, side=1, line=3.8, cex=base_cex*0.8, adj=0)
  if (!is.null(panel)) mtext(panel, side=3, line=1.5, adj=0, font=2, cex=base_cex*1.5)
  if (!is.null(filename)) {
    dev.off()
    cat("Saved:", filename, "\n")
  }
  invisible(ri)
}

################################################################################
################################## Figure 3A ###################################
################################################################################
sec_mosq <- list(
  section("Population replacement: wMel Wolbachia", "darkblue",
          c("Wolbachia wMel replacement – DENV", "Wolbachia wMel replacement – CHIKV",
            "Wolbachia wMel replacement – ZIKV"),
          c("DENV", "CHIKV", "ZIKV")),
  section("Population replacement: wAlbB Wolbachia", "skyblue4",
          "Wolbachia wAlbB replacement – DENV", "DENV"),
  section("Population suppression (IIT-SIT): wAlbB Wolbachia", "cyan4",
          "Wolbachia wAlbB IIT-SIT suppression – DENV", "DENV"),
  section("Population suppression: sterile insect technique", "cyan4",
          "Sterile insect technique suppression", "DENV")
)

################################################################################
################################## Figure 3B ###################################
################################################################################
viruses <- c("DENV", "CHIKV", "ZIKV", "ABV")
sec_other <- lapply(
  c("Autodissemination", "Environmental management",
    "Insecticide-treated materials", "Indoor residual spraying", "Integrated programs",
    "Larvivorous predators", "Mass trapping",
    "Spatial repellents"),
  function(x) section(x, "firebrick3", paste0(x, " -- ", viruses), viruses))

is_mosq <- dat$intervention_class %in% c("Wolbachia", "Sterile insect technique")

################################################################################
############################### Combine figures ################################
################################################################################
di_A2 <- make_forest(dat[dat$in_primary & is_mosq, ],
                     "forest_A_mosqmod_design.pdf", "Mosquito modification interventions", sec_mosq)

di_B2 <- make_forest(dat[dat$in_primary & !is_mosq, ],
                     "forest_B_other.pdf", "Other Aedes-borne disease interventions", sec_other)

# main figure: A and B on one page
dA <- dat[dat$in_primary & !dat$drop_rob & is_mosq, ]
dB <- dat[dat$in_primary & !dat$drop_rob & !is_mosq, ]
titleA <- "Mosquito modification interventions"
titleB <- "Other Aedes-borne disease control interventions"

hA <- suppressWarnings(make_forest(dA, NULL, titleA, sec_mosq, measure = TRUE))
hB <- suppressWarnings(make_forest(dB, NULL, titleB, sec_other, measure = TRUE))

pdf("fig3_combined.pdf", width = 13, height = hA + hB, pointsize = 14)
layout(matrix(1:2), heights = c(hA, hB))
make_forest(dA, NULL, titleA, sec_mosq, panel = "A")
make_forest(dB, NULL, titleB, sec_other, panel = "B")
dev.off()

################################################################################
############################### Funnel & Egger's ###############################
################################################################################
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
