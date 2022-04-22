library(future)
# check the current active plan

# change the current plan to access parallelization
plan("multiprocess", workers = 4)
plan()


library(Seurat)
pbmc <- readRDS("./pbmc3k_final.rds")

# Enable parallelization
plan("multiprocess", workers = 4)
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)


##################################################################
#Benchmarks
timing.comparisons <-
  data.frame(fxn = character(),
             time = numeric(),
             strategy = character())

##################################################################

plan("sequential")
start <- Sys.time()
pbmc <-
  ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "ScaleData",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "sequential"
    )
  )

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "FindMarkers",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "sequential"
    )
  )

##################################################################
plan("multiprocess", workers = 2)
start <- Sys.time()
pbmc <-
  ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "ScaleData-2cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-2cpu"
    )
  )

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "FindMarkers-2cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-2cpu"
    )
  )
##################################################################
plan("multiprocess", workers = 3)
start <- Sys.time()
pbmc <-
  ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "ScaleData-3cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-3cpu"
    )
  )

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "FindMarkers-3cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-3cpu"
    )
  )

##################################################################
plan("multiprocess", workers = 4)
start <- Sys.time()
pbmc <-
  ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "ScaleData-4cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-4cpu"
    )
  )

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <-
  rbind(
    timing.comparisons,
    data.frame(
      fxn = "FindMarkers-4cpu",
      time = as.numeric(end -
                          start, units = "secs"),
      strategy = "multiprocess-4cpu"
    )
  )
##################################################################
timing.comparisons$fxn <- factor(timing.comparisons$fxn, levels = timing.comparisons$fxn)

library(ggplot2)
library(cowplot)
ggplot(timing.comparisons, aes(fxn, time)) + 
  geom_bar(aes(fill = strategy), 
           stat = "identity", 
           position = "dodge") +
  ylab("Time(s)") + 
  xlab("Function") + 
  theme_cowplot()+ 
  theme(axis.text.x=element_text(size=10))





