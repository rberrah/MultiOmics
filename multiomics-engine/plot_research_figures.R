#!/usr/bin/env Rscript
# MultiOmics display-equivalent publication figures (base R + jsonlite).
# Plot the SAME display-ready values and sample order that were shown in-browser.
# This reproduces plotted data, not the Svelte vector geometry pixel-for-pixel.
# Run from an extracted research ZIP:
# Rscript figures_r/generate_figures.R provenance/analysis_complete.json figures_r/output
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript generate_figures.R ANALYSIS_JSON OUTPUT_DIRECTORY")
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("R package jsonlite is required")
source_path <- normalizePath(args[[1]],mustWork=TRUE)
output <- args[[2]]
dir.create(output,showWarnings=FALSE,recursive=TRUE)
record <- jsonlite::fromJSON(source_path,simplifyVector=FALSE)
visuals <- record$visualizations
if(is.null(visuals)) stop("This analysis contains no archived display-ready visualization data")
palette <- grDevices::colorRampPalette(c("#3567B2","#F7F7F7","#CC3E3E"))(201L)
color_for <- function(value, cap=3) {
  if(is.null(value) || length(value)!=1L || !is.finite(suppressWarnings(as.numeric(value))))
    return("#E0E0E0")
  index <- as.integer(round(100+100*max(-1,min(1,as.numeric(value)/max(cap,1e-9)))))+1L
  palette[max(1L,min(201L,index))]
}
num <- function(x) {
  if(is.null(x)||length(x)!=1L)return(NA_real_)
  z <- suppressWarnings(as.numeric(x))
  if(length(z)!=1L || !is.finite(z)) NA_real_ else z
}
chr <- function(x,otherwise="") if(is.null(x)||!length(x))otherwise else as.character(x[[1L]])
# Visual consistency contract, independent of statistical inference.
legend_values <- c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5)
utils::write.csv(data.frame(z_score=legend_values,
  fill_hex=vapply(legend_values,color_for,character(1),cap=2.5)),
  file.path(output,"heatmap_color_legend.csv"),row.names=FALSE)
outputs <- list()
with_devices <- function(name,width,height,draw,kind,details) {
  for(extension in c("pdf","svg")) {
    path <- file.path(output,paste0(name,".",extension))
    if(extension=="svg"&&!capabilities("cairo")) {
      cat("SKIP SVG (R was built without Cairo): ",name,"\n")
      next
    }
    if(extension=="pdf")grDevices::pdf(path,width=width,height=height,useDingbats=FALSE)
    else grDevices::svg(path,width=width,height=height,onefile=TRUE)
    error <- tryCatch({draw();NULL},error=function(e)e)
    grDevices::dev.off()
    if(!is.null(error))stop("Figure "+name+" failed: "+conditionMessage(error))
    outputs[[length(outputs)+1L]] <<- list(filename=basename(path),kind=kind,
      generated_by="R_base_graphics",comparison_to_browser="scientifically_similar_not_pixel_identical",
      plotted_data_source=details,
      bytes=as.integer(file.info(path)$size))
  }
}
heatmaps <- visuals$heatmaps
if(!is.null(heatmaps))for(layer in names(heatmaps)) {
  h <- heatmaps[[layer]]
  if(is.null(h$rows)||!length(h$rows)||is.null(h$samples)||!length(h$samples))next
  features <- vapply(h$rows,function(r)chr(r$feature),character(1))
  samples <- vapply(h$samples,function(r)chr(r$sampleId),character(1))
  if(any(!nzchar(features))||any(!nzchar(samples)))stop("Missing heatmap identifiers")
  values <- do.call(rbind,lapply(h$rows,function(r) {
    if(length(r$z)!=length(samples))stop("Archived heatmap row/sample mismatch")
    vapply(r$z,num,numeric(1))
  }))
  if(!identical(dim(values),c(length(features),length(samples))))stop("Heatmap dimension mismatch")
  # Exact plotted cells: users can verify every color against this numeric CSV.
  plot_values <- data.frame(feature_id=features,values,check.names=FALSE)
  names(plot_values) <- c("feature_id",samples)
  utils::write.csv(plot_values,file.path(output,paste0("heatmap_",layer,"_plot_values.csv")),
                   row.names=FALSE,na="NA")
  local({
    values0<-values; feat0<-features; sam0<-samples; layer0<-layer
    with_devices(paste0("heatmap_",layer0),
       max(8,4+0.19*length(sam0)),max(6,2+0.22*length(feat0)),function(){
        op<-par(mar=c(6,8,3,2));on.exit(par(op))
        nr<-nrow(values0);nc<-ncol(values0)
        plot(c(0.5,nc+0.5),c(0.5,nr+0.5),type="n",axes=FALSE,
             xlab="",ylab="",main=paste("MultiOmics",layer0,"— row z-score"))
        for(i in seq_len(nr))for(j in seq_len(nc))
          rect(j-.5,nr-i+.5,j+.5,nr-i+1.5,
               col=color_for(values0[i,j],2.5),border="#FFFFFF")
        axis(1,at=seq_len(nc),labels=sam0,las=2,cex.axis=.7,tick=FALSE)
        axis(2,at=seq_len(nr),labels=rev(feat0),las=2,cex.axis=.7,tick=FALSE)
        box(col="#BBBBBB")
        mtext("Grey = missing value; blue = lower, red = higher (within feature)",
              side=1,line=4,cex=.7)
       },"heatmap",paste0("visualizations.heatmaps.",layer0,".rows[].z"))
  })
}
# Circular effects: same ranked feature effects and q-values, distinct layer scales.
meta <- visuals$metabologram
if(!is.null(meta)) {
  transcript<-meta$transcriptomics
  metabolites<-meta$metabolomics
  if(length(transcript)||length(metabolites)){
    plot_effects <- do.call(rbind,lapply(c("transcriptomics","metabolomics"),function(layer){
      l <- meta[[layer]]
      if(!length(l))return(NULL)
      data.frame(omic=rep(layer,length(l)),
        feature=vapply(l,function(e)chr(e$feature),character(1)),
        effect=vapply(l,function(e)num(e$effect),numeric(1)),
        effect_scale=vapply(l,function(e)chr(e$effectScale),character(1)),
        q_value=vapply(l,function(e)num(e$qValue),numeric(1)),
        stringsAsFactors=FALSE)
    }))
    if(!is.null(plot_effects))
      utils::write.csv(plot_effects,file.path(output,"metabologram_plot_values.csv"),
                       row.names=FALSE,na="NA")
    draw_half <- function(entries,start,end,radius1,radius2,title) {
      if(!length(entries))return(invisible(NULL))
      effects <- vapply(entries,function(e)num(e$effect),numeric(1))
      scales <- vapply(entries,function(e)chr(e$effectScale),character(1))
      cap <- if(all(scales=="log2"))3 else max(1e-9,abs(effects[is.finite(effects)]),na.rm=TRUE)
      n<-length(entries)
      for(i in seq_len(n)){
        a<-seq(start+(end-start)*(i-1)/n,start+(end-start)*i/n,length.out=20)
        xx<-c(radius2*cos(a),rev(radius1*cos(a)))
        yy<-c(radius2*sin(a),rev(radius1*sin(a)))
        polygon(xx,yy,col=color_for(effects[i],cap),border="#FFFFFF",lwd=.7)
        angle<-(start+(end-start)*(i-.5)/n)
        label<-chr(entries[[i]]$labelFr,chr(entries[[i]]$feature))
        if(nchar(label)>15)label<-paste0(substr(label,1,13),"…")
        text(1.2*cos(angle),1.2*sin(angle),labels=label,cex=.47,
             adj=if(cos(angle)<0)c(1,0.5) else c(0,0.5))
      }
      text(0,if(start==pi/2)0.28 else -0.28,title,cex=.8,font=2)
    }
    with_devices("metabologram",11,9,function(){
      op<-par(mar=c(2,2,4,2));on.exit(par(op))
      plot(c(-1.75,1.75),c(-1.5,1.5),type="n",asp=1,axes=FALSE,xlab="",ylab="",
           main="MultiOmics — descriptive circular effects (not pathway significance)")
      draw_half(transcript,pi/2,3*pi/2,.58,.97,"RNA")
      draw_half(metabolites,-pi/2,pi/2,.58,.97,"MET")
      symbols(0,0,circles=.48,inches=FALSE,fg="#AAAAAA",bg="#FFFFFF",add=TRUE)
      text(0,0,"Effects\n(not flux)",cex=.8)
      mtext("Separate effect scales by layer. Grey = unknown. Labels use original features.",side=1,cex=.7)
    },"metabologram","visualizations.metabologram[transcriptomics,metabolomics]")
  }
}
central <- visuals$centralCarbon
if(!is.null(central) && length(central$metabolites)){
  nodes <- c(central$metabolites,central$enzymes)
  nmap <- setNames(nodes,vapply(nodes,function(n)chr(n$id),character(1)))
  layer_limit<-function(entries) {
    measured<-Filter(function(n)is.finite(num(n$measurement$effect)),entries)
    if(!length(measured))return(3)
    scales<-vapply(measured,function(n)chr(n$measurement$effectScale),character(1))
    if(all(scales=="log2"))return(3)
    max(1e-9,vapply(measured,function(n)abs(num(n$measurement$effect)),numeric(1)))
  }
  met_limit<-layer_limit(central$metabolites)
  rna_limit<-layer_limit(central$enzymes)
  node_table <- do.call(rbind,lapply(nodes,function(n){
    m <- n$measurement
    data.frame(id=chr(n$id),label=chr(n$label),
      x=num(n$x),y=num(n$y),measurement_status=chr(m$status,"unmeasured"),
      effect_scale=chr(m$effectScale),
      effect=if(identical(chr(m$status),"measured"))num(m$effect) else NA_real_)
  }))
  utils::write.csv(node_table,file.path(output,"central_carbon_plot_values.csv"),
                   row.names=FALSE,na="NA")
  with_devices("central_carbon",12,8,function(){
    op<-par(mar=c(2,2,3,2));on.exit(par(op))
    plot(c(0,920),c(610,0),ylim=c(610,0),type="n",axes=FALSE,xlab="",ylab="",
         main="Central carbon — descriptive annotation (not metabolic flux)")
    for(e in central$edges) {
      if(length(e)!=2L)next
      from<-nmap[[chr(e[[1L]])]];to<-nmap[[chr(e[[2L]])]]
      if(is.null(from)||is.null(to))next
      arrows(num(from$x),num(from$y),num(to$x),num(to$y),length=.05,col="#B5B9C0")
    }
    for(n in nodes){
      m<-n$measurement
      effect<-if(identical(chr(m$status),"measured"))num(m$effect) else NA_real_
      cap<-if(n$id %in% vapply(central$metabolites,function(k)chr(k$id),character(1)))met_limit else rna_limit
      fill<-if(is.finite(effect))color_for(effect,cap) else "#DCDCDC"
      symbols(num(n$x),num(n$y),circles=18,inches=FALSE,add=TRUE,
              fg="#5B6471",bg=fill)
      text(num(n$x),num(n$y),labels=chr(n$label),cex=.6)
    }
    mtext("Blue/red = measured effects on their declared scale (per layer). Grey = unavailable.",
          side=1,cex=.75)
  },"central_carbon_map","visualizations.centralCarbon[nodes,edges]")
}
if(!length(outputs))stop("No supported figure with archived display-ready data")
manifest<-list(schema_version="1.0.0",figure_mode="R_publication_similar",
  does_not_refit_model=TRUE,pixel_identity_to_Svelte=FALSE,
  data_source="provenance/analysis_complete.json",
  visual_contract=list(
    heatmap_scale="row z-score",heatmap_color_limits=c(-2.5,2.5),
    palette=list(low="#3567B2",neutral="#F7F7F7",high="#CC3E3E",
                 missing_heatmap="#E0E0E0",missing_metabolic_map="#DCDCDC"),
    heatmap_order="original site sample and ranked-feature order",
    effect_scales="retained per omics layer: log2 or source-declared units",
    missing_values="NA remains missing; never coerced to zero",
    no_statistical_refit=TRUE),
  figures=outputs)
writeLines(jsonlite::toJSON(manifest,auto_unbox=TRUE,pretty=TRUE),
           file.path(output,"figure_manifest.json"))
cat("PASS: R created",length(outputs),"figure files from archived display data at",output,"\n")
