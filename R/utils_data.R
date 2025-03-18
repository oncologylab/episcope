.make_bed = function(bed, op_dir = tempdir(), up = 2500, down = 2500, tss = "center", for_profile = FALSE){
  #bwtool tool requires only three columns

  tss = match.arg(arg = tss, choices = c("start", "end", "center"))

  if(!dir.exists(paths = op_dir)){
    dir.create(path = op_dir, showWarnings = FALSE, recursive = TRUE)
  }

  temp_op_bed = tempfile(pattern = "profileplot", tmpdir = op_dir, fileext = ".bed")

  if(is.data.frame(bed)){
    bed = data.table::as.data.table(x = bed)
    #data.table::setDT(x = bed)
    colnames(bed)[1:3] = c("chr", "start", "end")
    bed[, chr := as.character(chr)]
    bed[, start := as.numeric(as.character(start))]
    bed[, end := as.numeric(as.character(end))]
  }else if(file.exists(bed)){
    bed = data.table::fread(file = bed, select = list(character = 1, numeric = c(2, 3)), col.names = c("chr", "start", "end"))
    bed = bed[,.(chr, start, end)]
  }

  if(!for_profile){
    if(tss == "center"){
      bed[, focal_point := as.integer(apply(bed[,2:3], 1, mean))]
      bed[, bed_start := focal_point-up]
      bed[, bed_end := focal_point+down]
    }else if(tss == "start"){
      bed[, bed_start := start-up]
      bed[, bed_end := start+down]
    }else{
      bed[, bed_start := end-up]
      bed[, bed_end := end+down]
    }
    bed = bed[,.(chr, bed_start, bed_end)]
    data.table::setkey(x = bed, chr, bed_start, bed_end)
  }

  data.table::fwrite(x = bed, file = temp_op_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  return(temp_op_bed)
}
