#'@title Write a standard mzXML.
#'
#'@description
#'\code{write.mzXML} is a copy from the caMassClass package which is no longer actively on CRAN
#'
#'@details
#'The main task of \code{\link{read.mzXML}} and write.mzXML functions is to extract and save scan data of mzXML files. 
#'In addition attempt is made to keep all other sections of mzXML file as unparsed XML code, so the data can be extracted latter or saved into new mzXML files.
#'Those unparsed sections are stored as XML text.
#'
#'@param mzXML The generic mzXML object
#'@param filename A mzXML filename to export to.
#'@param precision Either '32' or '64' byte.
#'
#'@return
#'Nothing. mzXML object is exported to filename. If it was imported using \code{\link{read.mzXML}} it should contain all previous fields.
#'
#'@export
#'
#'@import bitops
#'
#'@importFrom utils packageDescription
#'
write.mzXML = function(mzXML, filename, precision=c('32', '64'))
{
#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
  # helper functions
  Paste <- function(...) paste(..., sep="", collapse="")
  
  fprintf <- function(fp, level, ..., append=TRUE)
  {
    x = paste(..., sep="")
    if (length(x)==0 || is.null(x)) return(NULL)
    spaces = if (level>0) Paste(rep("  ", level)) else ""
    x = gsub("'", "\"", x)
    cat(spaces, x, file=fp, sep="")
    NULL
  }
  
  # [Modification_by_JaLi:] 
  # this is a copy from the caTools function as this package is about to be archived (07/2018)
  base64encode <- function (x, size = NA, endian = .Platform$endian) 
  {
    if ((typeof(x) != "character") & (typeof(x) != "raw")) 
      x = writeBin(x, raw(), size = size, endian = endian)
    if ((typeof(x) == "character") & (typeof(x) != "raw")) {
      nlen <- nchar(x)
      x = writeBin(x, raw(), size = size, endian = endian)
      length(x) <- nlen
    }
    x = as.integer(x)
    ndByte = length(x)
    nBlock = ceiling(ndByte/3)
    neByte = 4 * nBlock
    if (ndByte < 3 * nBlock) 
      x[(ndByte + 1):(3 * nBlock)] = 0
    dim(x) = c(3, nBlock)
    y = matrix(as.integer(0), 4, nBlock)
    y[1, ] = bitops::bitShiftR(x[1, ], 2)
    y[2, ] = bitops::bitOr(bitops::bitShiftL(x[1, ], 4), bitops::bitShiftR(x[2, ], 4))
    y[3, ] = bitops::bitOr(bitops::bitShiftL(x[2, ], 2), bitops::bitShiftR(x[3, ], 6))
    y[4, ] = x[3, ]
    y = bitops::bitAnd(y, 63)
    alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
    alpha = strsplit(alpha, NULL)[[1]]
    z = alpha[y + 1]
    npbytes = 3 * nBlock - ndByte
    if (npbytes > 0) 
      z[(neByte - npbytes + 1):neByte] = "="
    z = paste(z, collapse = "")
    return(z)
  }
  # done with local functions
  
  precision = match.arg(precision)
  if (!is.character(filename)) stop("read.mzXML: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string
  fp  = file(filename, "w")
  
  if (is.null(mzXML))
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");
  if (attr(mzXML, "class")!="mzXML")
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");
  
  #------------------------------------
  # Fill-in required sections if empty
  #------------------------------------
  if (is.null(mzXML$header)) {
    Str = "http://sashimi.sourceforge.net/schema_revision/mzXML_2.1"
    mzXML$header = Paste( "<mzXML xmlns='",Str,"'\n  ",
                          "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n  ",
                          "xsi:schemaLocation='", Str, " ", Str, "/mzXML_idx_2.1.xsd'>\n")
  }
  if (is.null(mzXML$parentFile)) {
    mzXML$parentFile = Paste( "    <parentFile filename='file://unknown' ",
                              "fileType='RAWData' fileSha1='0000000000000000000000000000000000000000'/>\n")
  }
  if (is.null(mzXML$dataProcessing)) {
    Version = utils::packageDescription("caMassClass")$Version
    Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    mzXML$dataProcessing = Paste("    <dataProcessing>\n",
                                 "      <software type='processing' name='cran.r-project.org/caMassClass' ",
                                 "version='",Version,"' completionTime='",Time,"'/>\n    </dataProcessing>")
  }
  
  #-----------------------------
  # Write beggining of file
  #-----------------------------
  fprintf(fp, 0, "<?xml version='1.0' encoding='ISO-8859-1'?>\n", append=FALSE)
  fprintf(fp, 0, mzXML$header)
  fprintf(fp, 1, "<msRun scanCount='",length(mzXML$scan),"'>\n")
  fprintf(fp, 0, mzXML$parentFile)
  fprintf(fp, 0, mzXML$msInstrument)
  fprintf(fp, 0, mzXML$dataProcessing)
  fprintf(fp, 0, mzXML$separation)
  fprintf(fp, 0, mzXML$spotting)
  
  #---------------------------------
  # Write scan Section
  #---------------------------------
  indexScan = "  <index name='scan'>\n"
  n   = length(mzXML$scan)
  Num = integer(n)
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$num
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$parentNum
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$msLevel
  Num = 1-diff(c(Num,1)) # number of </scan> after each scan
  size = (if (precision=="32") 4 else 8)
  for(i in 1:n) if (!is.null(mzXML$scan[[i]])) {
    indexScan = Paste(indexScan, "    <offset id='",mzXML$scan[[i]]$num,"'>",
                      seek(fp),"</offset>\n")
    mass  = mzXML$scan[[i]]$mass
    peaks = mzXML$scan[[i]]$peaks
    stopifnot(length(mass)==length(peaks))
    ScanHeader = Paste("<scan num='"   , mzXML$scan[[i]]$num,
                       "' msLevel='"   , mzXML$scan[[i]]$msLevel, 
                       "' peaksCount='", length(peaks),
                       "'  lowMz='"     , min(mass),
                       "'  highMz='"    , max(mass),
                       "'", mzXML$scan[[i]]$scanAttr, ">\n")
    fprintf(fp, 2, ScanHeader)
    fprintf(fp, 3, mzXML$scan[[i]]$scanOrigin)
    fprintf(fp, 0, mzXML$scan[[i]]$precursorMz)
    fprintf(fp, 0, mzXML$scan[[i]]$maldi)
    p = as.vector(rbind(mass,peaks))
    fprintf(fp, 3, Paste("<peaks precision='",precision, "' byteOrder='network' pairOrder='m/z-int'>", base64encode(p, endian="big", size=size), "</peaks>\n"))
    fprintf(fp, 0, mzXML$scan[[i]]$nameValue)
    if(Num[i]) for (j in 1:Num[i]) fprintf(fp, 2, "</scan>\n")
  }
  indexScan = Paste(indexScan, "  </index>\n")
  
  
  #---------------------------------
  # Write end of file
  #---------------------------------
  fprintf(fp, 1, "</msRun>\n")
  n = seek(fp)
  fprintf(fp, 0, indexScan)
  fprintf(fp, 1, "<indexOffset>",n,"</indexOffset>\n")
  cat("  <sha1>", file=fp, sep="")
  n = seek(fp)
  close(fp)
  sha1 = digest::digest(filename, algo="sha1", file=TRUE, length=n)
  cat(sha1, "</sha1>\n</mzXML>\n", file=filename, append=TRUE, sep="")
  invisible(NULL)
}
