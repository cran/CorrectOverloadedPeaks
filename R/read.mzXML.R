#'@title Read a standard mzXML
#'
#'@description
#'\code{read.mzXML} is a copy from the caMassClass package which is no longer actively on CRAN.
#'
#'@details
#'The main task of read.mzXML and \code{\link{write.mzXML}} functions is to extract and save scan data of mzXML files. 
#'In addition attempt is made to keep all other sections of mzXML file as unparsed XML code, so the data can be extracted latter or saved into new mzXML files.
#'Those unparsed sections are stored as XML text.
#'
#'@param filename A mzXML file as exported by vendor software.
#'
#'@return
#'A generic R object of class mzXML.
#'
#'@export
#'
read.mzXML = function(filename)
{
#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
  Paste = function(...) paste(..., sep="", collapse="")
  
  strtrunc = function(Str,Sub) {
    lp = attr(regexpr(paste(".*",Sub,sep=""),Str),'match.length')
    return( substring(Str, 1, lp) )
    #y = unlist(strsplit(Str,Sub)) # other way of doing it
    #return( paste(y[-length(y)], sub,  sep="", collapse="") )
  }
  
  fregexpr = function(pattern, filename)
  { # similar to gregexpr but operating on files not strings
    buf.size=1024
    n  = file.info(filename)$size
    pos = NULL
    fp = file(filename, "rb")
    for (d in seq(1,n,by=buf.size)) {
      m = if (n-d>buf.size) buf.size else n-d
      p = gregexpr(pattern, readChar(fp, m))[[1]]
      if(p[1]>0) pos=c(pos, p+d-1)
    }
    close(fp)
    if (is.null(pos)) pos=-1
    return (pos)
  }
  
  #-------------------------------
  # define XML handler function
  #-------------------------------
  mzXMLhandlers <- function()
  {
    #---------------------------------------------------------
    # local variables
    #---------------------------------------------------------
    obj      = new.mzXML() # create new mzXML object
    iScan    = 0
    ParentID = vector(mode="integer")
    sha1     = vector(mode="list", length=2) # optional - element - sha-1 sums
    sha1[1] <- sha1[2] <- 0
    # Optional attributes that might come with a scan that will be stored
    OptScanAttr = c("polarity", "scanType", "centroided", "deisotoped", 
                    "chargeDeconvoluted", "retentionTime", "ionisationEnergy", 
                    "collisionEnergy", "cidGasPressure", "totIonCurrent") 
    
    #-------------------------------
    # local functions
    #-------------------------------
    ToString = function(x, indent = 1)
    { # converts content of a node to a string
      if (is.null(x)) return(NULL);
      spaces = if (indent>0) Paste(rep("  ", indent)) else ""
      Name = XML::xmlName(x, TRUE)
      val  = XML::xmlValue(x)
      if (Name=="text") return( Paste(spaces, val, "\n") )
      if (!is.null(XML::xmlAttrs(x))) {
        att = paste(names(XML::xmlAttrs(x)), paste("\"", XML::xmlAttrs(x),
                                              "\"", sep = ""), sep = "=", collapse = " ")
        att = paste(" ", att, sep="")
      } else att = ""
      chl = ""
      for (i in XML::xmlChildren(x)) chl = Paste(chl, ToString(i, indent+1))
      if (chl=="") Str = Paste(spaces, "<" , Name, att, "/>\n")
      else Str = Paste(spaces, "<" , Name, att, ">\n", chl, spaces, "</", Name, ">\n")
      return(Str)
    }
    
    CatNodes = function(x,Name, indent = 2)
    { # concatinate strings of several nodes
      Str=NULL
      for (y in XML::xmlElementsByTagName(x, Name))
        Str = paste(Str, ToString(y,indent), sep="")
      return(Str)
    }
    
    read.mzXML.scan = function(x)
    { # process scan section of mzXML file
      if (is.null(x)) return(NULL)
      if (XML::xmlName(x) != "scan") return(NULL)
      scanOrigin <- precursorMz <- nameValue <- maldi <- mass <- peaks <- NULL
      att         = XML::xmlAttrs(x)
      num         = as.integer(att["num"])
      msLevel     = as.integer(att["msLevel"])
      peaksCount  = as.integer(att["peaksCount"]) # Total number of m/z-intensity pairs in the scan            
      msk         = names(att) %in% OptScanAttr
      if (sum(msk)==0) {
        scanAttr = ""
      } else {
        scanAttr = paste( names(att[msk]), paste("\"", att[msk], "\"", sep = ""), sep = "=", collapse = " ")
        scanAttr = paste(" ", scanAttr, sep="")
      }
      maldi       = ToString(x[["maldi"]])
      scanOrigin  = CatNodes(x, "scanOrigin", 3)
      nameValue   = CatNodes(x, "nameValue", 3)
      precursorMz = CatNodes(x, "precursorMz", 3)
      precursorMz = gsub("\n      " , " " , precursorMz)
      for (y in XML::xmlElementsByTagName(x, "scan")) ParentID[as.integer(XML::xmlAttrs(y)["num"])] <<- num
      y           = x[["peaks"]]
      att         = XML::xmlAttrs(y)
      peaks       = XML::xmlValue(y) # This is the actual data encoded using base64
      precision   = att["precision"] # nr of bits used by each component (32 or 64)
      byteOrder   = att["byteOrder"] # Byte order of the encoded binary information (must be network)
      # [Modification_by_JaLi:] ProteoWizard exported "contentType" instead of "pairOrder" field
      if ("contentType" %in% names(att)) pairOrder   = att["contentType"] # Order of the m/z intensity pairs (must be m/z-int)
      if ("pairOrder" %in% names(att)) pairOrder   = att["pairOrder"] # Order of the m/z intensity pairs (must be m/z-int)
      endian      = ifelse(byteOrder=="network", "big", "little")
      if (precision=="32") size=4
      else if(precision=="64") size=8
      else stop("read.mzXML.scan: incorrect precision attribute of peaks field")
      if (is.na(pairOrder) | pairOrder!="m/z-int")
        warning("read.mzXML.scan: incorrect pairOrder attribute of peaks field")
      if (peaksCount>0) {
        p = caTools::base64decode(peaks, "double", endian=endian, size=size)
        np = length(p) %/% 2
        if (np != peaksCount)
          warning("read.mzXML.scan: incorrect 'peakCount' attribute of 'peaks' field: expected ",
                  peaksCount, ", found ", np, "  ",(3*((nchar(peaks)*size)/4))/2, " (scan #",num,")")
        dim(p)=c(2, np)
        mass =p[1,]
        peaks=p[2,]
      }
      #x$children=NULL; # needed to capture the header
      #header <<- toString(x)
      return( list(mass=mass, peaks=peaks, num=num, parentNum=num,
                   msLevel=msLevel, scanAttr=scanAttr, maldi=maldi,
                   scanOrigin=scanOrigin, precursorMz=precursorMz, nameValue=nameValue) )
    }
    
    #---------------------------------------------------------
    # the instructions how to parse each section of mzXML file
    #---------------------------------------------------------
    list(
      mzXML  = function(x, ...) {
        y = x[["sha1"]]
        sha1[1]    <<- if (!is.null(y)) XML::xmlValue(y) else 0
        x$children =  NULL
        obj$header <<- toString(x,terminate=FALSE)
        NULL
      },
      
      msRun = function(x, ...) {
        y = x[["sha1"]]
        sha1[2]            <<- if (!is.null(y)) XML::xmlValue(y) else 0
        obj$msInstrument   <<- ToString(x[["msInstrument"]],2)
        obj$separation     <<- ToString(x[["separation"]],2)
        obj$spotting       <<- ToString(x[["spotting"]],2)
        obj$parentFile     <<- CatNodes(x,"parentFile")
        obj$dataProcessing <<- CatNodes(x,"dataProcessing")
        NULL
      },
      
      scan  = function(x, ...) {
        iScan <<- iScan+1
        obj$scan[[iScan]] <<- read.mzXML.scan(x)
        x$children=NULL
        x
      },
      
      data = function() {
        if (is.null(obj$header)) NULL
        else list(mzXML=obj, ParentID=ParentID, sha1=sha1)
      }
    ) #end of list of handler functions
  } # done with local functions
  
  #---------------------------------
  # begining of read.mzXML function
  #---------------------------------
  if (!is.character(filename)) stop("read.mzXML: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string
  
  sha1File = digest::digest(filename, algo="sha1", file=TRUE)
  x = XML::xmlTreeParse(file=filename, handlers=mzXMLhandlers(),
                        addAttributeNamespaces=TRUE) $ data()
  if (is.null(x)) # is this file a mzXML file ?
    stop("read.mzXML: This is not mzXML file");
  mzXML    = x$mzXML
  sha1Read = x$sha1
  
  # sort scans into correct order; find parent numbers of recursive nodes
  n = length(mzXML$scan)
  NumID = integer(n)
  for (i in 1:n) {
    NumID[i] = mzXML$scan[[i]]$num
    mzXML$scan[[i]]$scanOrigin = paste("<scanOrigin parentFileID='",sha1File,
                                       "' num='",NumID[i],"'/>\n", sep="");
  }
  mzXML$scan = mzXML$scan[ order(NumID) ]
  for (i in 1:n)
    if (!is.na(x$ParentID[i])) mzXML$scan[[i]]$parentNum = x$ParentID[i]
  else x$ParentID[i] = mzXML$scan[[i]]$parentNum
  mzXML$scan = mzXML$scan[ order(x$ParentID) ]
  
  ## read sha1 section
  n = sum(as.integer(lapply(sha1Read, is.character))) # how many sha1 were found
  if( n>0 ) {
    ## sha1 - sha-1 sum for this file (from the beginning of the file up to
    ## (and including) the opening tag of sha1
    if (is.null(sha1Read[[1]])) sha1Read[[1]]=sha1Read[[2]]
    sha1Pos = fregexpr("<sha1>", filename) + 6 # 6 = length("<sha1>")
    for(i in n) { # multiple sha1 sections are possible
      sha1Calc = digest::digest(filename, algo="sha1", file=TRUE, length=sha1Pos[i]-1)
      if (sha1Read[[i]]!=sha1Calc)
        warning("Stored and calculated Sha-1 sums do not match (stored '",
                sha1Read[[i]],"'; calculated '", sha1Calc,"')")
    }
  }
  
  # strip mzXML terminator from header section
  mzXML$header = gsub("/>", ">\n", mzXML$header)
  mzXML$header = gsub("^ +", "", mzXML$header)
  # Remove incorrect "-quotes inserted in 2.10.0
  mzXML$header = gsub("[\u0093\u0094\u201C\u201D]", '"', mzXML$header)
  # add info about parent file (the file we just read)
  mzXML$parentFile = Paste(mzXML$parentFile, "    <parentFile filename='file://",
                           filename, "' fileType='processedData' fileSha1='", sha1File, "'/>\n")
  return( mzXML )
}