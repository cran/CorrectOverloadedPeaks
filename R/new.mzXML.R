new.mzXML = function()
{
  object = list(
    header         = NULL, # required - list    - Path to all the ancestor files (up to the native acquisition file) used to generate the current XML instance document.
    parentFile     = NULL, # required - list    - Path to all the ancestor files (up to the native acquisition file) used to generate the current XML instance document.
    dataProcessing = NULL, # required - list    - Description of any manipulation (from the first conversion to mzXML format until the creation of the current mzXML instance document) applied to the data.
    msInstrument   = NULL, # optional - element - General information about the MS instrument.
    separation     = NULL, # optional - element - Information about the separation technique, if any, used right before the acquisition.
    spotting       = NULL, # optional - element - Acquisition independent properties of a MALDI experiment.
    scan           = vector(mode="list")
  )
  class(object) <- "mzXML"
  return(object)
}