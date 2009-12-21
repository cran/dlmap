`summary.dlmap` <- 
function(object, ...)
{
  cat(" Summary of input data: \n")
  summary(object$input, ...)
  cat(" Summary of final results: \n")
  print(object$Summary)
}
 
