

#' A gate for the set of all yeast cells
#'
#' Typically set in FSC.A by SSC.A space to exclued any debris
#' @format formal class polygonGate

'yeastGate'

#' A gate for the set of all diploid singlet yeast cells
#'
#' Typically set in FSC.A by FSC.H space
#' Diploids are typically 5um x 6um ellipsoids while haploids are typically
#' 4um x 4um spheroids. As a result, diploids are longer and you get a
#' larger 'area/volume'.
#'
#' @format formal class polygonGate
'dipsingletGate'

#' A gate for the set of all diploid doublets
#'
#' @format formal class polygonGate
'dipdoubletGate'

#' A gate for the set of all haploid singlets
#' @format formal class polygonGate
'hapsingletGate'

#' A gate for the set of all haploid doublets
#' @format formal class polygonGate
'hapdoubletGate'
