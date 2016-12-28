#' A gate for the set of all yeast cells
#'
#' Typically set in FSC.A by SSC.A space to exclued any debris
#'     FSC.A  SSC.A
#'    400000  10000
#'     50000  15000
#'     40000  80000
#'   1250000 600000
#'   2000000 500000
#' @format formal class polygonGate
#'
'yeastGate'

#' A gate for the set of all diploid singlet yeast cells
#'
#' Typically set in FSC.A by FSC.H space
#' Diploids are typically 5um x 6um ellipsoids while haploids are typically
#' 4um x 4um spheroids. As a result, diploids are longer and you get a
#' larger 'area/volume' FSC.A
#'     FSC.A   FSC.H
#'     80000   80000
#'     70000  100000
#'    800000 1000000
#'   1250000 1250000
#'   1000000 1000000
#' @format formal class polygonGate
#'
'dipsingletGate'
