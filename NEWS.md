# flowTime 1.3

- `summarizeFlow` now summarizes all data channels by default
    - now implements the dplyr-based `meanMedianSD` function
    - seems to be a little bit faster 
- `tidyFlow` is a new function that generates a rectangular, tidy 
dataset containing all collected parameters for each cell, including the time 
at which the cell was measured
- default arguments for `steadyState` have been updated to be most general 
