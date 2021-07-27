# flowTime 1.15.1

- bugfix in `summarizeFlow` which was causing numeric columns to be converted to characters

# flowTime 1.15

- `summarizeFlow` now summarizes all data channels by default
    - now implements the dplyr-based `meanMedianSD` function
    - seems to be a little bit faster 
- `tidyFlow` is a new function that generates a rectangular, tidy 
dataset containing all collected parameters for each cell, including the time 
at which the cell was measured
- default arguments for `steadyState` have been updated to be most general 
