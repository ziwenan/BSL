# BSL 2.0.0

* Second major release. Please note some arguments have been renamed.
* Added a new method ("semiBSL") in the main function "bsl".
* Now bsl returns an S4 class object instead of S3. The plot method gives an option for whether ggplot2 or R graphics (the default) will be used.
* The simulation function and summary statistics function do not need to have a list as the second argument now.
* The new bsl function uses "parallelArgs" instead of "parallel_packages", where users can specify all arguments supported by "foreach".
* Code update for three examples: MA(2), multivariate G & K, and cell biology.
