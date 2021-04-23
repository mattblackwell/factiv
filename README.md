# factiv: Analyzing Factorial Experiments with Noncompliance

`factiv` is an R package to estimate causal effects in factorial experiments with noncompliance on multiple factors. It implements the methodology of [Blackwell (2017)](https://mattblackwell.org/files/papers/joint-iv.pdf) and [Blackwell and Pashley (2020)][BP2020]. 

To install the package, you can use the `remotes` package:
    
    remotes::install_github("mattblackwell/factiv")


There are two main functions in the package, `iv_factorial`, which  performs superpopulation-based inference, and `iv_finite_factorial`, which performs finite-population inference for the factorial effects using the Fieller method described in [Blackwell and Pashley (2020)][BP2020]. You can use these functions in the following manner:

```r
library(factiv)
data(newhaven)

superpop <- iv_factorial(turnout_98 ~ inperson + phone | inperson_rand + 
  phone_rand, data = newhaven)

finite <- iv_finite_factorial(formula = turnout_98 ~ inperson + phone |
  inperson_rand + phone_rand, data = newhaven)
```


[BP2020]: https://mattblackwell.org/files/papers/factorial-iv.pdf
