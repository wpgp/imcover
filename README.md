
<!-- README.md is generated from README.Rmd. Please edit that file -->

# imcover: Immunisation coverage modelling in R

<!-- badges: start -->
<!-- badges: end -->

The goal of the `imcover` package is to provide access to data download,
processing and modelling tools and to support a Bayesian statistical
modelling approach to generate national estimates of immunization
coverage from multiple time series of data on coverage. `imcover` is
built as part of the open-source statistical computing and modelling
language, `R` (<https://cran.r-project.org/>). This package is designed
to support broadly replicable and reproducible analyses of immunization
coverage.

## Statistical modelling software

The core of `imcover` is the functionality to fit a Bayesian statistical
model of multiple time series. The sources of coverage data (in this
example administrative, official and surveys) are taken as multiple,
partial estimates of the true, unobserved immunization coverage in a
country. A Bayesian estimation approach allows us to incorporate these
multiple datasets, place prior beliefs on which sources are more
reliable, share information between countries, and to quantify
uncertainty in our estimate of the latent immunization coverage.

`imcover` provides an interface to [`Stan`
(https://mc-stan.org/)](https://mc-stan.org/) for statistical
computation. This means that, in addition to `imcover`, many of the
tools for assessing model performance and visualizing results from
`Stan` will work for `imcover` results. However, users must have `Stan`
installed and linked with `R` in order to use `imcover`

To install `Stan`, please follow the instructions for your operating
system described here:
<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

## Installing `imcover`

The `imcover` package is not yet on CRAN and can be installed from
Github.

Because this is still in a private repository, you must first create an
authentication token. To do this, log in to github.com and got to the
Developer settings: <https://github.com/setting/tokens>.

1.  Select the button to “Generate new token” at the top of the page.
2.  Optionally type a note to remind yourself what this token is.
3.  Set an expiration date.
4.  Check the box for **repo** scope to give access to private
    repositories
5.  Copy your token somewhere safe and don’t share it! You will need it
    for the next step.

Now, in `R` you can install the package. In the following command insert
your token as a string:

``` r
devtools::install_github('wpgp/imcover@main', auth_token = 'YOUR TOKEN')
```

The build process takes some time because it is compiling the C++ code
for the `Stan` models. It may also ask you to install some additional
dependencies.

## Contributions

Feedback and contributions are welcome. Please raise or respond to an
issue, or create a new branch to develop a feature/modification and
submit a pull request.
