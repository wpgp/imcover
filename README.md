
<!-- README.md is generated from README.Rmd. Please edit that file -->

# imcover

<!-- badges: start -->
<!-- badges: end -->

## Testing package of immunisation coverage modelling in R.

The `imcover` package is not yet on CRAN and can be installed from
Github.

Because this is still in a private repository, you must first create an
authentication token. To do this, log in to github.com and got to the
Developer settings: <https://github.com/setting/tokens>.

1.  Select the button to “Generate new token” at the top of the page.
2.  Optionally type a note to remind yourself what this token is.
3.  Set an expiration date
4.  Check the box for **repo** scope to give access to private
    repositories
5.  Copy your token somewhere safe and don’t share it! You will need it
    for the next step

Now, in `R` you can install the package. In the following command insert
your token as a string:

``` r
devtools::install_github('wpgp/imcover@main', auth_token = 'YOUR TOKEN')
```

The build process takes some time because it is compiling the C++ code
for the `Stan` models. It may also ask you to install some additional
dependencies.
