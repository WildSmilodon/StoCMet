# StoCMet

StoCMet is a Fortran implementation of the **stochastic collocation method** for non-intrusive uncertainty quantification. It generates the parameter combinations at which an existing deterministic model must be evaluated and can then use quadrature weights to estimate the expected value, variance, standard deviation, skewness, and kurtosis of the model output.

The deterministic solver is not coupled to StoCMet: run it separately for every row in the generated run list, then arrange its results in one of the supported post-processing formats.

## Features

- Independent uniform random variables using tensor-product Gauss–Legendre quadrature
- Independent normal random variables using tensor-product Gauss–Hermite quadrature
- Independent uniform random variables using a Smolyak/Clenshaw–Curtis sparse grid
- Plain-text run lists that can be consumed by any deterministic solver
- General section/column post-processing (`brain`) and specialized post-processors for the included `heatEx`, `lio`, and `bioheat` examples

## Requirements

- GNU Fortran (`gfortran`)
- GNU Make

The Makefile uses GNU Fortran-specific options, including `-fcray-pointer` and `-fdefault-real-8`.

## Build

From the repository root:

```sh
mkdir -p bin
make -C src
```

This creates `bin/StoCMet`. To remove the executable and object/module files:

```sh
make -C src clean
```

`make -C src install` copies the executable to `~/bin/StoCMet-1.0-o`; create `~/bin` first if it does not exist.

## Quick start

StoCMet reads a file named `scm.inp` from its current working directory. The smallest useful workflow is:

```sh
mkdir example
cd example
```

Create `scm.inp`:

```text
# Number of random variables
1
# Number of collocation points per variable
3
# Distribution
uniform
# name  minimum  maximum
conductivity  0.8  1.2
# Post-processing mode
no
```

Run StoCMet from that directory:

```sh
../bin/StoCMet
```

The generated `scm-runList.txt` contains the run ID and parameter value for each deterministic simulation:

```text
 # number of random variables = 1, number of runs = 3
 # ID conductivity
 1  0.8450806661517034
 2  1.0000000000000000
 3  1.1549193338482966
```

Evaluate the deterministic model once per data row, preserving the run order. If the input selects a post-processing mode other than `no`, place those results in the format expected by that mode and run StoCMet again.

## Input format

The non-comment records in `scm.inp` are positional:

```text
number_of_random_variables
number_of_collocation_points
distribution
one_variable_definition_per_random_variable
post_processing_mode
mode_specific_options, if any
```

Blank lines and lines containing `#` are ignored. Because the parser skips any line that contains `#`, comments must be on their own lines; inline comments are not supported.

### Distributions

`uniform` uses one variable definition per line:

```text
name  minimum  maximum
```

`normal` uses:

```text
name  mean  standard_deviation
```

Both methods create a full tensor-product grid, so the number of deterministic runs is

```text
number_of_runs = number_of_collocation_points ^ number_of_random_variables
```

The quadrature tables currently support these point counts:

| Distribution | Supported point counts |
| --- | --- |
| `uniform` | 2–9, 16, 32, 48, 64 |
| `normal` | 2, 3, 4, 5, 7, 9, 11 |

`sparseuniform` uses the same `name minimum maximum` variable definitions as `uniform`, but interprets the collocation-point field as the Smolyak grid level. It is generally preferable to the full tensor grid when the number of random variables makes `nCP^nRV` too expensive.

Random variables are treated as independent in all three methods.

### Post-processing modes

The record after the variable definitions selects what happens after the run list is written:

| Mode | Purpose | Main output |
| --- | --- | --- |
| `no` (or any unrecognized value) | Generate the run list only | `scm-runList.txt` |
| `brain` | Process selected columns and line ranges in numbered run directories | `ExpV-*`, `Vari-*`, `SDev-*`, `Skew-*`, `Kurt-*` |
| `heatEx` | Process columns from a heat-exchanger `Results.dat` file | `scm-results.dat` |
| `lio` | Process the expected `rptidNNNN.out` and `juridNNNN.out` files | `scm-results.dat` |
| `bioheat` | Process `RUN_ID/all_T_values.txt` files | Five `scm-results-*.dat` files |

These post-processors expect specific directory layouts and file contents. Start with the matching configuration under `run/` and adapt its paths and variable definitions:

- `run/brain/scm.inp` (set its current `no` mode to `brain` to enable post-processing)
- `run/heatEx/scm.inp`
- `run/lio/scm.inp`
- `run/bioheat/scm.inp`

All relative paths are resolved from the directory in which StoCMet is launched. Output files are also written there and may be overwritten on subsequent runs.

## Statistical output

For every processed quantity, StoCMet computes quadrature-weighted estimates of:

- expected value
- variance
- standard deviation
- skewness
- kurtosis (the standardized fourth central moment, not excess kurtosis)

The accuracy depends on the selected quadrature order or sparse-grid level and on how smoothly the deterministic response varies with the uncertain inputs.

## Repository layout

```text
src/             Fortran sources and Makefile
bin/             Compiled executable (after building)
build/           Object and module files (after building)
run/             Example configurations and result data
LICENSE          GNU General Public License v3
```

## License

StoCMet is distributed under the [GNU General Public License v3](LICENSE).

## Author

Jure Ravnik — <jure.ravnik@um.si>
