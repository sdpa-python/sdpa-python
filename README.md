# SDPA for Python

`sdpa-python` is a Python 3 wrapper for SDPA (SemiDefinite Programming Algorithm). SDPA is a software package for solving general SDPs based on primal-dual interior-point methods with the HRVW/KSH/M search direction [1].

This package is a Python 3 port of SDPA-P, the Python 2 based wrapper originally written by **Kenta KATO** provided at the [official SDPA website](http://sdpa.sourceforge.net/download.html). This repository aims to provide Python 3 support for SDPA. It will be accompanied with documentation to allow users to link SDPA against modern BLAS libraries (i.e. Intel MKL, Apple Accelerate and NVIDIA cuBLAS).

For installation instructions, please see the [documentation website](https://sdpa-python.github.io).

## History

SDPA was officially developed between 1995 and 2012 by **Makoto Yamashita, Katsuki Fujisawa, Masakazu Kojima, Mituhiro Fukuda, Kazuhiro Kobayashi, Kazuhide Nakata and Maho Nakata** [1] [2] [3]. The [official SDPA website](http://sdpa.sourceforge.net/download.html) contains an unmaintained version of SDPA.

Owing to it's implementation that uses LAPACK for numerical linear algebra for dense matrix computation [1], it can be linked against modern LAPACK implementations, providing performance gains on a variety of architectures.

## References

If you are using SDPA for Python in your research, please cite SDPA by citing the following papers and book chapters. The BibTex has been provided in `CITATIONS.bib`.

[1] Makoto Yamashita, Katsuki Fujisawa and Masakazu Kojima, "Implementation and evaluation of SDPA 6.0 (Semidefinite Programming Algorithm 6.0)," *Optimization Methods and Software, vol. 18, no. 4, pp. 491–505*, 2003, doi: [10.1080/1055678031000118482](https://doi.org/10.1080/1055678031000118482).

[2] Makoto Yamashita, Katsuki Fujisawa, Kazuhide Nakata, Maho Nakata, Mituhiro Fukuda, Kazuhiro Kobayashi, and Kazushige Goto, "A high-performance software package for semidefinite programs: SDPA 7," *[Research Report B-460](http://www.optimization-online.org/DB_HTML/2010/01/2531.html) Dept. of Mathematical and Computing Science, Tokyo Institute of Technology, Tokyo, Japan, September*, 2010.

[3] Makoto Yamashita, Katsuki Fujisawa, Mituhiro Fukuda, Kazuhiro Kobayashi, Kazuhide Nakata and Maho Nakata, “Latest Developments in the SDPA Family for Solving Large-Scale SDPs,” in *Handbook on Semidefinite, Conic and Polynomial Optimization, M. F. Anjos and J. B. Lasserre, Eds. Boston, MA: Springer US*, 2012, pp. 687–713. doi: [10.1007/978-1-4614-0769-0_24](https://doi.org/10.1007/978-1-4614-0769-0_24).
