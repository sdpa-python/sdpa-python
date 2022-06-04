# SDPA for Python

SDPA for Python (`sdpa-python`) is a Python 3 wrapper for SDPA (SemiDefinite Programming Algorithm). SDPA is a software package for solving general SDPs based on primal-dual interior-point methods with the HRVW/KSH/M search direction [1].

This package is a Python 3 port of SDPAP, the Python 2 based wrapper originally written by **Kenta Kato** provided at the [official SDPA website](http://sdpa.sourceforge.net/download.html). This repository aims to provide Python 3 support for SDPA.

SDPA for Python can be installed by

```bash
pip install sdpa-python
```

For usage documentation or to build from source, please see the [documentation website](https://sdpa-python.github.io).

Owing to its implementation that uses BLAS/LAPACK for numerical linear algebra for dense matrix computation [1], it can potentially be linked against hardware specific BLAS implementations, providing acceleration on various architectures. Currently the binaries available on PyPI are built using OpenBLAS for Windows and Linux, and Apple Accelerate for macOS.

## History

SDPA was officially developed between 1995 and 2012 by **Makoto Yamashita, Katsuki Fujisawa, Masakazu Kojima, Mituhiro Fukuda, Kazuhiro Kobayashi, Kazuhide Nakata, Maho Nakata and Kazushige Goto** [1] [2] [3]. The [official SDPA website](http://sdpa.sourceforge.net/download.html) contains an unmaintained version of SDPA.

SDPAP was written by **Kenta Kato** as a Python 2 interface for SDPA. The [official SDPA website](http://sdpa.sourceforge.net/download.html) also contains an unmaintained version of SDPAP.

This package is a Python 3 port of SDPAP.

## References

If you are using SDPA for Python in your research, please cite SDPA by citing the following papers and book chapters. The [BibTex of the below](https://github.com/sdpa-python/sdpa-python/blob/main/CITATIONS.bib) has been included in the repository.

[1] Makoto Yamashita, Katsuki Fujisawa and Masakazu Kojima, "Implementation and evaluation of SDPA 6.0 (Semidefinite Programming Algorithm 6.0)," *Optimization Methods and Software, vol. 18, no. 4, pp. 491–505*, 2003, doi: [10.1080/1055678031000118482](https://doi.org/10.1080/1055678031000118482).

[2] Makoto Yamashita, Katsuki Fujisawa, Kazuhide Nakata, Maho Nakata, Mituhiro Fukuda, Kazuhiro Kobayashi, and Kazushige Goto, "A high-performance software package for semidefinite programs: SDPA 7," *[Research Report B-460](http://www.optimization-online.org/DB_HTML/2010/01/2531.html) Dept. of Mathematical and Computing Science, Tokyo Institute of Technology, Tokyo, Japan, September*, 2010.

[3] Makoto Yamashita, Katsuki Fujisawa, Mituhiro Fukuda, Kazuhiro Kobayashi, Kazuhide Nakata and Maho Nakata, “Latest Developments in the SDPA Family for Solving Large-Scale SDPs,” in *Handbook on Semidefinite, Conic and Polynomial Optimization, M. F. Anjos and J. B. Lasserre, Eds. Boston, MA: Springer US*, 2012, pp. 687–713. doi: [10.1007/978-1-4614-0769-0_24](https://doi.org/10.1007/978-1-4614-0769-0_24).
