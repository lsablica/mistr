<div align="center">
  <img src="vignettes/figures/sticker.png" alt="mistr Logo" width="200"/>

  # mistr

  [![R](https://img.shields.io/badge/R-%23E67E22.svg?&logo=R&logoColor=white)](https://www.r-project.org/)
  [![CRAN Status](https://www.r-pkg.org/badges/version/mistr)](https://cran.r-project.org/package=mistr)
  [![License: GPL-3.0](https://img.shields.io/badge/License-GPL%203.0-blue.svg)](https://opensource.org/licenses/GPL-3.0)

  *A computational framework for mixture and composite distributions.*

  [Key Features](#key-features) •
  [Project Overview](#project-overview) •
  [Citation](#citation)
</div>

---

## Key Features  

📦 **Flexible Distribution Modeling**
- Comprehensive support for standard, composite, and mixture distributions.
- Built-in symbolic differentiation engine for transformations and inverse transformations, simplifying complex operations on random variables.

📊 **Advanced Statistical Functions**
- Seamlessly evaluates PDFs, CDFs, quantiles, and random samples.
- Includes monotonic transformations and mixtures with custom weights.

🧮 **Efficient Mathematical Framework**
- Automatic construction of required equations for transformed random variables.
- Supports truncated and hierarchical distributions for specialized modeling.

📖 **Extensive Documentation**
- Detailed vignettes and examples, including tutorials on risk measures, parameter estimation, and real-world applications.

## Project Overview

`mistr` is a cutting-edge R package that enables users to define, manipulate, and analyze univariate probability distributions. Whether modeling complex financial data or designing custom distributions, `mistr` provides the tools needed for a wide range of statistical tasks.

Highlights include:
- **Composite and Mixture Models**: Easily construct and analyze hierarchical and truncated distributions.
- **Symbolic Differentiation**: Built-in engine to derive and manage equations for transformations and inversions of random variables.
- **Versatility**: From academic research to real-world applications, `mistr` simplifies handling distributions.

### Applications
- **Financial Mathematics**: Model heavy-tailed data for asset returns, risk management, and loss modeling.
- **Actuarial Science**: Fit and evaluate composite models for insurance claim distributions.
- **Custom Distributions**: Create and analyze tailored distributions for experimental data.

## Citation

If you use `mistr` in your research, please cite:

```bibtex
@article{mistr2020,
  title = "{mistr: A Computational Framework for Mixture and Composite Distributions}",
  author = {Lukas Sablica and Kurt Hornik},
  journal = "{The R Journal}",
  year = "2020",
  volume = "12",
  number = "1",
  pages = "283--299",
  doi = "10.32614/RJ-2020-003",
  url = "https://journal.r-project.org/archive/2020/RJ-2020-003/index.html"
}
```



