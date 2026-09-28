```markdown
---
title: 'Fit-abundace: A Python package for fit chemical abundance gradients'
tags:
  - Python
  - astronomy
  - abundance gradients
  - spiral galaxies
authors:
  - name: André Felie de Siqueira Cardoso
    orcid: 0000-0003-1097-3247
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Author Oscar Cavichia
    orcid: 0000-0002-7103-8036
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: PPGCosmo, Núcleo Cosmo - Ufes, CCE, Universidade Federal do Espírito Santo, Av. Fernando Ferrari, 540, 29075-910 Vitória, ES, Brazil
   index: 1
   ror: 00hx57361
 - name: Instituto de Física e Química, Universidade Federal de Itajubá, Av. BPS, 1303, 37500-903 Itajubá-MG, Brazil
   index: 2
date: 28 September 2026
bibliography: paper.bib

# Optional fields for papers that are part of a joint submission.
# For example, submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#
# If you are not making a joint submission you should remove these lines.
#
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

`fit_abundance` is a Python package for the automated analysis of radial
oxygen and nitrogen abundance gradients in spiral galaxies using spectroscopic
observations of H II regions. It combines abundance determination, H II-region
selection, gradient fitting, statistical model selection, and diagnostic
visualization in a single pipeline.

# Statement of need

Chemical abundance gradients provide important information about the chemical
evolution of galaxies, but their analysis involves several sequential data-processing
and modeling steps. `fit_abundance` provides a useful implementation of this pipeline,
reducing the need for repeated ad hoc analysis codes.

# State of the field                                                                                                                  

Several astronomical software packages provide tools for individual tasks such
as spectral analysis, statistical fitting, or astronomical data processing.
`fit_abundance` focuses specifically on the automated analysis of radial
oxygen abundance gradients from H II region spectroscopy, integrating these steps
into a single analysis pipeline.

# Software design

The package is organized into modules for galactocentric distance calculation,
abundance determination, H II region selection, abundance calibrators, gradient models,
and visualization. The main `fit_final` function combines these components into an
automated workflow while allowing the user to specify calibrators, selection criteria,
and fitting options.

# Research impact statement

`fit_abundance` was developed from the computational methodology used to analyze oxygen
abundance gradients in 147 spiral galaxies in Cardoso et al. (2025) [@Cardoso2025].
The software was subsequently used to perform additional scientific analyses, providing a
reusable implementation of the methodology for future studies of galaxy chemical abundance gradients.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

```
