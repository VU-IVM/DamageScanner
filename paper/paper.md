---
title: 'DamageScanner: A Python package for natural hazard risk assessments'
tags:
  - Python
  - natural hazards
  - damage
  - exposure
  - risk
authors:
  - name: Elco E. Koks
    corresponding: true # (This is how to denote the corresponding author)
    orcid: 0000-0002-4953-4527
    affiliation: 1 # (Multiple affiliations must be quoted)
  - given-names: Hans
    dropping-particle: de
    surname: Moel
    orcid: 0000-0002-6826-1974
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Sadhana Nirandjan
    orcid: 0000-0002-2967-7782
    affiliation: 1
  - given-names: Jens
    dropping-particle: de
    surname: Bruijn
    orcid: 0000-0003-3961-6382
    affiliation: "1,2"
affiliations:
 - name: Water & Climate Risk, Institute for Environmental Studies, Vrije Universiteit Amsterdam, The Netherlands
   index: 1
 - name: International Institute for Applied Systems Analysis (IIASA), Laxenburg, Austria
   index: 2
date: 20 March 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`DamageScanner` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `DamageScanner` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `DamageScanner` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`DamageScanner` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `DamageScanner` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

