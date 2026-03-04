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
date: 1 March 2026
bibliography: paper.bib

---
# Summary

Due to a global increase in the severity and frequency of a wide-range of natural hazards, easy-to-use tools are required to properly understand the potential impacts from such events. This could either be a local rainfall event, to a large-scale earthquake event. While the magnitude and intensity of such events differ greatly, the computational workflow to estimate those impact do not. As such, the `DamageScanner` has been developed a simply and ready-to-use framework to allow exposure, damage and risk assesments either through user-defined input, or simply using public data as a starting point. 

# Statement of need

`DamageScanner` is a Python package to perform exposure, damage and risk assessments for natural hazards and climate extremes. It allows users to easily perform a risk assessment for both vector and raster data through a single computational framework.     

The `DamageScanner` is designed to be used by both researchers, risk modellers, and by students in introductory courses on geospatial analysis and natural hazard risk modelling. It has already been used in a number of scientific publications [@Koksetal2019 ; @KoksHaer2020], European projects [@MIRACA:2026] and integrated into other python packages (e.g., `GEB`).   

# State of the field                                                                                                                  

Several tools exist for natural hazard risk assessments, of which the `climada` framework [@Aznar2019; @Bresch2021] is most well-known. The `climada` framework is a widely-used Python ecosystem to perform large-scale risk assessment for various hazards. It allows for a multitude use cases and contains a wide range of submodules to tailor specific needs within the natural hazard risk domain [@Muhlhofer2023 ; @Riedel2024]. All of which is maintained by dedicated softwareengineers. While the `climada` framework is commendable and is used widely, the large ecosystem and the many dependencies within the several submodules also make it less flexible for integration within our own use cases. 

As such, the `DamageScanner` was built rather than contributing to existing projects for several reasons. Firstly, the `DamageScanner` has been developed organically to tailor specific use cases within a multitude of projects. For example, the original `DamageScanner` was built to perform raster-based damage asssessment, mostly using land use information, in combination with flood hazard maps [@DeMoel]. 

Secondly, with the rise of publicly-available object-based data, in particular due to OpenStreetMap (OSM), a need arose to efficiently perform object-based damage calculations [@VanGinkel2023; @Koksetal2019]. These object-based damage and risk assessment are, however, much more computational intensive. The `DamageScanner` fills a gap hhere, to provide an open-access approach to easily perform such risk assessment on top of object information.

Thirdly, geospatial information of natural hazards and climate extremes are stored in a range of different formats. Ranging from GeoTiffs and netCDFs, to vector-based representations of earthquake intensities and flood depths. To model the impacts across those hazards within the same modelling framework, a tool was needed that can readily read and use different storage formats. 

# Software design

The `DamageScanner`'s core design philosophy is focused on simplicity. 

# Research impact statement

The `DamageScanner` has demonstrated significant research impact and grown both its user base
and contributor community since its initial release in 2019. 

While the `DamageScanner` started as a tool primarily to support the core developer's
research, it has expanded organically to support a range of applications 


# AI usage disclosure

Generative AI tools were used to develop the test environment for the `DamageScanner`.

# Acknowledgements
This project has received funding from the European Union’s Horizon Europe research and innovation programme (Grant Agreement No. 101093854). 

# References
