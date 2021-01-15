# CO LTE multi-layer radiative transfer model
### Astrophysical method for describing and modelling the complex line profiles of CO molecules

In this work, we investigated the features of the analysis of complex CO line profiles in giant molecular clouds (GMCs). A technique has been developed to use several emission lines of the CO molecule to build a model and determine GMCs physical parameters within the local thermodynamic equilibrium framework. The technique includes clumps extraction using the GAUSSCLUMP algorithm and constructing a multilayer radiation transfer model for clumps using optimization and Monte Carlo methods. As an example, the technique was applied to analyze the large-scale mapping of the S231-S235 star formation complex in four different CO lines.

## 1. Install nessesary instruments

In order to go throught the analysis we need several instruments.

1. Perl (installed by default in must Linux-based OS) with the following modules:
    - List::MoreUtils (install using `cpan install List::MoreUtils`)
    - List::Util (install using `cpan install List::Util`)
    - Number::Closest (install using `cpan install Number::Closest`)
    - PDL::LiteF (install using `cpan install PDL::LiteF`)
    - PDL::NiceSlice (install using `cpan install PDL::NiceSlice`)
    - PDL::Stats (install using `cpan install PDL::Stats`)
2. Python (version 3.7 or later) with the following modules:
    - Numpy (install using command `pip install numpy`)
    - scipy.optimize (install using command `pip install scipy`)
    - scipy.signal 
    - emcee (install using command `pip install emcee`)
    - corner (install using command `pip install corner`)
3. ATNF Miriad package for spectra extraction
    See following link for details on installation: https://www.atnf.csiro.au/computing/software/miriad/INSTALL.html
4. STARLINK software package for GAUSSCLUMP extraction
    Get the lastest Starlink release at http://starlink.eao.hawaii.edu/starlink/Releases


## 1. Prepare data

Assume we have several CO data cubes for the same region
