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
    * See following link for details on installation: https://www.atnf.csiro.au/computing/software/miriad/INSTALL.html
4. STARLINK software package for GAUSSCLUMP extraction
    * Get the lastest Starlink release at http://starlink.eao.hawaii.edu/starlink/Releases


## 1. Prepare data

Assume we have several CO data cubes for the same region. Firstly we need to make cubes comparable to each other.

1. Convert cubes from fits to Miriad dataset using command ``fits in=Data1.fits out=Data1 op=xyin``
2. If nessesary, convolve cubes to the same beam size using command ``convol map=Data1 out=Data1.conv beam=13.2`
3. If some cubes are in galactic coordinates, then convert them to equatorial using command ``regrid in=Data1.conv out=Data1.conv.regrid options=galeqsw,offset``

## 2. Extract clumps

We assume that 13CO(1-0) line will be used for clump extraction. In order to extract clumps, the following command of Starlink should be used:

```
convert
fits2ndf 13CO.fits 13CO
cupid
findclumps in=13CO  out=13CO_model outcat=13CO_cat rms=0.18 method=gaussclumps config=^gaussclumps.cfg WCSPAR=true LOGFILE=catalog.txt
ndf2fits 13CO_model 13CO_model.fits
```
The content of gaussclumps.cfg is following (should be modified according to data):
```
GaussClumps.FwhmBeam=2.0 # beam size of cube in pixels
GaussClumps.ExtraCols=1 # Output nessesary columns
GaussClumps.NPad=100 # Specifies a termination criterion for the GaussClumps algorithm.  The algorithm will terminate when "Npad" consecutive clumps have been fitted all of which have peak values less than the threshold value specified by the "Thresh" parameter, or when one of the other termination criteria is met. [10] 
GaussClumps.Thresh=1.0 # Gives the minimum peak amplitude of clumps to be fitted by theGaussClumps algorithm (see alsoGaussClumps.NPad). The supplied value is multipled by the RMS noise level before being used. [2.0] 
GaussClumps.VeloRes=2.0 ## The velocity resolution of the instrument, in channels. The velocity FWHM of each clump is not allowed to be smaller than this value. Only used for 3D data. [2.0] 
```
As the result, we will have catalog of clumps in file catalog.txt and model of clumps in file 13CO_model.fits. The catalog should be converted to any table processor program (like MS Excel, Google Sheets, etc.) to work with it.
