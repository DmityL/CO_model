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
2. In nessesary, add rest frequency information using command ``puthd in=Data1/restfreq value=
3. If nessesary, convolve cubes to the same beam size using command ``convol map=Data1 out=Data1.conv beam=13.2`
4. If some cubes are in galactic coordinates, then convert them to equatorial using command ``regrid in=Data1.conv out=Data1.conv.regrid options=galeqsw,offset``

Here is the table containing rest frequency F0 (in Hz) and value of T0 (= h*nu/c) of different CO lines
|Transition|F0(CO)|F0(13CO)|F0(C17O)|F0(C18O)|T0(CO)|T0(13CO)|
|-----|-----|-----|-----|-----|-----|-----|
|1-0|115.2712018E+09|110.2013543E+09|112.3592837E+09|109.7821734E+09|5.5|5.3|
|2-1|230.5380000E+09|220.3986841E+09|224.7143850E+09|219.5603541E+09|11.1|10.6|
|3-2|345.7959899E+09|330.5879652E+09|337.0611298E+09|329.3305525E+09|16.6|15.9|
|4-3|461.0407682E+09|440.7651735E+09|449.3953412E+09|439.0887658E+09|22.1|21.2|
|5-4|576.2679305E+09|550.9262850E+09|561.7127845E+09|548.8310055E+09|27.7|26.4|
|6-5|691.4730763E+09|661.0672766E+09|674.0093443E+09|658.5532782E+09|33.2|31.7|
|7-6|806.6518060E+09|771.1841255E+09|786.2808166E+09|768.2515933E+09|38.7|37.0|
|8-7|921.7997000E+09|881.2728093E+09|898.5230217E+09|877.9219553E+09|44.2|42.3|

The convolution size should be computed as Beam = sqrt(Beam_final^2 - Beam_source^2)


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

The Gaussclump catalog is looks like this:
```
|Index|Peak1|Peak2|Peak3|Cen1|Cen2|Cen3|Size1|Size2|Size3|Sum|Peak|Volume|GCMEANPEAK|GCFWHM1|GCFWHM2|GCFWHM3|GCVELGRAD1|GCVELGRAD2|GCANGLE|GCBG|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|1|173.6809749|2.862433769|-18947.63233|173.6801833|2.863415146|-18915.47601|41.78595814|68.5669739|689.714063|9436.09794|18.24441286|3.59E+08|11.18156155|4.956228399|7.996271926|12.6606003|-0.389197949|-0.325293626|-11.25550824|3.186878383|
|2|173.6309618|2.887428652|-20807.14745|173.6278437|2.885465973|-20849.01349|75.25239955|56.41409135|610.708283|10992.94493|16.82952066|4.29E+08|21.76989453|9.325969233|5.606977234|11.49370661|-0.071603437|-0.349490468|-28.16718025|2.608133951|
|3|173.7184992|2.693687483|-16689.64967|173.7184102|2.695273252|-16629.26003|67.78094557|80.11623379|627.5875703|12811.11741|14.49698834|5.51E+08|30.75038375|6.91981433|9.937106321|11.70617921|-0.389524152|0.077169482|-31.76902994|4.204375308|
|4|173.6184711|2.774927344|-21338.43749|173.6168261|2.773805867|-21355.42696|61.29618221|54.99242473|701.1990071|9513.704235|15.55134412|3.98E+08|40.37107352|7.904827149|5.574627009|13.10262736|-0.333176555|0.162672734|-36.46933829|1.274985791|
|5|173.6372203|2.818679302|-19744.56738|173.6347192|2.817034972|-19691.52315|57.1025503|68.07299746|706.3527784|10171.1941|13.79874127|4.58E+08|48.89812688|6.643001077|8.007835062|12.83113146|0.62717653|0.218486793|15.49746208|0.901523622|
|6|173.1684387|2.356122299|-19611.74487|173.1654713|2.358375097|-19506.44966|94.50408122|71.09663136|976.7823119|19568.50905|11.15808744|9.94E+08|55.82002788|7.922738548|11.16532352|17.69790103|0.092286987|0.726093703|108.4669192|5.609638543|
|7|173.7810138|2.668693433|-16424.00466|173.7843963|2.666735238|-16408.70335|60.94708659|41.1504897|604.2729987|5338.051414|12.06671065|2.73E+08|70.68626515|7.256585774|5.060505622|11.73578738|0.039881272|-0.111926251|8.674529555|1.724292739|
```
We need to add extra columns: coordinates RA/Dec in HMS/DMS format. I prefer to use TOPCAT software for coodinate conversion: http://www.star.bris.ac.uk/~mbt/topcat/

## 3. Extract spectra

For each clump peak position that were found using GAUSSCLUMP we have to extract the spectra from each data cube. Assume that we create subfolders for each data cube: `_FCRAO_12CO_spectra`, `_FCRAO_13CO_spectra`, `_SMT_12CO_spectra`, `_SMT_13CO_spectra`. 

Then spectra for each clump should be extracted using following commands (coordinates and clump number should be replaced using GAUSSCLUMP catalog):
```mbspect in=12CO_data coord=05:41:28.74,+35:48:56.5 log=_FCRAO_12CO_spectra/Clump001.dat device=_FCRAO_13CO_spectra/Clump001.ps/PS options=measure xaxis=velo > _FCRAO_13CO_spectra/Clump001.fit```
We will now have extracted spectra and linewidth estimation for each line.

After we extract the spectra for each clump, we should combine different lines for each clump to single file. That can be done using *spect_comb.pl* script from this repo. The source code should be modified to fit your data. The script will automatically compute RMS for each spectra using first 70 channels of each file. The number of channels to cumpute RMS should be also modified in the file *spect_comb.pl*. We create folder *_spect_comb* to store the combined spectra.
```perl spect_comb.pl```
The resulting files ClumpNNN.dat looks like this (first column - velocity, second - intensity, last column - RMS estimation):
```
|-34.98225403|1.69511E+00|1.39041037401241|
|-34.85528946|-3.24238E-02|1.39041037401241|
|-34.72832489|-2.70249E+00|1.39041037401241|
|-34.60136032|1.10365E+00|1.39041037401241|
|-34.47439575|-1.53760E+00|1.39041037401241|
```
Note that spectra are beind combined using velocity shift specified in file *spect_comb.pl* ($dv = 35). If velocity shift will be too small, then data will be broken. 
The important parameter is velocity inverval for emission in line 67 of *spect_comb.pl* ($v > -35 and $v < 0). All data points outside of this interval will not be included to the combined spectra. Thus if you have emission line at ~ -20 km/s, then selecting velocity inverval -35>v>0 and velocity shift dv = 35 is good way to go.

## 4. Create an initial estimate of model parameters

In order to get an initial estimate for each clump we need following values for each clump: peak values of 12CO and 13CO lines, linewidth of 13CO line. The values are being extracted from the spectra using `scan_spectra.pl` and `scan_fit.pl` utils.  

Before we start extracting we need to create the simplified catalog of clumps that we name `clumps_cat.csv` with following content:
```
Clump	Peak1	Peak2	V
Clump001	173.681	2.862	-18.95
Clump002	173.631	2.887	-20.81
Clump003	173.718	2.694	-16.69
```

The `scan_spectra.pl` tool looks for ClumpNNN.dat files in the specific folder (should be specified in the source code) and using the peak velocities of each clumps from the Clumps_l_b_v.csv file extract the value of spectra intensity at the specific velocity. Actually it extracts three nearest points at specific velocity and returns the average of three points. Using this tools we extract intensity of each clump at peak velocity in several lines: 13CO and 12CO. The resulting files for each CO line looks like this:
```
Clump	Tpeak
Clump001	28.25
Clump002	27.92
Clump003	14.55
Clump004	27.87
Clump005	34.28
Clump006	13.53
Clump007	10.45
```
The `scan_fit.pl` tool scan for fit files *.fit that comes from the previous spectra extraction step. It looks for the following line: ``#MNW Line width:`` and push the found value to the CSV file for each clump, thus extracting the estimation of the linewidth. The resulting file looks like this:
```
Clump	FWHM
Clump001	1.927
Clump002	1.855
Clump003	2.412
Clump004	2.278
Clump005	1.765
Clump006	3.46
Clump007	1.644
```
After executing these tools we will have following files: _FCRAO_12CO_spectra_Tpeak.csv, _FCRAO_13CO_spectra_Tpeak.csv, _FCRAO_13CO_spectra_FWHM.csv. These files should be combined in any table processor program (easy to do with copy and paste). The resulting file
