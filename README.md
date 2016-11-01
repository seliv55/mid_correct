![Logo](text4217.png)

# MIDcor
Version: 1.0
## Short Description

“R”-program that corrects 13C mass isotopomers spectra of metabolites for natural occurring isotopes and peaks overlapping

## Description

MIDcor is an “R”-program that performs a primary analysis of isotopic isomers (isotopomers) distribution obtained by Gas Cromatography coupled with Mass Spectrometry (GCMS). The aim of this analysis is to have a correct distribution of isotopes originated from substrates that are artificially enriched with specific isotopes (usually 13C). To this end the program performs a correction for natural occurring isotopes and also correction for “impurities” of the assay media that give peaks overlapping with the spectra of analyzed labeled metabolites. This program offers two ways of corrections of “impurities” resulted from overlapping the assayed mass isotopomer distribution with peaks produced either by unknown metabolites in the media, or by different fragments produced by the assayed metabolite. 

This program offers two ways of corrections of “impurities” resulted from overlapping the assayed mass isotopomer distribution with peaks produced either by unknown metabolites in the media, or by different fragments produced by the assayed metabolite.

## Key features

- primary processing of 13C mass isotopomer data obtained with GCMS

## Functionality

- Preprocessing
- Statistical Analysis
- Workflows

## Approaches

- Isotopic Labelling Analysis / 13C
    
## Instrument Data Types

- MS

## Data Analysis

- correction for H+ loss produced by electron impact, natural occurring isotopes, peaks overlapping that depends on isotopic composition of metabolites

## Tool Authors

- Vitaly Selivanov (Universitat de Barcelona)

## Container Contributors

- [Pablo Moreno](EBI)

## Website

- N/A

## Git Repository

- https://github.com/seliv55/mid_correct

## Installation

-  As independent program. MIDcor itself does not require installation. Standing in the MIDcor directory enter in R environment with the command:
  
'''  R '''
  
 read the necessary functions:
  
''' source("lib.R")
  
source("midcor.R")'''
  
  
## Usage Instructions

1. Read a particular input file with data designed for analysis and analyze them using the functions available in Midcor:
> correct("filename")
        
Here "filename" is the name of the input file contaiming raw GCMS data
The input data is a file with raw mas/charge (m/z) distribution provided directly by a GCMS machine. An example of input data file is shown below. Comments are included between * *, they should not be present in a real input file. The files "GluC2C4b" and "GluC2C5" provide the examples. 

*below is an example of an input file *

carbons_total 18            * Total number of carbons in the derivated fragment *

fragment 4                     *<number of carbons in the assayed metabolite (or its fragment) >*

silicium 3                       *<Total number of Si atoms in the derivated fragment>*

sulfur   0                       *<Total number of S atoms in the derivated fragment.>*

m/z                  	417	418	419	420	421	422	423 

R29_01.D             	911	144256	52288	24832	9517	2647	1037 

R29_02.D             	929	127872	46416	21664	7421	2285	688 

R29_03.D             	1039	135296	49192	22760	7967	2202	681 

B29_01.D             	782	109024	39368	18592	7084	1894	674 

B29_02.D             	1059	134528	49200	22968	7643	2027	644 

B29_03.D             	1177	151808	54456	25152	8434	2396	773 

B10n8_01.D             	867	105080	40672	29472	22152	68248	21936 

B10n8_02.D             	1294	148416	56784	41776	29216	95936	31608 

B10n8_03.D             	1311	150528	57464	42104	28776	95184	30776 


* whatever number of lines, where the first word is the name of the sample. First three characters of this name, if they are the same, correspond to various injections of the same sample (). Next two characterizes should be the same for various samples corresponding to the same conditions of incubation. The subsequent values correspond to the given m/z in the GCMS spectrum.
Below is the last line, recognized by “fin” in the beginning.   *

finR24_02.D            	1381	101072	41633	33811	22779	53581	17014

* end of an input file *


The function “correct” reads the user provided GCMS data file, normalizes the input GCMS data, so that the sum of m/z values equals to 1, corrects the obtained distribution for the presence of naturally occurring isotopes, and saves the obtained corrected distributions, calculated for each line separately, in the output file, which name is the name of input file with “_c” added at the end. This part of analysis is saved under the title:
 "*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***"
 
 Then it sums the corresponding intensities in the injections referred to the same biological sample, performs the above described operations with respect to such sums, and shows the results under the title:
 *** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **

Then it calculates the difference D between the corrected distribution and theoretically expected one for the reference unlabeled sample, and  uses the obtained difference to ultimately correct for impurity the samples containing metabolites of living cells in the presence of labeled substrates, finds the means and standard deviation values betweem biological samples referred to the same conditions, and shows the results under the titles:
*** Statistics, samples fully corrected **
and
*** Correction factor: **

Examples of input data can be found in the files “GluC2C4b”, “GluC2C5” and the corresponding output data are in the files “GluC2C4b_c”, “GluC2C5_c”.



## Publications
- “MIDcor”, an R-program for deciphering mass interferences in mass spectra of metabolites enriched in stable isotopes. Submitted to BMC bioinformatics.


