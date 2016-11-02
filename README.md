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

Analysis the raw mass spectra presented in the file "filename" using the functions available in MIDcor:

> correct("filename",samb=1,samf=3,cndb=4,cndf=9)
        
Here "filename" is the name of the input file contaiming raw mass spectra. The parameters samb,samf,cndb,cndf are the positions of the initial and final characters in the row names designating the biological sample and conditions correspondingly
Here is an example of an input file. Comments are included between * *, they should not be present in a real input file.

carbons_total 18            * Total number of carbons in the derivated fragment *

fragment 4                     *<number of carbons in the assayed metabolite (or its fragment) >*

silicium 3                       *<Total number of Si atoms in the derivated fragment>*

sulfur   0                       *<Total number of S atoms in the derivated fragment.>*

m/z                  	417	418	419	420	421	422	423 

commer_01.D            	911	144256	52288	24832	9517	2647	1037

commer_02.D            	929	127872	46416	21664	7421	2285	688

commer_03.D            	1039	135296	49192	22760	7967	2202	681

natural_01.D           	782	109024	39368	18592	7084	1894	674

natural_02.D           	1059	134528	49200	22968	7643	2027	644

natural_03.D           	1177	151808	54456	25152	8434	2396	773

B10n8_01.D             	867	105080	40672	29472	22152	68248	21936

B10_02.D             	1294	148416	56784	41776	29216	95936	31608

B10_03.D             	1311	150528	57464	42104	28776	95184	30776

B11n8_01.D             	673	67576	25768	19064	16002	44320	14340

B11_02.D             	1113	108120	41208	30208	25192	71648	23056

B11_03.D             	1056	98024	37736	27880	22280	65240	21328

B12n8_01.D             	1091	250304	95240	74224	48360	180288	57240

B12_02.D             	918	117448	45728	36096	28784	95920	30728

B12_03.D             	1540	194560	75856	58952	44088	152704	48896

B4_c801.D              	1296	239744	87064	43832	16085	4421	1499

B4_02.D              	1468	283968	103032	52080	16744	4777	1367

B4_03.D              	1598	300224	108576	54832	17880	5120	1375

B5_c801.D              	1277	212160	78016	38984	13210	3783	1105

B5_02.D              	1895	331392	120048	60376	17800	5044	1326

B5_03.D              	1835	297664	108208	55768	17232	4738	1194

B6_c801.D              	989	146240	53128	26744	12221	3627	1254

B6_02.D              	516	23536	8442	4540	5035	1705	809

B6_03.D              	515	2997	1204	719	4041	1435	657

R10n8_01.D             	856	110320	43184	34032	22680	66200	21704

R10_02.D             	1197	163072	65632	50376	31136	98704	31080

R10_03.D             	1421	175552	69144	53384	32440	103688	32432

R11n8_01.D             	790	108040	42136	32840	21728	63776	20312

R11_02.D             	1214	169984	67200	52344	33120	100352	32424

R11_03.D             	1207	161856	63624	49032	30968	96616	30776

R12n8_01.D             	791	111560	44456	33976	21560	64456	20120

R12_02.D             	957	141440	56584	43480	29032	84056	26656

R12_03.D             	1077	154112	60488	47360	31152	92824	29616

R4_c801.D              	242	189	293	169	150	92	87

R4_02.D              	446	255	302	176	162	88	85

First four lines provide information necessary for the correction for natural occurring isotopes. The fifth line indicates the m/z range for the spectra starting from M-1 isotope. The next lines consist of a word representing a name of a given sample and the intensities of peaks corresponding to the indicated m/z values. The last line is recognized by “fin” in the beginning.

First of all, MIDcor corrects the spectra for proton lost due to electronic impact, defined by the ratio (M-1)/M. Then it normalizes the intensities in corrects each row and corrects for the natural occurrence of isotopes. It writes the results in "filename_c" under the title "*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***"

Then MIDcor sums the corresponding intensities presented in multiple injections from the same biological samples (technical replicates). In the example provided it recognizes an individual biological sample by three first characters. The same first three characters indicate that the corresponding injections are from the same sample. Respectively, the parameters samb and samf should have values 1 and 3. In principle, samples can be characterized by characters in different positions, and, respectively the parameters samb and samf can have different values. After summing the technical replicates, it corrects the sums for natural occurring isotopes and writes the result under the title "*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **".

After summation, the technical replicates the set of spectra contains only different biological replicates. Then MIDcor corrects them for peaks overlapping and finds the mean values and standard deviations among the replicates referring to the same conditions of incubation. The program recognizes the same conditions by the same characters that indicate the conditions. In the example provided it is the next two characters. In this case, the parameters cndb and cndf should have values 4 and 5 correspondingly. It shows the results under the titles:
*** Statistics, samples fully corrected **

To correct the peaks overlapping, MIDcor calculates the difference D between the measured mass spectra and the ones theoretically for the unlabeled samples. The unlabeled samples are recognized by the substring "commer" (for commercial preparation in minimal media containing only the substance and chemicals for its derivatization) or "natural" (containing the unlabeled substance in full media of cell incubation before the incubation starts).

Midcor performs the correction for peaks overlapping in labeled samples by one of two ways; either (i) subtract the obtained D without any change from the normalized measured spectra, or (ii) change the obtained D before its subtraction. The change includes shift D by the number of 13C labels obtained from artificially 13C enriched substrates and taking the part of D equal to the portion of corresponding mass isotopomer. To decide which of the two ways should be taken MIDcor compares the two D-vectors obtained for "commercial" or "natural" unlabeled samples. If the max value of D obtained for "natural" media is more than 0.5% larger than the one obtained in minimal media than MIDcor chooses the way (i) with D obtained for the "natural" conditions. In opposite cases, it chooses the way (ii).

Finally, it writes the D obtained for the "natural" conditions in the output file under the title

*** Correction factor: **

If a user decides to check whether some "isotopic effect" takes place, he/she should provide the set of measured intensities for an "artificially" labeled sample with labeling pattern a priori known (it could be a commercial preparation). The name of such a sample should contain the substring "fitf". Midcor fits the known 13C mass isotopomer pattern thus determining the factor of "isotopic effect"


## Publications
- “MIDcor”, an R-program for deciphering mass interferences in mass spectra of metabolites enriched in stable isotopes. Submitted to BMC bioinformatics.


