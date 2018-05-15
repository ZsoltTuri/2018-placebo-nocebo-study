# Evidence for Cognitive Placebo and Nocebo Effects in Healthy Individuals

<!---
[![DOI](https://zenodo.org/badge/19634/ihrke/2016-placebo-tdcs-study.svg)](https://zenodo.org/badge/latestdoi/19634/ihrke/2016-placebo-tdcs-study)

This repository contains data and analyses for the paper "Evidence for Cognitive Placebo and Nocebo Effects
in Healthy Individuals".

 If you want to use this data/analysis in a research publication,
please cite [our paper](http://link).


Turi, Z., Bjørkedal, E., Gunkel, L., Antal, A., Paulus, W. & Mittner, M. (yyyy).
Evidence for Cognitive Placebo and Nocebo Effects in Healthy Individuals. Journal. 

~~~{bibtex}
@article{turi_placebo_yyyy,
  author={Turi, Z., Bjørkedal, E., Gunkel, L., Antal, A., Paulus, W. and Mittner, M.},
  title={Evidence for Cognitive Placebo and Nocebo Effects in Healthy Individuals},
  year=yyyy,
  journal={Journal},
  volume=v,
  number=n,
  doi=d
}
~~~
 -->

## Requirements

Analysis are coded in [R](http://r-project.org) and [stan](http://mc-stan.org). Quite a lot R-packages and the [Stan](http://mc-stan.org) are required. It is easiest to set up the
R-packages using [conda](https://www.continuum.io/downloads).  We
provide an `environment.yml` file which allows to set up R with all
needed packages with very few commands.

1. download `anaconda` or `miniconda`

 e.g. on linux:
 ~~~{bash}
 wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
 bash Miniconda-latest-Linux-x86_64.sh
 conda update conda
 conda install anaconda-client anaconda-build conda-build
 ~~~

2. clone this repository
 ~~~{bash}
 git clone
 https://github.com/ZsoltTuri/2018-placebo-nocebo-study
 cd 2018-placebo-nocebo-study
 ~~~

3. re-create the environment used for creating these analyses:
 ~~~{bash}
 conda env create
 ~~~

If you are not willing to do this, you will have to look at the
[environment.yml](./environment.yml) file to check all
dependencies.

## Setup

This repository uses the
[ProjectTemplate](http://projecttemplate.net/) directory layout. It
also provides an `environment.yml` which allows to set up R with all
needed packages with very few commands.

## Data

Raw data is located in `data/raw` and is provided in `.csv` format.

The `.R` scripts located in `data` load the raw files into `R`
workspace under the name of the `R`-file (without the `.R` extension).

<!---**NOTE**: there are also pre-processed exports of all the variables discussed next; those are located under [data/export](data/export). These files have been created by the script [src/export_data.R](src/export_data.R).  -->

The data is structured as follows <!--- (refer to [the paper](http://www.link) for
details). -->

<!---
### Demographic
stored in variable `demographic`
~~~
> summary(demographic)
  participant             group         age                        subj   
 1      : 5   n_cond         :16   Min.   :18.00   natural_history_1 : 1  
 2      : 5   p_cond         :16   1st Qu.:22.00   natural_history_10: 1  
 3      : 5   n_cntrl        :16   Median :24.50   natural_history_11: 1  
 4      : 5   p_cntrl        :16   Mean   :24.54   natural_history_12: 1  
 5      : 5   natural_history:16   3rd Qu.:27.00   natural_history_13: 1  
 6      : 5                        Max.   :38.00   natural_history_14: 1  
 (Other):50   
~~~
Variables are coded as follows:

- `participant` - participant identification number from 1 till 16 in each group
- `group`  - group identification; 
	- n_cond: nocebo group
	- p_cond: placebo group
	- n_cntrl: nocebo control group
	- p_cntrl: placebo control group
- `age` - age of the participants in the time of the participation
- `subj` - 64 unique, group-independent participant identification number  -->

### Subjectively reported expected and perceived cognitive performance

stored in variable `subj.outcomes`

~~~
> summary(subj.outcomes)

  participant     group      expected          perceived             effect             type            subj   
 1      : 4   n_cond :16   Length:64          Length:64          nocebo :32   conditioning:32   n_cntrl_1 : 1  
 2      : 4   p_cond :16   Class :character   Class :character   placebo:32   control     :32   n_cntrl_10: 1  
 3      : 4   n_cntrl:16   Mode  :character   Mode  :character                                  n_cntrl_11: 1  
 4      : 4   p_cntrl:16                                                                        n_cntrl_12: 1  
 5      : 4                                                                                     n_cntrl_13: 1  
 6      : 4                                                                                     n_cntrl_14: 1  
 (Other):40  
~~~

Variables are coded as follows:

- `participant` - participant identification number from 1 till 16 in each group
- `group`  - group identification; 
	- n_cond: nocebo group
	- p_cond: placebo group
	- n_cntrl: nocebo control group
	- p_cntrl: placebo control group
- `expected` - expected direction of change in cognitive performance
    - decline means that subjects expected to get worse, 
    - neutral means no expected change, 
    - improve means subjects expected to have better performance
- `perceived` - perceived direction of change in cognitive performance
	- decline means that subjects perceived to be worse, 
    - neutral means no perceived change, 
    - improve means subjects perceived to be better
- `effect` - whether the group belongs to nocebo or placebo effect
    - nocebo: n_cond and n_cntrl
    - placebo: p_cond and p_cntrl
- `type` - type of manipulation
	- conditioning: n_cond and p_cond
	- control: n_cntrl and p_cntrl
- `subj` - 64 unique, group-independent participant identification number


### Data from the reward-based learning task (learning phase)

data from the five different groups are stored in `d`

~~~
> summary(d)
  participant        group      symbl_pair symbl_position    accuracy      reaction_time        trial            reward      
 1      : 2400   n_cntrl:7680   AB:12800   Min.   :1.0    Min.   :0.0000   Min.   :0.0333   Min.   :  1.00   Min.   :0.0000  
 2      : 2400   n_cond :7680   CD:12800   1st Qu.:1.0    1st Qu.:1.0000   1st Qu.:0.6830   1st Qu.: 60.75   1st Qu.:0.0000  
 3      : 2400   nhg    :7680   EF:12800   Median :1.5    Median :1.0000   Median :0.8331   Median :120.50   Median :1.0000  
 4      : 2400   p_cntrl:7680              Mean   :1.5    Mean   :0.7861   Mean   :0.8665   Mean   :120.50   Mean   :0.6165  
 5      : 2400   p_cond :7680              3rd Qu.:2.0    3rd Qu.:1.0000   3rd Qu.:1.0330   3rd Qu.:180.25   3rd Qu.:1.0000  
 6      : 2400                             Max.   :2.0    Max.   :1.0000   Max.   :1.6662   Max.   :240.00   Max.   :1.0000  
 (Other):24000                                                             NA's   :343                                       
 day               subj      
 1:19200   n_cntrl_1 :  480  
 2:19200   n_cntrl_10:  480  
           n_cntrl_11:  480  
           n_cntrl_12:  480  
           n_cntrl_13:  480  
           n_cntrl_14:  480  
           (Other)   :35520 
~~~

Variables are coded as follows:

- `participant` - participant identification number between 1 and 16 in each group
- `group`  - group identification; 
	- n_cond: nocebo group
	- p_cond: placebo group
	- n_cntrl: nocebo control group
	- p_cntrl: placebo control group
- `symbol_pair`  - pair type (AB,CD,EF) with (80/20, 70/30 or 60/40 % reward contingency)
- `accuracy` - accuracy: 1 correct, 0 incorrect
- `reaction_time` - reaction time in s
- `trial` - trial number between 1 and 240
- `reward` - 1: reward was received, 0: no reward
- `day` - 1: first day, 2: second day
- `subj` - 64 unique, group-independent participant identification number
 

## Analyses

All analyses are located in `src/`. To run the scripts, you need to
have the `ProjecTemplate` package and various other packages
installed.

The first two lines in each file
~~~{R}
library(ProjectTemplate)
load.project()
~~~
convert the raw data into a more convenient format by

1. running the `data/<dataset>.R` file
2. running the preprocessing scripts in `munge`
3. loading the convenience functions in `lib`
