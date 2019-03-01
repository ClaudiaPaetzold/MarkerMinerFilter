# MarkerMinerFilter
filter [MarkerMiner](https://bitbucket.org/srikarchamala/markerminer.git) (Chamala et al., 2015) output for exon length and variability to create fasta files with probes

### Prerequisites

This script runs exclusively on Python 2.x 

In addition, the script requires BioPython to be installed. 

BioPython can be installed 

via pip
```
pip install biopython
```

or with conda
```
conda install -c anaconda biopython
```
The output file of the MarkerMiner software has to be extracted before using the MinerOutpurFilter.

### Installing

download the script from this site. In a console, navigate to the folder you saved it in (or make it executable, if you prefer) and run it like so:
```
python MinerOutputFilter.py [Parameters] [Arguments]
```

## Usage

when employing the MinerOutputFilter you have will can filter the output of the 
```
usage: MinerOutputFilter.py [-h] --reference
                            {Tcacao,Rcommunis,Bdistachyon,Mesculenta,Osativa,Gmax,Sbicolor,Cpapaya,Ptrichocarpa,Vvinifera,Zmays,Fvesca,Mtrunculata,Mdomestica,Alyrata,Athaliana}
                            [--minlength MINLENGTH] [--maxlength MAXLENGTH]
                            [--maxvar MAXVAR] [--minvar MINVAR] [--cleanup]
                            MarkerMinerOutputDir

filter MarkerMiner (Chamala et al., 2015) output for exon length and variability to create fasta files with probes

positional arguments:
  MarkerMinerOutputDir  full path to extracted MarkerMiner output directory
optional arguments:
  -h, --help            show this help message and exit
  --reference {Tcacao,Rcommunis,Bdistachyon,Mesculenta,Osativa,Gmax,Sbicolor,Cpapaya,Ptrichocarpa,Vvinifera,Zmays,Fvesca,Mtrunculata,Mdomestica,Alyrata,Athaliana}
                        reference used in MarkerMiner analysis.
  --minlength MINLENGTH
                        minimum required exon length. Default: 120
  --maxlength MAXLENGTH
                        maximum allowed exon length. Default: 1000
  --maxvar MAXVAR       maximum allowed exon variability (discarding
                        reference). Default: 10
  --minvar MINVAR       minimum required exon variability (discarding
                        reference). Default: 0
  --cleanup             if selected, all intermediates files will be deleted
```

