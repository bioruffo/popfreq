# popfreq
Show population frequency tables from Ensembl data

# What is popfreq #

This script builds haplotype frequency tables from population data obtained from Ensembl. It is useful to determine ancestral and derived alleles.


1 . Requisites:  

 * Python 3.5+ (I suggest [Anaconda Python](https://www.continuum.io/downloads))  
 * The Pandas library

2. Downloading data from Ensembl:

 * Search Ensembl for a variant e.g. [rs699](https://www.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=1:230709548-230710548;v=rs699;vdb=variation;vf=19167044)
 * Click on "Population genetics"
 * On the 1000 genomes table, population (line) "ALL", column Genotypes, click "Show"
 * On the full genotypes table, hover on the spreadsheet icon on the top right and click "download whole table".

# Using popfreq #

As first commit, this script is edited in place to alter rs parameters. This will be improved on later.

# LICENSE #

popfreq is offered under the MIT License.
Please read it here: https://opensource.org/licenses/MIT
