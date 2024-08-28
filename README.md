# Flatworm CpG methylation exploration

## Aim:
There is no real consensus of CpG methylation presence in flatworm genomes <sup>[1](#myfootnote1),[2](#myfootnote2),[3](#myfootnote3),[4](#myfootnote4)</sup>.

Using PacBio HiFi read kinetics I attempt to find evidence for CpG methylation in flatworms.

## Repo structure:

I divided the project into two subparts. 

1. Indirect evidence for mCpG presence could be the presence of CpG islands. Hence a CpG island analysis of flatworm and other genomes is found in the `cpg_island/` directory.
2. Using kinetic information extracted during the sequencing and deep learning models I called CpG modifications in the _S.mediterranea_ genome. The analyses for this sub-project can be found in the `cytosine_methylation/` directory.

## Some VERY SUPERFICIAL conclusions:

* CpG methylation might be present at very low levels (see: `cytosine_methylation/outputs/cpg_methlyation-perc_genome-wide_histogram.pdf`). 

* The methylation percentage around gene bodies is low, especially in the promoter regions (see: `cytosine_methylation/outputs/ssmed_5mC_profile-mean.pdf`).

* As the PacBio reads originate from whole adults information is lost from different cell lineages. Might be some more interesting things present in specific cell lines. Or not.


## Bibliography:

<a name="myfootnote1">1</a>: Geyer, Kathrin K., Iain W. Chalmers, Neil MacKintosh, Julie E. Hirst, Rory Geoghegan, Mathieu Badets, Peter M. Brophy, Klaus Brehm, and Karl F. Hoffmann. 2013. “__Cytosine Methylation Is a Conserved Epigenetic Feature Found throughout the Phylum Platyhelminthes.__” BMC Genomics 14(1):462. doi: 10.1186/1471-2164-14-462.

<a name="myfootnote2">2</a>: Geyer, Kathrin K., Carlos M. Rodríguez López, Iain W. Chalmers, Sabrina E. Munshi, Martha Truscott, James Heald, Mike J. Wilkinson, and Karl F. Hoffmann. 2011. “__Cytosine Methylation Regulates Oviposition in the Pathogenic Blood Fluke Schistosoma Mansoni.__” Nature Communications 2(1):424. doi: 10.1038/ncomms1433.

<a name="myfootnote3">3</a>: Jaber-Hijazi, Farah, Priscilla J. K. P. Lo, Yuliana Mihaylova, Jeremy M. Foster, Jack S. Benner, Belen Tejada Romero, Chen Chen, Sunir Malla, Jordi Solana, Alexey Ruzov, and A. Aziz Aboobaker. 2013. “__Planarian MBD2/3 Is Required for Adult Stem Cell Pluripotency Independently of DNA Methylation.__” Developmental Biology 384(1):141–53. doi: ht<span>tp://</span>/10.1016/j.ydbio.2013.09.020.

<a name="myfootnote4">4</a>: Raddatz, Günter, Paloma M. Guzzardo, Nelly Olova, Marcelo Rosado Fantappié, Markus Rampp, Matthias Schaefer, Wolf Reik, Gregory J. Hannon, and Frank Lyko. 2013. “__Dnmt2-Dependent Methylomes Lack Defined DNA Methylation Patterns.__” Proceedings of the National Academy of Sciences 110(21):8627–31. doi: 10.1073/pnas.1306723110.
