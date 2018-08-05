# VEGFC
Python3 script for a phylogenetic analysis of VEGF-C amino acid sequences using BioPython [https://biopython.org], t_coffee [http://www.tcoffee.org/Projects/tcoffee], PhyML [http://www.atgc-montpellier.fr/phyml] and the ETE toolkit [http://etetoolkit.org]. Its main goal is not to get the phylogenetic relationship of the species correct (since that has been achieved multiple times elsewhere by using larger datasets), but to identify whether there are distinct subgroups within the VEGF-C orthologs and whether and how they cluster to the different vertebrate clades. It produces a graphical output in SVG format, which can be viewed and edited using e.g. Inkscape [https://inkscape.org/en] or viewed with any SVG-capable bowser like Firefox or Chrome. The animal SVG silouette images are from [http://phylopic.org] (if public domain or CC0) or own creations.


## Requirements
This script has been developed to run under Linux and was specifically tested ONLY under Ubuntu 16.04 and 18.04. It has the following requirements:

### UBUNTU PACKAGES
It requires the installation of the following Ubuntu packages:

* python3-setuptools
* clustalw
* mafft
* dialign-tx
* poa
* probcos
* muscle
* kalign
* amap-align
* proda
* prank
* t-coffee
* phyml

>sudo apt -y --show-progress install inkscape python3-setuptools python3-pyqt4 python3-pyqt4.qtopengl python3-pip autoconf t-coffee clustalw mafft dialign-tx poa probcons muscle kalign amap-align proda prank phyml t-coffee imagemagick build-essential libblas-dev liblapack-dev zlib1g-dev libcairo2-dev libcurl4-openssl-dev python3-numpy python3-lxml python3-six

For some reason the Ubuntu package names some of the alignment executables differently and some links need to be created in order for t-coffee to find them:

>sudo ln -s /usr/bin/dialign-tx /bin/dialign-t
>sudo ln -s /usr/bin/clustalw /bin/clustalw2

### MANUAL INSTALLS

In addition, the pcma executable need to be manually downloaded from [http://prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz] and compiled because it is not available from the Ubuntu repository:

>tar -xvzf pcma.tar.gz
>cd pcma
>make
>sudo cp pcma /bin

### OPTIONAL INSTALLS

The multicore-enabled version of phyml (phyml-mpi) is not available as a precompiled Ubuntu package and needs to be installed manually, but they single-core version works as well (is just slower). The command to execute the multicore version is:
>mpirun -n 4 phyml-mpi -i " + PHYLIP_ALIGNED_TRIMMED_CODED +  " -d aa -b -1

If this script is run on a (headless) server, the xvfb package is required since the ete3 package requires the presence of x.org.

### PYTHON MODULES

The following Python modules need to be installed:

* biopython
* ete3

>sudo pip3 install biopython ete3
