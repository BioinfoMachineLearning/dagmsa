# DAGMSA
DAGMSA: Directed Acyclic Graph-based Multiple Sequence Alignments (MSA) for Protein Multimer Structure Prediction

This software is designed to generate multiple sequence alignments (MSAs) of multimeric protein complexes in the Modeller .pir and Alphafold2-multimer (AF2_M) formats.


###                  Downloading             

Just type the following command to download the software.


`git clone https://github.com/BioinfoMachineLearning/dagmsa`

You can also use wget to download the software using the following commands:

```
mkdir dagmsa
cd dagmsa
wget https://github.com/BioinfoMachineLearning/dagmsa
```

Once the download is complete, follow the steps in the installation section to install and run the software.


###                  Installation             

This software was developed and tested using python, perl and other dependent software. Please make sure the following softwares versions are installed:

    (1) Python 3.6 or Above

    (2) Modeller 10.0 or above (https://salilab.org/modeller/download_installation.html)

    (3) Perl


The following python dependencies (packages) are needed for the code to execute:

    (1) H5py 2.9.0


Besides the python packages, In order to generate features for the deep learning predictor, the following software tools are necessary:

    (1) HH-suite-3.0 available at: https://github.com/soedinglab/hh-suite

    (2) JackHMMER/HMMER-3.1 available at: http://hmmer.org/download.html

    (3) Latest HH-suite searchable protein database like UniRef30_2020_06: available at: http://gwdu111.gwdg.de/~compbiol/uniclust/2020_06/

    (4) Modeller searchable PDB template database list (https://salilab.org/modeller/downloads/pdball.pir.gz)

### Usage
Use any tool either Modeller or HHsuite to generate the MSA of the individual monomers. Then put them in a single folder and add the file names into an input text file then run the following command:

```
#Usage
sh runDAGMSA_pir.sh <input_monomer_msa_file_directory> <list_file> <output_directory>

or 

sh runDAGMSA_hhb.sh <input_monomer_msa_file_directory> <list_file> <output_directory>

#Example:
sh runDAGMSA_pir.sh <input_monomer_msa_file_directory> <list_file> <output_directory>
 
```

