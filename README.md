# **MaCRe: Macro-Complex Reconstructor**

MaCRe is the Macro-Complex Reconstructor program for the SBI/PYT final project. The integers of this group who developed the program are Ricardo Fong Zazueta, Carles Navarro Ramírez and Ignacio Navas Camacho.

In this paper you will find an extend description of what MaCRe does and how does it works.


## **Index**

- [MaCRe](#macre)
- [Theoretical background](#theoretical-background)
- [MaCRe installation and requisites](#macre-installation-and-requisites)
- [MaCRe as package](#macre-as-package)
- [MaCRe tutorial](#macre-tutorial)
- [MaCRe algorithm description](#macre-algorithm-description)
- [Examples](#examples)
- [References](#references) 


## **MaCRe**

MacRe or Macro-Complex Reconstructor is a bioinformatic program which aims to reconstruct a biological macro-complex given a set of predefined paired interactions pdb files of the different chains of the complex. This program is written in Python following a superimposition strategy to build the macro-complexes, using the module Bio.PDB from Biopython as the main tool for performing the reconstruction. 


## **Theoretical background**

Proteins are essential macromolecules for life, they have structural and methabolic funtions for which they need to interact with other proteins or bio-molecules (Garcia-Garcia, J. et al. 2012) like DNA, RNA or other proteins for instance. Therefore, knowing the sturcuture of these biological macro-comlexes is essential for the understanding on how cells work (Lensink M.F., et al, 2016), and how proteins carry out their functions. This knowledge is of  paramount importance to further understand the aetiology of health and disease and for the development of medical aplications such as drug discovery. The structure of these macro-complexes is resolved by experimental procedures such as NMR spectrocospy or X-ray christalography. This information is deposited in public repositories such as the RCSB Protein Data Bank in the form of pdb files containing the information of the tridimentional configuration of structures, chains, residues and atoms of such complexes. Given this latter fact it is important to develop computational methods to carry out complex reconstruction from pdb files in an efficient manner, as the one here developed. Furthermore, based on information of protein primary structure, and its tridimentional conformation it is possible to explore the quality of the computational reconstructions by analizing the energetic architecture of the folds of the protein through programs like prosa2003 that will be here presented as a method to assess the reconstructions of the MaCRe program.   

Protein-protein (PP) interaction prediction is a key field to adress protein function and drug dicovery (Tuncbag, N., et al. 2011), given that this relationship is at the core of the understanding of polipeptide performance. To this end many computational strategies have been developed, such as template based modeling, PP docking and hybrid integrative modeling. Template based modeling is based on the principle of the 3D structure is conserved through evolution (Fiser, A., et al, 2010). PP docking addresses this problem, through mathematical methods like fast fourier transform based on the 3D coordinates of two christallized proteins that interact and through Fourier correlation builds a model (Comeau, S. Rr, et al, 2004). Docking analytical tools have been divided in rigid-body and energy minimization docking methods. These approaches are capable of determining if a given set of polipeptides are interacting, based on a mathemtical modeling approach of the surface of the proteins involved in the analyses. Hybrids methos use computational techniques and advanced experimental methos such as X-ray free-electron laser single- particle imaging, and electron cryo-tomography (Srivastava, A., et al, 2020).


## **MaCRe installation and requisites**

In order to use MaCRe, Python3 must be installed in your computer:
```
$ sudo apt install -y python3-pip
```
Moreover, the following Python modules has to be installed:
- BioPython
```
$ pip3 install biopython
```
- Argparse
```
$ pip3 install argparse
```
For further analysis using MaCRe other external programs need to be installed:
- Prosa
- Chimera


## **MaCRe as package**

The functions and arguments developed for this proyect have been organized in modules. This modules, named as the files arguments.py and functions.py are inside the folder named modules. In this folder the constructor file __init__.py is set up, and the main folder contains the file setup.py. This two files make the program able to be installed in the local machine as a package distribution. 
In order to install the package in the local machine the following command must be run from MaCRe's main folder:
```
$ python3 setup.py install
```
After running this command all the functions of the program can be imported through the synthax: 
```
from modules import [name of the function]
```
or importing the whole module as:
```
impot modules
```

## **MaCRe tutorial**

MaCRe is a recursive bioinformatic program that uses an iterative code to reconstruct biological macro-complexes. In order to use MaCRe from the terminal if it has not been installed, the code must follow this way, making sure that you are in the folder where ComplexRecostruction.py is located:
```
python3 ComplexReconstruction.py -o [output_name] [other_arguments]
```
The arguments of MaCRe can be divided in mandatory and optional:

#### **Mandatory arguments**

- **-o** or **- -out**: store the name of the output file that will be created as final result of the macro-complex reconstruction. This option is used as follows: ```-o example```; this code will create an output file called example.pdb.

#### **Optional arguments**

- **-i** or **- -input**: stores the path where the .pdb files that are going to be analized are located. By default, if no path is defined, the program will use the current working directory and execute with the .pdb files stored in it. To use this argument the code is ```-i path/to/the/files```.

- **-v** or **- -verbose**: by default this option is desactivated. If the value is set to 'True' it activates the possibilty of see the procedures that are being taken places as the program runs. To activate this option the following will be added ```-v```.

- **-f** or **- -force**: by defaul the value of this argument is 'False'. This way, if the output directory arleady exist the program stops and warn this. If the value is set to 'True' even if the output directory exist the program will overwrite it. 

- **-c** or **- -core**: when this option is activated, the program will use the chain selected from the input directory as the core chain to reconstruc the macro-complex. By deafault this option is not activated, and the chain used as core chain is the one with more interactions. To use this the code will be ```-c core_chain.pdb``` where core_chain.pdb is the name of the pdb file which will be used as core molecule.

- -**-s** or **- -stechiometry**: when this argument is not called the program will add as many chains to the complex as chains provided in the interactions. However, if stechiometry is set there are two options: if it is a number, this will be the number of shains added to the complex; if it is a file.txt (.txt format is mandatory) it will indicate how many chains each protein contributes. The file must be written like UNIPROT_id: n_chains. An example to activate this argument is `-s stechiometry.txt`.

- **-a** or **- -analysis**: if this option is activated, a prosa analysis of the reconstructed structure will be run. when used, the identifier of the protein is used in order to download the original cristalography pdb stucture from the RCSB database. A new folder named PDB is built and the downloaded structure is stored there. The reconstructed structure is stored in the folder FinalComplex/Structures, and in the folder FinalComplex/Analysis the <id>prosa.cmd file is stored along with the pdf image of both the model and target structure prosa energy profiles. Finally, chimera is launched with both model and target 3d representations. 

- **-rnuc** or **- -remove_nucleotides**: if this option is activated once the macro-complex is reconstructed the atoms which do not belong to an amino acid are removed from the final pdb output. This option is useful if we want to do further analysis with Prosa, which does not recognise atoms from nucleotides. To activate this argument the code that needs to be written is ```-rnuc```.

- **-rhet** or **- -remove_heteroatoms**: if this option is activated once the macro-complex is reconstructed if there exists hetero atoms will be removed from the final pdb output. Heteroatoms cause problems when superimposing tructures in ICM. To activate this option the code that needs to be written is ```-rhet```.
  

## **MaCRe algorithm description** 

Macro-Complex Reconstructor is a bioinformatic program written in Python3, which uses a iterative code to build biological macro-molecules out of paired sequences. Here you will fine a description of the procedure that MaCRe follows to achive his goal. 

The function that will generate the macro-complex is **get_macro_complex_iterative()**. However, before usin it, it is necessary to create the arguments that this function will use, the core strucutre and the list of pdbfiles that interact with this core chain.

MaCRe starts with **get_files**, reading all the pdb files found in the default o defined path `-i`. Between those files one will be selected to act as the core chain, around which the rest of the macro-molecule will be constructed. For this pourpose **get_structures()** will be executed, analizing the pdb files found in the `-i` directory and returnig their structures in a list. Once the list with all the pdb structures is generated, the function **define_core_chain()** will check the list for the most frequent structure and will return it in bio.pdb format. This method of choosing the most common chain as the core chain has the objective to speed up the code, reducing the numbers of iterations that will be needed later in the code. However, the user has the possibility of choosing which pdb file will act as the core structure by using the argument `-c examle.pdb`, if the user uses this option he/she can use the name he wants and will not intefere with the rest of the program, and this core structure must be in the same directory as the other structures..

Once the core structure is defined and the list with all the strucures is created the reconstruction can beging calling out **get_macro_complex_iterative()**. MaCRe will start by choosing which are the structures that will first be added to the core strucure using **get_model_to_add()**. This function will iterate trough all the strucutres(test models) and using **get_best_RMSD()** will calculate the best RMSD (and therefore superimposition) between the core strucuture and the test model (the test model needs to have one of the core model chain in ist structure in order to be superimposed with it). After this, with the test model and its RMSD score, **get_models_to_ad()** will create a dictionary with the best models to add ordered by their RMSD from the lowest score (best fited) to the largest.

The next step that MaCRe will follow is to iterate through the dicctionary with the ordered test models adding them to the core structure aplying the corresponding possition matrix. Moreover, the program will check for each chain of the test models if was already added to the core structure. This is done by looking at the clashes produced when the superimposer matrix is applyed. If the chain has more than 30% of it structure with clashes it is considered as it has already been aded to the core structure. On the other hand, if the chain is not considered to already exist in the model it will be added to it. 

The process of obtaining new models to add is repeated untill all the structures have benn added to the core structure forming the macrocomplex or the maximium iterations are achieved (50). Once the macro-complex is reconstructed it will be saved to a new folder "FinalComplex/Structures/" with the name that the user selected in the output argument `-o example`. 

Moreover, the user has some other options to build his complex. If he adds the optional argument **--stechiometry** he can choose between two options:

  *Enter an integer*: which will be the maximum number of chains added to the complex, if the complex reaches this number of chains the algorithm will stop. And, if the user, for example, selects a number of chains that the algorithm is not able to reach by using the **get_best_RMSD()**, which by default is considered to be when the program has done 5 iterations, the program will use the function **get_all_RMSD()** which, instead of looking for the best RMSD of each interaction, stores all the superimpositions and returns a list with the best 100, mixing all the chains. In this way, if a chain is not possible to be added because with the superimposition that the program did at the beggning there are clashes, the program will use other superimpositions of the structure until it can fitt it in the complex. If the program reaches the maximum number of iterations the program stops even if not all the chains have been added.
  
  *Enter a file*: This file will be a file (.txt) with the uniprot code of each protein that the user wants to add to the structure and the number of times it has to be added. Each line will be a protein and the format is the following : "UNIPROT_CODE: N". The program will process this file and store the stechiometry as a hash, then it will do a similar process as when the user enters an integer, but this time each time a chain is added the program stores in a dictionary the number of times each protein has been added. After that, the program compares if the protein has been added the number of times the user wants and, if so, this protein will not be added anymore. If the program does 5 iterations and there is some chain still to be added, it will use the same process as explained for wehn an integer is provided, using the **get_all_RMSD()** function. 

critical poinst:

limitations of the program:


## **Examples** 

### **Example 1: 1gzx**

The protein "1gzx" is the oxy T state haemoglobin from *Homo sapiens*, a tetramer which function is to transport oxygen. The reconstruction of the model whas generated by:
```
python3 ComplexReconstruction.py -i /path/to/folder/1gzx -rnuc -o 1gzx
```
The arguments used in this script are: -i, this flag points the program to the folder where the binary interactions files are stored;
-rnuc, this flag remove the nucleotides from the reconstruction if they exist so we can perform a Prosa analysis; -o, this flag sets the name of the output reconstructed pdb file which will be stored in the /path/to/folder/FinalComplex/Structures

Using Prosa the energy profile is ploted. In red is found the enrgy of the reconstructed macromolecule while in black the profile of the model. Whith a win size of 50, most of the energies are below cero, except for the final part. 
```
read pdb 1gzx.pdb Model
read pdb 1gzxReconstruction.pdb reconstruction
analyse energy *
winsize * 50
shift reconstruction 5
color comb reconstruction red
color comb Model black
graph title 1gzx protein reconstruction energy analysis
plot
```
Using ICM it is possible to superimpose the models and calculate the RMSD. In this case the RMSD value is 2,3181; this is below 10, therefore the recostructed model conformation is near the native conformation. This Prosa analysis could have been done also using the argument -a.
```
read pdb "1gzx.pdb"
read pdb "1gzxReconstruction.pdb"
superimpose a_1. a_2. align 
Rmsd(a_1. a_2. align)
```
The visualizatoin of the complex is made with chimera. MaCRe launches chimera to visualize both model and reconstruction structure if the -a flag is set in the command

<img src="https://github.com/11carlesnavarro/PythonProject/blob/main/images/1zgxcompare-pdb.png"  height="300"><img src="https://github.com/11carlesnavarro/PythonProject/blob/main/images/1gzxColor-1.png"  height="300">

### **Example 2: 5fj8**

The complex "5fj8" is a Cryo-EM structure of yeast RNA polymerase III elongation complex at 3. 9 A from *Saccharomyces cerevisiae* whith transcription function. The reconstruction of the model was generated by:
```
python3 ComplexReconstruction.py  -i /path/to/folder/5fj8 -rnuc -rhet -o 5fj8
```
The arguments used in this script are: -i, this flag points the program to the folder where the binary interactions files are stored; -rnuc, this flag remove the nucleotides from the reconstruction if they exist so we can perform a Prosa analysis; -rhet, this flag remove the eteroatoms, which might cause problems if ICM is used with the reconstructed molecule; -o, this flag sets the name of the output reconstructed pdb file which will be stored in the /path/to/folder/FinalComplex/Structures

Using Prosa the energy profile is ploted. In red is found the enrgy of the reconstructed macromolecule while in black the profile of the model. Whith a win size of 50, most of the energies are below cero, except for some peaks. This Prosa analysis could have been done also using the argument -a.
```
read pdb 5fj8.pdb Model
read pdb 5fj8Reconstruction.pdb reconstruction
analyse energy *
winsize * 50
shift reconstruction 20
color comb reconstruction red
color comb Model black
graph title 5fj8 protein reconstruction energy analysis
plot
```
Using ICM it is possible to superimpose the models and calculate the RMSD. In this case the RMSD value is 74,4199; this is over 10, therefore the recostructed model conformation is far from the native conformation. In this case to calculate the RMSD was necessary to remove the hetero atoms (HETATM) from the reconstructed pdb using the argument -rhet.
```
read pdb "5fj8.pdb"
read pdb "5fj8Reconstruction.pdb"
superimpose a_1. a_2. align 
Rmsd(a_1. a_2. align)
```
The visualizatoin of the complex is made with chimera.

<img src="https://github.com/11carlesnavarro/PythonProject/blob/main/images/5fj8compared-pdb.png"  height="300"><img src="https://github.com/11carlesnavarro/PythonProject/blob/main/images/5fj8Color-1.png"  height="300">

### **Example 3: 2O61**

The complex 2O61 is a macrocomplex formed by NFkB, IRF7, IRF3 bound to the interferon b-enhancer  whith transcription function. The reconstruction of the model whas generated by:
```
python3 ComplexReconstruction.py -i /path/to/folder/2O61 -s stechiometry.txt -c infb_dna.pdb -o out -rhet
```
The visualization of this complex was amde by Chimera, comparing the reconstructed macro-molecule with the original 2o61.pdb file.

<img src="https://github.com/11carlesnavarro/PythonProject/blob/main/images/2o61comparison.png"  height="300">

### **Further analysis** 

 Based on this information, obtained by analizing the atomic conformation of the provided files it is possible to take further decitions on the modeling of the reconstructed structures. One way to follow the analysis is to fix some of the energetic strains in the energetic profiles through the use of Modeller software. Still, we have to take into account that the prosa energy profiles obtained in this analysis were c
arried out against the original structures, which have the same profiles as our reconstructions, so it is
not likely for other proteins more distantly related to have better energy profiles than the ones here dep
icted.
  On the other hand, we have seen that the complex reconstruction of 2O61 have been successful, but the
 best way to follow the analysis would be to analize by docking each pair of relationships between the str
uctures involved. As we have seen, the christalographic data in the repositories is not always harmonic wh
en analizing its energy profile, so it is desireable to explore the possible tridimentional relationships
between these structures.

## **References**

- Garcia-Garcia, J. et al. Networks of ProteinProtein Interactions: From Uncertainty to Molecular Details. Mol Inform 31, 342-362, doi:10.1002/minf.201200005 (2012).

- Russell, R. B. et al. A structural perspective on protein-protein interactions. Curr Opin Struct Biol 14, 313-324, doi:10.1016/j.sbi.2004.04.006S0959440X04000739 [pii] (2004).

- Tuncbag, N., Gursoy, A., Nussinov, R. & Keskin, O. Predicting protein-protein interactions on a proteome scale by matching evolutionary and structural similarities at interfaces using PRISM. Nat Protoc 6, 1341-1354, doi:nprot.2011.367 [pii] 10.1038/nprot.2011.367 (2011).

- Lensink, M. F. et al. Prediction of homo- and hetero-protein complexes by protein docking and template-based modeling: a CASP-CAPRI experiment. Proteins, doi:10.1002/prot.25007 (2016).

- Comeau, S. R., Gatchell, D. W., Vajda, S., & Camacho, C. J. (2004). ClusPro: a fully automated algorithm for protein–protein docking. Nucleic acids research, 32(suppl_2), W96-W99.

- Fiser, A. (2010). Template-based protein structure modeling. Computational biology, 73-94.

- Srivastava, A., Tiwari, S. P., Miyashita, O., & Tama, F. (2020). Integrative/hybrid modeling approaches for studying biomolecules. Journal of molecular biology, 432(9), 2846-2860.
