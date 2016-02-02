## BROMOC Tutorial ##

This tutorial will guide you to learn how to prepare the necessary files to run BROMOC simulations.
The case example of alpha-hemolysin and single stranded DNA olygomer (40 Adenine nucleotides) is presented in this tutorial. 
However, BROMOC suite allows you to employ any kind of ion channel.  

In order to get a solid understanding on how BROMOC works, consult the original papers shown below. If you use BROMOC suite, please, cite us.

- Pablo M. De Biase, Carlos J. Solano, Suren Markosyan, Luke Czapla and Sergei Noskov. BROMOC-D: Brownian Dynamics/Monte-Carlo Program Suite to Study Ion and DNA Permeation in Nanopores. Journal of Chemical Theory and Computation (2012) 
 
- Pablo M. De Biase, Suren Markosyan and Sergei Noskov. Microsecond Simulations of DNA and Ion Transport in Nanopores with Noverl Ion-Ion and Ion-Nucleotides Effective Potentials. Journal of Computational Chemistry (2014) 

For complete a description of BROMOC usage, please refer to the "Documentation" section provided with the code tarball (freely available to download from http://noskovlab.com/downloads/bromoc.tar.bz2)
 
The software needed to follow this tutorial:
 - BROMOC (http://www.noskovlab.com)
 - EFPOT/ECP (http://www.noslovlab.com)
 - CHARMM (http://www.charmm.org)
 - VMD (http://www.ks.uiuc.edu/Research/vmd/)

Recommended OS:
 - Linux

## Extract this tutorial
        $ tar xvf bromoc-tutorial.tar.bz2
        $ cd bromoc-tutorial
        $ vi TUTORIAL

## Download and place Effective Potentials for DNA-KCl into efpot folder
        $ wget http://noskovlab.com/downloads/ep-ecp-dna-kcl.tar.bz2 # do this inside bromoc-tutorial folder
	$ tar xvf ep-ecp-dna-kcl.tar.bz2
	$ mv ep-ecp-dna-kcl efpot  #rename folder to efpot

## Install BROMOC (in Linux)
#       * GNU or Intel compilers must be installed
	$ wget http://noskovlab.com/downloads/bromoc.tar.bz2
	$ tar xvf bromoc.tar.bz2
        $ cd bromoc
        $ ./install       # for GNU compilers or
        $ ./install intel # for intel compilers
	$ export PATH=$PWD/bin:$PATH
	$ echo 'export PATH='$PWD'/bin:$PATH' >> ~/.bashrc


## Prepare protein structure of interest. In this example: staphylococcal alpha-hemolysin ##
	Get the structure of staphylococcal alpha-hemolysin, a heptameric transmembrane pore. 
        The pdb file named 7AHL can be downloaded from the Protein Data Bank (www.rcsb.org) (Direct link: http://dx.doi.org/10.2210/pdb7ahl/pdb). 
	Make you sure to adjust the position and orientation according to your needs. Use "transcrd" from BROMOC tools (it needs NOT Extended crd). 
        Also add hydrogens and remove waters from the crystal structure (if needed). 

        For this tutorial we have already included pre set-up structure "wt.pdb" in folder "crd"

 	Navigate to "crd" folder
	$ cd crd
    
        Convert pdb to crd (crd is CHARMM coordinates format) using "pdb2crd" from BROMOC tools
        $ pdb2crd wt.pdb wt.crd n

## Prepare Charges, PB and Stern Radii files using CHARMM ##

        In "crd" folder...
	
	Use "chr_pb_stern.inp" script to compute the following coordinate files using CHARMM:
        $ charmm < chr_pb_stern.inp > chr_pb_stern.out &	
        Note that this process may take several hours

        Output files:
	Charges "wt_chr.cor"
	PB Radii "wt_pb.cor"
	Stern Radii "wt_stern.cor" 

## Compute reaction field, static field and repulsive field maps using pb-pnp ##

	Navigate to "maps" folder
        $ cd ../maps
	Make yourself familiar with the contents of this directory

	Maps will have extension ".pbeq"

	To generate reaction field map do:
        $ pb-pnp < rf.inp > rf.out 
        Input: "rf.inp" 
        Output:"rf.pbeq","rf.out"
        
	To generate Static field with explicit ions map do: 
        $ pb-pnp < st.inp > st.out 
	Input: "st.inp"
	Output: "st.pbeq","st.out"

	To generate Static field with implicit ions map do:
        $ pb-pnp < st_noion.inp > st_noion.out 
	Input: "st_noion.inp"
	Output: "st_noion.pbeq","st_noion.out"
        
	To generate Repulsive field map do:
        $ pb-pnp < rep.inp > rep.out 
        Input: "rep.inp"
	Output "rep.pbeq","rep.out"

        Note that these processes may take several hours.

## Generate ss-DNA and set it up to fit alpha-hemolysin pore ##

	Navigate to "dna" directory
        $ cd ../dna
	Make yourself familiar with the contents of this directory

	Execute BROMOC and run "generate_dna.inp" script to generate a single stranded poly cytosine oligomer (40 bases).
        The DNA will be stretched applying forces at both ends and compressed laterally by using a narrow phantom pore.
        $ bromoc < generate_dna.inp > generate_dna.out
        Output: "dna-0000.bco","dna-0000.xyz","generate_dna.out" 
	Final structures are generated in two formats ".bco" (BROMOC coordinates) and ".xyz" (Standard XYZ)

## Generate a random ss-DNA structure doing a "heatup" protocol ##
        In this section we are going to do brownian dynamics of the previous structure
        at high temperature (330 K) and high diffusivity (0.1) to generate a random structure
        that will be used later to do simulations.
        $ bromoc < heatup_dna.inp > heatup_dna.out
        Ouput: "dna.bco", "heatup_dna.out"

## Run current measurement in BROMOC ##
        
        In this example we are going to use the previously heated-up DNA. First we are going to fill the system with ions using 
        GCMC (Gran Canonical Monte Carlo) and MCM (Monte Carlo Move). Then we will do a short equilibration using BD (Brownian Dynamics).
        Finally, a 0.1 us production BD run and measure the ion current. In this example we use effective potentials contained in 
        folder "efpot" (prepared above). DNA-KCl effective potential set were developed by Noskov Lab. Check our website: http://www.noskovlab.com
 
	First, navigate to "testrun" directory
        $ cd ../testrun
	Make yourself familiar with the contents of this directory

	Check "run-0000.inp" input file
	- Make sure to have the proper system box size (should exactly match the box size you used to calculate the maps).
	- Make sure to have the proper ionic concentration and the corresponding electrochemical potentials.
	- Make sure to have the proper voltage (it should exactly match the one you use to calculate the "st.pbeq").
	- Double check the paths for each supplementary file (i.e. effective potentials, maps and dna).
	
	Please consult the "Documentation" section for the details on syntaxes used in this input script.

	The output files you get are the following:
		"run-first-0000.bco" A file that contains the initial BROMOC coordinates
		"run-first-0000.xyz" A file that contains the initial coordinates in xyz format
		"run-0000.btr" A BROMOC trajectory file
		"run-0000.cur" A file that contains information on calculated ion current
		"run-last-0000.bco" A file that contains the final BROMOC coordinates
		"run-last-0000.xyz" A file that contains the final coordinates in xyz format

        To run the example:
        $ bromoc < run-0000.inp > run-0000.out

        Note that this process may take several hours.

## Other specific examples ##

	In addition to current measurement, you can also run some additional tasks using BROMOC suite.
	In this tutorial we have included several examples of additional BROMOC features that user can leverage off.

	## Analytical potential ##

	The interaction between two particles are described through two kinds of interaction pairwise potentials in BROMOC: 
        An analytical form and a discrete form. Previously, we used the discrete form which we called "Effective Potentials".
	The analytical form is composed by three analytical terms: a Coulombic, Lennard-Jones and Short-Range Potential of 
	Mean Force term. In this example we are going to use the latter.

	Navigate to "examples/analytical-potential" directory
	$ cd ../examples/analytical-potential

	Execute "run-0000.inp" input.
	$ bromoc < run-0000.inp > run-0000.out

	## Monte Carlo move ##

	This option allows the movement of mobile particles using a low-computational cost (fast) technique by moving ions 
	based on Metropolis Monte-Carlo (MMC) algorithm. MCM generates a particle configuration consistent with Boltzmann 
	probability distribution at the given temperature and in the working system. MCM, instead of computing the microscopic 
	forces, obtains the energy change for a given random move and accept it or discard it according to the MMC criteria.

	Navigate to "examples/monte-carlo-move" directory
	$ cd ../monte-carlo-move

	Execute "mc-0000.inp" input.
	$ bromoc < mc-0000.inp > mc-0000.out

	## Simulations without DNA or open pore simulations ##

	This example allows user to run BROMOC simulations using only a protein channel (no DNA included)

	Navigate to "examples/openpore" directory
	$ cd ../openpore

        Execute "run-0000.inp" input.
	$ bromoc < run-0000.inp > run-0000.out

	## Parallel simulations with immobilized DNA ##

	In this example, we are going immobilize a ss-DNA and run simulations on that system.
	This this example we are going to perform a parallel simulation. The more ions the 
	system has the better parallel scaling is.

	Navigate to "examples/parallel-bd-fixeddna" directory
	$ cd ../parallel-bd-fixeddna

	Execute "bd-0000.inp" input with bromoc-para using 4 processors
	$ OMP_NUM_THREADS=4 bromoc-para < bd-0000.inp > bd-0000.out

	## Simulations with DNA embeded in a phantom pore with implicit ions ##
	
	A membrane can be represented by an analytical repulsive potential (polynomial) applied for a switch region between bulk 
	and membrane regions so that the defined membrane region has a potential to represent an impermeable membrane. In addition, 
	a cylindrical pore with a given radius can also be setup for each ion type. A similar analytical repulsive potential is 
	applied for a switch region between the pore bulk and membrane regions with a potential. 
	
	This example performs simulations using "phantom" cylindrical pore. In this example we are not going to use ions
	explicitly. So ions will be considered in a implicit manner using a Debye-Huckel approximation.

	Navigate to "examples/phantom-pore-dna" directory
	cd ../phantom-pore-dna

	Execute "run-0000.inp" input with BROMOC.
	$ bromoc < run-0000.inp > run-0000.out

	## Simulations with DNA embeded in a phantom pore and in the presence of salt solution (KCl) ##

	This example shows how to perform simulations using phantom pore, explicit ions (KCl) and Analytical Potentials.

	Navigate to "examples/phantom-pore-dna-ions" directory
	cd ../phantom-pore-dna-ions

	Execute "run-0000.inp" input with BROMOC.
	$ bromoc < run-0000.inp > run-0000.out

	## Plot (static field, reaction field and repulsion field) maps ##

	This script plots reaction field, static field and repulsive field maps

	Navigate to "examples/plot-maps" directory
	$ cd ../plot-maps
	Execute "plotmap.inp" input with BROMOC.
	$ bromoc < plotmap.inp > plotmap.out &

	## Generate ions for a simulation box ##

	This example populates a system with ions.

	Navigate to "examples/populate-with-ions" directory
	$ cd ../populate-with-ions
	
	Run "bromoc-ions-gener.inp" BROMOC input
	$ bromoc < bromoc-ions-gener.inp > bromoc-ions-gener.out &

	## Generate a trajectory with a Spherical System and Compute RDF
	In this example a Brownian Dynamics is produced with DNA and KCl ions in a Spherical System.
	Radial Distribution Function is computed using rdf-bromoc from the resulted trajectory (BROMOC Tools). 

	Navigate to "examples/sph-rdf" directory
	$ cd ../sph-rdf

	Generate ions and equilibrate system. Run "bd-ions-gener.inp" input ang generate "bd-ions-gener.bco"
	$ bromoc < bd-ions-gener.inp > bd-ions-gener.out

	Run Brownian Dynamics and generate trajectory "bd-0001.btr"
	$ bromoc < bd-0001.inp > bd-0001.out

	Compute RDF from the output trajectory using Spherical System Radius 37.22 Ang
	$ rdf-bromoc bd-0001.btr 37.2213301228 0.1
	
	## Visualizing a trajectory ##
	
	Navigate to "testrun" directory
	$ cd ../../testrun

	Convert btr (BROMOC trajectory) to xtc (GROMACS trajectory)
	$ btrconv run-0000.btr 1 y xtc "" 3
	
	Visualize with VMD the trajectory and protein map
	$ vmd -f run-0000.gro run-0000.xtc -f ../maps/rep.pbeq



