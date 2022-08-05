# pyDAEDALUSX
The purpose of this software is to design nucleic acid-scaffolded wireframe origami with double duplex edges. With a PLY file as geometry input and with or without a TXT file as the scaffold sequence input, the algorithm calculates scaffold routing along the edges of the geometry and outputs the staple sequences required to fold the scaffold nucleic acid into the target geometry.

There are two ways to submit jobs to the algorithm:

 - Python batch scripts
 - Graphical user interface (GUI) *(Windows only)*

Both require the program to be started as a server in the background.

## Installation

**pyDAEDALUSX download**
Installation time is negligible, as no installation is necessary beyond downloading the required files. Simply clone the GitHub repository or download and unzip the files into your preferred installation location. 

The pyDAEDALUSX repository contains these instructions, an example batch submission script, a "PLY_Files" subfolder with example geometry file inputs, and a "pyDAEDALUSX" subfolder containing the program files. 

A GUI application for Windows users is available upon request.

 **System requirements:** To run pyDAEDALUSX, you minimally need to have Python 2 installed. We recommend setting up a **virtual environment with Python v2.7** in Anaconda (instructions for doing so can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)).
Briefly, in terminal or Anaconda Prompt (recommended for Windows users), type:

    conda create –name python2 python=2.7

to create the environment, named "python2" in this example.
To activate the environment:

    conda activate python2
   
Your command prompt prefix should then look like this:

       (python2) $

Several Python packages are required:

 - networkx==1.11 
 - numpy>=1.11.1
 - matplotlib>=1.5.1
 - scipy>=0.18.0
 - click>=6.6
 - mpmath>=0.19
 - mock>=2.0.0
 - tqdm>=4.9.0
 - flake8>=3.2.0

Use conda to install the required packages in the active python2 virtual environment, e.g.:

    conda install network==1.11 click mpmath mock tqdm flake8

   > We have tested pyDAEDALUSX on Windows 10 and Windows 11 desktops and laptops, and a MacOS 10.15 (Catalina) laptop.

## Running the program
To operate pyDAEDALUSX, you must first start the backend python process, and then submit design jobs either through the GUI or through a python script.

### Start the backend server for pyDAEDALUSX
In terminal or Anaconda Prompt, with your virtual environment activated and operating Python v2.7 as described above, navigate to the pyDAEDALUSX repository folder, wherever you installed it. For example:

    cd C:/users/username/pyDAEDALUSX
    
 Start the backend process with Python by calling the "DAEDserve.py" script in the pyDAEDALUSX subfolder:

    python pyDAEDALUSX/DAEDserve.py
  
  >Note: output folders for each design will be created in the working directory from which you enter the above command. If you want these folders to be output elsewhere, navigate to the desired output directory first, and then specify the full pathname to ./pyDAEDALUSX/pyDAEDALUSX/DAEDserve.py when you start the backend process with Python.

Leave this running in the background (you can minimize the command prompt window) while you submit design jobs using one of the options below.
  
### Submit design jobs

#### Inputs
There are several required inputs and one optional input to generate a 3D nucleic acid-scaffolded DX wireframe origami design with pyDAEDALUSX:

 - Project name (a name for your design)
	 *This will be the name of the output folder created.*
 - Helical form (A-form or B-form)
	 *Note: select A-form if your scaffold and/or staples will be RNA.*
 - Helical turns (minimum edge length, must be an integer)
	 *The minimum edge length will be (# of helical turns) x (bp/turn).*
	 *A-form helices have 11 bp/turn. B-form helices have approximately 10.5 bp/turn.*

	* *Note: Designs using A-form helices must have a minimum edge length of at least 4 helical turns (44 bp). Designs using B-form helices may have a minimum edge length as low as 3 helical turns (31 bp).*
- Target geometry (PLY file)
	*The PLY format is a common Computer Aided Design (CAD) file format; read about it [here](https://en.wikipedia.org/wiki/PLY_(file_format)). Many sample geometries are included in PLY format in the subfolder "PLY_Files" in the repository.*

- *Optional*: Scaffold sequence (TXT file containing only the sequence. A, C, T, G, U are all permissible).
	*If you do not provide a scaffold sequence file, by not selecting a file in the GUI or by specifying "M13.txt" as the scaffold sequence input when submitting a job with Python, the program will use M13mp18 phage sequence for scaffold lengths less than 7,249 nt and will generate random sequence for larger scaffold lengths.*

Longer edge lengths and geometries with more edges will require longer processing times. Design of an A-form regular tetrahedron ("01_tetrahedron.ply" in the PLY_Files subfolder) with 4 helical turns per edge takes about 3s on a typical desktop running Windows 11.
	
#### Option 1: Command Line Interface (CLI)

In a seperate terminal or Anaconda Prompt window from the one running the backend server, start Python:

    python

Then, within Python, import the xmlrpc.client package and define the server parameters and multicall function as follows:

    >>> import xmlrpc.client
    >>> proxy = xmlrpc.client.ServerProxy("http://localhost:4242/")
    >>> multicall = xmlrpc.client.MultiCall(proxy)

Then use the specify the inputs in the order listed above in `multicall.calc()`, and submit the design job to the backend server by running `multicall()`. For example:

    >>> multicall.calc('TestProject','Aform', 4, 'PLY_Files/01_tetrahedron.ply', 'M13.txt')
    >>> result = multicall()

#### Option 2: Batch jobs with a Python script

Similar to the CLI job submission above. You can submit multiple jobs by running the `multicall.calc()` function multiple times with different input parameters before running `result = multicall()`. This can be done either in command prompt running Python as above, or by writing and running a Python script.

An example batch job submission script, "example_submission_script.py" is included in the pyDAEDALUSX repository. In this example, the script loops through all geometry files in the directory "PLY_Files" and designs an A-form 3D wireframe origami object with minimum 4 helical turns per edge for each geometry, using default scaffold sequence.
  
 The script can be run with Python in a command prompt window (after navigating to the appropriate working directory containing your script), e.g.:
 
    python example_submission_script.py

 Or from your favorite Python editor.

#### Option 3: (Windows only) Graphical User Interface

If you have a Windows OS, we can send you a folder containing a GUI upon request.
To use the GUI, run within the shared "pyDAEDALUSX-win32-x64" folder, run the "pyDAEDALUSX.exe" file. A GUI window should appear. Provide the inputs (described above) in the appropriate locations. The "Select PLY" and "Select Sequence" buttons will open a File Explorer for you to select the appropriate file. Do not select a sequence file if you wish to use the program default scaffold sequence.

Click the "submit" button to submit the job to your backend server and generate the design. At the bottom of the GUI window, the word "Processing.." will show until the job is complete, at which time it will display "Fail!" (*see troubleshooting section below*) or "Success!" 

### Outputs

For each design job submitted, a folder titled with your input Project Name will be created in the working directory from which you called the backend server DAEDserve.py with Python. Important outputs include:
- A **.csv** file containing the sequences of all staples necessary to fold the design, as well as the scaffold sequence used.
- A **.cndo** file that describes the predicted nucleotide positions and topology of the output design, using 3DNA notation (for more information, see [here](https://cando-dna-origami.org/cndo-file-converter/)).
- An atomic model of the approximate predicted 3D structure, as **.pdb** files (several model type variants are output). These can be viewed in UCSF Chimera, pyMOL, or similar atomic model viewer, and can be used as inputs into molecular dynamics simulations.
	- [design name and type]_[date]**.pdb**: Defined as a single model composed of many chains. Each strand (scaffold and staples) has its own chain ID. The scaffold is Chain A.
	*For large designs with more staple strands than available alphanumeric chain IDs, chain IDs will be repeated.*
	- [design name and type]_[date]-**multimodel.pdb**: Defined as one model composed of many sub-models. Each strand (scaffold and staples) has its own sub-model number (#0.1, #0.2, etc.). The scaffold is sub-model #0.1.
	- [design name and type]_[date]-**segid.pdb**: Defined as a single model with one chain (pseudo-connected between strands). Each strand is defined as a segment of the chain.
	> Note: By default, atomic models for A-form designs are generated with RNA scaffold and DNA staples.
	
	> Note: For large designs that exceed the PDB file format limit of 99,999 atoms, atom numbers 100,000 + are designated with a hybrid base-36 encoding. These files are still viewable in UCSF Chimera, and typically work with common MD software.
	
- **.png** files showing the distribution of edge lengths, as well as whether some edges needed to be rounded to the nearest 10.5 or 11 bp (slightly distorting the ratio of edge lengths in the input geometry file to ensure all edges are composed of an integer number of full helical turns)

## Troubleshooting

If a design fails (identifiable by a “Fail!” note at the bottom of the GUI, or by missing output files, or by missing strand complements in the .pdb file), make sure:
- Your input scaffold sequence is long enough for the geometry and minimum edge length you have specified.  
	- *Try submitting the job with the same inputs but no scaffold sequence (so that the default will be used). If the design is successful, note the length of the default scaffold sequence used (in the “staples_[…].csv” file in the output folder, the final sequence entry corresponds to the scaffold sequence).*
- You have used at least 3 helical turns as the minimum edge length for B-form designs, or at least 4 helical turns as the minimum edge length for A-form designs.
- Your PLY file is formatted correctly (http://paulbourke.net/dataformats/ply/).  
	- *Try opening the .ply file in a 3D Viewer.*
	- *Try generating a design with one of the included sample PLY file inputs.*

