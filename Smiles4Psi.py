#!/usr/bin/env python
# coding: utf-8

#
# @BEGIN LICENSE
#
# Smiles4Psi: Automated conformational search and Psi4 conformer inputs
#             generation from SMILES strings
#
# Copyright (c) 2021 
# Carlos H. Borca, Satyen Dhamakar, and Michael A. Webb
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Smiles4Psi.
#
# Smiles4Psi is free software; you can redistribute it and/or modify
# it under the tesms of the GNU Lesser General Public License as
# published by the Free Software Foundation, version 3.
#
# Smiles4Psi is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Smiles4Psi; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA.
#
# @END LICENSE
#

__author__  = "Carlos H. Borca, Satyen Dhamakar"
__version__ = "0.0.4"
__credits__ = ["Carlos H. Borca", "Satyen Dhamakar", "Michael A. Webb"]
__email__   = "cborca@princeton.edu"
__license__ = "LGPLv3"


# Module imports.
import argparse
import os
import subprocess
import sys

# Setting up pretty print for debug printouts.
import pprint
pp = pprint.PrettyPrinter(indent=1)

# Non-standard module imports.
try:
    from openbabel import openbabel # NOTE: Requires the Open Babel python module available through python.

except ImportError:
    openBabelImportError = """
    +--------------------------------+
    | Open Babel Module Import Error | Unable to import the Open Babel module.
    +--------------------------------+
    
    Is Open Babel installed in your system?

    Create a conda environment with the dependencies required by this code.
    First download a version of Anaconda, i.e. Miniconda:
    
      $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    Execute the installer and follow on-screen instructions:
    
      $ bash Miniconda3-latest-Linux-x86_64.sh
    
    Then create a conda environment using the following command:
    
      $ conda create -n S4P psi4=1.3.2 psi4-rt=1.3 python=3.7 openbabel=3.1.0 psutil=5.8.0 -c psi4 -c conda-forge
    
    Once the environment is created created, restart your shell and activate the environment:
    
      $ conda activate S4P

    And execute this program again in the active environment.
    """
    print(openBabelImportError)
    sys.exit()

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def createParser():

    # For the argument parser. An example of usuage and a description of what this script does.
    exas =  "python {} 'CCO' -st Ethanol -ff xyz -ct c.xyz -nc 100 ".format(sys.argv[0])
    exas += "-ab EtOH -fn 'Alcohol' -vb 1 -oe True"
    docs =  "Smiles4Psi4: Automated conformational search with a genetic algorithm interfaced "
    docs += "with a Psi4-input generator for conformers that requires only SMILES as input."
    epis =  "Execution of this code requires the `openbabel` module. It may also need `psutil`."

    # Help strings.
    hsm = "SMILES string input, i.e. 'CCO'"
    hst = "molecular 3D-structure file name, i.e. 'Ethanol' (default: 'Molecule')"
    hff = "molecular 3D-structure file format, i.e. 'xyz' (default: 'xyz')"
    hct = "conformers trajectory full file name, i.e. 'c.xyz' (default: 'Confrmrs.xyz')"
    hnc = "Target number of conformers, i.e. 300 (default: 100)"
    hsf = "Type of scoring function for conformational search, i.e. 'energy' (default: 'rmsd')"
    hab = "Molecule name abbreviation for file names and variables, i.e. 'EtOH' (default: 'MOL')"
    hfn = "Full molecule name for text labels, i.e. 'Ethylic alcohol'"
    hvb = "Level of verbosity: 0 = Error, 1 = Warnings, 2 = Information, 3 = Debug. (default: 0)"
    hpi = "Psi4 psithon input file name without extension, i.e. 'Et-B97D3-aTZ'"
    hmm = "Psi4 memory in gigabytes, i.e. '2'"
    hmt = "Psi4 method specification string, i.e. 'B97-D3/aug-cc-pVTZ'"
    hsc = "Psi4 SCF algorithm type, i.e. 'DIRECT' (default: DF)"
    hmp = "Psi4 MP2 algorithm type, i.e. 'CONV' (default: DF)"
    hcc = "Psi4 CC algorithm type, i.e. 'DF' (default: DF)"
    hfc = "Psi4 freeze core electrons boolean (default: True)"
    hec = "Psi4 energy convergence threshold, i.e. 10 (default: 8)"
    hdc = "Psi4 density convergence threshold, i.e. 10 (default: 8)"
    hoe = "Psi4 one-electron properties bolean (default: True)"
    hcb = "Psi4 cube-based properties bolean (default: False)"

    # Create parser.
    parser = argparse.ArgumentParser(description=docs, epilog=epis, usage=exas)

    # Positional Arguments
    parser.add_argument("sm", metavar="SMILES", type=str, help=hsm)

    # No positional arguments have been specified. Optional arguments the script might expect.
    parser.add_argument("-st", "--structn", dest="st", default="Molecule",     help=hst)
    parser.add_argument("-ff", "--3dformt", dest="ff", default="xyz",          help=hff)
    parser.add_argument("-ct", "--cfilenm", dest="ct", default="Confrmrs.xyz", help=hct)
    parser.add_argument("-nc", "--nconfrs", dest="nc", default=100,            help=hnc)
    parser.add_argument("-sf", "--scorefn", dest="sf", default="rmsd",         help=hsf)
    parser.add_argument("-ab", "--molabrv", dest="ab", default="MOL",          help=hab)
    parser.add_argument("-fn", "--fulmonm", dest="fn", default=None,           help=hfn)
    parser.add_argument("-vb", "--verbose", dest="vb", default=0,              help=hvb)
    parser.add_argument("-pi", "--psi4inp", dest="pi", default=None,           help=hpi)
    parser.add_argument("-mm", "--psi4mem", dest="mm", default=None,           help=hmm)
    parser.add_argument("-mt", "--psi4mtd", dest="mt", default=None,           help=hmt)
    parser.add_argument("-sc", "--psi4scf", dest="sc", default="DF",           help=hsc)
    parser.add_argument("-mp", "--psi4mp2", dest="mp", default="DF",           help=hmp)
    parser.add_argument("-cc", "--psi4cca", dest="cc", default="DF",           help=hcc)
    parser.add_argument("-fc", "--psi4fce", dest="fc", default="True",         help=hfc)
    parser.add_argument("-ec", "--psi4ecv", dest="ec", default=8,              help=hec)
    parser.add_argument("-dc", "--psi4dcv", dest="dc", default=8,              help=hdc)
    parser.add_argument("-oe", "--psi4oep", dest="oe", default="True",         help=hoe)
    parser.add_argument("-cb", "--psi4cbp", dest="cb", default="False",        help=hcb)

    return parser
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def processArgs(args):
    """
    Take the parsed arguments and return usable objects in a dictionary.

    Parameters
    ----------
    args : class 'argparse.Namespace'
        An object containing the runtime arguments.

    Returns
    -------
    keywords : dict
        Set of runtime variables 
    """

    keywords = {} # A dictionary for variable objects.

    keywords['SMILES']               = args.sm
    keywords['3DStructFileName']     = args.st
    keywords['3DStructFileFormat']   = args.ff
    keywords['ConfTrajFileName']     = args.ct
    keywords['NumberOfConformers']   = int(args.nc)
    keywords['ScoringFunction']      = args.sf
    keywords['MoleculeNameAbbrev']   = args.ab
    keywords['MoleculeNameFull']     = args.fn if args.fn else args.st
    keywords['Verbosity']            = int(args.vb)
    keywords['Psi4']                 = {}
    keywords['Psi4']['InputName']    = args.pi if args.pi else args.ab
    keywords['Psi4']['Memory']       = args.mm
    keywords['Psi4']['Method']       = args.mt
    keywords['Psi4']['SCFType']      = args.sc
    keywords['Psi4']['MP2Type']      = args.mp
    keywords['Psi4']['CCType']       = args.cc
    keywords['Psi4']['FreezeCore']   = bool(False if (args.fc.lower() == "false") else True)
    keywords['Psi4']['EConv']        = int(args.ec)
    keywords['Psi4']['DConv']        = int(args.dc)
    keywords['Psi4']['OneElecProps'] = bool(False if (args.oe.lower() == "false") else True)
    keywords['Psi4']['CubeProps']    = bool(False if (args.cb.lower() == "false") else True)

    return keywords
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def printHeader(keywords):
    """
    Prints welcome message and additional execution details.

    Parameters
    ----------
    keywords : dict
        Set of runtime variables 

    Examples
    --------
    >>> printHeader(keywords)
    """

    # Warnings.
    if (keywords['Verbosity'] > 0):
        print("\n{}\n".format(36*"-="))
        print("{}{}".format(31*" ", "Smiles4Psi"))

        subtitle = """
Automated conformational search with a genetic algorithm interfaced with
a Psi4-input generator for conformers that requires only SMILES as input

        -------------------------------------------------------
         Carlos H. Borca, Satyen Dhamakar, and Michael A. Webb
           Department of Chemical and Biological Engineering
               School of Engineering and Applied Science
                          Princeton University
        -------------------------------------------------------"""

        print("{}".format(subtitle))
        print("\n{}\n".format(36*"-="))

    # Info.
    if (keywords['Verbosity'] > 1):
        print("{}".format(72*"="))
        print("Smiles4Psi runtime variables configuration (Keywords dictionary):")
        print("{}".format(72*"-"))
        pp.pprint(keywords)
        print("{}".format(72*"="))
        print("")

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def smilesToMolFile(smiles, molFileName='Molecule', fileFormat='xyz', verbose=0):
    """
    From a SMILES string use Open Babel to write a file with a molecule's 3D
    structure specified by a list of elements and their XYZ coordinates.

    Upon execution, Open Babel will output a file with the molecule's
    structure.

    Parameters
    ----------
    smiles : str
        A SMILES string for the molecule.
    molFileName : str
        Name of the OpenBabel output file with the molecular structure.
    fileFormat : str
        Format of the OpenBabel output.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    molFile : file
        File with a molecular 3D structure.

    Examples
    --------
    >>> smilesToMolFile('C=C(C)C(=O)OCCCC', molFileName='BMA', fileFormat='xyz', verbose=2)
    1 molecule converted
    """

    # Option `--gen3D` specifies the generation of XYZ coordinates.
    # Option `--append 'MW'` appends the molecular weight to the description of the molecule.
    openBabelCommand = "obabel -:'{}' -O {}.{} --gen3D --append 'MW'".format(smiles, molFileName, fileFormat)

    print("\n3D-structure generation:\n{}".format(openBabelCommand)) if (verbose > 1) else None # Info.
    print("\n{} Output from Open Babel {}".format(24*"-", 24*"-")) if (verbose > 2) else None # Debug.

    # Execute the command on the shell and wait until finished. An output file will be produced.
    # NOTE: The following line requires the executable of OpenBabel to be accessible from the shell.
    # TODO: Make sure that the file is output to the correct location!!
    # TODO: How is chirality going to be handled?
    if (verbose > 2):
        openBabelProcess = subprocess.Popen(openBabelCommand,
                                            shell=True)
    else:
        openBabelProcess = subprocess.Popen(openBabelCommand,
                                            shell=True,
                                            stderr=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL)
    openBabelProcess.wait()

    print("{}".format(72*"-")) if (verbose > 2) else None # Debug.

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def molFileToConfTraj(smiles,
                      molFileName='Molecule',
                      fileFormat='xyz',
                      confOutTrajName='Confrmrs.xyz',
                      numConfrmrs=3,
                      scoreFunction='rmsd',
                      verbose=0):
    """
    From a file with a molecule's 3D structure specified by a list of
    elements and their XYZ coordinates, run a conformational search and
    obtain an XYZ trajectory with conformer structures.

    Parameters
    ----------
    molFileName : str
        Name of the OpenBabel output file with the molecular structure.
    fileFormat : str
        Format of the OpenBabel output.
    confOutTrajName : str
        Name of the OpenBabel output trajectory file with conformers.
    numConfrmrs : str
        Number of requested conformers.
    scoreFunction : str
        Type of scoring function used by Open Babel.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    confOutTrajFile : file
        Trajectory-like file with conformers produced by the genetic algorithm search.

    Examples
    --------
    >>> molFileToConfTraj(smiles, molFileName='Molecule', fileFormat='xyz', scoreFunction='rmsd', verbose=2)
    1 molecule converted
    """

    # Option `--conformer` requests a conformational search based on the genetic algorithm.
    # Option `--nconf` specifies the number of conformers to be obtained.
    # Option `--score` specifies the scoring function employed to rank conformers.
    # Option `--writeconformers` outputs conformer structures to a trajectory-like file.
    openBabelCommand = "obabel {}.{} -O {} --conformer --nconf {} --score {} --writeconformers".format(molFileName,
                                                                                                       fileFormat,
                                                                                                       confOutTrajName,
                                                                                                       numConfrmrs,
                                                                                                       scoreFunction)

    print("\nGA conformational search:\n{}".format(openBabelCommand)) if (verbose > 1) else None # Info.
    print("\n{} Output from Open Babel {}".format(24*"-", 24*"-")) if (verbose > 2) else None # Debug.

    # Execute the command on the shell and wait until finished. An output file will be produced.
    # TODO: Make sure that the file is output to the correct location!!
    if (verbose > 2):
        openBabelProcess = subprocess.Popen(openBabelCommand,
                                            shell=True)
    else:
        openBabelProcess = subprocess.Popen(openBabelCommand,
                                            shell=True,
                                            stderr=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL)
    openBabelProcess.wait()

    print("{}".format(72*"-")) if (verbose > 2) else None # Debug.

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def getChargeFromSmiles(smiles, verbose=0):
    """
    From a SMILES string use OpenBabel to determine the molecule's charge.

    Parameters
    ----------
    smiles : str
        A SMILES string for the molecule.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    charge : int
        Molecular charge.

    Examples
    --------
    >>> getChargeFromSmiles('C=C(C)C(=O)OCCCC', verbose=3)
    1
    """

    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles)
    charge = mol.GetTotalCharge()

    #print("\nMolecular charge = {} a.u.".format(charge)) if (verbose > 2) else None # Debug.

    return charge
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def getMultipFromSmiles(smiles, verbose=0):
    """
    From a SMILES string use OpenBabel to determine the molecule's 
    multiplicity.

    Parameters
    ----------
    smiles : str
        A SMILES string for the molecule.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    multiplicity : int
        Molecular multiplicity.

    Examples
    --------
    >>> getMultipFromSmiles('C=C(C)C(=O)OCCCC', verbose=3)
    0
    """

    obConversion = openbabel.OBConversion()
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smiles)
    multiplicity = mol.GetTotalSpinMultiplicity()

    #print("\nMolecular spin multiplicity = {} a.u.".format(multiplicity)) if (verbose > 2) else None # Debug.

    return multiplicity
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def splitTrajToDict(smiles, confOutTrajName='Confrmrs.xyz', molNameAbv='MOL', molNameFull='Molecule', verbose=0):
    """
    From an XYZ trajectory with conformer structures, create a dictionary
    containing the data for each conformer.

    Parameters
    ----------
    smiles : str
        A SMILES string for the molecule.
    confOutTrajName : str
        Name of the OpenBabel output trajectory file with conformers.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    conformers : dict
        Dictionary with conformers structures and information.

    Examples
    --------
    >>> splitTrajToDict('CC', 'EtConf.xyz', verbose=0)

    """

    conformers = {}
    trajLines = []

    with open(confOutTrajName, 'r') as tF:
        for l in tF:
            trajLines.append(l)

    numAtoms = int(trajLines[0])
    confIdx = -1 # Indexing of conformers start at 0.

    for idx, line in enumerate(trajLines):

        if (idx%(numAtoms + 2) == 0):
            confIdx += 1

            conformers[confIdx] = {}                                                      # Conformer's dictionary.
            conformers[confIdx]['MoleculeName']    = molNameFull                          # Conformer's molecule name.
            conformers[confIdx]['MoleculeAbv']     = molNameAbv                           # Molecule name abbreviation.
            conformers[confIdx]['Index']           = confIdx                              # Conformer's index.
            conformers[confIdx]['SMILES']          = smiles                               # Conformer's SMILES string.
            conformers[confIdx]['NumberOfAtoms']   = numAtoms                             # Total number of atoms.
            conformers[confIdx]['Charge']          = getChargeFromSmiles(smiles, verbose) # Molecular charge.
            conformers[confIdx]['Multiplicity']    = getMultipFromSmiles(smiles, verbose) # Molecular spin multiplicity.
            conformers[confIdx]['MolecularWeight'] = None                                 # Molecular weight (below).
            conformers[confIdx]['Atoms']           = {}                                   # Dictionary conformer atoms.

        atomIdx = (idx - 2)%(numAtoms + 2)                                                # First two lines are extras

        if (idx%(numAtoms + 2) == 1):
            conformers[confIdx]['MolecularWeight'] = float(line[:-1])

        if (len(line.split()) >= 3):

            # A dictionary for the element and coordinates of each atom.
            atom = {}

            elemsAndCoords  = line.split()
            atom['e'] = elemsAndCoords[0]
            atom['x']       = float(elemsAndCoords[1])
            atom['y']       = float(elemsAndCoords[2])
            atom['z']       = float(elemsAndCoords[3])

            # Append the new atom to the conformer dictionary.
            conformers[confIdx]['Atoms'][atomIdx] = atom

    # Info.
    if (verbose > 1):
        print("\nConformers dictionary now contains {} entries".format(len(conformers)))

    return conformers
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def getAvailMemFromPsutil(verbose=0):
    """
    Imports the non-standard psutil module to return the amount of available
    system memory in bytes.

    Parameters
    ----------
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    memInBytes : int
        Current available system memory in bytes.

    Examples
    --------
    >>> getAvailMemFromPsutil(verbose=3)
    Current available system memory (RAM) = 3359.28 MB
    """

    try:
        import psutil

    except ImportError:
        psutilImportError = """
    +----------------------------+
    | PSUtil Module Import Error | Unable to import the PSUtil module.
    +----------------------------+

    Is PSUtil installed in your system?
    Create a conda environment with the dependencies required by this code.
    First download a version of Anaconda, i.e. Miniconda:

      $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    Execute the installer and follow on-screen instructions:

      $ bash Miniconda3-latest-Linux-x86_64.sh

    Then create a conda environment using the following command:

      $ conda create -n S4P psi4=1.3.2 psi4-rt=1.3 python=3.7 openbabel=3.1.0 psutil=5.8.0 -c psi4 -c conda-forge

    Once the environment is created created, restart your shell and activate the environment:

      $ conda activate S4P

    And execute this program again in the active environment.
    """
        print(psutilImportError)
        sys.exit()

    svmem      = psutil.virtual_memory()
    memInBytes = svmem.available

    print("\nAvailable memory (RAM) = {:.2f} MB".format(memInBytes/(1024*1024))) if (verbose > 1) else None # Info.

    return memInBytes
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def conformerToPsithon(conformer,
                       totalNumConf,
                       memStr,
                       psi4,
                       verbose=0):
    """
    Writes out a new folder in the current folder containing a Psithon input
    file for one provided conformer. Multiple options are user-accessible.

    Parameters
    ----------
    conformer : dict
        Dictionary with the structure and information for one conformer.
    totalNumConf : int
        Total number of structures in the dictionary of conformers.
    memStr : str
        String containing the description of the Psi4 memory setup.
    psi4 : dict
        Dictionary with Psi4 setup information.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    p4InFolder/p4InFile : file
        Psithon input file for the conformer.

    Examples
    --------
    >>> conformerToPsithon(conformer, totalNumConf, psi4)
    """

    # Title and credits
    p4In =  "# Psi4 input file produced by Smiles4Psi\n\n"
    p4In += "# Credits: Carlos H. Borca, Satyen Dhamakar, and Michael A. Webb\n"
    p4In += "#          Department of Chemical and Biological Engineering\n"
    p4In += "#          School of Engineering and Applied Science\n"
    p4In += "#          Princeton University\n"
    p4In += "#          Princeton, New Jersey, USA\n\n"

    # Molecule details
    p4InFileName = "{}-{}.in".format(psi4['InputName'], str(conformer['Index']).zfill(len(str(totalNumConf))))

    p4In += "# Psithon input name         = {}\n"      .format(p4InFileName                )
    p4In += "# Molecule full name         = {}\n"      .format(conformer['MoleculeName']   )
    p4In += "# Molecule name abbreviation = {}\n"      .format(conformer['MoleculeAbv']    )
    p4In += "# SMILES descriptor          = {}\n"      .format(conformer['SMILES']         )
    p4In += "# Number of atoms            = {}\n"      .format(conformer['NumberOfAtoms']  )
    p4In += "# Molecular weight           = {} g/mol\n".format(conformer['MolecularWeight'])
    p4In += "# Molecular charge           = {} a.u.\n" .format(conformer['Charge']         )
    p4In += "# Molecular spin multiplicty = {}\n"      .format(conformer['Multiplicity']   )
    p4In += "# Conformer index            = {}\n"      .format(conformer['Index']          )

    # Memory
    p4In += "\nmemory {}\n".format(memStr)

    # Molecule name (string used must be acceptable for variable creation in Python)
    p4MolName = "{}{}".format(conformer['MoleculeAbv'], str(conformer['Index']).zfill(len(str(totalNumConf))))

    normalChars = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                   "T", "U", "V", "W", "X", "Y", "Z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "_"]

    for i in p4MolName.upper():
        if not (i in normalChars):
            p4MolName = p4MolName.replace(i, '')

    p4In += "\nmolecule {} {{\n".format(p4MolName)

    # Charge and multiplicity
    p4In += "{} {}\n".format(conformer['Charge'], conformer['Multiplicity'])

    # Elements and coordinates
    for key in conformer['Atoms'].keys():
        p4In += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(conformer['Atoms'][key]['e'],
                                                              conformer['Atoms'][key]['x'],
                                                              conformer['Atoms'][key]['y'],
                                                              conformer['Atoms'][key]['z'])
    p4In += "}\n\n"

    # Creation of output directories # TODO: Check that this works for multiple conformers.
    if (psi4['CubeProps'] == True):
        p4In += "# Get current active molecule object, and its name.\n"
        p4In += "# Then, use that name to create adirectory for CUBE files.\n"
        p4In += "activemol  = psi4.get_active_molecule()\n"
        p4In += "moldirname = psi4.core.Molecule.name(activemol)\n"
        p4In += "actmolcwd  = os.getcwd()\n"
        p4In += "outpropdir = os.path.join(actmolcwd, moldirname)\n\n"
        p4In += "try:\n"
        p4In += "    os.mkdir(outpropdir)\n\n"
        p4In += "except FileExistsError:\n"
        p4In += "    pass\n\n"

    # Control settings
    p4In += "# Psi4 control settings:\n"
    p4In += "set {\n"
    p4In += "  scf_type          = {}\n".format(psi4['SCFType'])
    p4In += "  mp2_type          = {}\n".format(psi4['MP2Type'])
    p4In += "  cc_type           = {}\n".format(psi4['CCType'])
    p4In += "  freeze_core       = {}\n".format(psi4['FreezeCore'])
    p4In += "  e_convergence     = {}\n".format(psi4['EConv'])
    p4In += "  d_convergence     = {}\n".format(psi4['EConv'])

    if (psi4['CubeProps'] == True):
        p4In += "  cubeprop_filepath = $outpropdir\n"
        p4In += "  cubeprop_tasks    = ['frontier_orbitals', 'density',\n"
        p4In += "                       'ESP', 'dual_descriptor']\n"

    p4In += "}\n\n"

    # Execution statement
    p4In += "# Run a single-point energy calculation.\n"
    p4In += "# Return total energy and wavefunction object.\n"
    p4In += "E, wfn = energy('{}', return_wfn=True)\n\n".format(psi4['Method'])

    # Evaluation of one-electron properties
    if (psi4['OneElecProps'] == True):
        p4In += "# Print ground-state one-electron molecular properties\n"
        p4In += "oeprop(wfn, 'DIPOLE', 'QUADRUPOLE', 'MULTIPOLE(4)', 'MO_EXTENTS',\n"
        p4In += "       'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'WIBERG_LOWDIN_INDICES',\n"
        p4In += "       'MAYER_INDICES', 'NO_OCCUPATIONS', title='{}')\n\n".format(conformer['MoleculeName'])

    # All variables printout
    p4In += "# Print all output variables.\n"
    p4In += "print_variables()\n"

    # Directives to write CUBE files (last in case it fails)
    if (psi4['CubeProps'] == True):
        p4In += "\n# Write CUBE files as specified by the cubeprop_tasks settings.\n"
        p4In += "cubeprop(wfn)"

    # Debug.
    if (verbose > 2):
        print("{}\n".format(72*"."))
        print("{}".format(p4In))
        print("{}".format(72*"."))

    # Write Psithon input to a new directory
    owd = os.getcwd()
    p4InFolder = psi4['InputName']

    try:
        os.mkdir(p4InFolder)

    except FileExistsError:
        pass

    os.chdir(p4InFolder)

    with open(p4InFileName, 'w') as p4InFile:
        for line in p4In:
            p4InFile.write(line)

    os.chdir(owd)

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def buildPsithonsFromConfDict(conformers, psi4, verbose=0):
    """
    Calls the conformerToPsithon() function to build Psithon input files for
    every conformer stored in a provided conformers dictionary.

    Parameters
    ----------
    conformers : dict
        Dictionary with conformers structures and information.
    psi4 : dict
        Dictionary with Psi4 setup information.
    verbose : int
        Adjusts the level of detail of the printouts.

    Returns
    -------
    p4InFolder : directory
        A new directory in the current location with Psithon inputs.

    Examples
    --------
    >>> buildPsithonsFromConfDict(conformers)
    """

    # Determine available memory if not specified by the user.
    if not psi4['Memory']:
        import math
        memInBytes = getAvailMemFromPsutil(verbose)
        memInGigaB = memInBytes/(1024*1024*1024)
        memStr = "{} GB".format(math.floor(memInGigaB))

    else:
        memStr = "{:.1f} GB".format(float(psi4['Memory']))

    # For style purposes
    print("") if (verbose > 1) else None # Info.
    print("{}".format(72*"-")) if (verbose > 2) else None # Debug.

    # Generate a psithon input for each conformer.
    for key in conformers.keys():
        print("") if (verbose > 2) else None # Debug.
        print("Generating Psi4 psithon input for conformer {}".format(key)) if (verbose > 1) else None # Info.

        if (verbose > 2): # Debug.
            print("{}".format(72*"."))
            pp.pprint(conformers[key])

        conformerToPsithon(conformers[key], len(conformers), memStr, psi4, verbose)

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def printFooter(keywords):
    """
    Prints exit message and additional execution details.

    Parameters
    ----------
    keywords : dict
        Set of runtime variables 

    Examples
    --------
    >>> printFooter(keywords)
    """

    print("") if (keywords['Verbosity'] > 1) else None # Info.

    # Warnings.
    if (keywords['Verbosity'] > 0):
        print("{}------------------------------------".format(18*" "))
        print("{} Execution terminated successfully! ".format(18*" "))
        print("{}------------------------------------".format(18*" "))
        print("\n{}\n".format(36*"-="))
        print("{}{}".format(21*" ", "Thank you for using Smiles4Psi"))
        print("\n{}\n".format(36*"-="))

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def main(argv):
    """
    This program automates the generation quantum chemical calculation input
    files for Psi4. The input files contain geometries extracted from a
    genetic algorithm conformational search executed with Open Babel.

    In a nutshell, what it does is:

    (1) Starting from the SMILES description of a molecule it uses Open
    Babel to construct the molecule's 3D-structure.

    (2) From the 3D-structure, it uses Open Babel to run a genetic algorithm
    conformational search and obtain a given number of geometries for the
    molecule in a trajectory-like file.

    (3) It then extracts the geometry of each conformer from the trajectory-
    like file and stores it in a python dictionary object along with
    complementary information about the molecule.

    (4) Using the geometries and additional finormation of each conformer it
    generates Psi4 Psithon input files with a given QM setup.

    Parameters
    ----------
    argv : list
        Arguments passed after the script call during terminal execution.

    Returns
    -------
    molFile : file
        File with a molecular 3D structure.
    confOutTrajFile : file
        Trajectory-like file with conformers from the GA search.
    p4InFolder/p4InFiles : directory/files
        Psithon files for every conformer contained in a new directory.

    Examples
    --------
    >>> main(sys.argv[1:])
    """

    # Check arguments given:
    print("No arguments passed. Check documentation: `python Smiles4Psi --help`") if (not argv) else ""
    
    # Call parser creator to process options introduced by the user at execution.
    parser = createParser()
    args = parser.parse_args()

    # Process parsed arguments into usable objects.
    keywords = processArgs(args)

    # Print the program's header.
    printHeader(keywords)

    # Using Open Babel, pass the SMILES string to write a file with the molecular structure.
    smilesToMolFile(keywords['SMILES'],
                    keywords['3DStructFileName'],
                    keywords['3DStructFileFormat'],
                    keywords['Verbosity'])

    # Use the molecular structure file to run a Genetic Algorithm Conformational search using Open Babel.
    molFileToConfTraj(keywords['SMILES'],
                      keywords['3DStructFileName'],
                      keywords['3DStructFileFormat'],
                      keywords['ConfTrajFileName'],
                      keywords['NumberOfConformers'],
                      keywords['ScoringFunction'],
                      keywords['Verbosity'])

    # From the resulting trajectory-like file, extract all the conformers information to a dictionary.
    conformers = splitTrajToDict(keywords['SMILES'],
                                 keywords['ConfTrajFileName'],
                                 keywords['MoleculeNameAbbrev'],
                                 keywords['MoleculeNameFull'],
                                 keywords['Verbosity'])

    # Build and write Psithon inputs from the conformers dictionary.
    buildPsithonsFromConfDict(conformers, 
                              keywords['Psi4'],
                              keywords['Verbosity'])

    # Print the program's footer.
    printFooter(keywords)

    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

# Main code execution
if __name__ == "__main__":
    main(sys.argv[1:])
