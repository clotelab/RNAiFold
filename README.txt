RNAiFold3.1 - Clotelab - Boston College

/******************************************************************************
 *   Copyright (C) 2014  Juan Antonio Garcia Martin, Peter Clote ,Ivan Dotu   *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/

COMPILING RNAiFold 3.1:

Requirements:
	- Software is developed in OR-Tools (https://github.com/google/or-tools), which is required for compiling RNAiFold.
	- A version of compiled Vienna Package libraries for Linux-amd64 is included. Download and compile the latest vienna libraries for your SO and copy libRNA.a file (from lib directory) into RNAiFold ViennaRNA directory. (This version is tested for ViennaRNA 2.1.9)
	- A version of compiled RNAstructure libraries for Linux-amd64 and OSX is included. Download (http://rna.urmc.rochester.edu/RNAstructure.html) Compile the latest RNAstructure libraries for your SO and libHybridRNA.so and libRNA.so to RNAiFold RNAstructure directory. (This version is tested for RNAstructure 5.7)
Compilation step:
1. Download (https://github.com/google/or-tools) and install OR-Tools following the instructions in the installation guide included.
2. The expected default directory for OR-Tools installation is $(HOME)/OR-Tools-src/or-tools-read-only, if you install OR-Tools in a different directory change the value of ORTOOLS_DIR into Makefile 
3. Assign the variable ORTOOLS_SRC in Makefile to the directory of your current OR-Tools installation
4. make

One single executable will be generated:
	RNAiFold

USAGE: 

Only a target secondary structure is required for running RNAiFold, either as an input parameter or into the input file

Only a target secondary structure is required, which can be provided either as a parameter by command line flags or inside an input file.
   ./RNAiFold -RNAscdstr <SECONDARY_STRUCTURE>
   ./RNAiFold -InputFile <INPUT_FILE>

INPUT parameters:
   Design constraints can be provided as command line flags or inside the input file using the appropriate label preceeded by the "pound" symbol ("#") and writing the value in the next line. 
   This are the valid labels/flags organized by sections.

   Input file.
     -InputFile: (string) If provided, any parameters will be read from the given file. As explained parameters must be preceeded by the pound ('#') symbol
         Input File format:
           Input file must contain a valid secondary structure, all the other fields are optional, RNAiFold input file format is:
             > Fasta comment
             Target structure
             Sequence constraints
             # Parameter
             Parameter value

   Target structure(s) and folding temperature(s)
     -RNAscdstr: (string) Input RNA secondary structure(s). To indicate multiple target structures they must be specified in the same line separated by the pipe '|' symbol.
     -temp: (string) Target folding temperature(s) in celsius. DEFAULT:37 To indicate multiple folding temperatures for the corresponding number of target structures write them in the same line separated by commas.

   Number of solutions and search time limit:
     -TimeLimit: (int64) Search time limit in seconds. DEFAULT: 600 (10  minutes)
     -MAXsol: (int32) Input MAX solutions. DEFAULT:8

   Sequence and structure compatibility:
     -RNAseqcon: (string) Input RNA sequence constraints. DEFAULT:
     -RNAcompstr: (string) Compatible RNA secondary structure. DEFAULT:
     -IncompBP: (string) List of incompatible base pairs. DEFAULT:

   Thermodynamic parameters:
     -RNAlibrary: (string) RNA library for folding and computing energy values: Allowed models are: Vienna (ViennaRNA package) RNAstructure (Mathews' Lab). Default RNA library is ViennaRNA (Vienna). DEFAULT:Vienna
     -EnergyModel: (string) Energy model for ViennaRNa library. Allowed models are: 2004 (Turner '04), 1999 (Turner '99), 2007 (Andronescu '07) . Default energy model is Turner '04 (2004). DEFAULT:2004
     -dangles: (int32) Dangling treatment. DEFAULT:2

   Objective functions:
     -MFEstructure: (int32) Sequence MFE structure(s) must be the target structure(s). DEFAULT:1
     -MinimizeEnsDef: (int32) Minimize ensemble defect for the target structure: print only sequences whose ensemble defect for the target structure is equal or lower than the previous one. DEFAULT:0
     -MinimizeMFE: (int32) Minimize free energy of the MFE structure: print only sequences whose free energy when folded into their MFE structure is equal or lower than the previous one. DEFAULT:0

   Local constraints:
     -HelixCstrs: (string) Helix local constraints. DEFAULT:
     -LocalCstrs: (string) Local constraints. DEFAULT:

   Amino acid constraints:
     -AAseqcon: (string) List of amino acid sequence constraints. DEFAULT:
     -AAsimilCstr: (int32) Blosum similarity threshold (-4 to 4)(default 5 forces aminoacids to be equal to those specified) or allow similar amino acids based on a specific classification (6- Classification 1, 7-Classification 2 (see manual). DEFAULT:5
     -AAstartPos: (string) List of starting positions for amino acid constraints. DEFAULT:
     -AAtarget: (string) List of target amino acid sequence for blosum similarity score maximization. DEFAULT:

   Large neighborhood search:
     -LNS: (int32) Activate Large Neighborhood Search. DEFAULT:0
     -LNSrestartTime: (int32) Maximum consecutive restarts in LNS without changes on fixed positions. DEFAULT:15
     -LNStimeMultiplier: (int32) Multiplier for restart time in milliseconds (Time=length*multiplier). DEFAULT:0
     -LNSunchangedRestarts: (int32) Maximum consecutive restarts in LNS without changes on fixed positions. DEFAULT:10
     -MaxBlosumScore: (int32) Maximize blosum score of target sequences (default disabled). DEFAULT:0

   Limit free energy of the structure:
     -MinMFE: (double) Minimum free energy allowed for the sequence folded into its minimum free energy structure. DEFAULT:1000
     -MaxMFE: (double) Maximum free energy allowed for the sequence folded into its minimum free energy structure. DEFAULT:1000


   Output parameters:
     -ShowMeasures: (int32) Show structural diversity measures. DEFAULT:1
     -ShowHelices: (int32) Show helix identifiers for local constraints. This flag is only to visualize heclices, no search will be performed.


   Nucleotide composition:
     -minAU: (int32) Minimum number of AU base pairs. DEFAULT:0
     -maxAU: (int32) Maximum number of AU base pairs. DEFAULT:-1
     -minGC: (int32) Minimum  number of GC base pairs. DEFAULT:0
     -maxGC: (int32) Maximum number of GC base pairs. DEFAULT:-1
     -minGU: (int32) Minimum  number of GU base pairs. DEFAULT:0
     -maxGU: (int32) Maximum number of GU base pairs. DEFAULT:-1

     -minA: (string) List of minimum number of As in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -maxA: (string) List of maximum number of As in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -minC: (string) List of minimum number of Cs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -maxC: (string) List of maximum number of Cs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -minG: (string) List of minimum number of Gs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -maxG: (string) List of maximum number of Gs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -minU: (string) List of minimum number of Us in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -maxU: (string) List of maximum number of Us in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:

     -consA: (string) List of maximum number of consecutive As in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -consC: (string) List of maximum number of consecutive Cs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -consG: (string) List of maximum number of consecutive Gs in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:
     -consU: (string) List of maximum number of consecutive Us in the full sequence (N) or in a specific range (N StartPos EndPos). DEFAULT:

     -minGCcont: (double) Minimum GC content. DEFAULT:0
     -maxGCcont: (double) Maximum GC content. DEFAULT:100

   Search parameters and heuristics:
     -RandomAssignment: (int32) Activate random value heuristic. DEFAULT:0
     -BPthreshold: (int32) Probablility of selecting the next BP assignment in value heuristic. DEFAULT:85
     -UPthreshold: (int32) Probablility of selecting the next UP assignment in value heuristic. DEFAULT:100
     -VarHeuristic: (int32) Variable heuristic 1-Helices bottom to top 2-In to out . DEFAULT:2
     -HelixHeuristic: (int32) Helix ordering heuristic for the search 1-Simple overlap 2-Base pair overlap 3-Total overlap (default 2). DEFAULT:2
     -IncludeDangles: (int32) Include dangling positions when creating helices. DEFAULT:1

OUTPUT:
  Three possible types of results can be returned:
    - Solution found: For each solution found the following information is displayed.
        Sequence.
        GC content and the number of base pairs of each type (strong, weak and wobble).
        Amino acid sequence and blosum score if blosum maximization is active.
        Free energy of the structure in kcal/mol
        Additional measures
    - No solution found: If search time is reached and no solution has been found within this time limit.
    - No possible solution: If the target structure (with specified constraints) has no solution and the time limit has not been reached.

EXAMPLES: 
	The following command line example illuestrates the use of RNAiFold with input flags (no input file) corresponding for the parameters for tRNA.fas input file:
    ./RNAiFold -RNAscdstr '(((((.(..((((.........)))).(((((.......))))).....(((((.......)))))).))))).' -RNAseqcon 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCUGGCUCG' -MAXsol 5 -dangles 2 -minGC 10 -maxGC 20 -LNS 1

    Several sample input files are included, each one illustrating a different design usign the features included in RNAiFold.
    
	File: tRNA.fas
	Command: ./RNAiFold -InputFile examples/tRNA.fas
	Description: tRNA structure with sequence constraints, solutions must have a minimum of 10 and a maximum of 20 GC base pairs. Large neighborhood search strategy is applied in the search.
	

	File: cofold.fas
	Command: ./RNAiFold -InputFile examples/cofold.fas
	Description: Hybridization, design two RNA chains whose hybidized MFE structure is the target
	

	File: partialTarget.fas
	Command: ./RNAiFold -InputFile examples/partialTarget.fas
	Description: Partial target, where some positions in are not constrained to be paired or unpaired, include sequence constraints using IUPAC codes and incompatible base pair connstraints


	File: Overlapping-HIV1-FSS.fas.fas
	Command: ./RNAiFold -InputFile examples/Overlapping-HIV1-FSS.fas.fas
	Description: This example corresponds to fremeshift stimulating signal in HIV-1 the overlapping region of Gag and Gag-Pol polyproteins
	             Generates sequences that code for amina acid sequences with blosum similarity >= 2 to the given amino acid sequences in two different coding frames (0 and +1). 
	             

	File: SECIS.fas.fas
	Command: ./RNAiFold -InputFile examples/SECIS.fas.fas
	Description: Design of selenocystein codon insertion in AJF73661.1 gene.
	             Generate a sequence that optimizes BLOSUM62 similarity to the given amino acid sequence from AJF73661.1 gene containing a selenocystein codon at position 1.

	File: Lambda_thermo.fas
	Command: ./RNAiFold -InputFile examples/Lambda_thermo.fas
	Description: Example of design for multiple target structure and folding tempertures. MFE structure of the solutions sequences must be target structure 1 at 32ºC and target structure 2 at 62ºC.
	             Large neighborhood search strategy is applied in the search.
	
	File: AAoptimization.fas
	Command: ./RNAiFold -InputFile examples/AAoptimization.fas
	Description: Amino acid optimization example. Solution sequences fold into the trget MFE structure and optimize blosum similarity to IFSSLPGLVPKGKE amino acid sequence.

Additional help and manual:
    For additional help please visit http://bioinformatics.bc.edu/clotelab/RNAiFold/index.php?tab=manual


