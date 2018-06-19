1. USAGE
--------------
To run the Software Package from the command line using default parameters,
Run the run_predict_binding.py script as:
```shell
./run_predict_binding.py XXXX.pdb
```


2. INPUT
--------------
The three dimensional coordinates of the atomic positions of the protein
complex should be given in a PDB format.
Place the PDB files in the 'input' directory of the package.


3. OUTPUT
--------------
After a successful run of the program, the corresponding output files are
generated in the 'output' directory.
For each complex, two text files will be generated -
1. XXXX_scores.out - Contains the prediction and other metrics like z-score
2. XXXX_respairs.out - Contains the residue pairs comprising the interface

In case the Alanine Scanning flag is provided, an additional file will be generated,
3. XXXX_alascan.out - Contains the residues which upon mutation might destabilize the interaction


4. EXAMPLES
--------------
Examples for the input and output files are provided in the 'examples' directory


5. ADDITIONAL PARAMETERS
--------------
Optional flags that can be turned on -
USAGE: ./run_predict_binding.py XXXX.pdb -d 8.0 -t mm

1. Distance threshold for interaction definition: -d, --cutoff, choices = {4.0, 6.0, 8.0}
	NOTE: Default value = 4.0

2. Types of interacting atoms between residues: -t, --intertype, 
    choices = {'all', 'mm', 'ms', 'ss'}
	NOTE: Default value = ss
	....all - Interactions between all atoms
	....mm  - Interactions between main chain atoms
	....ms  - Interactions between main chain atoms and side chain atoms
	....ss  - Interactions between side chain atoms

3. To specify your own 20x20 matrix for residue-pair scores: -cp, --custom_pot, file
	NOTE: The custom potential should be in the format of a .csv file with the residue_pair in
		  first column and the corresponding score in the second column.
		  Example file provided in the 'examples' directory (custom_pot.csv)

4. To perform Alanine Scanning with the default potential: -alscn, --alascan, 
   choices = {'0', 1', '2', '3'}
	NOTE: Default value = '0'
		  0 - Do not perform Alanine Scanning
		  1 - Mutate all the residues on the interface to Alanine
		  2 - Mutate all the residues on the interface to every other residue
		  3 - Mutate a given list of residues to every other residue
		  (residues to be provided in a separate file)

	*NOTE: In case mode '3' is selected, a separate file containing the residues
		   to be mutated should be provided in the 'input' directory. An example
		   file is given in the 'examples' directory (4pza_alscan_input.txt).

5. In case of multimers, specify particular interfaces, -p1, -p2, --protein_1, --protein_2

6. To specify a name for the output files, -o, --outfile

This help can be accessed on the command line using:
./run_predict_binding.py -h or ./run_predict_binding.py --help


6. SOURCE CODE
----------------
All the scripts for the package can be found in the 'scripts' directory.


7. LICENSE
----------------
Please refer to COPYING.txt for the full license or COPYING_LESSER.txt for a condensed version.
