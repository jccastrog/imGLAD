#This is an example file on how to run imGLAD

A simple run of imGLAD can be:

   ```bash
   $> ./fitModel.py -t _example/example_targetGenome.fa -sp 'Escherichia coli'
   ```
Followed by:

   ```bash
   $> ./probEstimate.py -t _example/example_targetGenome.fa -m _example/example_positiveSample.tbl -p _example/example_parameters.txt -l single
   ```
The files listed in [_example](./) are

[example_genomeList.txt](./example_genomeList.txt):	A list with genome accesion numbers can be passed as option -l to fitModel.py

[example_positiveSample.tbl](./example_positiveSample.tbl):	BLAST alignment of a positive metagenomic sample against the target genome

[example_negativeSample.tbl](./example_negativeSample.tbl):	BLAST alignment of a negative metagenomic sample against the target genome

[example_parameters.txt](./example_parameters.txt):	Parameters generated with fitModel.py for an E.coli genome

[example_targetGeneome.fa](./example_targetGeneome.fa):	A reference E. coli genome str. O157:H7
 
