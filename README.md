# Uncertain Butterfly Counting

AUTHOR: Alexander Zhou

INSTITUTION: Hong Kong University of Science and Technology

E-MAIL: atzhou@cse.ust.hk

Please E-mail me with any queries about the code or the subsequent paper.

We provide the IMDB dataset, which we have assigned each edge an existential probabilities, in the file format used for our implementation. Further information regarding the IMDB dataset may be found at https://www.cise.ufl.edu/research/sparse/matrices/Pajek/IMDB.html

# Arguments

argv[0] = ./iubfc

argv[1] = Process Type

	*	0 : UBFC (Baseline)
	*	1 : IUBFC (Improved) - VP
	*	2 : IUBFC (Improved) - EP
	*	3 : UBS-Vertex
	*	4 : UBS-Edge
	*	5 : PES-Vertex
	*	6 : PES-Edge
	*	13 : 500 UBS-Vertex Runs
	*	14 : 500 UBS-Edge Runs
	*	15 : 500 PES-Vertex Runs
	*	16 : 500 PES-Edge Runs

argv[2] = Threshold Probability

	*	Between 0.0 and 1.0
	
argv[3] = ID File

	*	Containing all IDs, non-repeating
	
argv[4] = Edge File

	*	tsv formal (ID1	ID2 Probability)

argv[5] = # of samples (optional)

	*	Required for processes 13, 14, 15 and 16
	*	If no argument, defaults to 1000

argv[6] = Output File (optional)
	
	*	Required for processes 13, 14, 15 and 16
	
Examples:

Run Baseline with threshold 0.7

	*	./iubfc 0 0.7 dataID.txt dataEdge.txt

Run UBS-Vertex 500 times at threshold 0.5, each batch containing 10,000 samples
	
	*	./iubfc 13 0.5 dataID.txt dataEdge.txt 10000 dataOut.txt
