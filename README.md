# UncertainButterflyCounting

AUTHOR: Alexander Zhou

INSTITUTION: Hong Kong University of Science and Technology

E-MAIL: atzhou@cse.ust.hk

Please E-mail me with any queries about the code or the subsequent paper but given when the mail is sent there is a non-zero (very likely) chance I have already forgotten everything about the code.

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
	*	13 : Multiple UBS-Vertex Runs
	*	14 : Multiple UBS-Edge Runs
	*	15 : Multiple PES-Vertex Runs
	*	16 : Multiple PES-Edge Runs

argv[2] = Threshold Probability

argv[3] = ID File

argv[4] = Edge File

argv[5] = # of samples (optional)

	*	Required for processes 13, 14, 15 and 16
	*	If no argument, defaults to 1000

argv[6] = Output File (optional)
	
	*	Required for processes 13, 14, 15 and 16
	*	Required for processes 