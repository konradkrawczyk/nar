Natural Antibody Reference.

A. EXPLANATION.

These files serve to perform hydrophobic annotation of surfaces and compare them to the scores obtained from a comprehensive model of human antibody repertoire and therapeutic antibodies.

The two scores included are:

1. Linear score- in human terms it is the measure of how big is the hydrophobic surface when all its constituents are counted.

In maths terms it is:

SUM[for each residue R in CDR+vicinity] ASA(R)*HYDROPHOBICITY(R,S)

Where ASA(R) is the surface accessible area (absolute terms).

Where HYDROPHOBICITY(R,S) is the normalized hydrophobicity score according to scheme S (of which there are several:

#0: Gravy - Kyte & Doolttle (1982)
#1: Wimley & White (1996)
#2: Hessa et al. (2005)
#3: Eisenberg & McLachlan (1986)
#4: Black & Mould (1990)

NB: CDR vicinity is currently defined as 4.0 away from the CDRs proper (Chothia CDRs)

2. Patch score - in human terms it is the measure of how big the hydrophobic surface is, allowing the proximal residues to reinforce the score (it favors hydrophobic residues clumped up rather than scattered on the surface).

In mathematical terms:

Sum Coulombic(R1,R2) | Residue R2 in CDR+vicinity, Residue R2 in CDR+vicinity, R1!=R2, where d(R1,R2)< 7.5A

where Coulombic(R1,R2) = [HYDROPHOBICITY(R1,S)*HYDROPHOBICITY(R2,S)]/(d*d)

where d(R1,R2) is the distance in angstroms of two closest heavy atoms of residues R1, R2.

Coulombic has NOTHING to do with charge -- it is just a shorthand to understand whats going on:

if the two residues are both high hydrophobicity score -- this is rewarded.
If the two residues are far apart: this is penalized.

B. USAGE

The pipeline is to firstly annotate a structure and then to 

It is advisable to feed it Chothia numbered sequence first!

1. Annotate:

General syntax: python Annotator.py [hydrpohobic annotation(see above)] [input_file] [ab chains] [output_file]
Specific example: python Annotator.py 0 examples/cristian.pdb IG results.txt

2. Display results.

General Syntax: python PlotResults.py [results_file] [linear/patch] [hydrophobicity parameter]
Specific example: python PlotResults.py results.txt linear 0

Plotting will display your query antibody on the backdrop of the natural antibodies and therapeutic ones.

C. TODOs

Allow for other schemes other than Chothia, Chothia CDRs.
Comment better.
Deal with Surface accessibility software failing every so often.
Add ASA to the patch score?
R-include the PROPKA charge information in the calculations -- currently they can just be confusing.
