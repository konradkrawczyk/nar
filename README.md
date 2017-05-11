# nar


These files serve to perform hydrophobic annotation of surfaces.

The two scores included are:

1. Linear score- in human terms it is the measure of how big is the hydrophobic surface when all its constituents are counted.

In maths terms it is:

SUM[for each residue R in CDR] ASA(R)*HYDROPHOBICITY(R,S)

Where ASA(R) is the surface accessible area (absolute terms).

Where HYDROPHOBICITY(R,S) is the normalized hydrophobicity score according to scheme S (of which there are several:

#0: Gravy - Kyte & Doolttle (1982)
#1: Wimley & White (1996)
#2: Hessa et al. (2005)
#3: Eisenberg & McLachlan (1986)
#4: Black & Mould (1990)

2. Patch score - in human terms it is the measure of how big the hydrophobic surface is, allowing the proximal residues to reinforce the score (it favors hydrophobic residues clumped up rather than scattered on the surface).


