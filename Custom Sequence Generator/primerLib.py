import string
import math
import Levenshtein

def get_melting_temp(enthalpy, entropy, selfcomp, gccont, plen):
	R=1.9872 #cal/mol/K
#	R=3.327
	molarConc=0.0000002 #mol/L oligo
	nacl=0.05 #mol/L NaCl
	cfac=((4.29*gccont-3.95)*math.log(nacl)+0.94*(math.log(nacl))**2)*10**-5 # monovalent concentration factor
	if (selfcomp): # if self-complementary
		x = 1.0;
	else:
		x = 4.0;
	TmB = enthalpy * 1000 / (entropy + R*math.log(molarConc / x)) - 273.15;
	Tm=1/((1/TmB)+cfac);
	Tm=(Tm+273.15)*0.982 - 273.15; # SR
#	Tm=Tm*0.96
#	Tm=Tm*0.947
	return Tm;

def seqCheck(templ, prim): # templ = template sequence, prim = primer to check
	j=0; # initialize counter
	lp=len(prim); # length of primer
	txt3=prim[::-1] # reverse strand
	temLen=int(len(templ)) # length of template is half of input because it is doubled for rollover
	txt4=''; # initialize reverse complementary strand
	for i in range(0, len(txt3)):
		if txt3[i]=='G':
			txt4+='C';
		elif txt3[i]=='C':
			txt4+='G';
		elif txt3[i]=='A':
			txt4+='T';
		elif txt3[i]=='T':
			txt4+='A';
	txt3='' # clear
	txt3=txt4; # set to be reverse complementary strand
	txt4='' # clear
	for i in range(0, len(templ)-lp-2):
		if Levenshtein.ratio(templ[i:i+lp+2], prim) > 0.83:
			j=j+1
	for i in range(0, len(templ)-lp-2):
		if Levenshtein.ratio(templ[i:i+lp+2], txt3) > 0.83:
			j=j+1
	if j>14: # what is this?
		return 0;
	else:
		return 1;

def get_enthalpy(seq): # seq = sequence
	nnenthalpy=0.0; # initialize nearest-neighbor enthalpy, kcal/mol
	i=0 # initialize counter for ith nucleotide
	enthalpy = 0.2; # initialization constant for enthalpy, kcal/mol, from INSERT REF HERE
	for m in seq: # for each nt in the sequence
		i=i+1 # get next nucleotide
		if ((i==1 or i==len(seq)) and (m=='A' or m=='T')):
			enthalpy=enthalpy+2.2 # From INSERT REF HERE, penalty for starting/ending A/T
		if (i<len(seq)): # sliding window to get nearest neighbors
			n=seq[i] # next nucleotide is nearest neighbor
			if (m == 'A'): # From INSERT REF HERE, nearest neighbor terms
				if (n == 'A'):
					nnenthalpy = -7.6; # kcal/mmol
				elif (n == 'T'):
					nnenthalpy = -7.2;
				elif (n == 'C'):
					nnenthalpy = -8.4;
				elif (n == 'G'):
					nnenthalpy = -7.8;
			elif (m == 'T'):
				if (n == 'A'):
					nnenthalpy = -7.2;
				elif (n == 'T'):
					nnenthalpy = -7.6;
				elif (n == 'C'):
					nnenthalpy = -8.2;
				elif (n == 'G'):
					nnenthalpy = -8.5;
			elif (m == 'C'):
				if (n == 'A'):
					nnenthalpy = -8.5;
				elif (n == 'T'):
					nnenthalpy = -7.8;
				elif (n == 'C'):
					nnenthalpy = -8.0;
				elif (n == 'G'):
					nnenthalpy = -10.6;
			elif (m == 'G'):
				if (n == 'A'):
					nnenthalpy = -8.2;
				elif (n == 'T'):
					nnenthalpy = -8.4;
				elif (n == 'C'):
					nnenthalpy = -9.8;
				elif (n == 'G'):
					nnenthalpy = -8.0;
		else:
			nnenthalpy=0.0
		enthalpy = enthalpy + nnenthalpy # add nearest-neighbor enthalpies to initialization, kcal/mol
	return enthalpy; # kcal/mol

def get_entropy(seq): # seq = sequence
	R=1.9872 #cal/mol/K
	molarConc=0.0004 #mol/L oligo
	nacl=0.05 #mol/L NaCl
#	cfac=((4.29*GCcontent(seq)-3.95)*math.log(nacl)+0.94*(math.log(nacl))**2)*10**-5
#	mg=0.002
#	aa=0.0000392
#	ba=-0.00000911
#	ca=0.0000626
#	da=0.0000142
#	ea=-0.000482
#	fa=0.000525
#	ga=0.0000831
#	mgfac=aa+ba*math.log(mg)+GCcontent(seq)*(ca+da*math.log(mg))+(ea+fa*math.log(mg)+ga*(math.log(mg))**2)/(2*(len(seq)-1))
	nnentropy=0.0; # initialize nearest-neighbor entropy, cal/mol/K, from http://www.math.utah.edu/~palais/pcr/papers/SL2004.pdf
	i=0; # initialize counter for ith nucleotide
	entropy = -5.7; # initialization constant for entropy, cal/mol/K, from INSERT REF HERE
	for m in seq: # for each nt in the sequence
		i=i+1 # get next nucleotide
		if ((i==1 or i==len(seq)) and (m=='A' or m=='T')):
			entropy=entropy+6.9; # pentaly for starting/ending A/T
		if (i<len(seq)): # sliding window to get nearest neighbors
			n=seq[i] # next nucleotide is nearest neighbor
			if (m == 'A'):
				if (n == 'A'): # NN parameters from INSERT REF HERE
					nnentropy = -21.3; #cal/mol/K
				elif (n == 'T'):
					nnentropy = -20.4;
				elif (n == 'C'):
					nnentropy = -22.4;
				elif (n == 'G'):
					nnentropy = -21.0;
			elif (m == 'T'):
				if (n == 'A'):
					nnentropy = -21.3;
				elif (n == 'T'):
					nnentropy = -21.3;
				elif (n == 'C'):
					nnentropy = -22.2;
				elif (n == 'G'):
					nnentropy = -22.7;
			elif (m == 'C'):
				if (n == 'A'):
					nnentropy = -22.7;
				elif (n == 'T'):
					nnentropy = -21.0;
				elif (n == 'C'):
					nnentropy = -19.9;
				elif (n == 'G'):
					nnentropy = -27.2;
			elif (m == 'G'):
				if (n == 'A'):
					nnentropy = -22.2;
				elif (n == 'T'):
					nnentropy = -22.4;
				elif (n == 'C'):
					nnentropy = -24.4;
				elif (n == 'G'):
					nnentropy = -19.9;
		else:
			nnentropy=0.0
		entropy = entropy + nnentropy # add nearest-neighbor entropies to initialization entropy, cal/mol/K
#	entropy=entropy+get_enthalpy(seq)*mgfac*1000
	nacl=nacl+(math.sqrt(2.0-0.2))/1000; # https://www.ncbi.nlm.nih.gov/pubmed/11673362
#	entropy=entropy+0.368*len(seq)*math.log(nacl)/2.0 http://www.math.utah.edu/~palais/pcr/papers/SL2004.pdf
	entropy=entropy+0.368*len(seq)*math.log(nacl)/2.0 # add NaCl correction
	return entropy; #cal/mol/K

def check_self_comp(seq): # check if self-complementary
	checker = 0; # initialize as false
	seqRC=[] # initialize reverse complementary sequence
	seqRev=list(reversed(seq)); # reverse sequence
	seqOri=list(reversed(seqRev)); # original sequence
	for i in seqRev: # generate reverse complementary sequence
		if i=='A':
			seqRC+='T'
		elif i=='T':
			seqRC+='A'
		elif i=='C':
			seqRC+='G'
		elif i=='G':
			seqRC+='C'
	if seqOri==seqRC: # if original is same as reverse complement
		checker=1; # set to true
	return checker;

def GCcontent(seq):
	content = (seq.count('G')+seq.count('C'))/len(seq); # the total number of G and C divided by total length
	return content;

