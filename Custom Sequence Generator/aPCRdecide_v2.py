#!/usr/local/bin/python3.4
# Copyright (C) Tyson Shepherd, Sakul Ratanalert, Adam Tao, LCBB, MIT  name of author
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
import math
from primerLib import *
import os, sys
import hashlib
import random
import sqlite3
#
#
#
def primer_tm(temp, length): # inputs: temp = template sequence, length = desired length amplified sequence
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	primers=[] # initialize storage
	primers.append([])
	primers.append([])
	primers.append([])
	primers.append([])
	testOne = ''; # initialize placeholder for first primer
	testTwo = ''; # initialize placeholder for second primer
	count = 0; # hanging variable?
	for p in range(18,23): # for first primer length [18,23)
		for q in range(0, int(len(temp)-length)): # for starting positions of the amplified length from 0 to end of template (template input sequence is duplicated to help with rollover)
			testOne=temp[q:p+q]; # set first primer to test
			Tm1 = get_melting_temp(get_enthalpy(testOne), get_entropy(testOne), check_self_comp(testOne), GCcontent(testOne), len(testOne)); # get Tm of testOne
			if ((Tm1 >= 53.5) and (Tm1 < 57.5)): # if Tm is [54,57]
				if (((testOne[p-1]=='A') or (testOne[p-1]=='T')) and (GCcontent(testOne[:9]) > GCcontent(testOne[len(testOne)-9:]))): # if primer ends in A/T and GC content of first 9 nt is higher than last 9 nt
					if (((GCcontent(testOne) >= 0.4) and (GCcontent(testOne) <= 0.6)) and ((testOne[0]=='C') or (testOne[0]=='G'))): # if primer overall GC content is [0.4,0.6] and primer starts with C/G
						for t in range(18,23): # for second primer length [18,23)
							testTwo=temp[q+length-t:q+length]; # set second primer to test, ending the desired length away
							Tm2 = get_melting_temp(get_enthalpy(testTwo), get_entropy(testTwo), check_self_comp(testTwo), GCcontent(testTwo), len(testTwo)); # get Tm of testTwo
							Tm1_int = round(Tm1);
							Tm2_int = round(Tm2);
							if ((Tm2_int >= Tm1_int+1) and (Tm2_int <= Tm1_int+3)): # if Tm of testTwo is [1,3] higher than testOne as ints
								if ((GCcontent(testTwo) >= 0.4) and (GCcontent(testTwo) <= 0.6)): # if second primer overall GC content is [0.4,0.6]
									if (((testTwo[t - 1] == 'C') or (testTwo[t - 1] == 'G')) and ((testTwo[0] == 'C') or (testTwo[0] == 'G'))): # if testTwo starts and ends with C/G
										if not (testOne in primers[0]): # store the primers if second primer has not been stored already (what about first?)
											primers[0].append(testOne)
											testTwoB = "".join(complement.get(base, base) for base in reversed(testTwo))
											primers[1].append(testTwoB)
											primers[2].append(str(Tm1)[:2])
#											primers[3].append(str(GCcontent(temp[q:q+length])*100)[:4])
#											primers[3].append(str(temp[q:q+length]))
											primers[3].append(str(Tm2)[:2]);
	return primers
#
# Process form to python
#
fInfor=open("./customSeq.txt","r")
fInrev=open("./customSeq_rev.txt", "r")
tmplA=[]
tmpSeq=''
for line in fInfor:
	tmpSeq+=line.rstrip().lstrip()
tmplA.append(tmpSeq);
tmpSeq=''
for line in fInrev:
	tmpSeq+=line.rstrip().lstrip()
tmplA.append(tmpSeq);
fInfor.close();
fInrev.close();
#fInLfor=open("./lambda.txt", "r")
#fInLrev=open("./lambda_rev.txt", "r")
#tmplB=[]
#tmpSeq=''
#for line in fInLfor:
#        tmpSeq+=line.rstrip().lstrip()
#tmplB.append(tmpSeq);
#tmpSeq=''
#for line in fInLrev:
#        tmpSeq+=line.rstrip().lstrip()
#tmplB.append(tmpSeq);
#fInLfor.close();
#fInLrev.close();
fOutTmp=open("priTst.txt", "w")
fOutAno=open("register.txt", "w")
fOutTwo=open("priRev.txt", "w")
conn=sqlite3.connect('aPrimerCustomDB.db')
c=conn.cursor()
# Change range to the length of the DNA fragment you want to check
for i in range(900, 1000):
	targetLen=i;
	for templDNA in tmplA:
		primers=primer_tm(templDNA, targetLen);
		for j in range(len(primers[0])):
			tmpUID=str(hashlib.md5(str(primers[0][j]+primers[1][j]).encode('utf-8')).hexdigest())
			exState="INSERT INTO primers VALUES ('"+tmpUID+"',"+str(i)+",'"+primers[0][j]+"','"+primers[1][j]+"',"+str(primers[2][j])+","+str(primers[3][j])+",'Custom','','')"
#			print(exState);
			c.execute(exState)
			fOutTmp.write(">"+tmpUID+"\n");
			fOutTmp.write(primers[0][j]+"\n\n");
			fOutTwo.write(">"+tmpUID+"\n");
			fOutTwo.write(primers[1][j]+"\n\n");
			fOutAno.write(tmpUID+"\n");
fOutTmp.close();
fOutAno.close();
fOutTwo.close();
conn.commit();
c.close()
