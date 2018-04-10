#!/usr/local/bin/python3.4
# Copyright (C) Tyson Shepherd, Sakul Ratanalert, Adam Tao, LCBB, MIT
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
import sqlite3
occ=[]
occRev=[]
regi=[]
fIn=open("blastOut.txt", "r")
for line in fIn:
	occ.append(line.split()[0]);
fIn.close()
fIn=open("blastOutRev.txt", "r")
for line in fIn:
	occRev.append(line.split()[0]);
fIn.close()
fIn=open("register.txt", "r")
for line in fIn:
	regi.append(line.rstrip())
fIn.close()
#fOut=open("currate.txt", "w")
conn=sqlite3.connect('aPrimerCustomDB.db')
c=conn.cursor()
for instan in regi:
	if occ.count(instan)>1:
		tmpTxt="DELETE FROM primers WHERE uid='"+instan+"'"
		c.execute(tmpTxt);
	elif occRev.count(instan)>1:
		tmpTxt="DELETE FROM primers WHERE uid='"+instan+"'"
		c.execute(tmpTxt);
conn.commit()
c.close()
