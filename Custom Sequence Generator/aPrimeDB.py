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
import sqlite3

conn = sqlite3.connect('aPrimerCustomDB.db');
c = conn.cursor()

# Create table
c.execute('''CREATE TABLE primers
             (uid text, length integer, forPrime text, revPrime text, forTm integer, revTm integer, template text, prodSeq text, closeValid text)''')

# Insert primers
conn.commit()
conn.close()
