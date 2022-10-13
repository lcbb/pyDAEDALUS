import xmlrpc.client
import os
#
# connect to localhost, port 4242 running the RPC server
#
proxy = xmlrpc.client.ServerProxy("http://localhost:4242/")
multicall = xmlrpc.client.MultiCall(proxy)
#
# calc: function on the RPC server to route the origami
# calc(projectName, helicalForm, helicalTurns, plyFile, sequenceFile):
#	projectName: short descriptive name, name of the generated folder
#                   (write "M13.txt" to use default sequence)
#	helicalForm: 'Bform', 'Aform', 'Hybrid', 'Twisted'
#	helicalTurns: multiplied by 11 (Aforms) or 10.5 (Bform) 
#			    to get minimum edge length
#	plyFile: Path to ply file
#	sequenceFile: Path to sequence file
#

ply_dir = "C:/pyDAEDALUSX/PLY_Files" # path to directory containing input geometries

for geom in os.listdir(ply_dir):
    ply = os.path.join(ply_dir, geom) # path to specific geometry file
	name = geom # project name
    multicall.calc(name,"Aform",4,ply,"M13.txt");

result = multicall()
