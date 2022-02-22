import xmlrpc.client
#
# connect to localhost, port 4242 running the RPC server
#
proxy = xmlrpc.client.ServerProxy("http://localhost:4242/")
multicall = xmlrpc.client.MultiCall(proxy)
#
# calc: function on the RPC server to route the origami
# calc(projectName, helicalForm, helicalTurns, plyFile, sequenceFile):
#	projectName: short descriptive name, name of the generated folder
#	helicalForm: 'Bform', 'Aform', 'Hybrid', 'Twisted'
#	helicalTurns: multiplied by 11 (Aforms) or 10.5 (Bform)
#	plyFile: Path to ply file
#	sequenceFile: Path to sequence file
#
multicall.calc("TestTet66","Aform",6,"./tet.ply","./M13.txt");
result = multicall()
multicall.calc("TestOct66","Aform",6,"./oct.ply","./M13.txt");
result = multicall()
multicall.calc("TestOct44","Aform",4,"./oct.ply","./M13.txt");
result = multicall()
multicall.calc("TestPBip66","Aform",6,"./pbip.ply","./M13.txt");
result = multicall()
