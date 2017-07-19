# Sumarize BASTA output trees 

import sys
import os
import random
import math
import argparse
import re

import numpy as np
#import networkx as nx
import matplotlib.pyplot as plt
from graph_tool.all import *


parser = argparse.ArgumentParser()
parser.add_argument('--inputF',"-i", help='input file containing the sampled trees in BEAST (usually with extension .trees).',default="")
parser.add_argument('--outputF',"-o", help='output file containing the inferred transmissions.',default="")
parser.add_argument('--burnin',"-b", help='percentage of trees to discard (default 20).', type=int, default=20)
parser.add_argument('--minValue',"-m", help='minimum value for which to add edges to the plot (default 0.1).', type=float, default=0.1)
parser.add_argument('--edgeThickness',"-e", help='maximum thickness of edges (default max 10).', type=float, default=10.)
parser.add_argument('--vertexFont',"-f", help='vertex font size (default 15).', type=int, default=15)
parser.add_argument('--vertexColor',"-vC", help='vertex color scale (default \"Greens\").', default="Greens")
parser.add_argument('--edgeColor',"-eC", help='edge color scale (default \"Reds\").', default="Reds")
parser.add_argument('--outputSize',"-s", help='output figure size (default 900).', type=int, default=900)
parser.add_argument('--format',"-fmt", help='format of output plots. (default \"pdf\", but can be any of \"auto\", \"ps\", \"pdf\", \"svg\", and \"png\").', default="pdf")
parser.add_argument('--noPlotDirect', '-nPD', dest='plotDirect', action='store_false')
parser.set_defaults(plotDirect=True)
parser.add_argument('--noPlotIndirect', '-nPI', dest='plotIndirect', action='store_false')
parser.set_defaults(plotIndirect=True)
parser.add_argument('--noRootInfo', '-nRI', dest='plotRoot', action='store_false')
parser.set_defaults(plotRoot=True)

args = parser.parse_args()


if(args.inputF=="" or args.outputF==""):
	print("Error, input and output files must be specified with -i and -o options.")
	exit()


## extract information from metadata
def extractInfo(string):
	index=0
	while string[index]!="[":
		index+=1
	index2=len(string)-1
	while string[index2]!="]":
		index2-=1
	subString=string[index+1:index2]
	listTraits=subString.split(",")
	host="-1"
	numT="-1"
	for i in range(len(listTraits)):
		trait=listTraits[i].split("=")[0]
		if trait=="&host":
			host=listTraits[i].split("=")[1]
		if trait=="numTransmissions":
			numT=listTraits[i].split("=")[1]
	if host=="-1" or numT=="-1":
		print("Traits in tree are not recognised: could not find host or number of transmission events along branch!")
		print(string)
		exit()
	return [host,int(numT)]



## Split tree into 2 subtrees
def splitTree(tree):
	if tree[0]!="(" or tree[-1]!=")":
		print("Tree does not strart or end in parenthesis")
		print(tree)
		exit()
	openCount=0
	index=0
	while tree[index]!="," or openCount!=1:
		if tree[index]=="(" or tree[index]=="[":
			openCount+=1
		elif tree[index]==")" or tree[index]=="]":
			openCount-=1
		index+=1
		if index==len(tree):
			return [tree]
	
	return [tree[1:index],tree[index+1:-1]]
	

## Separate metadata and rest of the tree
def metaData(tree):
	#print tree
	index=len(tree)-1
	while tree[index]!=")":
		index-=1
		if index==0:
			return [tree]
	
	return [tree[0:index+1],tree[index+1:]]
	
## Find new hosts in the tree
def recurFindHosts(tree,hosts):
	mt=metaData(tree)
	if len(mt)>1:
		host=extractInfo(mt[1])[0]
	else:
		host=extractInfo(mt[0])[0]
	if (host!="Unsampled") and (not (host in hosts )):
		hosts.append(host)
		#print host+" added to host list"
	if len(mt)>1:
		splitT=splitTree(mt[0])
		recurFindHosts(splitT[0],hosts)
		recurFindHosts(splitT[1],hosts)


## perform one recursion on one subtree
def handleTree(mt,root,parentHost,numTransParent,directTrans,indirectTrans,metAlready,metAlreadyInd,origins):
	#print mt
	#print origins
	#print root
	#print parentHost
	#print numTransParent
	if root[1]==0: #no change
		#print origins
		#print "\n\n"
		if len(mt)>1:
			recurTransm(mt[0],parentHost,numTransParent,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
			
	elif root[0]=="Unsampled": #going to an unsampled node
		#print origins
		#print "\n\n"
		if len(mt)>1:
			recurTransm(mt[0],parentHost,numTransParent+root[1],directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
				
	elif parentHost=="Unsampled": #coming from an unsampled node
		if root[0]!="Unsampled":
			if not (root[0] in origins.keys()):
				origins[root[0]]="Unsampled"
		#print origins
		#print "\n\n"
		if len(mt)>1:
			recurTransm(mt[0],root[0],0,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
		
	elif parentHost!=root[0]: #from one host into a different one, both sampled
		if (root[1]+numTransParent)==1: #direct transmission
			if ((not parentHost in metAlready.keys()) or (not root[0] in metAlready[parentHost])) and ((not parentHost in metAlreadyInd.keys()) or (not root[0] in metAlreadyInd[parentHost])):
				if parentHost in metAlready.keys():
					metAlready[parentHost].append(root[0])
				else:
					metAlready[parentHost]=[root[0]]
				directTrans[parentHost][root[0]]+=1
			if root[0] in origins.keys():
				if origins[root[0]]!=parentHost and origins[root[0]]!="Unsampled":
					origins[root[0]]="doubleOrigin"
				else:
					origins[root[0]]=parentHost
			else:
				origins[root[0]]=parentHost
		elif (root[1]+numTransParent)>1: #indirect transmission
			if ((not parentHost in metAlready.keys()) or (not root[0] in metAlready[parentHost])) and ((not parentHost in metAlreadyInd.keys()) or (not root[0] in metAlreadyInd[parentHost])):
				if parentHost in metAlreadyInd.keys():
					metAlreadyInd[parentHost].append(root[0])
				else:
					metAlreadyInd[parentHost]=[root[0]]
				indirectTrans[parentHost][root[0]].append(root[1]+numTransParent)
			if not (root[0] in origins.keys()):
				origins[root[0]]="Unsampled"
		else:
			print("there is a problem, this should not be 0")
			print(root[1]+numTransParent)
			exit()
		#print origins
		#print "\n\n"
		if len(mt)>1:
			recurTransm(mt[0],root[0],0,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
	else: #from one host to itself
		if root[1]+numTransParent==1: #direct transmission?
			print("there is a problem, this should not be 1")
			print(root[1]+numTransParent)
			exit()
		if ((not parentHost in metAlready.keys()) or (not root[0] in metAlready[parentHost])) and ((not parentHost in metAlreadyInd.keys()) or (not root[0] in metAlreadyInd[parentHost])):
			if parentHost in metAlreadyInd.keys():
				metAlreadyInd[parentHost].append(root[0])
			else:
				metAlreadyInd[parentHost]=[root[0]]
			indirectTrans[parentHost][root[0]].append(root[1]+numTransParent)
		if not (root[0] in origins.keys()):
			origins[root[0]]="Unsampled"
		#print origins
		#print "\n\n"
		if len(mt)>1:
			recurTransm(mt[0],root[0],0,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
		
	


# Update counts of transmission samples
def recurTransm(tree,parentHost,numTransParent,directTrans,indirectTrans,metAlready,metAlreadyInd,origins):
	splitT=splitTree(tree)
	#subtree1
	mt=metaData(splitT[0])
	root=extractInfo(mt[len(mt)-1])
	handleTree(mt,root,parentHost,numTransParent,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
	#subtree2
	mt=metaData(splitT[1])
	root=extractInfo(mt[len(mt)-1])
	handleTree(mt,root,parentHost,numTransParent,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)
	
	
	
## Read file to find burnin
inpF=open(args.inputF)
line="\n"
while len(line.split())<1 or line.split()[0]!="tree":
	line=inpF.readline()
	if line=="":
		print("Incorrect input file, is this a BEAST2 trees output file?")
		exit()
numTrees=0
normalL=len(line.split())
while len(line.split())==normalL and line.split()[0]=="tree":
	numTrees+=1
	line=inpF.readline()

burned=(float(args.burnin)/100)*numTrees
print("The first "+str(int(burned))+" trees out of "+str(numTrees)+" will be discarded as burnin.")

## re-Read file to find trees and collect values
inpF=open(args.inputF)
line="\n"
while len(line.split())<1 or line.split()[0]!="tree":
	line=inpF.readline()
	if line=="":
		print("Incorrect input file, is this a BEAST2 trees output file?")
		exit()
hosts=[]
directTrans={}
indirectTrans={}
totOrigins={}
roots={}
roots["Unsampled"]=0
numTrees=0
normalL=len(line.split())
while len(line.split())==normalL and line.split()[0]=="tree":

	if numTrees<int(burned):
		numTrees+=1
		line=inpF.readline()
		continue
	tree=line.split()[3]

	recurFindHosts(tree,hosts)
	
	if len(directTrans.keys())<len(hosts):
		for i in range(len(hosts)):
			if not (hosts[i] in directTrans.keys()):
				directTrans[hosts[i]]={}
				indirectTrans[hosts[i]]={}
				roots[hosts[i]]=0
				for j in range(len(hosts)):
					directTrans[hosts[i]][hosts[j]]=0
					indirectTrans[hosts[i]][hosts[j]]=[]
					
	mt=metaData(tree)
	root=extractInfo(mt[1])
	roots[root[0]]+=1
	metAlready={}
	metAlreadyInd={}
	origins={}
	if root[0]!="Unsampled":
		origins[root[0]]="Unsampled"
	recurTransm(mt[0],root[0],0,directTrans,indirectTrans,metAlready,metAlreadyInd,origins)#what happens if the root is sampled????????????????????????????????????????????????????????????????
	#print "final origins:"
	#print origins
	#print "\n\n\n\n"
	#exit()
	numTrees+=1
	line=inpF.readline()
	for host in origins.keys():
		if host in totOrigins.keys():
			if origins[host] in totOrigins[host].keys():
				totOrigins[host][origins[host]]+=1
			else:
				totOrigins[host][origins[host]]=1
		else:
			totOrigins[host]={}
			totOrigins[host][origins[host]]=1

numTrees=numTrees-int(burned)

#print("Hosts: ")
#print(hosts)
#print("\n\n Direct transmissions: ")
#print(directTrans)
#print("\n\n Indirect transmissions: ")
#print(indirectTrans)
#print("\n\n Roots: ")
#print(roots)
#print("\n\n nummTrees: ")
#print(numTrees)


#Record inferred network in text file
#hosts
outF=open(args.outputF+"_network.txt","w")
outF.write("Hosts: ")
for i in range(len(hosts)):
	outF.write(hosts[i]+", ")
outF.write("\n\n\n")
#roots
outF.write("Probabilities of being root: ")
for i in range(len(hosts)):
	outF.write(hosts[i]+" "+str(float(roots[hosts[i]])/numTrees)+", ")
outF.write("\n\n\n")
#direct transmissions
outF.write("Probabilities direct transmission: \n\n")
for i in range(len(hosts)):
	outF.write("From host "+hosts[i]+" to : \n")
	for j in range(len(hosts)):
		if j!= i and (directTrans[hosts[i]][hosts[j]])!=0:
			outF.write(hosts[j]+" "+str(float(directTrans[hosts[i]][hosts[j]])/numTrees)+", ")
	outF.write("\n\n")
outF.write("\n\n")
#indirect transmissions
outF.write("Probabilities indirect transmission: \n")
for i in range(len(hosts)):
	outF.write("From host "+hosts[i]+" to : \n")
	for j in range(len(hosts)):
		if j!= i and (len(indirectTrans[hosts[i]][hosts[j]]))!=0:
			outF.write(hosts[j]+" "+str(float(len(indirectTrans[hosts[i]][hosts[j]]))/numTrees)+", ")
	outF.write("\n\n")
outF.write("\n\n")
#origins
outF.write("Probabilities of direct transmittor to each sampled host: \n\n")
for i in range(len(hosts)):
	outF.write("To host "+hosts[i]+" from : \n")
	if hosts[i] in totOrigins.keys():
		for j in range(len(hosts)):
			if j!= i and (hosts[j] in totOrigins[hosts[i]].keys()):
				outF.write(hosts[j]+" "+str(float(totOrigins[hosts[i]][hosts[j]])/numTrees)+", ")
		if "Unsampled" in totOrigins[hosts[i]].keys():
			outF.write("Unsampled"+" "+str(float(totOrigins[hosts[i]]["Unsampled"])/numTrees)+", ")
		if "doubleOrigin" in totOrigins[hosts[i]].keys():
			outF.write("doubleOrigin"+" "+str(float(totOrigins[hosts[i]]["doubleOrigin"])/numTrees)+", ")
		outF.write("\n\n")
outF.write("\n\n")
outF.close()
print("\n\n"+"File "+args.outputF+"_network.txt containing output information successfully created!\n\n")




#Make graph for direct transmissions with graph_tools
plotD=args.plotDirect
if(plotD):
	minV=float(args.minValue)
	g = Graph()
	verteces=[]
	directEdges=[]
	for h1 in hosts:
		verteces.append(g.add_vertex())
	for h1 in range(len(hosts)):
		for h2 in range(len(hosts)):
			if (float(directTrans[hosts[h1]][hosts[h2]])/numTrees)>minV:
				directEdges.append(g.add_edge(verteces[h1], verteces[h2]))
				
				
	#Add information to graph
	eThick=float(args.edgeThickness)
	
	plotR=args.plotRoot
	rootProb = g.new_vertex_property("double")
	if plotR:
		rootProbText = g.new_vertex_property("string")
	for h1 in range(len(hosts)):
		rootProb[g.vertex(h1)]=float(roots[hosts[h1]])/numTrees
		if plotR:
			rootProbText[g.vertex(h1)]=("%.2f" % (float(roots[hosts[h1]])/numTrees))
	g.vertex_properties["root probability"] = rootProb
	if plotR:
		g.vertex_properties["root probability, text"] = rootProbText
	vertName = g.new_vertex_property("string")
	for h1 in range(len(hosts)):
		if plotR:
			vertName[g.vertex(h1)]=hosts[h1]+" "+("%.2f" % (float(roots[hosts[h1]])/numTrees))
		else:
			vertName[g.vertex(h1)]=hosts[h1]
	g.vertex_properties["host name"] = vertName
	
	transProb = g.new_edge_property("double")
	transProbText = g.new_edge_property("string")
	#print hosts
	for h1 in range(len(hosts)):
		for h2 in range(len(hosts)):
			#print hosts[h1]
			#print hosts[h2]
			#print float(directTrans[hosts[h1]][hosts[h2]])/numTrees
			if (float(directTrans[hosts[h1]][hosts[h2]])/numTrees)>minV:
				#print "sufficiently large"
				transProb[g.edge(h1,h2)] = (float(directTrans[hosts[h1]][hosts[h2]])/numTrees)*eThick +2.0
				transProbText[g.edge(h1,h2)] = ("%.2f" % ((float(directTrans[hosts[h1]][hosts[h2]])/numTrees)))
	g.edge_properties["direct transmission probability"] = transProb
	g.edge_properties["direct transmission probability, text"] = transProbText	
				
	#plot direct  transmission graph
	vFont=int(args.vertexFont)
	vColor=args.vertexColor
	eColor=args.edgeColor
	oSize=args.outputSize
	oformat=args.format
	#print("plotting")
	#print transProb
	graph_draw(g, vertex_text=vertName, vertex_font_size=vFont, output_size=(oSize, oSize), edge_pen_width=transProb, output=args.outputF+"_direct_transmissions."+oformat, vertex_fill_color=rootProb, vcmap=plt.get_cmap(vColor), vertex_pen_width=2.0, bg_color=[1., 1., 1., 1.], edge_text=transProbText, edge_color=transProb, ecmap=plt.get_cmap(eColor), fmt=oformat)
	print("\n\n"+"File "+args.outputF+"_direct_transmissions."+oformat+" containing graph of direct transmission successfully created!\n\n")	
	



#Make graph for indirect transmissions with graph_tools
plotI=args.plotIndirect
if(plotI):
	minV=float(args.minValue)
	g = Graph()
	verteces=[]
	directEdges=[]
	for h1 in hosts:
		verteces.append(g.add_vertex())
	for h1 in range(len(hosts)):
		for h2 in range(len(hosts)):
			if (float(len(indirectTrans[hosts[h1]][hosts[h2]])+directTrans[hosts[h1]][hosts[h2]])/numTrees)>minV:
				directEdges.append(g.add_edge(verteces[h1], verteces[h2]))
				
				
	#Add information to graph
	eThick=float(args.edgeThickness)
	
	rootProb = g.new_vertex_property("double")
	for h1 in range(len(hosts)):
		rootProb[g.vertex(h1)]=float(roots[hosts[h1]])/numTrees
	g.vertex_properties["root probability"] = rootProb
	vertName = g.new_vertex_property("string")
	for h1 in range(len(hosts)):
		vertName[g.vertex(h1)]=hosts[h1]
	g.vertex_properties["host name"] = vertName
	
	transProb = g.new_edge_property("double")
	transProbText = g.new_edge_property("string")
	#print hosts
	for h1 in range(len(hosts)):
		for h2 in range(len(hosts)):
			#print hosts[h1]
			#print hosts[h2]
			#print float(directTrans[hosts[h1]][hosts[h2]]+len(indirectTrans[hosts[h1]][hosts[h2]]))/numTrees
			if (float(directTrans[hosts[h1]][hosts[h2]]+len(indirectTrans[hosts[h1]][hosts[h2]]))/numTrees)>minV:
				#print "sufficiently large"
				transProb[g.edge(h1,h2)] = (float(directTrans[hosts[h1]][hosts[h2]]+len(indirectTrans[hosts[h1]][hosts[h2]]))/numTrees)*eThick +2.0
				transProbText[g.edge(h1,h2)] = ("%.2f" % ((float(directTrans[hosts[h1]][hosts[h2]]+len(indirectTrans[hosts[h1]][hosts[h2]]))/numTrees)))
	g.edge_properties["direct + indirect transmission probability"] = transProb
	g.edge_properties["direct + indirect transmission probability, text"] = transProbText	
				
	#plot indirect  transmission graph
	vFont=int(args.vertexFont)
	vColor=args.vertexColor
	eColor=args.edgeColor
	oSize=args.outputSize
	oformat=args.format
	
	#print("plotting again")
	#print transProb
	graph_draw(g, vertex_text=vertName, vertex_font_size=vFont, output_size=(oSize, oSize), edge_pen_width=transProb, output=args.outputF+"_indirect_transmissions."+oformat, vertex_fill_color=rootProb, vcmap=plt.get_cmap(vColor), vertex_pen_width=2.0, bg_color=[1., 1., 1., 1.], edge_text=transProbText, edge_color=transProb, ecmap=plt.get_cmap(eColor), fmt=oformat)	
	print("\n\n"+"File "+args.outputF+"_indirect_transmissions."+oformat+" containing graph of indirect transmission successfully created!\n\n")
	
	
exit()



