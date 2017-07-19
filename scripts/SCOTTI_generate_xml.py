import sys
import os
import math
import argparse
#import random
#import re
#import os.path
#import math

#does it represent an integer?
def IsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

#does it represent an integer?
def IsFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

#make tree string into format for mugration
def removeStates(treeString):
	l=len(treeString)
	newString=""
	i=0
	while i < l :
		if treeString[i]=="[":
			while treeString[i]!="]":
				i+=1
			i+=1
		newString+=treeString[i]
		i+=1
	return newString		

#scale branch lengths 
def scaleTreeString(treeString,scale):
	l=len(treeString)
	newString=""
	i=0
	while i < l :
		if treeString[i]==":":
			newString+=":"
			numb=""
			i+=1
			while treeString[i]!="," and treeString[i]!=")" and treeString[i]!=";":
				numb+=treeString[i]
				i+=1
			newString+=str(float(numb)*scale)
		newString+=treeString[i]
		i+=1
	return newString


parser = argparse.ArgumentParser()
#parser.add_argument('--template',"-t", help='the file SCOTTI_template.xml with path. It will be used as template to generate the xml file.',default=os.path.dirname(os.path.abspath(__file__))+"SCOTTI_template.xml")
parser.add_argument('--fasta',"-f", help='Input fasta alignment. If not specified, looks in th working directory for SCOTTI_aligment.fas',default=os.getcwd()+"SCOTTI_alignment.fas")
parser.add_argument('--output',"-o", help='output file that will contain the newly created SCOTTI xml to be run with BEAST2.',default=os.getcwd()+"SCOTTI_analysis")
parser.add_argument('--overwrite',"-ov", action='store_true', help="Overwrite output files.")#, dest='overW'
parser.set_defaults(overwrite=False)
parser.add_argument('--dates',"-d", help='Input csv file with sampling dates associated to each sample name (same name as fasta file). If not specified, looks for SCOTTI_dates.csv',default=os.getcwd()+"SCOTTI_dates.csv")
parser.add_argument('--hosts',"-ho", help='Input csv file with sampling hosts associated to each sample name (same name as fasta file). If not specified, looks for SCOTTI_hosts.csv',default=os.getcwd()+"SCOTTI_hosts.csv")
parser.add_argument('--hostTimes',"-ht", help='Input csv file with the earliest times in which hosts are infectable, and latest time when hosts are infective. Use same host names as those used for the sampling hosts file. If not specified, looks for SCOTTI_hostTimes.csv',default=os.getcwd()+"SCOTTI_hostTimes.csv")
parser.add_argument('--maxHosts',"-m", help='Maximum number of hosts (between sampled, non-sampled, and non-observed) allowed.', type=int, default=-1)
parser.add_argument('--unlimLife',"-u", action='store_true', help="non-sampled, non-observed hosts have unlimited life span (unlimited contribution in time to the outbreak). By default non-sampled non-observed hosts have limited duration. Usually asymptomaic patients have limited contributions, while to model environmental contamination it's better to use unlimited hosts. Unlimited life span can also decrease computational demand.")#, dest='limL', action='store_true'
parser.set_defaults(unlimLife=False)
parser.add_argument('--penalizeMigration',"-p", action='store_true', help="penalize lineages that are still in a deme after deme closure. By default, don't (original version).")#, dest='limL', action='store_true'
parser.set_defaults(penalizeMigration=False)
parser.add_argument('--numIter',"-n", help='Number of iterations of the MCMC. Default=100,000,000.', type=int, default=100000000)
parser.add_argument('--mutationModel',"-mu", help='String represnting the mutation model to be used. Possibilities are \"HKY\" or \"JC\". Further alternaties can be specified by modifying the output xml manually.', default="HKY")
parser.add_argument('--fixedAs', "-fA", type=int, help = 'Number of sites fixed for A in the genome and not included in the alignment.',default=0)
parser.add_argument('--fixedCs', "-fC", type=int, help = 'Number of sites fixed for C in the genome and not included in the alignment.',default=0)
parser.add_argument('--fixedGs', "-fG", type=int, help = 'Number of sites fixed for G in the genome and not included in the alignment.',default=0)
parser.add_argument('--fixedTs', "-fT", type=int, help = 'Number of sites fixed for T in the genome and not included in the alignment.',default=0)
parser.add_argument('--tracelog',"-l", help='Logging frequency in the trace file. Default=numIter/10000.', type=int, default=-1)
parser.add_argument('--screenlog',"-sl", help='Logging frequency to the screen. Default=numIter/5000.', type=int, default=-1)
parser.add_argument('--treelog',"-tl", help='Logging frequency in the tree file. Default=numIter/1000.', type=int, default=-1)



args = parser.parse_args()


#if !(os.path.isfile(args.template)):
#    print "Error: template file SCOTTI_template.xml not found, please specify it correctly in input using the option -t .\n"
#    exit()
if (not (os.path.isfile(args.fasta))):
    print("Error: input fasta file not found, please specify it correctly in input using the option -f .\n")
    exit()
if (os.path.isfile(args.output + ".xml")) and (not args.overwrite):
    print("Output xml file already exists, to allow overwriting of output file use option -ov, otherwise specify a different output file.\n")
    exit()
if not (os.path.isfile(args.dates)):
    print("Error: input dates csv file not found, please specify it correctly in input using the option -d .\n")
    exit()
if not (os.path.isfile(args.hosts)):
    print("Error: input hosts csv file not found, please specify it correctly in input using the option -s .\n")
    exit()
if not (os.path.isfile(args.hostTimes)):
    print("Error: input hosts csv file not found, please specify it correctly in input using the option -i .\n")
    exit()
if (args.numIter<=0):
    print("Error: the number of MCMC iterations specified is negative or null.\n")
    exit()
elif (args.numIter<=1000):
    print("Very small number of MCMC iterations, might want to specifiy something bigger.\n")
if (args.mutationModel!="HKY" and args.mutationModel!="JC"):
    print("Mutation model not supported. Please specify something among \"HKY\" or \"JC\".\n")
    exit()
if (args.fixedAs<0 or args.fixedCs<0 or args.fixedGs<0 or args.fixedTs<0):
    print("Error: the number of fixed sites for each base type must be positive.\n")
    exit()

#read fasta file
sqFile=open(args.fasta)
seqs={}
line=sqFile.readline()
linesplit=line.split()
while line!="" and (len(linesplit)<1 or linesplit[0][0]!=">"):
    line=sqFile.readline()
    linesplit=line.split()
if line=="":
    print("Wrong format in fasta file\n")
    exit()
while line!="":
    if len(linesplit[0])>1:
        name=linesplit[0].replace(">","")
    else:
        name=linesplit[1]
    seq=""
    line=sqFile.readline()
    linesplit=line.split()
    while line!="" and (len(linesplit)>0  and linesplit[0][0]!=">") and line!="\n":
        seq+=linesplit[0]
        line=sqFile.readline()
        linesplit=line.split()
    if name in seqs.keys():
        print("Error, multiple entries with same name "+name+" in fasta file. Only first part of the name line is used for the name.")
        exit()
    seqs[name]=seq
    while line!="" and  (len(linesplit)<1 or linesplit[0][0]!=">"):
        line=sqFile.readline()
        linesplit=line.split()
sqFile.close()

#read in file with host-sample information
hFile=open(args.hosts)
hosts={}
line=hFile.readline()
linesplit=("".join(line.split())).split(",")  
while line!="" and len(linesplit)>1: 
    sam=linesplit[0]
    hos=linesplit[1]
    if not (sam in seqs.keys()):
        print("Error: sample "+sam+" found in sample-host file but not found in the fasta file.")
        exit()
    if sam in hosts.keys():
        print("Error, multiple entries with same name "+sam+" in hosts file. ")
        exit()
    hosts[sam]=hos
    line=hFile.readline()
    linesplit=("".join(line.split())).split(",")
    while line=="\n":
        line=hFile.readline()
        linesplit=("".join(line.split())).split(",")
hFile.close()

#read in sampling dates
dFile=open(args.dates)
dates={}
line=dFile.readline()
linesplit=("".join(line.split())).split(",")  
while line!="" and len(linesplit)>1:  
    sam=linesplit[0]
    dat=linesplit[1]
    if not (IsFloat(dat)):
        print("Error: date "+str(dat)+" for sample "+sam+" cannot be converted to a number")
        exit()
    if not (sam in seqs.keys()):
        print("Error: sample "+sam+" found in sampling dates file but not found in the fasta file.")
        exit()
    if sam in dates.keys():
        print("Error, multiple entries with same name "+sam+" in dates file. ")
        exit()
    dates[sam]=float(dat)
    line=dFile.readline()
    linesplit=("".join(line.split())).split(",")
    while line=="\n":
        line=dFile.readline()
        linesplit=("".join(line.split())).split(",")
dFile.close()
    
#read in file with host times information
hFile=open(args.hostTimes)
hostT={}
line=hFile.readline()
linesplit=("".join(line.split())).split(",")  
while line!="" and len(linesplit)>2: 
    hos=linesplit[0]
    dat=float(linesplit[1])
    dat2=float(linesplit[2])
    if (not (IsFloat(dat))) or (not (IsFloat(dat2))):
        print("Error: date "+str(dat)+" or "+str(dat2)+" for sample "+sam+" cannot be converted to a number")
        exit()
    #if not (hos in seqs.keys()):
    #    print "Error: sample "+hos+" found in host times file but not found in the fasta file."
    #    exit()
    if hos in hostT.keys():
        print("Error, multiple entries with same name "+hos+" in hosts-times file. ")
        exit()
    hostT[hos]=[dat, dat2]
    line=hFile.readline()
    linesplit=("".join(line.split())).split(",")
    while line=="\n":
        line=hFile.readline()
        linesplit=("".join(line.split())).split(",")
hFile.close()


    
#check that the details are correct
for s in seqs.keys():
    if not (s in dates.keys()):
        print("Error, no date specified for sample "+s)
        exit()
    if not (s in hosts.keys()):
        print("Error, no host specified for sample "+s)
        exit()
for s in hosts.keys():
    if not (s in seqs.keys()):
        print("Error, no sequence specified for sample "+s)
        exit()
    if not (s in dates.keys()):
        print("Error, no date specified for sample "+s)
        exit()
for s in dates.keys():
    if not (s in seqs.keys()):
        print("Error, no sequence specified for sample "+s)
        exit()
    if not (s in hosts.keys()):
        print("Error, no host specified for sample "+s)
        exit()
for s in hosts.keys():
    if not (hosts[s] in hostT.keys()):
        print("Warning: no dates specified for host "+hosts[s]+". I will assume this host is available throughout")
        hostT[hosts[s]]=[-1000000000000,1000000000000]


for h in hostT.keys():
    found=0
    for s in hosts.keys():
        if hosts[s]==h:
            found=1
    if found==0:
        print("Warning: no sample specified for host "+h+". I will assume this host could have been infected or not. If you are sure this host was infected, then please include a non-informative sample (all \"N\"s) from it.")

  
# logging frequency
tracelog  = args.tracelog  if args.tracelog  > 0 else int(args.numIter/10000)
screenlog = args.screenlog if args.screenlog > 0 else int(args.numIter/5000)
treelog   = args.treelog   if args.treelog   > 0 else int(args.numIter/1000)
    
    
    
xml=open(args.output+".xml","w")
xml.write("<beast version=\'2.0\' namespace=\'beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util\'>\n\n"+"  <data id=\"alignmentVar\" dataType=\"nucleotide\">\n")
for s in seqs.keys():
    xml.write("    <sequence taxon=\'"+s+"\' value=\'"+seqs[s]+"\'/>\n")
xml.write("  </data>\n")
st=""
if args.fixedAs>0 and args.fixedCs>0 and args.fixedGs>0 and args.fixedTs>0:
    st="constantSiteWeights=\'"+str(args.fixedAs)+" "+str(args.fixedCs)+" "+str(args.fixedGs)+" "+str(args.fixedTs)+"\'"
xml.write("  <data id=\'alignment\' spec=\'FilteredAlignment\' filter=\'-\' data=\'@alignmentVar\' "+st+"/>\n\n")

#hosts
xml.write("   <typeTraitSet spec=\'TraitSetWithLimits\' id=\'typeTraitSet\' traitname=\"type\" value=\"")
hL=list(hosts.keys())
for s in hL:
    xml.write(s+"="+hosts[s])
    if s!=hL[len(hL)-1]:
        xml.write(",")
xml.write("\">\n     <taxa spec=\'TaxonSet\' alignment=\'@alignment\'/>\n   </typeTraitSet>\n\n")

#dates
xml.write("   <timeTraitSet spec=\'TraitSetWithLimits\' id=\'timeTraitSet\' traitname=\"date-forward\" value=\"")
dL=list(dates.keys())
for s in range(len(dL)):
    xml.write(dL[s]+"="+str(dates[dL[s]]))
    if (s!=(len(dL)-1)):
        xml.write(",")
xml.write("\">\n     <taxa spec=\'TaxonSet\' alignment=\'@alignment\'/>\n   </timeTraitSet>\n\n")

#end infectability
xml.write("   <endInfectionsTraitSet spec=\'InfectionTraitSet\' id=\'endInfectionsTraitSet\' type=\'end\' traitname=\"date-forward\" value=\"")
dL=list(hostT.keys())
for s in dL:
    xml.write(s+"="+str(hostT[s][1]))
    if s!=dL[len(dL)-1]:
        xml.write(",")
xml.write("\">\n 	  <taxaDates idref=\'timeTraitSet\'/>\n 	  <taxaTypes idref=\'typeTraitSet\'/>\n   </endInfectionsTraitSet>\n\n")

#start infectability
xml.write("   <startInfectionsTraitSet spec=\'InfectionTraitSet\' id=\'startInfectionsTraitSet\' type=\'start\' traitname=\"date-forward\" value=\"")
dL=list(hostT.keys())
for s in dL:
    xml.write(s+"="+str(hostT[s][0]))
    if s!=dL[len(dL)-1]:
        xml.write(",")
xml.write("\">\n 	  <taxaDates idref=\'timeTraitSet\'/>\n 	  <taxaTypes idref=\'typeTraitSet\'/>\n   </startInfectionsTraitSet>\n\n")

#check that times fit
for s in seqs.keys():
    h=hosts[s]
    if dates[s]<=hostT[hosts[s]][0]:
        print("Error, sampling time "+str(dates[s])+" of sample "+s+" must be after start of infectability time "+str(hostT[hosts[s]][0])+" of host "+h)
        exit()
    if dates[s]>=hostT[hosts[s]][1]:
        print("Error, sampling time "+str(dates[s])+" of sample "+s+" must be before end of infectiousness time "+str(hostT[hosts[s]][0])+" of host "+h)
        exit()

#include mutation model
xml.write("   <siteModel spec=\"SiteModel\" id=\"siteModel\">\n     <mutationRate spec=\'RealParameter\' id=\"mutationRate\" value=\"1.0\"/>\n     <substModel spec=\"")
if args.mutationModel=="HKY":
   xml.write("HKY\">\n       <kappa spec=\'RealParameter\' id=\"hky.kappa\" value=\"1.0\"/>\n")
elif args.mutationModel=="JC":
   xml.write("jc69\">\n")
elif args.mutationModel=="GTR":
   xml.write("GTR\" rateAC=\"@rateAC\" rateAG=\"@rateAG\" rateAT=\"@rateAT\" rateCG=\"@rateCG\" rateGT=\"@rateGT\" >\n       <parameter estimate=\'false\' name=\"rateCT\" lower=\"0.0\" value=\"1.0\" id=\"rateCTfixed\"/>\n")
else:
    print("Error, mutation model "+args.mutationModel+" not recognized.")
    exit()
xml.write("       <frequencies estimate=\"false\" spec=\'Frequencies\'>\n	 <frequencies spec=\'RealParameter\' id=\""+args.mutationModel+".freq\" value=\"0.25 0.25 0.25 0.25\"/>\n       </frequencies>\n     </substModel>\n   </siteModel>\n\n")
              
#migration model
nHos=len(hostT.keys())
meanL=0.0
for h in hostT.keys():
    meanL+=hostT[h][1]-hostT[h][0]
meanL=meanL/len(hostT.keys())
ma=args.maxHosts
if ma<nHos:
    ma=nHos+5
    print("Maximum number of hosts specified is less than the minimum number of hosts. Setting it to "+str(ma))
if ma<(nHos+2):
    print("Warning: the maximum number of hosts is less than "+str(nHos+2)+". If there are problems in the execution of the BEAST2 analysis, try a higher values of maxHosts.")
start=int((nHos+ma)/2)
if start<(nHos+2) and ma>=(nHos+2):
    start=nHos+2
xml.write("   <migrationModelUniform spec=\'MigrationModelUniform\' id=\'migModel\' minDemes=\'"+str(nHos)+"\'>\n     <rate spec=\'RealParameter\' value=\""+str(1.0/(meanL*nHos))+"\" dimension=\"1\" id=\"rate\"/>\n     <popSize spec=\'RealParameter\' value=\"1.0\" dimension=\"1\" id=\"popSize\"/>\n     <numDemes spec=\'IntegerParameter\' value=\""+str(start)+"\" dimension=\"1\" id=\"numDemes\" lower=\'"+str(nHos)+"\' upper=\'"+str(ma)+"\'/>\n     <trait idref=\'startInfectionsTraitSet\' type=\'start\'/>\n     <trait idref=\'endInfectionsTraitSet\' type=\'end\'/>\n   </migrationModelUniform>\n\n")            

#priors
xml.write("   <input spec=\'CompoundDistribution\' id=\'parameterPriors\'>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@mutationRate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"-3.0\" S=\"6.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@popSize\">\n       <distr spec=\"LogNormalDistributionModel\"  M=\"0.0\" S=\"6.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\""+str(math.log(1.0/(meanL*nHos)))+"\" S=\"0.2\"/>\n     </distribution>\n")
if args.mutationModel=="HKY":
    xml.write("     <distribution spec=\'beast.math.distributions.Prior\' x=\"@hky.kappa\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n   </input>\n\n")
elif args.mutationModel=="JC":
   xml.write("   </input>\n\n")
elif args.mutationModel=="GTR":
    xml.write("     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAC\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAG\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAT\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateCG\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateGT\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n   </input>\n\n")      
              
#Likelihoods
xml.write("<input spec=\'TreeLikelihood\' id=\"treeLikelihood\">\n     <data idref=\"alignment\"/>\n     <tree idref=\"tree\"/>\n     <siteModel idref=\'siteModel\'/>\n   </input>\n\n   <input spec=\'StructuredCoalescentTreeDensityNew\' id=\'treePrior\' limitedLifespan=\'"+str(args.unlimLife)+"\' penalizeMigration=\'"+str(args.penalizeMigration)+"\'>\n     <multiTypeTreeConcise idref=\"tree\"/>\n     <migrationModelUniform idref=\"migModel\"/>\n   </input>\n\n   <run spec=\"MCMC\" id=\"mcmc\" chainLength=\""+str(args.numIter)+"\" storeEvery=\"10000\">\n\n")

#migration model
xml.write("     <init spec=\'StructuredCoalescentMultiTypeTreeConcise\' id=\'tree\' nTypes=\""+str(start)+"\">\n         <migrationModelUniform spec=\'MigrationModelUniform\' minDemes=\'"+str(nHos)+"\'>\n             <rate spec=\'RealParameter\' value=\""+str(1.0/(meanL*nHos))+"\" dimension=\"1\"/>\n             <popSize spec=\'RealParameter\' value=\"1.0\" dimension=\"1\"/>\n             <numDemes spec=\'IntegerParameter\' value=\""+str(start)+"\" dimension=\"1\"/>\n         </migrationModelUniform>\n         <trait idref=\'typeTraitSet\'/>\n         <trait idref=\'timeTraitSet\'/>\n     </init>\n\n")              
              
#state nodes              
xml.write("     <state>\n       <stateNode idref=\"tree\"/>\n       <stateNode idref=\"rate\"/>\n       <stateNode idref=\"popSize\"/>\n       <stateNode idref=\"numDemes\"/>\n       <stateNode idref=\"mutationRate\"/>\n")
if args.mutationModel=="HKY":
    xml.write("       <stateNode idref=\"hky.kappa\"/>\n")
elif args.mutationModel=="GTR":              
    xml.write("       <stateNode idref=\"rateAC\"/>\n       <stateNode idref=\"rateAG\"/>\n       <stateNode idref=\"rateAT\"/>\n       <stateNode idref=\"rateCG\"/>\n       <stateNode idref=\"rateGT\"/>\n")
xml.write("       <stateNode idref=\""+args.mutationModel+".freq\"/>\n     </state>\n\n")     
              
#compound distribution
xml.write("     <distribution spec=\'CompoundDistribution\' id=\'posterior\'>\n       <distribution idref=\"treeLikelihood\"/>\n       <distribution idref=\'treePrior\'/>\n       <distribution idref=\"parameterPriors\"/>\n     </distribution>\n\n")
              
#operators
xml.write("     <operator spec=\'ScaleOperator\' id=\'RateScaler\' parameter=\"@rate\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"ScaleOperator\" id=\"PopSizeScaler\" parameter=\"@popSize\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"ShortRangeUniformOperator\" id=\"NumDemesScaler\" parameter=\"@numDemes\" weight=\"1\"/>\n     <operator spec=\"ScaleOperator\" id=\"muRateScaler\" parameter=\"@mutationRate\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\" parameter=\"@"+args.mutationModel+".freq\" delta=\"0.01\" weight=\"0.1\"/>\n     <operator id=\'treeScaler.t\' spec=\'ScaleOperator\' scaleFactor=\"0.5\" weight=\"3\" tree=\"@tree\"/>\n     <operator id=\'treeRootScaler.t\' spec=\'ScaleOperator\' scaleFactor=\"0.5\" weight=\"3\" tree=\"@tree\" rootOnly=\'true\'/>\n     <operator id=\'UniformOperator.t\' spec=\'Uniform\' weight=\"30\" tree=\"@tree\"/>\n     <operator id=\'SubtreeSlide.t\' spec=\'SubtreeSlide\' weight=\"15\" gaussian=\"true\" size=\"1.0\" tree=\"@tree\"/>\n     <operator id=\'narrow.t\' spec=\'Exchange\' isNarrow=\'true\' weight=\"15\" tree=\"@tree\"/>\n     <operator id=\'wide.t\' spec=\'Exchange\' isNarrow=\'false\' weight=\"3\" tree=\"@tree\"/>\n     <operator id=\'WilsonBalding.t\' spec=\'WilsonBalding\' weight=\"3\" tree=\"@tree\"/>\n")
if args.mutationModel=="HKY":
    xml.write("     <operator spec=\'ScaleOperator\' id=\'kappaScaler\' parameter=\"@hky.kappa\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n")
elif args.mutationModel=="GTR": 
    xml.write("     <operator spec=\'ScaleOperator\' id=\'ACScaler\' parameter=\"@rateAC\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'AGScaler\' parameter=\"@rateAG\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'ATScaler\' parameter=\"@rateAT\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'CGScaler\' parameter=\"@rateCG\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'GTScaler\' parameter=\"@rateGT\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n")
xml.write("\n")

#loggers
xml.write("     <logger logEvery=\""+str(tracelog)+"\" fileName=\""+args.output+".log\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>       <log idref=\"treeLikelihood\"/>\n       <log idref=\"migModel\"/>\n       <log idref=\"mutationRate\"/>\n       <log spec=\'TreeHeightLogger\' tree=\'@tree\'/>\n       <log idref=\""+args.mutationModel+".freq\"/>\n")
if args.mutationModel=="HKY":
    xml.write("       <log idref=\"hky.kappa\"/>\n")
elif args.mutationModel=="GTR":               
    xml.write("       <log idref=\"rateAC\"/>\n       <log idref=\"rateAG\"/>\n       <log idref=\"rateAT\"/>\n       <log idref=\"rateCG\"/>\n       <log idref=\"rateGT\"/>\n")
xml.write("     </logger>\n")
     
xml.write("     <logger logEvery=\""+str(treelog)+"\" fileName=\""+args.output+".trees\" mode=\"tree\">\n       <log idref=\'treePrior\'/>\n     </logger>\n")

xml.write("     <logger logEvery=\""+str(screenlog)+"\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>\n       <log idref=\"treeLikelihood\"/>\n       <log spec=\'TreeHeightLogger\' tree=\'@tree\'/>\n       <log idref=\"migModel\"/>\n       <log idref=\"mutationRate\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@treePrior\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@posterior\"/>\n     </logger>\n\n   </run>\n </beast>\n\n")




            
      
              
              
              
              
              
              
                  
      
			
exit()





		