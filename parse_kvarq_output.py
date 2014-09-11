"""
Paul McAdam
4th August 2014
usage parse_kvarq_output.py <collated scan table>
"""


import sys
import collections

kvarq_table=sys.argv[1]



def returnPhylo(phylo):
	if phylo.startswith("?"):
		return "Low Coverage"
	elif len(phylo.split("//"))!=1:
		return "Mixed Sample"
	elif phylo.startswith("lineage"):
		split=phylo.split('/')
		if len(split)==2:
			return "Beijing"
		else:
			split=phylo.split()
			return "Lineage " + split[1]

def returnResistance(res):
#	strain_resistant_to=set([])
	strain_resistant_to={"rifampicin":0,"isoniazid":0,"streptomycin":0,"ethambutol":0,"fluoroquinolones":0,"Kanamycin/Amikacin":0}
	if res==[]:
		return "Susceptible"
	else:
		res=res.strip("[]")
		res=res.replace("u\'","")
		res=res.replace("u'","")
		res=res.replace("'","")
		res=res.replace("\"","")
		res=res.replace("]","")
		res=res.replace("[","")
		resistances=res.split(",")
#		print resistances
		for resistance in resistances:
			resistance=resistance.split("resistance")[0].strip()
			if resistance.lower() in strain_resistant_to:
				strain_resistant_to[resistance.lower()]=1
#			strain_resistant_to.add(resistance)
	od = collections.OrderedDict(sorted(strain_resistant_to.items()))
	return od

def resistanceLoci(res):
	snp_in_locus={"inhA":0,"katG":0,"rpoB":0,"rpoA":0,"rpoC":0,"gyrA":0,"gyrB":0,"rpsL":0,"rrs":0,"embB":0,}
	if res==[]:
		return snp_in_locus
	else:
		for locus in snp_in_locus:
			if locus in res:
				snp_in_locus[locus]=1
	od = collections.OrderedDict(sorted(snp_in_locus.items()))
	return od


out_file=open(sys.argv[1].split(".")[0]+"_parsed.csv","w")

with open(kvarq_table,'r') as in_file:
	out_file.write("Isolate,Lineage,")
	resistant_to=returnResistance("HEADER")
	for x in resistant_to:
		out_file.write(x+",")
	snp_loci=resistanceLoci("HEADER")
	for y in snp_loci:
		out_file.write(y+",")
	out_file.write("\n")
	lines=in_file.readlines()
	for line in lines[1:]:
		isolate, lineage, resistance=line.split('\t')[0], line.split('\t')[1], line.split('\t')[2]
		isolate=isolate.replace("/home/UNIMELB/paulm-u/TB_sequencing_data/all_data/","")
		isolate=isolate.split("_")[0]
		if isolate.startswith("02MTB")==False:
			isolate="02MTB"+((4-len(isolate))*"0")+isolate
#		returnResistance(resistance)
#		print line
		out_file.write(isolate+","+returnPhylo(lineage)+",")
		resistant_to=returnResistance(resistance)
		for x in resistant_to:
			out_file.write(str(resistant_to[x])+",")
		snp_loci=resistanceLoci(resistance)
		for y in snp_loci:
			out_file.write(str(snp_loci[y])+",")
		out_file.write("\n")

out_file.close()


"""
To add:
positional SNPs identified
"""

streptomycin={}
ethambutol={}
isoniazid={}
fluoroquinolone={}
rifampicin={}
kanamycin={}

