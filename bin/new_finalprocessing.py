#!/usr/bin/python
import sys
import re

def translate(seq):
      
	table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	} 
	protein =""
	if len(seq)%3 == 0:
		for i in range(0, len(seq), 3):
			codon = seq[i:i + 3]
			if "N" in codon:
				protein+="?"
			else :
				protein+= table[codon]
	else:
		for j in range(0, len(seq)-len(seq)%3, 3):
			codon = seq[j:j + 3]
			if "N" in codon:
				protein+="?"
			else:
				protein+= table[codon]
	return protein

#######################################################################################################################################################################################
##### INSECTOR - STEP 11 - Modify the ends for start and stop codons ##################################################################################################################
#######################################################################################################################################################################################

dict_genseq={}
value=""
key=None
count = 0
with open(str(sys.argv[1]), 'r') as f:
	for line in f:
		line2=line.strip('\n')
		#print line2	
		if line2[0]==">":
			#count += 1
			if key: 
				dict_genseq[key] = value
			key=line2.split(' ')[0]
			#print key
			#a=a+1
			value=""
		else:
			value=value + line2[:]
			
	
	dict_genseq[key] = value
	#print len(dict_genseq)
	#print count
	
sorted_dict_genseq= sorted(dict_genseq.items())
final_hit_seq={}
dna=""
start_code=""
stop_code=""
key1=None
value1=""
with open (str(sys.argv[2]),"r") as fi:
	for line in fi:
		line2=line.strip('\n')	
		if line2[0]==">":
			if key1: 
				final_hit_seq[key1] = value1
			key1=line2
			value1=""
		else:
			value1=value1 + line2[:]
			
	
	final_hit_seq[key1] = value1

sorted_final_hit_seq=sorted(final_hit_seq.items())

f_table=open(str(sys.argv[4])+".final_table.txt", "w+")
f_protein=open(str(sys.argv[4])+".final_proteins.pep", "w+")
f_starminus_protein=open(str(sys.argv[4])+".ORs.starRemoved.pep", "w+")
f_gff=open(str(sys.argv[4])+".final_gff_file.gff", "w+")
strand=""
with open (str(sys.argv[3]), "r") as fil:
	for line in fil:
		line2=line.strip("\n").split("\t")
		scaf=str(line2[0])
		start_exon=int(line2[7].split(" ")[0])
		end_exon=int(line2[7].split(" ")[-1])
		start_code=""
		stop_code=""
		modify_start=0
		modify_end=0
		#print "iiiiiiiiiiiiiii", end_exon, scaf
		for k, v in dict_genseq.iteritems() :
			header=k
			if (scaf in header):
				if line2[3]=="forward" and int(line2[4])>=350:
					if "N" in v[end_exon-3 :end_exon+60]: 
						nt=v[end_exon-3:end_exon+60].find('N')
						stop_dna= translate(v[end_exon-3:end_exon+nt].upper())
					else:
						stop_dna= translate(v[end_exon-3:end_exon+60].upper())
					if "*" in stop_dna:
						stop_code=stop_dna[1:stop_dna.find('*')]
						last_exon_boundary=end_exon+(stop_dna.find('*'))*3
						line2[7]=line2[7].replace(line2[7].split(" ")[-1], str(last_exon_boundary),1)
						line2[2]=str(last_exon_boundary)
						modify_end=1

					if(v[start_exon-1:start_exon+2].upper()!="ATG"):
						if "N" in v[start_exon-61:start_exon-1]: 
							nt1=v[start_exon-61:start_exon-1].rfind('N')
							start_dna= translate(v[start_exon-61+nt1+1:start_exon-1].upper())
						else:									
							start_dna=  translate(v[start_exon-61:start_exon-1].upper())
						if "M" in start_dna:
							start_code=start_dna[start_dna.rfind('M'):]
							first_exon_boundary=start_exon-(len(start_dna)-start_dna.rfind('M'))*3
							line2[7]=line2[7].replace(line2[7].split(" ")[0], str(first_exon_boundary),1)
							line2[1]=str(first_exon_boundary)
							modify_start=1
						if start_dna.rfind('M')<start_dna.rfind('*'):
							start_code=""	
							first_exon_boundary=start_exon	
							line2[7]=line2[7].replace(line2[7].split(" ")[0], str(first_exon_boundary),1)
							line2[1]=str(first_exon_boundary)								
							modify_start=0	

						first_exon_length=int(line2[7].split(" ")[1]) - int(line2[7].split(" ")[0])
						# if modify_start!=1 and first_exon_length>30: 
						# 	length_of_trim=

					
				elif line2[3]!="forward" and int(line2[4])>=350:
					if(v[end_exon-3:end_exon].upper()!="CAT"):
						if "N" in v[end_exon:end_exon+60]: 
							nt=v[end_exon:end_exon+60].find('N')
							start_dna= translate(v[end_exon:end_exon+nt].upper().replace('A','1').replace('T','2').replace('G','3').replace('C','4').replace('1','T').replace('2','A').replace('3','C').replace('4','G')[::-1])
						else:
							start_dna= translate(v[end_exon:end_exon+60].upper().replace('A','1').replace('T','2').replace('G','3').replace('C','4').replace('1','T').replace('2','A').replace('3','C').replace('4','G')[::-1])
						if "M" in start_dna:
							start_code=start_dna[start_dna.rfind('M'):]	
							first_exon_boundary=end_exon+(len(start_dna)-start_dna.rfind('M'))*3	
							line2[7]=line2[7].replace(line2[7].split(" ")[-1], str(first_exon_boundary),1)
							line2[2]=str(first_exon_boundary)
							modify_start=1							 
						if start_dna.rfind('M')<start_dna.rfind('*'):
							start_code=""
							first_exon_boundary=end_exon
							line2[7]=line2[7].replace(line2[7].split(" ")[-1], str(first_exon_boundary),1)
							line2[2]=str(first_exon_boundary)							
							modify_start=0

						last_exon_length=int(line2[7].split(" ")[-1]) - int(line2[7].split(" ")[-2])

					if "N" in v[start_exon-61:start_exon+2]: 
						nt1=v[start_exon-61:start_exon+2].rfind('N')
							
						stop_dna= translate(v[start_exon-61+nt1+1:start_exon+2].upper().replace('A','1').replace('T','2').replace('G','3').replace('C','4').replace('1','T').replace('2','A').replace('3','C').replace('4','G')[::-1])
					else:
						stop_dna=  translate(v[start_exon-61:start_exon+2].upper().replace('A','1').replace('T','2').replace('G','3').replace('C','4').replace('1','T').replace('2','A').replace('3','C').replace('4','G')[::-1])
													
					if "*" in stop_dna:
						stop_code=stop_dna[1:stop_dna.find('*')] 
						last_exon_boundary=start_exon-(stop_dna.find('*'))*3
						line2[7]=line2[7].replace(line2[7].split(" ")[0], str(last_exon_boundary),1)
						line2[1]=str(last_exon_boundary)
						modify_end=1
				else:
					print "Partial", line2[0]	

##############################first_exon_boundary is for the first exon (that is supposed to start with M) while last_exon_boundary is for the last exon (which is supposed to contain *)################

		for k1, v1 in final_hit_seq.iteritems():	
			x=str(scaf)+"_"+str(start_exon)+"-"+str(end_exon)
			x_reverse=str(scaf)+"_"+str(end_exon)+"-"+str(start_exon)	########## Added on 26thJune2020 for cases where the exon boundary is reversed
			if ((x in k1) or (x_reverse in k1)):
				header=k1.split("_")
				header_start_end=k1.split("-")[0].split("_")[-1]+"-"+k1.split("-")[1].split("_")[0]
				temp=header.index(header_start_end)
				header[0:temp]=['_'.join(header[0:temp])]
				if line2[3]=="forward" and int(line2[4])>=350:
					if modify_start==1:
						header[1]=header[1].replace(header[1].split("-")[0], str(first_exon_boundary))
						final_hit_seq[k1]=start_code+v1
						header[5]="StartCodonPresent"

					if v1.find("M")!=-1 and v1.find("M")<10 and v1[0]!="M":	
						length_of_trim=v1.find("M") +1
						if modify_start!=1 and length_of_trim<11 and first_exon_length>30:
							final_hit_seq[k1]=v1[length_of_trim-1:]
							first_exon_boundary=start_exon+(length_of_trim*3)
							line2[7]=line2[7].replace(line2[7].split(" ")[0], str(first_exon_boundary),1)
							line2[1]=str(first_exon_boundary)
							header[1]=header[1].replace(header[1].split("-")[0], str(first_exon_boundary))
							header[5]="StartCodonPresent"


					if modify_end==1:	
						header[1]=header[1].replace(header[1].split("-")[1], str(last_exon_boundary))
						final_hit_seq[k1]+=stop_code

				else:
					if modify_start==1:
						header[1]=header[1].replace(header[1].split("-")[1], str(first_exon_boundary))
						final_hit_seq[k1]=start_code+v1
						header[5]="StartCodonPresent"

					if v1.find("M")!=-1 and v1.find("M")<10 and int(line2[4])>=350 and v1[0]!="M": 
						length_of_trim=v1.find("M") +1
						if modify_start!=1 and length_of_trim<11 and last_exon_length>30:
							final_hit_seq[k1]=v1[length_of_trim-1:]
							first_exon_boundary=end_exon-(length_of_trim*3)
							line2[7]=line2[7].replace(line2[7].split(" ")[-1], str(first_exon_boundary),1)
							line2[2]=str(first_exon_boundary)
							header[1]=header[1].replace(header[1].split("-")[1], str(first_exon_boundary))
							header[5]="StartCodonPresent"


					if modify_end==1:
						header[1]=header[1].replace(header[1].split("-")[0], str(last_exon_boundary))				
						final_hit_seq[k1]+=stop_code
				
				length=len(str(final_hit_seq[k1]))
				line2[4]=length
				if (length>=int(sys.argv[5])):
					header[3]="Complete"				
					line2[5]="Complete"
				#print "hello", final_hit_seq[k1], modify_start, modify_end, header[1], line2[4]

				#header=str(header)
				#print "header", header

				f_protein.writelines("%s_%s_%s_%s_%s_%s\n%s\n" %(header[0], header[1], header[2], header[3], header[4], header[5], final_hit_seq[k1]))	
				final_hit_seq[k1]=final_hit_seq[k1].upper().replace('*', '').replace('B', '').replace('J', '').replace('O', '').replace('Z', '').replace('U', '').replace('X','')
				f_starminus_protein.writelines("%s_%s_%s_%s_%s_%s\n%s\n" %(header[0], header[1], header[2], header[3], header[4], header[5], final_hit_seq[k1]))

		#print line2[0]
		f_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (line2[0], line2[1], line2[2], line2[3], line2[4], line2[5], line2[6], line2[7], line2[8]))	
		if line2[3]=="forward":
			strand="+"
		else:
			strand="-"

		f_gff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s_%s_%s_%s_%s_%s\n" % (line2[0], "GWS", "gene", line2[1], line2[2], ".", strand, ".", "ID=gene_"+line2[0], header[1], header[2], header[3], header[4], header[5],))		
		f_gff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s_%s_%s_%s_%s_%s;%s\n" % (line2[0], "GWS", "mRNA", line2[1], line2[2], ".", strand, ".", "ID=mRNA_"+line2[0], header[1], header[2], header[3], header[4], header[5], "Parent=gene_"+line2[0]+"_"+header[1]+"_"+header[2]+"_"+header[3]+"_"+header[4]+"_"+header[5],))	

		each_exon=line2[7].split(" ")
		#print len(each_exon)
		for i in range (0, len(each_exon), 2):
			s=each_exon[i]
			e=each_exon[i+1]
			f_gff.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s_%s_%s_%s_%s_%s;%s\n" % (line2[0], "GWS", "exon", s, e, ".", strand, ".", "ID=exon_"+line2[0], header[1], header[2], header[3], header[4], header[5], "Parent=mRNA_"+line2[0]+"_"+header[1]+"_"+header[2]+"_"+header[3]+"_"+header[4]+"_"+header[5],))		


f_table.close()
f_protein.close()
f_gff.close()
f_starminus_protein.close()

