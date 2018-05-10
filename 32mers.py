
import os
#directory=input("directory: ") 
#can also hard code this as directory="<path>", make sure ends in forward slash
directory="/Users/avahoffman/Documents/CSU/Research/python_misc/poop/"
for filename in os.listdir(directory):
	if filename.endswith(".txt"):
		print(os.path.join(directory, filename))
		with open(directory+"/"+filename) as open_file, open(directory+"/"+filename+"results.txt","w") as new_file:
			kmer_dic={}
			for line in open_file:
				if line.startswith(">"):
					eachgene=open_file.readline()
					charcount=0
					for char in eachgene:
						charend=charcount+32
						kmer=eachgene[charcount:charend]
						if len([c for c in kmer if c.isalpha()]) == 32:
							key=kmer
							if key in kmer_dic:
								kmer_dic[key]+=1
							else:
								kmer_dic[key]=1
						charcount+=1
			new_file.write("32mer"+"\t"+"count"+"\n")			
			for key in kmer_dic:
				new_file.write(str(key)+"\t"+str(kmer_dic[key])+"\n")
