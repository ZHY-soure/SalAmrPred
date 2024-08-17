import os,csv,sys,time,subprocess

# get all fasta files
directory = str(input("directory: "))
directory = directory.replace("\\","/")
if not directory.endswith("/"):
	directory += "/"

raw_fna_files = []
for filename in os.listdir(directory):
	if filename.endswith(".fna"):
		raw_fna_files.append(filename.split(".fna")[0])

blast_dir = ""
identity_n = 70

n_clus = 0
cluster_info = {}
cluster_db = {}
ref_cluster = {}

def makedb(in_name,db_name):
	makeblastdb = blast_dir+"makeblastdb"
	cmd = [makeblastdb,"-in",in_name,"-out",db_name,"-dbtype","nucl"]
	subprocess.check_output(cmd)
	#subprocess.check_output([makeblastdb+" -in "+in_name+" -out "+db_name+" -dbtype nucl"])

def os_blast(query_name,db_name,out_name,outfmt="6 std slen"):
	# the directory of blastn
	blastn = blast_dir+"blastn"
	cmd = [blastn,"-query",query_name,"-db",db_name,"-out",out_name,"-outfmt",outfmt,\
	"-evalue","1e-20","-perc_identity",str(identity_n),"-num_threads","36"]
	subprocess.check_output(cmd)

def xml_blast(query_name,db_name,out_name,outfmt="5"):
	# the directory of blastn
	blastn = blast_dir+"blastn"
	cmd = [blastn,"-query",query_name,"-db",db_name,"-out",out_name,"-outfmt",outfmt,"-num_threads","36"]
	subprocess.check_output(cmd)

def gather_rows(file):
	new_content = []
	with open(file,"r",encoding="utf-8") as f:
		content = list(f.readlines())
		seq_index = []
		for i in range(len(content)):
			try:
				if ">" == content[i][0]:
					seq_index.append(i)
			except:
				pass
		seq_index.append(len(content))
		for i in range(len(seq_index)-1):
			new_content.append(content[seq_index[i]])
			seq = ""
			line_start = seq_index[i]+1
			line_end = seq_index[i+1]
			for m in range(line_start,line_end):
				seq += content[m].strip("\n")
			new_content.append(seq)
	return new_content	

def find_seq(genome_content,gene):
	return_seq = ""
	for i in range(len(genome_content)):
		if ">"+gene == genome_content[i].strip("\n").split(" ")[0]:
			return_seq = genome_content[i+1]
	if return_seq == "":
			print(gene," not found")
	return return_seq

def self_cluster(file,n_cluster):

	cluster = {}
	genome_content = gather_rows(file)

	makedb(file,"self_db")
	os_blast(file,"self_db","self_cluster.xls")
	with open("self_cluster.xls","r",encoding="utf-8") as csvfile:
		plots = list(csv.reader(csvfile, delimiter="	"))
		for row in plots:
			
			if row[0] not in cluster and row[1] not in cluster:
				n_cluster += 1
				cluster[row[0]] = "cluster"+str(n_cluster)
				cluster[row[1]] = "cluster"+str(n_cluster)
				cluster_db["cluster"+str(n_cluster)] = set()
				if row[3] >= row[12]:
					ref_cluster["cluster"+str(n_cluster)] = row[0]
				else:
					ref_cluster["cluster"+str(n_cluster)] = row[1]
			elif row[0] in cluster and row[1] not in cluster:
				cluster[row[1]] = cluster[row[0]]
				if row[3] >= row[12]:
					ref_cluster[cluster[row[0]]] = row[0]
				else:
					ref_cluster[cluster[row[0]]] = row[1]
			elif row[0] not in cluster and row[1] in cluster:
				cluster[row[0]] = cluster[row[1]]
				if row[3] >= row[12]:
					ref_cluster[cluster[row[1]]] = row[0]
				else:
					ref_cluster[cluster[row[1]]] = row[1]

	print("total cluster number:",n_cluster)
	for gene in cluster:
		cluster_db[cluster[gene]].add(gene)
		if cluster[gene] not in cluster_info:
			cluster_info[cluster[gene]] = find_seq(genome_content,gene)

	with open("cluster.fasta","w",encoding="utf-8") as f:
		for clus in cluster_info:
			f.write(">"+clus+"\n"+cluster_info[clus]+"\n")
	return n_cluster

def update_cluster(file):

	genome_content = gather_rows(file)
	clusterd = []

	makedb("cluster.fasta","cluster")
	os_blast(file,"cluster","cluster.xls")
	with open("cluster.xls","r",encoding="utf-8") as csvfile:
		plots = list(csv.reader(csvfile, delimiter="	"))
		for row in plots:
			if float(row[2]) >= 70.0 and float(row[12])*1.3>=float(row[3])>=float(row[12])*0.7:
				cluster_db[row[1]].add(row[0])
				clusterd.append(">"+row[0])
				if row[3] > row[12]:
					ref_cluster[row[1]] = row[0]

	f_absent = open("absent.fasta","w",encoding="utf-8")
	for i in range(len(genome_content)):
		try:
			if ">" == genome_content[i][0] and genome_content[i].strip("\n").split(" ")[0] not in clusterd:
				f_absent.write(genome_content[i]+genome_content[i+1]+"\n")
		except:
			pass
	f_absent.close()

def process_print(list_info,item,interval=20):
	if (list_info.index(item)+1)%interval == 0 and (list_info.index(item)+1) != len(list_info):
		print(list_info.index(item)+1,"/",len(list_info))
	if (list_info.index(item)+1) == len(list_info):
		print(list_info.index(item)+1,"/",len(list_info))


# get all cluster
if not os.path.exists("cluster.cls"):
	n_clus = self_cluster(directory+raw_fna_files[0]+".fna",n_clus)
	for name in raw_fna_files[1:]:
		process_print(raw_fna_files,name,1)
		update_cluster(directory+name+".fna")
		try:
			n_clus = self_cluster("absent.fasta",n_clus)
		except:
			print("no absent gene")
			pass

	with open("cluster.cls","w",encoding="utf-8") as f:
		for clus in cluster_db:
			f.write(">"+clus+"\n")
			for gene in cluster_db[clus]:
				f.write(gene+"\n")

	with open("ref.cls","w",encoding="utf-8") as f:
		for clus in ref_cluster:
			f.write(clus+"	"+ref_cluster[clus]+"\n")
else:
	# get core SNPs
	with open("cluster.cls","r",encoding="utf-8") as f:
		line = " "
		while line:
			line = f.readline()
			try:
				if ">" == line[0]:
					clus = line.split("\n")[0][1:]
					cluster_db[clus] = set()
				elif len(line) != 0:
					cluster_db[clus].add(line.strip("\n"))
			except:
				pass

	with open("ref.cls","r",encoding="utf-8") as f:
		line = " "
		while line:
			line = f.readline()
			try:
				ref_cluster[line.split("	")[0]] = line.split("	")[1].strip("\n")
			except:
				pass

core_genes = {}
acc_genes = {}

def find_mid(line,start,end):
	if start in line and end in line:
		mid_string = line.split(start)[1].split(end)[0]
	return mid_string

gene_db = {}
for name in raw_fna_files:
	process_print(raw_fna_files,name,50)
	content = gather_rows(directory+name+".fna")
	for i in range(len(content)):
		try:
			if ">" == content[i][0]:
				gene_name = (content[i].strip("\n").split(" ")[0])[1:]
				gene_db[gene_name] = [content[i+1],name]
		except:
			pass

for clus in cluster_db:
	name_sum = set()
	for gene in cluster_db[clus]:
		name_sum.add(gene_db[gene][1])
	if len(name_sum) >= len(raw_fna_files)*0.99:
		core_genes[clus] = cluster_db[clus]
	else:
		acc_genes[clus] = cluster_db[clus]

print("numer of total accessory clusters:",len(acc_genes))
print("numer of total core clusters:",len(core_genes))

def SNP():
	snp_from = []
	snp_db = {}
	for name in raw_fna_files:
		snp_db[name] = []

	count = 0
	for clus in core_genes:

		has_genome = set()
		for gene in core_genes[clus]:
			has_genome.add(gene_db[gene][1])

		genome = raw_fna_files[:]
		count += 1
		print("extracting SNPs ",count,"/",len(core_genes))

		ref_name = ref_cluster[clus]
		ref_seq = gene_db[ref_name][0]

		f_core = open("core.fasta","w",encoding="utf-8")
		for gene_name in core_genes[clus]:
			f_core.write(">"+gene_name+"\n")
			f_core.write(gene_db[gene_name][0]+"\n")
		f_core.close()

		f_ref = open("ref.fasta","w",encoding="utf-8")
		f_ref.write(">"+ref_name+"\n")
		f_ref.write(ref_seq+"\n")
		f_ref.close()

		makedb("core.fasta","self_db")
		xml_blast("ref.fasta","self_db","coreSNP",outfmt="5")

		snp_cont = {}
		with open("coreSNP","r",encoding="utf-8") as f:
			gene_length = 0
			qseq = ""
			hseq = ""
			total_gap = 0
			content = list(f.readlines())
			for i in range(len(content)):
				if "<BlastOutput_query-len>" in content[i]:
					gene_length = int(find_mid(content[i],'<BlastOutput_query-len>',"</BlastOutput_query-len>"))
				if "<Hit_def>" in content[i]:
					gene_name = find_mid(content[i],"<Hit_def>","</Hit_def>").\
					split(" ")[0].replace("&gt;",">").replace("&apos;","'").replace("&amp;","&")
					geno = gene_db[gene_name][1]
					try:
						genome.remove(geno)
					except:
						print(geno," has been already removed")
						pass
					snp_cont[geno] = {}
				if "<Hsp_qseq>" in content[i]:
					qseq = find_mid(content[i],"<Hsp_qseq>","</Hsp_qseq>")
					hseq = find_mid(content[i+1],"<Hsp_hseq>","</Hsp_hseq>")
					total_gap = 0
				for i in range(len(qseq)):
					if qseq[i] != hseq[i]:
						if qseq != "-":
							snp_cont[geno][i-total_gap] = hseq[i]
						else:
							total_gap += 1
							x = 1
							while qseq[i-x] == "-":
								x += 1
							snp_cont[geno][(i-total_gap)+x*0.01] = hseq[i]
		
		all_snp_pos = []
		for name in snp_cont:
			for pos in snp_cont[name].keys():
				if pos not in all_snp_pos:
					all_snp_pos.append(pos)
		all_snp_pos = sorted(all_snp_pos)

		if all_snp_pos != []:
			for pos in all_snp_pos:
				snp_from.append(str(pos)+"	"+ref_name)

			for name in raw_fna_files:
				for pos in all_snp_pos:
					if name in snp_cont:
						if pos not in snp_cont[name]:
							try:
								snp_cont[name][pos] = ref_seq[pos]
							except:
								snp_cont[name][pos] = "-"
					else:
						try:
							snp_cont[name][pos] = "-"
						except:
							snp_cont[name] = {}
							snp_cont[name][pos] = "-"

			for name in snp_cont:
				for pos in all_snp_pos:
					snp_db[name].append(snp_cont[name][pos])

		print("number of SNPs:",len(snp_db[raw_fna_files[0]]),len(snp_db[raw_fna_files[11]]),"\n")

	letter = {"A":"1","C":"2","G":"3","T":"4"}
	f_snp = open("SNP.csv","w",encoding="utf-8")
	for name in snp_db:
		f_snp.write(name.split("_ar")[0])
		for snp in snp_db[name]:
			if snp in letter:
				f_snp.write(","+letter[snp])
			else:
				f_snp.write(",0")
		f_snp.write("\n")
	f_snp.close()

	f_snp_from = open("SNPfrom.txt","w",encoding="utf-8")
	for item in snp_from:
		f_snp_from.write(item+"\n")
	f_snp_from.close()

def ACC():
	acc_from = []
	acc_db = {}
	for name in raw_fna_files:
		acc_db[name] = []

	for clus in acc_genes:
		accg = list(acc_genes.keys())
		process_print(accg,clus,1000)
		acc_count = {}
		acc_from.append(list(acc_genes[clus])[0])

		for gene_name in acc_genes[clus]:
			geno = gene_db[gene_name][1]
			try:
				acc_count[geno] = 1
			except:
				acc_count[geno] = 1

		for name in raw_fna_files:
			try:
				acc_db[name].append(acc_count[name])
			except:
				acc_db[name].append(0)

	f_acc = open("ACC.csv","w",encoding="utf-8")
	for name in acc_db:
		f_acc.write(name.split("_ar")[0])
		for acc in acc_db[name]:
			f_acc.write(","+str(acc))
		f_acc.write("\n")
	f_acc.close()

	f_acc_from = open("ACCfrom.txt","w",encoding="utf-8")
	for item in acc_from:
		f_acc_from.write(item+"\n")
	f_acc_from.close()

if __name__ == '__main__':
	SNP()
	ACC()