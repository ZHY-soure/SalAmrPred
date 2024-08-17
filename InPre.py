import numpy as np
import os,csv,sys

directory = str(input("directory: "))
directory = directory.replace("\\","/")
if not directory.endswith("/"):
	directory += "/"
	
try:
	raw_fna_files = []
	for filename in os.listdir(directory):
		if filename.endswith(".fna"):
			raw_fna_files.append(filename.split(".fna")[0])
except:
	pass

raw_fna_files = sorted(raw_fna_files,key = lambda x:float(x))

def x_get_snp(one_hot=True):

	print("loading x dataset...")
	if not os.path.exists("snp.npy"):

		n_col = 0
		with open("SNP.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			n_col += len(plots[0])-1
		print(n_col)

		if one_hot:
			X_train = np.zeros((len(raw_fna_files),n_col,4),dtype="int8")
			letter = {"0":[0,0,0,0],"1":[1,0,0,0],"2":[0,1,0,0],"3":[0,0,1,0],"4":[0,0,0,1],"5":[1,1,1,1]}
		else:
			X_train = np.zeros((len(raw_fna_files),n_col),dtype="int8")
			letter = {"0":0,"1":1,"2":2,"3":3,"4":4,"5":5}

		with open("SNP.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			for m in range(len(plots)):
				if m%100 == 0:
					print(m)
				pos = raw_fna_files.index(plots[m][0])
				for n in range(len(plots[m])-1):
					X_train[pos][n] = letter[plots[m][n+1]]

		np.save("snp.npy",X_train)
	else:
		X_train = np.load("snp.npy")
	print("shape of x:",X_train.shape)
	return X_train

def x_get_acc(one_hot=True):

	print("loading x dataset...")
	if not os.path.exists("acc.npy"):

		n_col = 0
		with open("ACC.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			n_col += len(plots[0])-1
		print(n_col)

		if one_hot:
			X_train = np.zeros((len(raw_fna_files),n_col,4),dtype="int8")
			accer = {"0":[0,0,0,0],"1":[1,1,1,1]}
		else:
			X_train = np.zeros((len(raw_fna_files),n_col),dtype="int8")
			accer = {"0":0,"1":5}

		with open("ACC.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			for m in range(len(plots)):
				if m%100 == 0:
					print(m)
				pos = raw_fna_files.index(plots[m][0])
				for n in range(len(plots[m])-1):
					X_train[pos][n] = accer[plots[m][n+1]]

		np.save("acc.npy",X_train)
	else:
		X_train = np.load("acc.npy")

	print("shape of x:",X_train.shape)
	return X_train

def x_get_all(one_hot=True):

	print("loading x dataset...")
	if not os.path.exists("all.npy"):

		n_col1 = 0
		with open("SNP.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			n_col1 += len(plots[0])-1

		n_col2 = 0
		with open("ACC.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			n_col2 += len(plots[0])-1
		print(n_col1+n_col2)

		if one_hot:
			X_train = np.zeros((len(raw_fna_files),n_col1+n_col2,4),dtype="int8")

		if one_hot:
			letter = {"0":[0,0,0,0],"1":[1,0,0,0],"2":[0,1,0,0],"3":[0,0,1,0],"4":[0,0,0,1],"5":[1,1,1,1]}
			accer = {"0":[0,0,0,0],"1":[1,1,1,1]}
		else:
			letter = {"0":0,"1":1,"2":2,"3":3,"4":4,"5":5}
			accer = {"0":0,"1":5}

		with open("SNP.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			for m in range(len(plots)):
				if m%100 == 0:
					print(m)
				pos = raw_fna_files.index(plots[m][0])
				for n in range(len(plots[m])-1):
					X_train[pos][n] = letter[plots[m][n+1]]

		with open("ACC.csv","r") as csvfile:
			plots = list(csv.reader(csvfile, delimiter=","))
			for m in range(len(plots)):
				if m%100 == 0:
					print(m)
				pos = raw_fna_files.index(plots[m][0])
				for n in range(len(plots[m])-1):
					X_train[pos][n+n_col1] = accer[plots[m][n+1]]

		np.save("all.npy",X_train)
	else:
		X_train = np.load("all.npy")

	print("shape of x:",X_train.shape)
	return X_train

def x_get_seq():

	def process_print(list_info,item,interval=20):
		if (list_info.index(item)+1)%interval == 0 and (list_info.index(item)+1) != len(list_info):
			print(list_info.index(item)+1,"/",len(list_info))
		if (list_info.index(item)+1) == len(list_info):
			print(list_info.index(item)+1,"/",len(list_info))

	def gather_rows(file):
		new_content = []
		with open(file,"r") as f:
			content = list(f.readlines())
			seq_index = [i for i in range(len(content)) if ">" in content[i]]
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

	cluster_db = {}
	max_clus = 0
	with open("cluster.cls","r") as f:
		print("reading cluster database...")
		line = " "
		while line:
			try:
				line = f.readline()
			except:
				break
			if len(line) != 0:
				if ">" == line[0]:
					clus = int(line.strip("\n").split("cluster")[1])
					max_clus += 1
				else:
					cluster_db[line.strip("\n")] = int(clus)
	print("number of cluster:",max_clus)

	print("loading x dataset...")
	if not os.path.exists("snp.npy"):

		max_len = 0
		letter = {"A":[1,0,0,0],"C":[0,1,0,0],"G":[0,0,1,0],"T":[0,0,0,1]}

		for name in raw_fna_files:
			process_print(raw_fna_files,name,100)
			content = gather_rows(directory+name+"_ar.fna")
			for i in range(len(content)):
				if ">" == content[i][0]:
					max_len = max(max_len,len(content[i+1]))
		print("total length of raw seq:",max_len)

		X_train = list(np.zeros((len(raw_fna_files),max_clus,max_len,4),dtype="int8"))
		for name in raw_fna_files:
			index = raw_fna_files.index(name)
			process_print(raw_fna_files,name,100)
			content = gather_rows(directory+name+"_ar.fna")
			for i in range(len(content)):
				if ">" == content[i][0]:
					now_clus = cluster_db[content[i].split(" [")[0][1:]]
					for x in range(len(content[i+1])):
						try:
							X_train[index][now_clus][x] = letter[content[i+1][x].upper()]
						except:
							print(content[i+1][x])


		X_train = np.array(X_train,dtype="int8")
		np.save("snp.npy",X_train)
	else:
		X_train = np.load("snp.npy")
	print("shape of x:",X_train.shape)
	return X_train

def y_value_get_one_hot(row_num, res_sus=False, no_int=True):

	genome_class = {}
	phenotype = ["susce","inter","resis"]
	value = {}
	drug = ""
	res = 0
	sus = 0
	y_catego = []
	with open(directory+"matrix.csv","r") as csvfile:
		plots = list(csv.reader(csvfile, delimiter=","))
		drug = plots[0][row_num].replace("/","_")
		drug_mic = {"amoxicillin/clavulanic acid":[8,16,32],
			"ampicillin":[8,16,32],
			"azithromycin":[8,-1,16],
			"cefoxitin":[8,16,32], #tou bao xi ding
			"ceftiofur":[4,-1,8], #tou bao sai fu
			"ceftriaxone":[1,2,4], #tou bao qu song
			"chloramphenicol":[8,16,32],
			"ciprofloxacin":[0.0625,-1,1],
			"gentamicin":[4,8,16],
			"kanamycin":[16,32,64],
			"nalidixic acid":[16,-1,32],
			"streptomycin":[16,-1,32],
			"sulfisoxazole":[64,128,256],
			"tetracycline":[4,8,16],
			"trimethoprim/sulfamethoxazole":[2,-1,4]
			}
		sus = drug_mic[plots[0][row_num]][0]
		res = drug_mic[plots[0][row_num]][-1]
		print(drug+":")
		for row in plots[1:]:

			if row[row_num] != "-1":

				if not res_sus:
					try:
						value[row[0]] = float(row[row_num])

						try:
							genome_class[float(row[row_num])].append(raw_fna_files.index(row[0]))
						except:
							genome_class[float(row[row_num])] = [raw_fna_files.index(row[0])]

						if float(row[row_num]) not in y_catego:
							y_catego.append(float(row[row_num]))
					except:
						pass

					y_catego = sorted(y_catego,key = lambda x:float(x))

				elif res_sus:

					if not no_int:
						if row[row_num] not in phenotype:
							if float(row[row_num]) >= res:
								value[row[0]] = 2
								try:
									genome_class["resis"].append(raw_fna_files.index(row[0]))
								except:
									genome_class["resis"] = [raw_fna_files.index(row[0])]
							elif float(row[row_num]) <= sus:
								value[row[0]] = 0
								try:
									genome_class["susce"].append(raw_fna_files.index(row[0]))
								except:
									genome_class["susce"] = [raw_fna_files.index(row[0])]
							else:
								value[row[0]] = 1
								try:
									genome_class["inter"].append(raw_fna_files.index(row[0]))
								except:
									genome_class["inter"] = [raw_fna_files.index(row[0])]

						elif row[row_num] in phenotype:

							value[row[0]] = phenotype.index(row[row_num])
							try:
								genome_class[row[row_num]].append(raw_fna_files.index(row[0]))
							except:
								genome_class[row[row_num]] = [raw_fna_files.index(row[0])]

						y_catego = [0,1,2]

					else:
						if row[row_num] not in phenotype:
							if float(row[row_num]) >= res:
								value[row[0]] = 1
								try:
									genome_class["resis"].append(raw_fna_files.index(row[0]))
								except:
									genome_class["resis"] = [raw_fna_files.index(row[0])]
							elif float(row[row_num]) <= sus:
								value[row[0]] = 0
								try:
									genome_class["susce"].append(raw_fna_files.index(row[0]))
								except:
									genome_class["susce"] = [raw_fna_files.index(row[0])]

						elif row[row_num] in phenotype:

							if row[row_num] == "susce":
								value[row[0]] = 0
							elif row[row_num] == "resis":
								value[row[0]] = 1

							if row[row_num] != "inter":
								try:
									genome_class[row[row_num]].append(raw_fna_files.index(row[0]))
								except:
									genome_class[row[row_num]] = [raw_fna_files.index(row[0])]

						y_catego = [0,1]

	y = []
	count = 0
	for name in raw_fna_files:
		append = [0 for i in range(len(y_catego))]
		try:
			append[y_catego.index(value[name])] += 1
			count += 1
			y.append(append)
		except:
			y.append(append)
			pass
	print("add:",count)
	y = np.array(y)
	print("shape of label of each genome: ",y.shape)
	return y,y_catego,genome_class



"""
functions above: data matrice generation

functions below: train and test data split, and other graph functions

"""

def train_split(x_all, y_all, genome_class, p_train=0.8, train_seed=1,\
 row_num=1, replace_status=False, upper=False, split=False):

	ranges = {}
	if bool(split):
		if split != -1:
			with open("xgb_mic_range.txt","r") as f:
				content = list(f.readlines())
				row = 0
				for line in content:
					if "	" in line:
						row = line.split("	")[0]
						ranges[row] = line.split("	")[1].split(",")
				ranges = [int(each) for each in ranges[str(row_num)]]
			x_all = x_all[:,ranges[:split],:]
		else:
			x_all = x_all[:,[0],:]

	np.random.seed(train_seed)

	all_train = []
	all_test = []

	train_class = {}

	print("original:")
	for each in genome_class:
		print(each,len(genome_class[each]))

	print("class - train - test:")

	for each in genome_class:
		n_train = int(len(genome_class[each])*p_train)
		n_test = len(genome_class[each])-n_train

		this_train = list(np.random.choice(genome_class[each],n_train,replace=replace_status))
		this_test = [item for item in genome_class[each] if item not in this_train]

		all_train += this_train
		all_test += this_test

		train_class[each] = this_train

		print(each,"-",n_train,"-",n_test)

	all_train = np.array(all_train)
	all_train = all_train[np.random.permutation(len(all_train))]
	all_test = np.array(all_test)

	print("all train number: ",len(all_train))
	print("all test number: ",len(all_test))

	x_train = x_all[all_train]
	y_train = y_all[all_train]
	x_test = x_all[all_test]
	y_test = y_all[all_test]

	return x_train, x_test, y_train, y_test, train_class

def train_split_res(x_all, y_all, genome_class, row_num=1, p_train=0.8, train_seed=1, replace_status=False, upper=True, split=False):

	n_susce_train = 0
	n_susce_test = 0
	n_resis_train = 0
	n_resis_test = 0

	print("original:")
	for each in genome_class:
		print(each,len(genome_class[each]))

	print("class-train-test:")

	n_susce_train = int(len(genome_class["susce"])*p_train)
	n_susce_test = len(genome_class["susce"])-n_susce_train
	n_resis_train = int(len(genome_class["resis"])*p_train)
	n_resis_test = len(genome_class["resis"])-n_resis_train

	np.random.seed(train_seed)

	sus_train = list(np.random.choice(genome_class["susce"],n_susce_train,replace=replace_status))
	res_train = list(np.random.choice(genome_class["resis"],n_resis_train,replace=replace_status))

	if upper:
		n_increase = max(n_susce_train,n_resis_train)/(2*min(n_susce_train,n_resis_train))

		if n_increase > 1:
			if n_susce_train >= n_resis_train:
				n_resis_train = int(n_resis_train*n_increase)
				res_train = list(np.random.choice(res_train,n_resis_train,replace=True))
				

	print("susce-"+str(n_susce_train)+"-"+str(n_susce_test))
	print("resis-"+str(n_resis_train)+"-"+str(n_resis_test))

	all_genome = genome_class["susce"]+genome_class["resis"]
	random_train = np.array(sus_train+res_train)
	remain = np.array([item for item in all_genome if item not in random_train])

	random_train = random_train[np.random.permutation(len(random_train))]

	print("all train number: ",len(random_train))
	print("all test number: ",len(remain))

	x_train = x_all[random_train]
	y_train = y_all[random_train]
	x_test = x_all[remain]
	y_test = y_all[remain]


	return x_train, x_test, y_train, y_test, sus_train, res_train

def bagging_generate(x_train, y_train, bagging_seed=1, bagging_num=100):

	np.random.seed(bagging_seed)
	bagging = np.random.choice(range(len(x_train)),bagging_num,replace=True)

	x_train = x_train[bagging]
	y_train = y_train[bagging]
	
	return x_train, y_train

def bagging_generate_down(x_all, y_all, sus_train, res_train, bagging_seed=1, p_large=1.0):

	np.random.seed(bagging_seed)
	
	bagging = np.array(\
			list(np.random.choice(sus_train,min(100,len(sus_train)),replace=False))\
			+list(np.random.choice(res_train,min(100,len(res_train)),replace=False)))

	x_train = x_all[bagging]
	y_train = y_all[bagging]

	return x_train, y_train

def bagging_generate_down_int(x_all, y_all, sus_train, int_train, res_train, bagging_seed=1, p_large=2.5):

	np.random.seed(bagging_seed)

	n_sus = len(sus_train)
	n_int = len(int_train)
	n_res = len(res_train)

	n_sus_bagging = 0
	n_int_bagging = 0
	n_res_bagging = 0

	if n_sus >= int(n_int*p_large) and n_res >= int(n_int*p_large):
		if n_int >= 50:
			n_sus_bagging = int(n_int*p_large)
			n_int_bagging = n_int
			n_res_bagging = int(n_int*p_large)
		else:
			if n_sus > int(n_res*p_large):
				n_sus_bagging = int(n_res*p_large)
				n_int_bagging = n_int
				n_res_bagging = n_res
			elif n_res > int(n_sus*p_large):
				n_sus_bagging = n_sus
				n_int_bagging = n_int
				n_res_bagging = int(n_sus*p_large)
			else:
				n_sus_bagging = n_sus
				n_int_bagging = n_int
				n_res_bagging = n_res

	elif n_int > n_res:
		n_sus_bagging = int(n_int*p_large)
		n_int_bagging = n_int
		n_res_bagging = n_res

	else:
		print("bad distribution")

	print("bagging-train:")
	print("susce-"+str(n_sus_bagging))
	print("inter-"+str(n_int_bagging))
	print("resis-"+str(n_res_bagging))
	
	bagging = np.array(\
			list(np.random.choice(sus_train,n_sus_bagging,replace=False))\
			+list(np.random.choice(int_train,n_int_bagging,replace=False))\
			+list(np.random.choice(res_train,n_res_bagging,replace=False)))

	bagging = bagging[np.random.permutation(len(bagging))]

	x_train = x_all[bagging]
	y_train = y_all[bagging]

	return x_train, y_train


def output_heatmap(model, query_layer_name, img, phenotype, dim=2):
	"""Get the heatmap for image.

	Args:
		   model: keras model.
		   query_layer_name: name of query layer in the model.
		   img: processed input image.

	Returns:
		   heatmap: heatmap.
	"""
	# predict the image class
	preds = model.predict(img,batch_size=10)
	# get the last conv layer
	query_layer_name = model.get_layer(query_layer_name)

	all_heatmap = []
	for m in range(0,len(img)):
		if (m+1)%10 == 0:
			print(m+1,"/",len(img))
		# find the class index
		index = np.argmax(preds[m])
		# This is the entry in the prediction vector
		target_output = model.output[:, index]

		# compute the gradient of the output feature map with this target class
		grads = K.gradients(target_output, query_layer_name.output)[0]

		# mean the gradient over a specific feature map channel
		if dim == 2:
			pooled_grads = K.mean(grads, axis=(0, 1, 2))
		# this function returns the output of query_layer_name and grads 
		# given the input picture
		iterate = K.function([model.input], [pooled_grads, query_layer_name.output[0]])
		pooled_grads_value, conv_layer_output_value = iterate([img])

		# We multiply each channel in the feature map array
		# by "how important this channel is" with regard to the target class

		if dim == 2:
			for i in range(conv_layer_output_value.shape[-1]):
				conv_layer_output_value[:,:, i] *= pooled_grads_value[i]

		# The channel-wise mean of the resulting feature map
		# is our heatmap of class activation
		heatmap = np.mean(conv_layer_output_value, axis=-1)
		heatmap = np.maximum(heatmap,0)

		all_heatmap.append(heatmap)

	for heatmap in all_heatmap:
		try:
			total_map += heatmap
		except:
			total_map = heatmap

	if np.max(total_map) != 0:
		total_map /= np.max(total_map)

	return total_map

def weights_write(f_name, total_map):

	f_weights = open(f_name,"w")

	map_in_row = total_map.reshape(-1)
	n_max = len(map_in_row[map_in_row>0.1])
	max_indices = np.argsort(-map_in_row)[:n_max]
	max_from_clus = {}
	for index in max_indices:
		row = index//total_map.shape[1]
		col = index%total_map.shape[1]
		try:
			max_from_clus[row].append(col)
		except:
			max_from_clus[row] = [col]
	max_from = sorted(max_from_clus.items(),key=lambda x:len(x[1]),reverse=True)
	
	for each in max_from[:]:
		print("cluster:",each[0],"points num:",len(each[1]))
		f_weights.write("cluster:"+str(each[0])+" points num:"+str(len(each[1]))+"\n")

	f_weights.close()