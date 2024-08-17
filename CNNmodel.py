import datetime,sys,os
import InPre as ip
import numpy as np
import matplotlib.pyplot as plt

import tensorflow as tf
from keras import backend as K
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Embedding, Input, Reshape
from keras.layers import Conv1D, MaxPooling1D, AveragePooling2D
from keras.models import load_model

from sklearn.metrics import accuracy_score, auc, roc_curve, confusion_matrix, f1_score
from scipy import interp
from itertools import cycle


x_all = ip.x_get_all()

pic_dir = os.getcwd()+"/figures/"
if not os.path.exists(pic_dir):
	os.mkdir(pic_dir)

n_length = x_all.shape[1]

test_row = [i for i in range(1,16)]

f_metrics = open("cnn_metrics.txt","w")
for p in test_row:

	f_metrics.write("row_"+str(p))
	y_all,y_catego,genome_class = ip.y_value_get_one_hot(p,res_sus=False)
	n_outputs = len(y_catego)
	y_labels = y_catego

	x_train, x_test, y_train, y_test, train_class\
	 = ip.train_split(x_all, y_all, genome_class, p_train=0.8, train_seed=1, upper=False)
		
	starttime = datetime.datetime.now()

	if not os.path.exists("cnn_model_"+str(p)+".h5"):

		inputs = Input(shape=(n_length,4))
		x = Conv1D(8, 3, activation='relu')(inputs)
		x = Conv1D(8, 3, activation='relu')(x)
		x = MaxPooling1D(3,3)(x)
		x = Dropout(0.25)(x)
		x = Conv1D(16, 3, activation='relu')(x)
		x = Conv1D(16, 3, activation='relu')(x)
		x = MaxPooling1D(3,3)(x)
		x = Dropout(0.25)(x)
		x = Conv1D(32, 3, activation='relu')(x)
		x = Conv1D(32, 3, activation='relu')(x)
		x = MaxPooling1D(3,3)(x)
		x = Dropout(0.25)(x)
		x = Flatten()(x)
		x = Dense(512,activation='relu')(x)
		x = Dropout(0.5)(x)
		x = Dense(64,activation='relu')(x)
		x = Dropout(0.5)(x)
		preditions = Dense(n_outputs, activation='softmax')(x)

		model = Model(inputs=inputs,outputs=preditions)

		print(model.summary())

		model.compile(loss='categorical_crossentropy',
					  optimizer='adam',
					  metrics=['accuracy'])

		model.fit(x_train, y_train, verbose=2, batch_size=10, epochs=20)
		model.save("cnn_model_"+str(p)+".h5")

	else:
		model = load_model("cnn_model_"+str(p)+".h5")

	print("predicting...")
	y_pred_prob = model.predict(x_test,batch_size=10)

	K.clear_session()
	tf.reset_default_graph()

	y_pred = np.argmax(y_pred_prob,axis=1)
	y_true = np.argmax(y_test,axis=1)

	# for prob, pred, true in zip(y_pred_prob,y_pred,y_true):
	# 	print(prob,pred,true)

	xlocations = np.array(range(n_outputs))
	matrix = confusion_matrix(y_true=y_true,y_pred=y_pred,labels=xlocations)

	accuracy = accuracy_score(y_true,y_pred)
	f1 = f1_score(y_true,y_pred,average="macro")
	f_metrics.write(","+str(accuracy)+","+str(f1))
	endtime = datetime.datetime.now()
	print("accuracy:",accuracy)
	print("f1_score:",f1)
	print("run time:",(endtime-starttime).seconds,"\n")

	fpr = dict()
	tpr = dict()
	roc_auc = dict()
	for i in range(len(y_catego)):
		fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_pred_prob[:, i])
		roc_auc[i] = auc(fpr[i], tpr[i])

	fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_pred_prob.ravel())
	roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

	all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(y_catego))]))

	# Then interpolate all ROC curves at this points
	mean_tpr = np.zeros_like(all_fpr)
	for i in range(len(y_catego)):
		mean_tpr += interp(all_fpr, fpr[i], tpr[i])

	# Finally average it and compute AUC
	mean_tpr /= len(y_catego)

	fpr["macro"] = all_fpr
	tpr["macro"] = mean_tpr
	roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

	# Plot all ROC curves
	plt.figure()
	lw = 1
	plt.plot(fpr["micro"], tpr["micro"],
			 label='micro-average ROC curve (area = {0:0.2f})'.format(roc_auc["micro"]),
			 color='deeppink', linestyle=':', linewidth=2)

	plt.plot(fpr["macro"], tpr["macro"],
			 label='macro-average ROC curve (area = {0:0.2f})'.format(roc_auc["macro"]),
			 color='navy', linestyle=':', linewidth=2)

	f_metrics.write(","+str(roc_auc["micro"]))
	f_metrics.write(","+str(roc_auc["macro"]))
	colors = cycle(['aqua', 'darkorange', 'cornflowerblue'])
	for i, color in zip(range(len(y_catego)), colors):
		f_metrics.write(","+str(roc_auc[i]))
		plt.plot(fpr[i], tpr[i], color=color, lw=lw,
				 label=('ROC curve of %s (area = {1:0.2f})'% y_catego[i]).format(i, roc_auc[i]))
	f_metrics.write("\n")

	plt.plot([0, 1], [0, 1], 'k--', lw=lw)
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.0])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right")
	plt.savefig(pic_dir+str(p)+"_cnn_roc.png",format="PNG",dpi=300,bbox_inches="tight")
	plt.close()
	
	fig, ax = plt.subplots()
	ax.matshow(matrix,cmap=plt.cm.gray_r)
	plt.xticks(xlocations, y_labels, rotation=40)
	plt.yticks(xlocations, y_labels)
	for i in range(len(matrix)):
		for j in range(len(matrix)):
			text = ax.text(j, i, matrix[i, j],ha="center", va="center", color="r")
	plt.savefig(pic_dir+str(p)+"_cnn_confmat.png",format="PNG",dpi=300,bbox_inches="tight")
	plt.close()

f_metrics.close()