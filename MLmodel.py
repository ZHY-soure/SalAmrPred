import datetime,os
import InPre as ip
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier

from sklearn.metrics import accuracy_score, auc, roc_curve, confusion_matrix, f1_score
from scipy import interp
from itertools import cycle
import pickle

name = [
		"RandomForestClassifier",
		"XGBClassifier",
		]

classifiers = [
				RandomForestClassifier(n_jobs=42,n_estimators=100),
				XGBClassifier(n_jobs=42,n_estimators=100),
				]

pic_dir = os.getcwd()+"/figures/"
if not os.path.exists(pic_dir):
	os.mkdir(pic_dir)

np.random.seed(1)
x_all = ip.x_get_all()
n_pos = x_all.shape[1] # number of posistions of x data
x_all = x_all.reshape((-1,n_pos*4)) # reshape into single dimension

test_row = [i for i in range(1,16)]

f_metrics = open("ml_metrics.txt","w")
f_weights = open("ml_feature_importance.txt","w")
f_range = open("ml_range.txt","w")

test_row = range(1,16)
for p in test_row:

	f_range.write(str(p)+"	")

	y_all,y_catego,genome_class = ip.y_value_get_one_hot(p,res_sus=False)
	n_outputs = len(y_catego)
	y_labels = y_catego

	x_train, x_test, y_train, y_test, train_class\
	 = ip.train_split(x_all, y_all, genome_class, p_train=0.8, train_seed=1, upper=False)

	for x in range(len(classifiers)):

		starttime = datetime.datetime.now()
		print(name[x]+":")
		f_metrics.write(name[x]+"_row_"+str(p))

		if not os.path.exists(name[x]+"_model_"+str(p)):
			clf = OneVsRestClassifier(classifiers[x])
			clf.fit(x_train,y_train)
			pickle.dump(clf, open(name[x]+"_model_"+str(p), 'wb'))
		else:
			clf = pickle.load(open(name[x]+"_model_"+str(p), 'rb'))

		for q in range(n_outputs):
			f_weights.write(name[x]+"_row_"+str(p)+"_class_"+str(y_catego[q])+":\n")
			try:
				importances = clf.estimators_[q].feature_importances_
				res_contribution = [[i,0] for i in range((len(importances)//4))]
				for i in range(len(importances)):
					res_contribution[i//4][1] += importances[i]
				res_contribution = sorted(res_contribution,key=lambda x:float(x[1]),reverse=True)
				for each in res_contribution[:500]:
					f_range.write(str(each[0])+",")
				f_range.write("-1\n")

				with open("allfrom.txt","r") as f:
					content = list(f.readlines())
					for i in range(100):
						f_weights.write(str(i)+"	"+str(res_contribution[i][1])+\
							"	"+str(res_contribution[i][0])+"	"+content[res_contribution[i][0]])
			except:
				print("wrong weights analysis of ",y_catego[q])

		print("predicting...")
		y_pred_prob = clf.predict_proba(x_test)
		y_pred = np.argmax(y_pred_prob,axis=1)
		y_true = np.argmax(y_test,axis=1)

		# for prob, pred, true in zip(y_pred_prob,y_pred,y_true):
		# 	print(prob,pred,true)

		xlocations = np.array(range(len(y_catego)))
		matrix = confusion_matrix(y_true=y_true,y_pred=y_pred,labels=xlocations)

		accuracy = accuracy_score(y_true,y_pred)
		f1 = f1_score(y_true,y_pred,average="macro")
		f_metrics.write(","+str(accuracy)+","+str(f1))
		endtime = datetime.datetime.now()
		print("accuracy:",accuracy)
		print("f1:",f1)
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

		plt.plot([0, 1], [0, 1], 'k--', lw=lw)
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.0])
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.title('Some extension of Receiver operating characteristic to multi-class')
		plt.legend(loc="lower right")
		plt.savefig(pic_dir+name[x]+"_"+str(p)+"_roc.png",format="PNG",dpi=300,bbox_inches="tight")
		plt.close()
		
		fig, ax = plt.subplots()
		ax.matshow(matrix,cmap=plt.cm.gray_r)
		plt.xticks(xlocations, y_labels, rotation=40)
		plt.yticks(xlocations, y_labels)
		for i in range(len(matrix)):
			for j in range(len(matrix)):
				text = ax.text(j, i, matrix[i, j],ha="center", va="center", color="r")
		plt.savefig(pic_dir+name[x]+"_"+str(p)+"_confmat.png",format="PNG",dpi=300,bbox_inches="tight")
		plt.close()

		f_metrics.write("\n")
		f_weights.write("\n")

f_metrics.close()
f_weights.close()
f_range.close()