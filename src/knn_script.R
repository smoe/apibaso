
library(class)

k.min<-2
k.max<-7

# Takes a directory as first argument and the name of a
# domain as a second. It returns a vector describing
# the success with the three parameters
# 
# all      correctness of majority of k Neighbours
# certain  correctness of those cases that had 100% agreement
# coverage the number of sequences these reached

analyse.directory<-function(d,dom) {


	#
	# Reading in data
	#

	cat("Reading in main data set.\n")

	m.baso.baso<-read.table(paste(d,
		"/valid_basolateral.patterns_on_basolateral.sequences_",dom,
		"_main.distribution.tsv",sep=""))
	m.api.baso<-read.table(paste(d,
		"/valid_apical.patterns_on_basolateral.sequences_",dom,
		"_main.distribution.tsv",sep=""))
	m.api.api<-read.table(paste(d,
		"/valid_apical.patterns_on_apical.sequences_",dom,
		"_main.distribution.tsv",sep=""))
	m.baso.api<-read.table(paste(d,
		"/valid_basolateral.patterns_on_apical.sequences_",
		dom,"_main.distribution.tsv",sep=""))

	cat("Reading in test data set.\n")

	t.baso.baso<-read.table(paste(d,
		"/valid_basolateral.patterns_on_basolateral.sequences_",dom,
		"_test.distribution.tsv",sep=""))
	t.api.baso<-read.table(paste(d,
		"/valid_apical.patterns_on_basolateral.sequences_",dom,
		"_test.distribution.tsv",sep=""))
	t.api.api<-read.table(paste(d,
		"/valid_apical.patterns_on_apical.sequences_",dom,
		"_test.distribution.tsv",sep=""))
	t.baso.api<-read.table(paste(d,
		"/valid_basolateral.patterns_on_apical.sequences_",dom,
		"_test.distribution.tsv",sep=""))

	t.api<-cbind(t.api.api,t.api.baso)
	t.baso<-cbind(t.baso.api,t.baso.baso)
	t.t<-t(rbind(t.api,t.baso))


	cat("Assembling and preparing data matrix for main and test.\n")

	m.api<-cbind(m.api.api,m.api.baso)
	m.baso<-cbind(m.baso.api,m.baso.baso)
	m<-rbind(m.api,m.baso)
	m.t<-t(m)


	cat("Construction of classification vectors.\n")

	cl<-c(rep("a",ncol(m.api.api)),rep("b",ncol(m.baso.baso)))
	t.cl<-c(rep("a",ncol(t.api.api)),rep("b",ncol(t.baso.baso)))

	if (FALSE) {

		# INTERNAL EVALUATION

		maxindex<-71
		m.t.selection<-sample(1:maxindex,maxindex*9/10,replace=FALSE)
		m.t.train<-m.t[m.t.selection,]
		m.t.test<-m.t[(1:maxindex)[-m.t.selection],]        

		#cl.train<-cl[m.t.selection]
		cl.test<-cl[(1:maxindex)[-m.t.selection]]
		cl.submit<-factor(c(cl.train,cl.test))

		m.knn<-knn(m.t.train,m.t.test,cl.train,k=3,prob=TRUE)  

		r<-sapply(m.knn,function(X){levels(m.knn)[X];})
		r.eval<-sum(r==cl.test)*100/length(cl.test)

		r.pos<-attr(m.knn,"prob")==1
		r.pos.r<-sapply(m.knn[r.pos],function(X){levels(m.knn)[X];})
		r.eval.r<-sum(r.pos.r==cl.test[r.pos])*100/length(cl.test[r.pos])


		t.knn<-knn(m.t.train,t.t,cl.train,k=3,prob=TRUE) 
		t.r<-sapply(t.knn,function(X){levels(t.knn)[X];})
		t.r.eval<-sum(t.r==t.cl)*100/length(t.cl)

		t.r.pos<-attr(t.knn,"prob")==1
		t.r.pos.r<-sapply(t.knn[t.r.pos],function(X){levels(t.knn)[X];})
		t.r.eval.r<-sum(t.r.pos.r==t.cl[t.r.pos])*100/length(t.cl[t.r.pos])

	}

	cat("Performing validation at different values for ...\n")

	r<-rep(-1,(1+(k.max-k.min))*3)
	dim(r)<-c(3,(1+(k.max-k.min)))
	rownames(r)<-c("all","certain","coverage")
	colnames(r)<-paste("k=",(k.min:k.max),sep="")

	for (my.k in k.min:k.max) {

		all.knn<-knn(m.t,t.t,cl,k=my.k,prob=TRUE)
		all.r<-sapply(all.knn,function(X){levels(all.knn)[X];})
		all.r.eval<-sum(all.r==t.cl)*100/length(t.cl)

		cat("k=",my.k," ::: ",all.r.eval,"% correct on all ",
			length(t.cl)," test cases.\n",sep="")

		r["all",k.min-my.k+1]<-all.r.eval

		#significance.levels<-c(07.75,0.9,0.95,1)
		significance.levels<-c(1)
		for (significance in significance.levels) {
			all.r.pos<-attr(all.knn,"prob")>=significance
			all.r.pos.r<-sapply(all.knn[all.r.pos],function(X){levels(all.knn)[X];})
			all.r.eval.r<-sum(all.r.pos.r==t.cl[all.r.pos])*100/length(t.cl[all.r.pos])
			cat("k=",my.k," ::: ",all.r.eval.r,"% correct and ",
				significance*100,"% certain on ",sum(all.r.pos),"\n")
			r["certain",k.min-my.k+1]<-all.r.eval.r
			r["coverage",k.min-my.k+1]<-sum(all.r.pos)
		}
		cat("\n")
	}
	return(r)
}

data.sets<-Sys.glob("08_Datas_Ramrath/data_set_*")
domains<-c("complete","in","tm","out")

if (0 == length(data.sets)) stop("Did not find any files to analyse.")

analyse.directory.set<-function(ds,dom) {
	r.of.dir<-sapply(data.sets,analyse.directory,dom)
	dim(r.of.dir)<-c(3,length(data.sets),k.max-k.min+1)
	dimnames(r.of.dir) <- list(
		c("all","certain","coverage"),
		data.sets,
		paste("k=",(k.min:k.max),sep="")
	)
	return(r.of.dir)
}

results<-list()
for(dom in domains) {
	cat("\n\nDomain: ",dom,"\n\n")
	results[[dom]]<-analyse.directory.set(data.sets,dom)
}

for(dom in domains) {
	#fname<-paste("~/boxplots_",dom,".png",sep="")
	fname<-paste("~/boxplots_",dom,".pdf",sep="")
	pdf(fname)
	boxplot(as.data.frame(results[[dom]]["all",,c("k=2","k=3","k=4","k=5")]),
		#boxwex = 0.25, at = 1:3 - 0.2,
		col="yellow",main=paste("% correct in domain '",dom,"'")
			#,sub="% correct"
		)
	dev.off()
}


# retrieving data on coverage

for(dom in domains) {

	cat(paste("Investigating domain",dom,"\n"))

	m<-results[[dom]][c("certain","coverage"),,k=3]
	#print(m)
	v<-is.nan(m[1,])
	print(v)
	cat(paste("Nan:",sum(v),"  decent:",sum(!v),"\n"))
	print(m[,!v])
	stop
}
