#####################################
#		Amplitudes			#
#####################################


amps<-function(x){
  
  dens<-density(x)
  dd<-data.frame(dens$x,dens$y)
  
  maxima.ind<-which(diff(sign(diff(dd[,2])))==-2)+1
  minima.ind<-which(diff(sign(diff(dd[,2])))==2)+1
  colnames(dd)<-c("Location (x)","Amplitude (y)")
  peaks<-as.matrix(dd[maxima.ind,])
  antimode<-as.matrix(dd[minima.ind,])
  ls<-list(peaks, antimode)
  names(ls)<-c("Peaks", "Antimode")
  return(ls)
}

#####################################
#	Bimodality Coefficient		#
#####################################



bimodality_coefficient<-function(x, finite=TRUE,...){
  if(finite==TRUE){
    G=skewness(x,finite)
    sample.excess.kurtosis=kurtosis(x,finite)
    K=sample.excess.kurtosis
    n=length(x)
    B=((G^2)+1)/(K+ ((3*((n-1)^2))/((n-2)*(n-3))))
  }
  else{
    G=skewness(x,FALSE)
    K=kurtosis(x,FALSE)
    B=((G^2)+1)/(K)
  }
  return(B)
}


#####################################
#	 	    Modes			#
#####################################

modes<-function(data,type=1,digits="NULL",nmore="NULL"){
  if (type=="1"){
    tabs<-rle(sort(as.integer(data)))
    tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
    nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
                           tabs[2,][tabs[2,]==max(tabs[2,])]))
    
  }
  if (type=="2"){
    tabs<-rle(sort(signif(data,digits)))
    tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
    nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
                           tabs[2,][tabs[2,]==max(tabs[2,])]))
    
  }	
  if(!nmore=="NULL" & is.numeric(nmore) & type=="1"|type=="2"){
    for (i in 1:nmore){
      nmode<-cbind(nmode,as.matrix(rbind(tabs[1,][tabs[,2]==max(
        tabs[,2][tabs[,2]<nmode[2,i]])],tabs[2,][tabs[,2]==max(
          tabs[,2][tabs[,2]<nmode[2,i]])])))
      
    }
  }
  
  if (type=="3"){
    tabs<-rle(sort(as.vector(as.factor(data))))
    tabs<-t(as.matrix(cbind(tabs$values,tabs$lengths)))
    nmode<-as.matrix(rbind(tabs[1,][tabs[2,]==max(tabs[2,])],
                           tabs[2,][tabs[2,]==max(tabs[2,])]))
    
    if(!nmore=="NULL" & is.numeric(nmore) & type=="3"){
      
      j<-dim(nmode)[2]
      for (i in j+1:nmore){
        tmp<-nth_highest(tabs[2,],i)
        tmp2<-tabs[1,][tabs[2,]==tmp]
        nmode<-cbind(nmode,rbind(tmp2,tmp))
      }}}
  
  if (any(nmode[2,]==1)==TRUE){warning ("A single observation
	 is being observed as a mode.
	Double check the class or inspect the data.
	Alternatively, you may have specified 'nmore' too many times 
	for this data.")}
  
  rownames(nmode)<-c("Value","Length")
  nmode<-as.matrix(nmode[,!duplicated(nmode[1,])])
  return(nmode)
}

#####################################
#		  Skewness			#
#####################################




skewness<-function(x, finite=TRUE){
  n=length(x)
  S=(1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
  if(finite==FALSE){
    S=S
  }else{
    S=S*(sqrt(n*(n-1)))/(n-2)
  }
  return(S)	
}



#####################################
#		  Kurtosis			#
#####################################



kurtosis<-function(x, finite){
  n=length(x)
  K=(1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2) - 3
  if(finite==FALSE){
    K=K
  }
  else{
    K=((n-1)*((n+1)*K - 3*(n-1))/((n-2)*(n-3))) +3
  }
  return(K)	
}



#####################################
#		Amplitudes			#
#####################################


amps<-function(x){
  
  dens<-density(x)
  dd<-data.frame(dens$x,dens$y)
  
  maxima.ind<-which(diff(sign(diff(dd[,2])))==-2)+1
  minima.ind<-which(diff(sign(diff(dd[,2])))==2)+1
  colnames(dd)<-c("Location (x)","Amplitude (y)")
  peaks<-as.matrix(dd[maxima.ind,])
  antimode<-as.matrix(dd[minima.ind,])
  ls<-list(peaks, antimode)
  names(ls)<-c("Peaks", "Antimode")
  return(ls)
}



#####################################
#	N-th Highest Value		#
#####################################



nth_highest<-function(x,k=1){
  x<-as.vector(x)
  N<-length(x)
  sort(x,partial=N-k+1)[N-k+1]
}

##	Parametric Functions



#########################
#	Ashman's D		#
#########################


Ashmans_D<- function(mu1, mu2, sd1, sd2,...){
  D= sqrt(2)*(abs(mu1-mu2))/sqrt((sd1^2)+(sd1^2))
  return(D)
}


#####################################
#	Bimodality Separation		#
#####################################


bimodality_separation<-function(mu1,mu2,sd1,sd2,...){
  S=0.5*(mu1-mu2)/(sd1+sd2)
  return(S)
}
