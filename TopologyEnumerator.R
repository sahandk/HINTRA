##  FUNCTIONS  ##

BFS_check<-function(parvec){
  queue=c(0)
  flag=T
  for(index in 1:(length(parvec)+1)){
    if(index>length(queue)){#disconnected
      flag=F
      break
    }
    children=which(parvec==queue[index])
    if(length(children>0))
      if(length(intersect(children,queue))>0){#loop
        flag=F
        break
      }else{
        queue=c(queue,children)
      }
  }
  
  flag
}

Tree_Prior<-function(parvec){
  degs = tabulate(parvec+1)
  if(length(degs)< length(parvec)+1)
    degs=c(degs,rep(0,length(parvec)+1-length(degs)))
  
  lnum=0
  denum=0
  for(d in degs){
    lnum=lnum+lgamma(d+0.1)
    denum=denum+d+0.1
  }
  
  return(lnum - lgamma(denum))
}


##  MAIN CODE  ##

K=3 # The number of mutations

results=matrix(0,nrow = (K+1)^(K-1),ncol = K)
pv=rep(0,K)
counter=0
while(T){
  print(t(pv))
  #check
  if(BFS_check(pv)){
    print("accepted")
    counter=counter+1
    results[counter,]=pv
  }else
    print("rejected")
  
  #next
  flag=F
  for(i in K:1)
    if(pv[i]<K){
      flag=T
      break
    }
  if(flag){
    pv[i]=pv[i]+1
    if(i < K)
      for(j in (i+1):K)
        pv[j]=0
  }else
    break
}

priors=apply(results,1,Tree_Prior)
ord=order(priors,decreasing = T)

results=results[ord,]
nozeros=apply(results,1,function(x) length(which(x==0)))

write.table(results[nozeros==1,],paste0("./ParentVectors_K-",K,".txt"),quote = F,col.names = F,row.names = F,sep = "\t")
