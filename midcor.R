#source("lib.R")
#source("midcor.R")
#correct("GluC2C4b","var")
correct<-function(fname,mdcor="con"){
  md<-substr(mdcor,1,2);
# read experimental data
 fn<-file.path(fname);
  info<-read.table(fn, nrows=4);  nC<-info$V2[1]; nfrg<-info$V2[2]; nSi<-info$V2[3]; nS<-info$V2[4];
      gcms<-read.table(fn, skip=5);
           coln<-length(gcms);  nln<-nrow(gcms); nmass<-coln-2;
# theoretic distribution:
          mmlab<-mtr(nfrg,nmass,nC,nSi,nS);
# normalization
    ef<-sum(gcms[2])/sum(gcms[3]);
         gcmsn<-norm(gcms,ef);
# mass fractions
   fr<-mdistr(nfrg,gcmsn,mmlab,nln);# write mass fractions without correction:
 fn1<-paste(fn,"_c",sep="");
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# Sum the data for various injections from the same plate:
  lst=suminj(gcms,nln,coln); gcms=lst[[1]]; nln=lst[[2]];
# normalization
         gcmsn<-norm(gcms,ef);
# mass fractions
      fr<-mdistr(nfrg,gcmsn,mmlab,nln);
 write("\n*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **",fn1,append=TRUE)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# correction
           corr<-numeric(nmass);
       corr<-gcmsn[1,3:(2+nmass)]-mmlab[1,];
  if(md=="va") { for(ii in 1:9) {tmp<-gcmsn;  
#    for(j in 1:1)  for(k in 1:(nfrg+1))  for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
    for(j in 1:nln) { k=1; for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
               k=2; for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k])*0.573;#
               k=3; for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k])*0.;#
               for(k in 4:(nfrg+1)) for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k])*0.62;#
    }
      fr<-mdistr(nfrg,tmp,mmlab,nln);
        }}  else {tmp<-gcmsn; for(j in 1:nln)  tmp[j,3:(2+nmass)]<-tmp[j,3:(2+nmass)]-corr;
     fr<-mdistr(nfrg,tmp,mmlab,nln);}
# statistics in row data     
  lst=stat(tmp,nln,nfrg); gcmsn2=lst[[1]]; gcmsn2$V2=NULL; 
# write("\n*** Summed injections, statistics for raw data corrected for peaks overlap **",fn1,append=TRUE)
#  write.table(format(gcmsn2,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
#     statistics in final distribution of 13C isotopomers
lst=stat(fr,nln,nfrg); fr=lst[[1]]; nln=lst[[2]]; len=length(fr)-1;
# write data:
       razn<-gcmsn[1,];  razn[1]<-as.factor("correction"); for(i in 1:nmass) razn[i+1]<-corr[i];
       fr<-fr[1:len]
 write("\n*** Statistics, samples fully corrected **",fn1,append=TRUE)
    write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F); 
 write("*** Correction factor: **",fn1,append=TRUE)
     write.table(format(razn,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
   frmod=remsd(fr,"**sd**")
   gcmsn2=remsd(gcmsn2,"**sd**")
     gcms=alchi(frmod,gcmsn2,nmass,nfrg,mmlab);
 write("\n*** chi: **",fn1,append=TRUE)
     write.table(format(gcms,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
     mid=frmod[2,]; chi=0; nma=which(mid==max(mid[2:len])); fac=0.98;
     rada=gcmsn2[2,]
     for(i in 1:3){
   lst=ci99(mid,rada,mmlab,nma,fac,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,4,2,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,6,4,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,6,2,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,6,5,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,4,5,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   lst=ch2par(mid,5,2,0.001,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
     }
     razn=mid; razn[2:len]=mid[2:len]-frmod[2,2:len];
     cat(" diff: ",as.character(format(razn,digits=4)),"\n")
     
     cat("\n diff: ",as.character(format(razn,digits=4)),"\n",file=fn1,append=TRUE)
     cat(" mid: ",as.character(format(mid,digits=4)),"\n",file=fn1,append=TRUE)
}
