#source("lib.R")
#source("midcor.R")
#correct("GluC2C4b","var")
case2<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2
    for(j in 1:nln) { k=1; for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
               for(k in 2:(nfrg+1)) for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k])*ff;#
    }
      fr<-mdistr(nfrg,tmp,mmlab,nln); return(fr)
}

corrm<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2;
     k=1; for(i in 1:(nmass+1-k)) tmp[i+1+k]<-tmp[i+1+k]-corr[i]*(fr[1+k]);
   for(k in 2:(nfrg+1)) for(i in 1:(nmass+1-k)) tmp[i+1+k]<-tmp[i+1+k]-corr[i]*(fr[1+k])*ff;#
       return(mdistr(nfrg,tmp,mmlab,1))
}

fitf<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2;
  dpred=numeric(nfrg); dpred[2]=1.; xisq=0.; modf=1.01; ff1=ff
    fr<-corrm(tmp,mmlab,corr,ff,fr,nfrg);   xisq0=sum((fr[2:(1+nfrg)]-dpred)**2)
   iw=0;   while(iw<2){ ff1=ff1*modf; print(paste("iw=",iw))
    fr<-corrm(tmp,mmlab,corr,ff,fr,nfrg);   xisq=sum((fr[2:(1+nfrg)]-dpred)**2)
                if(xisq<xisq0) {ff=ff1; }
                 else {modf=1./modf;  iw=iw+1}
      }
      
       return(list(ff,xisq,fr))
}

correct<-function(fname,injb=11,injf=12,plab=4,plaf=9){#fname is the name of file with raw data;
# injb,injf,plab,plaf are the positions of the initial and final characters in the row name
# designating the biological sample and conditions correspondingly
 fn<-file.path(fname);# read experimental data
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
  lst<-suminj(gcms,nln-1,coln,injb,injf); gcms=lst[[1]]; nln=lst[[2]];
# normalization
         gcmsn<-norm(gcms,ef);
# mass fractions
      fr<-mdistr(nfrg,gcmsn,mmlab,nln);
 write("\n*** Summed injections for each plate, corrected only for natural 13C, 29,30Si, 33,34S **",fn1,append=TRUE)
  write.table(format(fr,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# correction
           corr<-numeric(nmass); corr1=numeric(nmass); icomm=0;
           patt=grep('Cold',fr[,1]);
      if(length(patt)>0) corr<-gcmsn[patt[1],3:(2+nmass)]-mmlab[1,]; print(corr)
           patt=grep('natur',fr[,1]);
      if(length(patt)>0) corr1<-gcmsn[patt[1],3:(2+nmass)]-mmlab[1,];
           patt=grep('fitf',fr[,1]); abc=fr[patt[1],]; ff=1.
      if(length(patt)>0) {ff=0.63; abc=fitf(tmp=gcmsn[patt[1],],mmlab,corr,ff,abc,nfrg); ff=abc[[1]]; print(abc)}
       mdmax=max(corr1-corr);print(paste("max=",mdmax," ff=",ff))
       if((mdmax<0.005)&&(mdmax > -0.005)) md="va"
        else md="co";
  if(md=="va") { for(ii in 1:9) {tmp<-gcmsn;
      fr<-case2(gcmsn,mmlab,corr,ff,fr,nfrg);
        }}  else {tmp<-gcmsn; for(j in 1:nln)  tmp[j,3:(2+nmass)]<-tmp[j,3:(2+nmass)]-corr;
     fr<-mdistr(nfrg,tmp,mmlab,nln);}
# statistics in row data     
  lst=stat(tmp,nln,nfrg,plab,plaf); gcmsn2=lst[[1]]; gcmsn2$V2=NULL; 
# write("\n*** Summed injections, statistics for raw data corrected for peaks overlap **",fn1,append=TRUE)
#  write.table(format(gcmsn2,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
#     statistics in final distribution of 13C isotopomers
lst=stat(fr,nln,nfrg,plab,plaf); fr=lst[[1]]; nln=lst[[2]]; len=length(fr)-1; #CF changed lst=stat(fr,nln,nfrg,4,9); fr=lst[[1]]; nln=lst[[2]]; len=length(fr)-1;
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
# write("\n*** chi: **",fn1,append=TRUE)
#     write.table(format(gcms,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
 write("\n*** 99% CI: **",fn1,append=TRUE)
     
     mid=frmod[2,];  fac=0.98; fac1=0.9;     rada=gcmsn2[2,]
  #confin(mid,rada,fac,fac1,mmlab,nmass,nfrg,11,fn1)
  #   mid=frmod[3,];  rada=gcmsn2[3,]
  #confin(mid,rada,fac,fac1,mmlab,nmass,nfrg,11,fn1)
}
