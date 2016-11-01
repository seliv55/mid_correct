binom<-function(n){ #calculation of binomial cefficients
 bf<-numeric(n); bf[1]<-1;
 if(n>1){ for(i in 2:n) bf[i]<-bf[i-1]*i;  bc<-numeric(n); 
   for(i in 1:(n-1)) bc[i]=bf[n]/bf[i]/bf[n-i];}
   else  bc<-1;
     bc }

carbons<-function(n, bc) {#natural distribution of 13C according to the number of carbons in whole fragment
  pc12<-0.989; pc13<-0.011;
    mc<-numeric(n+1); mc[1]<-pc12^n;
 if(n>1) {for(i in 1:(n-1)) mc[i+1]<-bc[i]*pc12^(n-i)*pc13^i; mc[n+1]<-pc13^n;}
   else {if(n==1) mc[2]<-pc13;}
       mc }
    
silicium<-function(mc){#correction of natural mass distribution accounting for Si
  pSi0<-0.9223; pSi1<-0.0467; pSi2<-0.031;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pSi0;
     m[2]<-mc[2]*pSi0+mc[1]*pSi1;
  for(i in 3:n) m[i]<-mc[i]*pSi0+mc[i-1]*pSi1+mc[i-2]*pSi2;
     m[n+1]<-mc[n]*pSi1+mc[n-1]*pSi2;
       m[n+2]<-mc[n]*pSi2;
         m }
      
sulfur<-function(mc){#correction of natural mass distribution accounting for S
  pS0<-0.9504; pS1<-0.0075; pS2<-0.0421;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pS0;
     m[2]<-mc[2]*pS0+mc[1]*pS1;
  for(i in 3:n) m[i]<-mc[i]*pS0+mc[i-1]*pS1+mc[i-2]*pS2;
     m[n+1]<-mc[n]*pS1+mc[n-1]*pS2;
       m[n+2]<-mc[n]*pS2;
         m }
      
tdist<-function(nc,nsi,ns){ #final natural mass distribution
   bcc<-binom(nc);
    m<-carbons(nc, bcc);
 if(nsi>0) for(i in 1:nsi)  m<-silicium(m);
 if(ns>0) for(i in 1:ns)  m<-sulfur(m);
      m }
 
norm<-function(msd,f1){#normalization of GCMS data accounting for the loss of protons
        colon<-length(msd);# f1<-f1[1];
    for (i in 3:(colon-1)) msd[i]<-(msd[i]+msd[i]*f1-msd[i+1]*f1);
   scol<- msd[3]; for(i in 4:(colon-1)) scol<-scol+msd[i];
  msd[2]<-msd[2]*0; 
 for (i in 3:(colon-1)) msd[i]<-(msd[i]/scol); #normalization
  scol<- msd[3]; for(i in 4:(colon-1)) scol<-scol+msd[i];#=1
    msd[colon]<-scol; 
      msd }
   
mtr<-function(nfrg,nmass,nC,nSi,nS){ #various numbers of labeled 13C with tails of natural distributions
    mt<-tdist(nC,nSi,nS); 
     lem<-length(mt); if(lem<nmass) for(i in (lem+1):nmass) mt[i]<-0.0;
      mt<-mt[1:nmass];
       mt<-mt/sum(mt);
        mm<-numeric(nmass*(nfrg+1));          
         attr(mm,"dim")<-c(nfrg+1,nmass);
          mm[1,]<-mt;
       for(i in 1:nfrg) { mt<-tdist(nC-i,nSi,nS); lem<-length(mt);
         if(lem<(nmass-i)) {for(k in (lem+1):(nmass-i)) mt[k]<-0.0;}
             mt<-mt[1:(nmass-i)];  mt<-mt/sum(mt);
               for(j in 1:length(mt)) mm[i+1,j+i]<-mt[j];}
              mm}
              
mdistr<-function(nreal,msd,mm,nln){ #label incorporation             
          fr<-msd; colon<-length(fr); for(i in 2:colon) fr[i]<-0;
 for(j in 1:(nreal+1)) {fr[j+1]<-msd[j+2]/mm[j,j];
   for(i in j:(colon-2)) msd[i+2]<-msd[i+2]-fr[j+1]*mm[j,i];}
# normalization:
  fr[colon]=fr[2];
        for(i in 3:(colon-2)) {fr[colon]<- fr[colon]+fr[i];}
#     for(i in 2:(colon-2)) fr[i]<- fr[i]/fr[colon];
      fr }

stat<-function(dfr,nln,nfrg,first,last){            i<-1; lfr<-length(dfr);
  dfr1=data.frame();
  while(i<nln){ cnlab<-substr(as.character(dfr[i,1]),first,last); k<-i;
    while(substr(as.character(dfr[k+1,1]),first,last)==cnlab) k<-k+1; 
      if(k>i) {sredn<-dfr[i,];
                 for(j in 2:lfr) sredn[j]<-mean(dfr[i:k,j]);
                 for(j in 2:lfr) if(sredn[j]<0.) sredn[j]<-0.;
#                 srsum<-sum(sredn[2:(2+nfrg)]);
#               sredn[2:(2+nfrg)]<-sredn[2:(2+nfrg)]/srsum;
               otkl<-dfr[1,]; otkl[1]<-as.factor("**sd**");
                 for(j in 2:length(otkl)) otkl[j]<-sd(dfr[i:k,j]);
        dfr1<-rbind(dfr1,sredn,otkl);
                  i<-k;}
       else{dfr1=rbind(dfr1,dfr[i,]);}
          i<-i+1;}
          return(list(dfr1,nln))
} 

suminj<-function(dfr,nln,colon,first,last){         i<-1;
  while(i<nln-1){ cnlab<-substr(as.character(dfr[i,1]),first,last); 
    while(substr(as.character(dfr[i+1,1]),first,last)==cnlab) {
      dfr[i,2:colon]<-dfr[i,2:colon]+dfr[i+1,2:colon];
        dfr<-dfr[-(i+1),]; nln<-nln-1;}
             i<-i+1;}
             return(list(dfr,nln))
}
remsd<-function(mid,lbl){ numl=nrow(mid);
    k=1
        for(i in 1:(numl)){
         if(as.character(mid[k,1])==lbl) {mid=mid[-k,]}
          else {k=k+1}
        }
    return(mid)
}

chirow<-function(dfr,draw,nmi,dlfrg,dteor){     len=length(draw);
        mid=numeric(nmi); tmp<-draw;
         for(j in 1:(dlfrg+1)) mid=mid+dfr[1,j+1]*dteor[j,];
          tmp[2:(2+dlfrg)]= ((draw[2:(2+dlfrg)]-mid[1:(1+dlfrg)])/0.005)^2;
          tmp[len]=sqrt(sum(tmp[2:(2+dlfrg)]));
    return(tmp)
}

alchi<-function(dfr,draw,nmi,dlfrg,dteor){     tmp<-draw; numl=nrow(draw);
            for(i in 1:(numl)){ tmp[i,]=chirow(dfr[i,],draw[i,],nmi,dlfrg,dteor)};
    return(tmp)
}

chandis<-function(mdis,pvar,pfix,fac){
   len=length(mdis); suma=0;
    for(i in 2:len) if(i!=pfix) suma=suma+mdis[i];
     pnew=fac*mdis[pvar]; delta=mdis[pvar]-pnew
  mdis[pvar]=pnew;   sumb=suma-delta
    for(i in 2:len) if(i!=pfix) mdis[i]=mdis[i]*suma/sumb;
   return(mdis)
}

ci99<-function(mid,rada,nadi,pvar,fac,nmass,nfrg,chi){
      len=length(rada)
     while(chi<6.63){ if(mid[pvar]>0.99){break}; if(mid[pvar]<0.001){break}
#     for(i in 1:5){ if(mid[pvar]>0.99){break}; if(mid[pvar]<0.001){break}
     mid=chandis(mid,pvar,99,fac);
     tmp=chirow(mid,rada,nmass,nfrg,nadi); chi=tmp[len]
    }
       return(list(mid,chi))
}

dechi<-function(mid0,pvar,pfix,fac,rada,nadi,nmass,nfrg,chi0){
     chi=chi0; len=length(rada); mid=mid0
     for(i in 1:2){      if((mid[pvar]>0.99)|(mid[pvar]<0.005)) {break}
      mid=chandis(mid,pvar,pfix,fac)
      tmp=chirow(mid,rada,nmass,nfrg,nadi); chi=tmp[len];
        if(chi>chi0) {fac=1./fac}
         else { mid0=mid; chi0=chi;}
     }
       return(list(mid0,chi0))
}

halfci<-function(mid,rada,fac,fac1,mmlab,nmass,nfrg, ncyc){
     chi=0; len=length(mid); nma=which(mid==max(mid[2:len]));
     for(i in 1:ncyc){
   lst=ci99(mid,rada,mmlab,nma,fac,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]
   for(j in 2:(nfrg+2)) if(j!=nma){lst=dechi(mid,j,nma,fac1,rada,mmlab,nmass,nfrg,chi); mid=lst[[1]]; chi=lst[[2]]}
     }
     return(lst)
}
confin<-function(mid0,rada,fac,fac1,mmlab,nmass,nfrg, ncyc,fn1){
  lst=halfci(mid0,rada,fac,fac1,mmlab,nmass,nfrg,11)
     cat(" mid: ",as.character(format(lst[[1]],digits=4))," chi=",as.character(lst[[2]]),"\n",file=fn1,append=TRUE)
    fac=1/fac;
  lst=halfci(mid0,rada,fac,fac1,mmlab,nmass,nfrg,11)
     cat(" mid: ",as.character(format(lst[[1]],digits=4))," chi=",as.character(lst[[2]]),"\n",file=fn1,append=TRUE)
}

