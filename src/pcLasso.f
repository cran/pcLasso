c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine pclasso(no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam,ulam,  th
     *r,maxit,verbose,ao,ia,kin,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),x(no,ni),w(no),aa(ni,*),ulam(nlam),ao(nx,nl
     *am)
      integer ia(nx),kin(nlam),mg(ng+1),verbose                         
      double precision, dimension (:), allocatable :: a,r,sx,s,uu,wt    
      integer, dimension (:), allocatable :: mm                         
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(sx(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(s(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(uu(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(wt(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      wt=no*w/sum(w)                                                    
      k=1                                                               
      a=0.0                                                             
      uu=0                                                              
      mm=0                                                              
      s=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      r=y                                                               
      sw=sum(w)                                                         
      do 10011 j=1,ni                                                   
      sx(j)=sum(wt*x(:,j)**2)                                           
10011 continue                                                          
      continue                                                          
      ysq=sum(w*y**2)/sw                                                
      rsq0=1.0                                                          
      do 10021 m=1,nlam                                                 
      alm=ulam(m)                                                       
      if(verbose.eq.1) call dblepr("lambda=",-1,alm,1);                 
      continue                                                          
10031 continue                                                          
      nlp=nlp+1                                                         
      do 10041 kg=1,ng                                                  
      lg=mg(kg+1)-mg(kg)                                                
      do 10051 k=1,lg                                                   
      j=k+mg(kg)-1                                                      
      gj=dot_product(wt*r,x(:,j))                                       
      aj=a(j)                                                           
      u=gj+sx(j)*aj                                                     
      aajj=aa(j,j-mg(kg)+1)                                             
      s(j)=uu(j)-aajj*aj                                                
      u=u-theta*s(j)                                                    
      v=abs(u)-alm                                                      
      a(j)=0.0                                                          
      if(v.gt.0.0) a(j)=sign(v,u)/(sx(j)+theta*aajj)                    
      if(a(j).eq.aj)goto 10051                                          
      if(mm(j) .ne. 0)goto 10071                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10052                                           
      mm(j)=nin                                                         
      ia(nin)=j                                                         
10071 continue                                                          
      del=a(j)-aj                                                       
      r=r-del*x(:,j)                                                    
      do 10081 l=1,lg                                                   
      jl=l+mg(kg)-1                                                     
      uu(jl)=uu(jl)+del*aa(jl,k)                                        
10081 continue                                                          
      continue                                                          
10051 continue                                                          
10052 continue                                                          
10041 continue                                                          
      continue                                                          
      if(nin.gt.nx)goto 10032                                           
      if(nlp .le. maxit)goto 10101                                      
      jerr=-m                                                           
      return                                                            
10101 continue                                                          
      rsq=sum(w*r**2)/(ysq*sw)                                          
      if(abs(rsq0-rsq).lt.thr)goto 10032                                
      rsq0=rsq                                                          
      goto 10031                                                        
10032 continue                                                          
      if(nin .le. nx)goto 10121                                         
      jerr=-10000-m                                                     
      goto 10022                                                        
10121 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      me=0                                                              
      do 10131 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
10131 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10022                                            
10021 continue                                                          
10022 continue                                                          
      deallocate(a,r,mm,sx,s,uu)                                        
      return                                                            
      end                                                               
      subroutine uncomp(ni,ao,ia,kin,a)                                 
      implicit double precision(a-h,o-z)                                
      double precision ao(*),a(ni)                                      
      integer ia(*)                                                     
      a=0.0                                                             
      if(kin.gt.0) a(ia(1:kin))=ao(1:kin)                               
      return                                                            
      end                                                               
c NOT CURRENTLY USED; f=0 maybe should be f=a0
      subroutine modval(ao,ia,kin,n,x,f)                                
      implicit double precision(a-h,o-z)                                
      double precision ao(kin),x(n,*),f(n)                              
      integer ia(kin)                                                   
      f=0                                                              
      if(kin.le.0) return                                               
      do 10141 i=1,n                                                    
      f(i)=f(i)+dot_product(ao(1:kin),x(i,ia(1:kin)))                   
10141 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine logpclasso(no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam,ulam, 
     * thr,maxit,verbose,a0,ao,ia,kin,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),x(no,ni),w(no),a0(nlam),aa(ni,*),ulam(nlam)
     *,ao(nx,nlam)
      integer ia(nx),kin(nlam),mg(ng+1),verbose                         
      double precision, dimension (:), allocatable :: a,r,sx,s,uu,wt,   
     *zz,pr                                                             
      double precision, dimension (:), allocatable :: eta,eta0,warg     
      integer, dimension (:), allocatable :: mm                         
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(sx(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(s(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(uu(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(wt(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(zz(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(pr(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(eta0(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(eta(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(warg(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      xminw=.0001                                                       
      del2thr=.01                                                       
      wt=no*w/sum(w)                                                    
      uu=0                                                              
      mm=0                                                              
      s=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      r=y                                                               
      sw=sum(w)                                                         
      do 10151 j=1,ni                                                   
      sx(j)=sum(wt*x(:,j)**2)                                           
10151 continue                                                          
      continue                                                          
      ysq=sum(w*y**2)/sw                                                
      az=0.0                                                            
      a=0.0                                                             
      it=0                                                              
      nlp=0                                                             
      nth=nlp                                                           
      kp=nth                                                            
      kth=kp                                                            
      r=y                                                               
      warg=w                                                            
      zz=y                                                              
      eta=0                                                             
      pr=1/(1+exp(-eta))                                                
      w=0.25                                                            
      w=warg*w                                                          
      w=no*w/sum(w)                                                     
      sw=sum(w)                                                         
      a0=0                                                              
      mlam=0                                                            
      rsq0=1.0                                                          
      do 10161 m=1,nlam                                                 
      alm=ulam(m)                                                       
      if(verbose.eq.1) call dblepr("lambda=",-1,alm,1);                 
      continue                                                          
10171 continue                                                          
      continue                                                          
10181 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      aj=az                                                             
      del=dot_product(w,r)/sw                                           
      az=az+del                                                         
      r=r-del                                                           
      do 10191 kg=1,ng                                                  
      lg=mg(kg+1)-mg(kg)                                                
      do 10201 k=1,lg                                                   
      j=k+mg(kg)-1                                                      
      gj=dot_product(wt*r,x(:,j))                                       
      aj=a(j)                                                           
      u=gj+sx(j)*aj                                                     
      aajj=aa(j,j-mg(kg)+1)                                             
      s(j)=uu(j)-aajj*aj                                                
      u=u-theta*s(j)                                                    
      v=abs(u)-alm                                                      
      a(j)=0.0                                                          
      if(v.gt.0.0) a(j)=sign(v,u)/(sx(j)+theta*aajj)                    
      if(a(j).eq.aj)goto 10201                                          
      if(mm(j) .ne. 0)goto 10221                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10202                                           
      mm(j)=nin                                                         
      ia(nin)=j                                                         
10221 continue                                                          
      del=a(j)-aj                                                       
      r=r-del*x(:,j)                                                    
      do 10231 l=1,lg                                                   
      jl=l+mg(kg)-1                                                     
      uu(jl)=uu(jl)+del*aa(jl,k)                                        
10231 continue                                                          
      continue                                                          
10201 continue                                                          
10202 continue                                                          
10191 continue                                                          
      continue                                                          
      if(nin.gt.nx)goto 10182                                           
      if(nlp .le. maxit)goto 10251                                      
      jerr=-m                                                           
      return                                                            
10251 continue                                                          
      rsq=sum(w*r**2)/(ysq*sw)                                          
      if(abs(rsq0-rsq).lt.thr)goto 10182                                
      rsq0=rsq                                                          
      goto 10181                                                        
10182 continue                                                          
      eta0=eta                                                          
      eta=zz-r                                                          
      del2=sum(abs(eta-eta0))/no                                        
      pr=1/(1+exp(-eta))                                                
      w=0.25                                                            
      zz=eta+(y-pr)/w                                                   
      w=warg*w                                                          
      w=no*w/sum(w)                                                     
      sw=sum(w)                                                         
      a0(m)=az                                                          
      r=zz-eta                                                          
      if(del2.lt.del2thr)goto 10172                                     
      goto 10171                                                        
10172 continue                                                          
      if(nin .le. nx)goto 10271                                         
      jerr=-10000-m                                                     
      goto 10162                                                        
10271 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      me=0                                                              
      do 10281 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
10281 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10162                                            
10161 continue                                                          
10162 continue                                                          
      deallocate(a,r,mm,sx,s,uu)                                        
      return                                                            
      end                                                               
