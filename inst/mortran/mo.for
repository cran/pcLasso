c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine pclasso(no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam,ulam,  th    177 
     *r,maxit,verbose,ao,ia,kin,nlp,jerr)
      implicit double precision(a-h,o-z)                                    178
      double precision y(no),x(no,ni),w(no),aa(ni,*),ulam(nlam),ao(nx,nl    179 
     *am)
      integer ia(nx),kin(nlam),mg(ng+1),verbose                             180
      double precision, dimension (:), allocatable :: a,r,sx,s,uu,wt            
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           184
      if(jerr.ne.0) return                                                  185
      allocate(r(1:no),stat=jerr)                                           185
      if(jerr.ne.0) return                                                  186
      allocate(mm(1:ni),stat=jerr)                                          186
      if(jerr.ne.0) return                                                  187
      allocate(sx(1:ni),stat=jerr)                                          187
      if(jerr.ne.0) return                                                  188
      allocate(s(1:ni),stat=jerr)                                           188
      if(jerr.ne.0) return                                                  189
      allocate(uu(1:ni),stat=jerr)                                          189
      if(jerr.ne.0) return                                                  190
      allocate(wt(1:no),stat=jerr)                                          190
      if(jerr.ne.0) return                                                  192
      wt=no*w/sum(w)                                                        195
      k=1                                                                   195
      a=0.0                                                                 195
      uu=0                                                                  195
      mm=0                                                                  195
      s=0.0                                                                 195
      nlp=0                                                                 195
      nin=nlp                                                               195
      r=y                                                                   195
      sw=sum(w)                                                             196
10010 do 10011 j=1,ni                                                       196
      sx(j)=sum(wt*x(:,j)**2)                                               196
10011 continue                                                              196
10012 continue                                                              196
      ysq=sum(w*y**2)/sw                                                    198
      rsq0=1.0                                                              199
10020 do 10021 m=1,nlam                                                     199
      alm=ulam(m)                                                           200
      if(verbose.eq.1) call dblepr("lambda=",-1,alm,1);                         
10030 continue                                                              203
10031 continue                                                              203
      nlp=nlp+1                                                             204
10040 do 10041 kg=1,ng                                                      204
      lg=mg(kg+1)-mg(kg)                                                    205
10050 do 10051 k=1,lg                                                       205
      j=k+mg(kg)-1                                                          205
      gj=dot_product(wt*r,x(:,j))                                           206
      aj=a(j)                                                               206
      u=gj+sx(j)*aj                                                         206
      aajj=aa(j,j-mg(kg)+1)                                                 208
      s(j)=uu(j)-aajj*aj                                                    209
      u=u-theta*s(j)                                                        209
      v=abs(u)-alm                                                          209
      a(j)=0.0                                                              210
      if(v.gt.0.0) a(j)=sign(v,u)/(sx(j)+theta*aajj)                        212
      if(a(j).eq.aj)goto 10051                                              214
      if(mm(j) .ne. 0)goto 10071                                            214
      nin=nin+1                                                             214
      if(nin.gt.nx)goto 10052                                               215
      mm(j)=nin                                                             215
      ia(nin)=j                                                             216
10071 continue                                                              217
      del=a(j)-aj                                                           217
      r=r-del*x(:,j)                                                        218
10080 do 10081 l=1,lg                                                       218
      jl=l+mg(kg)-1                                                         218
      uu(jl)=uu(jl)+del*aa(jl,k)                                            218
10081 continue                                                              220
10082 continue                                                              220
10051 continue                                                              221
10052 continue                                                              221
10041 continue                                                              222
10042 continue                                                              222
      if(nin.gt.nx)goto 10032                                               223
      if(nlp .le. maxit)goto 10101                                          223
      jerr=-m                                                               223
      return                                                                223
10101 continue                                                              224
      rsq=sum(w*r**2)/(ysq*sw)                                              225
      if(abs(rsq0-rsq).lt.thr)goto 10032                                    226
      rsq0=rsq                                                              228
      goto 10031                                                            229
10032 continue                                                              229
      if(nin .le. nx)goto 10121                                             229
      jerr=-10000-m                                                         229
      goto 10022                                                            229
10121 continue                                                              230
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 230
      kin(m)=nin                                                            231
      me=0                                                                  231
10130 do 10131 j=1,nin                                                      231
      if(ao(j,m).ne.0.0) me=me+1                                            231
10131 continue                                                              231
10132 continue                                                              231
      if(me.gt.ne)goto 10022                                                233
10021 continue                                                              234
10022 continue                                                              234
      deallocate(a,r,mm,sx,s,uu)                                            235
      return                                                                236
      end                                                                   237
      subroutine uncomp(ni,ao,ia,kin,a)                                     238
      implicit double precision(a-h,o-z)                                    239
      double precision ao(*),a(ni)                                          239
      integer ia(*)                                                         240
      a=0.0                                                                 240
      if(kin.gt.0) a(ia(1:kin))=ao(1:kin)                                   241
      return                                                                242
      end                                                                   243
      subroutine modval(ao,ia,kin,n,x,f)                                    244
      implicit double precision(a-h,o-z)                                    245
      double precision ao(kin),x(n,*),f(n)                                  245
      integer ia(kin)                                                       246
      f=a0                                                                  246
      if(kin.le.0) return                                                   247
10140 do 10141 i=1,n                                                        247
      f(i)=f(i)+dot_product(ao(1:kin),x(i,ia(1:kin)))                       247
10141 continue                                                              248
10142 continue                                                              248
      return                                                                249
      end                                                                   341
      subroutine logpclasso(no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam,ulam,     343 
     * thr,maxit,verbose,a0,ao,ia,kin,nlp,jerr)
      implicit double precision(a-h,o-z)                                    344
      double precision y(no),x(no,ni),w(no),a0(nlam),aa(ni,*),ulam(nlam)    345 
     *,ao(nx,nlam)
      integer ia(nx),kin(nlam),mg(ng+1),verbose                             346
      double precision, dimension (:), allocatable :: a,r,sx,s,uu,wt,           
     *zz,pr                                                                     
      double precision, dimension (:), allocatable :: eta,eta0,warg             
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           352
      if(jerr.ne.0) return                                                  353
      allocate(r(1:no),stat=jerr)                                           353
      if(jerr.ne.0) return                                                  354
      allocate(mm(1:ni),stat=jerr)                                          354
      if(jerr.ne.0) return                                                  355
      allocate(sx(1:ni),stat=jerr)                                          355
      if(jerr.ne.0) return                                                  356
      allocate(s(1:ni),stat=jerr)                                           356
      if(jerr.ne.0) return                                                  357
      allocate(uu(1:ni),stat=jerr)                                          357
      if(jerr.ne.0) return                                                  358
      allocate(wt(1:no),stat=jerr)                                          358
      if(jerr.ne.0) return                                                  359
      allocate(zz(1:no),stat=jerr)                                          359
      if(jerr.ne.0) return                                                  360
      allocate(pr(1:no),stat=jerr)                                          360
      if(jerr.ne.0) return                                                  361
      allocate(eta0(1:no),stat=jerr)                                        361
      if(jerr.ne.0) return                                                  362
      allocate(eta(1:no),stat=jerr)                                         362
      if(jerr.ne.0) return                                                  363
      allocate(warg(1:no),stat=jerr)                                        363
      if(jerr.ne.0) return                                                  365
      xminw=.0001                                                           366
      del2thr=.01                                                           369
      wt=no*w/sum(w)                                                        371
      uu=0                                                                  371
      mm=0                                                                  371
      s=0.0                                                                 371
      nlp=0                                                                 371
      nin=nlp                                                               371
      r=y                                                                   371
      sw=sum(w)                                                             372
10150 do 10151 j=1,ni                                                       372
      sx(j)=sum(wt*x(:,j)**2)                                               372
10151 continue                                                              372
10152 continue                                                              372
      ysq=sum(w*y**2)/sw                                                    374
      az=0.0                                                                374
      a=0.0                                                                 374
      it=0                                                                  374
      nlp=0                                                                 374
      nth=nlp                                                               374
      kp=nth                                                                374
      kth=kp                                                                374
      r=y                                                                   375
      warg=w                                                                375
      zz=y                                                                  375
      eta=0                                                                 375
      pr=1/(1+exp(-eta))                                                    377
      w=0.25                                                                379
      w=warg*w                                                              380
      w=no*w/sum(w)                                                         380
      sw=sum(w)                                                             380
      a0=0                                                                  380
      mlam=0                                                                382
      rsq0=1.0                                                              383
10160 do 10161 m=1,nlam                                                     383
      alm=ulam(m)                                                           385
      if(verbose.eq.1) call dblepr("lambda=",-1,alm,1);                         
10170 continue                                                              389
10171 continue                                                              389
10180 continue                                                              390
10181 continue                                                              390
      nlp=nlp+1                                                             390
      dlx=0.0                                                               391
      aj=az                                                                 391
      del=dot_product(w,r)/sw                                               392
      az=az+del                                                             392
      r=r-del                                                               393
10190 do 10191 kg=1,ng                                                      393
      lg=mg(kg+1)-mg(kg)                                                    394
10200 do 10201 k=1,lg                                                       394
      j=k+mg(kg)-1                                                          394
      gj=dot_product(wt*r,x(:,j))                                           395
      aj=a(j)                                                               395
      u=gj+sx(j)*aj                                                         395
      aajj=aa(j,j-mg(kg)+1)                                                 397
      s(j)=uu(j)-aajj*aj                                                    398
      u=u-theta*s(j)                                                        398
      v=abs(u)-alm                                                          398
      a(j)=0.0                                                              399
      if(v.gt.0.0) a(j)=sign(v,u)/(sx(j)+theta*aajj)                        401
      if(a(j).eq.aj)goto 10201                                              403
      if(mm(j) .ne. 0)goto 10221                                            403
      nin=nin+1                                                             403
      if(nin.gt.nx)goto 10202                                               404
      mm(j)=nin                                                             404
      ia(nin)=j                                                             405
10221 continue                                                              406
      del=a(j)-aj                                                           406
      r=r-del*x(:,j)                                                        407
10230 do 10231 l=1,lg                                                       407
      jl=l+mg(kg)-1                                                         407
      uu(jl)=uu(jl)+del*aa(jl,k)                                            407
10231 continue                                                              408
10232 continue                                                              408
10201 continue                                                              409
10202 continue                                                              409
10191 continue                                                              412
10192 continue                                                              412
      if(nin.gt.nx)goto 10182                                               413
      if(nlp .le. maxit)goto 10251                                          413
      jerr=-m                                                               413
      return                                                                413
10251 continue                                                              414
      rsq=sum(w*r**2)/(ysq*sw)                                              416
      if(abs(rsq0-rsq).lt.thr)goto 10182                                    417
      rsq0=rsq                                                              418
      goto 10181                                                            422
10182 continue                                                              422
      eta0=eta                                                              422
      eta=zz-r                                                              423
      del2=sum(abs(eta-eta0))/no                                            424
      pr=1/(1+exp(-eta))                                                    427
      w=0.25                                                                430
      zz=eta+(y-pr)/w                                                       433
      w=warg*w                                                              433
      w=no*w/sum(w)                                                         433
      sw=sum(w)                                                             434
      a0(m)=az                                                              434
      r=zz-eta                                                              437
      if(del2.lt.del2thr)goto 10172                                         437
      goto 10171                                                            442
10172 continue                                                              442
      if(nin .le. nx)goto 10271                                             442
      jerr=-10000-m                                                         442
      goto 10162                                                            442
10271 continue                                                              443
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 443
      kin(m)=nin                                                            444
      me=0                                                                  444
10280 do 10281 j=1,nin                                                      444
      if(ao(j,m).ne.0.0) me=me+1                                            444
10281 continue                                                              444
10282 continue                                                              444
      if(me.gt.ne)goto 10162                                                446
10161 continue                                                              448
10162 continue                                                              448
      deallocate(a,r,mm,sx,s,uu)                                            449
      return                                                                450
      end                                                                   455
