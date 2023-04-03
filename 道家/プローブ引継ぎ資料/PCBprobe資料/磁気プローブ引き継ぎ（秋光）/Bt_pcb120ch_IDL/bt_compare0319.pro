;+
; :Author: Moe Akimitsu
;-1次元解析
Pro Bt_compare0319,shot,pf_ref=pf_ref,tf_ref=tf_ref
  ;  shot=588
  ;

  ;stop

  device,decomposed=0
  !P.background=16777215L
  !P.color=0
    loadct,39

  pathname="o:\work\KENKYU\magnetic_probe\PCB120ch"
  cd, pathname
  
  ;cd, "o:\work\KENKYU\magnetic_probe\magnetic"
  ;cd, "C:\Users\moeru\OneDrive - The University of Tokyo\work\KENKYU\magnetic_probe\magnetic"
  ;filename="181115.csv"
;  log=read_csv(filename)
;
;  Dac=log.field2
;  PFset=log.field4
;  TFset=log.field5
;  Gasset=fix(log.field6)
;  rgwshot=log.field1
;
;  ts=[3.4,2.5,3.0,2.2]
;  ps=[39,34]
;  tref=[513,500,507,497]
;  ;tref=[513,500,507]
;  pref=[559,501]
  n_frame=50.
  t0=460;800.;750;800.;840.;800.;
  t_frame=t0+findgen(n_frame)*1.0

  ;Vac_ref=371
  Z_resol=5.
  R_resol=24.
  n_grid=60.
  myu0 = 4.*!pi*1.e-7
  coeff=read_ascii("RC_ANS_120.txt")
  coeff=coeff.field1
  direction=read_ascii("direction0120.txt")
  direction=direction.field1
  rpos=read_ascii("position_r.txt")
  rpos=rpos.field1
  zpos=fltarr(120)
  zpos5=[-0.042,-0.021,0.,0.021,0.042]
  for i=0,4 do zpos[24*i:24*i+23]=zpos5[i]
  
   zout=findgen(n_grid+1)/n_grid*(max(zpos)-min(zpos))+min(zpos);[-0.205,-0.155,-0.105,-0.055,-0.025,-0.005,0.005,0.025,0.055,0.105,0.155,0.205]
  rout=findgen(n_grid+1)/n_grid*(max(rpos)-min(rpos))+min(rpos);[0.09,0.11,0.15,0.19,0.23]
  Bt_2D=fltarr([n_elements(zout),n_elements(rout),n_elements(t_frame)])

;  for tscan=0,n_elements(tref)-1 do begin
;
;    for pscan=0,n_elements(pref)-1 do begin
;      ; ;TF_ref=568
      ; TF_ref=592
      ; PF_ref=559
      ;PF_ref=579


;      PF_ref=pref[pscan]
;      TF_ref=tref[tscan]
;      pla_No=where((dac*pfset*tfset*gasset*log.field3 ne 0.)*(pfset eq ps[pscan])*(tfset*10 eq ts[tscan]*10))
      ;stop
   ;   if total(pla_No) ne -1 then begin
      ;  for nn=0,n_elements(pla_no)-1 do begin
          ;shot=dac[pla_no[nn]]
          print,shot;,pfset,tfset
          ;stop
          ;TF_ref=568
          ;rpos[6]=rpos[6]-0.005

          ;for
          ;  ;++++++++++++++++++++++++osiro+++++++++++++++++++++++++++++++++++++
          ;   filename=dialog_pickfile(path="X:\results\ts-3u\",filter="*.rgw")
          ;  filename="X:\results\ts-3u\181117\181117"+strcompress(rgwshot/100,/remove_all)+strcompress(rgwshot/10,/remove_all)+strcompress(rgwshot mod 10,/remove_all)+".rgw"
          ;  data=read_ascii(filename,data_start=1) & data=data.field01(1:*,0:*)
          ;  otime=data[0,*]
          ;  Itfc= -129192*12 * reform(data[3,*])
          ;
          ;   bt_cal=myu0*Itfc[i]*1.e3/(2.*!pi*reform(position[1,*]))
          ;  ;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          cd,"C:\Program Files\MDSplus\idl"
          mdsconnect, '192.168.1.140' ;mdsconnect, 'IP address'
          CH_i=[1:128]; integrator CH
          chord1=[1:31:2]
          chord2=[2:32:2]
          chord3=[33:61:2]
          chord4=[34:62:2]
          chord5=[65:128]
          CH_d=[chord1,chord2,chord3,63,chord4,64,chord5] ;these are to make CH orders of digitizer
          mdsopen,'a038',PF_ref ; mdsopen,'treename',shot

          n= mdsvalue("_n=.AI:CH001")
          B_pf=fltarr([n_elements(n),128]);empty array
          ; coeff[30:44]=coeff[30:44]*1.8
          ;coeff[45:59]=coeff[30:44]*2

          for i=1,128 do begin ;shot
            fname='.AI:CH'+strcompress(ch_d[i-1]/100,/remove_all)+strcompress((ch_d[i-1] mod 100)/10,/remove_all)+strcompress((ch_d[i-1] mod 100) mod 10,/remove_all); you have to change the CH00*
            x= mdsvalue("_x="+fname)
            x=x-mean(x(0:40))
            B_pf[*,i-1]=x*coeff[i-1]
          end
          mdsclose,'a038',PF_ref
          mdsopen,'a038',TF_ref ; mdsopen,'treename',shot

          n= mdsvalue("_n=.AI:CH001")
          B_tf=fltarr([n_elements(n),128]);empty array
          ; coeff[30:44]=coeff[30:44]*1.8
          ;coeff[45:59]=coeff[30:44]*2

          for i=1,128 do begin ;shot
             fname='.AI:CH'+strcompress(ch_d[i-1]/100,/remove_all)+strcompress((ch_d[i-1] mod 100)/10,/remove_all)+strcompress((ch_d[i-1] mod 100) mod 10,/remove_all); you have to change the CH00*
     x= mdsvalue("_x="+fname)
            x=x-mean(x(0:40))
            B_tf[*,i-1]=x*coeff[i-1]
          end
          mdsclose,'a038',TF_ref
          mdsopen,'a038',shot ; mdsopen,'treename',shot

          n= mdsvalue("_n=.AI:CH001")
          Bt=fltarr([n_elements(n),128]);empty array
          ; coeff[30:44]=coeff[30:44]*1.8
          ;coeff[45:59]=coeff[30:44]*2

          for i=1,128 do begin ;shot
           fname='.AI:CH'+strcompress(ch_d[i-1]/100,/remove_all)+strcompress((ch_d[i-1] mod 100)/10,/remove_all)+strcompress((ch_d[i-1] mod 100) mod 10,/remove_all); you have to change the CH00*
       x= mdsvalue("_x="+fname)
            x=x-mean(x(0:40))
            Bt[*,i-1]=x*coeff[i-1]
          endfor
          mdsclose,'a038',shot
          mdsdisconnect
          ; Bt[*,59]=-bt[*,59]

          ;

          ;

          ;  stop
          ;  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ;B_pf=bt
          ;save,filename='pfref0816.sav',b_pf
          ;  restore,'C:\Users\Moe Akimitsu\IDLWorkspace\Default\pfref0816.sav';222?
          ;  restore,'C:\Users\Moe Akimitsu\IDLWorkspace\Default\tfref0816.sav';223
          bt=bt-b_pf
          ;+++++++++++++++++++++++++++++++++++++++++
          ;plot,z,r,psym=1 ,/isotropic
          time=findgen(n_elements(x))-500
         ; minus_ind=[187,92,108,82,87,52,[25:48],[97:113],[115:117],119,120];,[20:24],[44:48],[68:72],[92:96],[116:119]]
          neg1=[20,36,111,113,104]
          ;neg2=[18,23,32,84,87,99,109,78,48,52,67];[18,20,23,32,36,87,89,99,103,104,109,111,113];,117,77,78,52]
          neg2=[18,23,99,109,103,67,97,100,114,117,118,81,78];[18,20,23,32,36,87,89,99,103,104,109,111,113];,117,77,78,52],81,87,89
          change1=[1,2,4,14,30,41,46,52,81,82,97,105,109,112,114,118,78,87]
          change2=[7,33,72,15,76,96,127,128,121,122,123,124,125,126,108,119,71,91]
          CH_i=[1:128]; integrator CH

          index=[0:120]
          index[neg1-1]=-1
          index[neg2-1]=-1
          index[change2-1]=-1
          index=where(index ne -1)
         
         ; Bt[*,minus_ind]=-Bt[*,minus_ind]
         for i=0,n_elements(change1)-1 do bt[*,change1[i]-1]=bt[*,change2[i]-1]
           for i=0,n_elements(change1)-1 do b_tf[*,change1[i]-1]=b_tf[*,change2[i]-1]
            minus_ind=where(total(b_tf[60:400,*],1) le 0)
          b_tf[*,minus_ind]=-b_tf[*,minus_ind]
          bt[*,minus_ind]=-bt[*,minus_ind]; 極性反転
          bt[*,change2-1]=bt[*,change2-1]*0+1000.;交換したやつを表示させないため
          bt[*,neg1-1]=bt[*,neg1-1]*0+1000.;もともとないチャンネルも表示させない
          ;index=where(abs( total(bt[500:1000,*],1)-mean(total(bt[500:1000,*],1))) le  0.8*mean(total(bt[500:1500,*],1)))
          ;index=[1,3,[5:13],[15:17],19,21,22,[24:29],[31:34],37,38,40,[42:45],[47:51],53,54,56]-1
          ;
          ;bt=-bt
          ;b_tf=-b_tf

;          index=fltarr(120)
;          index0=-1+[3,[5:13],[15:17],19,21,22,[24:29],[32:34],37,38,40,[42:45],[47:49],51,53,54,56,57,[68:70],[72:120]]
;          index=-1+[3,[5:13],[15:17],19,21,22,[24:29],[32:34],37,38,40,[42:45],[47:49],51,53,54,56,57];,[61:66],[68:75]]
;          err=total(bt[1500:2000,*],1)/total(bt[1500:2000,*],1)
;          err2=err(index)
;          res=createboxplotdata(err2,mean_values=means,outlier_values=outliers)
;          if keyword_set(outliers) eq 1 then begin
;            outliers=reform(outliers[1,*])
;            for out=0,n_elements(outliers)-1 do index[where(index eq outliers[out])]=-1
;            index=index[where(index ne -1)]
;          endif
;          ; stop
;          ;index=[0:74]
;          index2=[0:119]
         ; index[index]=index2[index]
         
          ;   window,1 & !p.multi=0
          ;   plot,total(bt[500:1500,*],1)
          ;   oplot,total(bt[500:1500,index],1)
          ;  stop
          ;  window,10 ,xsize=1000,ysize=600
          ;  !p.multi=[0,9,6]
          ;rbtf=b_tf*0
          ; for i=0,74 do rbtf[*,i]=rpos[i]*bt[*,i]
          ;
          bk=[100:900]
          ; a=[0.,1.]*(2.*!pi)/(myu0*Itfc[bk]*1.e3)
          a=[0,0.04]
          ;  rb_tf_med=rbtf[*,index0]
          ;  rb_tf_med=median(total(rb_tf_med[bk,*],1))
          ;  rb_tf_all=rbtf[*,index]
          ;  rb_tf_all=total(rb_tf_all[bk,*],1)
;          coeffcal=fltarr([n_elements(index)])
;          ;  coeff= rb_tf_med/rb_tf_all
;          rb_tf=fltarr([n_elements(b_tf[*,0])])
;          ;  for i=0,n_elements(b_tf[*,0])-1 Do rb_tf[i]=reform(median(rpos[index0]*b_tf[i,index0],/even))
          ;  for i=0,n_elements(index)-1 do begin
          ;    coeff2=rb_tf[*]/b_tf[*,index[i]]/rpos[index[i]]
          ;
          ;    ; plot,time,coeff2,title=strcompress(index[i]),xr=[0,2000],yr=[0.5,1.5]
          ;    ; oplot,[-1000,3000],[1,1]*median(coeff2[400:1400]),color=200
          ;    bt[*,index[i]]=bt[*,index[i]]*median(coeff2[bk],/even)
          ;    b_tf[*,index[i]]=b_tf[*,index[i]]*median(coeff2[bk],/even)
          ;    coeffcal[i]=median(coeff2[bk],/even)
          ;  endfor

          ;bt=bt/b_tf*
          ;
          ; stop
          ;goto, kousei
          ; stop
          ;for i=0,n_elements(coeff)-1 do begin
          ;  bt[*,index[i]]=bt[*,index[i]]*coeff[i]
          ;  b_tf[*,index[i]]=b_tf[*,index[i]]*coeff[i]
          ;endfor

          ;  Bt_calc=matrix_multiply(myu0*reform(Itfc)*1.e3,(1/(2.*!pi*rpos)))

          ;insitu=reform(mean(Bt_calc[100:300,*],dimension=1)/mean(bt[600:800,*],dimension=1))
          ; for i=0,74 do Bt[*,i]=bt[*,i]*insitu[i]
          ; rpos[45:59]=rpos[15:29]+0.03
          ;rpos=rpos[index]

          ;  itfc=smooth(itfc,3)
          ;
;          bt=smooth(bt,[3,0])
;          b_tf=smooth(b_tf,[3,0])

          ;  stop

          ;t_frame=500+[450,465,466,467,468,469,470,471,473,475]
          ;  window,31,xsize=1200,ysize=700
          ;  !P.multi=[0,n_frame,5]


          ; fit校正
          ;  for k=0,4 do begin
          ;    array=[0:14]+15*k
          ;    array_Ch=where(index[array] gt 0.0)+15*k
          ;;    stop
          ;;   result=fltarr([2,40])
          ;;   for i=0,n_elements(result[0,*])-1 do  result[*,i] = reform(COMFIT(rpos[array_CH],b_tf[i+Bk+500,array_CH], A,/HYPERBOLIC))
          ;;;  stop
          ;;  result=median(result,dimension=2)
          ; ; result = reform(COMFIT(rpos[array_CH],b_tf[Bk+500,array_CH], A,/HYPERBOLIC))
          ; ; print, result
          ;; ::::
          ;;   rpos[array_ch]=rpos[array_ch]+result[0]/result[1]
          ;;   coeff[array_ch]= 1/(result[1]*rpos[array_ch])/b_tf[350+500,array_ch]
          ;;      for i=0,n_elements(bt[*,0])-1 do begin
          ;;        bt[i,array_ch]=bt[i,array_ch]*coeff[array_ch]
          ;;        b_tf[i,array_ch]=b_tf[i,array_ch]*coeff[array_ch]
          ;;      endfor
          ;;;;
          ;;stop
          ;  ; rout=findgen(n_grid+1)/n_grid*(max(rpos)-min(rpos))+min(rpos);[0.09,0.11,0.15,0.19,0.23]
          ; endfor
       ;   logtitle='_'+strcompress(ps[pscan],/remove_all)+'_'+strcompress(fix(ts[tscan]*1000),/remove_all)+'_'+strcompress(gasset[pla_no[nn]],/remove_all)
          ;cd, "o:\work\KENKYU\magnetic_probe\magnetic\"
          ;cd, "\\obi-wan\htanabe_backup\ts3u\IDL\makimitsu\tara-san"
          logtitle='_'+strcompress(shot,/remove_all)
          cd, "C:\Users\moeru\OneDrive - The University of Tokyo\work\KENKYU\magnetic_probe\magnetic\img\"
          FILE_MKDIR,'shot'+strcompress(shot,/remove_all)+logtitle
          cd, 'shot'+strcompress(shot,/remove_all)+logtitle


          FILE_MKDIR,'Bt'
          ;  FILE_MKDIR,'pala'
          ;   FILE_MKDIR,'pla_b'

      xsize=1200
ysize=700
for w=0,4 do begin
  window,w+11, xsize=xsize,ysize=ysize,title='probe'+strcompress(w+1,/remove_all)
  !p.multi=[0,8,3]
  for i=24*w+0,24*w+23 do begin
    title='ch'+strcompress(i+1,/remove_all)
    !p.color=0
    if where(index eq i) eq -1 then begin
      title=title+'-neg'
      !p.color=100
    endif

    plot,(bt-b_tf)[350:600,i],title=title,xst=1,charsize=0.7,Xtickname=REPLICATE('',60),Ytickname=REPLICATE('',60),yr=[-0.005,0.04]
    ; oplot,probe_b[*,i]*0,color=100
    ;oplot,bt_calc[300:700,i],color=200
    ;                plot,probe_b[*,i],title=title,yr=[-0.01,0.10],xst=1,charsize=0.7,Xtickname=REPLICATE('',60),Ytickname=REPLICATE('',60)
    ;                ; oplot,probe_b[*,i]*0,color=100
    ;                oplot,probe_b[*,i]*0.,color=200
  endfor
endfor

          ;WRITE_JPEG , "probe_ch.jpg",TVRD(/TRUE),/TRUE
           stop

          page=n_frame/25

          for k=0,4 do begin
            ;array=[0:23]+24*k
            array_Ch=index[where((index ge 0+24*k)*(index le 23+24*k ))]

            window,k ,xsize=1200,ysize=700

            !P.multi=[0,5,5]
            for l=0,page-1 do begin
              for j=0,n_frame/page-1  do begin
                i=25*l+j
                a=[0.,1.];*(2.*!pi)/(myu0*Itfc[t_frame[i]-500]*1.e3)

                ;   result = reform(COMFIT(rpos[array_CH],b_tf[t_frame[i],array_CH], A,/HYPERBOLIC))
                result = reform(COMFIT(rpos[index],b_tf[t_frame[i],index], A,/HYPERBOLIC))
                ;  print,result
                b_fit=1/(result[0]+result[1]*rout)

                plot,rpos[array_CH],b_tf[t_frame[i],array_CH],psym=1,xst=1,yst=1,yr=[0,0.500],title=strcompress(t_frame[i])
                oplot, rpos[array_CH],bt[t_frame[i],array_CH],psym=1,thick=2,color=100
                ;oplot,rpos[array_CH],bt_calc[(t_frame[i]-400)+500,array_CH],psym=5,color=200
                ; oplot,rpos[array_CH],bt_calc[t_frame[i]-500,array_CH],psym=5,color=200
                oplot,rout,b_fit,color=50
                ;stop
                WRITE_JPEG , "1D_bt"+strcompress(k+1,/remove_all)+'_'+strcompress(l+1,/remove_all)+".jpg",TVRD(/TRUE),/TRUE
              endfor
            endfor
            ;stop
            ;;;;
            ;
            ;kはアレイの本数　1本目から五本目、lはウィンドウ数
            ;
            ;
            ;
            ;
            ;
            ;;;;;;
            window,k+5,xsize=1200,ysize=700
            !P.multi=[0,5,5]
            for l=0,page-1 do begin
              for j=0,n_frame/page-1  do begin
                i=25*l+j
                ;      plot,rpos[array_CH],(b_tf[t_frame[i],array_CH])*rpos[array_CH],psym=1,xst=1,yst=1,yr=[-50,50],thick=2,title=strcompress(t_frame[i]-500.0)
                ;      oplot,rpos[array_CH],(bt[t_frame[i],array_CH])*rpos[array_CH],psym=1,color=50
                plot,rpos[array_CH],(bt[t_frame[i],array_CH]-b_tf[t_frame[i],array_CH])*rpos[array_CH],psym=1,xst=1,yst=1,yr=[-0.010,0.010],thick=2,title=strcompress(t_frame[i])
               
                oplot,[0,0.3],[0,0],psym=0,color=50
                ; WRITE_JPEG , "1D_rbt"+strcompress(k+1,/remove_all)+'_'+strcompress(l+1,/remove_all)+".jpg",TVRD(/TRUE),/TRUE
              endfor
            endfor
          endfor
          ;  stop




          For i=0,n_elements(t_frame)-1 do begin
            X=zpos[index]
            Y=rpos[index]
            ; Z=(bt[t_frame[i],index]-b_tf[t_frame[i],index]);*rpos
            Z=(bt-b_tf)[t_frame[i],index]*rpos
            ; Z=bt[t_frame[i],index];/rpos; 裏技
            ;Triangulate,X,Y,Tr,B1

            FOR k=0, N_ELEMENTS(tr)/3-1 DO t_index = [tr[*,k], tr[0,k]]
            ; Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout,Extra=B1)
            ; Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout)
            ; Bt_2D[*,*,i]=GRID_TPS(x, y, z, NGRID=[n_grid+1, n_grid+1], START=[0,0], DELTA=[1,1])
            Bt_2D[*,*,i]=griddata(X,Y,Z,/RADIAL_BASIS_FUNCTION, FUNCTION_TYPE=3,smoothing=0.005,/grid,xout=zout,yout=rout );smoothing=0.02,
          endfor

          ;bt_2D=smooth(Bt_2D,[3,2,1])
          rout2=fltarr([n_grid+1,n_grid+1])+1
          for ii=0,n_grid do rout2[*,ii]=rout2[*,ii]*rout[ii]
          cd, "Bt"
          window,20,xsize=600,ysize=800,title='BT'
          !p.multi=0
          ;   !P.multi=[0,8,3]
          for i=0,n_frame-1 do begin
            contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol])*rout2,zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i])+"s"),charsize=1.;,zr=[-40,40];*max(abs(Bt_2D));,xtitle="Z [m]",ytitle="R [m]",
            ;   oplot,zpos,rpos,psym=1,color=0
            contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol])*rout2,zout,rout,nlevels=128,/overplot,color=0

            ; contour,Bt_2D[*,*,i],zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i]-500)+"s"),charsize=1.5;,yr=[0.17,0.28];,xtitle="Z [m]",ytitle="R [m]";,zr=[-1.,1.]*max(abs(Bt_2D[*,*,i]))
            oplot,zpos[index],rpos[index],psym=1,color=0
            ;contour,Bt_2D[*,*,i],zout,rout,nlevels=120,/overplot,color=0
            WRITE_JPEG , strcompress("time"+string(fix(t_frame[i]))+"s.jpg",/remove_all),TVRD(/TRUE),/TRUE
          endfor
          save, filename='shot'+strcompress(shot,/remove_all)+'_Btdata.sav',Bt,Bt_2D,zpos,rpos,index,zout,rout,t_frame
          ;stop

          ;  For i=0,n_elements(t_frame)-1 do begin
          ;    X=zpos[index]
          ;    Y=rpos[index]
          ;    Z=(bt[t_frame[i],index]-b_tf[t_frame[i],index])*rpos
          ;    ;Z=bt[t_frame[i],index]*rpos
          ;    Triangulate,X,Y,Tr,B1
          ;
          ;    FOR k=0, N_ELEMENTS(tr)/3-1 DO t_index = [tr[*,k], tr[0,k]]
          ;    ;Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout,Extra=B1)
          ;    Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout)
          ;    ; Bt_2D[*,*,i]=GRID_TPS(x, y, z, NGRID=[n_grid+1, n_grid+1], START=[0,0], DELTA=[1,1])
          ;  endfor
          ;
          ;  ;bt_2D=smooth(Bt_2D,[3,2,1])
          ;  cd, "../pala"
          ;  window,20,xsize=600,ysize=800,title='plasma_rBt'
          ;  !p.multi=0
          ;  ;   !P.multi=[0,8,3]
          ;  for i=0,n_frame-1 do begin
          ;    contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol]),zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i]-500)+"s"),charsize=1.,zr=[-5,5];*max(abs(Bt_2D));,xtitle="Z [m]",ytitle="R [m]",
          ;    ;   oplot,zpos,rpos,psym=1,color=0
          ;    contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol]),zout,rout,nlevels=128,/overplot,color=0
          ;
          ;    ; contour,Bt_2D[*,*,i],zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i]-500)+"s"),charsize=1.5;,yr=[0.17,0.28];,xtitle="Z [m]",ytitle="R [m]";,zr=[-1.,1.]*max(abs(Bt_2D[*,*,i]))
          ;    ;oplot,zpos[index],rpos[index],psym=1,color=0
          ;    ;contour,Bt_2D[*,*,i],zout,rout,nlevels=120,/overplot,color=0
          ;    WRITE_JPEG , "time"+strcompress(fix(t_frame[i]-500),/remove_all)+"s.jpg",TVRD(/TRUE),/TRUE
          ;  endfor
          ;
          ;  For i=0,n_elements(t_frame)-1 do begin
          ;    X=zpos[index]
          ;    Y=rpos[index]
          ;    Z=(bt[t_frame[i],index]-b_tf[t_frame[i],index])
          ;    ;Z=bt[t_frame[i],index]*rpos
          ;    Triangulate,X,Y,Tr,B1
          ;
          ;    FOR k=0, N_ELEMENTS(tr)/3-1 DO t_index = [tr[*,k], tr[0,k]]
          ;    ;Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout,Extra=B1)
          ;    Bt_2D[*,*,i]=Trigrid(X,Y,Z,Tr,xout=zout,yout=rout)
          ;    ; Bt_2D[*,*,i]=GRID_TPS(x, y, z, NGRID=[n_grid+1, n_grid+1], START=[0,0], DELTA=[1,1])
          ;  endfor
          ;
          ;
          ;  ;bt_2D=smooth(Bt_2D,[3,2,1])
          ;  cd, "..\pla_b"
          ;  window,20,xsize=600,ysize=800,title='plasmaBt'
          ;  !p.multi=0
          ;  ;   !P.multi=[0,8,3]
          ;  for i=0,n_frame-1 do begin
          ;    contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol]),zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i]-500)+"s"),charsize=1.,zr=[-20,20];*max(abs(Bt_2D));,xtitle="Z [m]",ytitle="R [m]",
          ;    ;   oplot,zpos,rpos,psym=1,color=0
          ;    contour,smooth(Bt_2D[*,*,i],[n_grid/z_resol,n_grid/r_resol]),zout,rout,nlevels=128,/overplot,color=0
          ;
          ;    ; contour,Bt_2D[*,*,i],zout,rout,nlevels=256,/fill,/isotropic,xst=1,yst=1,zst=1,title=strcompress("time ="+string(t_frame[i]-500)+"s"),charsize=1.5;,yr=[0.17,0.28];,xtitle="Z [m]",ytitle="R [m]";,zr=[-1.,1.]*max(abs(Bt_2D[*,*,i]))
          ;    ;oplot,zpos[index],rpos[index],psym=1,color=0
          ;    ;contour,Bt_2D[*,*,i],zout,rout,nlevels=120,/overplot,color=0
          ;    WRITE_JPEG , "time"+strcompress(fix(t_frame[i]-500),/remove_all)+"s.jpg",TVRD(/TRUE),/TRUE
          ;  endfor
          ;stop

;        endfor
;        ;kousei:stop
;      endif
;
;    endfor
;    ;stop
;  endfor
  ;stop
End
