PRO Multi288CH
!P.multi=0 & device,decomposed=0 & loadct,39 &!P.color=255 & !P.background=0
a= findgen(17) * (!PI*2.0/16.0) & USERSYM, cos(A), sin(A), /FILL
tmax=100
EM_max=1.e8
lambda0=480.6;656.3;471.3;656.3;480.6;587.562;486.133;656.3;486.133;468.57;486.133;nm
mass=39.95;4.;12.01;39.95;1.
resolution=-1.420387E-12*lambda0^3 - 2.156031E-09*lambda0^2 + 1.250038E-06*lambda0 + 3.830769E-03;0.0037714; -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
cd, "C:\Users\Moe Akimitsu\onedriveb\OneDrive - The University of Tokyo\work\2019JSPF\Doppler"
z=read_ascii("z.txt") & z=reform(z.field1)*1.e-3
filename=dialog_pickfile(path="\\FOURIER\md0\Doppler\288CH\20190608",filter="*.asc")
d=read_ascii(filename[0]) & d=transpose(d.field0001(1:*,0:*))
bg=read_Ascii("\\192.168.1.140\md0\makimitsu\bg_gain3000.asc") & bg=transpose(bg.field0001(1:*,0:*)) & D=D-Bg
calib=read_ascii("Ar.calibration.0604.txt",data_start=1) & calib=calib.field1
CH=reform(calib[0,*])
Center=reform(calib[1,*])
Smile=reform(calib[2,*])
relative=reform(calib[4,*])
instru=reform(calib[5,*]) 

separation=where((CH-1) mod 16 eq 1)-1 & separation[3]=separation[3]+1
p=read_Ascii("r.txt") & p=p.field1*1.e-3
edge=0.3
Ti_instru=1.69e8*mass*(2.*resolution*instru*sqrt(2.*alog(2.))/lambda0)^2

lambda=dblarr([1024,n_elements(CH)])
x=dindgen(1024)
for i=0,n_elements(CH)-1 do begin
lambda[*,i]=(x-smile[i])*resolution+lambda0;+0.13
endfor
spectra=fltarr([1024,n_elements(CH)])
!p.MULTI=[0,1,2]
window,0,xsize=800,ysize=1000
loadct,1
interval=7
contour,d,x,x,/fill,nlevels=16,xst=1,yst=1,zst=1,zr=[min(d),(max(d)-min(d))*0.1+min(d)],yr=[min(center)-interval-5,max(center)+interval+5],position=[0.1,0.4,0.95,0.95]
loadct,39
for i=0,n_elements(CH)-1 do begin
oplot,smile[i]+[-40,40],center[i]*[1,1],color=250 & oplot,smile[i]*[1.,1.],center[i]+[-7,7],color=100 &
oplot,smile[i]+[-40,40],center[i]*[1,1]+interval,color=200,linestyle=2 & oplot,smile[i]+[-40,40],center[i]*[1,1]-interval,color=200,linestyle=2
endfor

xr=[410,540]
lambdaA=(x-(xr[0]+xr[1])/2.)*resolution+lambda0 & lambdaB=lambdaA & lambdaA=lambdaA(xr[0]:xr[1])
ybin=TOTAL(d,2)
PLOT,x,ybin,XST=1,YST=1,PSYM=-1,position=[0.1,0.05,0.95,0.4]
oplot,x(xr[0]:xr[1]),gaussfit(x[xr[0]:xr[1]],ybin[xr[0]:xr[1]],coeff,nterm=5),color=250,thick=2
A=[coeff[0],200,coeff[2],coeff[0],350,coeff[2],coeff[0],500,coeff[2],coeff[0],650,coeff[2],coeff[0],800,coeff[2],coeff[3],0,0]
fita=A*0.+1.
;stop
;**************************
;*Line-integrated analysis*
;**************************

passive_Ti=CH*0.
passive_Timax=CH*0.
passive_Timin=CH*0.
passive_Em=CH*0.

!P.multi=[0,4,4]
window,xsize=1800,ysize=900
for i=0,n_elements(CH)-1 do begin
if CH[i] mod 16 eq 0 then window, CH[i]/16+2,xsize=1800,ysize=900
  for j=0,n_elements(X)-1 do begin
    spectra[j,i]=spectra[j,i]+total(d(j:j,center[i]-interval:center[i]+interval))*relative[i]         ;各チャンネルのガウシアン信号を積分
  endfor
  input=reform(spectra[*,i]) 
  fit=gaussfit(x(round(smile[i])-75:round(smile[i])+75),input(round(smile[i])-75:round(smile[i])+75),coeff,nterms=4,sigma=sigma)
  passive_Ti[i]=1.69e8*mass*(2.*resolution*(coeff[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  passive_Timax[i]=1.69e8*mass*(2.*resolution*(coeff[2]+sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  passive_Timin[i]=1.69e8*mass*(2.*resolution*(coeff[2]-sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  plot,x(round(smile[i])-75:round(smile[i])+75),input(round(smile[i])-75:round(smile[i])+75),xst=1,yst=1,title="CH::" + strcompress(CH[i]),charsize=2,psym=1
  oplot,x(round(smile[i])-75:round(smile[i])+75),fit,color=250
  passive_Em[i]=resolution*total(input(round(smile[i])-75:round(smile[i])+75)-min(smooth(reform(input(round(smile[i])-75:round(smile[i])+75)),20)))
  spectra[*,i]=spectra[*,i]-min(smooth(reform(input(round(smile[i])-75:round(smile[i])+75)),20))
 ;stop
endfor

window,31 & !P.multi=0
plot,CH,passive_Ti,psym=8,xst=1,yst=1,yr=[0,max(passive_Ti)]
errplot,CH,passive_Timin,passive_Timax
for i=0,n_elements(separation)-1 do oplot,CH[separation[i]]*[1,1],[-10000,10000],color=150

Ti2D=fltarr([n_elements(p),n_elements(z)])
Em2D=fltarr([n_elements(p),n_elements(z)])
for i=0,n_elements(z)-2 do begin
  if i mod 2 eq 0 then Ti2D[*,i]=spline(CH(separation[i]:separation[i+1])-CH(separation[i]),passive_Ti(separation[i]:separation[i+1]),indgen(16))
  if i mod 2 eq 1 then Ti2D[*,i]=spline(CH(separation[i]:separation[i+1])-CH(separation[i]),passive_Ti(separation[i]:separation[i+1]),indgen(16)) 
  if i mod 2 eq 0 then Em2D[*,i]=spline(CH(separation[i]:separation[i+1])-CH(separation[i]),passive_Em(separation[i]:separation[i+1]),indgen(16))
  if i mod 2 eq 1 then Em2D[*,i]=spline(CH(separation[i]:separation[i+1])-CH(separation[i]),passive_Em(separation[i]:separation[i+1]),indgen(16))
if i mod 2 eq 1 then Ti2D[*,i]=reverse(Ti2D[*,i]) 
if i mod 2 eq 1 then Em2D[*,i]=reverse(Em2D[*,i])

endfor

Ti2D[*,n_elements(z)-1]=spline(CH(separation[i]:*)-CH(separation[i]),passive_Ti(separation[i]:*),indgen(16)) & Ti2D[*,n_elements(z)-1]=reverse(Ti2D[*,n_elements(z)-1])
Em2D[*,n_elements(z)-1]=spline(CH(separation[i]:*)-CH(separation[i]),passive_Em(separation[i]:*),indgen(16)) & Em2D[*,n_elements(z)-1]=reverse(Em2D[*,n_elements(z)-1])

window,31,ysize=800 & !P.multi=[0,1,2]
contour,transpose(Ti2D)/max(Ti2D),z,p,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1]
loadct,1 & contour,transpose(Em2D)/(max(Em2D)*0.75),z,p,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1];max(Em2D)*0.75]
loadct,39
!P.multi=0
;stop
;******************
;*Abel-ininversion*
;******************
num=250
yy=findgen(num+1)/num*(edge-min(p))+min(p)
spectra_interp=fltarr([n_elements(lambdaA),n_elements(yy),n_elements(z)])

for i=0,n_elements(z)-2 do begin
window,i;,xsize=400,ysize=300
  if i mod 2 eq 0 then result=trigrid_interpor_for_r_lambda([[reform(spectra(0:*,separation[i]:separation[i+1]-1))],[fltarr(n_elements(lambdaB))]],[[reform(lambda[0:*,separation[i]:separation[i+1]-1])],[lambdaB]],[p(CH(separation[i]:separation[i+1]-1)-CH(separation[i])),edge],lambdaA,yy)
  if i mod 2 eq 1 then result=trigrid_interpor_for_r_lambda([[fltarr(n_elements(lambdaB))],[reform(spectra(0:*,separation[i]:separation[i+1]-1))]],[[lambdaB],[reform(lambda[0:*,separation[i]:separation[i+1]-1])]],[edge,reverse(p(CH(separation[i]:separation[i+1]-1)-CH(separation[i])))],lambdaA,yy)
  spectra_interp[*,*,i]=result.z
endfor
window,17
result=trigrid_interpor_for_r_lambda([[fltarr(n_elements(lambdaB))],[reform(spectra(0:*,separation[17]:*))]],[[lambdaB],[reform(lambda[0:*,separation[17]:*])]],[edge,reverse(p(CH(separation[17]:N_ELEMENTS(ch)-1)-CH(separation[17])))],lambdaA,yy)
spectra_interp[*,*,17]=result.z
!p.multi=[0,6,6]
for i=0,num do begin
  if i mod 36 eq 0 then window, fix(i/36),xsize=1700,ysize=900
  plot,lambdaA,reform(spectra_interp[*,i,8]),xst=1,yst=1,title="(R,Z) ="+strcompress(yy[i],/remove_all),charsize=1.5
  oplot,lambdaA,gaussfit(lambdaA,reform(spectra_interp[*,i,8]),coeff,nterms=3,sigma=sigma),color=250,thick=2
endfor 

Local_spectra=fltarr([n_elements(lambdaA),num+1,n_elements(z)])
spectra_interp=smooth(spectra_interp,[5,num/16,1])
dy=yy[1]-yy[0]

for k=0,n_elements(z)-1 do begin
  for l=0,n_elements(lambdaA)-1 do begin
  derivative=deriv(yy,reform(spectra_interp[l,*,k]))
    for i=0,num do begin
      for j=i,num-1 do local_spectra[l,i,k]=local_spectra[l,i,k]-1./!pi*derivative[j]*alog(yy[j+1]*(1.+sqrt(1.-yy[i]^2/yy[j+1]^2))/(yy[j]*(1.+sqrt(1.-yy[i]^2/yy[j]^2))));Balandin's Abel inversion
    endfor
  endfor
endfor
Local_Spectra=smooth(Local_Spectra,[5,num/16,1])

emission=fltarr([num+1,n_elements(z)])
Ti_2D=fltarr([num+1,n_elements(z)])
Ti_max=fltarr([num+1,n_elements(z)])
Ti_min=fltarr([num+1,n_elements(z)])
Ti_instru2=total(Ti_instru(separation[0]:separation[1]))/(separation[1]-separation[0])
for i=1,n_elements(z)-2 do Ti_instru2=[Ti_instru2,total(Ti_instru(separation[i]:separation[i+1]))/(separation[i+1]-separation[i])]
Ti_instru2=[Ti_instru2,total(Ti_instru(separation[i]:n_elements(CH)-1))/(n_elements(CH)-1-separation[i])]

!P.multi=[0,6,6]
for i=0,num do begin
;if i mod 36 eq 0 then window, fix(i/36),xsize=1700,ysize=900
  for j=0,n_elements(z)-1 do begin
  input=reform(Local_spectra[*,i,j]) ;&  if input lt 0 then input[*] = -input[*] ;-min(smooth(reform(Local_spectra[*,i,j]),20))
  fit = gaussfit(lambdaA,input,coeff,nterms=3,sigma=sigma)
 ; stop
; plot,lambdaA,input,xst=1,yst=1,title="(R,Z) ="+strcompress(yy[i],/remove_all)+","+strcompress(z[j],/remove_all),charsize=1.5
; oplot,lambdaA,fit,color=250,thick=2
  Ti_2D[i,j]=1.69e8*MASS*(2.*coeff[2]*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  Ti_max[i,j]=1.69e8*MASS*(2.*(abs(coeff[2])+3.*abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  Ti_min[i,j]=1.69e8*MASS*(2.*(abs(coeff[2])-3.*abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  emission[i,j]=total(input*resolution);coeff[0]
  checker=float(abs(coeff[1]-lambda0) lt 0.3)*(coeff[0] gt 0)*(emission[i,j] gt 100)*float(abs(Ti_max[i,j]-Ti_min[i,j]) lt Ti_2D[i,j]);*float(abs(Ti_max[i,j]-Ti_min[i,j]) lt Ti_2D[i,j]+Ti_instru2[j]);*(emission_local[i,j] gt EM_max*0.1)
  Ti_2D[i,j]=Ti_2D[i,j]*checker & Ti_max[i,j]=Ti_max[i,j]*checker & Ti_min[i,j]=Ti_min[i,j]*checker
  endfor
endfor

;remove "NAN" grid point
for i=0,n_elements(z)-1 do begin
  for j=0,n_elements(yy)-1 do begin
    if finite(Ti_2D[j,i]) ne 1 then begin
    Ti_2D[j,i]=0.
    Ti_max[j,i]=0.
    Ti_min[j,i]=0.
    print,i,j   
    endif 
  endfor
endfor


start=where((yy gt p[0])) & start=start[0]
;Ti_2D=Ti_2D*float(emission_local gt 0)
;Ti_2D=Ti_2D*float(emission_local gt EM_max*0.04)
;for i=0,n_elements(z)-1 do Ti_2D[*,i]=Ti_2D[*,i]*float(abs(Ti_max[*,i]-Ti_min[*,i]) lt Ti_2D[*,i]+Ti_instru2[i])
checker=Ti_2D le 0
for i=start,n_elements(yy)-2 do begin
    for j=0,n_elements(z)-1 do begin
      if checker[i,j] then begin
        Ti_2D[i,j]=sqrt((Ti_2D[i+1,j]+Ti_2D[i-1,j])/2.)
        Ti_max[i,j]=(Ti_max[i+1,j]+Ti_max[i-1,j])/2.  
        Ti_min[i,j]=(Ti_min[i+1,j]+Ti_min[i-1,j])/2.        
      endif
    endfor  
endfor
loadct,39
window,26,ysize=800 & !P.multi=[0,1,2]
contour,Ti_2D/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1],xr=[min(p),max(p)]
loadct,1 &contour,emission/max(emission*0.75),yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1],xr=[min(p),max(p)];max(Em2D)*0.75]

window,27,ysize=500,xsize=1400
!P.multi=[0,4,1]
loadct,39 & contour,transpose(Ti2D),z,p,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,tmax],xtitle="Z[m]",ytitle="R[m]",title="Ti(Projection)";,/isotropic
plot,Ti2D[*,0],p,xst=1,yst=1,xr=[0,tmax*2],/nodata,xtitle="Ti [eV]",ytitle="R[m]",title="Ti(Projection)" & for i=0,5 do oplot,Ti2D[*,i],p,color=i*50,psym=-8
loadct,39 & contour,transpose(Ti_2D),z,yy,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,tmax],xtitle="Z[m]",ytitle="R[m]",title="Ti(Local)",yr=[min(p),max(p)];,/isotropic
plot,Ti_2D[*,0],yy,xst=1,yst=1,xr=[0,tmax*2],/nodata,xtitle="Ti [eV",ytitle="R[m]",title="Ti(Local)",yr=[min(p),max(p)] & for i=0,17 do oplot,Ti_2D[*,i],yy,color=i*15,psym=-8

window,28,xsize=1800,ysize=700 & loadct,5
!P.multi=[0,4,5]
for i=0,n_elements(z)-1 do begin
contour, Local_spectra[*,*,i],lambdaA,YY,/fill,nlevels=32,xst=1,yst=1,zst=1,zr=[0,max(local_spectra)]
oplot,lambda0*[1,1],[0,1000],linestyle=2,color=125,thick=2
endfor

;data=MAG_CAO(EF=-150,m_num=4,n_num=4,n_frame=16,dateshot=190605010)
restore,"\\192.168.1.140\md0\makimitsu\fordoppler2174.sav"
window,30,xsize=1200,ysize=600 & !P.multi=[0,2,1] & loadct,39
contour,transpose(Ti_2D)/tmax,z,yy,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="Z [m]",ytitle="R [m]",title="Ti (Local)",color=255,background=0;yr=[0.15,max(p)],xr=[-0.025,0.025],
;for i=0,n_elements(z)-1 do oplot,fltarr(n_elements(p))+Z[i],p,psym=8,color=250
;contour,data.A[13:28,*,1]/max(abs(data.A[13:28,*,1])),data.b[13:28],data.c,/overplot,nlevels=50,color=0,thick=2.0
;contour,-data.A[13:28,*,1]/max(abs(data.A[13:28,*,1])),data.b[13:28],data.c,/overplot,nlevels=50,color=0,thick=2.0
contour,(psi[*,*,12]/max(abs(psi))+1)/2,xout,yout,/overplot,nlevels=100
loadct,5 & contour,transpose(emission)/max(emission*0.75),z,yy,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="Z [m]",ytitle="R [m]",title="Emission (Local)",color=255,background=0;,yr=[0.15,max(p)],xr=[-0.025,0.025]
contour,(psi[*,*,12]/max(abs(psi))+1)/2,xout,yout,/overplot,nlevels=100
loadct,39;& for i=0,n_elements(z)-1 do oplot,fltarr(n_elements(p))+Z[i],p,psym=8,color=150

;contour,data.A[13:28,*,1]/max(abs(data.A[13:28,*,1])),data.b[13:28],data.c,/overplot,nlevels=50,color=0,thick=2.0
;contour,-data.A[13:28,*,1]/max(abs(data.A[13:28,*,1])),data.b[13:28],data.c,/overplot,nlevels=50,color=0,thick=2.0
stop

contour,smooth(Ti_2D,[num/16,2])/tmax,yy,z,/fill,nlevels=256,/isotropic,title="Ti::"+filename,xst=1,yst=1,zst=1,zr=[0,1],xtitle="R [m]",ytitle="Z [m]",xr=[min(p),max(p)],charsize=1.0,color=0,background=16777215L;,color=16777215L
for i=0,n_elements(z)-1 do oplot,p,fltarr(n_elements(p))+Z[i],psym=8,color=250
contour,transpose(data.A[62:143,2:194,10])/max(abs(data.A[62:143,2:194,10])),data.c[2:194],data.b[62:143],/overplot,nlevels=30,color=0,thick=2.0
contour,-transpose(data.A[62:143,2:194,10])/max(abs(data.A[62:143,2:194,10])),data.c[2:194],data.b[62:143],/overplot,nlevels=30,color=0,thick=2.0

;stop
data=mrd_analysis(150.)
window,31 ,ysize=1000,xsize=700 & !P.multi=0 
contour,transpose(data.A[*,*,80])/max(abs(data.A[*,*,80])),data.c,data.b,nlevels=30,xst=1,yst=1,zst=1,color=0,xtitle="R [m]",ytitle="Z [m]",thick=2.0,background=16777215L
contour,-transpose(data.A[*,*,80])/max(abs(data.A[*,*,80])),data.c,data.b,nlevels=30,/overplot,color=0,thick=2.0

colorbar=[[findgen(256)],[findgen(256)]]
x=findgen(256)/256*tmax
y=[0,1]
window,ysize=400
contour,colorbar,x,y,/fill,nlevel=256,xst=1,yst=7,zst=1,title="Ti[eV]",position=[0.1,0.85,0.9,0.9],color=0,background=255
;cb=colorbar(target=smooth(Ti_2D,[num/16,2])/tmax,orientation=0
;colorbar =[[findgen(255)-127.)/127.],[findgen(255)-127.)/127.]]
;contour,colorbar,(findgen(255)-127.)/127.,[0,1],/fill,nlevels=32,/noerase,position=[0.1,0.85,0.9,0.9],title="Ti::"+filename,color=white,xst=1,yst=4
stop
;
;stop
;save,filename=filename+".sav",Ti_2D,yy,z,emission_local
;kannkeinaityo
;save,filename=strmid(filename,0,48)+".sav",Ti_2D,yy,z,emission_local
;screen_output,filename+".png"
stop
END
