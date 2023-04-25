PRO IDL_2Dcontour
device,decomposed=0
loadct,39
test=dist(256)
;contour,test
;contour,test,nlevels=32
contour,test,nlevels=256,/fill,xstyle=1,ystyle=1,zstyle=1
contour,test,nlevels=32,/overplot
END