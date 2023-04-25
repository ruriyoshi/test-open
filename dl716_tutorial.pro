PRO DL716_tutorial
filename=dialog_pickfile(path="K:\results\ts-3u\170907",filter="*.rgw")
data=read_ascii(filename,data_start=1) & data=data.field01(1:*,0:*)
time=data[0,*]
N_CH=15
window,10 & !P.multi=[0,4,4]
for i=0,15 do plot,time,data[i,*],xst=1,yst=1,xtitle="time [us]",ytitle="signal [V]",charsize=2
END