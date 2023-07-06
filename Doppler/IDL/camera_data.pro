Function camera_data, img=img
d=read_ascii(filename[0]) & d=transpose(d.field0001(1:*,0:*))
  for i=0,n_elements(CH)-1 do begin
    lambda[*,i]=(x-smile[i])*resolution+lambda0;+0.13
  endfor
End