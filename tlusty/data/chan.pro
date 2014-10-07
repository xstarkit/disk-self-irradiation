pro chan,iat,ionmax
;
ion=indgen(ionmax)+1
print,ion
for i=0,ionmax-1 do begin
file=iat+strtrim(string(ion(i)),2)
a='/bin/cp -f '+file+' '+file+'_old'
print,a
spawn,a
a='sed ''s/ 1 15 4 0/ 1 15 99 0/''
a=a+' '+file+'_old >! '+file
print,a
spawn,a
endfor
;
return
end
