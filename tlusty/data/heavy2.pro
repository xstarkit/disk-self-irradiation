pro heavy2
;
ene=[ 173.93, 330.391, 511.8,   783.3,  1018.,  1273.8,  1671.792, $
     1928.462,9645.005,10986.876]
;
emg=[  61.671, 121.268, 646.41, 881.1, 1139.4, 1504.3, 1814.3, 2144.7, $
     2645.2,  2964.4, 14210.261,15829.951]
;
esi=[65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4, $
     2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63]
;
es= [83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4, $
     3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663, $
     28182.535]
;
ear=[127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08, $
     3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4, $
     33237.173,35699.936]
;
eca=[49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6, $
     1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860., $
     8770.,9338.,41366.,44177.41]
;
efe=[63.737,130.563,247.22,442.,605.,799.,1008.,1218.38, $
     1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82, $
     10180.,10985.,11850.,12708.,13620.,14510.,15797., $
     16500.,71203.,74829.6]
;
eni=[61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560., $
     1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606., $
     4896.2,12430.,13290.,14160.,15280.,16220.,17190., $
     18510.,19351.,82984.,86909.4]
;
gli=[2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.]
gfe=[2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1., $
     10.,21.,28.,25.,6.,25.,30.,25.]
gni=[2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1., $
     10.,21.,28.,25.,6.,25.,28.,21.,10.,21.]
;
cons=0.1239529*3.28805e15/13.595
;cons=0.1239529
;
ene=ene*cons
emg=emg*cons
ensi=esi*cons
es=es*cons
ear=ear*cons
eca=eca*cons
efe=efe*cons
eni=eni*cons
;
a1='****** Levels'
;
; 1.55935E+16      2.    2  'c4   2Se 1'  0   0.  0
;
lab1=['ne','mg','si','s','ar','ca','fe','ni']
lab0=['neon','magnesium','silicon','sulfur','argon','calcium','iron','nickel']
at=[10,12,14,16,18,20,26,28]
;
lb1="    3  ' "
lb2="    '  0   0.  0"
;
close,4
openw,4,'dat.tmp'
for i=0,7 do begin
  if i eq 0 then e=ene $
  else if i eq 1 then e=emg $
  else if i eq 2 then e=esi $
  else if i eq 3 then e=es $
  else if i eq 4 then e=ear $
  else if i eq 5 then e=eca $
  else if i eq 6 then e=efe $
  else if i eq 7 then e=eni 
  if i le 5 then g=reverse(gli(0:n_elements(e)-1)) $
  else if i eq 6 then g=reverse(gfe) $
  else if i eq 7 then g=reverse(gni) 
;
  for j=0,n_elements(e)-2 do begin
    lab2=lab1(i)+strtrim(string(j+1),2)
    print,format='(e12.5,f8.0,a9,a5,a16)',e(j),g(j),lb1,lab2,lb2
    openw,3,'tmp'
    printf,3,format='(a13)',a1
    printf,3,format='(e12.5,f8.0,a9,a5,a16)',e(j),g(j),lb1,lab2,lb2
    close,3
    a='cat tmp /stis9/blaes/atoms/pho*/heavies/'+lab0(i)+strtrim(string(j+1),2) 
    a=a+' >! '+lab2
    spawn,a
    file=lab2
     a='/bin/cp -f '+file+' '+file+'_old'
;     print,a
     spawn,a
     a='sed ''s/ 1 15 4 0/ 1 15 99 0/''
     a=a+' '+file+'_old >! '+file
;     print,a
     spawn,a
;
     lab3=strtrim(lab2,2)
     if strlen(lab3) eq 2 then lab3='  '+lab3
     if strlen(lab3) eq 3 then lab3=' '+lab3
     lab3=lab3+"' '/stis19/hubeny/omer2/data/"+lab2+"'"
     printf,4,format='(2i4,2a)',at(i),j,"    1     0     0     0   '",lab3
  endfor 
     lab2=lab1(i)+strtrim(string(n_elements(e)),2)
     lab3=strtrim(lab2,2)
     if strlen(lab3) eq 2 then lab3='  '+lab3
     if strlen(lab3) eq 3 then lab3=' '+lab3
     lab3=lab3+"'     ' '"
 printf,4,format='(2i4,2a)',at(i),n_elements(e)-1,$
           "    1     1     0     0   '",lab3
endfor
close,4
;
return
end
  

