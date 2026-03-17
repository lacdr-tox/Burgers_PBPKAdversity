Kgut_plasma = 6.48;		 	     #gut/plasma partition coefficient
Kliver_plasma = 5.17;		 	   #Liver/blood partition coefficient  
Kkidney_plasma = 3.21;		   #kidney/blood partition coefficient
Kfat_plasma = 22.68;	       #Fat/blood partition coefficient  
Krestbody_plasma = 3.25;     #Rest of the body/blood partition coefficient	assumed same as fat

fu = 0.01;   
kbp = 0.6;  

QurineC = 0.001
kfeces = 0.3

kgutabs = 10

if(Species == "Human"){
  fugut = 0.0967
  VmaxMC = 999.98
  KmM = 9.865
  
}else{
  fugut = 0.01377208
  VmaxMC = 507.8474727
  KmM = 6.536582192
}

if (gender =="female"){
  VmaxMC = VmaxMC/2
}

# not relevant for KET
Vehrc = 0
Kehr = 1
kbile = 0

Trc = 0
Tmc = 0
Kt = 1

