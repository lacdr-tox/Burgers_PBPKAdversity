Kgut_plasma = 4.10;		 	    #Liver/blood partition coefficient
Kliver_plasma = 3.31;		 	  #Liver/blood partition coefficient
Kkidney_plasma = 2.11;		  #kidney/blood partition coefficient
Kfat_plasma = 0.10;	        #Fat/blood partition coefficient
Krestbody_plasma = 2.14;    #Rest of the body/blood partition coefficient

fu = 0.003
fugut = 0.043322627

kbp = 0.605

QurineC = 1
kfeces = 0.0001

kgutabs = 10
if(Species == "Human"){
  VmaxMC = 7.8087
  KmM = 1.42E-05
}else{
  VmaxMC = 312.8816
  KmM = 0.083
}

# not relevant for DIC
Vehrc = 0
Kehr = 1
kbile = 0

Trc = 0
Tmc = 0
Kt = 1
