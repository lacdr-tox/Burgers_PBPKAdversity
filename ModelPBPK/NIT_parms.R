Kgut_plasma = 0.622;		 	   #Liver/blood partition coefficient
Kliver_plasma = 0.651;		 	 #Liver/blood partition coefficient  
Kkidney_plasma = 0.671;		   #kidney/blood partition coefficient
Kfat_plasma = 0.159;	       #Fat/blood partition coefficient  (fitted)
Krestbody_plasma = 0.423;    #Rest of the body/blood partition coefficient	(fitted)

fu = 0.42; 
kbp = 0.76;
kgutabs = 0.5275
kfeces = 0.0187;
Kt = 0.059
kbile = 0.063

if(Species == "Human"){
  fugut = fu
  QurineC = 4.98
  
  Vehrc = 0.00182
  Kehr = 0.0042
  
  Trc = 0.00826
  Tmc = 0.0498
  
  
  VmaxMC = 204.32
  KmM = 5.83
}else{
  fugut = 0.006

  QurineC = 20.36
  
  Vehrc = 0.52
  Kehr = 0.017

  Trc = 2.365
  Tmc = 14.26
  
  VmaxMC = 90.97
  KmM = 11.95
}


