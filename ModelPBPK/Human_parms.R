BW = 70

#constant Fraction of blood flows to organs (blood flow rate)
QCC = 4.8;                  #Total Cardiac blood output (L/h/kg)  
FQliver = 0.25; 		        #Fraction cardiac output going to liver 
FQkidney = 0.19;  		      #Fraction cardiac output going to kidney  
FQfat =  0.05;  		        #Fraction cardiac output going to fat  
FQgut = 0.196;   			      #Fraction cardiac output going to gut

Qfilterate = 7.5 ;				  #Glomerulation filteration rate L/hr equivalent to 125 ml/min

#constant organ volume as a fraction of total body weight
Fliver = 0.0257; 		        #Fraction liver volume
Fkidney = 0.0044; 			    #Fraction kidney volume  
Ffilterate= 0.00044;		    #Fraction tubules volume 
Ffat = 0.213;               #fractional volume of fat
Fgut = 0.0165;   				    #fractional volume of gut
Fplasma = 0.0428  	        #fractional volume of plasma

if(gender == "female"){
  BW = 58
  Ffat = 0.327
}

if(elderly){
  Qfilterate = 3.6 #Noronha 60 ml/min
  QCC = QCC*(0.99^40) #Butler
  Fliver = 0.016 #Butler
  Fkidney = 0.75*Fkidney #Drenth-van Maanen
  Ffat = ifelse(gender == "male", 0.35, 0.49) #Forbes
}

# simulating severe GFR
if(severe_GFR) {
  Qfilterate = 1.2 # 20 ml/min
}
