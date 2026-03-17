BW = 0.25

#constant Fraction of blood flows to organs (blood flow rate)
QCC = 17.76;                  #Total Cardiac blood output (L/h/kg)  
FQliver = 0.174; 		          #Fraction cardiac output going to liver 
FQkidney = 0.141;  		        #Fraction cardiac output going to kidney  
FQfat =  0.07;  		          #Fraction cardiac output going to fat  
FQgut = 0.1;   			          #Fraction cardiac output going to gut

Qfilterate= 0.0786 ;				    #Glomerulation filteration rate L/hr equivalent to 120 ml/min

#constant organ volume as a fraction of total body weight
Fliver = 0.0366; 		          #Fraction liver volume
Fkidney = 0.0073; 			      #Fraction kidney volume 
Ffilterate=0.00073;           #10 percent of kidney volume
Ffat = 0.0721;                #fractional volume of fat  
Fgut = 0.027;   				      #fractional volume of gut
Fplasma = 0.0321  	          #fractional volume of plasma

if(gender == "female"){
  BW = BW*0.8
}