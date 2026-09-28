// 6/18/2014: Hexc, Olesya: Making major upgrade of the detector contruction 
//            by using template vectors for describing atm layers.
//
// 09/05/2014: Hexc, Olesya: We are continuing modified the code and removing some of the G4cout's
// 11/06/2014: Hexc, Olesya: We fixed a bug in the constrcutor, i.e., atmType and denPer were NOT initialized
//            properly.
// 11/10/2014: Hexc, Olesya: Fixed a bug in the Layer character array.
// 12/04/2014: Hexc, Olesya: Change the material type to water for the earth.
// 04/06/2015: Hexc, Olesya: update the molecular composition of atmospheric air
//                          http://en.wikipedia.org/wiki/Atmosphere_of_Earth
//  
#include "ECRSDetectorConstruction.hh"
#include "ECRSDetectorMessenger.hh"
#include "ECRSAtmosphereSD.hh"
#include "ECRSMagneticField.hh"

//#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
//#include "PhysicalConstants.h"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

using namespace CLHEP;

ECRSDetectorConstruction::ECRSDetectorConstruction()
  : NbOfLayers(0), Si(NULL), space(NULL), elN(NULL), elO(NULL), elAr(NULL), earth(NULL),
   elC(NULL),N2(NULL),O2(NULL),Ar(NULL),CO2(NULL),solidUniverse(NULL),
   logicUniverse(NULL),physiUniverse(NULL),ECRS_sphere(NULL),
    ECRS_log(NULL),ECRS_phys(NULL), atmosphereSD(0), atmType(0), denPer(1.0) 
{
  NbOfLayers  = NUMBER_LAYER;
  mult = 1.0;
  multDensity = 1.0;
  visibility = true;
  externalMag = true;
  //externalMag = false;
    
  // Set detector sizes
  // set the radius of the Earth
  earthRadius = 6.3712e6*m; 
  oceanDepth = 11.e3*m;  

  universeSize = 3*earthRadius;   // world size is three times Earth's radius

  upperATMHeight = 1.e07*cm;
  airLayerThickness = 1.0e05*cm;

  G4cout << "Basic runconfiguration: " << " atmType = " << atmType << "density scale factor = " << denPer << G4endl;

  //  G4cout << "Build ECRS Detector Messenger ... " << G4endl;
  ECRSDetector = new ECRSDetectorMessenger(this);  
  //G4cout << "Build ECRS Detector Messenger ... Done!!! " << G4endl;
  
  if (  externalMag == false )
    {
      G4cout << "External Magnetic Field is off " <<  externalMag << G4endl;
    }
  
  if ( externalMag == true){
    theMagneticField=new ECRSMagneticField();
    G4cout << "External Magnetic Field is on " <<  externalMag << G4endl;
  }

}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

ECRSDetectorConstruction::~ECRSDetectorConstruction()
{
  delete ECRSDetector;
  delete theMagneticField;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

G4VPhysicalVolume* ECRSDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructUniverse();
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::DefineMaterials()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  
  //Define the ECRS to be made of solid silicone
  a = 28.086*g/mole;
  density = 2.4*g/cm3;
  //  density *= multDensity;
  Si = new G4Material(name="Silicone",z=14., a, density);
  
  // G4cout << " Multy density " << multDensity << G4endl;
  
  //Define elements and gasses to compose atmosphere
  a = 14.007*g/mole;
  elN = new G4Element(name="Nitrogen",symbol="N",z=7.,a);
  
  a = 15.999*g/mole;
  elO = new G4Element(name="Oxygen",symbol="O",z=8.,a);
  
  a = 39.948*g/mole;
  elAr = new G4Element(name="Argon",symbol="Ar",z=18.,a);
  
  a = 12.011*g/mole;
  elC = new G4Element(name="Carbon",symbol="C",z=6.,a);
  
  density = 0.001251*g/cm3;
  //  density *= multDensity;
  N2 = new G4Material(name="n_gas",density,ncomponents=1);
  N2->AddElement(elN,natoms=2);
  
  density = 0.001429*g/cm3;
  //  density *= multDensity;
  O2 = new G4Material(name="o_gas",density,ncomponents=1);
  O2->AddElement(elO,natoms=2);
  
  density = 0.001784*g/cm3;
  //  density *= multDensity;
  Ar = new G4Material(name="Ar_gas",density,ncomponents=1);
  Ar->AddElement(elAr,natoms=1);
  
  density = 0.001965*g/cm3;
  //  density *= multDensity;
  CO2 = new G4Material(name="carbon_dioxide",density,ncomponents=2);
  CO2->AddElement(elC,natoms=1);
  CO2->AddElement(elO,natoms=2);
  
  // Define air layer material
  
  // Need to define the air-density for each airlayer, i.e., fill the LayerDensity vector
  // Maximum number of layers should be less than 99.  Otherwise, one need to increase the array from 10 to 11 
  // for the Layer array.

  char Layer[11];
  // G4Material* AirLayer[NUMBER_LAYER];
  G4double height, densityValue, temperatureValue, pressureValue;
 
  for (int iLayer = 0; iLayer<NbOfLayers; iLayer++)
    {
      // Make volume
      sprintf(Layer, "AirLayer%02d", iLayer);
      height = upperATMHeight - iLayer * airLayerThickness;
      densityValue = GetAirDensity(height,0);
      //Xiaohang added following two lines to get temperature and pressure
      temperatureValue = GetAirDensity(height,2)+273.1;
      pressureValue = GetAirDensity(height,1)*1000;
      
      //G4cout << "Height: " << height << "   Density: " << densityValue 
      //     << "   Pressure: " << pressureValue 
      //     << "   Temperature: " << temperatureValue<< G4endl;

      //fileOut->fout  << "Height: " << height << "   Density: " << densityValue << G4endl;
      //AirLayer[iLayer] = new G4Material(Layer, density= densityValue*mg/cm3, ncomponents=2);
      //AirLayer[iLayer] = new G4Material(Layer, ncomponents=2, kStateGas, temperatureValue *kelvin, pressureValue*pascal);

      G4Material *AirLayer = new G4Material(Layer, density= densityValue*mg/cm3, 
					ncomponents=4, kStateGas, 
					temperatureValue *kelvin, 
					pressureValue*pascal);
      //      AirLayer->AddElement(elN, fractionmass=0.7);
      //      AirLayer->AddElement(elO, fractionmass=0.3);
      AirLayer->AddMaterial(N2, fractionmass=0.780870);  
      AirLayer->AddMaterial(O2, fractionmass=0.209476);
      AirLayer->AddMaterial(Ar, fractionmass=0.009340);
      AirLayer->AddMaterial(CO2,fractionmass=0.000314); // 1962 version of US atm standard 
      airLayer_Matt.push_back(AirLayer);
    }

  //Make the vacuum of space
  density = universe_mean_density;  //included from PhysicalConstants.h
  pressure = 1.0E-19*pascal;
  temperature = 2.74*kelvin;
  space = new G4Material(name="space",z=1.,a=1.01*g/mole,
			 density,kStateGas,temperature,pressure);
  
  //Make the Univers
  density = .0012250*g/cm3;  //included from PhysicalConstants.h
  pressure = 1.01325E5*pascal;
  temperature = 288.150*kelvin;
  earth = new G4Material(name="earth",density,ncomponents=4,
			 kStateGas,temperature,pressure);
  earth->AddMaterial(N2,fractionmass=75.521*perCent);
  earth->AddMaterial(O2,fractionmass=23.143*perCent);
  earth->AddMaterial(Ar,fractionmass=1.288*perCent);
  earth->AddMaterial(CO2,fractionmass=0.048*perCent);

  // Water is defined from NIST material database
  G4NistManager *man = G4NistManager::Instance();
  H2O = man->FindOrBuildMaterial("G4_WATER");
  
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* ECRSDetectorConstruction::ConstructUniverse()
{
  // THE UNIVERSE (i.e. World)

  // Define the angles for spheres of atmosphere and ECRS
  G4double startAnglePhi = 0.0*deg;
  G4double spanningAnglePhi = 360.0*deg;
  G4double startAngleTheta = 0.0*deg;
  G4double spanningAngleTheta = 180.0*deg; //90

  solidUniverse = new G4Sphere("World", 0.0*m, universeSize,
			       startAnglePhi, spanningAnglePhi, 
			       startAngleTheta, spanningAngleTheta);

  logicUniverse = new G4LogicalVolume(solidUniverse, 
				      space,             // material 
				      "World",
				      0,
				      0,
				      0);
  
  physiUniverse = new G4PVPlacement(0,                    // no rotation
				    G4ThreeVector(),      // put at (0,0,0)
				    "World",           // name
				    logicUniverse, 
				    0,                    // no mother
				    false, 
				    0);                   // copy number
  
  //                                 
  // Atmosphere Layers
  //
  
  char Layer[8];
  G4double atmLayerHeight = upperATMHeight, atmLayerHeight_prev;
  for (int iLayer = 0; iLayer<NbOfLayers; iLayer++)
    {
      // Make volume
      sprintf(Layer, "Layer%02d", iLayer);   // Maximum number of layer should be less than 99!!!
      atmLayerHeight  = upperATMHeight - iLayer*airLayerThickness;
      //G4cout << "************ atmLayerHeight " << atmLayerHeight << G4endl;
      airLayer.push_back(new G4Sphere(Layer, 0.0*m, earthRadius+atmLayerHeight,
			       startAnglePhi, spanningAnglePhi, 
			       startAngleTheta, spanningAngleTheta));
      
      // Define logical volume for the volume
      //sprintf(Layer_logic, "Logic_Layer%02d", iLayer);
      airLayer_logic.push_back(new G4LogicalVolume(airLayer[iLayer], airLayer_Matt[iLayer], Layer));
      
      if (iLayer == 0)
	{
	  airLayer_phys.push_back(new G4PVPlacement(0,	                     // no rotation
						    G4ThreeVector(),         // at (0,0,0)
						    Layer,		     // its name
						    airLayer_logic[iLayer],  // its logical volume     
						    physiUniverse,	     // its mother  volume
						    false,		     // no boolean operation
						    0));		     // copy number
	  atmLayerHeight_prev = atmLayerHeight;
	}
      else 
	{
	  airLayer_phys.push_back(new G4PVPlacement(0,	                     // no rotation
						    G4ThreeVector(),	     // at (0,0,0)
						    Layer,		     // its name
						    //	airLayer_logic[iLayer],  //its logical volume     
						    //	airLayer_phys[iLayer-1], //its mother  volume
						    airLayer_logic.at(iLayer),	 //its logical volume	     
						    airLayer_phys.at(iLayer-1),	 //its mother  volume
						    false,		     //no boolean operation
						    0));		     //copy number
	  atmLayerHeight_prev = atmLayerHeight;
	}
    } 

  // Set the universe color attribute
  G4VisAttributes* simpleWorldVisAtt= new G4VisAttributes(G4Colour(1., 0., 0.));   // red color
  simpleWorldVisAtt->SetVisibility(false);
  simpleWorldVisAtt->SetForceWireframe(false);
  //simpleWorldVisAtt->SetForceSolid(false);

  logicUniverse->SetVisAttributes(simpleWorldVisAtt);

  // Define atm layer color attributes
  G4VisAttributes* simpleLayerVisAtt[NUMBER_LAYER];
  
  for (int iLayer = 0; iLayer < NbOfLayers; iLayer++)
    {
      // simpleLayerVisAtt[iLayer] = new G4VisAttributes(G4Colour(1.0/(iLayer+1),1.0/(iLayer+1),1.0/(iLayer+1)));
      simpleLayerVisAtt[iLayer] = new G4VisAttributes(G4Colour(0.1, 0.1, (0.2 + 0.01*iLayer)));
      simpleLayerVisAtt[iLayer]->SetVisibility(false);
      //simpleLayerVisAtt[iLayer]->SetForceSolid(false);
      //simpleLayerVisAtt[iLayer]->SetForceWireframe(true);
      (airLayer_logic[iLayer])->SetVisAttributes(simpleLayerVisAtt[iLayer]);
    }
    
  // THE EARTH Water Layer
  G4double outerRadiusECRS = mult*earthRadius;
  G4double innerRadiusECRS = outerRadiusECRS - oceanDepth;
  ECRS_sphere = new G4Sphere("ECRS_earth",
			     innerRadiusECRS,outerRadiusECRS,
			     startAnglePhi,spanningAnglePhi,
			     startAngleTheta,spanningAngleTheta);
  ECRS_log = new G4LogicalVolume(ECRS_sphere,H2O,"ECRS_earth");   // changed to H2O from air
  ECRS_phys = new G4PVPlacement(0,               // rotation
				G4ThreeVector(), // position
				"ECRS_earth",          // name
				ECRS_log,        // logic name 
				// physiUniverse,   // its mother logic volume
				airLayer_phys.at(NUMBER_LAYER - 1),  // its mother logic volume
				false,
				0);
  
  G4VisAttributes* ECRSAtt = new G4VisAttributes(G4Colour(0., 0., 1.));  // blue color
  ECRSAtt->SetVisibility(true);
  //ECRSAtt->SetForceSolid(true);
  ECRSAtt->SetForceWireframe(true);
  ECRS_log->SetVisAttributes(ECRSAtt);
  
  // Set the atmosphere as the Sensitive Detector
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if (!atmosphereSD)
    {
      atmosphereSD = new ECRSAtmosphereSD("AtmoSD",this);
      SDman->AddNewDetector(atmosphereSD);
    }

  PrintUniverseParameters();
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  return physiUniverse;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::PrintUniverseParameters()
{
  G4cout << G4endl;
  G4cout << "Universe radius: ";
  G4cout << G4BestUnit(universeSize,"Length") << G4endl;
  G4cout << "Earth's radius: ";
  G4cout << G4BestUnit(mult*earthRadius,"Length") << G4endl;
}
  
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::SetMultiplier(G4double newValue)
{
  mult = newValue;
  G4cout << "You must use the command: /geometry/update before beamOn";
  G4cout << G4endl;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::SetDensityMult(G4double newValue)
{
  multDensity = newValue;
  G4cout << "You must use the command: /geometry/update before beamOn";
  G4cout << G4endl;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::SetVisFlag(G4bool newValue)
{
  /*
    Atm37Att->SetVisibility(newValue);
    atm37_log->SetVisAttributes(Atm11Att);
  */

  G4cout << "You may have to use: /vis/viewer/refresh to see changes"
	  << " in visualization.";
  G4cout << G4endl;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructUniverse());
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::UpdateMaterials()
{
  DefineMaterials();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructUniverse());
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void ECRSDetectorConstruction::PrintAllParameters()
{
  PrintUniverseParameters();
  /*
  G4cout << G4endl;
  G4cout << "Atmospheric Layer:\tHeight Above Surface:\tHeight From Center:\n";
  G4cout << "Atm1\t\t\t" << G4BestUnit(mult*atmHeight[1],"Length") << "\t\t\t";
  G4cout << G4BestUnit(mult*(atmHeight[1]+earthRadius),"Length") << G4endl;
  */
}



//
// Here is the ATM model for determining the air density, provided by Kanishka.
// This model is modified to calculate density, temperature and pressure for input height, modified by Xiaohang
//

G4double ECRSDetectorConstruction::GetAirDensity(G4double hin, G4int option)
{
  G4double temp, press, density, h;
  G4int opt = option; //chose output parameters (temperature 2/pressure 1/density 0)
  h = 0.001*hin;  // convert to meters
  if(atmType==0){ // general atm model
    //fileOut->fout  <<" General atm model is selected " << G4endl;
    if(h>25000){   
      temp = -131.21+(0.00299*h);
      press = 2.488*( pow( (temp+273.1)/216.6 , -11.388) );
      density = (press/( 0.2869*(temp+273.1) ));
      //G4cout << " In function: density = " << density << G4endl;
      //return density;
    } 
    else if(h>11000 && h<=25000){   
      //temp = -56.46;
      //lower-stratosphere temperature change
      //if(h>=14000 && h<=18000) temp = -56.46+(216.64*(0+0.00));
      //else temp = -56.46;
      //temp = -56.46+(216.64*(0.));
      temp = -56.46;
      press = 22.65*(exp(1.73-(0.000157*h)));
      density = (press/( 0.2869*(temp+273.1) ));
      //G4cout << " In function: density = " << density << G4endl;
      //return density;
    } 
    else if (h<=11000){   
      temp = 15.04-(0.00649*h);
      press = 101.29*( pow( (temp+273.1)/288.08 , 5.256) );
      //density = (press/( 0.2869*(temp+273.1) ));
      density = (press/( 0.2869*(temp+273.1) )) * denPer;
      //G4cout << " In function: density = " << density << G4endl;
      //return density;
    } 
    else {
            // 9/28/2026: replaced exit(1) with G4Exception for clean termination.
      G4ExceptionDescription ed;
      ed << "Atmospheric layer height " << h << " is outside the valid range.";
      G4Exception("ECRSDetectorConstruction::atmModel", "ECRS_ATM001",
                  FatalException, ed);
      return 0.;
    }

    if(opt==0) return density;
    else if(opt==1) return press;
    else if(opt==2) return temp;
    else {
            // 9/28/2026: replaced exit(1) with G4Exception for clean termination.
      G4ExceptionDescription ed;
      ed << "Invalid option " << opt << " requested from the atmosphere model.";
      G4Exception("ECRSDetectorConstruction::atmModel", "ECRS_ATM002",
                  FatalException, ed);
      return 0.;
    }
    
  }

  if(atmType==1){ //summer time atm model
    // fileOut->fout <<"Summer time atm model is selected " << G4endl;
       if(h>25000){   
      temp = -97.835+(0.00182493*h);
      press = 2.86767*( pow( (temp+273.1)/218.146, -15.2137) );
      density = press/( 0.2869*(temp+273.1) );
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else if(h>11000 && h<=25000){   
      if(h>11000 && h<=16250) temp = 14.3842 - (0.00546857*h);
      if(h>16250 && h<=25000) temp = -112.07 + (0.00246027*h);
      press = 19.8803*(exp(2.03009-(0.000163279*h)));
      density = press/( 0.2869*(temp+273.1) );
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else if (h<=11000){   
      temp = 27.8357-(0.00608262*h);
      press = 100.094*( pow( (temp+273.1)/299.941, 5.66264) );
      density = press/( 0.2869*(temp+273.1) ); 
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else {
            // 9/28/2026: replaced exit(1) with G4Exception for clean termination.
      G4ExceptionDescription ed;
      ed << "Atmospheric layer height " << h << " is outside the valid range.";
      G4Exception("ECRSDetectorConstruction::atmModel", "ECRS_ATM001",
                  FatalException, ed);
      return 0.;
    }
  }
  
  
  if(atmType==2){ //winter time atm model
    //fileOut->fout  <<"Winter time atm model is selected " << G4endl;
    if(h>25000){   
      temp = -108.101 +(0.00205126*h);
      press = 2.88776*( pow( (temp+273.1)/206.093, -8.62931) );
      density = press/( 0.2869*(temp+273.1) );
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else if(h>11000 && h<=25000){   
      if(h>11000 && h<=18250) temp = -24.6769 - (0.00257295*h);
      if(h>18250 && h<=25000) temp = -103.906 + (0.00191275*h);
      press = 19.3072*(exp(1.97427-(0.000162203*h)));
      density = press/( 0.2869*(temp+273.1) );
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else if (h<=11000){   
      temp = 14.2545-(0.00583544*h);
      press = 109.005*( pow( (temp+273.1)/293.321, 5.29241) );
      density = press/( 0.2869*(temp+273.1) );
      G4cout << " In function: density = " << density << G4endl;
      return density;
    } else {
            // 9/28/2026: replaced exit(1) with G4Exception for clean termination.
      G4ExceptionDescription ed;
      ed << "Atmospheric layer height " << h << " is outside the valid range.";
      G4Exception("ECRSDetectorConstruction::atmModel", "ECRS_ATM001",
                  FatalException, ed);
      return 0.;
    }
  }
}
