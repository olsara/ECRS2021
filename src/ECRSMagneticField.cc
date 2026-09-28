// Updated on 5/29/2014, hexc, olesya
// Updated on 8/29/2014, hexc, olesya
//   Followed the MAGCOSMagneticField.cc code to implement the update IGRF table
//   We still need to understand why we have to implement a break statement in the igrf read 
// Updated on 7/23/2015, hexc, olesya
//   Corrected the year-range of field coefficients to the proper year-max.
#include "ECRSMagneticField.hh"
#include "ECRSFieldMessenger.hh"
#include "globals.hh"
#include "geomdefs.hh"
#include "magneto_fsubroutine_def.hh"
#include "time.h"
#include "G4ios.hh"
#include "fstream"
#include "G4MagIntegratorStepper.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "SpaceCoordinateConvertor.hh"
#include "ECRSEquationOfMotion.hh"
#include "ECRSUnits.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ProcessTable.hh"
#include "BSIntegrator.hh"
#include "BSEquationInG4.hh"

#include <iostream>
#include <sstream>   // 9/28/2026: replaced deprecated <strstream> (removed in C++20)

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////////////
//
ECRSMagneticField::ECRSMagneticField()
{ 
  //messenger
  theFieldMessenger = new ECRSFieldMessenger(this); 
  
  //Parameters
  
  Pdyn=2.;
  Dst=0.;
  Imfy=0.;
  Imfz=1.*nanotesla;
  G1=1.;
  G2=0.;
  nm_igrf=13;
  
  times_of_data=0;
  Pdyn_data =0;
  Dst_data=0;
  Imf_gsm_y_data=0;
  Imf_gsm_z_data=0;
  g1_tsy_param_data=0;
  g2_tsy_param_data=0;
  n_tsy_data=0;
  
  //Igrf table
  if (getenv("IGRF_TABLE")){
    ReadIgrfTable(getenv("IGRF_TABLE"));
    G4cout << " Done with reading the IGRF table! " << G4endl;
  }
  else {
    G4cout<<"The environment variable IGRF_TABLE has not been defined !!!"<<std::endl;
  }
  
  //Set pointers on IGRF table for the SpaceCoordinateConvertor
  SpaceCoordinateConvertor* theConvertor = SpaceCoordinateConvertor::getInstance();
  theConvertor->Setp_hh_nm(&hh_nm);
  theConvertor->Setp_gg_nm(&gg_nm);
  
  // Set the start date and compute IGRF coefficient 
  TimeOfB=0.; 
  StartDate = DateAndTime(2000,1,1);
  G4cout << "The start date for computing IGRF coefficient: 1/1/2000 (as a default)" << G4endl;
  ReferenceDate = StartDate;
  ComputeIgrfCoefAccordingToTime();
  theConvertor->SetReferenceDate(StartDate);
 
  //dipole parameter
  SetEccentricDipoleParameterFromIGRF();
  GetInternalField = &ECRSMagneticField::GetIGRF;
  GetExternalField = &ECRSMagneticField::GetTSY89;
  SelectedOutsideMagnetosphere = &ECRSMagneticField::IGRFOutsideMagnetosphere;
 
  External=false;
  Internal=true;
  iopt=1; 
  
  //EquationOfMotion
  fEquationOfMotion = new ECRSEquationOfMotion(this);
   
  //Integration methods

  fChordFinder = 0;
  fStepper = 0;
  DefaultDeltaChord = .01*Re;
  DeltaChord = DefaultDeltaChord;
  DefaultDeltaIntersection = .001*re;
  DefaultEpsilon = 1.e-6;
  DefaultBSMaxStep =.5*re;
    
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()
	                        ->SetLargestAcceptableStep(10.*Re);
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                              ->SetDetectorField(this);
    
  ResetIntegrationParameters(); 
    
  ConsiderDipoleShift=false;

}
////////////////////////////////////////////////////////////////////////
//
ECRSMagneticField::~ECRSMagneticField()
{ if (fEquationOfMotion) delete fEquationOfMotion;
  if (fChordFinder) delete fChordFinder;
  if (fStepper) delete fStepper;
}

///////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::GetFieldValue( const G4double yTrack[7],
				       G4double B[3]         ) const 
{ 
  G4ThreeVector pos=G4ThreeVector(yTrack[0],yTrack[1],yTrack[2]);
  
  G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);
  if (Internal) Bfield = (this->*GetInternalField)(pos);
  if (External) Bfield=Bfield+ (this->*GetExternalField)(pos);
  
  B[0]=Bfield.x();
  B[1]=Bfield.y();
  B[2]=Bfield.z();    
  
  return ;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector ECRSMagneticField::GetFieldValue
(const G4ThreeVector GEOposition) const 
{ 
  G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);
  if (Internal) Bfield = (this->*GetInternalField)(GEOposition);
  if (External) Bfield=Bfield+ (this->*GetExternalField)(GEOposition);
  return Bfield;
}			       
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector ECRSMagneticField::FindPositionOnMagnetopause(G4double theta, 
							    G4double xgsm, 
							    G4double precision) const
{ 
  G4ThreeVector pos1 = G4ThreeVector(xgsm,0.,0.);
  
  SpaceCoordinateConvertor* theCoordinateConvertor
    = SpaceCoordinateConvertor::getInstance();
  theCoordinateConvertor->SetSystemInAndOut("GSM","GEO");
  
  if  (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos1)))
    return G4ThreeVector (0.,0.,0.);
  
  G4ThreeVector pos2 =G4ThreeVector(xgsm,40.*re*std::sin(theta),40.*re*std::cos(theta));
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
    pos3= (pos1+pos2)/2.;
    if (OutsideMagnetosphere(theCoordinateConvertor->Transform(pos3))) pos2=pos3;
    else pos1=pos3;
  }
  return (pos1+pos2)/2.;
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector ECRSMagneticField::FindStandOffPosition(G4double precision ) const
{ 
  G4ThreeVector pos1 = G4ThreeVector(0.,0.,0.)*Re;
  G4ThreeVector pos2 = G4ThreeVector(40.,0.,0.)*Re;
  G4ThreeVector pos3;  
  while ( (pos2-pos1).mag()> precision){
    pos3= (pos1+pos2)/2.;
    if (OutsideMagnetosphere(pos3)) pos2=pos3;
    else pos1=pos3;
  }
  return pos1;	
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetTimeOfB(G4double val)
{ 
  TimeOfB=val;
  G4double JulDate = StartDate.JulianDate() + ( val / (24. * 3600.) );
  /* G4cout<<StartDate.min<<std::endl;
     G4cout<<StartDate.sec<<std::endl;*/
  ReferenceDate = DateAndTime(JulDate);
  // G4cout<<ReferenceDate.month<<" Reference date"<<std::endl;
  /*G4cout<<ReferenceDate.min<<std::endl;
    G4cout<<ReferenceDate.sec<<std::endl;*/
  ComputeTSYParameter(val);
  
  // new Igrf Coeff 
  ComputeIgrfCoefAccordingToTime();
  
  // set Reference date for Space coordinate convertor 
  SpaceCoordinateConvertor::getInstance()->SetReferenceDate(ReferenceDate);
  
  // new dipole parameters 
  SetEccentricDipoleParameterFromIGRF();
}

////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetStartDate(G4int Year, G4int Month,G4int Day,
				     G4int Hour, G4int Minute,G4int Second )
{
  // Set the start date
  StartDate = DateAndTime(Year,Month,Day,Hour,Minute,Second);
  G4cout << "Set the field start date to (Y, M, D, H, M, S): " << Year << "  " << Month;
  G4cout << "  " << Day << "  " << Hour << "  " << Minute << "  " << Second << G4endl;
  double val =TimeOfB;
  SetTimeOfB(val);
}		

////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetIopt(  G4int val)
{
  iopt=val;
  if (val <1) iopt=1;
  if (val >7) iopt=7;
  return ;
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetInternalField(  G4String aString)
{ 
  Internal=true;
  if (aString == "DIPOLE") 
    GetInternalField=&ECRSMagneticField::GetGEODipole;
  else if (aString == "IGRF")
    GetInternalField=&ECRSMagneticField::GetIGRF; 
  else if (aString == "NOFIELD")
    Internal=false;     
}
///////////////////////////////////////////////////////////////////////////////

void ECRSMagneticField::SetExternalField(  G4String aString)
{External=false;
  if (aString == "TSY89"){
    External=true;
    GetExternalField=&ECRSMagneticField::GetTSY89;
    SetMagnetopauseModel("TSY89");
  }
  else if (aString == "TSYBOB89"){ 
    External=true;
    GetExternalField=&ECRSMagneticField::GetTSYBOB89;
    SetMagnetopauseModel("TSY89");
  } 
  else if (aString == "TSY96"){
    External=true;
    GetExternalField=&ECRSMagneticField::GetTSY96;
    SetMagnetopauseModel("TSY96");
  }
  else if (aString == "TSY2001"){
    External=true;
    GetExternalField=&ECRSMagneticField::GetTSY2001;
    SetMagnetopauseModel("TSY2001");
  }
  else  SetMagnetopauseModel("NOFIELD");
  
}
////////////////////////////////////////////////////////////////////////////////

void  ECRSMagneticField::SetMagnetopauseModel( G4String aString) 
{// G4cout<<aString<<std::endl; 
  if (aString == "TSY2001")
    SelectedOutsideMagnetosphere
      =&ECRSMagneticField::TSY2001OutsideMagnetosphere;
  
  else if (aString == "TSY96")
    SelectedOutsideMagnetosphere
      =&ECRSMagneticField::TSY96OutsideMagnetosphere; 		  
  
  else  if (aString == "TSY89")
    SelectedOutsideMagnetosphere
      =&ECRSMagneticField::TSY89OutsideMagnetosphere; 
  
  else 	
    SelectedOutsideMagnetosphere
                  =&ECRSMagneticField::IGRFOutsideMagnetosphere; 	    
  
}
////////////////////////////////////////////////////////////////////////////////

void ECRSMagneticField::SetDipoleB0(G4double aVal)
{ 
  DipoleB0 = aVal;
  bdip_.b0dip = aVal/nanotesla;
}

////////////////////////////////////////////////////////////////////////////////

void ECRSMagneticField::SetTiltedDipoleParameterFromIGRF()
{
  SpaceCoordinateConvertor* theConvertor = 
    SpaceCoordinateConvertor::getInstance();
  
  DipoleB0=theConvertor->GetGeoDipoleB0();
  G4ThreeVector DipoleAxis=theConvertor->GetGeoDipoleAxisInGEO();
  DipoleTheta=DipoleAxis.theta();
  DipolePhi=DipoleAxis.phi();
  bdip_.b0dip=DipoleB0/nanotesla;
  DipolePS=theConvertor->GetTiltAngle();
  DipoleSchift=G4ThreeVector(0.,0.,0.);
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetEccentricDipoleParameterFromIGRF()
{ 
  SpaceCoordinateConvertor* theConvertor = 
    SpaceCoordinateConvertor::getInstance();
  DipoleB0=theConvertor->GetGeoDipoleB0();
  G4ThreeVector DipoleAxis=theConvertor->GetGeoDipoleAxisInGEO();
  DipoleTheta=DipoleAxis.theta();
  DipolePhi=DipoleAxis.phi();
  bdip_.b0dip=DipoleB0/nanotesla;
  DipolePS=theConvertor->GetTiltAngle();
  DipoleSchift=theConvertor->GetGeoDipoleSchiftInGEO();
 
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetEpsilon(G4double aVal)
{ //precision for G4Integration methods
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMinimumEpsilonStep(aVal);
  G4TransportationManager::GetTransportationManager()
                                  ->GetPropagatorInField()
                                  ->SetMaximumEpsilonStep(aVal*1.00001);
  //precision for Bulirsh Stoer method
  BSIntegrator::getInstance()->SetEps(aVal);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetDeltaChord(G4double aVal)
{ if (fChordFinder) fChordFinder->SetDeltaChord(aVal);
  DeltaChord=aVal;
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetBSMaxStep(G4double aVal)
{ BSIntegrator::getInstance()->SetHmax(aVal);
}

////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetDeltaIntersection(G4double aVal)
{ ///////Delta intersection for G4Integration method
  G4TransportationManager::GetTransportationManager()->GetFieldManager()
                                                 ->SetDeltaIntersection(aVal);

 ///////Delta intersection Bulirsh method
 BSEquationInG4* theBSEquation = 
           ((BSEquationInG4*) BSIntegrator::getInstance()
                                ->GetEquationTobeIntegrated());
 if (theBSEquation) theBSEquation->SetCrossingDelta(aVal);

}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::ResetIntegrationParameters()
{ SetStepper("CashKarpRKF45");
  SetDeltaChord(DefaultDeltaChord);
  SetDeltaIntersection(DefaultDeltaIntersection);
  SetEpsilon(DefaultEpsilon);
  SetBSMaxStep(DefaultBSMaxStep);
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SelectBulirshStoerMethod()
{ 
  // Commented out by hexc  5/13/2009
  // G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",true);
  // G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",false);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SelectG4IntegrationMethods()
{ G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",false);
  G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",true);
 
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::ReadTSYParameter(G4String nameFile)
{
  if (times_of_data) delete times_of_data;
  if (Pdyn_data) delete Pdyn_data;
  if ( Dst_data) delete  Dst_data;
  if ( Imf_gsm_y_data) delete  Imf_gsm_y_data;
  if ( Imf_gsm_z_data) delete  Imf_gsm_z_data;
  if ( g1_tsy_param_data) delete  g1_tsy_param_data;
  if ( g2_tsy_param_data) delete  g2_tsy_param_data;
 
  n_tsy_data=0;
  std::fstream File_Input(nameFile, (std::ios::binary | std::ios::in));
  File_Input>>n_tsy_data;
 
  if (n_tsy_data <= 0){
  	G4cout<<"the file does not exist or is not well written"<<G4endl;
  } 
  else{ 
  	times_of_data= new G4double[n_tsy_data];
   	Pdyn_data= new G4double[n_tsy_data];
   	Dst_data= new G4double[n_tsy_data];
   	Imf_gsm_y_data= new G4double[n_tsy_data];
   	Imf_gsm_z_data= new G4double[n_tsy_data];
   	g1_tsy_param_data= new G4double[n_tsy_data];
   	g2_tsy_param_data= new G4double[n_tsy_data];
 
   	G4int nyear;
   	G4int nmonth;
   	G4int nday;
   	G4int nhour;
   	G4int nminute;
   	G4int nsecond=0;
   	File_Input>>nyear>>nmonth>>nday>>nhour>>nminute;
   	SetStartDate(nyear, nmonth, nday, nhour, nminute, nsecond); 
 
   	for (G4int i=0;i<n_tsy_data;i++){
   		File_Input>>times_of_data[i]
	      			>>Pdyn_data [i]
              			>> Dst_data[i]
              			>> Imf_gsm_y_data[i]
              			>> Imf_gsm_z_data[i]
              			>> g1_tsy_param_data[i]
              			>>g2_tsy_param_data[i];
   	}
  } 
  File_Input.close(); 
  double val =TimeOfB;
  SetTimeOfB(val);
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::ReadIgrfTable(G4String nameFile)
{ 
  n_tsy_data=0;
  std::fstream File_Input(nameFile,  std::ios::in);
 
  //clear vectors
  igrf_year.clear();
  h_nm.clear();
  g_nm.clear();
  dh_nm.clear();
  dg_nm.clear();
  hh_nm.clear();
  gg_nm.clear();
  recurrence_coef.clear();

  //Process the first three lines 
  char Line[300];
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
  File_Input.getline(Line,300);
  G4cout << Line << G4endl;
 
  char ch;
  std::stringstream astream;   // 9/28/2026: was std::strstream (deprecated)
  File_Input.get(ch);
  while (ch != '\n') {
    astream<<ch;
    File_Input.get(ch);
  }

  G4String str1;
  astream>>str1>>str1>>str1;
  if (str1!="m"){
	// 9/28/2026: use G4Exception instead of exit(0) so the run terminates
	// cleanly (destructors, MT worker synchronisation) with a clear message.
	G4ExceptionDescription ed;
	ed << "IGRF_TABLE file '" << nameFile << "' is not formatted as expected "
	   << "(expected an 'm' column header, found '" << str1 << "').";
	G4Exception("ECRSMagneticField::ReadIgrfTable", "ECRS_IGRF001",
		    FatalException, ed);
  }
  G4double year=0.;
  while (!astream.eof() && year >= 0.){
    year=-999999.;
    astream>>year;
    if (year!=-999999.){ 
      igrf_year.push_back(year);
    }
  }
  //number of column
  
  G4int nyear = igrf_year.size();
  G4cout << "Number of years read in: " << nyear << G4endl;
  nmax=0;
  G4int n,m;
  
  // initialisation h_nm,g_nm;  
  
  for (int i=0; i <nyear ; i++){
    h_nm.push_back(std::vector<float> ());
    g_nm.push_back(std::vector<float> ());
  }
  
  G4int nStop = 0;

  while (!File_Input.eof()){
    if (nStop++ > 200) break;  // temp stuff
    str1="bibab";
    n=-999;
    m=-999;
    File_Input>>str1>>n>>m;
    if (str1 != "bibab" && n != -999 && m != -999) {
      // G4cout<<str1<<'\t'<<n<<'\t'<<m<<'\t';
      if (n>nmax){
	nmax=n;
	for (int i=0; i <nyear ; i++){
	  h_nm[i].insert(h_nm[i].end(),nmax+1,0.);
	  g_nm[i].insert(g_nm[i].end(),nmax+1,0.);
	}
	dg_nm.insert(dg_nm.end(),nmax+1,0.);
	dh_nm.insert(dh_nm.end(),nmax+1,0.);
	hh_nm.insert(hh_nm.end(),nmax+1,0.);
	gg_nm.insert(gg_nm.end(),nmax+1,0.);	    
      }
      G4int index = ((n+2)  * (n-1)) /2 + m;
      //G4cout<<"index "<<index<<std::endl;
      G4double value;
      int i_start=0;
      
      // The following section was for reading old IGRF data table
      /*
      if (n >10) {    // For old igrf table, no data after n>10
	i_start = 20;
	for (int i=0;i<i_start;i++){	
	  h_nm[i][index]  = 0.;
	  g_nm[i][index]  = 0.;
	}
      }				
      */

      for (int i= i_start; i <nyear ; i++){
	File_Input>>value;
	//G4cout<<value<<'\t';
	if ( (str1 == "h") && (i < nyear-1) ) h_nm[i][index] = value;
	if ( (str1 == "g") && (i < nyear-1) ) g_nm[i][index] = value;
	if ( (str1 == "h") && (i == nyear-1) ) dh_nm[index] = value;
	if ( (str1 == "g") && (i == nyear-1) ) dg_nm[index] = value;     	
      }
      
    }
  }

  G4cout << " Table reading is finished!!!!! " << G4endl;
  
  //Here we calculate the coefficient for the recurrence equation used to compute 
  //the associated Legendre function in the Gauss normalsation. 
  // In this normalisation the equation of recurrence is
  // P_n_m(x) = x * p_n-1_m -(n+m)*(n-m)* p_n-2_m / ((2n-3)*(2n-1)) (1)
  // The coefficient that we calculate here is (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1))
  // The qeaution (1) is a transformation of the classical  equation of recurrence of lengendre assocaited function
  // (n-m) * P_n_m = x  * p_n-1_m *(2n-1) -(n + m -1) * p_m_n-2  (2)
  // by using the Gauss normalisation (2) transform in (1)
  // (1) is better for computing purpose because you have just to compute 
  //the coefficient  (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) at initialisation an multiplied it by
  //   p_m_n-2  when computinhg p_n_m
  // In  case of (2) you have to compute three vector of cooefficient at initilaisation you have to
  // make three multiplications instead of one 
  

  for (n=1;n<nmax + 1;n++)
    for (int m= 0; m < n+1; m++) 
      recurrence_coef.push_back(double ((n-1+m)*(n-1-m)) / double ((2*n-1)*(2*n-3)) );
  File_Input.close();

}

////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::ComputeIgrfCoefAccordingToTime()
{
  G4int nyear = igrf_year.size();
  G4int n=(nmax+3)*nmax/2;
  
 // G4cout<<"number of year is "<<nyear<<endl;
 
  if (ReferenceDate.year < igrf_year[0]  || 
                     ReferenceDate.year > igrf_year[nyear-1]){
     	G4cout<<"The year should be in the period "<<int(igrf_year[0])
                                   <<"-"<<int(igrf_year[nyear-1]) <<std::endl;
      
      	G4cout<<"The coefficient will be set according to year "<<
       	                                   int(igrf_year[nyear-1])<<std::endl;
      	for (int i=0;i<n;i++){
       		gg_nm[i]=g_nm[nyear-1][i];
        	hh_nm[i]=h_nm[nyear-1][i];}
  }
  else {
     	G4int index=0;
     
      	while (ReferenceDate.year >= (igrf_year)[index] && index < nyear) index++;
      
     	// G4cout <<"Selected index"<<index<<std::endl;
      	if (index<nyear){
        	DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
	 	DateAndTime Date2=DateAndTime(int (igrf_year[index]), 1,1);
	 	G4double delta_day12=Date2.DifferenceInDays(Date1);
	 	G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
	 	G4double f2= delta_day / delta_day12;
         	G4double f1=1.-f2;
       		for (int i=0;i<n;i++){
          		gg_nm[i]=g_nm[index-1][i]*f1 + g_nm[index][i]*f2;
           		hh_nm[i]=h_nm[index-1][i]*f1 + h_nm[index][i]*f2;
		}
      	}
      	else {
        	DateAndTime Date1=DateAndTime(int (igrf_year[index-1]), 1,1);
	 	G4double delta_day=ReferenceDate.DifferenceInDays(Date1);
	 	G4double dt=delta_day / 365.25;
         	for (int i=0;i<n;i++){
          		gg_nm[i]=g_nm[nyear-1][i] + dt*dg_nm[i];
           		hh_nm[i]=h_nm[nyear-1][i] + dt*dh_nm[i];
		}
        }
  }
     
  //The coefficient are  Schmidt quasi normalised
  // For using the recuurence relation in Gauss normalisation
  // we have to convert the coefficient in Gauss normalisation
  // the g_nm and h_nm coefficient become
  G4double s=1.;
  for (n=1;n<nmax+1;n++){
  	G4int mn = ((n + 1) * n) / 2  -1;
       	s*=double(2*n-1)/double(n);
       	hh_nm[mn] *= s;
       	gg_nm[mn] *= s;
       	G4double p=s;
       	for (int m=1;m<n+1;m++){
        	G4double aa=1.;
	 	if (m == 1) aa=2.;
	 	p *= std::sqrt(aa * double(n - m + 1) / double (n + m) );
	 	int nmn =  mn + m ;
	 	hh_nm[nmn] *= p;
         	gg_nm[nmn] *= p;
	}
  }   				   			     
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::PrintBfield(G4ThreeVector geo_pos) const
{ SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 
		   
  G4ThreeVector Bfield =G4ThreeVector(0.,0.,0.);
  if (Internal) Bfield = (this->*GetInternalField)(geo_pos);
  if (External) Bfield=Bfield+ (this->*GetExternalField)(geo_pos);
  Bfield /= nanotesla;
  G4ThreeVector Bgsm = theConvertor->Transform(Bfield, "GEO", "GSM");
  G4ThreeVector Bmag = theConvertor->Transform(Bfield, "GEO", "MAG");
  G4ThreeVector Bsm = theConvertor->Transform(Bfield, "GEO", "SM");
 
  G4cout<<"Bfield at primary position" <<std::endl; 		   
  G4cout<<"Bfield GEO: "<<Bfield<<" "<<Bfield.mag()<<std::endl;
  G4double altitude, longitude, latitude;
  theConvertor->ComputeGEOIDCoordinatesFromGEOPosition(geo_pos,
 						     altitude,
						     longitude,latitude);
  //Northward, eastward, dowmward component
  //---------------------------------------
 
  G4ThreeVector vertical = G4ThreeVector(1.,0.,0.);
  vertical.setRThetaPhi(1.,90.*degree-latitude,longitude);
  G4ThreeVector east_dir =G4ThreeVector(-std::sin(longitude),std::cos(longitude),0. );
  G4ThreeVector north_dir = -east_dir.cross(vertical);
  G4double Bdownward = -Bfield.dot(vertical);
  G4double Bnorthward = Bfield.dot(north_dir);
  G4double Beastward = Bfield.dot(east_dir);
 
  G4cout<<"Bfield north, east, down: "<<Bnorthward
                       <<" "<<Beastward<<" "<<Bdownward<<std::endl;
  G4cout<<"Bfield GSM: "<<Bgsm<<" "<<Bgsm.mag()<<std::endl;
  G4cout<<"Bfield SM: "<<Bsm<<" "<<Bsm.mag()<<std::endl;
  G4cout<<"Bfield MAG: "<<Bmag<<" "<<Bmag.mag()<<std::endl; 
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::PrintStormParameter()
{ G4cout<<"iopt "<<iopt<<G4endl;
  G4cout<<"Pdyn "<<Pdyn<<G4endl;
  G4cout<<"Dst "<<Dst/nanotesla<<G4endl;
  G4cout<<"Imfy "<<Imfy/nanotesla<<G4endl;
  G4cout<<"Imfz "<<Imfz/nanotesla<<G4endl;
  G4cout<<"G1 "<<G1<<G4endl;
  G4cout<<"G2 "<<G2<<G4endl;
  G4cout<<"PS "<<DipolePS/degree<<G4endl;
}
////////////////////////////////////////////////////////////////////////////////

G4ThreeVector ECRSMagneticField::GetIGRFFortran(G4ThreeVector pos) const 
{ G4double r=pos.mag()/re;
  G4double t=pos.theta();
  G4double p=pos.phi();
  float Br=1.;
  float Bp=1.;
  float Bt=1.;
  G4int nm=nm_igrf;
  //year_igrf=1980;
   
  G4double ndays_per_y  = DateAndTime(ReferenceDate.year+1,1,1)
                                  .DifferenceInDays(DateAndTime(ReferenceDate.year,1,1));
   
  G4double fraction_of_year =
                     ReferenceDate.DifferenceInDays(
		                    DateAndTime(ReferenceDate.year,1,1))/ndays_per_y;
  float yea=float(ReferenceDate.year + fraction_of_year)  ;
  float r1 =float(r);
  float t1 =float(t);
  float p1 =float(p);

   
  if (nm > 10) nm=10;
  
  // Modified by hexc 5/19/2009 in order to match the fort fct name
  igrf_to_cc__(&yea,&nm,&r1,&t1,&p1,&Br,&Bt,&Bp);
   
  G4double Bunit=nanotesla;
  G4double st=sin(t);
  G4double ct=cos(t);
  G4double sp=sin(p);
  G4double cp=cos(p);
  G4double Be=double(Br)*st+double(Bt)*ct;
  G4double Bx=Bunit*(Be*cp-double(Bp)*sp);
  G4double By=Bunit*(Be*sp+double(Bp)*cp);
  G4double Bz=Bunit*(double(Br)*ct-double(Bt)*st);
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////

G4ThreeVector ECRSMagneticField::GetIGRF(G4ThreeVector pos) const 
{ 
  // This method computes the igrf magnetic field
  // For this purpose associated legendre function are calculated by the recurrence equation
  // P_n_m(x) = x * p_n-1_m(x) - recurrence_coef(n,m) * p_n-2_m(x)  (1)
  // where recurrence_coef = (n+m-1)*(n-m-1)/ ((2n-3)*(2n-1)) and x= cos_theta
  // By deriving  equation (1) by theta we obtain
  // dp_n_m/dtheta = -sin_theta * p_n-1_m + cos_theta * dp_n-1_m/dtheta  -  
  //                                                        recurrence_coef * dp_n-2_m/dtheta
  // and dp_n_m/dtheta can  also be computed by recurrence. 
  // we should also now p_m_m = (sin (theta))^m


  float r=pos.mag()/re;
  float theta=pos.theta();
  if (theta < 0.01*degree) theta =0.01*degree;
  else if (theta > 179.99*degree) theta =179.99*degree;
  float phi=pos.phi();
   
  float pp = 1. / r;
  float p =  pp * pp;
 
   
   
   
  float cos_phi = std::cos(phi);
  float two_cos_phi = cos_phi * 2.;   // 9/28/2026: reuse cos_phi (was recomputing std::cos(phi))
  float sin_phi = std::sin(phi);
  float cos_theta = std::cos(theta);
  float sin_theta = std::sin(theta);
   
   
   
  // a_n = (1 /r) ^ (n+1)
  // da_n_dth = -(a_n1) * n+1 
 
  // 9/28/2026: Use fixed-size stack arrays instead of std::vector to avoid a
  // heap allocation (and vector growth) on every field evaluation - GetIGRF is
  // called from GetFieldValue on the integration hot path. nm_igrf is <= 13.
  const G4int kMaxHarm = 16;
  float a[kMaxHarm];
  float da_dr[kMaxHarm];
  G4int n;
  G4int old_nm_igrf =nm_igrf;
  if (ReferenceDate.year <2000. && nm_igrf >10) old_nm_igrf =10; 	
  if (old_nm_igrf > kMaxHarm) old_nm_igrf = kMaxHarm;
   
   
  for ( n=1;n<old_nm_igrf+1;n++){
    	p *= pp;
     	a[n-1] = p;
     	da_dr[n-1] = -p * float (n+1);
  }
    
  // case where m=0, p_1_0 =1.
  //---------------------------
    
  // n=1 
  float p_nm=cos_theta;
  float dp_nm_dth=-sin_theta;
    
  float Br=-da_dr[0]*gg_nm[0]*p_nm;
  float Bt=-a[0]*gg_nm[0]*dp_nm_dth;
  float Bp=0.;
    
  float p_nm1=p_nm; //corresponds to P for n-1 and m
  float p_nm2=1.;   //corresponds to P for n-2 and m
  
  float dp_nm1_dth=dp_nm_dth;
  float dp_nm2_dth=0.;
    
  // n>2
    
 for (n = 2 ; n < old_nm_igrf +1 ; n++ ){
    	G4int ii= ((n+2) * (n-1) )/2; 
      	p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
      	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                       - recurrence_coef[ii]*dp_nm2_dth;
      	Br-=da_dr[n-1]*gg_nm[ii]*p_nm;
      	Bt-=a[n-1]*gg_nm[ii]*dp_nm_dth; 
      	p_nm2=p_nm1;
      	p_nm1=p_nm;
      	dp_nm2_dth=dp_nm1_dth;
      	dp_nm1_dth=dp_nm_dth;
  } 
    
  // continue for m>0
  //---------  
    
  float cos_mphi2=cos_phi;
  float cos_mphi1=1.;
  float cos_mphi=cos_phi;
  float sin_mphi2=-sin_phi;
  float sin_mphi1=0.;
  float sin_mphi=sin_phi;
    
    
  float p_mm=1.;
  float dp_mm_dth=1.; 
  float Bpm=0.;
  float Bp_hh,Bp_gg,Br_hh,Br_gg,Bt_hh,Bt_gg;
  for (int m = 1 ; m < old_nm_igrf +1 ; m++ ){
    	//Initialisation
     	//------------------
     
       	Bp_hh=0.;
       	Bp_gg=0.;
       	Br_hh=0.;
       	Br_gg=0.;
       	Bt_hh=0.;
       	Bt_gg=0.;
      
     	//compute cos_mphi and sin_mphi with
     	// recurrence formula from numerical recipies in C++ p 184 
      	cos_mphi= two_cos_phi*cos_mphi1 -cos_mphi2;
      	sin_mphi= two_cos_phi*sin_mphi1 -sin_mphi2;
      
      
      	cos_mphi2=cos_mphi1;
      	sin_mphi2=sin_mphi1;
      	cos_mphi1=cos_mphi;
      	sin_mphi1=sin_mphi;
    
      	dp_mm_dth= float(m) * p_mm *cos_theta;
      	p_mm *= sin_theta;
      	Bpm=0;
     	// case where n=m
      	n=m;
      	G4int ii= ((n+3) * (n) )/2 -1; 
      	p_nm = p_mm;
      	dp_nm_dth = dp_mm_dth;
      
      	float ggg = gg_nm[ii]*p_nm;
      	float hhh = hh_nm[ii]*p_nm;
      
      	Br_gg-=da_dr[n-1]*ggg;
      	Br_hh-=da_dr[n-1]*hhh;
      
      	Bt_gg-=a[n-1]*gg_nm[ii]*dp_nm_dth;
      	Bt_hh-=a[n-1]*hh_nm[ii]*dp_nm_dth;
      
      	//Bt-=a[n-1]*gh*dp_nm_dth;
      
      	//Bpm-=a[n-1]*dgh_dphi*p_nm; 
      
      	Bp_gg-=a[n-1]*ggg;
      	Bp_hh-=a[n-1]*hhh; 
      
      
      	p_nm1=p_nm; //corresponds to P for n-1 and m
      	p_nm2=0.;  //corresponds to P for n-2 and m
      	dp_nm1_dth=dp_nm_dth;
     	 dp_nm2_dth=0.;
        
     
        for ( n = m + 1 ; n < old_nm_igrf +1 ; n++ ){
        	G4int ii= ((n+2) * (n-1) )/2 + m; 
         	p_nm = cos_theta *p_nm1 - recurrence_coef[ii]*p_nm2;
         	dp_nm_dth = -sin_theta * p_nm1 +  cos_theta *dp_nm1_dth 
                                                 - recurrence_coef[ii]*dp_nm2_dth;
	 
		 /*gh = gg_nm[ii]*cos_mphi + hh_nm[ii]*sin_mphi;
         	dgh_dphi = hh_nm[ii]*cos_mphi - gg_nm[ii]*sin_mphi;
        
	 	Br-=da_dr[n-1]*gh*p_nm;
         	Bt-=a[n-1]*gh*dp_nm_dth;
         	Bpm-=a[n-1]*dgh_dphi*p_nm;*/ 
	 
	 
	 	float ggg = gg_nm[ii]*p_nm;
         	float hhh = hh_nm[ii]*p_nm;
      
         	Br_gg-=da_dr[n-1]*ggg;
        	 Br_hh-=da_dr[n-1]*hhh;
      
         	Bt_gg-=a[n-1]*gg_nm[ii]*dp_nm_dth;
         	Bt_hh-=a[n-1]*hh_nm[ii]*dp_nm_dth;
      
         	Bp_gg-=a[n-1]*ggg;
         	Bp_hh-=a[n-1]*hhh; 

	 	p_nm2=p_nm1;
	 	p_nm1=p_nm;
	 	dp_nm2_dth=dp_nm1_dth;
	 	dp_nm1_dth=dp_nm_dth;
	}
	Br+=Br_gg*cos_mphi + Br_hh*sin_mphi;
	Bt+=Bt_gg*cos_mphi + Bt_hh*sin_mphi;
	Bp+=float(m)*(Bp_hh*cos_mphi - Bp_gg*sin_mphi);

  } 
  
  if (sin_theta >= 1.e-8) Bp/=sin_theta; 
  else if (cos_theta<0.) Bp=-Bp;
   
  G4double Be= Br*sin_theta+ Bt*cos_theta;
  G4double Bx= Be*cos_phi- Bp*sin_phi;
  G4double By= Be*sin_phi+ Bp*cos_phi;
  G4double Bz= Br*cos_theta - Bt*sin_theta;
  if (r < 0.01){   
   	Bx=.00000001*tesla;
    	By=0.;
    	Bz=0.;
  }
    
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz)*nanotesla;
  return Bfield;
}
/////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector ECRSMagneticField::GetGEODipole(G4ThreeVector pos) const
{ G4double xgeo,ygeo,zgeo; 
  xgeo=pos.x()/re;
  ygeo=pos.y()/re;
  zgeo=pos.z()/re;
   
  G4ThreeVector SchiftedPosition = 
  G4ThreeVector(xgeo,ygeo,zgeo)-(DipoleSchift/re); 
   
  G4double sth=sin(DipoleTheta);
  G4double cth=cos(DipoleTheta);
  G4double sphi=sin(DipolePhi);
  G4double cphi=cos(DipolePhi);
   
  //coordinates in DipoleSystem
  G4double xd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*cth
                -sth*SchiftedPosition.z();
  G4double yd=SchiftedPosition.y()*cphi-SchiftedPosition.x()*sphi;
  G4double zd=(SchiftedPosition.x()*cphi+SchiftedPosition.y()*sphi)*sth
                +cth*SchiftedPosition.z();
   
   
  //component in dipole system		
  G4double Bxd,Byd,Bzd;
  G4double r2=xd*xd +yd *yd +zd*zd; 
  G4double r5=std::pow(r2,2.5);
  G4double factor=DipoleB0/r5;
  Bxd=-3*factor*xd*zd;
  Byd=-3*factor*yd*zd;
  Bzd=factor*(r2-3*zd*zd);
   
   
  //back to Geo coordinate
  G4double Bx,By,Bz;
  Bx=(cth*cphi*Bxd - sphi*Byd  + sth*cphi*Bzd);
  By=(cth*sphi*Bxd + cphi*Byd  + sth*sphi*Bzd);
  Bz=(-sth*Bxd  + cth*Bzd);
   
  G4ThreeVector Bfield=G4ThreeVector(Bx,By,Bz);
  return Bfield;
}
////////////////////////////////////////////////////////////////////////////////
/////////////// External field models

G4ThreeVector ECRSMagneticField::GetTSY89( G4ThreeVector pos) const
{  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

   G4ThreeVector pos_geo=pos;
   if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
   //GEOtoGSM
   G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
   
   float xgsm=float(pos_gsm.x()/re);
   float ygsm=float(pos_gsm.y()/re);
   float zgsm=float(pos_gsm.z()/re);
  
 
   
   float parmod[10]; 
   float Bxgsm,Bygsm,Bzgsm;
   float ps=float(DipolePS);
   G4int ipar=iopt;
   // 9/11/2014, hexc:  added anothe "_" below
   t89c__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                      &Bxgsm, &Bygsm, &Bzgsm );
   
   G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
   Bgsm *= nanotesla;
   
   //GSMtoGEO
   return theConvertor->TransformGSMinGEO(Bgsm);
}
////////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector ECRSMagneticField::GetTSYBOB89( G4ThreeVector pos) const
{  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

   G4ThreeVector pos_geo=pos;
   if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
   //GEOtoGSM
   G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
   
   float xgsm=float(pos_gsm.x()/re);
   float ygsm=float(pos_gsm.y()/re);
   float zgsm=float(pos_gsm.z()/re);
  
 
   
   float parmod[10]; 
   parmod[1]=Dst;
   float Bxgsm,Bygsm,Bzgsm;
   float ps=float(DipolePS);
   G4int ipar=iopt;
   t89c_boberg__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                      &Bxgsm, &Bygsm, &Bzgsm );
   
   G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
   Bgsm *= nanotesla;
   
   //GSMtoGEO
   return theConvertor->TransformGSMinGEO(Bgsm);
}
////////////////////////////////////////////////////////////////////////////
G4ThreeVector ECRSMagneticField::GetTSY96(G4ThreeVector pos) const
{ 
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

  G4ThreeVector pos_geo=pos;
  if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  
  float xgsm=float(pos_gsm.x()/re);
  float ygsm=float(pos_gsm.y()/re);
  float zgsm=float(pos_gsm.z()/re);
  

  float parmod[10]; 
  parmod[0]=Pdyn;
  parmod[1]=Dst/nanotesla; 
  parmod[2]=Imfy/nanotesla;
  parmod[3]=Imfz/nanotesla;
   
   
  float Bxgsm,Bygsm,Bzgsm;
  float ps=float(DipolePS);
  G4int ipar=iopt;

  // changed by hexc:  from t96_01 to mcos_t96   5/19/2009
  mcos_t96_01__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                        &Bxgsm, &Bygsm, &Bzgsm );

  G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
  Bgsm *= nanotesla;
   
  //GSMtoGEO
  return theConvertor->TransformGSMinGEO(Bgsm);

}
////////////////////////////////////////////////////////////////////////////////

G4ThreeVector
   ECRSMagneticField::GetTSY2001( G4ThreeVector pos) const
{ 
 SpaceCoordinateConvertor* theConvertor =
                             SpaceCoordinateConvertor::getInstance(); 

   
 G4ThreeVector pos_geo=pos;
 
 // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
 //GEOtoGSM
 G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
 
 float xgsm=float(pos_gsm.x()/re);
 float ygsm=float(pos_gsm.y()/re);
 float zgsm=float(pos_gsm.z()/re);
  
// call tsy2001 fortran subroutine   
   
 float parmod[10]; 
 parmod[0]=Pdyn;
 parmod[1]=Dst/nanotesla; 
 parmod[2]=Imfy/nanotesla;
 parmod[3]=Imfz/nanotesla;
 parmod[4]=G1;
 parmod[5]=G2;
   
 float Bxgsm,Bygsm,Bzgsm;
 float ps=float(DipolePS);
   
 G4int ipar=iopt;
 t01_01__(&ipar,parmod,&ps,&xgsm,&ygsm,&zgsm,
                                        &Bxgsm, &Bygsm, &Bzgsm );
   
 G4ThreeVector Bgsm=
             G4ThreeVector(G4double(Bxgsm),G4double(Bygsm),G4double(Bzgsm));
 Bgsm *= nanotesla;
   
 //GSMtoGEO
 return theConvertor->TransformGSMinGEO(Bgsm);
}

////////////////////////////////////////////////////////////////////////////////
//-------------------------------Outside magnetosphere
bool ECRSMagneticField::OutsideMagnetosphere( G4ThreeVector pos) const
{ return  (this->*SelectedOutsideMagnetosphere)(pos);
}
////////////////////////////////////////////////////////////////////////////////
//
bool ECRSMagneticField::IGRFOutsideMagnetosphere( G4ThreeVector pos)  const
{ G4double r=pos.mag()/re;
  bool res=false;
  if (r > 25.) res =true;
  return res;
}
////////////////////////////////////////////////////////////////////////////////

bool ECRSMagneticField::TSY89OutsideMagnetosphere( G4ThreeVector pos) const
{ //model from Kobel PhD "Zu die magnetosph�rischen Effekten der kosmischen
  //  Strahlung"
  G4double a=-0.0545;
  G4double b_iopt[]={11.7,11.1,10.8,10.4,10.4,10.2,10.2};
  G4double f_iopt[]={20.,15.,6.67,10.,5.,6.,6.};
  
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos)/re;
  G4double rho2 = pos_gsm.y()*pos_gsm.y() + pos_gsm.z()*pos_gsm.z();
  G4bool res = (rho2>900 || pos_gsm.x() <-60);
  if (!res){ 
  	double  k_par = DipolePS/f_iopt[iopt-1];
  	G4double cos_k=std::cos(k_par);
  	G4double sin_k=std::sin(k_par);
  	G4RotationMatrix mat = G4RotationMatrix(G4ThreeVector(cos_k,0,sin_k),
  					  G4ThreeVector(0,1,0), 
					  G4ThreeVector(-sin_k,0,cos_k));
  	G4ThreeVector pos_rot = mat*pos_gsm;
	G4double rho_rot= pos_rot.y()*pos_rot.y() + pos_rot.z()*pos_rot.z();
	res = (pos_rot.x() > a*rho_rot +b_iopt[iopt]); 
  }
 						  
  return res;
}
///////////////////////////////////////////////////////////////////////////////
//
bool ECRSMagneticField::TSY96OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko96 code
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

  G4ThreeVector pos_geo=pos;
  
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  
  // magnetopause model
  G4double am0=70.;
  G4double s0=1.08;
  G4double x00=5.48;
  G4double dsig=.005;
  G4double xappa=std::pow(Pdyn/2.,0.14);
  G4double x0=x00/xappa;
  G4double am=am0/xappa;
  G4double rho2=ygsm*ygsm+zgsm*zgsm;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
///////////////////////////////////////////////////////////////////////////////
//
bool ECRSMagneticField::TSY2001OutsideMagnetosphere( G4ThreeVector pos) const 
{ 
  // taken form the fortran tsyganenko2001 code
  SpaceCoordinateConvertor* theConvertor =
                   SpaceCoordinateConvertor::getInstance(); 

  G4ThreeVector pos_geo=pos;
  // if (ConsiderDipoleShift) pos_geo-= DipoleSchift;  
  //GEOtoGSM
  G4ThreeVector pos_gsm =theConvertor->TransformGEOinGSM(pos_geo);
  

  G4double xgsm=pos_gsm.x()/re;
  G4double ygsm=pos_gsm.y()/re;
  G4double zgsm=pos_gsm.z()/re;
  G4double r= pos_gsm.mag()/re;
  
  G4double xss=xgsm;
  G4double zss=zgsm;
  
  G4double rh0=8.94335;
  G4double dd;
  G4double rh2= -5.2;
  G4double sps = std::sin(DipolePS);
  //change of coordinate due to the tail warping
  do{ G4double xsold=xss;
      G4double zsold=zss;
      G4double rh=rh0+rh2*(zss/r)*(zss/r);
      G4double sinpsas=sps/std::pow(1.+std::pow(r/rh,3.),1./3.);
      G4double cospsas=std::sqrt(1. - sinpsas*sinpsas);
      zss=xgsm*sinpsas+zgsm*cospsas;
      xss=xgsm*cospsas-zgsm*sinpsas;
      dd=std::abs(xss-xsold)+std::abs(zss-zsold);
  } while (dd > 1.e-6);
  
  // magnetopause model
  G4double a0_a=34.586;
  G4double s0=1.196;
  G4double a0_x0=3.4397;
  G4double dsig=.003;
  G4double xappa=std::pow(Pdyn/2.,0.15790);
  G4double x0=a0_x0/xappa;
  G4double am=a0_a/xappa;
  G4double rho2=ygsm*ygsm+zss*zss;
  G4double asq=am*am;
  G4double xmxm=am+xgsm-x0;
  if (xmxm < 0.) xmxm=0.;
  G4double axx0=xmxm*xmxm;
  G4double aro=asq+rho2;
  G4double sigma=std::sqrt((aro+axx0+std::sqrt((aro+axx0)*(aro+axx0)-4.*asq*axx0))/(2.*asq));
  if (sigma < s0-dsig && xgsm>-50.) return false;
  else return true;
}
////////////////////////////////////////////////////////////////////////////////

void ECRSMagneticField::ComputeTSYParameter(G4double t)
{ G4int index =-1;
  if (n_tsy_data>0){
  	for (G4int i=0;i<n_tsy_data-1;i++){
   		if ( (t >= times_of_data[i]) && (t <= times_of_data[i+1])){
       			index=i;
        		i=n_tsy_data-1;
		}
	}
  	if (index > -1){
    		G4double dx=(t-times_of_data[index])
                	/(times_of_data[index+1]-times_of_data[index]);
     		Pdyn=Pdyn_data[index] + dx * (Pdyn_data[index+1]-Pdyn_data[index]);
     		Dst=Dst_data[index] + dx * (Dst_data[index+1]-Dst_data[index]);
     		Dst*=nanotesla;
		Imfy=Imf_gsm_y_data[index] + dx * 
                      (Imf_gsm_y_data[index+1]-Imf_gsm_y_data[index]);
		Imfy*=nanotesla;      
     		Imfz=Imf_gsm_z_data[index] + dx * 
                      (Imf_gsm_z_data[index+1]-Imf_gsm_z_data[index]);
     		Imfz*=nanotesla;
		G1=g1_tsy_param_data[index] + dx * 
                       (g1_tsy_param_data[index+1]-g1_tsy_param_data[index]);
     		G2=g2_tsy_param_data[index] + dx * 
                      (g2_tsy_param_data[index+1]-g2_tsy_param_data[index]);     
  	}
  }  		
}
////////////////////////////////////////////////////////////////////////////////
//
void ECRSMagneticField::SetStepper(G4String aString)
{ if (fStepper) delete fStepper;
   
  if (aString == "ExplicitEuler"){
    	fStepper = new G4ExplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ExplicitEuler is called"<<G4endl;
  }
  else if  (aString == "ImplicitEuler"){
    	fStepper = new G4ImplicitEuler( fEquationOfMotion );
     	G4cout<<"G4ImplicitEuler is called"<<G4endl;
  }
  else if  (aString == "SimpleRunge"){
    	fStepper = new G4SimpleRunge( fEquationOfMotion );
     	G4cout<<"G4SimpleRunge is called"<<G4endl;
  } 
  else if  (aString == "ClassicalRK4"){
    	fStepper = new G4ClassicalRK4( fEquationOfMotion );
     	G4cout<<"G4ClassicalRK4 is called"<<G4endl;
  }  
  else if  (aString == "RKG3_Stepper"){
    	fStepper = new G4RKG3_Stepper( fEquationOfMotion );
     	G4cout<<"G4RKG3_Stepper is called"<<G4endl;
  }
  else if  (aString == "CashKarpRKF45"){
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }   
  else {
    	fStepper = new G4CashKarpRKF45( fEquationOfMotion );
     	G4cout<<"The selected stepper is not available"<<G4endl;
     	G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
  }
  
  if (fChordFinder) delete fChordFinder;
  G4double min_step= .000001*re;
  fChordFinder = new G4ChordFinder(this,min_step,fStepper);
  fChordFinder->SetDeltaChord(DeltaChord); 
  
  G4FieldManager* fieldMgr= 
  G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMgr->SetChordFinder(fChordFinder); 

}
