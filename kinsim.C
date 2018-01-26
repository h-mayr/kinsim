// Code to read in the kinematics from LISE++ and plot in ROOT with some random spread

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph"
#include "TRandom3.h"


#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

string convertInt( int number ) {

	stringstream ss;
	ss << number;
	return ss.str();

}

double projLab( double com, double Ap, double At, double Eb, double Ex ) {

	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);
	
//	if( tau > 1 ) double Th_max = TMath::ASin(tau*epsilon);	
		
	// y = tan(theta_lab)
	double y = TMath::Sin(com) / ( TMath::Cos(com) + tau*epsilon );

	double Th = TMath::ATan(y);
	if( Th < 0. ) Th += TMath::Pi();
	return TMath::RadToDeg()*Th;
	
}

double targLab( double com, double Ap, double At, double Eb, double Ex ) {

	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Sin(TMath::Pi()-com) / ( TMath::Cos(TMath::Pi()-com) + epsilon );

	double Th = TMath::ATan(y);
	if( Th < 0. ) Th += TMath::Pi();
	return TMath::RadToDeg()*Th;

}

double projCoM( double theta_lab, double Ap, double At, double Eb, double Ex ) {

	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Tan(theta_lab);
	// x = cos(com)
	double x = (-y*y*epsilon*tau + TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);

	double Th = TMath::ACos(x);
	if( Th < 0. ) Th += 2*TMath::Pi();
	return TMath::RadToDeg()*Th;

}

double targCoM( double theta_lab, double Ap, double At, double Eb, double Ex ) {

	double tau = Ap/At;
	double Eprime = Eb*Ap - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ap/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Tan(theta_lab);
	// x = cos(com)
	double x = (-y*y*epsilon*tau + TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);

	double Th = TMath::ACos(x);
	if( Th < 0. ) Th += 2*TMath::Pi();
	return TMath::RadToDeg()*Th;

}

void kinsim( string inputfile, int Zp, int Zt, double Ap, double At, double Eb /* MeV/u */,
					double Ex /* MeV */, double sigma = 20. /* MeV */, int Nevts = 1000000 ) {

	// Load library for spherical harmonics
//	gSystem->Load("libMathMore");
	
	// Open input file
	ifstream in;
	in.open(inputfile.c_str(), ios::in);
	if( !in.is_open() ) return;
	
	// Open output file
	string outname = convertInt(Ap) + "p_" + convertInt(At) + "t_" + convertInt((int)(Eb*1000));
	outname += "AkeV_sig" + convertInt((int)sigma) + "MeV.root";
	TFile *out = new TFile(outname.c_str(),"RECREATE");

	// Define parameters to read in
	double Thp, Tht, Ep, Et;
	
	// Define and initiate histograms and graphs to fill
	TGraph *gProj = new TGraph();
	TGraph *gTarg = new TGraph();
	gProj->SetTitle("Projectile kinematics;Laboratory angle [deg];Energy [MeV]");
	gTarg->SetTitle("Recoil kinematics;Laboratory angle [deg];Energy [MeV]");
	TGraph *gRuth;

	double stepSize = 1.0; // degrees
	string title = "Kinematics in the lab frame for A_{p}=" + convertInt(Ap) + " and A_{t}="; 
	title += convertInt(At) + " at " + convertInt((int)(Eb*1000)) + " keV/u";
	string title1 = title + " (projectile);Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab = new TH2F("kin_lab",title1.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title2 = title + " (projectile);Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab_p = new TH2F("kin_lab_p",title2.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title3 = title + " (recoil);Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab_t = new TH2F("kin_lab_t",title3.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title4 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com = new TH2F("kin_com",title4.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title5 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com_p = new TH2F("kin_com_p",title5.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title6 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com_t = new TH2F("kin_com_t",title6.c_str(),(int)(180./stepSize),0,180,250,0,1000);
	string title7 = title + ";Lab angle of recoil [deg];Lab angle of projectile [deg]";
	TH2F *lab_lab = new TH2F("lab_lab",title7.c_str(),(int)(180./stepSize),0,180,(int)(180./stepSize),0,180);
	
	// Define and initiate Rutherford distribution
	string eqnR = "1.44*((";
	eqnR += convertInt(Zp) + "*" + convertInt(Zt) + ")/" + convertInt((int)(Eb*Ap)) + ")**2";
	eqnR += "/(sin(x*pi/360.)**4)";
	TF1 *ruth = new TF1("ruth",eqnR.c_str(),1.0,180.0);
	gRuth = new TGraph(ruth);
	gRuth->SetTitle("Rutherford cross-section;Centre of mass angle [deg];d#sigma_{R}/d#Omega");
	
	// Define and initiate Coulex probability
	gClxp = new TGraph();
	gClxp->SetTitle("Coulex probability;Centre of mass angle [deg];P_{CE}");	
	gClxp->SetPoint(0, 0.0   ,0.0000);	gClxp->SetPoint(10, 10.0   ,0.0000);
	gClxp->SetPoint(1, 1.0   ,0.0000);	gClxp->SetPoint(11, 11.0   ,0.0000);
	gClxp->SetPoint(2, 2.0   ,0.0000);	gClxp->SetPoint(12, 12.0   ,0.0000);
	gClxp->SetPoint(3, 3.0   ,0.0000);	gClxp->SetPoint(13, 13.0   ,0.0000);
	gClxp->SetPoint(4, 4.0   ,0.0000);	gClxp->SetPoint(14, 14.0   ,0.0000);
	gClxp->SetPoint(5, 5.0   ,0.0000);	gClxp->SetPoint(15, 15.0   ,0.0000);
	gClxp->SetPoint(6, 6.0   ,0.0000);	gClxp->SetPoint(16, 16.0   ,0.0000);
	gClxp->SetPoint(7, 7.0   ,0.0000);	gClxp->SetPoint(17, 17.0   ,0.0000);
	gClxp->SetPoint(8, 8.0   ,0.0000);	gClxp->SetPoint(18, 18.0   ,0.0000);
	gClxp->SetPoint(9, 9.0   ,0.0000);	gClxp->SetPoint(19, 19.0   ,0.0000);
	gClxp->SetPoint(20, 20.0  ,0.0001);
	gClxp->SetPoint(21, 40.0  ,0.0046);
	gClxp->SetPoint(22, 60.0  ,0.0234);
	gClxp->SetPoint(23, 80.0  ,0.0550);
	gClxp->SetPoint(24, 100.0 ,0.0900);
	gClxp->SetPoint(25, 120.0 ,0.1198);
	gClxp->SetPoint(26, 140.0 ,0.1400);
	gClxp->SetPoint(27, 160.0 ,0.1507);
	gClxp->SetPoint(28, 180.0 ,0.1539);
	
	// Define and initiate Coulex cross-section
	gClx = new TGraph();
	gClx->SetTitle("Coulex cross section;Centre of mass angle [deg];d#sigma_{CE}/d#Omega");
	double P_CE, dsigma_R, dsigma_CE, ang;

	cout << "\t\tdsigma_CE = P_CE x dsigma_R\n";
	
	for( int k=0; k<200; k++ ) {
	
		ang = 1.0000001 + 179.*k/200.;
		dsigma_R = gRuth->Eval(ang,0,"S");
		P_CE = gClxp->Eval(ang,0,"S");
		if( P_CE < 5E-05 ) P_CE = 0;
		dsigma_CE = P_CE * dsigma_R;
		gClx->SetPoint( k, ang, dsigma_CE );
		
		if( k%8==0 )
			cout << "At " << ang << " deg\t" << dsigma_CE << " = " << P_CE << " x " << dsigma_R << endl;
		
	}

	// Read all lines in file and fill graphs
	int i = 1;
	while( in >> Thp >> Ep >> Tht >> Et ) {
		
		gProj->SetPoint(i, Thp, Ap*Ep);
		gTarg->SetPoint(i, Tht, At*Et);
		i++;
		
	}
	in.close();

	// Write graphs to file
	gProj->Write("gProj");
	gTarg->Write("gTarg");
	gRuth->Write("gRuth");
	gClxp->Write("gClxp");
	gClx->Write("gClx");

	// Some parameters needed for filling
	double com, p_lab, p_en, t_lab, t_en;
	TRandom3 rand;
	double w;  // weight of fill, depending on 
	
	// CD detector info
	double CDlow = 15.2; // Angular limits of the CD detector in the lab (degrees)
	double CDupp = 52.0;
	long long pEvts = 0, tEvts = 0;  // counter for number of projectiles and recoils in the CD
	
	// Loop over number of events
	for( int i=0; i<Nevts; i++ ){
	
		if( i%100000 == 0 ) cout << "\t" << i << "/" << Nevts << " - " << i*100./Nevts << "\%" << endl;
	
		com = 1.0000001 + 179.*rand.Rndm();
		w = gClx->Eval(com,0,"S");
		
		p_lab = projLab( com*TMath::DegToRad(), Ap, At, Eb, Ex );
		t_lab = targLab( com*TMath::DegToRad(), Ap, At, Eb, Ex );

		p_en =  gProj->Eval(p_lab);
		t_en =  gTarg->Eval(t_lab);
		p_en += rand.Gaus(0,sigma);
		t_en += rand.Gaus(0,sigma);
		
		if( p_lab >= CDlow && p_lab <= CDupp ) pEvts++;
		if( t_lab >= CDlow && t_lab <= CDupp ) tEvts++;
				
		lab_lab->Fill(p_lab,t_lab,w);
		kin_lab_p->Fill(p_lab,p_en,w);
		kin_lab_t->Fill(t_lab,t_en,w);
		kin_com_p->Fill(com*TMath::RadToDeg(),p_en,w);
		kin_com_t->Fill(com*TMath::RadToDeg(),t_en,w);
		
	}
	
	kin_lab->Add(kin_lab_p,kin_lab_t);
	kin_com->Add(kin_com_p,kin_com_t);
	
	cout << "Fraction of recoil events in CD range: " << (float)tEvts*100./(float)(pEvts+tEvts) << "\%\n";

	out->Write();
	//out->Close();

}