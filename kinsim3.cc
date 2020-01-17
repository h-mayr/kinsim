// Code to read in stopping powers from SRIM and plot in ROOT with some random spread
// Liam Gaffney (liam.gaffney@cern.ch) - Originally authord in September 2015
// Follow https://github.com/lpgaff/kinsim for updates and bug fixes
//
// Thanks to Nigel Warr, Kenzo Abrahams, Amar Boukhari, and others at Miniball
// for help with testing, improvements and bug fixes.

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
using namespace std;

TGraph *gSP[4];

string convertInt( int number ) {
	
	stringstream ss;
	ss << number;
	return ss.str();
	
}

string convertFloat( float number, int precision ) {
	
	stringstream ss;
	ss << setprecision( precision ) << number;
	return ss.str();
	
}

const string gElName[110] = {
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
    "Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",
    "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
    "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
    "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
    "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
    "Mt","Ds" };

bool stoppingpowers( int Zb, int Zt, double Ab, double At, string srim_dir, string opt ) {
	
	int index = 0;
    string srimfile = srim_dir + "/"; // prefix
	string title = "Stopping powers for ";
	
	// Beam or target like..?
	if( opt.substr(0,1) == "B" ) {
		
		srimfile += convertInt(Ab+0.5) + gElName[Zb-1];
		title += convertInt(Ab+0.5) + gElName[Zb-1];
		
	}
	
	else if( opt.substr(0,1) == "T" ) {
		
		srimfile += convertInt(At+0.5) + gElName[Zt-1];
		title += convertInt(At+0.5) + gElName[Zt-1];
		index++;
		
	}
	
    else {
		
        cout << "opt must equal BT, TT, BS or TS\n";
        return false;
		
	}
    
	// Target, contaminant or alumium dead layer..?
	if( opt.substr(1,1) == "T" ) {
		
		srimfile += "_" + convertInt(At+0.5) + gElName[Zt-1] + ".txt";
		title += " in " + convertInt(At+0.5) + gElName[Zt-1];
		title += ";Ion energy [MeV];Stopping power [MeV/(mg/cm^2)]";
		
	}
	
	else if( opt.substr(1,1) == "S" ) {
		
		srimfile += "_Si.txt";
		title += " in the Si dead layer";
		title += ";Ion energy [keV];Stopping power [MeV/mm]";
		index += 2;
		
	}
	
	else {
		
		cout << "opt must equal BT, TT, BS or TS\n";
		return false;
		
	}
	
    ifstream infile;
    infile.open( srimfile.c_str(), ios::in );
    
    if( !infile.is_open() ) {
        
        cout << "Cannot open " << srimfile << endl;
        return false;
        
    }
    
	gSP[index]->SetTitle( title.c_str() );
	
    string line, units, tmp_str;
    stringstream line_ss;
    bool endflag = false;
    double BEn, nucl, elec, total, tmp_dbl;
    int p = 0;
    
    // Test file format
    getline( infile, line );
    if( line.substr( 0, 5 ) == " ====" ) {
        
        while( line.substr( 0, 5 ) != "  ---" )
            getline( infile, line );
        
        getline( infile, line ); // read first line of data
        
    }
    
    while( !infile.eof() && !endflag ) {
        
        // Read in data
        line_ss.str("");
        line_ss << line;
        line_ss >> BEn >> units >> nucl >> elec >> tmp_dbl >> tmp_str >> tmp_dbl >> tmp_str;
        
        if( units == "eV" ) BEn *= 1E-6;
        else if( units == "keV" ) BEn *= 1E-3;
        else if( units == "MeV" ) BEn *= 1E0;
        else if( units == "GeV" ) BEn *= 1E3;
        
        total = nucl + elec ; // MeV / ( mg / cm^2 )
        
        gSP[index]->SetPoint( p, BEn, total );
        //cout << p << "\t" << BEn << "\t" << total << endl;
        
        // Get next line
        getline( infile, line );
        p++;
        
        // If we've reached the end, stop
        if( line.substr( 0, 9 ) == "---------" ) endflag = true;
        if( line.substr( 0, 9 ) == " Multiply" ) endflag = true;
        
    }
    
    TCanvas *c = new TCanvas();
    gSP[index]->Draw("A*");
    gSP[index]->GetXaxis()->SetTitleOffset(1.3);
    gSP[index]->GetYaxis()->SetTitleOffset(1.3);
    TGaxis::SetMaxDigits(3);
    string pdfname = srimfile.substr( 0, srimfile.find_last_of(".") ) + ".pdf";
    c->SetLogx();
    c->SaveAs( pdfname.c_str() );
	
	delete c;
    
    return true;
    
}

bool stoppingpowers( int Zb, int Zt, double Ab, double At, string srim_dir ) {
	
	bool success = true;
	
	for( int i = 0; i < 4; i++ )
		gSP[i] = new TGraph();
	
	success *= stoppingpowers( Zb, Zt, Ab, At, srim_dir, std::string("BT") );
	success *= stoppingpowers( Zb, Zt, Ab, At, srim_dir, std::string("TT") );
	success *= stoppingpowers( Zb, Zt, Ab, At, srim_dir, std::string("BS") );
	success *= stoppingpowers( Zb, Zt, Ab, At, srim_dir, std::string("TS") );
	
	return success;
	
}

double GetTh( double anno, double cd_dist ) {

    // Returns theta angle from ann strip number in radians */
    return ( atan( ( 9 + ( 15.5 - anno ) * 2 ) / cd_dist ) );

}

double projLab( double com, double Ab, double At, double Eb, double Ex ) {

	double tau = Ab/At;
	double Eprime = Eb*Ab - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ab/Eprime);
	
//	if( tau > 1 ) double Th_max = TMath::ASin(tau*epsilon);	
		
	// y = tan(theta_lab)
	double y = TMath::Sin(com) / ( TMath::Cos(com) + tau*epsilon );

	double Th = TMath::ATan(y);
	if( Th < 0. ) Th += TMath::Pi();
	return TMath::RadToDeg()*Th;
	
}

double targLab( double com, double Ab, double At, double Eb, double Ex ) {

	double tau = Ab/At;
	double Eprime = Eb*Ab - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ab/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Sin(TMath::Pi()-com) / ( TMath::Cos(TMath::Pi()-com) + epsilon );

	double Th = TMath::ATan(y);
	if( Th < 0. ) Th += TMath::Pi();
	return TMath::RadToDeg()*Th;

}

double projCoM( double theta_lab, double Ab, double At, double Eb, double Ex ) {

	double tau = Ab/At;
	double Eprime = Eb*Ab - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ab/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Tan(theta_lab);
	// x = cos(com)
	double x = (-y*y*epsilon*tau + TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);

    double Th;
    if( theta_lab < 0.5*TMath::Pi() ) Th = TMath::ACos(x);
    else Th = TMath::Pi() - TMath::ACos(x);
    if( Th < 0. ) Th += TMath::Pi();
    return TMath::RadToDeg()*Th;

}

double targCoM( double theta_lab, double Ab, double At, double Eb, double Ex ) {

	double tau = Ab/At;
	double Eprime = Eb*Ab - Ex*(1+tau);
	double epsilon = TMath::Sqrt(Eb*Ab/Eprime);

	// y = tan(theta_lab)
	double y = TMath::Tan(theta_lab);
	// x = cos(com)
	double x = (-y*y*epsilon*tau + TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);

	double Th = TMath::ACos(x);
	if( Th < 0. ) Th += TMath::Pi();
	return TMath::RadToDeg()*Th;

}

double projEn( double Ab, double At, double BEn, double Ex, double th_cm ) {
    
    double Eprime = BEn - ( Ex * ( 1 + (Ab/At) ) );
    double tau = (Ab/At) * TMath::Sqrt( BEn / Eprime );
    
    double Eproj = TMath::Power( At/(At+Ab), 2.0 );
    Eproj *= 1. + tau*tau + 2.*tau*TMath::Cos( th_cm );
    Eproj *= Eprime;
    
    return Eproj;
    
}

double targEn( double Ab, double At, double BEn, double Ex, double th_cm ) {
    
    double Eprime = BEn - ( Ex * ( 1 + (Ab/At) ) );
    double tau = (Ab/At) * TMath::Sqrt( BEn / Eprime );
    double epsilon = TMath::Sqrt( BEn / Eprime );
    
    double Etarg = (At*Ab) / TMath::Power( (At+Ab), 2.0 );
    Etarg *= 1. + epsilon*epsilon + 2.*epsilon*TMath::Cos( TMath::Pi() - th_cm );
    Etarg *= Eprime;
    
    return Etarg;

}

double GetELoss( float Ei, float dist, int opt, string combo ) {
	
	// Returns the energy loss at a given initial energy and distance travelled in the target or Si dead layer
	// Ei is the initial energy in MeV
	// dist is the distance travelled in the target in mg/cm2
	// opt = 0 calculates normal energy loss as particle moves through target (default)
	// opt = 1 calculates energy increase, i.e. tracing particle back to reaction point
	// combo = "BT", "TT", "BS" or "TS" for the beam in target, target in target,
	// beam in Si or target in Si, respectively.
	// Stopping power data is taken from SRIM the output files must be placed in the './srim/'
	// folder with the format 62Fe_109Ag.txt, 62Fe_Si.txt, 109Ag_109Ag.txt or 109Ag_Si.txt,
	// for combo = "BT", "TT", "BS" and "TS", repsectively.

	double dedx = 0;
	int Nmeshpoints = 20; // number of steps to take in integration
	double dx = dist/(double)Nmeshpoints;
	double E = Ei;
	
	for( int i = 0; i < Nmeshpoints; i++ ){
		
		if( E < 1. ) break; // when we fall below 1 MeV we assume maximum energy loss

		if( combo == "BT" ) dedx = gSP[0]->Eval(E);
		else if( combo == "TT" ) dedx = gSP[1]->Eval(E);
		else if( combo == "BS" ) dedx = gSP[2]->Eval(E);
		else if( combo == "TS" ) dedx = gSP[3]->Eval(E);
		
		if( opt == 1 )
			E += dedx*dx;
		
		else
			E -= dedx*dx;
		
	}
	
	
	if( opt == 0 ) return Ei - E;
	else return E - Ei;
	
}

double GetTEn( double Ab, double At, double Eb, double Ex, double TTh, double th_cm, double thick, double depth ) {

	// Returns energy of target after exiting the target

    if( TTh < 0.501*TMath::Pi() && TTh > 0.499*TMath::Pi() ) return 0.1; // stopped
    
    // energy at interaction point
    double Ereac = Ab*Eb - GetELoss( Ab*Eb, depth, 0, "BT" );
    
    // energy after reaction
    double Etarg = targEn( Ab, At, Ereac, Ex, th_cm );
    
    // energy loss over distance to exit of target
    double dist = TMath::Abs( (double)(thick-depth) / TMath::Cos( TTh ) );
    Etarg -= GetELoss( Etarg, dist, 0, "TT" );
    
    if( Etarg < 0. ) return 0.1; // recoil is stopped in target
	
	// Correct for dead layer loss
	dist = TMath::Abs( 0.0007 / TMath::Cos( TTh ) );
	Etarg -= GetELoss( Etarg, dist, 0, "TS" );

    return Etarg;
    
}

double GetBEn( double Ab, double At, double Eb, double Ex, double BTh, double th_cm, double thick, double depth ) {

	// Returns energy of target after exiting the target
    if( BTh < 0.501*TMath::Pi() && BTh > 0.499*TMath::Pi() ) return 0.1; // stopped
    
    // energy at interaction point
    double Ereac = Ab*Eb - GetELoss( Ab*Eb, depth, 0, "BT" );
    
    // energy after reaction
    double Eproj = projEn( Ab, At, Ereac, Ex, th_cm );
    
    // energy loss over distance to exit of target
    double dist = TMath::Abs( (double)(thick-depth) / TMath::Cos( BTh ) );
    Eproj -= GetELoss( Eproj, dist, 0, "BT" );
    
    if( Eproj < 0. ) return 0.1; // beam is stopped in target
    
	// Correct for dead layer loss
	dist = TMath::Abs( 0.0007 / TMath::Cos( BTh ) );
	Eproj -= GetELoss( Eproj, dist, 0, "BS" );

    return Eproj;
    
}

void kinsim3( int Zb, int Zt, double Ab, double At, double thick /* mg/cm^2 */, double Eb /* MeV/u */,
    double dEb = 0.1 /* MeV/u */, double Ex = 1.0 /* MeV */, double res = 0.6 /* % */,
	double cd_dist = 28.0 /* mm */, bool flat = false /* angular distribution? */,
	long Nevts = 1E6, string srim_dir = "../srim" ) {
		
    // Suppress some message from root
    gErrorIgnoreLevel = kWarning;
    
    // Check we have sensible elements
    if( Zb > 110 || Zt > 110 ) {
        cout << "Super heavy elements!" << endl;
        return;
    }
    
    // Setup stopping powers
	for( int i = 0; i < 4; i++ ) gSP[i] = new TGraph();
	if( !stoppingpowers( Zb, Zt, Ab, At, srim_dir ) )
        return;

	// Open output file
	string outname = convertInt(Ab+0.5) + gElName[Zb-1] + "_" + convertInt(At+0.5) + gElName[Zt-1] + "_";
    outname += convertFloat(thick,3) + "mg_" + convertFloat(Eb,3) + "MeVu_d";
	outname += convertFloat(dEb,3) + "MeVu_res" + convertFloat(res,1) + ".root";
	TFile *out = new TFile(outname.c_str(),"RECREATE");

	// Define and initiate histograms to fill
	double stepSize = 1.0; // degrees
    double cd_angles[17];
    for( int k=0; k<17; k++ )
        cd_angles[k] = GetTh( 15.5 - k, cd_dist ) * TMath::RadToDeg();

	string title = "Kinematics in the lab frame for " + convertInt(Ab+0.5) + gElName[Zb-1] + " on ";
	title += convertInt(At+0.5) + gElName[Zt-1] + " at " + convertFloat(Eb,3) + " MeV/u";
	string title1 = title + ";Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab = new TH2F("kin_lab",title1.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
	string title2 = title + " (beam);Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab_b = new TH2F("kin_lab_b",title2.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
	string title3 = title + " (recoil);Laboratory angle [deg];Energy [MeV]";
	TH2F *kin_lab_t = new TH2F("kin_lab_t",title3.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
    string title7 = title + ";Lab angle of recoil [deg];Lab angle of beam [deg]";
    TH2F *lab_lab = new TH2F("lab_lab",title7.c_str(),(int)(180./stepSize),0,180,(int)(180./stepSize),0,180);
    string title8 = title + ";Laboratory angle [deg];Energy [MeV]";
    TH2F *cd_sim = new TH2F("cd_sim",title8.c_str(),16,cd_angles,1000,0,1000);

    title = "Kinematics in the CoM frame for " + convertInt(Ab+0.5) + gElName[Zb-1] + " on ";
    title += convertInt(At+0.5) + gElName[Zt-1] + " at " + convertFloat(Eb,3) + " MeV/u";
    string title4 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com = new TH2F("kin_com",title4.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
	string title5 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com_b = new TH2F("kin_com_b",title5.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
	string title6 = title + ";Centre of mass angle [deg];Energy [MeV]";
	TH2F *kin_com_t = new TH2F("kin_com_t",title6.c_str(),(int)(180./stepSize),0,180,1000,0,1000);
	
	// Define and initiate Rutherford distribution
	string eqnR = "1.44*((";
	eqnR += convertInt(Zb) + "*" + convertInt(Zt) + ")/" + convertInt((int)(Eb*Ab)) + ")**2";
	eqnR += "/(sin(x*pi/360.)**4)";
	TF1 *ruth = new TF1("ruth",eqnR.c_str(),1.0,180.0);
	TGraph *gRuth = new TGraph(ruth);
	gRuth->SetTitle("Rutherford cross-section;Centre of mass angle [deg];d#sigma_{R}/d#Omega");
	
	// Define and initiate Coulex probability
	TGraph *gClxp = new TGraph();
	gClxp->SetTitle("Coulex probability;Centre of mass angle [deg];P_{CE}");	
	gClxp->SetPoint(0,  0.0   ,0.000000);
	gClxp->SetPoint(1,  5.0   ,0.000000);
	gClxp->SetPoint(2,  10.0  ,0.000001);
	gClxp->SetPoint(3,  16.0  ,0.000013);
	gClxp->SetPoint(4,  22.0  ,0.0001);
	gClxp->SetPoint(5,  28.0  ,0.0006);
	gClxp->SetPoint(6,  34.0  ,0.0020);
	gClxp->SetPoint(7,  40.0  ,0.0046);
	gClxp->SetPoint(8,  60.0  ,0.0234);
	gClxp->SetPoint(9,  80.0  ,0.0550);
	gClxp->SetPoint(10, 100.0 ,0.0900);
	gClxp->SetPoint(11, 120.0 ,0.1198);
	gClxp->SetPoint(12, 140.0 ,0.1400);
	gClxp->SetPoint(13, 160.0 ,0.1507);
	gClxp->SetPoint(14, 180.0 ,0.1539);
	
	// Define and initiate Coulex cross-section
    TH1F *hClx = new TH1F( "hClx", "hClx", 200, 0, 180 );
	TGraph *gClx = new TGraph();
	gClx->SetTitle("Coulex cross section;Centre of mass angle [deg];d#sigma_{CE}/d#Omega");
	double P_CE, dsigma_R, dsigma_CE, ang;

	for( int k=0; k<200; k++ ) {
	
		ang = 0.0000001 + 180.*k/200.;
		dsigma_R = gRuth->Eval( ang, 0, "S" );
		P_CE = gClxp->Eval( ang, 0, "S" );
		if( P_CE < 1E-06 ) P_CE = 0;
		dsigma_CE = P_CE * dsigma_R;
		gClx->SetPoint( k, ang, dsigma_CE );
        hClx->SetBinContent( k+1, dsigma_CE );

	}
    
	// Write graphs to file
	gRuth->Write("gRuth");
	gClxp->Write("gClxp");
	gClx->Write("gClx");

	// Some parameters needed for filling
	double com, b_lab, b_en, t_lab, t_en, depth, Eb_real;
	TRandom3 rand;
	
	// Loop over number of events
	for( int i=0; i<Nevts; i++ ){
	
        if( (i+1)%10000 == 0 ) {
            cout << "\t" << i+1 << "/" << Nevts << " - " << (int)((i+1)*100./Nevts) << "\%\r";
            cout.flush();
        }
	
		if( flat ) com = 180.0 * rand.Rndm(i);
		else com = hClx->GetRandom();
        depth = rand.Rndm(i) * thick;
        Eb_real = Eb + rand.Gaus( 0, dEb );
		
		b_lab = projLab( com*TMath::DegToRad(), Ab, At, Eb_real, Ex );
		t_lab = targLab( com*TMath::DegToRad(), Ab, At, Eb_real, Ex );

        b_en = GetBEn( Ab, At, Eb_real, Ex, b_lab*TMath::DegToRad(), com*TMath::DegToRad(), thick, depth );
		t_en = GetTEn( Ab, At, Eb_real, Ex, t_lab*TMath::DegToRad(), com*TMath::DegToRad(), thick, depth );
        b_en += rand.Gaus( 0, res*b_en*0.01 ); // detector resolution %
        t_en += rand.Gaus( 0, res*t_en*0.01 ); // detector resolution %
		
		lab_lab->Fill( b_lab, t_lab );
        kin_lab_b->Fill( b_lab, b_en );
        kin_lab_t->Fill( t_lab, t_en );
        cd_sim->Fill( b_lab, b_en );
        cd_sim->Fill( t_lab, t_en );
		kin_com_b->Fill( com, b_en );
		kin_com_t->Fill( com, t_en );
		
	}
	
    cout << endl;
	
	kin_lab->Add( kin_lab_b, kin_lab_t );
	kin_com->Add( kin_com_b, kin_com_t );
	
    string name;
    for( int i = 0; i < cd_sim->GetNbinsX(); i++ ) {
    
        name = "cd_sim_" + convertInt(i+1);
        cd_sim->ProjectionY( name.c_str(), i+1, i+1 );
        
    }
    
	out->Write();
	//out->Close();

}
