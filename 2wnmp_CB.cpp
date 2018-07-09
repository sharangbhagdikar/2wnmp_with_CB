#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <map>
#include <sstream>
using namespace std;

double randnu ();
double roundin (double);
int rand_l_lim(int comp_size);
int rand_w_lim(int comp_size);
int rand_tox_lim(int comp_size);
void setzero(double &val);//??
void carry_out_reaction();//??
int kil;
int end_of_str_traps;
double k31, k13, mul_fac;
double E1, E2v, E2c, E2g, R, SHW;
int num_compartment, domains_per_nm;
int j;
int nr1, nr2, nf1, nf2, ndz1, ndz2;
double q1, q2, q3;
int c1, c2, c3, c4;
int k1, k2, k3;
double somevar, somevar2, somevar3;
double alpha0, yemp, tau;
double kc1,kc2;
double temp13, temp11, temp12;
double width, length, height;
double delVot = 0,delVt_prev;
int w_lim, l_lim, z_lim;
double psi_str, FSiO2;
int emi, emiprev, previous_emi, previous_emiH2;
double str_Nit[10000];
double str_time[10000];
double t,t_prev,t_base;
int prev,curr;
int counter, f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16;
double r, r1, r2;
double kk1, kk2, kk3;
double alpha_k13_total;
double alpha_k31_total;
double temp_1, temp_2;
double err;

// First vector holds value, second holds cumulative probabiility
pair < vector<double>, vector <double> > E1_arr, E2v_arr, E2c_arr, E2g_arr, SHW_arr, R_arr;
// Stores chosen parameter values for reuse during recovery
vector <double> E1_store, E2v_store, E2c_store, E2g_store, SHW_store, R_store;

vector <double> gauss_arr, avg_vth_str, avg_vth_rec, time_arr;
vector <int> avg_traps_str, avg_traps_rec;


// Solves for surface potential
#define fpsi(VGS,VFB,psi,bot) VGS - VFB - psi - bot*sqrt( (2*q*ESi*E0*ND) * ( VT*exp(-psi/VT) + psi - VT + exp(-2*phif/VT)*(VT*exp(psi/VT) - psi - VT) ))/cox

// Gaussian
#define gaussian(x,mu,sigma) exp(- pow((x-mu)/sigma,2)/2)/(sqrt(2*pi)*sigma)



struct threeD {

	int x;
	int y;
	int z;

	threeD(int x, int y, int z) : x(x), y(y), z(z) {}
	threeD(const threeD &c) {x = c.x; y = c.y; z = c.z;}

	bool operator<(const threeD &o)  const {
		return ((x < o.x) || ((x == o.x) && (y < o.y)) || ((x == o.x) && (y == o.y) && (z < o.z)));
	}

    bool operator== (const threeD &o) const {
    return ((x == o.x) && (z == o.z) && (y == o.y));
    }

};

struct locate_threeD {
    int a;
    int b;
    int c;

    locate_threeD(int a, int b, int c) : a(a), b(b), c(c) {}

    bool operator () (const threeD &o) const{
        return (o.x==a && o.y==b && o.z==c);
    }
};



/**********************************************************************************/
/********************************    Device Data     ******************************/
/**********************************************************************************/

const double q = 1.6e-19; // Farad, charge on an electron
const double ESi = 11.8; // Relative permittivity of Silicon
const double ESiO2 = 3.9; // Relative permittivity of Silicon-di-oxide
const double E0 = 8.854e-14; // Absolute permittivity in farad/cm
const double ni0 = 1e10; // Intrensic carrier concentration at 300 K in /cm^3
const double Eg = 1.12; // Bandgap of Si at 300 K in eV
const double Kb = 1.380658e-23; // Boltzman constant in J/K

/**********************************************************************************/
/***************************    Changeable Parameters     *************************/
/**********************************************************************************/


    #define WIDTH  50
    #define LENGTH 20
    #define HEIGHT 1

    bool isVthneeded = 0;   // Whether delVth information is required or not
    double compartment_length = 1/2.0;
	double z_comp_length = HEIGHT/2.0;
    int sim_size = 50;  // Number of devices that need to be simulated


    double Temp = 130;  // Celsius
    double Vgs_dc = 2; //DC stess voltage in Volts

    double Vgr_dc = 0;                    //DC recovery voltage in Volts
	double Tox = HEIGHT*1e-7;             //cm
	double dev_area = WIDTH*LENGTH*1E-14; // cm2
	double ND = 1e18;       //  Channel Doping in atoms/cm^3
    double N0 = 5e20;       //  Density of gate insulator traps in /cm^3

    int num_defects = 50;       // Make ir N0*Volume

    double vfb = 0.4992;   // Flatband voltage in Volts

    double E1m = -0.3;     // eV
    double E1s = 0.05;        // eV

    double EV = 0;         // eV
    double EC = 1.12;       // eV

    double E2vm = 0;         // eV
    double E2vs = 0.05;        // eV

    double E2cm = EC;       //eV
    double E2cs = 0.05;     //eV

//  OVER HERE!!!
    double E2gm = 0;       //eV
    double E2gs = 0;     //eV

	double SHWm = 2;      // 5.9eV
	double SHWs = 0.4;       // eV

	double Rm = 0.5;       //0.52
	double Rs = 0;         // eV
    double Gamma = 9e-9; // eV cm/V

/*********************************************************************************/

	double VT1 = (Kb*300)/q;    // Thermal voltage Vt@ 300 K = 0.0259eV in eV
	double VT= (Kb*(273+Temp))/q;   // Thermal voltage Vt@ 273+Temp in eV
	double ni = ni0*(pow(static_cast<double>((273+Temp)/300),1.5))*exp(-Eg/2*(1/VT-1/VT1)); // in cm^-3
	double phif=VT*log(ND/ni);                // Bulk potential in V
	double beta = 1/VT;                       // in 1/eV ;34.7meV for 130 deg C

	double kt=Kb*(273+Temp);// Joules
	double cox=(ESiO2*E0)/Tox; // in F/cm^2
	double Ldp=sqrt((ESi*E0*kt)/(q*q*ND)); // Debye Length in cm

    double pi=22.0/7;
	double mp = 0.81*9.1e-31;                   // hole mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double me = 1.08*9.1e-31;                  // electron mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double vp=sqrt(8*kt/(pi*mp))*100;         // should be thermal velocity (check the formula) in cm/s
    double sp=(1e-15) ;                       // cross-section area in cm^2
    double h=6.62607004e-34;                  // in Joule*seconds
    double NV=2*(pow(static_cast<double>(2*pi*mp*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double NC=2*(pow(static_cast<double>(2*pi*me*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double p0=(pow(ni,2))/ND;                 // Bulk minority concentration  /cm^3

// Initial sites allocated randomly and Defect sites upon hole capture
vector <bool> init_sites(num_defects,1);

// Maps each trap to some integer value for convenience
map<threeD, int> threeD2int;
vector<threeD> int2threeD;
vector <int> for_traps, rev_traps;
//Lookup table for rXn rates in all cells
vector <double> k13_arr(num_defects,0), k31_arr(num_defects,0);
// Cumulative rXn rates
vector <double> cum_for(num_defects,0), cum_rev(num_defects,0);

void generate_rand_params()
{
    int index;
    int size_arr = E1_arr.first.size() - 1;
    double g1;

    if(E1s==0) E1 = E1m;

    else
    {
        g1 = randnu();

        index = lower_bound(E1_arr.second.begin(), E1_arr.second.end(), g1) - E1_arr.second.begin();
        E1 = abs(E1_arr.first[index]) < 1e-8 ? 0 : E1_arr.first[index];
    }

    if(E2vs==0) E2v = E2vm;
    else
    {
        g1 = randnu();

        index = lower_bound(E2v_arr.second.begin(), E2v_arr.second.end(), g1) - E2v_arr.second.begin();
        E2v = abs(E2v_arr.first[index]) < 1e-8 ? 0 : E2v_arr.first[index];
    }

    if(E2cs==0) E2c = E2cm;
    else
    {
        g1 = randnu();

        index = lower_bound(E2c_arr.second.begin(), E2c_arr.second.end(), g1) - E2c_arr.second.begin();
        E2c = abs(E2c_arr.first[index]) < 1e-8 ? 0 : E2c_arr.first[index];
    }

    if(SHWs==0) SHW = SHWm;
    else
    {
        g1 = randnu();

        index = lower_bound(SHW_arr.second.begin(), SHW_arr.second.end(), g1) - SHW_arr.second.begin();
        SHW = abs(SHW_arr.first[index]) < 1e-8 ? 0 : SHW_arr.first[index];
    }

    if(Rs==0) R = Rm;
    else
    {

        g1 = randnu();

        index = lower_bound(R_arr.second.begin(), R_arr.second.end(), g1) - R_arr.second.begin();
        R = abs(R_arr.first[index]) < 1e-8 ? 0 : R_arr.first[index];
    }

}

void ox_field(double VGb_str)
{

//    if(VGb_str == 1.5) FSiO2 = 10e6;
//    else if(VGb_str == 0) FSiO2 = 1e6;
	//Solving for surface potential
	int cnt=1;

	int maxiter = 2500;
	double lpsi;
	double rpsi;
    int b1;

	if(VGb_str >= vfb)
    {
        lpsi = 0;
        rpsi = 3*phif;
        b1 = 1;
    }

    else
    {
        lpsi = -2*phif;
        rpsi = 0;
        b1 = -1;
    }
    //cout<<"rpsi "<<rpsi<<endl;
    //cout<<"lpsi "<<lpsi<<endl;

	double fvall = fpsi(VGb_str,vfb,lpsi,b1);
    //cout<<"fvall "<<fvall<<endl;

	double fvalu = fpsi(VGb_str,vfb,rpsi,b1);
    //cout<<"fvalu "<<fvalu<<endl;
    //cout<<VGb_str<<endl;
    double newpsi;
	double tmp=0.0;
	double tol=1e-10;
	double fvalleft,fvalright;

	if (fvall*fvalu < 0.0)
	{
	    //cout<<"hello"<<endl;
		newpsi = (lpsi+rpsi)/2;
		do {
			fvalleft  = fpsi(VGb_str,vfb,lpsi,b1);
			fvalright = fpsi(VGb_str,vfb,newpsi,b1);

			if (fvalleft*fvalright < 0)
					rpsi=newpsi;
				else
					lpsi=newpsi;

            tmp = (lpsi+rpsi)/2;
            newpsi=tmp;
            cnt=cnt+1;
		}

		while(cnt < maxiter && fabs(fvalright) > tol );
		psi_str=newpsi;
	}

    //cout<<psi_str<<endl;

    //psi_str = 1.129843;
	//Surface Electric Field
	double FSi = sqrt(2)*(VT/Ldp)*sqrt((exp(-psi_str/VT)+(psi_str/VT)-1)+((exp(-2*phif/VT))*(exp(psi_str/VT)-(psi_str/VT)-1)));

	//Electric Field in SiO2
	FSiO2 = FSi*ESi/ESiO2;
	//cout<<FSiO2<<endl;
}

void rate(double xvecbytox)
{
    double k13v, k13c, k13g, k31v, k31c, k31g;
    double ps=p0*(exp(beta*psi_str)); // Surface minority carrier concentration near the interface which is inversion carrier density in  /cm^3
    double ns = ND*p0/ps;

    double pf_ps = ps*vp*sp;
    double pf_nv = NV*vp*sp;
    double pf_ns = ns*vp*sp;
    double pf_nc = NC*vp*sp;

    double delET = E1-EV;

    // CHECK SIGN!
    double epsTV = E2v - delET - Gamma*(xvecbytox + 0.5/z_lim)*FSiO2;   // Need to check sign
    double epsTC = epsTV - E2v + E2c;   // Need to check sign

    double eps12v = SHW/(pow((1 + R),2)) + R*epsTV/(1+R) + R*(pow(static_cast<double>(epsTV),2))/(4*SHW);
    double eps21v = eps12v - epsTV;
    double eps12c = SHW/(pow((1 + R),2)) + R*epsTC/(1+R) + R*(pow(static_cast<double>(epsTC),2))/(4*SHW);
    double eps21c = eps12c - epsTC;

    k13v = eps12v > 0 ? pf_ps*exp(-beta*eps12v) : pf_ps;
    k31v = eps21v > 0 ? pf_nv*exp(-beta*eps21v) : pf_nv;

    k13c = eps12c > 0 ? pf_nc*exp(-beta*eps12c) : pf_nc;
    k31c = eps21c > 0 ? pf_ns*exp(-beta*eps21c) : pf_ns;

    k13 = k13v + k13c;
    k31 = k31v + k31c;

}


int main(){


	srand (7998);
	//cout<<rand()<<endl;
	double sw_time_1 [] = {1e-6, 1e3};

    //****************Device Dimensions*********************//


	w_lim = (int) ceil(WIDTH/compartment_length);

	l_lim = (int) ceil(LENGTH/compartment_length);

	z_lim =  (int) ceil(HEIGHT/z_comp_length);

    //******************************************************//


	E1_store.resize(num_defects);
	E2v_store.resize(num_defects);
	E2c_store.resize(num_defects);
	R_store.resize(num_defects);
    SHW_store.resize(num_defects);


    // Generating custom filenames
    ostringstream file_name;
    string str_file,rec_file;
    file_name<<"E1_"<<E1m<<"_shw"<<SHWm<<"_R_"<<Rm<<"_"<<Temp<<"C"<<".csv";
    rec_file = "rec_" + file_name.str();
    str_file = "str_" + file_name.str();

	ofstream myfile1,myfile,myfile2,myfile3, myfile4, myfile5;

    //Stress and recovery files
	myfile1.open (rec_file.c_str(), ios::out | ios::trunc);
	myfile.open (str_file.c_str(), ios::out | ios::trunc);

    // Str and rec stat files
    myfile2.open ("Str_data.txt", ios::out | ios::trunc);
    //myfile3.open ("Rec_data.txt", ios::out | ios::trunc);


    if(isVthneeded)
    {
        // Vth files
        myfile4.open(("Vstr_" + file_name.str()).c_str(), ios::out | ios::trunc);
        myfile5.open(("Vrec_" + file_name.str()).c_str(), ios::out | ios::trunc);

    }



    int cn = 0,cnt;
	int flagger = 0;
    int d1 = 0;


/*******************************************************************************
 PDFs for all 4 Gaussian distributions

 Store probabilities of Gaussian distributed RV from mu-3*sigma to mu+3*sigma
 mu = 1 and sigma=1
 6*sigma/50 = 0.12
********************************************************************************/

    somevar = -2;
    somevar2 = 0;
    for(int i = 0; i < 50; i++)
    {
    somevar += 0.12;
    somevar2 += gaussian(somevar,1,1);
    gauss_arr.push_back(somevar2);
}

    for(int i = 0; i<50; i++)
    {
        gauss_arr[i] = gauss_arr[i]/somevar2; // Approximating the gaussian cdf
    }

    if(E1s != 0)
    {
        somevar = E1m-3*E1s;
        somevar3 = 6*E1s/50;

        for(int i = 0; i<50; i++)
        {
            somevar += somevar3;
            E1_arr.first.push_back(somevar);
            E1_arr.second.push_back(gauss_arr[i]);
        }
    }

    if(E2vs != 0)
    {
        somevar = E2vm-3*E2vs;
        somevar3 = 6*E2vs/50;

        for(int i = 0; i<50; i++)
        {
            somevar += somevar3;
            E2v_arr.first.push_back(somevar);
            E2v_arr.second.push_back(gauss_arr[i]);
        }
    }

    if(E2cs != 0)
    {
        somevar = E2cm-3*E2cs;
        somevar3 = 6*E2cs/50;

        for(int i = 0; i<50; i++)
        {
            somevar += somevar3;
            E2c_arr.first.push_back(somevar);
            E2c_arr.second.push_back(gauss_arr[i]);
        }
    }


    if(SHWs != 0)
    {
        somevar = SHWm-3*SHWs;
        somevar3 = 6*SHWs/50;

        for(int i = 0; i<50; i++)
        {
            somevar += somevar3;
            SHW_arr.first.push_back(somevar);
            SHW_arr.second.push_back(gauss_arr[i]);
        }
    }

    if(Rs != 0)
    {
        somevar = Rm-3*Rs;
        somevar3 = 6*Rs/50;

        for(int i = 0; i<50; i++)
        {
            somevar += somevar3;
            R_arr.first.push_back(somevar);
            R_arr.second.push_back(gauss_arr[i]);
        }
    }



    double t_step;
    int pt_dec = 10;           // #Points per decade required

    t_prev = sw_time_1[0];
    cn = 0;

/************************************************************************/
    //Create time vector with pt_dec points per decade spaced//

    //    **LINEAR**

    while (t_prev < sw_time_1[1] - 0.001)
    {
        t_step = t_prev*10.0/pt_dec;

        for(int j = 0; j < int(pt_dec*0.9); j++)    // Only 9/10 needed since starts with 1/10 not zero
        {
            myfile<<t_prev<<",";
            myfile1<<t_prev<<",";
            time_arr.push_back(t_prev);

            if(isVthneeded)
            {
                // Vth files
                myfile4<<t_prev<<",";
                myfile5<<t_prev<<",";
            }

            t_prev = t_prev + t_step;
            cn += 1;
        }

    }


   //    **LOGARITHM**

//    while (t_prev < sw_time_1[1] - 0.001)
//    {
//        t_step = pow(10, 1.0/pt_dec);
//
//        for(int j = 0; j < int(pt_dec*0.9); j++)    // Only 9/10 needed since starts with 1/10 not zero
//        {
//            myfile<<t_prev<<",";
//            myfile1<<t_prev<<",";
//            time_arr.push_back(t_prev);
//
//            // Vth files
//            myfile4<<t_prev<<",";
//            myfile5<<t_prev<<",";

//            t_prev = t_prev*t_step;
//            cn += 1;
//        }
//
//    }

    myfile<<t_prev<<endl;
    myfile1<<t_prev<<endl;
    time_arr.push_back(t_prev);

    if(isVthneeded)
    {
        // Vth files
        myfile4<<t_prev<<endl;
        myfile5<<t_prev<<endl;
        avg_vth_str.assign(cn+1,0);
        avg_vth_rec.assign(cn+1,0);
    }

    // Pre-assign cn+1 sizes: equal to size of time vector generated above
    avg_traps_str.assign(cn+1,0);
    avg_traps_rec.assign(cn+1,0);



for (int ii = 0; ii < sim_size; ii++)
{

    cnt = 0;

//**********************************************************************

	d1 =0;
    alpha_k13_total = 0;                   // Initial total forward propensity: from ground state to transport state
    alpha_k31_total = 0;                  // Initial reverse propensity:

//  Distribute defect sites spatially in the bulk

	for(d1= 0; d1<num_defects; d1++){
		c1 = rand_l_lim(l_lim);
		c2 = rand_w_lim(w_lim);
		c3 = rand_tox_lim(z_lim);

        threeD some3d(c1,c2,c3);

		map<threeD, int>::iterator it;
		it = threeD2int.find(some3d);

		if (it != threeD2int.end())
        {
			d1--;
			continue;
		}

		else
        {
			threeD2int.insert(make_pair(some3d,d1));
			int2threeD.push_back(some3d);
		}
	}



/*********************************************************************************/
/******************************     STRESS     ***********************************/
/*********************************************************************************/
	emi = 0;
    emiprev = 0;
    delVot = 0;
    delVt_prev = 0;
	t = sw_time_1[0];

	counter = 0;
	flagger =0;

    ox_field(Vgs_dc);

    cn = 0;

    // Assign random parameters to each trap/defect location
    // Calculate and store forward and reverse rXn rates for each defect
    for(int j = 0; j < num_defects; j++)
    {
        generate_rand_params();

        // Store these since will need during recovery
        E1_store[j] = E1;
        E2v_store[j] = E2v;
        E2c_store[j] = E2c;
        SHW_store[j] = SHW;
        R_store[j] = R;
        cn += 1;

        rate((double) (int2threeD[j].z)/z_lim);
        //cout<<int2threeD[j].z<<endl;

        k13_arr[j] = k13;
        k31_arr[j] = k31;

    }
    // Create the cumulative rXn propensity array
    partial_sum(k13_arr.begin(), k13_arr.end(), cum_for.begin());

    alpha_k13_total = cum_for.back();
    cout<<alpha_k13_total<<endl;


	while (t < sw_time_1[1] && flagger == 0){

		r1 = randnu();
		while (r1 == 0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total + alpha_k31_total;

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile << "Error"<<endl;
		}

        if(t+tau <= sw_time_1[1])   // Carry rXn only when 'next' time is within end time
        {
            r2 = randnu();
            while (r2 == 0){
                r2 = randnu();
            }
            r2 = r2*alpha0;

            carry_out_reaction();
        }

        if (t == sw_time_1[0])
        {
            t_prev = time_arr[0];
        }

            t = t + tau;
            while (t_prev < t && cnt < time_arr.size()-1)

            {

                myfile<<emiprev<<",";
                avg_traps_str[cnt] += emiprev;

                if(isVthneeded)
                {
                    myfile4<<delVt_prev<<",";
                    avg_vth_str[cnt] += delVt_prev;
                }
                cnt += 1;

                t_prev = time_arr[cnt];

            }

            if(t <= sw_time_1[1])
            {
                emiprev = emi;
                if(isVthneeded)
                {
                    delVt_prev = delVot;
                }

            }

	}
    myfile<<emi<<endl;
    avg_traps_str[cnt] += emi;

    if(isVthneeded)
    {
        myfile4<<delVot<<endl;
        avg_vth_str[cnt] += delVot;
    }

    cnt = 0;
    end_of_str_traps = emi;


    // Note down stats in Str_data file
	myfile2 << emi <<" out of "<< num_defects<<" defects broken in bulk"<<endl;

	for (int j = 0; j < num_defects; j++) {
        if(!init_sites[j])
        {
            myfile2<<int2threeD[j].x<<" "<<int2threeD[j].y<<" "<<int2threeD[j].z<<" : ";
            myfile2<<E1_store[j]<<" "<<SHW_store[j]<<" "<<R_store[j]<<endl;
        }
	}
	myfile2<<endl;


/*********************************************************************************/
/****************************     RECOVERY     ***********************************/
/*********************************************************************************/

	t = sw_time_1[0];

    ox_field(Vgr_dc);
    cn = 0;

    for(int j = 0; j < num_defects; j++)
    {
        E1 = E1_store[j];
        E2v = E2v_store[j];
        E2c = E2c_store[j];
        SHW = SHW_store[j];
        R = R_store[j];

        // Uodate forward and backwrd rates at possible defect sites
        rate((double) (int2threeD[j].z)/z_lim);
        k13_arr[j] = k13;
        k31_arr[j] = k31;

        if(init_sites[j])
        {
            cum_for[j] = j==0 ? k13 : cum_for[j-1] + k13;
            cum_rev[j] = j==0 ? 0 : cum_rev[j-1];
        }

        else
        {
            cum_for[j] = j==0 ? 0 : cum_for[j-1];
            cum_rev[j] = j==0 ? k31 : cum_rev[j-1] + k31;
        }
    }

    alpha_k13_total = cum_for.back();

//    for(int h = 0; h<num_defects; h++)
//    {
//        cout<<init_sites[h]<<" "<<k13_arr[h]<<" "<<cum_for[h]<<" "<<k31_arr[h]<<" "<<cum_rev[h]<<endl;
//    }
    alpha_k31_total = cum_rev.back();





    while (t < sw_time_1[1] && flagger == 0){

		r1 = randnu();
		while (r1 ==0){
			r1 = randnu();
		}

		alpha0 =  alpha_k13_total+ alpha_k31_total;

		tau = (1/alpha0)*log(1/r1);
		if (tau < 0){
			flagger = 1;
			myfile1 << "Error"<<endl;
            return 0;
		}

        if(t+tau <= sw_time_1[1])
        {
            r2 = randnu();
            while (r2 ==0){
                r2 = randnu();
            }
            r2 = r2*alpha0;

            carry_out_reaction();

        }

        if (t == sw_time_1[0])
        {
            t_prev = time_arr[0];
        }

            t = t + tau;

            while(t_prev < t && cnt < time_arr.size()-1)
            {
                myfile1<<emiprev<<",";
                avg_traps_rec[cnt] += emiprev;

                if(isVthneeded)
                {
                    myfile5<<delVt_prev<<",";
                    avg_vth_rec[cnt] += delVt_prev;
                }
                cnt += 1;

                t_prev = time_arr[cnt];
            }

            if(t < sw_time_1[1])
            {
                if(isVthneeded)
                {
                    delVt_prev = delVot;
                }

                emiprev = emi;
            }
    }

    myfile1<<emi<<endl;
    avg_traps_rec[cnt] += emi;



if(isVthneeded)
{
    myfile5<<delVot<<endl;
    avg_vth_rec[cnt] += delVot;
}


    init_sites.assign(num_defects,1);
    cum_for.assign(num_defects,0);
    cum_rev.assign(num_defects,0);
    threeD2int.clear();
    int2threeD.clear();

}


    //Printing out average values at end of file
    for(int h = 0; h < avg_traps_rec.size(); h++)
    {
        myfile1<<avg_traps_rec[h]*1.0/sim_size<<",";
        myfile<<avg_traps_str[h]*1.0/sim_size<<",";

        if(isVthneeded)
        {
            myfile4<<avg_vth_str[h]/sim_size<<",";
            myfile5<<avg_vth_rec[h]/sim_size<<",";
        }

    }

    myfile1.close();
    //myfile3.close();
    myfile.close();

if(isVthneeded)
{
    myfile2.close();
    myfile5.close();
    myfile4.close();
}

    //cout<<avg_traps_rec[82]<<endl;
	return 0;
}

double randnu(){
    return (rand() * 1.0) / (RAND_MAX);
}

int rand_l_lim(int l_lim){
	int randllim = (int) floor((randnu()*(l_lim)));
	if(randllim==l_lim) return l_lim-1;
	return randllim;
}

int rand_w_lim(int w_lim){
	int randwlim = (int) floor((randnu()*(w_lim)));
	if(randwlim==w_lim) return w_lim-1;
	return randwlim;
}

int rand_tox_lim(int z_lim){
	int randtoxlim = (int) floor((randnu()*(z_lim)));
	if(randtoxlim==z_lim) return z_lim-1;
	return randtoxlim;
}

double roundin(double d) {
	return floor(d + 0.5);
}

void setzero(double &val) {
	if (val > -1e-5 && val < 1e-5) {
		val = 0;
	}
}

void carry_out_reaction() {

	int l;
	double rate_change;
	//K31
	if (r2 < alpha_k31_total) {

            l =(int) (lower_bound(cum_rev.begin(), cum_rev.end(), r2) - cum_rev.begin());

            if(init_sites[l] == 1)
            {
                cout<<"error1"<<" "<<l<<endl;
                cout<<cum_rev[l]<<" ";
                for(int h = 0; h < num_defects; h++)
                {
                    cout<<cum_rev[h]<<endl;
                }
                cout<<"r2: "<<r2<<endl;
                cout<<"alpha_k31_total: "<<alpha_k31_total<<endl;
                while(init_sites[l])
                {
                    l--;
                }
                cout<<cum_rev[l]<<endl;
            }

            init_sites[l] = 1;
            rate_change = k31_arr[l];

            // Forward rXn propensity
            transform(cum_for.begin()+l, cum_for.end(), cum_for.begin()+l, bind2nd(plus<double>(),k13_arr[l]));

            // Reverse rXn propensity
            if(--emi != 0)
            {
                if(l==0 || cum_rev[l-1] == 0)
                {
                    while(init_sites[l])
                    {
                        cum_rev[l++] = 0;
                    }
                }
                transform(cum_rev.begin()+l, cum_rev.end(), cum_rev.begin()+l, bind2nd(minus<double>(), rate_change));
            }

            else cum_rev.assign(num_defects, 0);

            alpha_k31_total = cum_rev.back();
            alpha_k13_total = cum_for.back();

            if(isVthneeded)
            {
                if(emi==0) delVot = 0;
                else delVot -= 1000*q*Tox*(1 - (int2threeD[l].z + 0.5)/(1.0*z_lim))/(dev_area*ESiO2*E0);
            }

	}

	//K13
	else if (r2 < alpha0) {

            l = (int) (lower_bound(cum_for.begin(), cum_for.end(), r2 - alpha_k31_total) - cum_for.begin());

            if(init_sites[l] == 0)
            {
                cout<<"error2"<<" "<<l<<endl;
                cout<<cum_for[l]<<" ";
//                for(int h = 0; h < num_defects; h++)
//                {
//                    cout<<cum_for[h]<<endl;
//                }
//                cout<<r2<<endl;
//                cout<<alpha0<<endl;
                while(!init_sites[l])
                {
                    l--;
                }
                cout<<cum_for[l]<<endl;
            }

            init_sites[l] = 0;
            rate_change = k13_arr[l];

            // Reverse rXn propensity
            transform(cum_rev.begin()+l, cum_rev.end(), cum_rev.begin()+l, bind2nd(plus<double>(),k31_arr[l]));

            // Forward rXn propensity
            if(++emi!=num_defects)
            {
                // To correct errors due to numerical precision
                if(l == 0 || cum_for[l-1]==0)
                {
                    while(!init_sites[l])
                    {
                        cum_for[l++] = 0;
                    }
                }
                transform(cum_for.begin()+l, cum_for.end(), cum_for.begin()+l, bind2nd(minus<double>(),rate_change));
            }

            else cum_for.assign(num_defects,0);

            alpha_k31_total = cum_rev.back();
            alpha_k13_total = cum_for.back();

            if(isVthneeded)
            {
                delVot += 1000*q*Tox*(1 - (int2threeD[l].z + 0.5)/(1.0*z_lim))/(dev_area*ESiO2*E0);
            }

            //setzero(delVot);
	}
}
