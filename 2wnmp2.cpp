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
double k21, k12;
double mul_fac;
double E1, E2v, E2c, E2g, R, SHW;
int num_compartment, domains_per_nm;
int j;

double Vg;
int c1, c2, c3, c4;

double somevar, somevar2, somevar3;
double alpha0, yemp, tau;
double temp13, temp11, temp12;
double width, length, height;
double delVot = 0,delVt_prev;
int w_lim, l_lim, z_lim;
double psi, FSiO2;
int emi, emiprev, previous_emi, previous_emiH2;
double str_Nit[10000];
double str_time[10000];
double t,t_prev,t_base;
int prev,curr;
int counter;
double r1, r2;

double cum_total[2] = {};
double alpha_total[2] = {};

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

// Function object for adding/subtracting propensity from cumulative
struct update_prop {           // function object type:
    double k;
    update_prop (double k) : k(k) {}

  void operator()(double &som) {som += k;}
} ;

/**********************************************************************************/
/********************************    Device Data     ******************************/
/**********************************************************************************/

const double q = 1.6e-19; // Farad, charge on an electron
const double ESi = 11.8; // Relative permittivity of Silicon
const double ESiO2 = 3.9; // Relative permittivity of Silicon-di-oxide
const double E0 = 8.854e-14; // Absolute permittivity in farad/cm
const double ni0 = 1e10; // Intrensic carrier concentration at 300 K in /cm^3
const double Eg0 = 1.166; // Bandgap of Si at 300 K in eV
const double Kb = 1.380658e-23; // Boltzman constant in J/K

/**********************************************************************************/
/***************************    Changeable Parameters     *************************/
/**********************************************************************************/
    #define WIDTH  50
    #define LENGTH 20
    #define HEIGHT 1

    bool isVthneeded = 1;   // Whether delVth information is required or not
    bool isCB = 1, isGate = 0;
    bool lin1_log0 = 1;
    double compartment_length = 1/2.0;
    double z_comp_length = HEIGHT/5.0;
    int sim_size = 50;  // Number of devices that need to be simulated
    int pt_dec = 10; //Number of points per decade; Should be multiple of 10

    double Temp = 125;  // Celsius
    double Vgs_dc = 2.5; //DC stess voltage in Volts

    double Vgr_dc = 0;                    //DC recovery voltage in Volts
    double Tox = HEIGHT*1e-7;             //cm
    double dev_area = WIDTH*LENGTH*1E-14; // cm2
    double ND = 1e18;       //  Channel Doping in atoms/cm^3
    double N0 = 5e20;       //  Density of gate insulator traps in /cm^3

    int num_defects = 140;       // Make it N0*Volume

    double vfb = -0.532;   // Flatband voltage in Volts

    double E1m = -1.36;     // eV
    double E1s = 0;        // eV

    double EV = 0;         // eV

    double E2vm = EV;         // eV
    double E2vs = 0;          // eV

    double E2cs = 0;        //eV

//  OVER HERE!!!
    double E2gm = 0;     // eV
    double E2gs = 0;     // eV

    double NVg = 2e23;      //cm^-3
    double ng = 2e22;       //cm^-3

    double SHWm = 7.95;      // 5.9eV
    double SHWs = 0;       // eV

    double Rm = 2.59;       //0.52
    double Rs = 0;          // eV

/*********************************************************************************/
    double Eg = Eg0 - 4.73e-4*pow(Temp,2)/(Temp+636);    // Temperature dependent BG
	double EC = Eg + EV;    //eV
    double E2cm = EC;       //eV


	double VT1 = (Kb*300)/q;    // Thermal voltage Vt@ 300 K = 0.0259eV in eV
	double VT= (Kb*(273+Temp))/q;   // Thermal voltage Vt@ 273+Temp in eV
	double ni = ni0*(pow(static_cast<double>((273+Temp)/300),1.5))*exp(-Eg/2*(1/VT-1/VT1)); // in cm^-3	// T^1.5 dependence?
	double phif = VT*log(ND/ni);                // Bulk potential in V
	double beta = 1/VT;                       // in 1/eV ;34.7meV for 130 deg C

	double kt = Kb*(273+Temp);// Joules
	double cox = (ESiO2*E0)/Tox; // in F/cm^2
	double Ldp = sqrt((ESi*E0*kt)/(q*q*ND)); // Debye Length in cm

    double pi=22.0/7;
    double mp = 0.81*9.1e-31;                   // hole mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double me = 1.08*9.1e-31;                  // electron mass in kg (http://ecee.colorado.edu/~bart/book/effmass.htm)
    double vp = sqrt(8*kt/(mp*pi))*100;         // should be thermal velocity (check the formula) in cm/s
    double vn = sqrt(8*kt/(me*pi))*100;         // should be thermal velocity (check the formula) in cm/s
    double sp = 1e-19 ;                       // cross-section area in cm^2
    double h = 6.62607004e-34;                  // in Joule*seconds
    double NV=2*(pow(static_cast<double>(2*pi*mp*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double NC=2*(pow(static_cast<double>(2*pi*me*Kb*(273+Temp)/(pow(h,2))),1.5))*1e-6; // in /cm^3
    double p0=(pow(ni,2))/ND;                 // Bulk minority concentration  /cm^3

/*********************************************************************************/

// Initial sites allocated randomly and Defect sites upon hole capture
vector <int> state_1, state_2;
// Maps each trap to some integer value for convenience
map<threeD, int> threeD2int;
vector<threeD> int2threeD;
//Lookup table for rXn rates in all cells
vector <double> k12_arr(num_defects,0), k21_arr(k12_arr);
// Cumulative rXn rates
vector <double> cum_12(num_defects,0), cum_21;

void generate_rand_params()
{
    int index;
    int size_arr = E1_arr.first.size() - 1;
    double g1;

    if(E1s==0) E1 = E1m;

    else
    {
        g1 = randnu();

        index = lower_bound(gauss_arr.begin(), gauss_arr.end(), g1) - gauss_arr.begin();
        E1 = E1m - 3*E1s + 6*(index + 1)*E1s/50;
    }

    if(E2vs==0) E2v = E2vm;
    else
    {
        g1 = randnu();

        index = lower_bound(gauss_arr.begin(), gauss_arr.end(), g1) - gauss_arr.begin();
        E2v = E2vm - 3*E2vs + 6*(index+1)*E2vs/50;
    }

    if(E2cs==0) E2c = E2cm;
    else
    {
        g1 = randnu();

        index = lower_bound(gauss_arr.begin(), gauss_arr.end(), g1) - gauss_arr.begin();
        E2c = E2cm - 3*E2cs + 6*(index+1)*E2cs/50;
    }

    if(SHWs==0) SHW = SHWm;
    else
    {
        g1 = randnu();

        index = lower_bound(gauss_arr.begin(), gauss_arr.end(), g1) - gauss_arr.begin();
        SHW = SHWm - 3*SHWs + 6*(index+1)*SHWs/50;
    }

    if(Rs==0) R = Rm;
    else
    {

        g1 = randnu();

        index = lower_bound(gauss_arr.begin(), gauss_arr.end(), g1) - gauss_arr.begin();
        R = Rm - 3*Rs + 6*(index+1)*Rs/50;
    }

}

void ox_field(double VGb_str)
{


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
		psi=newpsi;
	}

    //cout<<psi<<endl;

    //psi = 1.129843;
	//Surface Electric Field
	double FSi = sqrt(2)*(VT/Ldp)*sqrt((exp(-psi/VT)+(psi/VT)-1)+((exp(-2*phif/VT))*(exp(psi/VT)-(psi/VT)-1)));

	//Electric Field in SiO2
	FSiO2 = FSi*ESi/ESiO2;
	//cout<<FSiO2<<endl;
}


void rate(double xvecbytox, double psi, double field)
{
    double k12v, k12c, k12g, k21v, k21c, k21g;
    double ps=p0*(exp(beta*psi)); // Surface minority carrier concentration near the interface which is inversion carrier density in  /cm^3
    double ns = ND*p0/ps;
//    cout<<"p0 "<<p0<<endl;
//
//    cout<<"ns "<<ns<<endl;
//    cout<<"ni "<<ni<<endl;
//    cout<<"ps "<<ps<<endl;

    double pf_ps = ps*vp*sp;
    double pf_nv = NV*vp*sp;
    double pf_ns = ns*vn*sp;
    double pf_nc = NC*vn*sp;
    double pf_Ng = NVg*vp*sp;
    double pf_ng = ng*vp*sp;

    double delET = E1-EV;

    double delE_channel_tun = E1 + (xvecbytox + 0.5/z_lim)*field*Tox;
    //cout<<"del_channel "<<delE_channel_tun<<endl;
    //Check below

    double epsTV = E2v - delE_channel_tun;   // Need to check sign

    //cout<<"epsTC "<<epsTC<<endl;
    //cout<<"epsTV "<<epsTV<<endl;

    //Check below

    // Channel Valence Band
    double eps12v = SHW/(pow((1 + R),2)) + R*epsTV/(1+R) + R*(pow(static_cast<double>(epsTV),2))/(4*SHW);
    double eps21v = eps12v - epsTV;

    k12v = eps12v > 0 ? pf_ps*exp(-beta*eps12v) : pf_ps;
    k21v = eps21v - epsTV > 0 ? pf_nv*exp(-beta*(eps21v - epsTV)) : pf_nv;

    if(isCB)
    {
        // Channel conduction Band
        double epsTC = epsTV - E2v + E2c;

        double eps12c = SHW/(pow((1 + R),2)) + R*epsTC/(1+R) + R*(pow(static_cast<double>(epsTC),2))/(4*SHW);
        double eps21c = eps12c - epsTC;

        k12c = eps12c > 0 ? pf_nc*exp(-beta*eps12c) : pf_nc;
        k21c = eps21c - epsTC > 0 ? pf_ns*exp(-beta*(eps21c - epsTC)) : pf_ns;
    }

    else
    {
        k12c = 0;
        k21c = 0;
    }

    if(isGate)
    {
        // Gate
        /******************************/
        //double delE_gate_tun = delET - (1 - xvecbytox - 0.5/z_lim)*field*Tox;
        double delE_gate_tun = E1 - (1 - xvecbytox - 0.5/z_lim)*field*Tox;
        //cout<<"del_gate "<<delE_gate_tun<<endl;
        /******************************/

        /******************************/
        //double epsTG = 0.5 - delE_gate_tun;
        double epsTG = Vg + vfb - delE_gate_tun;
        //cout<<"epsTG "<<epsTG<<endl;
        /******************************/

        double eps12g = SHW/(pow((1 + R),2)) + R*epsTG/(1+R) + R*(pow(static_cast<double>(epsTG),2))/(4*SHW);
        double eps21g = eps12g - epsTG;
        //cout<<"eps12g "<<eps12g<<endl;
        //cout<<"eps21g "<<eps21g<<endl;

        k12g = eps12g > 0 ? pf_Ng*exp(-beta*eps12g) : pf_Ng;
        k21g = eps21g + delE_gate_tun > 0 ? pf_ng*exp(-beta*(eps21g + delE_gate_tun)) : pf_ng;
    }
    else
    {
        k12g = 0;
        k21g = 0;
    }


    // Final rates
    k12 = k12v + k12c + k12g;
    k21 = k21v + k21c + k21g;

}


int main()
{
	srand (9981);
	//cout<<rand()<<endl;
	double sw_time_1 [] = {1e-6, 1e3};
    bool flagger;
    int cnt, index;
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


    // Generating filenames
    ofstream Trap_rec_file, Trap_str_file, Str_data, Rec_data, Vth_str_file, Vth_rec_file;
    stringstream file_name;
    string s;
    file_name<<"E1_"<<E1m<<"_shw"<<SHWm<<"_R_"<<Rm<<"_"<<Temp<<"C"<<".csv";

    Trap_str_file.open(("str_"+file_name.str()).c_str(), ios::out | ios::trunc);
    Trap_rec_file.open(("rec_"+file_name.str()).c_str(), ios::out | ios::trunc);
    Str_data.open ("Str_data.txt", ios::out | ios::trunc);

    if(isVthneeded)
    {
        Vth_str_file.open(("Vstr_" + file_name.str()).c_str(), ios::out | ios::trunc);
        Vth_rec_file.open(("Vrec_" + file_name.str()).c_str(), ios::out | ios::trunc);
    }



/*******************************************************************************
 CDF for all Gaussian distributions

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



/************************************************************************/
    //Create time vector with pt_dec points per decade spaced//
    double t_step;
    t_prev = sw_time_1[0];
    int cn = 0;

    //    **LINEAR**
    if(lin1_log0)
    {
        t_step = t_prev*10.0/pt_dec;
        while (t_prev < sw_time_1[1] - 0.001)
        {
            for(int j = 0; j < int(pt_dec*0.9); j++)    // Only 9/10 needed since starts with 1/10 not zero
            {
                Trap_str_file<<t_prev<<",";
                Trap_rec_file<<t_prev<<",";
                time_arr.push_back(t_prev);

                if(isVthneeded)
                {
                    // Vth files
                    Vth_str_file<<t_prev<<",";
                    Vth_rec_file<<t_prev<<",";
                }

                t_prev = t_prev + t_step;
            }
            t_step = 10*t_step;
        }
    }


    //    **LOGARITHM**
    else
    {
        t_step = pow(10, 1.0/pt_dec);
        while (t_prev < sw_time_1[1] - 0.001)
        {

            for(int j = 0; j < int(pt_dec*0.9); j++)    // Only 9/10 needed since starts with 1/10 not zero
            {
                Trap_str_file<<t_prev<<",";
                Trap_rec_file<<t_prev<<",";
                time_arr.push_back(t_prev);

                // Vth files
                Vth_str_file<<t_prev<<",";
                Vth_rec_file<<t_prev<<",";

                t_prev = t_prev*t_step;
            }
        }
    }

    Trap_str_file<<t_prev<<endl;
    Trap_rec_file<<t_prev<<endl;
    time_arr.push_back(t_prev);

    cn = time_arr.size();

    if(isVthneeded)
    {
        // Vth files
        Vth_str_file<<t_prev<<endl;
        Vth_rec_file<<t_prev<<endl;
        avg_vth_str.assign(cn,0);
        avg_vth_rec.assign(cn,0);
    }

    // Pre-assign cn+1 sizes: equal to size of time vector generated above
    avg_traps_str.assign(cn,0);
    avg_traps_rec.assign(cn,0);


// Initiialization of surface potential and electric fields for stress and recovery
    ox_field(Vgs_dc);
    double psi_str = psi, FSiO2_str = FSiO2;
    const double psi_str0 = psi;
    const double FSiO2_str0 = FSiO2;

    ox_field(Vgr_dc);
    double psi_rec = psi, FSiO2_rec = FSiO2;
    const double psi_rec0 = psi;
    const double FSiO2_rec0 = FSiO2;



    // SIMULATION LOOP
    for(int ii = 0; ii < sim_size; ii++)
    {

        //  Distribute defect sites spatially in the bulk
        for(int d1= 0; d1<num_defects; d1++){
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
                state_1.push_back(d1);
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
            Vg = Vgs_dc;
            counter = 0;
            flagger = 0;
            cn = 0;
            cnt = 0;

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

                rate((double) (int2threeD[j].z)/z_lim, psi_str, FSiO2_str);

                k12_arr[j] = k12;
                k21_arr[j] = k21;

            }
            // Create the cumulative rXn propensity array
            partial_sum(k12_arr.begin(), k12_arr.end(), cum_12.begin());

            alpha_total[0] = cum_12.back();
            cout<<alpha_total[0]<<endl;
//            for(auto iit = k1p2_arr.begin(); iit != k1p2_arr.end(); iit++)
//            {
//                cout<<*iit<<endl;
//            }
//            partial_sum(k1p2_arr.begin(), k1p2_arr.end(), cum_1p2.begin());
//            alpha_total[5] = cum_1p2.back();

            partial_sum(alpha_total, alpha_total + 2, cum_total);



            while (t < sw_time_1[1] && !flagger){

                r1 = randnu();
                while (r1 == 0){
                    r1 = randnu();
                }

                alpha0 =  cum_total[1];

                tau = (1/alpha0)*log(1/r1);
                if (tau < 0){
                    flagger = 1;
                    Trap_str_file << "Error"<<endl;
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

                        Trap_str_file<<emiprev<<",";
                        avg_traps_str[cnt] += emiprev;

                        if(isVthneeded)
                        {
                            Vth_str_file<<delVt_prev<<",";
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

            // One last entry is not filled in the while loop: cnt = size-1
            Trap_str_file<<emi<<endl;
            avg_traps_str[cnt] += emi;

            if(isVthneeded)
            {
                Vth_str_file<<delVot<<endl;
                avg_vth_str[cnt] += delVot;
            }

            cnt = 0;
            end_of_str_traps = emi;


            // Note down stats in Str_data file

            Str_data << emi <<" out of "<< num_defects<<" defects broken in bulk"<<endl;

            Str_data << "State 2: "<<state_2.size()<<endl;
            for (vector<int>::iterator it = state_2.begin(); it != state_2.end(); it++)
            {
                Str_data<<int2threeD[*it].x<<" "<<int2threeD[*it].y<<" "<<int2threeD[*it].z<<" : ";
                Str_data<<E1_store[*it]<<" "<<SHW_store[*it]<<" "<<R_store[*it]<<endl;
            }
            Str_data<<endl;

            Str_data << "State 1: "<<state_1.size()<<endl;
            for (vector<int>::iterator it = state_1.begin(); it != state_1.end(); it++)
            {
                Str_data<<int2threeD[*it].x<<" "<<int2threeD[*it].y<<" "<<int2threeD[*it].z<<" : ";
                Str_data<<E1_store[*it]<<" "<<SHW_store[*it]<<" "<<R_store[*it]<<endl;
            }
            Str_data<<endl;



        /*********************************************************************************/
        /****************************     RECOVERY     ***********************************/
        /*********************************************************************************/

            t = sw_time_1[0];
            Vg = Vgr_dc;
            cn = 0;
            cnt = 0;

            for(int j = 0; j < state_1.size(); j++)
            {
                index = state_1[j];

                E1 = E1_store[index];
                E2v = E2v_store[index];
                E2c = E2c_store[index];
                SHW = SHW_store[index];
                R = R_store[index];

                // Uodate forward and backwrd rates at possible defect sites
                rate((double) (int2threeD[index].z)/z_lim, psi_rec, FSiO2_rec);

                k12_arr[index] = k12;
                k21_arr[index] = k21;


                if(j == 0)
                {
                    cum_12[j] = k12;
                }
                else
                {
                    cum_12[j] = cum_12[j-1] + k12;

                }
            }
            alpha_total[0] = cum_12.back();


            for(int j = 0; j < state_2.size(); j++)
            {
                index = state_2[j];

                E1 = E1_store[index];
                E2v = E2v_store[index];
                E2c = E2c_store[index];
                SHW = SHW_store[index];
                R = R_store[index];

                // Uodate forward and backwrd rates at possible defect sites
                rate((double) (int2threeD[index].z)/z_lim, psi_rec, FSiO2_rec);

                k12_arr[index] = k12;
                k21_arr[index] = k21;

                if(j == 0)
                {
                    cum_21[j] = k21;
                }
                else
                {
                    cum_21[j] = cum_21[j-1] + k21;
                }
            }
            alpha_total[1] = cum_21.back();


            partial_sum(alpha_total, alpha_total + 2, cum_total);
//            for(int j = 0; j < 8; j++)
//            {
//                cout<<cum_total[j]<<endl;
//            }


           while (t < sw_time_1[1] && flagger == 0){

                r1 = randnu();
                while (r1 ==0){
                    r1 = randnu();
                }

                alpha0 =  cum_total[1];

                tau = (1/alpha0)*log(1/r1);
                if (tau < 0){
                    flagger = 1;
                    Trap_rec_file << "Error"<<endl;
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
                    Trap_rec_file<<emiprev<<",";
                    avg_traps_rec[cnt] += emiprev;

                    if(isVthneeded)
                    {
                        Vth_rec_file<<delVt_prev<<",";
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
            Trap_rec_file<<emi<<endl;
            avg_traps_rec[cnt] += emi;


            if(isVthneeded)
            {
                Vth_rec_file<<delVot<<endl;
                avg_vth_rec[cnt] += delVot;
            }




            cum_12.assign(num_defects,0);
            fill(cum_total, cum_total + 2, 0);
            fill(alpha_total, alpha_total + 2, 0);

            cum_21.clear();
            state_1.clear();
            state_2.clear();

            threeD2int.clear();
            int2threeD.clear();

            psi_str = psi_str0;
            psi_rec = psi_rec0;
            FSiO2_str = FSiO2_str0;
            FSiO2_rec = FSiO2_rec0;


    }
    //Printing out average values at end of file
    for(int h = 0; h < avg_traps_rec.size(); h++)
    {
        Trap_rec_file<<avg_traps_rec[h]*1.0/sim_size<<",";
        Trap_str_file<<avg_traps_str[h]*1.0/sim_size<<",";

        if(isVthneeded)
        {
            Vth_str_file<<avg_vth_str[h]/sim_size<<",";
            Vth_rec_file<<avg_vth_rec[h]/sim_size<<",";
        }

    }

    Trap_rec_file.close();
    //Rec_data.close();
    Trap_str_file.close();
    Str_data.close();

    if(isVthneeded)
    {
        Vth_rec_file.close();
        Vth_str_file.close();
    }


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

	int l, m;
	double change;
	vector<double>::iterator itcum;
	//
	if (r2 < cum_total[0]) {
            l =(int) (lower_bound(cum_12.begin(), cum_12.end(), r2) - cum_12.begin());
            m = state_1[l]; // Trap id

            // Update cumulative propensities for new state
            if(!state_2.empty())
            {
                cum_21.push_back(cum_21.back() + k21_arr[m]);
            }
            else
            {
                cum_21.push_back(k21_arr[m]);
            }
            alpha_total[1] = cum_21.back();
            // Add trap to state 2
            state_2.push_back(m);
            emi++;

            // Remove trap from state 1
            state_1.erase(state_1.begin() + l);
            // Subtract from cum prop
            if(!state_1.empty())
            {
                itcum = cum_12.erase(cum_12.begin() + l);
                for_each(itcum, cum_12.end(), update_prop(-k12_arr[m]));
                alpha_total[0] = cum_12.back();
            }

            else
            {
                cum_12.clear();
                alpha_total[0] = 0;
            }

            // Update total propensities for all 8 rXns
            partial_sum(alpha_total, alpha_total + 8, cum_total);

            if(isVthneeded)
            {
                delVot += 1000*q*Tox*(1 - (int2threeD[m].z + 0.5)/(1.0*z_lim))/(dev_area*ESiO2*E0);
            }

	}

	//
	else if (r2 < cum_total[1]) {


            l =(int) (lower_bound(cum_21.begin(), cum_21.end(), r2 - cum_total[0]) - cum_21.begin());
            m = state_2[l];    //Trap id


            // Update cumulative propensities for new state
            if(!state_1.empty())
            {
                cum_12.push_back(cum_12.back() + k12_arr[m]);
            }
            else
            {
                cum_12.push_back(k12_arr[m]);
            }
            alpha_total[0] = cum_12.back();
            // Add trap to state 1
            state_1.push_back(m);

            // Remove trap from state 1'
            state_2.erase(state_2.begin() + l);
            emi--;
            // Subtract from cum prop
            if(!state_2.empty())
            {
                itcum = cum_21.erase(cum_21.begin() + l);
                for_each(itcum, cum_21.end(), update_prop(-k21_arr[m]));
                alpha_total[1] = cum_21.back();
            }
            else
            {
                cum_21.clear();
                alpha_total[1] = 0;
            }

            // Update total propensities for all 8 rXns
            partial_sum(alpha_total, alpha_total + 2, cum_total);

            if(isVthneeded)
            {
                if(emi==0) delVot = 0;
                else delVot -= 1000*q*Tox*(1 - (int2threeD[m].z + 0.5)/(1.0*z_lim))/(dev_area*ESiO2*E0);
                //cout<<delVot<<endl;
            }




	}

}
