 


#ifndef _FUNCTIONS_H_INCLUDED_ //include guard
#define _FUNCTIONS_H_INCLUDED_

#include<iostream>
#include<vector>
#include<algorithm>
#include <fstream>
#include <string>
#include<sstream>
#include <ctime>
#include <cstdlib>
#include<time.h>
#include <cmath>
#include <stdio.h>
#include <numeric>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>




#define path1 "/xx/scripts/" /*folder for the codes*/

#define path2 "/xx/data/fmd_2001/"      /*FMD data*/
// #define path2b "/home/sl451/Dropbox/HWU/phylo_seir_fmd_inference/data/fmd_2001_partial_randscpt_2/set 1_10/"       /*FMD data*/
#define path4 "xx/mcmc samples/"      /*output folder*/

using namespace std;
//using namespace Rcpp;


typedef  vector < vector<int> > vec2int;
typedef  vector < vector<double> > vec2d;

typedef boost::mt19937 base_generator_type;

const double t_range = 7.0; // for update of t_i; t_i would range between [t_o - t_range, t_o + range] (see paper SI)
const double t_back =100.0; // assume a max length for the latent period 

const double rate_exp_prior = 0.001;// the rate of exp to be used as vague prior

const int ind_n_base_part =0;// =1 if the seq data is partial
const int n_base_part =1000; // the partial length used if ind_n_base_part =1


/* lower/upper bounds for parameters */

const double alpha_hi = 0.1; 
const double beta_hi = 50.0;
const double mu_lat_hi =50.0;

const double var_lat_lo = 0.1;
const double var_lat_hi = 30.0;

const double c_hi =100.0;
const double d_hi = 100.0;
const double k_1_hi = 20.0;
const double k_2_hi = 100.0;

const double mu_1_hi = 0.1;
const double mu_2_hi = 0.1;

const double p_ber_hi = 1.0;

/*
const double  mu_lat_p_1= 4.7; // first parameter for the prior dist of mu_lat
const double  mu_lat_p_2= 0.3; // second ...

const double  var_lat_p_1= 2.16; // first parameter for the prior dist of var_lat
const double  var_lat_p_2= 1.3; // second ...

const double  beta_p_1= 5.81; // first parameter for the prior dist of beta
const double  beta_p_2= 1.7; // second ...

const double  k_1_p_1= 0.02; // first parameter for the prior dist of k_1
const double  k_1_p_2= 0.002; // second ...

const double  c_p_1= 1.84; // first parameter for the prior dist of c
const double  c_p_2= 0.1; // second ...

const double  d_p_1= 2.4; // first parameter for the prior dist of d
const double  d_p_2= 0.5; // second ...*/


struct nt_struct { // structure of the sequence data

vec2int nt; // the sequence corresponds to the  infected-or-infecting event OR sampling of a farm
vec2d t_nt ; // the time corresponds to the sequence above
vector<int> current_size; // the number of above sequences recorded (in a farm)

vector<double> t_sample; // sampling time of sequences (in a farm)

vec2int infecting_list; // each elements contains the list of subjects being infected by a particular farm (sorted in order of infection)
vector<int> infecting_size; // number of subjects infected by the farm (if sampeled, it would be equal to current_size -2 as also excluding the first record corresponds to own infection)


};


struct para_key {
double alpha,beta,a,b,mu_lat,var_lat,c,d,k_1,k_2; //key model parameters (a,b are redundant)
double stb_1, stb_2;//reduandant, ignore them
double mu_1, mu_2; // mutation rates
double p_ber; // varation parameter p (for the master sequence)
};


//----------------


struct para_aux{ // other model quantities to be fixed/known
int n, n_seq, n_base; //n, population size (number of farms/sites); n_seq, max number of sequences for a farm (to reserve the capacity of the vector storing sequence data); n_base, number of bases for a sequence
string kernel_type; // type of spatial kernel used e.g exponential/power-law
double dimen_x,dimen_y; // ignore for now (was dimension of space for simulation)
double t_max,unassigned_time; // t_max is the upper limit of observation period; unassigned_time is an abritary extreme value e.g. 9e+10 , to indicate that an event does not happen e.g. no infection
int seed; // random seed used, may be ignored
};
//-----------------



struct epi_struct { // structure of the epidemiological data
 vector<int> k; // k=1,....n indexes for farms/sites
 vector<double> q, t_e, t_i,t_r; // ignore q for now
 vector<int> status; //ignore for now
 vector<double> coor_x,coor_y; // coordinates (lon and lat)
 
 vector<int> gp_stb; // redundant, ignore for now
 vector<double> stb; // redundant, ignore for now

 vector<int> infected_source; // sources of infections 

};

//------------------


struct lh_SQUARE { // the contribution to likelihood function from each individual is subdivided into f_U, f_E, f_I, f_EnI, f_R, f_InR, log_f_S, log_f_Snull (this is used so that we can only update necessary parts of the likelihood when doing MCMC updates, without re-calculating the whole likelihood); this trick may only work when you are working on a general SEIR spatial model (see the paper)
vector<long double> f_U, q_T, kt_sum_U; // f_U (and its "composition"), likelihood contribution from susceptibles
vector<long double> f_E, g_E,k_sum_E, h_E,q_E,kt_sum_E; // f_E(...), likelihood contribution of infected/exposed
vector<long double> f_I, f_EnI; // f_I, .. of those once become infectious; f_EnI, ...of those infected but not yet become infectious
vector<long double> f_R, f_InR; //f_R, ... of recovered; f_InR, .. infectious but not yet recovered
vector<long double> log_f_S, log_f_Snull;// log_f_S, log-likelihood  of sequence data; log-likelihood of the first sequence (i.e. background sequence, a variant from the master sequence)
};




//----------------------------

inline long double func_kernel (double, double, double, double, double, double, const string&); // function prototype for calculating kernel distance

double log_lh_func (lh_SQUARE, int);

inline double lh_snull(const vector<int>&, const vector<int>&, const double& , const int&);
inline double lh_snull_base(const int& , const int& , const double&); // compute the log pr a base for background

//--------------------------

inline void sample_snull (const vector<int>&, vector<int>&, const double&, const int&, gsl_rng* );//sample a seq for background

//----------------------------

inline void seq_propose_tri (vector<int>& ,  double&, const vector<int>& , const vector<int>& , const vector<int>& , const double& , const double& , const double& , const double& ,  const double& , const double& , const  int& , gsl_rng * );

inline void seq_backward_pr_tri (const vector<int>& ,  double& , const vector<int>& , const vector<int>& , const vector<int>& , const double& , const double& , const double& , const double& , const double& , const double& , const int& );

//--------------------------

inline void seq_propose_cond(vector<int>& , double&, const vector<int>&, const vector<int>& , const double& , const double& , const double&,  const double& , const double& ,int, gsl_rng*); // proposing a new sequence condtional on a future sequence at the direction of time change

inline void seq_propose_uncond(vector<int>& ,  double&, const vector<int>& , const double& , const double& , const double& ,  const double& , const double& , const double& , const double & , const double& , const double& , const int& , const int& , const int& , gsl_rng * );  // proposing a new sequence WITHOUT condtional on a future sequence at the direction of time change

//---------------------------

inline void seq_backward_pr_cond(const vector<int>& , double& , const vector<int>& , const vector<int>&, const double& , const double& ,  const double& ,   const double& , const double& ,int );

inline void seq_backward_pr_uncond(const vector<int>& ,  double&, const vector<int>& , const double&, const double& , const double& ,  const double& , const double& , int);
//--------------------------
double log_lh_base (int&, int&, double, double  , double, double);

double log_lh_seq (vector<int>&, vector<int>& , double , double  , double , double , int );

//--------------------------


class FUNC {

private:

double alpha_Clh;
double beta_Clh;
double a_Clh;
double b_Clh;
double mu_lat_Clh;
double var_lat_Clh;
double c_Clh;
double d_Clh;
double k_1_Clh;
double k_2_Clh;

double mu_1_Clh;
double mu_2_Clh;

double p_ber_Clh;

int n_Clh;
string kernel_type_Clh;
double dimen_x_Clh,dimen_y_Clh;
double t_max_Clh,unassigned_time_Clh;
int seed_Clh;

int n_seq_Clh, n_base_Clh;

vector< vector<double> > coordinate_Clh;

vector<int> xi_U_Clh, xi_E_Clh, xi_E_minus_Clh, xi_I_Clh, xi_R_Clh, xi_EnI_Clh, xi_InR_Clh;
vector<double> t_e_Clh, t_i_Clh, t_r_Clh;

vector<int> index_Clh;

// vector<int> con_seq_Clh;

vector<double> stb_Clh;

vector <int> infected_source_Clh;

public:

void set_para (para_key para_current_arg, para_aux para_other_arg, vector< vector<double> > coordinate_arg, vector<int> xi_U_arg, vector<int> xi_E_arg, vector<int> xi_E_minus_arg ,vector<int> xi_I_arg, vector<int> xi_R_arg, vector<int> xi_EnI_arg, vector<int> xi_InR_arg, vector<double> t_e_arg, vector<double> t_i_arg, vector<double> t_r_arg, vector<double> stb_arg, vector<int> index_arg, vector<int> infected_source_arg) {

alpha_Clh=para_current_arg.alpha;
beta_Clh=para_current_arg.beta;
a_Clh=para_current_arg.a;
b_Clh=para_current_arg.b;
mu_lat_Clh=para_current_arg.mu_lat;
var_lat_Clh=para_current_arg.var_lat;
c_Clh=para_current_arg.c;
d_Clh=para_current_arg.d;
k_1_Clh=para_current_arg.k_1;
k_2_Clh=para_current_arg.k_2;

mu_1_Clh=para_current_arg.mu_1;
mu_2_Clh=para_current_arg.mu_2;

p_ber_Clh = para_current_arg.p_ber;

n_Clh = para_other_arg.n;
kernel_type_Clh = para_other_arg.kernel_type;
dimen_x_Clh = para_other_arg.dimen_x;
dimen_y_Clh = para_other_arg.dimen_y;
t_max_Clh = para_other_arg.t_max;
unassigned_time_Clh = para_other_arg.unassigned_time;
seed_Clh = para_other_arg.seed;

n_base_Clh = para_other_arg.n_base;
n_seq_Clh = para_other_arg.n_seq;

//con_seq_Clh = con_seq;

coordinate_Clh = coordinate_arg;

xi_U_Clh = xi_U_arg;
xi_E_Clh = xi_E_arg;
xi_E_minus_Clh = xi_E_minus_arg;
xi_I_Clh = xi_I_arg;
xi_R_Clh = xi_R_arg;
xi_EnI_Clh = xi_EnI_arg;
xi_InR_Clh = xi_InR_arg;


t_e_Clh = t_e_arg;
t_i_Clh = t_i_arg;
t_r_Clh = t_r_arg;

index_Clh = index_arg;

stb_Clh = stb_arg;

infected_source_Clh = infected_source_arg;

}


void initialize_kernel_mat (vector< vector<double> >&, vector<double>&); // function prototype for initializing kernel distance  

void initialize_delta_mat (vector< vector<double> >&); // function prototype for initializing length of exposure time

void initialize_lh_square (lh_SQUARE&, vector< vector<double> >, vector< vector<double> >, vector<double>&, nt_struct&, vector<int>&);

};

//-----------------------------

class mcmc_UPDATE {

private:

// double alpha_CUPDATE;
// double beta_CUPDATE;
// double a_CUPDATE;
// double b_CUPDATE;
// double c_CUPDATE;
// double d_CUPDATE;
// double k_1_CUPDATE;
// double k_2_CUPDATE;

int n_CUPDATE;
string kernel_type_CUPDATE;
double dimen_x_CUPDATE,dimen_y_CUPDATE;
double t_max_CUPDATE,unassigned_time_CUPDATE;
int seed_CUPDATE;

int n_seq_CUPDATE, n_base_CUPDATE;

vector< vector<double> > coordinate_CUPDATE;

// vector<int> xi_U_CUPDATE, xi_E_CUPDATE, xi_E_minus_CUPDATE, xi_I_CUPDATE, xi_R_CUPDATE;
// vector<double> t_e_CUPDATE, t_i_CUPDATE, t_r_CUPDATE;

vector<int> index_CUPDATE;

// vector<int> con_seq_CUPDATE;

public:

void set_para (para_aux para_other_arg, vector< vector<double> > coordinate_arg) {

// alpha_CUPDATE=para_current_arg.alpha;
// beta_CUPDATE=para_current_arg.beta;
// a_CUPDATE=para_current_arg.a;
// b_CUPDATE=para_current_arg.b;
// c_CUPDATE=para_current_arg.c;
// d_CUPDATE=para_current_arg.d;
// k_1_CUPDATE=para_current_arg.k_1;
// k_2_CUPDATE=para_current_arg.k_2;

n_CUPDATE = para_other_arg.n;
kernel_type_CUPDATE = para_other_arg.kernel_type;
dimen_x_CUPDATE = para_other_arg.dimen_x;
dimen_y_CUPDATE = para_other_arg.dimen_y;
t_max_CUPDATE = para_other_arg.t_max;
unassigned_time_CUPDATE = para_other_arg.unassigned_time;
seed_CUPDATE = para_other_arg.seed;

coordinate_CUPDATE = coordinate_arg;

// con_seq_CUPDATE = con_seq;

n_base_CUPDATE = para_other_arg.n_base;
n_seq_CUPDATE = para_other_arg.n_seq;


/*
xi_U_CUPDATE = xi_U_arg;
xi_E_CUPDATE = xi_E_arg;
xi_E_minus_CUPDATE = xi_E_minus_arg;
xi_I_CUPDATE = xi_I_arg;
xi_R_CUPDATE = xi_R_arg;

t_e_CUPDATE = t_e_arg;
t_i_CUPDATE = t_i_arg;
t_r_CUPDATE = t_r_arg;

index_CUPDATE = index_arg;*/
}

public:

void alpha_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector <double>&, const vector<int>&, const vector<double>& , const vector<int>&, para_key&, int,gsl_rng* &);
void beta_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector <double>&, const vector<int>&, const vector<double>& , const vector<int>&,para_key&, int,gsl_rng* &);

void c_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, const vector<int>&, para_key& ,int, gsl_rng* &);
void d_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, const vector<int>&, para_key& ,int,gsl_rng* &);
void k_1_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , const vector<int>& ,vector<double>&, int,gsl_rng* & );
void k_2_update(lh_SQUARE& , double& , vector< vector<double> >&, const  vector< vector<double> >&, const vector<int>&, const vector<int>&,const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&, const  vector<int>&, para_key&,   const vector<double>& , const vector<int>& , vector<double>&, int,gsl_rng* & );

void mu_1_update(lh_SQUARE& , double& ,  const vector<int>& , para_key& , nt_struct& , int iter,gsl_rng* &);
void mu_2_update(lh_SQUARE& , double& ,  const vector<int>& , para_key& , nt_struct& , int iter,gsl_rng* &);

void p_ber_update (lh_SQUARE& , double& ,  const vector<int>& , para_key& , nt_struct& , vector<int>& ,vector<int>&, int iter,gsl_rng* &);


void mu_joint_update(lh_SQUARE& , double& ,  const vector<int>& , para_key& , nt_struct& , int iter);

void mu_exp_update(lh_SQUARE&, double&, const vector<int>&,  const vector<int>&,  const vector<double>&, const vector<double>&, const vector<int>&, para_key&, int,gsl_rng* &);


void mu_lat_update(lh_SQUARE&, double&, const vector<int>&,  const vector<int>&,  const vector<double>&, const vector<double>&, const vector<int>&, para_key&, int,gsl_rng* &);
void var_lat_update(lh_SQUARE&, double&, const vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, const vector<int>&, para_key&, int,gsl_rng* &);

void stb_1_update(lh_SQUARE& , double& , const vector<int>& , const vector<int>& , const vector<double>& , const vector<int>& , vector<double>& , para_key& , const vector<int>& , int );
void stb_2_update(lh_SQUARE& , double& , const vector<int>& , const vector<int>& , const vector<double>& , const vector<int>& ,  vector<double>& , para_key& , const vector<int>& , int );


void t_e_seq(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>& , vec2int&, vec2d&,  vec2int&, const vector<int>&, vector<int>&,vector<int>&, int&, int,gsl_rng* &); // jointly update the infection time and sequences

void index_first_seq(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>& , vec2int&, vec2d&,  vec2int&, const vector<int>&, vector<int>&, vector<int>&,int&, int,gsl_rng* &);


void t_i_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, const vector<int>&, const vector<int>&, const vector<int>&, const vector<int>&, const vector<int>&, const vector<int>&, const vector<int>&,  const vector<double>&, vector<double>&, const vector<double>&, const vector<double>&, const vector<int>&, const para_key&, const vector<double>&, const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>&, vec2int&, vec2d&,  vec2int&, const vector<int>&, int);

void t_e_seq_add_del(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, vector<int>&, const vector<double>&, vector<int>& , vec2int&, vec2d&, vec2int&, vector<int>&, vector<int>&,int); // add or del an infection without sample and t_i

void con_seq_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>& , vec2int&, vec2d&, vector<int>&, vector<int>&, int, gsl_rng*&); //

void seq_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>& , vec2int&, vec2d&, vector<int>&, vector<int>&, const int&, int, gsl_rng*&); //

void seq_n_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, const vector<int>& , vec2int&, vec2d&, vector<int>&, const int&, gsl_rng*&); // update the whole seq

void source_t_e_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, vector<int>&, const vector<double>&, vector<int>& , vec2int&, vec2d&, vec2int&, vector<int>&,vector<int>&, int&, vector<int>& , vector<int>&,int,gsl_rng* &); // update the infecting source of an infection

void source_t_e_update_V2(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, vector<int>&, const vector<double>&, vector<int>& , vec2int&, vec2d&, vec2int&, vector<int>&,vector<int>&, int&, vector<int>& ,int,gsl_rng* &); // do not update t_e compare to source_t_e_update 


void source_t_e_update_tri(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, vector<int>&, const vector<double>&, vector<int>& , vec2int&, vec2d&, vec2int&, vector<int>&,vector<int>&, int&, vector<int>& ,int); // update the infecting source of an infection. consider the "triangle" in proposing a new source and respective sequence


//

void source_update(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, vector<int>&, const vector<double>&, vector<int>& , vec2int&, vec2d&, vec2int&, vector<int>&,vector<int>&, int&, vector<int>& ,int); // update the infecting source of an infection


void t_e_replct(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& , const vector<double>&, const vector<int>&, const vector<double>&, int);

void t_e_add_del(lh_SQUARE&, double&, const vector< vector<double> >&, vector< vector<double> >&, vector<int>&, vector<int>&, vector<int>& , const vector<int>&, vector<int>&, const vector<double>&, const vector<double>&, vector<double>& , vector<int>&, const para_key&,  const vector<double>& ,const vector<double>&, const vector<int>&,int);

//

 };

//-----------------------------
void count_type_all(nt_struct &, vector<int>&, int&, int& , int&, int& ); // count number of unchanged, transition, transversion (whole dataset)

void count_type_seq (vector<int>& , vector<int> , int&, int&, int& , int& );

//---------------------------
void residual_gnl(const vector< vector<double> >&, const vector<int>&, const vector<int>&, const vector<int>&,  const vector<double>&, const vector<double>&, const vector<double>&, const para_key&, const para_aux&,  const vector<double>& , const vector<double>&, int);

void residual_kernel_path(const vector< vector<double> >&, const vector<int>&, const vector<int>&,  const vector<double>&, const vector<double>&, const vector<double>&, const para_key&, const para_aux&, const vector<double>& , const vector<double>&,  const vector<int>&, int); // do not need to impute link given transmission path

void residual_lat( const vector<int>&, const vector<int>&,  const vector<int>& , const vector<double>&, const vector<double>&, const para_key&, const para_aux&,  int);

//-----------------------------

void residual_nucle (const para_key& , const para_aux& , vector<int>& ,const vector<int>& , const vec2int& , const vec2d&, int);

void residual_sub_nucle (vector<double>& , const vector<int>&, const vector<int>& , const double& , const double& , double , double , const int&, base_generator_type&);

//-----------------------------

void IO_simpara(para_key&, para_aux&);
void IO_simdata(para_key, para_aux ,vector < vector<double> >& , epi_struct& , nt_struct&, vector<int>&,vector<int>&);
//void lh_plot (para_key, para_aux , vector< vector<double> >, vector<int>, vector<int> , vector<int> , vector<int> ,vector<int>, vector<int> ,vector<int>, vector<double> , vector<double> , vector<double>, vector<int>);

void initialize_mcmc(para_key& , para_aux&,  vector<int>&, vector<int>& , vector<int>&, vector<int>& , vector<int>&,  vector<int>& , vector<int>&, vector<int>&,  vector<double>& , vector<double>& , vector<double>&, vector<int>& ,  vector<double>& , vector<int>& , vector<int>& , vector < vector<double> >&, vector <double>&, vec2int&, vector<double>&, nt_struct&, vector<int>&);

void seq_initialize_pair(nt_struct&, int , int, double, int, para_key&); // update sequences of infectious-infected pair

void delete_seq_samples(para_key& , para_aux&,  vector<int>&, vector<int>& , vector<int>&, vector<int>& , vector<int>&,  vector<int>& , vector<int>&, vector<int>&,  vector<double>& , vector<double>& , vector<double>&, vector<int>& ,  vector<double>& , vector<int>& , vector<int>& , vector < vector<double> >&, vector <double>&, vec2int&, vector<double>&, nt_struct& );


#endif
