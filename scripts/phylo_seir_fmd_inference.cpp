 
#include "phylo_seir_fmd_inference_header.h"


/*
header file xx.h defines some global variables and function prototypes
xx_functions.cpp define functions in details
xx_IO.cpp define functions for inputting data */

/* key data structures (see header file for details): nt_struct (for holding sequence data); epi_struct (holding epi data); lh_SQUARE (holding likelihood which is sudivided into many blocks); para_key, key model parameters; para_aux, other unchanged model parameters e.g. type of spatial kernel

/*Note: likelihood which is sudivided into many blocks so that we can perform partial updates of the likelihood to save some computational time; but this trick may work only for a general SEIR framework used in the paper. You may have to modifiy the code in other cases */

/*key classes (see header file for details): I defined a class mcmc_UPDATE so that I can use many variables, e.g. n_base (look for n_base_CUPDATE), as private and do not need to pass them in the function arguments; see header file, mcmc_UPDATE; you of course can pass those variables in function arguments and do not use the class mcmc_UPDATE*/
/* Similarly for class Func */

/*key functions: t_e_seq, update infection time and seq jointly; source_t_e_update, update infection time, source of infection and sequence jointly (see the paper)

/* explanations for some function arguments:*/
		
// kernel_mat_current_arg : a n by n symmetric matrix (n=number of farms/sites) that stores the "spatial risk" between an element in a row and an element in a column i.e. the K(d_ij); unchanged in this update
// delta_mat_current_arg: a n by n matrix that stores the exposure duration of an element in row to an element in column (values can be 0)

// xi_U_arg: farms (with indexes k=1,2,...n; see definition of epi_struct below) remian uninfected throughout throughout the whole observation period
// xi_E_arg: farms once infected
// xi_E_minus_arg: farms once infected, excluding the first infected farm
// xi_I_arg: farms once become infectious
// xi_EnI_arg: farms infected but not yet infectious
// t_i_arg: the times of becoming infectious (farms have not once become infectious would have an entry with an arbitrary extreme value e.g. 9e+10)
// t_e_arg: exposure/infection times
// index_arg: farms with earliest infections (may be more than one if they were infected together); assumed known in this FMD analysiss
// para_current_arg: current values of model parameters, see definition of para_key below
// stb_arg, ignore for now, let their vector elements be 1 
// norm_const_current_arg: normalizing constant, you can devise your own way of normalizing the spatial kernel (you can set the values to be 1 if you do not want to normalize)
// con_seq: the master sequence G_M

//subject_proposed is the farm you wanna propose change//
// r_c, a random seed generator for gsl, look for gsl_xxx in the function. You can delete it if you use other random number generators. Current way might not be very efficient.


using namespace std;

int main (int argc, char *argv[]){

//RInside R(argc, argv);


ifstream myfile_in;
ofstream myfile_out; 

para_key para_true; // redundant (was for simulation)

para_aux para_other;

epi_struct epi_final;

nt_struct nt_data;

vector<int> index;

IO_simpara (para_true, para_other);  //Importing unchanged parameters e.g. n_base, n_seq

vector < vector<double> > coordinate(para_other.n,vector<double>(2));

vector<int> con_seq, con_seq_estm;

IO_simdata(para_true, para_other, coordinate, epi_final, nt_data, index, con_seq_estm); //Importing epidemiological data, estimated master sequence (con_seq_estm) . etc


vec2int sample_data; // 2-d vector contains the sampled sequences; non-sampled premises would have unexpected values
sample_data.resize(para_other.n);



string line, field;
int line_count = 0;

myfile_in.open((string(path2)+ string("seq_fmd.txt")).c_str(),ios::in); // for importing sequence data (see below)

while (getline(myfile_in,line)) {
	stringstream ss(line);
	//field_count=0;
	
	while (getline(ss, field, ',')) {
		stringstream fs (field);
		int t;
		fs >> t;
		sample_data[line_count].push_back(t);
		
		//field_count = field_count + 1;
	}
	line_count = line_count + 1 ;
}// end while for getline

myfile_in.close();


for (int i=0;i <=(para_other.n-1);i++){

	myfile_out.open((string(path4) +string("sample_data.txt")).c_str(),ios::app);
//	if (nt_data.t_sample.at(i)!=para_other.unassigned_time){
		for (int j=0;j<=((int)sample_data.at(i).size()-1);j++){
		int rem = (j+1)%para_other.n_base;
		if ((rem!=0) | (j==0)) myfile_out << sample_data[i][j] << ",";
		if ((rem==0) & (j!=0)) myfile_out<< sample_data[i][j] << " " << endl;	
		}
//	}
	myfile_out.close();
}


/*----------------------------*/

vector<int> xi_U, xi_E, xi_E_minus, xi_I, xi_R, xi_EnI, xi_EnIS, xi_InR; // indices sets indicating the individuals stay in S OR have gone through the other classes (E OR I OR R), and individuals hve gone through E but not I (EnI) and I but not R (InR)

vector<int> xi_beta_E; // vector contains the individuals with secondary infection

xi_U.reserve(para_other.n); //dynamically updated if necessary
xi_E.reserve(para_other.n);
xi_E_minus.reserve(para_other.n);
xi_I.reserve(para_other.n);
xi_R.reserve(para_other.n);
xi_EnI.reserve(para_other.n);
xi_EnIS.reserve(para_other.n);
xi_InR.reserve(para_other.n);
xi_beta_E.reserve(para_other.n);


for (int i=0; i<=(para_other.n-1);i++){
if (epi_final.t_e.at(i)==para_other.unassigned_time) xi_U.push_back(i);
if (epi_final.t_e.at(i)!=para_other.unassigned_time) xi_E.push_back(i);
if (epi_final.t_i.at(i)!=para_other.unassigned_time) xi_I.push_back(i);
if (epi_final.t_r.at(i)!=para_other.unassigned_time) xi_R.push_back(i);
}

xi_E_minus = xi_E;
for (int i=0; i<= (int)(index.size()-1);i++){
xi_E_minus.erase(find(xi_E_minus.begin(),xi_E_minus.end(),index.at(i)));
} // E set excluding index

xi_EnI = xi_E;
for (int i=0; i<= (int)(xi_I.size()-1);i++){
xi_EnI.erase(find(xi_EnI.begin(),xi_EnI.end(),xi_I.at(i)));
} // E set excluding I

xi_EnIS = xi_EnI;
for (int i=0; i<= (int)(xi_EnI.size()-1);i++){
	if (nt_data.t_sample.at(xi_EnI.at(i))!=para_other.unassigned_time){
	xi_EnIS.erase(find(xi_EnIS.begin(),xi_EnIS.end(),xi_EnI.at(i)));
	}
} // E set excluding I and sampled


xi_InR = xi_I;
for (int i=0; i<= (int)(xi_R.size()-1);i++){
xi_InR.erase(find(xi_InR.begin(),xi_InR.end(),xi_R.at(i)));
} // I set excluding R

for (int i=0; i<= (int)(para_other.n-1);i++){
if (( epi_final.infected_source.at(i)!=9999) & ( epi_final.infected_source.at(i)!=-99)) xi_beta_E.push_back(i);
}


/*----------------------------*/

lh_SQUARE lh_square; // the object contains the information of the likelihood contribution of each individual, and it is dynamic and changed during MCMC sampling

lh_square.f_U.assign(para_other.n,1.0);
lh_square.q_T.assign(para_other.n,0.0);
lh_square.kt_sum_U.assign(para_other.n,0.0);
lh_square.f_E.assign(para_other.n,1.0);
lh_square.g_E.assign(para_other.n,1.0);
lh_square.h_E.assign(para_other.n,1.0);
lh_square.k_sum_E.assign(para_other.n,0.0);
lh_square.q_E.assign(para_other.n,0.0);
lh_square.kt_sum_E.assign(para_other.n,0.0);
lh_square.f_I.assign(para_other.n,1.0);
lh_square.f_R.assign(para_other.n,1.0);
lh_square.f_EnI.assign(para_other.n,1.0);
lh_square.f_InR.assign(para_other.n,1.0);
lh_square.log_f_S.assign(para_other.n,0.0);
lh_square.log_f_Snull.assign(para_other.n,0.0);





/*-----------------------------------------------------Start of MCMC sampling------------------------------------------*/

para_key para_current = para_true; // not valid
vector<int> xi_I_current = xi_I;
vector<int> xi_U_current = xi_U;
vector<int> xi_E_current = xi_E;
vector<int> xi_E_minus_current = xi_E_minus;
vector<int> xi_R_current= xi_R;
vector<int> xi_EnI_current = xi_EnI;
vector<int> xi_EnIS_current = xi_EnIS;
vector <int> xi_InR_current =xi_InR;
vector<double> t_e_current = epi_final.t_e;
vector<double>t_i_current =  epi_final.t_i;
vector<double>t_r_current =  epi_final.t_r;
vector<int>index_current = index;
vector<int> xi_beta_E_current = xi_beta_E;

vector<int> con_seq_current = con_seq_estm;

vector<double> stb_current =  epi_final.stb;
vector<int> gp_stb_current =  epi_final.gp_stb;

vector<int> infected_source_current =  epi_final.infected_source;

nt_struct nt_data_current;
nt_data_current.t_sample = nt_data.t_sample;


vector<double> t_onset =  epi_final.t_i ;// used as a ref point to sample t_i; assumed to be given,e.g, the onset time  (in simulation, this may be taken to be the true t_i)

//Note: struct copy is fragile: nt_data_current=nt_data wont work!!//

//----------------------------------------//
//stb_current.assign(para_other.n,1.0); // when ignoring the stb
//------------------------------------//


lh_SQUARE lh_square_current ;
vector < vector<double> > kernel_mat_current(para_other.n, vector<double>(para_other.n)), delta_mat_current(para_other.n, vector<double>(para_other.n)); // a dynamic matrix contain the "kernel distance" 
vector <double> norm_const_current(para_other.n);

//---------------------------------------------------//

initialize_mcmc(para_current, para_other, xi_I_current, xi_U_current,xi_E_current, xi_E_minus_current, xi_R_current,   xi_EnI_current, xi_EnIS_current,  xi_InR_current,  t_e_current,  t_i_current, t_r_current,  index_current,   stb_current,  gp_stb_current,  infected_source_current, kernel_mat_current, norm_const_current,  sample_data, t_onset, nt_data_current, con_seq_current); // initialze the parameters/unobserved data for mcmc



//delete_seq_samples(para_current, para_other, xi_I_current, xi_U_current,xi_E_current, xi_E_minus_current, xi_R_current,   xi_EnI_current, xi_EnIS_current,  xi_InR_current,  t_e_current,  t_i_current, t_r_current,  index_current,   stb_current,  gp_stb_current,  infected_source_current, kernel_mat_current, norm_const_current,  sample_data, t_onset, nt_data_current); // delete all samples and sampling times, keep other config 

//para_current.k_2 = 0.1; // needed when no initializtion and wanna set a reasonable value for k_2 in wrong kernel

//---------------------------------------------------//


// for (int i=0;i <=(para_other.n-1);i++){
// 
// 	ostringstream convert;
// 	convert << i;
// 
// 	myfile_out.open((string(path4)+ string("initial_subject_").c_str() + string(convert.str()) +string("_t_nt.txt")).c_str(),ios::app);
// 	for (int j=0;j<=(nt_data_current.current_size.at(i)-1);j++){
// 	myfile_out <<nt_data_current.t_nt[i].at(j)<< endl;
// 	}
// 	myfile_out.close();
// 	
// 	myfile_out.open((string(path4)+ string("initial_subject_").c_str() + string(convert.str()) +string("_nt.txt")).c_str(),ios::app);
// 	for (int j=0;j<=(nt_data_current.current_size.at(i)*para_other.n_base-1);j++){
// 	int rem = (j+1)%para_other.n_base;
// 	if ((rem!=0) | (j==0)) myfile_out << nt_data_current.nt[i][j] << ",";
// 	if ((rem==0) & (j!=0)) myfile_out<< nt_data_current.nt[i][j] << " " << endl;
// 	}
// 	myfile_out.close();
// 
// }

//--

//////------------------------------------------------------------------------//////

lh_square_current.f_U.assign(para_other.n,1.0);
lh_square_current.q_T.assign(para_other.n,0.0);
lh_square_current.kt_sum_U.assign(para_other.n,0.0);
lh_square_current.f_E.assign(para_other.n,1.0);
lh_square_current.g_E.assign(para_other.n,1.0);
lh_square_current.h_E.assign(para_other.n,1.0);
lh_square_current.k_sum_E.assign(para_other.n,0.0);
lh_square_current.q_E.assign(para_other.n,0.0);
lh_square_current.kt_sum_E.assign(para_other.n,0.0);
lh_square_current.f_I.assign(para_other.n,1.0);
lh_square_current.f_R.assign(para_other.n,1.0);
lh_square_current.f_EnI.assign(para_other.n,1.0);
lh_square_current.f_InR.assign(para_other.n,1.0);
lh_square_current.log_f_S.assign(para_other.n,0.0);
lh_square_current.log_f_Snull.assign(para_other.n,0.0);

FUNC func_mcmc;

func_mcmc.set_para(para_current, para_other, coordinate, xi_U_current, xi_E_current, xi_E_minus_current, xi_I_current, xi_R_current, xi_EnI_current, xi_InR_current, t_e_current, t_i_current, t_r_current, stb_current, index_current, infected_source_current);
func_mcmc.initialize_kernel_mat(kernel_mat_current, norm_const_current); // initialize the kernel matrix
func_mcmc.initialize_delta_mat(delta_mat_current); // initialize the kernel matrix
func_mcmc.initialize_lh_square(lh_square_current, kernel_mat_current, delta_mat_current,norm_const_current, nt_data_current, con_seq_current); //initialize lh_square

myfile_out.open((string(path4)+string("initial_norm_const.txt")).c_str(),ios::app);
for (int i=0;i<=((int)para_other.n-1);i++){
myfile_out<<norm_const_current.at(i)<<endl;
}
myfile_out.close();

double log_lh_current = log_lh_func (lh_square_current, para_other.n); // initialization of log-likelihood value

myfile_out.open((string(path4)+string("initial_lh.txt")).c_str(),ios::app);
myfile_out<<log_lh_current<<endl; //NOTE: this must be defined, otherwise has to re-initialize some components
myfile_out.close();

//-------------------
int total_count_1_initial, total_count_2_initial, total_count_3_initial;
total_count_1_initial=total_count_2_initial=total_count_3_initial=0;

count_type_all(nt_data_current,  xi_E_current, para_other.n_base, total_count_1_initial, total_count_2_initial, total_count_3_initial);

myfile_out.open((string(path4)+string("count_type_initial.txt")).c_str(),ios::app);
myfile_out<<total_count_1_initial <<"," << total_count_2_initial <<"," << total_count_3_initial <<endl;
myfile_out.close();

//-------------

myfile_out.open((string(path4)+string("initial_f_I.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_I_current.size()-1);i++){
myfile_out<<lh_square_current.f_I.at(xi_I_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_EnI.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_EnI_current.size()-1);i++){
myfile_out<<lh_square_current.f_EnI.at(xi_EnI_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_R.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_R_current.size()-1);i++){
myfile_out<<lh_square_current.f_R.at(xi_R_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_InR.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_InR_current.size()-1);i++){
myfile_out<<lh_square_current.f_InR.at(xi_InR_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_kt_sum_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.kt_sum_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_g_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.g_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();myfile_out.open((string(path4)+string("initial_h_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.h_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_q_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.q_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_E.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out<<lh_square_current.f_E.at(xi_E_minus_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_f_U.txt")).c_str(),ios::app);
for (int i=0;i<=((int)xi_U_current.size()-1);i++){
myfile_out<<lh_square_current.f_U.at(xi_U_current.at(i))<<endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_log_f_S.txt")).c_str(),ios::app);
for (int i=0;i<=(para_other.n-1);i++){
myfile_out<<lh_square_current.log_f_S.at(i)<<endl;
}
myfile_out.close();

/*--------------------*/
ifstream myfile1_in, myfile2_in, myfile3_in,  myfile4_in, myfile5_in, myfile6_in, myfile7_in,myfile8_in, myfile9_in, myfile10_in;
ofstream myfile1_out, myfile2_out, myfile3_out,  myfile4_out, myfile5_out, myfile6_out, myfile7_out,myfile8_out, myfile9_out, myfile10_out;

mcmc_UPDATE mcmc_update;
mcmc_update.set_para (para_other, coordinate);

myfile1_out.open((string(path4)+string("parameters_current.txt")).c_str(),ios::app); // for outputting current mcmc values of model parameters
myfile1_out << "alpha" << "," << "beta" << "," << "mu_lat" << "," << "var_lat" << "," << "c" << "," << "d" << "," << "k_1" << "," << "k_2" <<  "," << "stb_1" << "," << "stb_2" <<"," << "mu_1" << "," << "mu_2" << "," << "p_ber" << endl;
myfile2_out.open((string(path4)+string("lh_current.txt")).c_str(),ios::app); // .. current likelihood value
myfile4_out.open((string(path4)+string("infected_source_current.txt")).c_str(),ios::app); // ..current imputed infected sources
myfile10_out.open((string(path4)+string("con_seq_current.txt")).c_str(),ios::app);// .. current imputed master seq


// myfile3_out.open((string(path4)+string("count_type_check.txt")).c_str(),ios::app); // redundant
// myfile5_out.open((string(path4)+string("cover_rate_current.txt")).c_str(),ios::app);

// myfile6_out.open((string(path4)+string("index_first_seq_current.txt")).c_str(),ios::app); // record the first seq of index
// myfile7_out.open((string(path4)+string("imported_sampled_first_seq_current.txt")).c_str(),ios::app);// record the first seq of imported cases (those with samples)
// myfile8_out.open((string(path4)+string("imported_not_sampled_first_seq_current.txt")).c_str(),ios::app);// record the first seq of imported cases (those with samples)

// myfile9_out.open((string(path4)+string("ind_cryp_single.txt")).c_str(),ios::app);// record the first seq of imported cases (those with samples)
// myfile9_out << "t_i" <<"," << "current_size" << endl;


// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T

const gsl_rng_type* T_c_unvs= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c_unvs = gsl_rng_alloc (T_c_unvs); // r is pointer points to an object with Type T
gsl_rng_set (r_c_unvs, 1); // set a seed


vector<int> list_update; // would contain the subjects (descended from subject_proposed below) whose FIRST sequence would be updated, with a sequential order (i.e., level-wise and time-wise) of updating (note: as each event needed to be updated corresponds to an infection event, it would be sufficient to update the first sequence of necessary subjects so as to update all downstream seq)
list_update.reserve(1);

// ofstream myfile_temp_out; 
// myfile_temp_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// 

int n_iter = 1000000; //number of iterations for MCMC
int n_freq = 10; // frequency to translate an infection time


for (int i=0;i<=(n_iter-1);i++){


//gsl_rng_set (r_c, -1000*(i+1)); // set a seed

mcmc_update.k_1_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, infected_source_current,  norm_const_current, i, r_c_unvs); // mcmc update for spatial kerenel parameter



//mcmc_update.k_2_update(lh_square_current, log_lh_current, kernel_mat_current, delta_mat_current, xi_U_current, xi_E_minus_current, xi_I_current, t_r_current,  t_i_current, t_e_current, index_current, para_current, stb_current, infected_source_current,  norm_const_current, i, r_c_unvs); 



mcmc_update.beta_update(lh_square_current,log_lh_current, xi_U_current, xi_E_minus_current, t_e_current, index_current, stb_current,infected_source_current, para_current, i, r_c_unvs); // .. secondary transmisison rate

mcmc_update.alpha_update(lh_square_current, log_lh_current, xi_U_current, xi_E_minus_current, t_e_current,index_current, stb_current,infected_source_current, para_current, i, r_c_unvs); // .. primary rate

mcmc_update.mu_lat_update(lh_square_current, log_lh_current, xi_I_current, xi_EnI_current, t_i_current, t_e_current, index_current, para_current, i, r_c_unvs); // ...mean of latent period
mcmc_update.var_lat_update(lh_square_current, log_lh_current, xi_I_current, xi_EnI_current,t_i_current, t_e_current, index_current, para_current, i, r_c_unvs); // .. variance of latent peiod

mcmc_update.c_update(lh_square_current,log_lh_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, index_current, para_current, i, r_c_unvs); // .. parameter for infectious period
//mcmc_update.d_update(lh_square_current,log_lh_current, xi_R_current, xi_InR_current,t_r_current, t_i_current, index_current, para_current, i, r_c_unvs);


mcmc_update.mu_1_update(lh_square_current,log_lh_current, xi_E_current, para_current, nt_data_current, i, r_c_unvs); // .. mutation rate
mcmc_update.mu_2_update(lh_square_current,log_lh_current, xi_E_current, para_current, nt_data_current, i, r_c_unvs);

mcmc_update.p_ber_update(lh_square_current,log_lh_current, xi_E_current, para_current, nt_data_current, infected_source_current, con_seq_current,i, r_c_unvs); //.. variation parameter for the Master sequence G_M

for (int k=0;k<=(para_other.n_base-1);k++){
	mcmc_update.con_seq_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt, xi_beta_E_current, con_seq_current, k, r_c_unvs); // update the master sequence
}

//-----

//--------------------//

for (int j=0;j<=((int)(n_freq/10.0) - 1);j++){ // update sequence without changing t_e, source

	int subject_proposed_2 = xi_E_current.at(gsl_rng_uniform_int (r_c_unvs, xi_E_current.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn

	for (int k=0;k<=(para_other.n_base-1);k++){
	mcmc_update.seq_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt, xi_beta_E_current, con_seq_current,subject_proposed_2, k, r_c_unvs);
	}                                                                   

}

//--------------------//

for (int j=0;j<=(n_freq - 1);j++){

	int subject_proposed = xi_E_current.at(gsl_rng_uniform_int (r_c_unvs, xi_E_current.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn

	mcmc_update.t_e_seq(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt, nt_data_current.infecting_list,  nt_data_current.infecting_size, xi_beta_E_current, con_seq_current, subject_proposed, i+1, r_c_unvs); // update t_e and sequence jointly

	mcmc_update.source_t_e_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt,  nt_data_current.infecting_list,  nt_data_current.infecting_size, xi_beta_E_current,  subject_proposed, list_update, con_seq_current, i+1, r_c_unvs); // update t_e, source, sequence, jointly

}



for (int jk=0;jk<=((int)index_current.size()-1);jk++){// update the  first seq of indexes
	for (int k=0;k<=(para_other.n_base-1);k++){
	mcmc_update.seq_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, t_r_current,  t_i_current,  t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt, xi_beta_E_current, con_seq_current, index_current.at(jk), k, r_c_unvs);
	}                                                                   

}


mcmc_update.t_i_update(lh_square_current, log_lh_current,kernel_mat_current,  delta_mat_current,xi_U_current, xi_E_current,  xi_E_minus_current, xi_I_current, xi_EnI_current, xi_R_current, xi_InR_current,  t_r_current,  t_i_current, t_onset, t_e_current, index_current, para_current, stb_current, norm_const_current, infected_source_current, nt_data_current.t_sample, nt_data_current.current_size,  nt_data_current.nt, nt_data_current.t_nt, nt_data_current.infecting_list,  nt_data_current.infecting_size,  i+1);


//--------------

myfile1_out << para_current.alpha << "," << para_current.beta << "," << para_current.mu_lat << "," << para_current.var_lat << "," << para_current.c << "," << para_current.d << "," << para_current.k_1 << "," <<  para_current.k_2 << "," << para_current.stb_1<< "," <<  para_current.stb_2 <<  "," << para_current.mu_1<< "," <<  para_current.mu_2 << "," << para_current.p_ber << endl;
//myfile1_out.close();

myfile2_out <<  log_lh_current << endl;


div_t div_iter_lh;
div_iter_lh = div (i,100);
if (div_iter_lh.rem==0){ // for outputting infected_source (only at every 100 iterations, which saves space..)

	for (int js=0;js<=(para_other.n-1);js++){
	int rem = (js+1)%para_other.n;
	if ((rem!=0) | (js==0)) myfile4_out<< infected_source_current.at(js) << ",";
	if ((rem==0) & (js!=0)) myfile4_out <<  infected_source_current.at(js) << " " << endl;
	}


}

//------------------

div_t div_iter_seq;
div_iter_seq = div (i,2000);

switch (div_iter_seq.rem==0){ // for outputting the master sequence (only at every 2000 iterations, which saves space..)


	case 1:{

		for (int js=0;js<=(para_other.n_base-1);js++){
			int rem = (js+1)%para_other.n_base;
			if ((rem!=0) | (js==0)) myfile10_out << con_seq_current.at(js)<< ",";
			if ((rem==0) & (js!=0)) myfile10_out <<  con_seq_current.at(js)  << " " << endl;
		}

	
	break;
	}
	
	default:{
	break;
	}

}


} // end of MCMC loop

gsl_rng_free(r_c_unvs);

myfile1_out.close();
myfile2_out.close();
// myfile3_out.close();
myfile4_out.close();
// myfile5_out.close();
// myfile6_out.close();
// myfile7_out.close();
// myfile8_out.close();
// myfile9_out.close();
myfile10_out.close();


//------------------


return(0);
}


