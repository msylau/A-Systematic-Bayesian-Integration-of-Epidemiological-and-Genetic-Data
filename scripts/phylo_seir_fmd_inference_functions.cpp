
#include "phylo_seir_fmd_inference_header.h"


inline long double func_kernel (double x_1 , double y_1, double x_2, double y_2, double k_1_arg, double k_2_arg, const string& kernel_type_arg){

//long double eucli_dist = sqrt( (x_1-x_2)*(x_1-x_2) + (y_1-y_2)*(y_1-y_2) );

double pi = 3.1415926535897;
double rad = pi/180.0;
double a1 = x_1 * rad;
double a2 = y_1 * rad;
double b1 = x_2 * rad;
double b2 = y_2 * rad;
double dlon = b2 - a2;
double dlat = b1 - a1;
double a = pow((sin(dlat/2.0)),2.0) + cos(a1) * cos(b1) * pow((sin(dlon/2.0)),2.0);
double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
double R = 6378.145;
long double eucli_dist  = R * c;

long double func_ker = exp((-k_1_arg)*eucli_dist); //exp
//long double func_ker = 1.0/pow(eucli_dist, k_2_arg); //power-law


return(func_ker);
}

/*-------------------------------------------*/

inline double lh_snull(const vector<int>& con_seq, const vector<int>& seq, const double& p_ber, const int& n_base){ // compute the log pr a seq for background

	double lh_snull=0.0;

	int m=0;

	for (int i=0;i<=(n_base-1);i++){
		switch(seq.at(i)!=con_seq.at(i)){
			case 1:{ // a change
				m = m +1;
			break;
			}
			case 0:{ // not a change
			break;
			}
		}
	}

	lh_snull = m*log(p_ber) + (n_base-m)*log(1-p_ber) + m*log(1.0/3.0); 
//	lh_snull = m*log(p_ber) + (n_base-m)*log(1-p_ber) ;

	return(lh_snull);
}

/*-------------------------------------------------------*/

inline double lh_snull_base(const int& con_base, const int& base, const double& p_ber){ // compute the log pr a base for background

	double lh_snull_base=0.0;

	int m=0;

		switch(base!=con_base){
			case 1:{ // a change
				m = m +1;
			break;
			}
			case 0:{ // not a change
			break;
			}
		}

	lh_snull_base = m*log(p_ber) + (1-m)*log(1-p_ber) + m*log(1.0/3.0);
	//lh_snull_base = m*log(p_ber) + (1-m)*log(1-p_ber) ;

	return(lh_snull_base);
}

/*-------------------------------------------------------*/

inline void sample_snull (const vector<int>& con_seq, vector<int>& seq_proposed, const double& p_ber, const int& n_base, gsl_rng* r){ //sample a seq for background
	
for (int j=0;j<=(n_base-1);j++){

	int ber_trail =  gsl_ran_bernoulli (r,p_ber); // return 1 if a change is to happen
	int base_proposed=0;// any output of 0 would indicate a mistake


	switch(ber_trail==1){
		case 1:{ // randomly choose one among other 3
			switch(con_seq.at(j)){
				case 1:{
					int type = gsl_rng_uniform_int (r, 3);
					switch(type){
						case 0:{
							base_proposed = 2;
						break;
						}
						case 1:{
							base_proposed = 3;
						break;
						}
						case 2:{
							base_proposed = 4;
						break;
						}
					}
				break;
				}
				case 2:{
					int type = gsl_rng_uniform_int (r, 3);
			
					switch(type){
						case 0:{
							base_proposed = 1;
						break;
						}
						case 1:{
							base_proposed = 3;
						break;
						}
						case 2:{
							base_proposed = 4;
						break;
						}
					}	
				break;
				}
				case 3:{
					int type = gsl_rng_uniform_int (r, 3);
					switch(type){
						case 0:{
							base_proposed = 1;
						break;
						}
						case 1:{
							base_proposed = 2;
						break;
						}
						case 2:{
							base_proposed = 4;
						break;
						}
					}	
				break;
				}
				case 4:{
					int type = gsl_rng_uniform_int (r, 3);
					switch(type){
						case 0:{
							base_proposed = 1;
						break;
						}
						case 1:{
							base_proposed = 2;
						break;
						}
						case 2:{
							base_proposed = 3;
						break;
						}
					}	
				break;
				}
			}

			seq_proposed.at(j) = base_proposed; 

		break;
		}
		case 0:{
			seq_proposed.at(j) = con_seq.at(j); // same as consensus seq
		break;
		}
	}
	}

}

/*-------------------------------------------------------*/


void seq_propose_cond(vector<int>& seq_proposed, double& log_pr_forward, const vector<int>&  nt_past_forward, const vector<int>& nt_future_forward, const double& t_proposed, const double& t_past, const double& t_future, const double& mu_1, const double& mu_2, int n_base, gsl_rng* r_c){

double T = fabs(t_future - t_past);
double dt = fabs(t_proposed - t_past);
double p = dt/T; // pr of assigning a base to be same as the corresponding one in the future sequence

// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("dt.txt")).c_str(),ios::app);
// 		myfile_out <<dt << "," << T << ","<< T - dt << endl;
// 		myfile_out.close();	


double P[2] = {1.0-p, p};

gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (sizeof(P)/sizeof(P[0]),P);

int m=0; // number of sites from nt_past were different from the corresponding sites from nt_future
int dn=0; // number of sites from nt_past become the same as the coresponding sites from nt_future among m

for (int i=0; i<=(n_base-1); i++){

	switch(nt_past_forward.at(i)==nt_future_forward.at(i)){

		case 0:{ // not the same

			m = m + 1;
			int bool_c = gsl_ran_discrete (r_c, g); //  1 = assign to be the same as future

// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("bool_c.txt")).c_str(),ios::app);
// 		myfile_out <<bool_c << endl;
// 		myfile_out.close();	


			switch(bool_c){
				case 1:{
					dn =dn +1;
					seq_proposed.at(i)=nt_future_forward.at(i);
				break;
				}
				case 0:{
					seq_proposed.at(i)=nt_past_forward.at(i);
				break;
				}
			}

		break;
		}

		case 1:{ // same
			seq_proposed.at(i)=nt_past_forward.at(i);
		break;
		}

	}

}

	//log_pr_forward = log((long double)pow(p, dn)) + log((long double)pow(1.0-p, m-dn));
	log_pr_forward =dn*log(p) + (m-dn)*log(1.0-p);

//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_pr_forward_cond.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  (long double)pow(p, dn) << ", "<< (long double)pow(1.0-p, m-dn) <<","<< log_pr_forward << endl;
// 		//myfile_mcmc_out <<  p << ", "<< dn <<","<< log_pr_forward << endl;
// 		//myfile_mcmc_out <<  dt << ", "<< T <<","<< log_pr_forward << endl;
// 		//myfile_mcmc_out <<  t_proposed << ", "<< t_past <<","<< log_pr_forward << endl;
// 
// 		myfile_mcmc_out.close();

	gsl_ran_discrete_free (g);

}									

/*-------------------------------------------*/

void seq_propose_uncond(vector<int>& seq_proposed, double& log_pr_forward,  const vector<int>& nt_past_forward, const double& t_proposed, const double& t_past, const double& t_future,  const double& mu_1, const double& mu_2, int n_base, gsl_rng * r_c){

int n_1, n_2, n_3; // number of unchanged, transition, two types of transversion (compare original nt and nt_proposed)
n_1=n_2=n_3= 0;

double p_1, p_2, p_3;

double abs_dt = fabs(t_future - t_past);
//double abs_dt = fabs(t_proposed - t_past);


p_1 = 0.25 + 0.25*exp(-4.0*mu_2*abs_dt) + 0.5*exp(-2.0*(mu_1+mu_2)*abs_dt) ; // pr of a base not changing
p_2 = 0.25 + 0.25*exp(-4.0*mu_2*abs_dt) - 0.5*exp(-2.0*(mu_1+mu_2)*abs_dt); // pr of a transition of a base
p_3 = 2.0*(0.25 - 0.25*exp(-4.0*mu_2*abs_dt));  // pr of a transversion (two possible events)


// p_1 = 0.25;
// p_2 = 0.25;
// p_3 = 0.5;


double P[3] = {p_1, p_2, p_3};

gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (sizeof(P)/sizeof(P[0]),P);

for (int j=0;j<=(n_base-1);j++){

	int type= gsl_ran_discrete (r_c, g) + 1; 


	switch(nt_past_forward.at(j)){

	case 1:{ // an A

		switch(type){
		case 1:{
		seq_proposed.at(j) = nt_past_forward.at(j);
		n_1 = n_1 + 1;
		break;
		}
		case 2:{
		seq_proposed.at(j) = 2; 
		n_2 = n_2 + 1;
		break;
		}		
		case 3:{
		n_3 = n_3 + 1;

		int type_trans = gsl_rng_uniform_int (r_c, 2);//  uniformly drawn from [0,2) to determine the exact type of transversion

			switch(type_trans){
				case 0:{
				seq_proposed.at(j) = 3;
				break;
				}
				case 1:{
				seq_proposed.at(j) = 4;
				break;
				}
			}	

		break;
		}			
		}
	break;
	}
	
	case 2:{ // a G
		switch(type){
		case 1:{
		seq_proposed.at(j) = nt_past_forward.at(j);
		n_1 = n_1 + 1;
		break;
		}
		case 2:{
		seq_proposed.at(j) = 1;
		n_2 = n_2 + 1;
		break;
		}		
		case 3:{
		n_3 = n_3 + 1;

		int type_trans = gsl_rng_uniform_int (r_c, 2);//  uniformly drawn from [0,2) to determine the exact type of transversion

			switch(type_trans){
				case 0:{
				seq_proposed.at(j) = 3;
				break;
				}
				case 1:{
				seq_proposed.at(j) = 4;
				break;
				}
			}	

		break;
		}				
		}
	break;
	}
	
	case 3:{ // a T
		switch(type){
		case 1:{
		seq_proposed.at(j) = nt_past_forward.at(j);
		n_1 = n_1 + 1;
		break;
		}
		case 2:{
		seq_proposed.at(j) = 4;
		n_2 = n_2 + 1;
		break;
		}		
		case 3:{
		n_3 = n_3 + 1;

		int type_trans = gsl_rng_uniform_int (r_c, 2);//  uniformly drawn from [0,2) to determine the exact type of transversion

			switch(type_trans){
				case 0:{
				seq_proposed.at(j) = 1;
				break;
				}
				case 1:{
				seq_proposed.at(j) = 2;
				break;
				}
			}	

		break;
		}			
		}
	break;
	}

	case 4:{ // a C
		switch(type){
		case 1:{
		seq_proposed.at(j) = nt_past_forward.at(j);
		n_1 = n_1 + 1;
		break;
		}
		case 2:{
		seq_proposed.at(j) = 3;
		n_2 = n_2 + 1;
		break;
		}		
		case 3:{
		n_3 = n_3 + 1;

		int type_trans = gsl_rng_uniform_int (r_c, 2);//  uniformly drawn from [0,2) to determine the exact type of transversion

			switch(type_trans){
				case 0:{
				seq_proposed.at(j) = 1;
				break;
				}
				case 1:{
				seq_proposed.at(j) = 2;
				break;
				}
			}	

		break;
		}				
		}
	break;
	}

	}

}
	
	
	//log_pr_forward = log((long double)pow(p_1, n_1)) + log((long double)pow(p_2, n_2)) + log((long double)pow(p_3, n_3));
	log_pr_forward = n_1*log(p_1) +n_2*log(p_2) + n_3*log(p_3);


//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_pr_forward_uncond.txt")).c_str(),ios::app);
// 		//myfile_mcmc_out <<  (long double)pow(p_1, n_1) << ", "<<  (long double)pow(p_2, n_2)<< ", "<<  (long double)pow(p_3, n_3)<< ","<< log_pr_forward << endl;
// 		//myfile_mcmc_out << n_1 << ", "<< n_2 << ", "<<  n_3 << ","<< log_pr_forward << endl;
// 		myfile_mcmc_out << p_1 << ", "<< p_2 << ", "<<  p_3 << ","<< p_4 <<","<< log_pr_forward << endl;
// 		//myfile_mcmc_out << p_1 << ", "<< n_1<< ","<< log_pr_forward << endl;
// 
// 		myfile_mcmc_out.close();


	gsl_ran_discrete_free (g);

}
/*-------------------------------------------*/

void seq_backward_pr_cond(const vector<int>& seq_proposed_backward, double& log_pr_backward, const vector<int>& nt_past_backward, const vector<int>& nt_future_backward, const double& t_proposed_backward, const double& t_past_backward,  const double& t_future_backward, const double& mu_1, const double& mu_2, int n_base){

double dt = fabs(t_proposed_backward - t_past_backward);
double T =  fabs(t_future_backward - t_past_backward);

double p = dt/T;

int m=0; // number of sites from nt_past were different from the corresponding sites from nt_future
int dn=0; // number of sites from nt_past become the same as the coresponding sites from nt_future among m

for (int i=0; i<=(n_base-1); i++){

	switch((nt_past_backward.at(i)==seq_proposed_backward.at(i)) & (seq_proposed_backward.at(i)==nt_future_backward.at(i))){

		case 0:{// not three bases are the same
			m = m +1;

			switch(seq_proposed_backward.at(i)==nt_future_backward.at(i)){
				case 1:{// it is the same as future base but not the same as past (i.e.,a change)
					dn = dn + 1;
				break;
				}
				case 0:{
					//do nothing
				break;
				}
			}

		break;
		}

		case 1:{
			// do nothing
		break;
		}

	}
}

//log_pr_backward = log((long double)pow(p, dn)) + log((long double)pow(1.0-p, m-dn));
log_pr_backward = dn*log(p) + (m-dn)*log(1.0-p);

//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_pr_backward_cond.txt")).c_str(),ios::app);
// 		//myfile_mcmc_out <<  (long double)pow(p, dn) << ", "<< (long double)pow(1.0-p, m-dn) <<","<< log_pr_forward << endl;
// 		myfile_mcmc_out <<  p << ", "<< dn <<","<< log_pr_backward << endl;
// 		//myfile_mcmc_out <<  dt << ", "<< T <<","<< log_pr_forward << endl;
// 		//myfile_mcmc_out <<  t_proposed << ", "<< t_past <<","<< log_pr_forward << endl;
// 
// 		myfile_mcmc_out.close();


}
/*-------------------------------------------*/

void seq_backward_pr_uncond(const vector<int>& seq_proposed_backward,  double& log_pr_backward, const vector<int>& nt_past_backward, const double& t_proposed_backward, const double& t_past_backward, const double& t_future_backward,  const double& mu_1, const double& mu_2, int n_base){


int count_1, count_2, count_3;

count_1=count_2=count_3=0;

double p_1, p_2, p_3;

double abs_dt = fabs(t_future_backward - t_past_backward);
//double abs_dt = fabs(t_proposed_backward - t_past_backward);


p_1 = 0.25 + 0.25*exp(-4.0*mu_2*abs_dt) + 0.5*exp(-2.0*(mu_1+mu_2)*abs_dt) ; // pr of a base not changing
p_2 = 0.25 + 0.25*exp(-4.0*mu_2*abs_dt) - 0.5*exp(-2.0*(mu_1+mu_2)*abs_dt); // pr of a transition of a base
p_3 = 2.0*(0.25 - 0.25*exp(-4.0*mu_2*abs_dt));  // pr of a transversion (two possible events)

// p_1 = 0.25;
// p_2 = 0.25;
// p_3 = 0.5;


for ( int i=0;i<=(n_base-1); i++){

	switch(abs(nt_past_backward.at(i)-seq_proposed_backward.at(i))){
	case 0:{
	count_1 = count_1 + 1;
	break;
	}
	case 1:{
		switch( ((nt_past_backward.at(i)==2) & (seq_proposed_backward.at(i)==3)) | ((nt_past_backward.at(i)==3) & (seq_proposed_backward.at(i)==2)) ){
		case 1:{
		count_3 = count_3 + 1;
		break;
		}
		case 0:{
		count_2 = count_2 + 1;
		break;
		}
		}
	break;
	}
	case 2:{
	count_3 = count_3 + 1;
	break;
	}
	case 3:{
	count_3 = count_3 + 1;
	break;
	}
	}

}

	//log_pr_backward = log((long double)pow(p_1, count_1)) + log((long double)pow(p_2, count_2)) + log((long double)pow(p_3, count_3));
	log_pr_backward = count_1*log(p_1) + count_2*log(p_2) + count_3*log(p_3) ; 


//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_pr_backward_uncond.txt")).c_str(),ios::app);
// 		//myfile_mcmc_out <<  (long double)pow(p_1, n_1) << ", "<<  (long double)pow(p_2, n_2)<< ", "<<  (long double)pow(p_3, n_3)<< ","<< log_pr_forward << endl;
// 		//myfile_mcmc_out << n_1 << ", "<< n_2 << ", "<<  n_3 << ","<< log_pr_forward << endl;
// 		myfile_mcmc_out << p_1 << ", "<< p_2 << ", "<<  p_3 << ","<< log_pr_backward << endl;
// 		//myfile_mcmc_out << p_1 << ", "<< n_1<< ","<< log_pr_forward << endl;
// 
// 		myfile_mcmc_out.close();


}

double log_lh_base (int& base_1, int& base_2, double t_1_arg, double t_2_arg , double mu_1_arg, double mu_2_arg){

double log_lh=-99.0;;

double dt = t_2_arg - t_1_arg;

double p_1 = 0.25 + 0.25*exp(-4.0*mu_2_arg*dt) + 0.5*exp(-2.0*(mu_1_arg+mu_2_arg)*dt) ; // pr of a base not changing
double p_2 = 0.25 + 0.25*exp(-4.0*mu_2_arg*dt) - 0.5*exp(-2.0*(mu_1_arg+mu_2_arg)*dt); // pr of a transition of a base
double p_3 =1.0*(0.25 - 0.25*exp(-4.0*mu_2_arg*dt));  // pr of a transversion (two possible events)
//double p_1 = 1.0 - p_2 - 2.0*p_3; // pr of a base not changing


	switch(abs(base_1-base_2)){
	case 0:{
	log_lh = log(p_1);
	break;
	}
	case 1:{
		switch( ((base_1==2) & (base_2==3)) | ((base_1==3) & (base_2==2)) ){
		case 1:{
		log_lh = log(p_3);
		break;
		}
		case 0:{
		log_lh = log(p_2);
		break;
		}
		}
	break;
	}
	case 2:{
		log_lh = log(p_3);
	break;
	}
	case 3:{
		log_lh = log(p_3);
	break;
	}
	}


//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_lh_base.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_lh << endl;
// 		myfile_mcmc_out.close();
// 

return(log_lh);

}

/*-------------------------------------------*/

double log_lh_seq (vector<int>& seq_1_arg, vector<int>& seq_2_arg, double t_1_arg, double t_2_arg , double mu_1_arg, double mu_2_arg, int n_base_arg){

double dt = t_2_arg - t_1_arg;


int count_1, count_2, count_3; // count_1=count of unchanged sites, ..transition.., transversion
count_1=count_2=count_3=0; 

//vector<int> type(n_base_arg);


for ( int i=0;i<=(n_base_arg-1); i++){

	switch(abs(seq_1_arg.at(i)-seq_2_arg.at(i))){
	case 0:{
	count_1 = count_1 + 1;
	//type.at(i) = 1;
	break;
	}
	case 1:{
		switch( ((seq_1_arg.at(i)==2) & (seq_2_arg.at(i)==3)) | ((seq_1_arg.at(i)==3) & (seq_2_arg.at(i)==2)) ){
		case 1:{
		count_3 = count_3 + 1;
		//type.at(i) = 3;
		break;
		}
		case 0:{
		count_2 = count_2 + 1;
		//type.at(i) = 2;
		break;
		}
		}
	break;
	}
	case 2:{
	count_3 = count_3 + 1;
	//type.at(i) = 3;
	break;
	}
	case 3:{
	count_3 = count_3 + 1;
	//type.at(i) = 3;
	break;
	}
	}

}


double p_1 = 0.25 + 0.25*exp(-4.0*mu_2_arg*dt) + 0.5*exp(-2.0*(mu_1_arg+mu_2_arg)*dt) ; // pr of a base not changing
double p_2 = 0.25 + 0.25*exp(-4.0*mu_2_arg*dt) - 0.5*exp(-2.0*(mu_1_arg+mu_2_arg)*dt); // pr of a transition of a base
double p_3 =1.0*(0.25 - 0.25*exp(-4.0*mu_2_arg*dt));  // pr of a transversion (two possible events)
//double p_1 = 1.0 - p_2 - 2.0*p_3; // pr of a base not changing

double log_lh = count_1*log(p_1) + count_2*log(p_2) + count_3*log(p_3);


// double P[3] = {p_1, p_2, 2.0*p_3};
// unsigned int C[3] = {count_1, count_2, count_3};
// double log_lh = gsl_ran_multinomial_lnpdf (3, P, C);


// 		int total = count_1 + count_2 + count_3;
// 		double p_total = p_1 + p_2 + 2.0*p_3;
// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("lh_seq.txt")).c_str(),ios::app);
// 		myfile_out <<log_lh<< ","<< total << "," << p_total << "," << count_1 <<"," << count_2 << "," << count_3 << "," << p_1 <<"," << p_2 << "," << 2.0*p_3 << endl;
// 		myfile_out.close();	

// 		myfile_out.open((string(path4)+string("test_2_type.txt")).c_str(),ios::app);
// 		for (int i=0;i<=(n_base_arg-1);i++){
// 		if (i<(n_base_arg-1))  myfile_out << type[i] << ",";
// 		if (i==(n_base_arg-1))   myfile_out << type[i] << " " << endl;
// 		}
// 		myfile_out.close();
// 	
// 		myfile_out.open((string(path4)+string("test_3_abs.txt")).c_str(),ios::app);
// 		for (int i=0;i<=(n_base_arg-1);i++){
// 		myfile_out << abs(seq_1_arg.at(i) - seq_2_arg.at(i)) << "," << fabs(seq_1_arg.at(i) - seq_2_arg.at(i)) << endl;
// 		}
// 		myfile_out.close();


//return(lh);
return(log_lh);

}

/*-------------------------------------------*/

double dtnorm(double x, double mean, double sd, double a){ // pdf of truncated normal with lower bound = a; upper bound =Inf

double num = (1.0/sd)*gsl_ran_ugaussian_pdf((x-mean)/sd);

double denom = 1.0 - gsl_cdf_ugaussian_P((a-mean)/sd);

double d = num/denom;

return(d);

}


/*-------------------------------------------*/
long double func_latent_pdf(double t , double mu_lat, double var_lat){

long double func_lat_pdf;

// // double a_lat = mu_lat*mu_lat/var_lat; // when use GAMMA latent
// // double b_lat = var_lat/mu_lat;

func_lat_pdf = gsl_ran_gamma_pdf(t, mu_lat, var_lat);

// // double a_lat = log(mu_lat*mu_lat/(sqrt(var_lat +mu_lat*mu_lat))); // when use lognormal latent
// // double b_lat = sqrt(log(var_lat/(mu_lat*mu_lat)+1));
// // func_lat_pdf = gsl_ran_lognormal_pdf(t, a_lat, b_lat);

// // double a_lat = mu_lat; // when use exponential latent
// //func_lat_pdf = gsl_ran_exponential_pdf(t, a_lat);

//func_lat_pdf = gsl_ran_exponential_pdf(t, mu_lat);


return(func_lat_pdf);
}
/*-------------------------------------------*/

inline long double func_latent_cdf(double t , double mu_lat, double var_lat){

long double func_lat_cdf;

// // double a_lat = mu_lat*mu_lat/var_lat; // when use GAMMA latent
// // double b_lat = var_lat/mu_lat;

func_lat_cdf = gsl_cdf_gamma_P(t, mu_lat, var_lat);

// // double a_lat = log(mu_lat*mu_lat/(sqrt(var_lat +mu_lat*mu_lat))); // when use lognormal latent
// // double b_lat = sqrt(log(var_lat/(mu_lat*mu_lat)+1));
// // func_lat_cdf = gsl_cdf_lognormal_P(t, a_lat, b_lat);

// // double a_lat = mu_lat; // when use exponential latent
// // func_lat_cdf = gsl_cdf_exponential_P(t, a_lat);

//func_lat_cdf = gsl_cdf_exponential_P(t, mu_lat);


return(func_lat_cdf);
}
/*-------------------------------------------*/

double log_lh_func (lh_SQUARE lh_square_arg, int n_arg) {

double log_lh_value =0.0;

for (int i=0; i<=(n_arg-1);i++){
log_lh_value = log_lh_value +  lh_square_arg.log_f_Snull.at(i) +  lh_square_arg.log_f_S.at(i) + log(lh_square_arg.f_U.at(i)) +log(lh_square_arg.f_E.at(i)) + log(lh_square_arg.f_I.at(i)) +log(lh_square_arg.f_R.at(i)) + log(lh_square_arg.f_EnI.at(i)) + log(lh_square_arg.f_InR.at(i));
}

return(log_lh_value);

}


/*------------------------------------------------*/
void FUNC::initialize_kernel_mat (vector< vector<double> >& kernel_mat_arg, vector<double>& norm_const_arg) {


for (int i=0;i<=(n_Clh-1);i++) {
 for (int j=0;j<=(n_Clh-1);j++) {
 if (i==j) kernel_mat_arg[i][j]=0.0;
 if (i<j) kernel_mat_arg[i][j] = func_kernel (coordinate_Clh[i][0],coordinate_Clh[i][1],coordinate_Clh[j][0],coordinate_Clh[j][1],k_1_Clh,k_2_Clh,kernel_type_Clh);
 if (i>j) kernel_mat_arg[i][j]=kernel_mat_arg[j][i];
 }
}

for (int j=0;j<=(n_Clh-1);j++) {
norm_const_arg.at(j) = 0.0;
 for (int i=0;(i<=(n_Clh-1)); i++) {
norm_const_arg.at(j) = norm_const_arg.at(j) +  kernel_mat_arg[i][j];
}
}

}

/*------------------------------------------------*/
void FUNC::initialize_delta_mat (vector< vector<double> >& delta_mat_arg){
 

if (xi_U_Clh.size()>=1){
for (int i=0;i<= (int)(xi_U_Clh.size()-1);i++){

    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
    
    switch (t_r_Clh.at(xi_I_Clh.at(j))>t_max_Clh) {
    case 1:{ // not yet recovered
    delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)] = t_max_Clh - t_i_Clh.at(xi_I_Clh.at(j));
    break;
    }
    case 0:{ // recovered
    delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)] = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
    break;
    }    
    }

    }
}
}

//----------//


if (xi_E_minus_Clh.size()>=1){
for (int i=0;i<= (int)(xi_E_minus_Clh.size()-1);i++){
   
    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
   
        if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_minus_Clh.at(i))) { 

        switch (t_r_Clh.at(xi_I_Clh.at(j))>=t_e_Clh.at(xi_E_minus_Clh.at(i))) {
        case 1:{ // not yet recovered at e_i
        delta_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)] = t_e_Clh.at(xi_E_minus_Clh.at(i)) - t_i_Clh.at(xi_I_Clh.at(j));
        break;
        }
        case 0:{ // recovered before e_i
        delta_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)] = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
        break;
        }    
        }

        } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 


    } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)

}
}

for (int i=0;i<= (int)(index_Clh.size()-1);i++){
for (int j=0;j<= (int) (n_Clh-1);j++){
delta_mat_arg[index_Clh.at(i)][j] = 0.0;
}
}

}
/*------------------------------------------------*/

void FUNC::initialize_lh_square (lh_SQUARE& lh_square_arg, vector< vector<double> > kernel_mat_arg, vector< vector<double> > delta_mat_arg, vector<double>& norm_const_arg, nt_struct& nt_data_arg, vector<int>& con_seq){


if (xi_U_Clh.size()>=1){
for (int i=0;i<= (int)(xi_U_Clh.size()-1);i++){

    for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
 
       double delta_t = delta_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)]; 
//     double delta_t = 0.0;
//     switch (t_r_Clh.at(xi_I_Clh.at(j))>t_max_Clh) {
//     case 1:{ // not yet recovered
//     delta_t = t_max_Clh - t_i_Clh.at(xi_I_Clh.at(j));
//     break;
//     }
//     case 0:{ // recovered
//     delta_t = t_r_Clh.at(xi_I_Clh.at(j)) - t_i_Clh.at(xi_I_Clh.at(j));
//     break;
//     }    
//     }

    lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)) = lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i)) + delta_t*kernel_mat_arg[xi_U_Clh.at(i)][xi_I_Clh.at(j)]/norm_const_arg.at(xi_I_Clh.at(j));
    }

lh_square_arg.q_T.at(xi_U_Clh.at(i)) = alpha_Clh*t_max_Clh + beta_Clh*stb_Clh.at(xi_U_Clh.at(i))*lh_square_arg.kt_sum_U.at(xi_U_Clh.at(i));

lh_square_arg.f_U.at(xi_U_Clh.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_arg.q_T.at(xi_U_Clh.at(i)),1.0);


}
}

//------------------------------------------//

for (int i=0;i<= (int)(xi_E_Clh.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_Clh.at(i);

	switch (nt_data_arg.current_size.at(k_E)>1) {

		case 1:{
	
	
	// 	vector<int> seq_1(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_Clh);
	// 	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+n_base_Clh, nt_data_arg.nt[k_E].begin()+2*n_base_Clh);
	
		for (int j=0;j<=(nt_data_arg.current_size.at(k_E)-2);j++){
	
		vector<int> seq_1(nt_data_arg.nt[k_E].begin()+j*(n_base_Clh), nt_data_arg.nt[k_E].begin()+(j+1)*(n_base_Clh));
		vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*(n_base_Clh), nt_data_arg.nt[k_E].begin()+(j+2)*(n_base_Clh));
		lh_square_arg.log_f_S.at(k_E) =lh_square_arg.log_f_S.at(k_E)+ log_lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_Clh, mu_2_Clh, n_base_Clh);
	
		}
	
		break;
		}

		case 0:{	
		break;
		}
	}

	//--
	switch (infected_source_Clh.at(k_E)==9999) {
		case 1:{
// 			lh_square_arg.log_f_Snull.at(k_E) = n_base_Clh*log(0.25); // assume background infection gives a random sequence from stationary dist of seq
			vector<int> seq(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_Clh); // the first seq
			lh_square_arg.log_f_Snull.at(k_E) = lh_snull(con_seq, seq, p_ber_Clh, n_base_Clh); // compute the log pr a seq for background

		break;
		}
			
		case 0:{
		break;
		}
	}

	//--

}

//------------------------------------------//


if (xi_E_minus_Clh.size()>=1){

for (int i=0;i<= (int)(xi_E_minus_Clh.size()-1);i++){
   
	for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++){
	
		if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_minus_Clh.at(i))) { 
	
		double delta_t = delta_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)]; 

// 		switch (t_r_Clh.at(xi_I_Clh.at(j))>=t_e_Clh.at(xi_E_minus_Clh.at(i))) {
// 		case 1:{ // not yet recovered at e_i
// 		lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) = lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) + kernel_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)]/norm_const_arg.at(xi_I_Clh.at(j)); // update k_sum_E
// 		break;
// 		}
// 		case 0:{ // recovered before e_i
// 		break;
// 		}    
// 		}
	
		lh_square_arg.kt_sum_E.at(xi_E_minus_Clh.at(i)) = lh_square_arg.kt_sum_E.at(xi_E_minus_Clh.at(i)) + delta_t*kernel_mat_arg[xi_E_minus_Clh.at(i)][xi_I_Clh.at(j)]/norm_const_arg.at(xi_I_Clh.at(j)); // update kt_sum_E
	
		} // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 
	
	
	} // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)

// lh_square_arg.g_E.at(xi_E_minus_Clh.at(i)) = alpha_Clh + beta_Clh*stb_Clh.at(xi_E_minus_Clh.at(i))*lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i));


	switch(infected_source_Clh.at(xi_E_minus_Clh.at(i))){
	case 9999:{ // by background
	lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) = 0.0; // update k_sum_E
	lh_square_arg.g_E.at(xi_E_minus_Clh.at(i)) = alpha_Clh;
	break;
	}

	// 	case -99:{
	// 	ofstream myfile_out; 
	// 	myfile_out.open((string(path4)+string("warnings.txt")).c_str(),ios::app);
	// 	myfile_out <<xi_E_minus_Clh.at(i)<<endl;
	// 	myfile_out.close();	
	// 	break;
	// 	}

	default :{ // not by background
	lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i)) = kernel_mat_arg[xi_E_minus_Clh.at(i)][infected_source_Clh.at(xi_E_minus_Clh.at(i))]/norm_const_arg.at(infected_source_Clh.at(xi_E_minus_Clh.at(i))); // update k_sum_E
	lh_square_arg.g_E.at(xi_E_minus_Clh.at(i)) = beta_Clh*stb_Clh.at(xi_E_minus_Clh.at(i))*lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i));
	break;
	}

	}


// 	ofstream myfile_out; 
// 	myfile_out.open((string(path4)+string("k_sum_E.txt")).c_str(),ios::app);
// 	myfile_out <<lh_square_arg.k_sum_E.at(xi_E_minus_Clh.at(i))<<endl;
// 	myfile_out.close();

lh_square_arg.q_E.at(xi_E_minus_Clh.at(i)) = alpha_Clh*t_e_Clh.at(xi_E_minus_Clh.at(i)) + beta_Clh*stb_Clh.at(xi_E_minus_Clh.at(i))*lh_square_arg.kt_sum_E.at(xi_E_minus_Clh.at(i));
lh_square_arg.h_E.at(xi_E_minus_Clh.at(i)) = gsl_ran_exponential_pdf(lh_square_arg.q_E.at(xi_E_minus_Clh.at(i)),1.0);

lh_square_arg.f_E.at(xi_E_minus_Clh.at(i)) = lh_square_arg.g_E.at(xi_E_minus_Clh.at(i))*lh_square_arg.h_E.at(xi_E_minus_Clh.at(i));

}

}

//----------//

if (xi_I_Clh.size()>=1){
for (int i=0;i<= (int)(xi_I_Clh.size()-1);i++){
//lh_square_arg.f_I.at(xi_I_Clh.at(i)) = gsl_ran_gamma_pdf(t_i_Clh.at(xi_I_Clh.at(i)) - t_e_Clh.at(xi_I_Clh.at(i)), a_Clh, b_Clh);
lh_square_arg.f_I.at(xi_I_Clh.at(i)) = func_latent_pdf(t_i_Clh.at(xi_I_Clh.at(i)) - t_e_Clh.at(xi_I_Clh.at(i)), mu_lat_Clh,var_lat_Clh);

}
}

//--------//

// if (xi_R_Clh.size()>=1){
// for (int i=0;i<= (int)(xi_R_Clh.size()-1);i++){
// lh_square_arg.f_R.at(xi_R_Clh.at(i)) = gsl_ran_weibull_pdf(t_r_Clh.at(xi_R_Clh.at(i)) - t_i_Clh.at(xi_R_Clh.at(i)), c_Clh, d_Clh);
// }
// }

if (xi_R_Clh.size()>=1){
for (int i=0;i<= (int)(xi_R_Clh.size()-1);i++){
lh_square_arg.f_R.at(xi_R_Clh.at(i)) = gsl_ran_exponential_pdf(t_r_Clh.at(xi_R_Clh.at(i)) - t_i_Clh.at(xi_R_Clh.at(i)), c_Clh);
}
}

//-------//

if (xi_EnI_Clh.size()>=1){
for (int i=0;i<= (int)(xi_EnI_Clh.size()-1);i++){
//lh_square_arg.f_EnI.at(xi_EnI_Clh.at(i)) = 1.0 -  gsl_cdf_gamma_P(t_max_Clh - t_e_Clh.at(xi_EnI_Clh.at(i)), a_Clh, b_Clh);
lh_square_arg.f_EnI.at(xi_EnI_Clh.at(i)) = 1.0 -  func_latent_cdf(t_max_Clh - t_e_Clh.at(xi_EnI_Clh.at(i)),mu_lat_Clh,var_lat_Clh);
}
}
//-------//

// if (xi_InR_Clh.size()>=1){
// for (int i=0;i<= (int)(xi_InR_Clh.size()-1);i++){
// lh_square_arg.f_InR.at(xi_InR_Clh.at(i)) = 1.0 -  gsl_cdf_weibull_P(t_max_Clh - t_i_Clh.at(xi_InR_Clh.at(i)), c_Clh, d_Clh);
// }
// }

if (xi_InR_Clh.size()>=1){
for (int i=0;i<= (int)(xi_InR_Clh.size()-1);i++){
lh_square_arg.f_InR.at(xi_InR_Clh.at(i)) = 1.0 -  gsl_cdf_exponential_P(t_max_Clh - t_i_Clh.at(xi_InR_Clh.at(i)), c_Clh);
}
}


}


/*------------------------------------------------*/

void mcmc_UPDATE::mu_joint_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg,  const vector<int>& xi_E_arg, para_key& para_current_arg, nt_struct& nt_data_arg, int iter){

double mu_1_proposed = 0.0;
double mu_2_proposed = 0.0;

double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

const gsl_rng_type* T_c= gsl_rng_ranlux;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed
mu_1_proposed = para_current_arg.mu_1 + 0.001*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

switch (mu_1_proposed<=0) {
case 1: {
mu_1_proposed = -mu_1_proposed; //reflection
break;
}
case 0: {
mu_1_proposed = mu_1_proposed;
break;
}
}

mu_2_proposed = mu_1_proposed;
//---

for (int i=0;i<= (int)(xi_E_arg.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_arg.at(i);

	switch (nt_data_arg.current_size.at(k_E)>1) {

	case 1:{

	log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(k_E); //subtract part of likelihood that would be updated below

	lh_square_modified.log_f_S.at(k_E) = 0.0;

// 	vector<int> seq_1(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_CUPDATE);
// 	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+2*n_base_CUPDATE);

	for (int j=0;j<=(nt_data_arg.current_size.at(k_E)-2);j++){

// 		switch(j==0){
// 
// 		case 1:{
// 		lh_square_modified.f_S.at(k_E) =lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		case 0:{
// 		seq_1 = seq_2;
// 		vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);
// 		lh_square_modified.f_S.at(k_E)  = lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		}

	vector<int> seq_1(nt_data_arg.nt[k_E].begin()+ j*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE);
	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);

 	lh_square_modified.log_f_S.at(k_E)  = lh_square_modified.log_f_S.at(k_E)+ log_lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, mu_2_proposed, n_base_CUPDATE);

 	}


	log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(k_E); 

	break;
	}

	default:{
	break;
	}

	}
}

// boost::math::normal_distribution <double> mu_1_prior(0.04,0.002);
// double prior_ratio = pdf(mu_1_prior, mu_1_proposed)/pdf(mu_1_prior, para_current_arg.mu_1);

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("prior_mu_1.txt")).c_str(),ios::app);
// myfile_mcmc_out << mu_1_proposed << "," << para_current_arg.mu_1 <<","<< prior_ratio << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_1= mu_1_proposed;
para_current_arg.mu_2= mu_2_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_1 = para_current_arg.mu_1;
para_current_arg.mu_2 = para_current_arg.mu_2;

break;

}

gsl_rng_free(r_c);


}
/*------------------------------------------------*/

void mcmc_UPDATE::mu_1_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg,  const vector<int>& xi_E_arg, para_key& para_current_arg, nt_struct& nt_data_arg, int iter,gsl_rng* & r_c ){

double mu_1_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_ranlux;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

//mu_1_proposed = para_current_arg.mu_1 + 0.0005*gsl_ran_gaussian(r_c,1.0);
mu_1_proposed = para_current_arg.mu_1 + 0.00005*gsl_ran_gaussian(r_c,1.0);

// gsl_rng_free(r_c);

// switch (mu_1_proposed<=0) {
// case 1: {
// mu_1_proposed = -mu_1_proposed; //reflection
// break;
// }
// case 0: {
// mu_1_proposed = mu_1_proposed;
// break;
// }
// }


//switch (mu_1_proposed<=0) {
switch( (mu_1_proposed<=0) | (mu_1_proposed>=mu_1_hi) ){

case 1: {
mu_1_proposed = para_current_arg.mu_1;
break;
}
case 0: {
mu_1_proposed = mu_1_proposed;
break;
}
}


//----------

// double mu_1_up = 0.005; 
// switch ((mu_1_proposed<=0) | (mu_1_proposed>=mu_1_up)) {
// case 1: {
// 	  if (mu_1_proposed<=0) mu_1_proposed = -mu_1_proposed; //reflection
// 	  if (mu_1_proposed>=mu_1_up) mu_1_proposed = mu_1_up - (mu_1_proposed - mu_1_up); //reflection
// break;
// }
// case 0: {
// mu_1_proposed = mu_1_proposed;
// break;
// }
// }

//----------

//mu_1_proposed = gsl_ran_flat(r_c, 0.0001, 0.004); // use a uniform prior for proposal

//----------

for (int i=0;i<= (int)(xi_E_arg.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_arg.at(i);

	switch (nt_data_arg.current_size.at(k_E)>1) {

	case 1:{

	log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(k_E); //subtract part of likelihood that would be updated below

	lh_square_modified.log_f_S.at(k_E) = 0.0;

// 	vector<int> seq_1(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_CUPDATE);
// 	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+2*n_base_CUPDATE);

	for (int j=0;j<=(nt_data_arg.current_size.at(k_E)-2);j++){

// 		switch(j==0){
// 
// 		case 1:{
// 		lh_square_modified.f_S.at(k_E) =lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		case 0:{
// 		seq_1 = seq_2;
// 		vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);
// 		lh_square_modified.f_S.at(k_E)  = lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		}

	vector<int> seq_1(nt_data_arg.nt[k_E].begin()+ j*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE);
	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);

 	lh_square_modified.log_f_S.at(k_E)  = lh_square_modified.log_f_S.at(k_E)+ log_lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);

 	}


	log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(k_E); 

	break;
	}

	default:{
	break;
	}

	}
}


//boost::math::gamma_distribution <double> mu_1_prior(1,0.003);
// boost::math::normal_distribution <double> mu_1_prior(0.002,0.0005);
// 
// double prior_ratio = pdf(mu_1_prior, mu_1_proposed)/pdf(mu_1_prior, para_current_arg.mu_1);
// 
// acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*prior_ratio);

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("prior_mu_1.txt")).c_str(),ios::app);
// myfile_mcmc_out << mu_1_proposed << "," << para_current_arg.mu_1 <<","<< prior_ratio << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_1= mu_1_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_1 = para_current_arg.mu_1;
break;

}

//gsl_rng_free(r_c);


}
/*------------------------------------------------*/


void mcmc_UPDATE::mu_2_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg,  const vector<int>& xi_E_arg, para_key& para_current_arg, nt_struct& nt_data_arg, int iter, gsl_rng* & r_c){

double mu_2_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_ranlux;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

//mu_2_proposed = para_current_arg.mu_2 + 0.00005*gsl_ran_gaussian(r_c,1.0);
mu_2_proposed = para_current_arg.mu_2 + 0.000005*gsl_ran_gaussian(r_c,1.0);

// gsl_rng_free(r_c);

// switch (mu_2_proposed<=0) {
// case 1: {
// mu_2_proposed = -mu_2_proposed; //reflection
// break;
// }
// case 0: {
// mu_2_proposed = mu_2_proposed;
// break;
// }
// }

//switch (mu_2_proposed<=0) {
switch( (mu_2_proposed<=0) | (mu_2_proposed>=mu_2_hi) ){

case 1: {
mu_2_proposed = para_current_arg.mu_2;
break;
}
case 0: {
mu_2_proposed = mu_2_proposed;
break;
}
}


//------------------------

// double mu_2_up = 0.005; 
// switch ((mu_2_proposed<=0) | (mu_2_proposed>=mu_2_up)) {
// case 1: {
// 	  if (mu_2_proposed<=0) mu_2_proposed = -mu_2_proposed; //reflection
// 	  if (mu_2_proposed>=mu_2_up) mu_2_proposed = mu_2_up - (mu_2_proposed - mu_2_up); //reflection
// break;
// }
// case 0: {
// mu_2_proposed = mu_2_proposed;
// break;
// }
// }

//-------------------------

//mu_2_proposed = gsl_ran_flat(r_c, 0.0001, 0.002); // use a uniform prior for proposal

//------------------

for (int i=0;i<= (int)(xi_E_arg.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_arg.at(i);

	switch (nt_data_arg.current_size.at(k_E)>1) {

	case 1:{

	log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(k_E); //subtract part of likelihood that would be updated below

	lh_square_modified.log_f_S.at(k_E) = 0.0;

// 	vector<int> seq_1(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_CUPDATE);
// 	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+2*n_base_CUPDATE);

	for (int j=0;j<=(nt_data_arg.current_size.at(k_E)-2);j++){

// 		switch(j==0){
// 
// 		case 1:{
// 		lh_square_modified.f_S.at(k_E) =lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		case 0:{
// 		seq_1 = seq_2;
// 		vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);
// 		lh_square_modified.f_S.at(k_E)  = lh_square_modified.f_S.at(k_E)* lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], mu_1_proposed, para_current_arg.mu_2, n_base_CUPDATE);
// 		break;
// 		}
// 
// 		}

	vector<int> seq_1(nt_data_arg.nt[k_E].begin()+ j*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE);
	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*n_base_CUPDATE, nt_data_arg.nt[k_E].begin()+(j+2)*n_base_CUPDATE);

 	lh_square_modified.log_f_S.at(k_E)  = lh_square_modified.log_f_S.at(k_E)+ log_lh_seq(seq_1, seq_2, nt_data_arg.t_nt[k_E][j], nt_data_arg.t_nt[k_E][j+1], para_current_arg.mu_1, mu_2_proposed,  n_base_CUPDATE);

 	}


	log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(k_E); 

	break;
	}

	default:{
	break;
	}

	}
}


//boost::math::gamma_distribution <double> mu_2_prior(1.0,0.003);
// boost::math::normal_distribution <double> mu_2_prior(0.0005,0.00005);
// 
// double prior_ratio = pdf(mu_2_prior, mu_2_proposed)/pdf(mu_2_prior, para_current_arg.mu_2);
// 
// acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*prior_ratio);


acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));



double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// myfile_mcmc_out << prior_ratio<< endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_2= mu_2_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_2 = para_current_arg.mu_2;
break;

}

//gsl_rng_free(r_c);


}
/*------------------------------------------------*/


void mcmc_UPDATE::p_ber_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg,  const vector<int>& xi_E_arg, para_key& para_current_arg, nt_struct& nt_data_arg, vector<int>& infected_source_current_arg, vector<int>& con_seq, int iter, gsl_rng* & r_c){

long double p_ber_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;


p_ber_proposed = para_current_arg.p_ber + 0.01*gsl_ran_gaussian(r_c,1.0);

//p_ber_proposed =  gsl_ran_beta(r_c, 1.0,10.0);
//p_ber_proposed = gsl_ran_flat(r_c, 0.0, 1.0);

switch( (p_ber_proposed<=0) | (p_ber_proposed>=p_ber_hi) ){

case 1: {
p_ber_proposed = para_current_arg.p_ber;
break;
}
case 0: {
p_ber_proposed = p_ber_proposed;
break;
}
}

double log_prior_y =0.0;
double log_prior_x=0.0;

// log_prior_y = log(gsl_ran_beta_pdf(p_ber_proposed, 1.0, 10.0));
// log_prior_x = log(gsl_ran_beta_pdf(para_current_arg.p_ber, 1.0, 10.0));

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// myfile_mcmc_out << exp(log_prior_y - log_prior_x)<< endl;
// myfile_mcmc_out.close();


//------------------------

for (int i=0;i<= (int)(xi_E_arg.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_arg.at(i);

	switch (infected_source_current_arg.at(k_E)==9999) {

		case 1:{//bg infection
	
			log_lh_modified = log_lh_modified - lh_square_modified.log_f_Snull.at(k_E); //subtract part of likelihood that would be updated below
		
			vector<int> seq(nt_data_arg.nt[k_E].begin(), nt_data_arg.nt[k_E].begin()+n_base_CUPDATE); // the first seq
			lh_square_modified.log_f_Snull.at(k_E) = lh_snull(con_seq, seq, p_ber_proposed, n_base_CUPDATE); // compute the log pr a seq for background

			log_lh_modified = log_lh_modified + lh_square_modified.log_f_Snull.at(k_E); 
	
		break;
		}
	
		default:{ // 2nd infection
		break;
		}

	}
}



acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) +  (log_prior_y - log_prior_x)));



double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// myfile_mcmc_out << prior_ratio<< endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.p_ber= p_ber_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.p_ber = para_current_arg.p_ber;
break;

}

//gsl_rng_free(r_c);


}
/*------------------------------------------------*/

void mcmc_UPDATE::alpha_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const vector<double>& stb_arg, const vector<int>& infected_source_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

double alpha_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

alpha_proposed = para_current_arg.alpha + 0.001*gsl_ran_gaussian(r_c,1.0);

// double log_alpha_proposed;
// log_alpha_proposed = log(para_current_arg.alpha) + 0.008*gsl_ran_gaussian(r_c,1.0);
// alpha_proposed = exp(log_alpha_proposed );


// gsl_rng_free(r_c);

//double up = 0.01;

// switch (alpha_proposed<=0) {
// 	case 1: {
// 	alpha_proposed = -alpha_proposed; //reflection
// 	break;
// 	}
// 	case 0: {
// 	alpha_proposed = alpha_proposed;
// 	break;
// 	}
// }


//switch (alpha_proposed<=0) {
switch ( (alpha_proposed<=0) | (alpha_proposed>=alpha_hi)) {

	case 1: {
	alpha_proposed = para_current_arg.alpha; //rejection
	break;
	}
	case 0: {
	alpha_proposed = alpha_proposed;
	break;
	}
}

double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(alpha_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.alpha, 1.0/rate_exp_prior));


if (xi_U_arg.empty()==0){
for (int i=0; i<=(int)(xi_U_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.q_T.at(xi_U_arg.at(i)) = alpha_proposed*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_E_minus_arg.empty()==0){
for (int i=0; i<= (int)(xi_E_minus_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

// lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_proposed + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
// 

	switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
	case 9999:{ // by background
	lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) ; // this unchanged as long as infectious source unchanged
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = alpha_proposed;
	break;
	}

	// 	case -99:{
	// 	ofstream myfile_out; 
	// 	myfile_out.open((string(path4)+string("warnings.txt")).c_str(),ios::app);
	// 	myfile_out <<xi_E_minus_arg.at(i)<<endl;
	// 	myfile_out.close();	
	// 	break;
	// 	}

	default :{ // not by background
	lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) ; // this unchanged as long as infectious source unchanged
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) ;
	break;
	}

	}

lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = alpha_proposed*t_e_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))) ;
}
}

acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) + (log_prior_y-log_prior_x)));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.alpha/alpha_proposed));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.alpha = alpha_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.alpha = para_current_arg.alpha;
break;

}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::beta_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const vector<double>& stb_arg, const vector<int>& infected_source_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double beta_proposed = 0.0;
double acp_pr = 0.0;

// double log_prior_x =0.0;
// double log_prior_y =0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

beta_proposed = para_current_arg.beta + 1.0*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

// switch (beta_proposed<=0) {
// case 1: {
// beta_proposed = -beta_proposed; //reflection
// break;
// }
// case 0: {
// beta_proposed = beta_proposed;
// break;
// }
// }

//switch (beta_proposed<=0) {
switch ((beta_proposed<=0) | (beta_proposed>=beta_hi)) {

case 1: {
beta_proposed =para_current_arg.beta;
break;
}
case 0: {
beta_proposed = beta_proposed;
break;
}
}



if (xi_U_arg.empty()==0){
for (int i=0; i<=(int)(xi_U_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + beta_proposed*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_E_minus_arg.empty()==0){
for (int i=0; i<= (int)(xi_E_minus_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

	// lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + beta_proposed*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));

	switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
	case 9999:{ // by background
	//lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) ; // this unchanged as long as infectious source unchanged
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) ;
	break;
	}

	// 	case -99:{
	// 	ofstream myfile_out; 
	// 	myfile_out.open((string(path4)+string("warnings.txt")).c_str(),ios::app);
	// 	myfile_out <<xi_E_minus_arg.at(i)<<endl;
	// 	myfile_out.close();	
	// 	break;
	// 	}

	default :{ // not by background
	//lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) ; // this unchanged as long as infectious source unchanged
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =  beta_proposed*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
 ;
	break;
	}

	}


lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i)) + beta_proposed*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))) ;
}
}


// boost::math::gamma_distribution <double> beta_prior(10,1);
// double prior_ratio = pdf(beta_prior, beta_proposed)/pdf(beta_prior, para_current_arg.beta);
// acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg))*prior_ratio);


acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_beta.txt")).c_str(),ios::app);
// myfile_mcmc_out  << acp_pr << "," <<log_prior_y << "," << log_prior_x <<endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.beta = beta_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.beta = para_current_arg.beta;
break;

}

//gsl_rng_free(r_c);

}


/*------------------------------------------------*/
void mcmc_UPDATE::mu_exp_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double mu_lat_proposed = 0.0;
double acp_pr = 0.0;
double log_prior_y = 0.0;
double log_prior_x = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;


mu_lat_proposed = para_current_arg.mu_lat + 1.0*gsl_ran_gaussian(r_c,1.0);


//switch (mu_lat_proposed<=0) {
switch ((mu_lat_proposed<=0) | (mu_lat_proposed>=mu_lat_hi)) {

case 1: {
mu_lat_proposed = para_current_arg.mu_lat;
break;
}
case 0: {
mu_lat_proposed = mu_lat_proposed;
break;
}
}


if (xi_I_arg.empty()==0){
for (int i=0; i<=(int)(xi_I_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_I.at(xi_I_arg.at(i)) = func_latent_pdf(t_i_arg.at(xi_I_arg.at(i)) - t_e_arg.at(xi_I_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_EnI_arg.empty()==0){
for (int i=0; i<=(int)(xi_EnI_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_EnI.at(xi_EnI_arg.at(i)) =  1.0 - func_latent_cdf(t_max_CUPDATE - t_e_arg.at(xi_EnI_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //add back part of likelihood that updated above
}
}


//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x)));



double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_mu_lat.txt")).c_str(),ios::app);
// myfile_mcmc_out << mu_lat_proposed << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_lat = mu_lat_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_lat = para_current_arg.mu_lat;
break;

}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::mu_lat_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double mu_lat_proposed = 0.0;
double acp_pr = 0.0;


lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

mu_lat_proposed = para_current_arg.mu_lat + 1.0*gsl_ran_gaussian(r_c,1.0);
//mu_lat_proposed = para_current_arg.mu_lat + gsl_ran_flat(r_c,-0.05,0.05);

// double log_mu_lat_proposed;
// log_mu_lat_proposed = log(para_current_arg.mu_lat) + 0.1*gsl_ran_gaussian(r_c,1.0);
// mu_lat_proposed = exp(log_mu_lat_proposed );

// gsl_rng_free(r_c);

// switch (mu_lat_proposed<=0) {
// case 1: {
// mu_lat_proposed = -mu_lat_proposed; //reflection
// break;
// }
// case 0: {
// mu_lat_proposed = mu_lat_proposed;
// break;
// }
// }

//switch (mu_lat_proposed<=0) {
switch ((mu_lat_proposed<=0) | (mu_lat_proposed>=mu_lat_hi)) {

case 1: {
mu_lat_proposed = para_current_arg.mu_lat;
break;
}
case 0: {
mu_lat_proposed = mu_lat_proposed;
break;
}
}

double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(mu_lat_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.mu_lat, 1.0/rate_exp_prior));

// double mu_up = 10.0;
// 
// switch ((mu_lat_proposed<=0) | (mu_lat_proposed>=mu_up)) {
// 	case 1: {
// 	
// 		switch (mu_lat_proposed<=0) {
// 		case 1: {
// 			mu_lat_proposed = -mu_lat_proposed; //reflection
// 		break;
// 		}
// 		case 0: {
// 			mu_lat_proposed = mu_up - (mu_lat_proposed - mu_up); //reflection
// 		break;
// 		}
// 		}
// 	break;
// 	}
// 
// 	case 0: {
// 		mu_lat_proposed = mu_lat_proposed;
// 	break;
// 	}
// 	}


if (xi_I_arg.empty()==0){
for (int i=0; i<=(int)(xi_I_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_I.at(xi_I_arg.at(i)) = func_latent_pdf(t_i_arg.at(xi_I_arg.at(i)) - t_e_arg.at(xi_I_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_EnI_arg.empty()==0){
for (int i=0; i<=(int)(xi_EnI_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_EnI.at(xi_EnI_arg.at(i)) =  1.0 - func_latent_cdf(t_max_CUPDATE - t_e_arg.at(xi_EnI_arg.at(i)), mu_lat_proposed, para_current_arg.var_lat);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //add back part of likelihood that updated above
}
}


//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x)));



double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_mu_lat.txt")).c_str(),ios::app);
// myfile_mcmc_out << mu_lat_proposed << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.mu_lat = mu_lat_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.mu_lat = para_current_arg.mu_lat;
break;

}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::var_lat_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double var_lat_proposed = 0.0;
double acp_pr = 0.0;


lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

var_lat_proposed = para_current_arg.var_lat + 0.2*gsl_ran_gaussian(r_c,1.0);
//var_lat_proposed = para_current_arg.var_lat + gsl_ran_flat(r_c,-0.03,0.03);


// double log_var_lat_proposed;
// log_var_lat_proposed = log(para_current_arg.var_lat) + 0.1*gsl_ran_gaussian(r_c,1.0);
// var_lat_proposed = exp(log_var_lat_proposed);

// switch (var_lat_proposed<=0.0) {
// case 1: {
// var_lat_proposed = -var_lat_proposed; //reflection
// break;
// }
// case 0: {
// 	var_lat_proposed = var_lat_proposed;
// break;
// }
// }

//switch (var_lat_proposed<=0.1) {
switch ((var_lat_proposed<=var_lat_lo) | (var_lat_proposed>=var_lat_hi)) {

case 1: {
var_lat_proposed = para_current_arg.var_lat ;
break;
}
case 0: {
	var_lat_proposed = var_lat_proposed;
break;
}
}

double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(var_lat_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.var_lat, 1.0/rate_exp_prior));


// double var_up = 10.0;
// 
// switch ((var_lat_proposed<=0) | (var_lat_proposed>=var_up)) {
// 	case 1: {
// 	
// 		switch (var_lat_proposed<=0) {
// 		case 1: {
// 			var_lat_proposed = -var_lat_proposed; //reflection
// 		break;
// 		}
// 		case 0: {
// 			var_lat_proposed = var_up - (var_lat_proposed - var_up); //reflection
// 		break;
// 		}
// 		}
// 	break;
// 	}
// 
// 	case 0: {
// 		var_lat_proposed = var_lat_proposed;
// 	break;
// 	}
// 	}


if (xi_I_arg.empty()==0){
for (int i=0; i<=(int)(xi_I_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_I.at(xi_I_arg.at(i)) = func_latent_pdf(t_i_arg.at(xi_I_arg.at(i)) - t_e_arg.at(xi_I_arg.at(i)), para_current_arg.mu_lat, var_lat_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(xi_I_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_EnI_arg.empty()==0){
for (int i=0; i<=(int)(xi_EnI_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //subtract part of likelihood that would be updated below

lh_square_modified.f_EnI.at(xi_EnI_arg.at(i)) =  1.0 - func_latent_cdf(t_max_CUPDATE - t_e_arg.at(xi_EnI_arg.at(i)), para_current_arg.mu_lat, var_lat_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(xi_EnI_arg.at(i))); //add back part of likelihood that updated above
}
}

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_var_lat.txt")).c_str(),ios::app);
// myfile_mcmc_out << acp_pr << "," << log_prior_y << "," << log_prior_x <<endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch((uniform_rv<=acp_pr) & (isnan(log_lh_modified)==0)){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.var_lat = var_lat_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.var_lat = para_current_arg.var_lat;
break;

}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::c_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<int>& index_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double c_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

c_proposed = para_current_arg.c + 1.0*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

// switch (c_proposed<=0) {
// case 1: {
// c_proposed = -c_proposed; //reflection
// break;
// }
// case 0: {
// c_proposed = c_proposed;
// break;
// }
// }


//switch (c_proposed<=0) {
switch ((c_proposed<=0 ) | (c_proposed>=c_hi)) {

case 1: {
c_proposed =para_current_arg.c;
break;
}
case 0: {
c_proposed = c_proposed;
break;
}
}

double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(c_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.c, 1.0/rate_exp_prior));

if (xi_R_arg.empty()==0){
for (int i=0; i<=(int)(xi_R_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //subtract part of likelihood that would be updated below

//lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_weibull_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), c_proposed, para_current_arg.d);
lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_exponential_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), c_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_InR_arg.empty()==0){
for (int i=0; i<=(int)(xi_InR_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //subtract part of likelihood that would be updated below

//lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_weibull_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)), c_proposed, para_current_arg.d);
lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)), c_proposed);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //add back part of likelihood that updated above
}
}

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x)));

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.c = c_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.c = para_current_arg.c;
break;

}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::d_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<int>& index_arg, para_key& para_current_arg, int iter, gsl_rng* & r_c){

double d_proposed = 0.0;
double acp_pr = 0.0;


lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed

d_proposed = para_current_arg.d + 1.0*gsl_ran_gaussian(r_c,1.0);
// gsl_rng_free(r_c);

// switch (d_proposed<=0) {
// case 1: {
// d_proposed = -d_proposed; //reflection
// break;
// }
// case 0: {
// d_proposed = d_proposed;
// break;
// }
// }

//switch (d_proposed<=0) {
switch ((d_proposed<=0) | (d_proposed>=d_hi)) {

case 1: {
d_proposed = para_current_arg.d;
break;
}
case 0: {
d_proposed = d_proposed;
break;
}
}


double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(d_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.d, 1.0/rate_exp_prior));

if (xi_R_arg.empty()==0){
for (int i=0; i<=(int)(xi_R_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //subtract part of likelihood that would be updated below

//lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_weibull_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), para_current_arg.c, d_proposed );
lh_square_modified.f_R.at(xi_R_arg.at(i)) = gsl_ran_exponential_pdf(t_r_arg.at(xi_R_arg.at(i)) - t_i_arg.at(xi_R_arg.at(i)), para_current_arg.c);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(xi_R_arg.at(i))); //add back part of likelihood that updated above
}
}

if (xi_InR_arg.empty()==0){
for (int i=0; i<=(int)(xi_InR_arg.size()-1);i++){

log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //subtract part of likelihood that would be updated below

//lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_weibull_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)),  para_current_arg.c, d_proposed);
lh_square_modified.f_InR.at(xi_InR_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(t_max_CUPDATE - t_i_arg.at(xi_InR_arg.at(i)),  para_current_arg.c);


log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(xi_InR_arg.at(i))); //add back part of likelihood that updated above
}
}

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.d = d_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.d = para_current_arg.d;
break;
}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::k_1_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , const vector<int>& infected_source_arg, vector<double>& norm_const_current_arg, int iter, gsl_rng* & r_c){

double k_1_proposed = 0.0;
double acp_pr = 0.0;



lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > kernel_mat_modified = kernel_mat_current_arg;
vector<double> norm_const_modified = norm_const_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed


//k_1_proposed = para_current_arg.k_1 + 1.0*gsl_ran_flat(r_c,0.0,0.01);
k_1_proposed = para_current_arg.k_1 + 0.05*gsl_ran_gaussian(r_c,1.0);


// gsl_rng_free(r_c);

// switch (k_1_proposed<=0) {
// case 1: {
// k_1_proposed = -k_1_proposed;//reflection
// break;
// }
// case 0: {
// k_1_proposed = k_1_proposed;
// break;
// }
// }


//switch (k_1_proposed<=0) {
switch ((k_1_proposed<=0) | (k_1_proposed>=k_1_hi)) {

case 1: {
k_1_proposed = para_current_arg.k_1;
break;
}
case 0: {
k_1_proposed = k_1_proposed;
break;
}
}


double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(k_1_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.k_1, 1.0/rate_exp_prior));

//----------
for (int i=0;i<=(n_CUPDATE-1);i++) {
 for (int j=0;j<=(n_CUPDATE-1);j++) {
 if (i==j) kernel_mat_modified[i][j]=0.0;
 if (i<j) kernel_mat_modified[i][j] = func_kernel (coordinate_CUPDATE[i][0],coordinate_CUPDATE[i][1],coordinate_CUPDATE[j][0],coordinate_CUPDATE[j][1],k_1_proposed,para_current_arg.k_2,kernel_type_CUPDATE);
 if (i>j) kernel_mat_modified[i][j]=kernel_mat_modified[j][i];
 }
}

for (int j=0;j<=(n_CUPDATE-1);j++) {
norm_const_modified.at(j)=0.0;
 for (int i=0;(i<=(n_CUPDATE-1));i++) {
norm_const_modified.at(j)= norm_const_modified.at(j) + kernel_mat_modified[i][j];
}
}

//----------

if (xi_U_arg.empty()==0){
for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

    log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
    lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

//     double delta_t = 0.0;
//     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
//     case 1:{ // not yet recovered
//     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
//     break;
//     }
//     case 0:{ // recovered
//     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
//     break;
//     }    
//     }

    double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

    lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + delta_t*kernel_mat_modified[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j));

      
    }

lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));

lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
}
}

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){

    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   
        if (t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_minus_arg.at(i))) { 

	double delta_t = delta_mat_current_arg[xi_E_minus_arg.at(i)][xi_I_arg.at(j)];

//         switch (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i))) {
//         case 1:{ // not yet recovered at e_i
//         lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) + kernel_mat_modified[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j)); // update k_sum_E
//         break;
//         }
//         case 0:{ // recovered before e_i
//         break;
//         }    
//         }
// 	  

        lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) + delta_t*kernel_mat_modified[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j)); // update kt_sum_E

        } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 

	//lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));

		switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
		case 9999:{ // by background
		//lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)); // update k_sum_E
		lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));
		break;
		}
	
		// 	case -99:{
		// 	ofstream myfile_out; 
		// 	myfile_out.open((string(path4)+string("warnings.txt")).c_str(),ios::app);
		// 	myfile_out <<xi_E_minus_arg.at(i)<<endl;
		// 	myfile_out.close();	
		// 	break;
		// 	}
	
		default :{ // not by background
		lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = kernel_mat_modified[xi_E_minus_arg.at(i)][infected_source_arg.at(xi_E_minus_arg.at(i))]/norm_const_modified.at(infected_source_arg.at(xi_E_minus_arg.at(i))); // update k_sum_E
		lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
		break;
		}
	
		}

    } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}
//----------

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x) ));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_k_1_proposed.txt")).c_str(),ios::app);
// myfile_mcmc_out <<   acp_pr << "," <<log_prior_y << "," << log_prior_x <<endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();



switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
kernel_mat_current_arg = kernel_mat_modified;
norm_const_current_arg = norm_const_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.k_1 = k_1_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.k_1 = para_current_arg.k_1;
break;
}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::k_2_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, vector< vector<double> >& kernel_mat_current_arg,  const vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, para_key& para_current_arg, const vector <double>& stb_arg , const vector<int>& infected_source_arg, vector<double>& norm_const_current_arg, int iter, gsl_rng* & r_c){

double k_2_proposed = 0.0;
double acp_pr = 0.0;



lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > kernel_mat_modified = kernel_mat_current_arg;
vector<double> norm_const_modified = norm_const_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c, 1000*iter); // set a seed


k_2_proposed = para_current_arg.k_2 + 1.0*gsl_ran_gaussian(r_c,1.0);

double log_prior_y =0.0;
double log_prior_x=0.0;

log_prior_y = log(gsl_ran_exponential_pdf(k_2_proposed, 1.0/rate_exp_prior));
log_prior_x = log(gsl_ran_exponential_pdf(para_current_arg.k_2, 1.0/rate_exp_prior));

// gsl_rng_free(r_c);

// switch (k_1_proposed<=0) {
// case 1: {
// k_1_proposed = -k_1_proposed;//reflection
// break;
// }
// case 0: {
// k_1_proposed = k_1_proposed;
// break;
// }
// }


//switch (k_1_proposed<=0) {
switch ((k_2_proposed<=0) | (k_2_proposed>=k_2_hi)) {

case 1: {
k_2_proposed = para_current_arg.k_2;
break;
}
case 0: {
k_2_proposed = k_2_proposed;
break;
}
}



// log_prior_x = log(dtnorm(para_current_arg.k_1, k_1_p_1, k_1_p_2, 0.0));
// log_prior_y = log(dtnorm(k_1_proposed, k_1_p_1, k_1_p_2, 0.0));

//----------
for (int i=0;i<=(n_CUPDATE-1);i++) {
 for (int j=0;j<=(n_CUPDATE-1);j++) {
 if (i==j) kernel_mat_modified[i][j]=0.0;
 if (i<j) kernel_mat_modified[i][j] = func_kernel (coordinate_CUPDATE[i][0],coordinate_CUPDATE[i][1],coordinate_CUPDATE[j][0],coordinate_CUPDATE[j][1],para_current_arg.k_1,k_2_proposed,kernel_type_CUPDATE);
 if (i>j) kernel_mat_modified[i][j]=kernel_mat_modified[j][i];
 }
}

for (int j=0;j<=(n_CUPDATE-1);j++) {
norm_const_modified.at(j)=0.0;
 for (int i=0;(i<=(n_CUPDATE-1));i++) {
norm_const_modified.at(j)= norm_const_modified.at(j) + kernel_mat_modified[i][j];
}
}

//----------

if (xi_U_arg.empty()==0){
for (int i=0;i<= (int)(xi_U_arg.size()-1);i++){

    log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
    
    lh_square_modified.kt_sum_U.at(xi_U_arg.at(i))  = 0.0;
    
    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){

//     double delta_t = 0.0;
//     switch (t_r_arg.at(xi_I_arg.at(j))>t_max_CUPDATE) {
//     case 1:{ // not yet recovered
//     delta_t = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
//     break;
//     }
//     case 0:{ // recovered
//     delta_t = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
//     break;
//     }    
//     }

    double delta_t = delta_mat_current_arg[xi_U_arg.at(i)][xi_I_arg.at(j)];

    lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(i)) + delta_t*kernel_mat_modified[xi_U_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j));

      
    }

lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(xi_U_arg.at(i))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));

lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
}
}

//----------

if (xi_E_minus_arg.empty()==0){
for (int i=0;i<= (int)(xi_E_minus_arg.size()-1);i++){

    log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below

    lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = 0.0;
    lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) =0.0;
    
    for (int j=0;j<= (int) (xi_I_arg.size()-1);j++){
   
        if (t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_minus_arg.at(i))) { 

	double delta_t = delta_mat_current_arg[xi_E_minus_arg.at(i)][xi_I_arg.at(j)];

//         switch (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i))) {
//         case 1:{ // not yet recovered at e_i
//         lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) + kernel_mat_modified[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j)); // update k_sum_E
//         break;
//         }
//         case 0:{ // recovered before e_i
//         break;
//         }    
//         }
// 	  

        lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i)) + delta_t*kernel_mat_modified[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/norm_const_modified.at(xi_I_arg.at(j)); // update kt_sum_E

        } // end of if (t_i_Clh.at(xi_I_Clh.at(j))<t_e_Clh.at(xi_E_Clh.at(i))) 

	//lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));

		switch(infected_source_arg.at(xi_E_minus_arg.at(i))){
		case 9999:{ // by background
		//lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)); // update k_sum_E
		lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i));
		break;
		}
	
		// 	case -99:{
		// 	ofstream myfile_out; 
		// 	myfile_out.open((string(path4)+string("warnings.txt")).c_str(),ios::app);
		// 	myfile_out <<xi_E_minus_arg.at(i)<<endl;
		// 	myfile_out.close();	
		// 	break;
		// 	}
	
		default :{ // not by background
		lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i)) = kernel_mat_modified[xi_E_minus_arg.at(i)][infected_source_arg.at(xi_E_minus_arg.at(i))]/norm_const_modified.at(infected_source_arg.at(xi_E_minus_arg.at(i))); // update k_sum_E
		lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) =para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
		break;
		}
	
		}

    } // end of  for (int j=0;j<= (int) (xi_I_Clh.size()-1);j++)


lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);

lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));

log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
}
}
//----------

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_prior_y - log_prior_x) ));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("acp_k_1_proposed.txt")).c_str(),ios::app);
// myfile_mcmc_out <<   acp_pr << "," <<log_prior_y << "," << log_prior_x <<endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();



switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
kernel_mat_current_arg = kernel_mat_modified;
norm_const_current_arg = norm_const_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.k_2 = k_2_proposed;
break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.k_2 = para_current_arg.k_2;
break;
}

//gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::stb_1_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, vector<double>& stb_arg, para_key& para_current_arg, const vector<int>& gp_stb_arg, int iter){

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

double stb_1_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

double step_stb_1 = 0.01;

stb_1_proposed = para_current_arg.stb_1 + step_stb_1*gsl_ran_gaussian(r_c,1.0);


switch (stb_1_proposed<=0) {
case 1: {
stb_1_proposed = -stb_1_proposed; //reflection
break;
}
case 0: {
stb_1_proposed = stb_1_proposed;
break;
}
}

if (xi_U_arg.empty()==0){
for (int i=0; i<=(int)(xi_U_arg.size()-1);i++){

	if (gp_stb_arg.at(xi_U_arg.at(i))==1) {
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
	
	lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_1_proposed*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
	lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
	}
}
}

if (xi_E_minus_arg.empty()==0){
for (int i=0; i<= (int)(xi_E_minus_arg.size()-1);i++){

	if (gp_stb_arg.at(xi_E_minus_arg.at(i))==1) {
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
	
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + para_current_arg.beta*stb_1_proposed*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
	
	lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) =  para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.beta*stb_1_proposed*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
	lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);
	
	lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))) ;
	}
}
}

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.alpha/alpha_proposed));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.stb_1 = stb_1_proposed;

for (int i=0; i<=((int) stb_arg.size()-1); i++){
	switch(gp_stb_arg.at(i)){
	case 1:{
	stb_arg.at(i) = stb_1_proposed;
	break;
	}
	default :{
	stb_arg.at(i) = stb_arg.at(i);
	}
	}
}

break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.stb_1 = para_current_arg.stb_1;
break;

}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::stb_2_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_minus_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, vector<double>& stb_arg, para_key& para_current_arg, const vector<int>& gp_stb_arg, int iter){
const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c, 1000*iter); // set a seed

double stb_2_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

double step_stb_2 = 0.01;

stb_2_proposed = para_current_arg.stb_2 + step_stb_2*gsl_ran_gaussian(r_c,1.0);


switch (stb_2_proposed<=0) {
case 1: {
stb_2_proposed = -stb_2_proposed; //reflection
break;
}
case 0: {
stb_2_proposed = stb_2_proposed;
break;
}
}

if (xi_U_arg.empty()==0){
for (int i=0; i<=(int)(xi_U_arg.size()-1);i++){

	if (gp_stb_arg.at(xi_U_arg.at(i))==2) {
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //subtract part of likelihood that would be updated below
	
	lh_square_modified.q_T.at(xi_U_arg.at(i)) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_2_proposed*lh_square_modified.kt_sum_U.at(xi_U_arg.at(i));
	lh_square_modified.f_U.at(xi_U_arg.at(i)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(i)),1.0);
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(i))); //add back part of likelihood that updated above
	}
}
}

if (xi_E_minus_arg.empty()==0){
for (int i=0; i<= (int)(xi_E_minus_arg.size()-1);i++){

	if (gp_stb_arg.at(xi_E_minus_arg.at(i))==2) {
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))); //subtract part of likelihood that would be updated below
	
	lh_square_modified.g_E.at(xi_E_minus_arg.at(i)) = para_current_arg.alpha + para_current_arg.beta*stb_2_proposed*lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(i));
	
	lh_square_modified.q_E.at(xi_E_minus_arg.at(i)) =  para_current_arg.alpha*t_e_arg.at(xi_E_minus_arg.at(i)) + para_current_arg.beta*stb_2_proposed*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(i));
	lh_square_modified.h_E.at(xi_E_minus_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(i)),1.0);
	
	lh_square_modified.f_E.at(xi_E_minus_arg.at(i)) = lh_square_modified.g_E.at(xi_E_minus_arg.at(i))*lh_square_modified.h_E.at(xi_E_minus_arg.at(i));
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(i))) ;
	}
}
}

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(para_current_arg.alpha/alpha_proposed));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("uniform_rv_acp.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("lh_change_.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_modified << "," << log_lh_current_arg << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
para_current_arg.stb_2 = stb_2_proposed;

for (int i=0; i<=((int) stb_arg.size()-1); i++){
	switch(gp_stb_arg.at(i)){
	case 2:{
	stb_arg.at(i) = stb_2_proposed;
	break;
	}
	default :{
	stb_arg.at(i) = stb_arg.at(i);
	}
	}
}

break;
}

case 0: {
}
lh_square_current_arg = lh_square_current_arg;
log_lh_current_arg = log_lh_current_arg;
para_current_arg.stb_2 = para_current_arg.stb_2;
break;

}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/
void mcmc_UPDATE::con_seq_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vector<int>&  xi_beta_E_arg, vector<int>& con_seq,  int iter, gsl_rng* & r_c){

double acp_pr = 0.0;

int position_proposed, base_proposed;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

position_proposed =iter;

//---

base_proposed=0;
int base_current = con_seq.at(position_proposed); 

switch(base_current){
	case 1:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 2;
			break;
			}
			case 1:{
				base_proposed = 3;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}
	break;
	}
	case 2:{
		int type = gsl_rng_uniform_int (r_c, 3);


		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 3;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}	
	break;
	}
	case 3:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 2;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}	
	break;
	}
	case 4:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 2;
			break;
			}
			case 2:{
				base_proposed = 3;
			break;
			}
		}	
	break;
	}
}
//---


for (int i=0;i<= (int)(xi_E_arg.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_arg.at(i);

	switch (infected_source_current_arg.at(k_E)==9999) {

		case 1:{//bg infection
		
			int base_k_E = nt_current_arg[k_E][position_proposed];

			log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(k_E); 
		
			double log_x =  lh_snull_base(con_seq.at(position_proposed),base_k_E , para_current_arg.p_ber);
			double log_y =  lh_snull_base(base_proposed, base_k_E, para_current_arg.p_ber);
		
			lh_square_modified.log_f_Snull.at(k_E) =  lh_square_modified.log_f_Snull.at(k_E)  - log_x + log_y;
		
			log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(k_E); 
			
		break;
		}
	
		default:{ // 2nd infection
		break;
		}

	}
}

//----

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("000.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_lh_current_arg<< "," << log_lh_modified<< endl;
// myfile_mcmc_out.close();


	switch(uniform_rv<=acp_pr){
	case 1: {
	
	lh_square_current_arg = lh_square_modified;
	log_lh_current_arg = log_lh_modified;
	con_seq.at(position_proposed) = base_proposed;	
	
	break;
	}
	
	case 0: {
	break;
	}
	}

}


/*------------------------------------------------*/

void mcmc_UPDATE::seq_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vector<int>&  xi_beta_E_arg, vector<int>& con_seq, const int& subject_proposed,  int iter, gsl_rng* & r_c){

//int subject_proposed ;
double acp_pr = 0.0;

double log_part_x_subject=0.0;
double log_part_y_subject=0.0;

double log_part_x_source=0.0;
double log_part_y_source =0.0;


int position_proposed, base_proposed;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

//vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
//vector<double> t_e_modified = t_e_arg;
//vector <int> index_modified = index_arg;
//vector <int> xi_E_minus_modified = xi_E_minus_arg;

//vector <int> xi_U_modified = xi_U_arg;
//vector <int> xi_E_modified = xi_E_arg;
//vector <int> xi_EnI_modified = xi_EnI_arg;

// const gsl_rng_type* T_c= gsl_rng_ranlux;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,(iter+1)*time(NULL)); // set a seed



//subject_proposed = xi_E_arg.at(gsl_rng_uniform_int (r_c, xi_E_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn
//subject_proposed = xi_E_arg.at(gsl_rng_uniform_int (r_c, xi_beta_E_arg.size())); // only change the one from secondary infection 

// //----- test---------//
// while((current_size_arg.at(subject_proposed)>1)| (infected_source_current_arg.at(subject_proposed)==9999)){
// subject_proposed = xi_E_arg.at(gsl_rng_uniform_int (r_c, xi_E_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn
// }
// //--------------------//

//subject_proposed =3;

//position_proposed = gsl_rng_uniform_int (r_c, n_base_CUPDATE);
position_proposed =iter;


//vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);

//vector<int> seq_proposed; 
//seq_proposed.assign(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE );


int subject_source = infected_source_current_arg.at(subject_proposed);


//---

//base_proposed = gsl_rng_uniform_int (r_c, 4) + 1; // a new base proposed

base_proposed=0;
int base_current = nt_current_arg[subject_proposed][position_proposed]; // always refers to the first sequence

switch(base_current){
	case 1:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 2;
			break;
			}
			case 1:{
				base_proposed = 3;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}
	break;
	}
	case 2:{
		int type = gsl_rng_uniform_int (r_c, 3);

// 		ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("type_test.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  type << endl;
// 		myfile_mcmc_out.close();


		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 3;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}	
	break;
	}
	case 3:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 2;
			break;
			}
			case 2:{
				base_proposed = 4;
			break;
			}
		}	
	break;
	}
	case 4:{
		int type = gsl_rng_uniform_int (r_c, 3);
		switch(type){
			case 0:{
				base_proposed = 1;
			break;
			}
			case 1:{
				base_proposed = 2;
			break;
			}
			case 2:{
				base_proposed = 3;
			break;
			}
		}	
	break;
	}
}
//---

int base_next_subject =0;

switch (current_size_arg.at(subject_proposed)>1) {

	case 1:{
		//--
		base_next_subject = nt_current_arg[subject_proposed][position_proposed + n_base_CUPDATE];

		log_part_x_subject = log_lh_base(base_current, base_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2);

		log_part_y_subject = log_lh_base(base_proposed, base_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2);

		//--
		lh_square_modified.log_f_S.at(subject_proposed) =  lh_square_modified.log_f_S.at(subject_proposed) - log_part_x_subject;
		log_lh_modified = log_lh_modified - log_part_x_subject;
		
		lh_square_modified.log_f_S.at(subject_proposed) =  lh_square_modified.log_f_S.at(subject_proposed) + log_part_y_subject;
		log_lh_modified = log_lh_modified + log_part_y_subject;


	break;
	}

	case 0:{
	break;
	}
}
//-------

int rank_source =-1;  //count the rank of the original t_e among t_nt_current_arg.at(subject_source) 
int base_before_source =0;
int base_next_source = 0;

switch(subject_source ){

case 9999:{ // by background
//t_proposed= gsl_ran_flat(r_c, 0.0, min( t_sample_arg.at(subject_proposed), min( t_i_arg.at(subject_proposed), t_max_CUPDATE)) );

	log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 

	double log_x =  lh_snull_base(con_seq.at(position_proposed), base_current, para_current_arg.p_ber);
	double log_y =  lh_snull_base(con_seq.at(position_proposed), base_proposed, para_current_arg.p_ber);

 	lh_square_modified.log_f_Snull.at(subject_proposed) =  lh_square_modified.log_f_Snull.at(subject_proposed)  - log_x + log_y;

	log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); 

break;
}

default :{ // not by background

//nt_modified_source = nt_current_arg.at(subject_source);
rank_source = distance( t_nt_current_arg.at(subject_source).begin(), find(t_nt_current_arg.at(subject_source).begin(), t_nt_current_arg.at(subject_source).end(), t_e_arg.at(subject_proposed)) );
//nt_modified_source.erase(nt_modified_source.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source

base_before_source =  nt_current_arg[subject_source][(rank_source-1)*n_base_CUPDATE + position_proposed];
 
	switch(current_size_arg.at(subject_source)>(rank_source+1)){
		case 1:{// there  is a valid base_next_source
			base_next_source =  nt_current_arg[subject_source][(rank_source+1)*n_base_CUPDATE + position_proposed];


			log_part_x_source = log_lh_base(base_before_source, base_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2) + log_lh_base(base_current, base_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2);

			log_part_y_source = log_lh_base(base_before_source, base_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2) + log_lh_base(base_proposed, base_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2);
	
	

		break;
		}
		
		case 0:{

			log_part_x_source = log_lh_base(base_before_source, base_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2);

			log_part_y_source = log_lh_base(base_before_source, base_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2);
	
			break;
		}
	}

lh_square_modified.log_f_S.at(subject_source) =  lh_square_modified.log_f_S.at(subject_source) - log_part_x_source;
log_lh_modified = log_lh_modified - log_part_x_source;


lh_square_modified.log_f_S.at(subject_source) =  lh_square_modified.log_f_S.at(subject_source) + log_part_y_source;
log_lh_modified = log_lh_modified + log_part_y_source;

break;
}

}

//------

// double log_part_y = log_part_y_subject + log_part_y_source;
// double log_part_x = log_part_x_subject + log_part_x_source;
//acp_pr = min(1.0,exp(log_part_y-log_part_x));

acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("22.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();


switch(uniform_rv<=acp_pr){
case 1: {

lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;
nt_current_arg[subject_proposed][position_proposed]= base_proposed;

	switch (subject_source){
		
		case 9999:{ // by background
		break;
		}
		
		default :{ // not by background
			nt_current_arg[subject_source][(rank_source)*n_base_CUPDATE + position_proposed] = base_proposed;
		}
	}
			

break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);

//--------------
		


		
		// myfile_mcmc_out.open((string(path4)+string("index_modified.txt")).c_str(),ios::app);
		// if (index_modified.empty()==0){
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out <<  index_modified.at(i) << "," <<  t_e_modified.at(index_modified.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg.txt")).c_str(),ios::app);
		// if (index_arg.empty()==0){
		// for (int i=0; i<=(int)(index_arg.size()-1); i++){
		// myfile_mcmc_out << index_arg.at(i) << "," << t_e_arg.at(index_arg.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("find_xi_E_mibus.txt")).c_str(),ios::app);
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out << 1*(find(xi_E_minus_modified.begin(), xi_E_minus_modified.end(),index_modified.at(i))==xi_E_minus_modified.end()) << endl; // should always equals to 1, i.e., new index has been excluded from xi_E_minus
		// }
		// myfile_mcmc_out.close();
		// 

		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_kt_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.kt_sum_E.at(index_arg.at(0)) << endl; // shouls always be 0
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_k_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.k_sum_E.at(index_arg.at(0)) << endl;  // shouls always be 0
		// myfile_mcmc_out.close();

// 		ofstream myfile_mcmc_out; 
// 
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_seq.txt")).c_str(),ios::app);
// 		//if (log_lh_current_arg==log_lh_modified) 
// 		//myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << "," << acp_pr << endl;
// 		myfile_mcmc_out <<  log_part_x  << "," <<log_part_y << "," << acp_pr << endl;
// 		myfile_mcmc_out.close();
// 
//  
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed.txt")).c_str(),ios::app); 		
// 		if ((log_lh_current_arg==log_lh_modified)) myfile_mcmc_out << subject_proposed  << ","<< position_proposed << ","  << base_current << "," << base_proposed << "," << acp_pr<<endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("uniform_rv.txt")).c_str(),ios::app);
// 		if (uniform_rv<=acp_pr) myfile_mcmc_out << uniform_rv << "," << acp_pr  << endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("acp_pr.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << min(1.0,exp(log_part_y-log_part_x)) << "," <<min(1.0,exp(log_lh_modified-log_lh_current_arg))  << endl;
// 		myfile_mcmc_out.close();
		//---------------

}


/*------------------------------------------------*/


void mcmc_UPDATE::t_i_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, const vector<int>& xi_U_arg, const vector<int>& xi_E_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, const vector<int>& xi_EnI_arg, const vector<int>& xi_R_arg, const vector<int>& xi_InR_arg, const vector<double>& t_r_arg, vector<double>& t_i_arg, const vector<double>& t_onset_arg, const vector<double>& t_e_arg, const vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg,  vec2int& infecting_list_current_arg, const vector<int>& infecting_size_current_arg, int iter){

double t_proposed; // new t_i to be proposed
double t_low, t_up;

double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_i_modified = t_i_arg;


const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter); // set a seed

int subject_proposed = xi_I_arg.at(gsl_rng_uniform_int (r_c, xi_I_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn
double t_o = t_onset_arg.at(subject_proposed); // this is assumed to be given,e.g, the onset time  (in simulation, this may be taken to be the true t_i)
//double t_range = 1.0; // prior belief; we assume the sampled t_i would range between [t_o - t_range, t_o + range]


switch(infecting_size_current_arg.at(subject_proposed)>=1){
	case 1:{
		double min_t = t_e_arg.at(infecting_list_current_arg[subject_proposed][0]); // the minimum t_e among all those exposures being infected by the subject; as infecting_list is sorted according to the order of infecting, the first one is the minimum

		t_low = max(t_o - t_range, t_e_arg.at(subject_proposed));
		t_up = min(t_o + t_range, min_t);

// 		t_low = t_e_arg.at(subject_proposed);
// 		t_up = min_t;

	break;
	}
	case 0:{

		t_low = max(t_o - t_range, t_e_arg.at(subject_proposed));
		t_up = min( t_max_CUPDATE, min(t_o + t_range, t_r_arg.at(subject_proposed)));

// 		t_low = t_e_arg.at(subject_proposed);
// 		t_up =  min(t_max_CUPDATE, t_r_arg.at(subject_proposed));

	break;
	}
}

t_proposed= gsl_ran_flat(r_c,t_low, t_up );


// t_proposed = t_i_arg.at(subject_proposed)+ 1.0*gsl_ran_gaussian(r_c,1.0); //random walk//
// 
// switch ((t_proposed<t_low) | (t_proposed>t_up)) {
// 
// case 1: {
// t_proposed = t_i_arg.at(subject_proposed);
// break;
// }
// case 0: {
// t_proposed = t_proposed;
// break;
// }
// }



t_i_modified.at(subject_proposed) = t_proposed;


double log_prior_y =0.0;
double log_prior_x=0.0;

// log_prior_y = log(gsl_ran_gamma_pdf(t_proposed, t_o/pow(0.5,2.0), pow(0.5,2.0)));
// log_prior_x = log(gsl_ran_gamma_pdf(t_i_arg.at(subject_proposed), t_o/pow(0.5,2.0), pow(0.5,2.0)));

//----------------------------------------------------------------------------------//

for (int j=0;j<=(int) (xi_U_arg.size()-1);j++){ 

log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(xi_U_arg.at(j))); //subtract part of likelihood that would be updated below


switch (t_r_arg.at(subject_proposed)>=t_max_CUPDATE) {
	case 1:{
		delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_max_CUPDATE  - t_proposed;
	break;
	}
	case 0:{
		delta_mat_modified[xi_U_arg.at(j)][subject_proposed] = t_r_arg.at(subject_proposed) - t_proposed;
	break;
	}
}	


lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) = lh_square_modified.kt_sum_U.at(xi_U_arg.at(j)) - delta_mat_current_arg[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) + delta_mat_modified[xi_U_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_U_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed)  ; //subtract the infectious challenge  due to the infectious subject chosen THEN add the updated one back

lh_square_modified.q_T.at(xi_U_arg.at(j)) = para_current_arg.alpha*stb_arg.at(xi_U_arg.at(j))* t_max_CUPDATE  + para_current_arg.beta*stb_arg.at(xi_U_arg.at(j))*lh_square_modified.kt_sum_U.at(xi_U_arg.at(j));

lh_square_modified.f_U.at(xi_U_arg.at(j)) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(xi_U_arg.at(j)),1.0);

log_lh_modified = log_lh_modified + log(lh_square_modified.f_U.at(xi_U_arg.at(j))); // add back the updated part of likelihood

}

//----------------------------------------------------------------------------------//

for (int j=0;j<=(int) (xi_E_minus_arg.size()-1);j++){ 
	
	if (xi_E_minus_arg.at(j)!=subject_proposed) {
	
		switch ((t_e_arg.at(xi_E_minus_arg.at(j))>t_proposed) ) {
		
			case 1:{
			
			log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(j))); //subtract part of likelihood that would be updated below

					switch (t_r_arg.at(subject_proposed)>=t_e_arg.at(xi_E_minus_arg.at(j))) {
						case 1:{
							delta_mat_modified[xi_E_minus_arg.at(j)][subject_proposed] = t_e_arg.at(xi_E_minus_arg.at(j)) - t_proposed;
						break;
						}
						case 0:{
							delta_mat_modified[xi_E_minus_arg.at(j)][subject_proposed] = t_r_arg.at(subject_proposed) - t_proposed;
						break;
						}
					}	
	
					
					switch((t_e_arg.at(xi_E_minus_arg.at(j))>t_i_arg.at(subject_proposed))){
					case 1:{ 
					//lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j))=lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j));
					lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j))= lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j)) - delta_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed)  + delta_mat_modified[xi_E_minus_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed)  ; 
					break;
					}
					case 0:{
					//lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j)) = lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j)) + kernel_mat_current_arg[xi_E_minus_modified.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) ;
					lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j))= lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j)) + delta_mat_modified[xi_E_minus_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) ; 
					break;
					}
					}
				
			lh_square_modified.q_E.at(xi_E_minus_arg.at(j))= para_current_arg.alpha*stb_arg.at(xi_E_minus_arg.at(j))*t_e_arg.at(xi_E_minus_arg.at(j)) + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(j))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j));
			lh_square_modified.g_E.at(xi_E_minus_arg.at(j))=   lh_square_modified.g_E.at(xi_E_minus_arg.at(j)); // unchanged as source does not change
			
			lh_square_modified.h_E.at(xi_E_minus_arg.at(j))= gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(j)),1.0);
			
			lh_square_modified.f_E.at(xi_E_minus_arg.at(j))= lh_square_modified.g_E.at(xi_E_minus_arg.at(j))*lh_square_modified.h_E.at(xi_E_minus_arg.at(j));
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(j))); // add back the updated part of likelihood
			
			break;
			
			}
		
			case 0:{
			
			log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(xi_E_minus_arg.at(j))); //subtract part of likelihood that would be updated below
			
			delta_mat_modified[xi_E_minus_arg.at(j)][subject_proposed] = 0.0;
			
					switch((t_e_arg.at(xi_E_minus_arg.at(j))>t_i_arg.at(subject_proposed))){
					case 1:{ 
					//lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j))= lh_square_modified.k_sum_E.at(xi_E_minus_arg.at(j)) - kernel_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) ;
					lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j))= lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j)) - delta_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]*kernel_mat_current_arg[xi_E_minus_arg.at(j)][subject_proposed]/norm_const_current_arg.at(subject_proposed) ;
					break;
					}
					case 0:{
					//lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j)) = lh_square_modified.k_sum_E.at(xi_E_minus_modified.at(j));
					lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j))= lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j));
					}
					}
					
			lh_square_modified.q_E.at(xi_E_minus_arg.at(j))= para_current_arg.alpha*stb_arg.at(xi_E_minus_arg.at(j))*t_e_arg.at(xi_E_minus_arg.at(j)) + para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(j))*lh_square_modified.kt_sum_E.at(xi_E_minus_arg.at(j));
			lh_square_modified.g_E.at(xi_E_minus_arg.at(j))=   lh_square_modified.g_E.at(xi_E_minus_arg.at(j));
			
			lh_square_modified.h_E.at(xi_E_minus_arg.at(j))= gsl_ran_exponential_pdf(lh_square_modified.q_E.at(xi_E_minus_arg.at(j)),1.0);
			
			lh_square_modified.f_E.at(xi_E_minus_arg.at(j))= lh_square_modified.g_E.at(xi_E_minus_arg.at(j))*lh_square_modified.h_E.at(xi_E_minus_arg.at(j));
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(xi_E_minus_arg.at(j))); // add back the updated part of likelihood
			
			break;
			
			}
		
		}
	
	}
}

//----------------------------------------------------------------------------------//

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_proposed - t_e_arg.at(subject_proposed), para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 


//----------

switch ( find(xi_R_arg.begin(), xi_R_arg.end(),subject_proposed) != (xi_R_arg.end()) ) { //return 1 when the subject is also in xi_R
case 1:{
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_R.at(subject_proposed)); //subtract part of likelihood that would be updated below
//	lh_square_modified.f_R.at(subject_proposed) = gsl_ran_weibull_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c, para_current_arg.d);
	lh_square_modified.f_R.at(subject_proposed) = gsl_ran_exponential_pdf(t_r_arg.at(subject_proposed) - t_proposed, para_current_arg.c);
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_R.at(subject_proposed)); 
	break;
	}
	case 0:{
	break;
	}
}

//----------
switch ( find(xi_InR_arg.begin(), xi_InR_arg.end(),subject_proposed) != (xi_InR_arg.end()) ) { //return 1 when the subject is also in xi_InR
	case 1:{
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_InR.at(subject_proposed)); //subtract part of likelihood that would be updated below
	//lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_weibull_P(t_max_CUPDATE - t_proposed, para_current_arg.c, para_current_arg.d);
	lh_square_modified.f_InR.at(subject_proposed) = 1.0 -  gsl_cdf_exponential_P(t_max_CUPDATE - t_proposed, para_current_arg.c);
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_InR.at(subject_proposed)); 
	
	break;
	}
	case 0:{
	break;
	}
}

//---------------
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg) +(log_prior_y-log_prior_x)));

// 		ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << acp_pr  << "," << t_i_arg.at(subject_proposed) << "," << t_proposed << endl;
// 		myfile_mcmc_out.close();

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
	case 1: {
		lh_square_current_arg = lh_square_modified;
		delta_mat_current_arg = delta_mat_modified;
		log_lh_current_arg = log_lh_modified;
		t_i_arg= t_i_modified;
	break;
	}
	
	case 0: {
	break;
	}
}

gsl_rng_free(r_c);

}

/*------------------------------------------------*/


void mcmc_UPDATE::t_e_seq(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg,  vec2int& infecting_list_current_arg, const vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, vector<int>& con_seq, int& subject_proposed, int iter,gsl_rng* & r_c){

//double t_back =10.0;

//int subject_proposed ;
double t_proposed = 0.0;
double t_low, t_up;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

// vector <int> xi_U_modified = xi_U_arg;
// vector <int> xi_E_modified = xi_E_arg;
// vector <int> xi_EnI_modified = xi_EnI_arg;

/*const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter); // set a see*/



vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

vector<int> nt_current_seq; // the orginal sequence of the subject which would be updated

int subject_source = infected_source_current_arg.at(subject_proposed);


//int rank_subject_x =distance( t_nt_current_arg.at(subject_proposed).begin(), find(t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_e_arg.at(subject_proposed)) ); //count the rank (distance from the first element) of the original t_e among t_nt_current_arg.at(subject_proposed) 
int rank_subject_x =0; // it is always zero as we are updating the sequence at its own infection

t_nt_modified_subject.erase(t_nt_modified_subject.begin() + rank_subject_x); // erase the original t_nt entry for subject_proposed

//---

nt_current_seq.assign(nt_modified_subject.begin()+n_base_CUPDATE*rank_subject_x , nt_modified_subject.begin()+n_base_CUPDATE*(rank_subject_x+1) ); //copy the original nt before erasing

nt_modified_subject.erase(nt_modified_subject.begin()+n_base_CUPDATE*rank_subject_x , nt_modified_subject.begin()+n_base_CUPDATE*(rank_subject_x+1) );  //erase the original nt entry for subject_proposed


//--

vector<int> nt_modified_source;
vector<double> t_nt_modified_source;

vector<int> infecting_list_modified_source;


int rank_source_x =-1;  //count the rank of the original t_e among t_nt_current_arg.at(subject_source) 

switch(subject_source ){

case 9999:{ // by background

t_up =  min( t_sample_arg.at(subject_proposed), min( t_i_arg.at(subject_proposed), t_max_CUPDATE));
t_low = max(0.0, t_up-t_back);
//t_low = max( t_e_arg.at(index_arg.at(0)), t_up-t_back);

// 	switch( t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
// 	
// 		case 1:{// with valid t_s
// 			double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 	
// 			switch(t_temp< t_sample_arg.at(subject_proposed)){
// 				case 1:{
// 					t_low =  max(0.0, t_temp);
// 				break;
// 				}
// 				case 0:{// should be unlikely if t_back is large enough
// 					double dt = t_temp -  t_sample_arg.at(subject_proposed);
// 					t_low =  max(0.0, t_sample_arg.at(subject_proposed) -  dt);
// 				break;
// 				}
// 			}
// 		break;
// 		}
// 	
// 		case 0:{ // no valid t_s
// 			t_low = max(0.0, t_up-t_back);
// 		break;
// 		}
// 	
// 	
// 	}

t_proposed= gsl_ran_flat(r_c,t_low, t_up );

//---

// 	t_up =  min( t_sample_arg.at(subject_proposed), min( t_i_arg.at(subject_proposed), t_max_CUPDATE));
// 	t_low = max(0.0, t_up-10.0);
// 
// 	t_proposed = t_e_arg.at(subject_proposed) + 0.1*gsl_ran_gaussian(r_c,1.0);
// 
// 	switch((t_proposed<t_low)| (t_proposed>t_up)){
// 		case 0:{
// 			//do nothing
// 		break;
// 		}
// 		case 1:{
// 			switch(t_proposed<t_low){
// 				case 0:{ // t_proposed>t_up
// 
// 					switch((t_proposed - t_up)>(t_up - t_low)){
// 
// 						case 0:{ 
// 						t_proposed = t_up -(t_proposed - t_up);
// 						break;
// 						}
// 
// 						case 1:{ 
// 						t_proposed = t_up - fmod(t_proposed - t_up, t_up - t_low);
// 						break;
// 						}
// 					}
// 
// 				break;
// 				}
// 
// 				case 1:{//t_proposed<t_low
// 
// 					switch(fabs(t_proposed - t_low)>(t_up - t_low)){
// 
// 						case 0:{ 
// 						t_proposed = t_low + fabs(t_proposed - t_low);
// 						break;
// 						}
// 
// 						case 1:{ 
// 						t_proposed = t_low + fmod(fabs(t_proposed - t_low), t_up - t_low);
// 						break;
// 						}
// 					}
// 
// 				break;
// 				}
// 			}
// 		break;
// 		}
// 	}

//--


break;
}

default :{ // not by background

nt_modified_source = nt_current_arg.at(subject_source);
t_nt_modified_source = t_nt_current_arg.at(subject_source);

t_up =   min( t_sample_arg.at(subject_proposed), min( min( t_i_arg.at(subject_proposed), t_r_arg.at(subject_source)),t_max_CUPDATE));
t_low = max(t_i_arg.at(subject_source), t_up -t_back );

// 	double t_temp_2 = min(t_sample_arg.at(subject_proposed), t_r_arg.at(subject_source));
// 	switch( t_temp_2!=unassigned_time_CUPDATE){
// 	
// 		case 1:{// with valid t_s (subject) or t_r(source)
// 			double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 			switch(t_temp< t_temp_2){
// 				case 1:{
// 					t_low =  max(t_i_arg.at(subject_source), t_temp);
// 				break;
// 				}
// 				case 0:{// should be unlikely if t_back is large enough
// 					double dt = t_temp -  t_temp_2;
// 					t_low =  max(t_i_arg.at(subject_source), t_temp_2 -  dt);
// 				break;
// 				}
// 			}
// 		break;
// 		}
// 	
// 		case 0:{ // no with valid t_s (subject) and t_r(source)
// 			t_low =  max(t_i_arg.at(subject_source), t_up - t_back);
// 		break;
// 		}
// 	}
	

t_proposed= gsl_ran_flat(r_c, t_low, t_up );

//-----

// 	t_up =   min( t_sample_arg.at(subject_proposed), min( min( t_i_arg.at(subject_proposed), t_r_arg.at(subject_source)),t_max_CUPDATE));
// 	t_low = max(t_i_arg.at(subject_source), t_up -10.0 );
// 
// 	t_proposed = t_e_arg.at(subject_proposed) + 0.1*gsl_ran_gaussian(r_c,1.0);
// 
// 	switch((t_proposed<t_low)| (t_proposed>t_up)){
// 		case 0:{
// 			//do nothing
// 		break;
// 		}
// 		case 1:{
// 			switch(t_proposed<t_low){
// 				case 0:{ // t_proposed>t_up
// 
// 					switch((t_proposed - t_up)>(t_up - t_low)){
// 
// 						case 0:{ 
// 						t_proposed = t_up -(t_proposed - t_up);
// 						break;
// 						}
// 
// 						case 1:{ 
// 						t_proposed = t_up - fmod(t_proposed - t_up, t_up - t_low);
// 						break;
// 						}
// 					}
// 
// 				break;
// 				}
// 
// 				case 1:{//t_proposed<t_low
// 
// 					switch(fabs(t_proposed - t_low)>(t_up - t_low)){
// 
// 						case 0:{ 
// 						t_proposed = t_low + fabs(t_proposed - t_low);
// 						break;
// 						}
// 
// 						case 1:{ 
// 						t_proposed = t_low + fmod(fabs(t_proposed - t_low), t_up - t_low);
// 						break;
// 						}
// 					}
// 
// 				break;
// 				}
// 			}
// 		break;
// 		}
// 	}

//--


rank_source_x = distance( t_nt_current_arg.at(subject_source).begin(), find(t_nt_current_arg.at(subject_source).begin(), t_nt_current_arg.at(subject_source).end(), t_e_arg.at(subject_proposed)) );


t_nt_modified_source.erase(t_nt_modified_source.begin() + rank_source_x); // erase the original t_nt entry for source
nt_modified_source.erase(nt_modified_source.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source

break;
}

}

t_e_modified.at(subject_proposed) = t_proposed;


//----------------------------------------------------------------------------------//

t_nt_modified_subject.push_back(t_proposed); // the unsorted t_nt with inserted t_proposed
sort( t_nt_modified_subject.begin(),  t_nt_modified_subject.end()); // the sorted t_nt with inserted t_proposed

//int rank_subject_y =distance( t_nt_modified_subject.begin(), find(t_nt_modified_subject.begin(), t_nt_modified_subject.end(), t_proposed) ); //count the NEW rank of the t_proposed  among t_nt_modifed.at(subject_proposed) 
int rank_subject_y = 0;

//nt_modified_subject.insert(nt_modified_subject.begin()+(rank_subject_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());  //insert  the new nt 

int rank_source_y =-1;  //count the NEW rank of the t_proposed  among t_nt_modifed.at(subject_source) 

switch(subject_source ){

case 9999:{ // by background
break;
}

default :{ // not by background

	
	//vector <double>  t_nt_source (t_nt_modified.at(subject_source)) ;
	
	t_nt_modified_source.push_back(t_proposed);
	
	sort( t_nt_modified_source.begin(),  t_nt_modified_source.end()); 

	rank_source_y = distance(t_nt_modified_source.begin(), find(t_nt_modified_source.begin(),t_nt_modified_source.end(), t_proposed) );

	//nt_modified_source.insert(nt_modified_source.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
	
	//------------

	infecting_list_modified_source = infecting_list_current_arg.at(subject_source);

	int rank_x = distance(infecting_list_modified_source.begin(), find(infecting_list_modified_source.begin(), infecting_list_modified_source.end(), subject_proposed));

	infecting_list_modified_source.erase(infecting_list_modified_source.begin()+rank_x);


	vector<double> t_y(infecting_size_current_arg.at(subject_source));
	for (int i=0;i<=(infecting_size_current_arg.at(subject_source)-1);i++){
	t_y.at(i) = t_e_modified.at(infecting_list_current_arg[subject_source][i]);
	}

	sort(t_y.begin(), t_y.end());

	int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_e_modified.at(subject_proposed)));
	infecting_list_modified_source.insert(infecting_list_modified_source.begin()+rank_y, subject_proposed);
	
	//--------------
	
	


break;
}

}



// 		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("t_e_proposed.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << subject_proposed <<   ","  << t_proposed <<","<< t_e_arg.at(subject_proposed) << "," <<  t_sample_arg.at(subject_proposed) - t_proposed  << endl;
// 		myfile_mcmc_out.close();

//---------------------------------------- proposing a new sequence & the proposal probability ----------------------------------------------------------//

vector<int> seq_proposed(n_base_CUPDATE);

double dt;

double t_past, t_future;

double log_pr_forward=0.0; // the log of proposal probability 

vector<int> nt_past_forward(n_base_CUPDATE); // the sequence at the nearest past (in the direction of time change) compared to the time of the proposed sequence; this might be or might not be the original sequence which gotta be replaced
vector<int> nt_future_forward(n_base_CUPDATE); // the sequence at the nearest future(in the direction of time change) compared to the time of the proposed sequence

dt = t_proposed - t_e_arg.at(subject_proposed); // the dimension of dt tells the direction of time change

//

double t_proposed_backward;

vector<int> seq_proposed_backward(n_base_CUPDATE);

double t_past_backward, t_future_backward;

double log_pr_backward=0.0; // the log of proposal probability 

vector<int> nt_past_backward(n_base_CUPDATE); 
vector<int> nt_future_backward(n_base_CUPDATE); 
//

switch(subject_source==9999){
	case 0:{ // NOT from background

		switch(current_size_arg.at(subject_proposed)>1){ // return 1 when the subject has more than one sequence available

			case 0:{ //  ONLY one sequence available for the subject
	
				switch(dt>=0){
					case 1:{ // propose the time to the right
	
						switch(rank_source_x==rank_source_y){
	
							case 1:{ // unchanged rank_source
								nt_past_forward = nt_current_seq;
								t_past = t_e_arg.at(subject_proposed);
								
								switch(rank_source_y==(current_size_arg.at(subject_source)-1)){ // see if the sequence is the last one
									case 1:{ // it is the last sequence
										t_future = t_proposed;
										seq_propose_uncond(seq_proposed, log_pr_forward, nt_past_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// no defined nt_future_forward
										
										//------------ backward direction------//
										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward,  nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // compute the backward proposal pr with future sequence
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // note: the argument takes t_future_backwards with value t_proposed_backward


										//------------------------------------------------//

									break;
									}
									case 0:{ // not the last sequence
										t_future = t_nt_modified_source.at(rank_source_y+1);
	
										nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+2)*n_base_CUPDATE); // the next sequence with rank=rank_source_y +1

										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// with defined nt_future_forward
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

										//------------------------------------------------//
										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);									
										//------------------------------------------------//
									
									break;
									}							
								}
		
							break;
							}
	
							case 0:{ //changed rank_source (rank_source_y will be > rank_source_x as the time is proposed to the right)
								nt_past_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE); // note: start point is not rank_source_y - 1 due to change of ranking of the sequence now before the proposed sequence
	
								t_past = t_nt_modified_source.at(rank_source_y-1);
	
	
								switch(rank_source_y==(current_size_arg.at(subject_source)-1)){ // see if the sequence is the last one
									case 1:{ // it is the last sequence
										t_future = t_proposed;
										seq_propose_uncond(seq_proposed,  log_pr_forward,nt_past_forward, t_proposed, t_past, t_future,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// no defined nt_future_forward

										//------------------------------------------------//

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//------------------------------------------------//
										
									break;
									}
									case 0:{ // not the last sequence
										t_future = t_nt_modified_source.at(rank_source_y+1);
	
										nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+2)*n_base_CUPDATE); // the next sequence with rank=rank_source_y +1; note: still rank_source_y+1 as the ranking of the sequence that still after proposed sequence deos not change

 										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// with defined nt_future_forward
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

										//------------------------------------------------//

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//------------------------------------------------//
									
									break;
									}							
								}
	
							break;
							}
						}
	
					break;
					}
	
					case 0:{  // propose the time to the left
	
						switch(rank_source_x==rank_source_y){
	
							case 1:{ // unchanged rank_source
	
								nt_past_forward = nt_current_seq;
								t_past = t_e_arg.at(subject_proposed);
	
								t_future = t_nt_modified_source.at(rank_source_y-1);
	
								nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE); // the previous sequence with rank=rank_source_y - 1 (it always exists as there is a least one sequence corresponds to the infection of the source itself

 								seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE,r_c);// with defined nt_future_forward
								//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

								//------------------------------------------------//

								switch(rank_source_x==(current_size_arg.at(subject_source)-1)){ // see if the (original) sequence was the last one

									case 0:{ // it was not the last sequence
										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;
		
										t_future_backward = t_nt_modified_source.at(rank_source_x+1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);
		
 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}

									case 1:{ // it was the last sequence

										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;
		
										t_future_backward = t_proposed_backward;

										seq_backward_pr_uncond(seq_proposed_backward,  log_pr_backward,nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // no future sequence

									break;
									}
								}
								//------------------------------------------------//

							break;
							}
	
							case 0:{ // changed rank_source (rank_source_y will be < rank_source_x as the time is proposed to the left)
	
								nt_past_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE); // note: start point is not rank_source_y - 1 due to change of ranking of the sequence now after the proposed sequence
	
								t_past = t_nt_modified_source.at(rank_source_y+1);
	
								t_future = t_nt_modified_source.at(rank_source_y-1);
	
								nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE); // the next sequence with rank=rank_source_y -1; note: still rank_source_y-1 as the ranking of the sequence that still before proposed sequence deos not change

 								seq_propose_cond(seq_proposed,  log_pr_forward,nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future,para_current_arg.mu_1, para_current_arg.mu_2,  n_base_CUPDATE, r_c);// with defined nt_future_forward
								//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

								//------------------------------------------------//

								switch(rank_source_x==(current_size_arg.at(subject_source)-1)){ // see if the (original) sequence was the last one

									case 0:{ // it was not the last sequence

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x+1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}

									case 1:{ // it was the last sequence

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_proposed_backward;

										seq_backward_pr_uncond(seq_proposed_backward,  log_pr_backward,nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // no future sequence

									break;
									}
								}
								//------------------------------------------------//
	
							break;
							}
	
						}
	
					break;
					}
				}

			break;
			}


			case 1:{ //  MORE than one sequence available for the subject
	
				switch(dt>=0){
					case 1:{ // propose the time to the right
	
						switch(rank_source_x==rank_source_y){
	
							case 1:{ // unchanged rank_source
								nt_past_forward = nt_current_seq;
								t_past = t_e_arg.at(subject_proposed);
								
								switch(rank_source_y==(current_size_arg.at(subject_source)-1)){ // see if the sequence is the last one (on source)
									case 1:{ // it is the last sequence

										//t_future = t_proposed;
										//seq_propose_uncond(seq_proposed, log_pr_forward, nt_past_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// no defined nt_future_forward

										t_future = t_nt_modified_subject.at(1);
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE); // the 2nd sequence; it cannot take over 2nd sequence as 2nd sequence only happens after time of becoming infectious where t_proposed cannot exceed

										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c); 
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

							
										//------------ backward direction------//
										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward,  nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // compute the backward proposal pr with future sequence
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); // note: the argument takes t_future_backwards with value t_proposed_backward


										//------------------------------------------------//

									break;
									}
									case 0:{ // not the last sequence

										switch(t_nt_modified_source.at(rank_source_y+1)<t_nt_modified_subject.at(1)){
											case 1:{
												t_future = t_nt_modified_source.at(rank_source_y+1);
												nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+2)*n_base_CUPDATE); // the next sequence with rank=rank_source_y +1
											break;
											}
											case 0:{
												t_future = t_nt_modified_subject.at(1);
												nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);		
											break;
											}
										}


										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// with defined nt_future_forward
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

										//------------------------------------------------//
										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);									
										//------------------------------------------------//
									
									break;
									}							
								}
		
							break;
							}
	
							case 0:{ //changed rank_source (rank_source_y will be > rank_source_x as the time is proposed to the right)
								nt_past_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE); // note: start point is not rank_source_y - 1 due to change of ranking of the sequence now before the proposed sequence
	
								t_past = t_nt_modified_source.at(rank_source_y-1);
	
	
								switch(rank_source_y==(current_size_arg.at(subject_source)-1)){ // see if the sequence is the last one
									case 1:{ // it is the last sequence

										//t_future = t_proposed;
										//seq_propose_uncond(seq_proposed,  log_pr_forward,nt_past_forward, t_proposed, t_past, t_future,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// no defined nt_future_forward

										t_future = t_nt_modified_subject.at(1);
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE); // the 2nd sequence; it cannot take over 2nd sequence as 2nd sequence only happens after time of becoming infectious where t_proposed cannot exceed

										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

											
										//------------------------------------------------//
										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//------------------------------------------------//
										
									break;
									}
									case 0:{ // not the last sequence

										switch(t_nt_modified_source.at(rank_source_y+1)<t_nt_modified_subject.at(1)){
											case 1:{
												t_future = t_nt_modified_source.at(rank_source_y+1);
												nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+2)*n_base_CUPDATE); // the next sequence with rank=rank_source_y +1
											break;
											}
											case 0:{
												t_future = t_nt_modified_subject.at(1);
												nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);		
											break;
											}
										}

 										seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// with defined nt_future_forward
										//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

										//------------------------------------------------//

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_source.at(rank_source_x-1);
										nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//------------------------------------------------//
									
									break;
									}							
								}
	
							break;
							}
						}
	
					break;
					}
	
					case 0:{  // propose the time to the left
	
						switch(rank_source_x==rank_source_y){
	
							case 1:{ // unchanged rank_source
	
								nt_past_forward = nt_current_seq;
								t_past = t_e_arg.at(subject_proposed);
	
								t_future = t_nt_modified_source.at(rank_source_y-1);
	
								nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE); // the previous sequence with rank=rank_source_y - 1 (it always exists as there is a least one sequence corresponds to the infection of the source itself

 								seq_propose_cond(seq_proposed,  log_pr_forward, nt_past_forward, nt_future_forward, t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE,r_c);// with defined nt_future_forward
								//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

								//------------------------------------------------//

								switch(rank_source_x==(current_size_arg.at(subject_source)-1)){ // see if the (original) sequence was the last one

									case 0:{ // it was not the last sequence

										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;
		

										switch(t_nt_modified_source.at(rank_source_x+1)<t_nt_modified_subject.at(1)){
											case 1:{
												t_future_backward = t_nt_modified_source.at(rank_source_x+1);
												nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);				
											break;
											}
											case 0:{
												t_future_backward = t_nt_modified_subject.at(1);
												nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);		
											break;
											}
										}


 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}

									case 1:{ // it was the last sequence

										t_past_backward = t_proposed;
										nt_past_backward = seq_proposed;
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;
		
										t_future_backward = t_nt_modified_subject.at(1);
										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);


									break;
									}
								}
								//------------------------------------------------//

							break;
							}
	
							case 0:{ // changed rank_source (rank_source_y will be < rank_source_x as the time is proposed to the left)
	
								nt_past_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y+1)*n_base_CUPDATE); // note: start point is not rank_source_y - 1 due to change of ranking of the sequence now after the proposed sequence
	
								t_past = t_nt_modified_source.at(rank_source_y+1);
	
								t_future = t_nt_modified_source.at(rank_source_y-1);
	
								nt_future_forward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_y)*n_base_CUPDATE); // the next sequence with rank=rank_source_y -1; note: still rank_source_y-1 as the ranking of the sequence that still before proposed sequence deos not change

 								seq_propose_cond(seq_proposed,  log_pr_forward,nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future,para_current_arg.mu_1, para_current_arg.mu_2,  n_base_CUPDATE, r_c);// with defined nt_future_forward
								//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_current_seq, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

								//------------------------------------------------//

								switch(rank_source_x==(current_size_arg.at(subject_source)-1)){ // see if the (original) sequence was the last one

									case 0:{ // it was not the last sequence

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										switch(t_nt_modified_source.at(rank_source_x+1)<t_nt_modified_subject.at(1)){
											case 1:{
												t_future_backward = t_nt_modified_source.at(rank_source_x+1);
												nt_future_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x+2)*n_base_CUPDATE);				
											break;
											}
											case 0:{
												t_future_backward = t_nt_modified_subject.at(1);
												nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);		
											break;
											}
										}


 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}

									case 1:{ // it was the last sequence

										t_past_backward = t_nt_modified_source.at(rank_source_x);
										nt_past_backward.assign(nt_current_arg.at(subject_source).begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_arg.at(subject_source).begin()+(rank_source_x)*n_base_CUPDATE);
										
										t_proposed_backward = t_e_arg.at(subject_proposed);
										seq_proposed_backward = nt_current_seq;

										t_future_backward = t_nt_modified_subject.at(1);
										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);

 										seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
										//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, seq_proposed,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);


									break;
									}
								}
								//------------------------------------------------//
	
							break;
							}
	
						}
	
					break;
					}
				}

			break;
			}

		}

	break;
	}

	case 1:{ //  from background

// 		sample_snull (con_seq_CUPDATE, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE, r_c); //sample a seq for background
// 		log_pr_forward = lh_snull(con_seq_CUPDATE, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);
// 
// 		log_pr_backward = lh_snull(con_seq_CUPDATE, nt_current_seq, para_current_arg.p_ber, n_base_CUPDATE);


		switch(current_size_arg.at(subject_proposed)>1){ // return 1 when the subject has more than one sequence available

			case 0:{ //  ONLY one sequence available for the subject

				//------------------------------------------------//

// 				t_past = t_e_arg.at(subject_proposed);
// 				nt_past_forward = nt_current_seq; 
// 				t_future = t_proposed;
// 				seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_future,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);// no defined nt_future_forward
// 				//------------------------------------------------//
// 
// 				t_past_backward = t_proposed;
// 				nt_past_backward =  seq_proposed;
// 				
// 				t_proposed_backward = t_e_arg.at(subject_proposed);
// 				seq_proposed_backward = nt_current_seq;
// 
// 				t_future_backward =t_proposed_backward;
// 
// 				seq_backward_pr_uncond(seq_proposed_backward,  log_pr_backward,nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); 

				//------------------------------------------------//

				sample_snull (con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE, r_c); //sample a seq for background
				log_pr_forward = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);
		
				log_pr_backward = lh_snull(con_seq, nt_current_seq, para_current_arg.p_ber, n_base_CUPDATE);

			break;
			}

			case 1:{ // MORE than one sequence available for the subject

				switch(dt>=0){
					case 1:{ // propose the time to the right
						t_past = t_e_arg.at(subject_proposed);
						nt_past_forward = nt_current_seq;
						
						t_future = t_nt_modified_subject.at(1);

	
						nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE); // the 2nd sequence; it cannot take over 2nd sequence as 2nd sequence only happens after time of becoming infectious where t_proposed cannot exceed
	
 						seq_propose_cond(seq_proposed,  log_pr_forward,nt_past_forward, nt_future_forward, t_proposed,  t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c); 
						//seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_proposed, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
						//------------------------------------------------//
		
						t_past_backward = t_proposed;
						nt_past_backward =  seq_proposed;
						
						t_proposed_backward = t_e_arg.at(subject_proposed);
						seq_proposed_backward = nt_current_seq;
		
						t_future_backward =t_proposed_backward;
		
						seq_backward_pr_uncond(seq_proposed_backward,  log_pr_backward,nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); 
		
		
						//------------------------------------------------//
	
					break;
					}
	
					case 0:{  // propose the time to the left 
						t_past = t_e_arg.at(subject_proposed);
						nt_past_forward = nt_current_seq;
						t_future = t_proposed;
						seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_future,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

						//------------------------------------------------//

						t_past_backward = t_proposed;
						nt_past_backward =  seq_proposed;
						
						t_proposed_backward = t_e_arg.at(subject_proposed);
						seq_proposed_backward = nt_current_seq;
	
						t_future_backward = t_nt_modified_subject.at(1);
						nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);
	
 						seq_backward_pr_cond(seq_proposed_backward, log_pr_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward, t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2,  n_base_CUPDATE);
						//seq_backward_pr_uncond(seq_proposed_backward, log_pr_backward, nt_past_backward,  t_proposed_backward, t_past_backward, t_proposed_backward,para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

						//------------------------------------------------//

					break;
					}
				}

			break;
			}

		}

	break;
	}
}

//-------------------------------------------end of proposing a new sequence ---------------------------------------------------------------------------------//



nt_modified_subject.insert(nt_modified_subject.begin()+(rank_subject_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());  //insert  the new nt 

switch(subject_source ){

case 9999:{ // by background
break;
}

default :{ // not by background

nt_modified_source.insert(nt_modified_source.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());

break;
}

}


//----------------------------------------------------------------------------------//

	switch (current_size_arg.at(subject_proposed)>1) {

	case 1:{

	log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below

	lh_square_modified.log_f_S.at(subject_proposed) = 0.0;

	for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){


	vector<int> seq_1(nt_modified_subject.begin()+j*(n_base_CUPDATE), nt_modified_subject.begin()+(j+1)*(n_base_CUPDATE));
	vector<int> seq_2(nt_modified_subject.begin()+(j+1)*(n_base_CUPDATE), nt_modified_subject.begin()+(j+2)*(n_base_CUPDATE));

	lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2, t_nt_modified_subject.at(j), t_nt_modified_subject.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

	}

	log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 

	break;
	}

	default:{
	break;
	}

	}

//----------------------------------------------------------------------------------//

switch (subject_source){
	
	case 9999:{ // by background

		log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 

		 lh_square_modified.log_f_Snull.at(subject_proposed)  = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);

		log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); 

	break;
	}
	
	default :{ // not by background
	
		
		switch (current_size_arg.at(subject_source)>1) {
	
		case 1:{
	
		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_source); //subtract part of likelihood that would be updated below
	
		lh_square_modified.log_f_S.at(subject_source) = 0.0;
	
		for (int j=0;j<=(current_size_arg.at(subject_source)-2);j++){
	
	
		vector<int> seq_1(nt_modified_source.begin()+j*(n_base_CUPDATE), nt_modified_source.begin()+(j+1)*(n_base_CUPDATE));
		vector<int> seq_2(nt_modified_source.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source.begin()+(j+2)*(n_base_CUPDATE));

		lh_square_modified.log_f_S.at(subject_source) =lh_square_modified.log_f_S.at(subject_source)+log_lh_seq(seq_1, seq_2, t_nt_modified_source.at(j), t_nt_modified_source.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
	
		}
	
		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_source); 
	
		break;
		}
	
		default:{
		break;
		}
	
		}
	break;
	}

}

//----------------------------------------------------------------------------------//

lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;

for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again

if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {

switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
case 1:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));

//lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j)) ;

break;
}
case 0:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
break;
}
}//end switch

lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

}
}

//----------
switch(infected_source_current_arg.at(subject_proposed)){

case 9999:{ // by background
//lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
break;
}

default :{ // not by background
lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][infected_source_current_arg.at(subject_proposed)]/norm_const_current_arg.at(infected_source_current_arg.at(subject_proposed)); // update k_sum_E
lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
break;
}

}


//---------------this part is needed when do not update index-----------------

// log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
// 
// lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
// // 		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
// lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
// lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
// 
// log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 

//--------------------------------
		
switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes

case 1:{

	switch(t_proposed<t_e_arg.at(index_arg.at(0))){ // original indexes would be replace by the chosen subject
	
	case 1:{

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("test_1_1.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
// myfile_mcmc_out.close();

	index_modified.clear();	
	index_modified.assign(1,subject_proposed);// replace index 
	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));

	for ( int i =0; i<= (int) (index_arg.size()-1); i++){
	xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
	}	

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
	
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;
	

	
	for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
	
	lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // this is not the subject_proposed
	lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
	lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
	lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
	}
	
	break;
	}
	
	
	case 0:{
	
	if (t_proposed==t_e_arg.at(index_arg.at(0))){ // addtion of one more index

	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_a.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();
	
		//if (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()){ // return 1 if proposed subject not one of the original indexes
	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
		index_modified.push_back(subject_proposed); // add index 
		xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.q_E.at(subject_proposed)=0.0;
		lh_square_modified.g_E.at(subject_proposed)=1.0;
		lh_square_modified.h_E.at(subject_proposed)=1.0;
		lh_square_modified.f_E.at(subject_proposed)=1.0;
	
		//}
	
	}
	
	if (t_proposed>t_e_arg.at(index_arg.at(0))){ // no shift of cases between xi_E and xi_E_minus



	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_b.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();

	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below


	
		lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
// 		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
		lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
	
		log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 


	} // end if t_proposs>t_e_arg.at()
	
	break;
	}
	
	}


break;
}


case 0: { // when chosen subject is one of the indexes



	index_modified.clear();

	int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
	double min_t = t_e_modified.at(first_min); // the minimum time of exposure
	
	int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time

// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("test_0.txt")).c_str(),ios::app);
// 	myfile_mcmc_out << first_min << "," << subject_proposed << "," <<  t_proposed << "," << num_min << endl;
// 	myfile_mcmc_out.close();
	
	switch (num_min>1) {
	case 1: {
	index_modified.reserve(n_CUPDATE);	
	for (int i=0; i<=(n_CUPDATE-1);i++){
	if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
	}
	break;
	}
	case 0:{
	index_modified.assign(1,first_min);
	break;
	}
	
	}

	xi_E_minus_modified = xi_E_arg;


	for (int i=0;i<= (int) (index_modified.size()-1); i++){

	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i))); // this subject would have to be removed from likelihood function ( new index might be orginally an index, but the the log(lh_square_modified.f_E.at(index_modified.at(i)) will be zero in this case)
	//lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.q_E.at(subject_proposed)=0.0;
	//lh_square_modified.g_E.at(subject_proposed)=1.0;
	//lh_square_modified.h_E.at(subject_proposed)=1.0;
	//lh_square_modified.f_E.at(subject_proposed)=1.0;

	lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.q_E.at(index_modified.at(i))=0.0;
	lh_square_modified.g_E.at(index_modified.at(i))=1.0;
	lh_square_modified.h_E.at(index_modified.at(i))=1.0;
	lh_square_modified.f_E.at(index_modified.at(i))=1.0;			
	}

	switch(find(index_modified.begin(),index_modified.end(),subject_proposed) ==index_modified.end() ){ //return 1 when the chosen  subject is NO longer an index
	case 1:{
/*
	ofstream myfile_mcmc_out; 
	myfile_mcmc_out.open((string(path4)+string("test_0_1.txt")).c_str(),ios::app);
	myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(index_modified.at(0)) << endl;
	myfile_mcmc_out.close();*/

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 	

	lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
	//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
	lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
	lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);

	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
	
	break;
	}
	case 0:{


	break;
	}
	}

break;
}

}
		
//--------------------//

switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_I.txt")).c_str(),ios::app);
// myfile_mcmc_out   <<  subject_proposed << "," << t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_I.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------

switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_EnI.txt")).c_str(),ios::app);
// myfile_mcmc_out  <<    subject_proposed << "," <<t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_EnI.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------

// switch(isfinite(exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward))){
// 	case 1:{
// 		acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward));
// 	break;
// 	}
// 
// 	case 0:{
// 		acp_pr =0.0;
// 	break;
// 	}
// }

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward));
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg)+(log_pr_backward-log_pr_forward)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
delta_mat_current_arg = delta_mat_modified;
log_lh_current_arg = log_lh_modified;
t_e_arg= t_e_modified;
index_arg = index_modified;
xi_E_minus_arg = xi_E_minus_modified;

// nt_current_arg = nt_modified;
// t_nt_current_arg = t_nt_modified;

nt_current_arg.at(subject_proposed) = nt_modified_subject;
t_nt_current_arg.at(subject_proposed) = t_nt_modified_subject;


	switch (subject_source){
		
		case 9999:{ // by background
		break;
		}
		
		default :{ // not by background
		nt_current_arg.at(subject_source) = nt_modified_source;
		t_nt_current_arg.at(subject_source) = t_nt_modified_source;
		infecting_list_current_arg.at(subject_source) = infecting_list_modified_source;		
		}
	}
			

break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);

//--------------
		



		// myfile_mcmc_out.open((string(path4)+string("find_xi_E_mibus.txt")).c_str(),ios::app);
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out << 1*(find(xi_E_minus_modified.begin(), xi_E_minus_modified.end(),index_modified.at(i))==xi_E_minus_modified.end()) << endl; // should always equals to 1, i.e., new index has been excluded from xi_E_minus
		// }
		// myfile_mcmc_out.close();
		// 

		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_kt_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.kt_sum_E.at(index_arg.at(0)) << endl; // shouls always be 0
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_k_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.k_sum_E.at(index_arg.at(0)) << endl;  // shouls always be 0
		// myfile_mcmc_out.close();

//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("index_t_e_seq.txt")).c_str(),ios::app);
// 		if (index_arg.empty()==0){
// 		for (int i=0; i<=(int)(index_arg.size()-1); i++){
// 		myfile_mcmc_out << index_arg.at(i) << "," << infected_source_current_arg.at(index_arg.at(i))  << endl;
// 		}
// 		}
// 		myfile_mcmc_out.close();

// 		myfile_mcmc_out.open((string(path4)+string("log_pr_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_pr_backward << "," << log_pr_forward << "," << exp(log_pr_backward-log_pr_forward) <<endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << "," << acp_pr << endl;
// 		myfile_mcmc_out.close();
// 
// 		
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_t_e_seq.txt")).c_str(),ios::app); 		
// 		myfile_mcmc_out << subject_proposed  << ","<< current_size_arg.at(subject_proposed) << "," << infected_source_current_arg.at(subject_proposed) <<endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("t_proposed_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(subject_proposed) << "," << t_e_arg.at(subject_proposed) << endl;
// 		myfile_mcmc_out.close();


// 		myfile_mcmc_out.open((string(path4)+string("p_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << p_1  << "," << p_2 << "," << p_3 << endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("m_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << m_1  << "," << m_2 << "," << m_3 << "," << (long double) pow( p_1,m_1)*(long double) pow( p_2, m_2)*(long double) pow( p_3,m_3)<< endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("n_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << n_1  << "," << n_2 << "," << n_3 << "," <<(long double)  pow( p_1,n_1)*(long double) pow( p_2, n_2)*(long double) pow( p_3,n_3)<< endl;
// 		myfile_mcmc_out.close();





		//---------------

}

/*------------------------------------------------*/
void mcmc_UPDATE::source_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vec2int& infecting_list_current_arg, vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, int& subject_proposed, vector<int>& list_update, int iter){

double acp_pr = 0.0;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

double log_pr_seq_forward=0.0; 
double log_pr_seq_backward=0.0;

double log_pr_ds_forward=0.0;
double log_pr_ds_backward=0.0;


double t_e_subject = t_e_arg.at(subject_proposed);

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector<int> current_size_modified = current_size_arg;

vec2int nt_modified = nt_current_arg;

vec2int infecting_list_modified= infecting_list_current_arg;
vector<int> infecting_size_modified = infecting_size_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter); // set a seed

//vector<int> infected_source_modified = infected_source_current_arg;

vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
//vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

int source_x = infected_source_current_arg.at(subject_proposed);

int rank_source_x; // the rank of removed sequence in old source (note: it has different meaning in function t_e_seq in which the source remains the same)
vector<int> nt_current_source_x;
vector<double> t_nt_current_source_x;
vector<int> nt_modified_source_x;
vector<double> t_nt_modified_source_x;

// vector<int> nt_subject_seq; // the orginal first sequence of the subject
// nt_subject_seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);


//-------------- propose a new source --------------//

int source_y;
int rank_source_y;
vector<int> nt_current_source_y;
vector<double> t_nt_current_source_y;
vector<int> nt_modified_source_y;
vector<double> t_nt_modified_source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

	switch( (t_i_arg.at(xi_I_arg.at(i))<t_e_subject) & (t_r_arg.at(xi_I_arg.at(i))>=t_e_subject)){

		case 1:{
			source_pool.push_back(xi_I_arg.at(i));	
		break;
		}
		case 0:{
		break;
		}
	}
}
	
// 		if (source_x!=9999){
// 		ofstream myfile_mcmc_out;
// 		myfile_mcmc_out.open((string(path4)+string("00_test.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << 1*(find(source_pool.begin(),source_pool.end(), source_x)==source_pool.end())  << endl;
// 		myfile_mcmc_out.close();
// 		}
// 
// 		if ((source_x!=9999) & ((find(source_pool.begin(),source_pool.end(), source_x)==source_pool.end())==1)){
// 		ofstream myfile_mcmc_out;
// 		myfile_mcmc_out.open((string(path4)+string("00_test2.txt")).c_str(),ios::out);
// 		for (int i=0; i<=(int) source_pool.size()-1;i++){
// 		myfile_mcmc_out << source_pool.at(i)  << endl;
// 		}
// 		myfile_mcmc_out << "source_x"<<endl;
// 		myfile_mcmc_out << source_x<<endl;
// 		myfile_mcmc_out << endl;	
// 		myfile_mcmc_out << t_i_arg.at(source_x) <<"," << t_r_arg.at(source_x) <<"," << t_e_subject << endl;
// 	
// 		myfile_mcmc_out.close();
// 		}



source_pool.insert(source_pool.begin(),9999);

int num_infectious = (int)source_pool.size();

//-----------------------------propose uniformly-------------------------------------------//

// //source_y = source_x;
// // while((source_y==source_x) & (num_infectious>1)){
// //source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_infectious)); // uniformly choose a new source (including bg)
// // }
// //infected_source_modified.at(subject_proposed) = source_y; 

//source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_infectious)); // uniformly choose a new source (including bg)


//-----propose according to infectious challenge--------------------//

vector<double> ic(num_infectious);
ic.at(0) = para_current_arg.alpha;

switch(num_infectious>=2){

	case 1:{ // with 2nd sources from pool

		for (int j=1;j<=(num_infectious-1);j++){
		ic.at(j)= para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_pool.at(j)]/norm_const_current_arg.at(source_pool.at(j)); // a new source will be proposed according to the infectious challenges
		}

		double *P=&ic.at(0); // convert vector to array
		gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
		int link= gsl_ran_discrete (r_c, g);
		gsl_ran_discrete_free (g);

		source_y = source_pool.at(link); // a new source
		log_pr_forward = log(ic.at(link));

		switch(source_x==9999){
			case 0:{
				double ic_source_x =  para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_x]/norm_const_current_arg.at(source_x);
				log_pr_backward =log(ic_source_x);
			break;
			}
			case 1:{
				double ic_source_x =  para_current_arg.alpha;
				log_pr_backward =log(ic_source_x);
			break;
			}
		}


		break;
	}

	case 0:{ // only primary source from pool

		source_y = 9999;
		log_pr_forward = log(para_current_arg.alpha);

		double ic_source_x =  para_current_arg.alpha;
		log_pr_backward =log(ic_source_x);

	break;
	}
}


//--------------//

//----------end of proposing a new source -----//


//----------------------------------------------------------------------------------------------------------------//

switch(source_y==source_x){

	case 0:{
		
		vector<int> seq_proposed(n_base_CUPDATE); // newly proposed sequence when source changes
		vector<int> nt_past_forward(n_base_CUPDATE);
		vector<int> nt_future_forward(n_base_CUPDATE);
		double t_proposed, t_past, t_future;
		
		vector<int> seq_proposed_backward(n_base_CUPDATE);
		vector<int> nt_past_backward(n_base_CUPDATE);
		vector<int> nt_future_backward(n_base_CUPDATE);
		double t_proposed_backward, t_past_backward, t_future_backward;
		
		//------------------------------------------------
		
		switch(source_y==9999){
		
			case 0:{// new 2nd infection

				infecting_size_modified.at(source_y) = infecting_size_modified.at(source_y) +1;
				
				vector<double> t_y(infecting_size_current_arg.at(source_y));
				for (int i=0;i<=(infecting_size_current_arg.at(source_y)-1);i++){
				t_y.at(i) = t_e_arg.at(infecting_list_current_arg[source_y][i]);
				}
				t_y.push_back(t_e_subject);
				sort(t_y.begin(), t_y.end());
			
				int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_e_subject));
				infecting_list_modified.at(source_y).insert(infecting_list_modified.at(source_y).begin()+rank_y, subject_proposed);
				
				//----------------------------------------------------//
	
				nt_current_source_y = nt_current_arg.at(source_y);
				t_nt_current_source_y = t_nt_current_arg.at(source_y);
		
				nt_modified_source_y = nt_current_source_y;
				t_nt_modified_source_y = t_nt_current_source_y;
		
				t_nt_modified_source_y.push_back(t_e_subject);
		
				sort( t_nt_modified_source_y.begin(),  t_nt_modified_source_y.end()); 
		
				rank_source_y = distance(t_nt_modified_source_y.begin(), find(t_nt_modified_source_y.begin(),t_nt_modified_source_y.end(), t_e_subject) );
		
				current_size_modified.at(source_y) = current_size_modified.at(source_y) + 1;
		
				//----------------------------------------------------------------------------------------------------------------//
				t_past = t_nt_modified_source_y.at(rank_source_y-1);
				nt_past_forward.assign(nt_current_source_y.begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE);
		
				t_proposed = t_e_subject;
		
				//switch(current_size_arg.at(subject_proposed)>1){
				switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with sample in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y (take closer seq as the future seq)
		
								//switch(t_nt_current_arg[subject_proposed][1]<t_nt_modified_source_y.at(rank_source_y +1)){
								switch(t_sample_arg.at(subject_proposed)<t_nt_modified_source_y.at(rank_source_y +1)){
		
									case 1:{// take the one from subject as future seq
										
										int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );
		
										t_future  = t_sample_arg.at(subject_proposed); 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);
										
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
		
									case 0:{ // take the one from source_y as future seq
										t_future = t_nt_modified_source_y.at(rank_source_y+1);
										nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
								}
		
							break;
							}
		
							case 0:{//  inserted seq will be last seq at rank_source_y  in source_y
		
										
										int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );
		
										t_future  = t_sample_arg.at(subject_proposed); 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);
										
										
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
		
		
					break;
					}
		
					case 0:{// with no sample in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y
								t_future = t_nt_modified_source_y.at(rank_source_y+1);
								nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
								seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);						
							break;
							}
		
							case 0:{// inserted seq will be last seq at rank_source_y  in source_y
								t_future = t_proposed;
								seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				//nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_subject_seq.begin(), nt_subject_seq.end());
				nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
		
				nt_modified.at(source_y) = nt_modified_source_y;
		
			break;
			}
		
			case 1:{// new bg infection

				//switch(current_size_arg.at(subject_proposed)>1){ 
				switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with sample in subject
						int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );

						t_past  = t_sample_arg.at(subject_proposed); //yes, t_past!
						nt_past_forward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);

						t_proposed = t_e_subject;
						t_future = t_proposed;

						seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);			
		
					break;
					}
		
					case 0:{// with no sample in subject
		
						log_pr_seq_forward = n_base_CUPDATE*log(0.25);
		
						for (int i=0; i<=(n_base_CUPDATE-1);i++){
						seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
						}
		
					break;
					}
				}
			break;
			}
		
		}
		
		//--------------------------------------------
		nt_modified_subject.erase(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE); 
		nt_modified_subject.insert(nt_modified_subject.begin(), seq_proposed.begin(), seq_proposed.end());
		nt_modified.at(subject_proposed) = nt_modified_subject;
		//---------------------------------------------
		
		switch(source_x==9999){
		
			case 0:{// was 2nd infection

				infecting_size_modified.at(source_x) = infecting_size_modified.at(source_x) - 1;

				int rank_x = distance(infecting_list_current_arg.at(source_x).begin(), find(infecting_list_current_arg.at(source_x).begin(), infecting_list_current_arg.at(source_x).end(), subject_proposed));
			
				infecting_list_modified.at(source_x).erase(infecting_list_modified.at(source_x).begin()+rank_x);

				//----------------------------------------------------------------------------------------------------------------//

				nt_current_source_x = nt_current_arg.at(source_x);
				t_nt_current_source_x = t_nt_current_arg.at(source_x);
		
				rank_source_x = distance( t_nt_current_source_x.begin(), find(t_nt_current_source_x.begin(), t_nt_current_source_x.end(), t_e_subject) ); 
		
				//----------------------------------------------------------------------------------------------------------------//
		
				t_past_backward = t_nt_current_source_x.at(rank_source_x-1);
				nt_past_backward.assign(nt_current_source_x.begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x)*n_base_CUPDATE);
		
				t_proposed_backward = t_e_subject;
				seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
				//switch(current_size_arg.at(subject_proposed)>1){ 
				switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
				
					case 1:{// with sample subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
		
							case 1:{// also NOT last seq at rank_source_x  in source_x (take closer seq as the future seq)
								//switch(t_nt_current_arg[subject_proposed][1]<t_nt_current_source_x.at(rank_source_x +1)){
								switch(t_sample_arg.at(subject_proposed)<t_nt_current_source_x.at(rank_source_x +1)){

									case 1:{// take the one from subject as future seq

										int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );
		
										t_future_backward  = t_sample_arg.at(subject_proposed); 
										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);
										
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
									break;
									}
									case 0:{// take the one from source_x as future seq
										t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
										nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
						
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
									break;
									}
								}
							break;
							}
		
							case 0:{ //last seq at rank_source_x  in source_x

								int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );

								t_future_backward  = t_sample_arg.at(subject_proposed); 
								nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);
												
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}	
						}
		
		
					break;
					}
					
					case 0:{ // with no sample in the subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
							case 1:{ // not last seq at rank_source_x  in source_x
								t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
								nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
								
							break;
							}
							case 0:{ // last seq rank_source_x  in source_x
								t_future_backward = t_proposed_backward;
		
								seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				nt_modified_source_x = nt_current_source_x;
				t_nt_modified_source_x = t_nt_current_source_x;
				
				t_nt_modified_source_x.erase(t_nt_modified_source_x.begin() + rank_source_x); // erase the original t_nt entry for source_x
				nt_modified_source_x.erase(nt_modified_source_x.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source_x.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source_x
				nt_modified.at(source_x) = nt_modified_source_x;

				current_size_modified.at(source_x) = current_size_modified.at(source_x) - 1;
		
			break;
			}
		
			case 1:{// was bg infection
		
				//switch(current_size_arg.at(subject_proposed)>1){
				switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){

					case 1:{// with sample in subject

						int rank_sample_subject = distance( t_nt_current_arg.at(subject_proposed).begin(), find( t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_sample_arg.at(subject_proposed)) );

						t_past_backward  = t_sample_arg.at(subject_proposed);  // note: it is "past"!
						nt_past_backward.assign(nt_current_arg.at(subject_proposed).begin()+ rank_sample_subject*n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ (rank_sample_subject+1)*n_base_CUPDATE);
							

						t_proposed_backward = t_e_subject;
						seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
						t_future_backward = t_proposed_backward;
		
						seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
		
					break;
					}
		
					case 0:{// with no sample in subject
						log_pr_seq_backward = n_base_CUPDATE*log(0.25); // note: this is the proposal density! not refers to the model!
					break;
					}
				}
		
			break;
			}
		
		}
		
		



	//------------- deal with change of likelihood due to change of source and sequences in subject_proposed)---------------//
		
		switch (current_size_arg.at(subject_proposed)>1) {
		
		case 1:{
		
		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below
		
		lh_square_modified.log_f_S.at(subject_proposed) = 0.0;
		
		for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){
		
		
		vector<int> seq_1(nt_modified.at(subject_proposed).begin()+j*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE));
		vector<int> seq_2(nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+2)*(n_base_CUPDATE));
		
		lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2, t_nt_current_arg[subject_proposed][j], t_nt_current_arg[subject_proposed][j+1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
		}
		
		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 
		
		break;
		}
		
		default:{
		break;
		}
		
		}
		
		//---
		switch(source_x){
		
			case 9999:{ // was background infection
		
				switch(source_y){
					case 9999:{ // new  bg infection 		
					break;
					}
		
					default:{ // new secondary infection
		
						log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 
		
						lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) = 0.0;
		
						log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); // in fact redudancy, but here for clarity
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
			break;
			}
			
			default :{ // was secondary infection
		
				switch(source_y){
					case 9999:{ // new  bg infection
		
						log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); // redudant as second term must be zero
		
		
						lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
						lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
						lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) =  n_base_CUPDATE*log(0.25);
		
						log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
			
		
				
					break;
					}
		
					default:{ // new secondary infection
			
						log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
		
						lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
		
						//--
		
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
		
						//---
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
		
			break;
			}
		
		}
		
		
		
	//------------- end of with change of likelihood (due to change of source and sequences in subject_proposed)---------------//
	
		//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward)*exp(log_pr_seq_backward-log_pr_seq_forward)*exp(log_pr_ds_backward-log_pr_ds_forward));

		acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+(log_pr_seq_backward-log_pr_seq_forward)+ (log_pr_ds_backward-log_pr_ds_forward)) );
		
		
		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		
		lh_square_current_arg = lh_square_modified;
		log_lh_current_arg = log_lh_modified;
		
		//nt_current_arg.at(subject_proposed) = nt_modified_subject;
		nt_current_arg = nt_modified;

		current_size_arg = current_size_modified;
		infected_source_current_arg.at(subject_proposed) =  source_y;
		//infected_source_current_arg = infected_source_modified;

		infecting_list_current_arg = infecting_list_modified;
		infecting_size_current_arg = infecting_size_modified;
		
			switch (source_x){
				
				case 9999:{ 
				break;
				}
				
				default :{ 
				//nt_current_arg.at(source_x) = nt_modified_source_x;
				t_nt_current_arg.at(source_x) = t_nt_modified_source_x;	
				break;	
				}
			}
		
			switch (source_y){
				
				case 9999:{ 
				break;
				}
				
				default :{ 
				//nt_current_arg.at(source_y) = nt_modified_source_y;
				t_nt_current_arg.at(source_y) = t_nt_modified_source_y;
				break;	
				}
			}

	break;
	}

	case 0: {
	break;
	}
}

gsl_rng_free(r_c);

break;
}

case 1:{ // source_y==source_x
break;
}
}

//--------------------------------------------------

ofstream myfile_mcmc_out; 

myfile_mcmc_out.open((string(path4)+string("subect_proposed_source_update.txt")).c_str(),ios::app);
if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<subject_proposed << ","<<  source_x <<","<< source_y << "," << infecting_size_current_arg.at(subject_proposed) << ","<< (int) list_update.size() << endl;
myfile_mcmc_out.close();	

myfile_mcmc_out.open((string(path4)+string("log_lh_change_source_update.txt")).c_str(),ios::app);
if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_lh_current_arg <<","<< log_lh_modified << "," << source_x <<","<< source_y <<","<< acp_pr << endl;
myfile_mcmc_out.close();	

myfile_mcmc_out.open((string(path4)+string("log_pr_forward_backward_source_update.txt")).c_str(),ios::app);
//if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) )
//myfile_mcmc_out << exp(log_lh_modified-log_lh_current_arg) <<"," << exp(log_pr_backward-log_pr_forward) << ","<<  exp(log_pr_seq_backward-log_pr_seq_forward) << "," << acp_pr << endl;
if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_pr_seq_backward <<","<< log_pr_seq_forward<< "," << log_pr_ds_backward <<","<< log_pr_ds_forward << endl;
myfile_mcmc_out.close();	


}

/*------------------------------------------------*/

void mcmc_UPDATE::source_t_e_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vec2int& infecting_list_current_arg, vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, int& subject_proposed, vector<int>& list_update,vector<int>& con_seq, int iter, gsl_rng* & r_c){

//double t_back =10.0;

double acp_pr = 0.0;
double t_low, t_up;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

double log_pr_t_e_forward=0.0; 
double log_pr_t_e_backward=0.0;

double log_pr_seq_forward=0.0; 
double log_pr_seq_backward=0.0;

double log_pr_ds_forward=0.0;
double log_pr_ds_backward=0.0;


double t_proposed;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;

vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

// vector <int> xi_U_modified = xi_U_arg;
// vector <int> xi_E_modified = xi_E_arg;
// vector <int> xi_EnI_modified = xi_EnI_arg;

vector<int> current_size_modified = current_size_arg;

vec2int nt_modified = nt_current_arg;
vec2d t_nt_modified = t_nt_current_arg;

vec2int infecting_list_modified= infecting_list_current_arg;
vector<int> infecting_size_modified = infecting_size_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,iter); // set a seed

//vector<int> infected_source_modified = infected_source_current_arg;

vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

//vector<int> nt_subject_seq; // the orginal first sequence of the subject
//nt_subject_seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);

int source_x = infected_source_current_arg.at(subject_proposed);

int rank_source_x; // the rank of removed sequence in old source (note: it has different meaning in function t_e_seq in which the source remains the same)
vector<int> nt_current_source_x;
vector<double> t_nt_current_source_x;
vector<int> nt_modified_source_x;
vector<double> t_nt_modified_source_x;



//-------------- propose a new source --------------//

int source_y;
int rank_source_y;
vector<int> nt_current_source_y;
vector<double> t_nt_current_source_y;
vector<int> nt_modified_source_y;
vector<double> t_nt_modified_source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject

double t_bound = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

	switch( t_i_arg.at(xi_I_arg.at(i))<t_bound){
//	switch( (t_i_arg.at(xi_I_arg.at(i))<t_bound) & ((t_bound  - t_r_arg.at(xi_I_arg.at(i)))<=t_back) ){

		case 1:{
			source_pool.push_back(xi_I_arg.at(i));	
		break;
		}
		case 0:{
		break;
		}
	}
}
	

source_pool.insert(source_pool.begin(),9999);

int num_infectious = (int)source_pool.size();

//-----------------------------propose uniformly-------------------------------------------//

//source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_infectious)); // uniformly choose a new source (including bg)

//-----propose according to infectious challenge--------------------//

vector<double> ic(num_infectious);
ic.at(0) = para_current_arg.alpha;

switch(num_infectious>=2){

	case 1:{ // with 2nd sources from pool

		for (int j=1;j<=(num_infectious-1);j++){
		ic.at(j)= para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_pool.at(j)]/norm_const_current_arg.at(source_pool.at(j)); // a new source will be proposed according to the infectious challenges
		}

		double *P=&ic.at(0); // convert vector to array
		gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
		int link= gsl_ran_discrete (r_c, g);
		gsl_ran_discrete_free (g);

		source_y = source_pool.at(link); // a new source
		log_pr_forward = log(ic.at(link));

		switch(source_x==9999){
			case 0:{
				double ic_source_x =  para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_x]/norm_const_current_arg.at(source_x);
				log_pr_backward =log(ic_source_x);
			break;
			}
			case 1:{
				double ic_source_x =  para_current_arg.alpha;
				log_pr_backward =log(ic_source_x);
			break;
			}
		}


		break;
	}

	case 0:{ // only primary source from pool

		source_y = 9999;
		log_pr_forward = log(para_current_arg.alpha);

		double ic_source_x =  para_current_arg.alpha;
		log_pr_backward =log(ic_source_x);

	break;
	}
}


//--------------//

//----------end of proposing a new source -----//

//----------------------------------------------------------------------------------------------------------------//

switch(source_y==source_x){

	case 0:{


		//-------- propose a new t_e---------------//
		switch(source_y){
		
		case 9999:{ // by background
		
			t_up = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(0.0, t_up-t_back);
		
// 			switch( t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 			
// 					switch(t_temp< t_sample_arg.at(subject_proposed)){
// 						case 1:{
// 							t_low =  max(0.0, t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_sample_arg.at(subject_proposed);
// 							t_low =  max(0.0, t_sample_arg.at(subject_proposed) -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no valid t_s
// 					t_low = max(0.0, t_up-t_back);
// 				break;
// 				}
// 			
// 			
// 			}
			
			t_proposed= gsl_ran_flat(r_c,t_low, t_up );
			
			log_pr_t_e_forward = log(1.0/(t_up-t_low));
		
		break;
		}
		
		default :{ // not by background
		
			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_y)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));			
			t_low = max(t_i_arg.at(source_y), t_up -t_back );

// 			double t_temp_2 = min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_y));
// 			switch( t_temp_2!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s (subject) or t_r(source)
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 					switch(t_temp< t_temp_2){
// 						case 1:{
// 							t_low =  max(t_i_arg.at(source_y), t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_temp_2;
// 							t_low =  max(t_i_arg.at(source_y), t_temp_2 -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no with valid t_s (subject) and t_r(source)
// 					t_low =  max(t_i_arg.at(source_y), t_up - t_back);
// 				break;
// 				}
// 			}

				
			t_proposed= gsl_ran_flat(r_c, t_low, t_up );

			log_pr_t_e_forward = log(1.0/(t_up-t_low));

		break;
		}
		
		}

		//------------//
		t_e_modified.at(subject_proposed) = t_proposed;
		//------------//

		switch(source_x){
		
		case 9999:{ // by background
		
			t_up = min( t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(0.0, t_up-t_back);

// 			switch( t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 			
// 					switch(t_temp< t_sample_arg.at(subject_proposed)){
// 						case 1:{
// 							t_low =  max(0.0, t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_sample_arg.at(subject_proposed);
// 							t_low =  max(0.0, t_sample_arg.at(subject_proposed) -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no valid t_s
// 					t_low = max(0.0, t_up-t_back);
// 				break;
// 				}
// 			
// 			
// 			}

						
			log_pr_t_e_backward = log(1.0/(t_up-t_low));
		
		break;
		}
		
		default :{ // not by background
		
			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_x)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(t_i_arg.at(source_x), t_up -t_back );

// 			double t_temp_2 = min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_x));	
// 			switch( t_temp_2!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s (subject) or t_r(source)
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 					switch(t_temp< t_temp_2){
// 						case 1:{
// 							t_low =  max(t_i_arg.at(source_x), t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_temp_2;
// 							t_low =  max(t_i_arg.at(source_x), t_temp_2 -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no with valid t_s (subject) and t_r(source)
// 					t_low =  max(t_i_arg.at(source_x), t_up - t_back);
// 				break;
// 				}
// 			}


			log_pr_t_e_backward = log(1.0/(t_up-t_low));

		break;
		}
		
		}

		//--------------------------------------//

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00test.txt")).c_str(),ios::app);
// myfile_mcmc_out << 0 << endl;
// myfile_mcmc_out.close();		
		
		vector<int> seq_proposed(n_base_CUPDATE); // newly proposed sequence when source changes
		vector<int> nt_past_forward(n_base_CUPDATE);
		vector<int> nt_future_forward(n_base_CUPDATE);
		double t_past, t_future;
		
		vector<int> seq_proposed_backward(n_base_CUPDATE);
		vector<int> nt_past_backward(n_base_CUPDATE);
		vector<int> nt_future_backward(n_base_CUPDATE);
		double t_proposed_backward, t_past_backward, t_future_backward;
		
		//------------------------------------------------
		
		switch(source_y==9999){
		
			case 0:{// new 2nd infection

				infecting_size_modified.at(source_y) = infecting_size_modified.at(source_y) +1;
				
				vector<double> t_y(infecting_size_current_arg.at(source_y));
				for (int i=0;i<=(infecting_size_current_arg.at(source_y)-1);i++){
				t_y.at(i) = t_e_arg.at(infecting_list_current_arg[source_y][i]);
				}
				t_y.push_back(t_proposed);
				sort(t_y.begin(), t_y.end());
			
				int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_proposed));
				infecting_list_modified.at(source_y).insert(infecting_list_modified.at(source_y).begin()+rank_y, subject_proposed);
				
				//----------------------------------------------------//
	
				nt_current_source_y = nt_current_arg.at(source_y);
				t_nt_current_source_y = t_nt_current_arg.at(source_y);
		
				nt_modified_source_y = nt_current_source_y;
				t_nt_modified_source_y = t_nt_current_source_y;
		
				t_nt_modified_source_y.push_back(t_proposed);
		
				sort( t_nt_modified_source_y.begin(),  t_nt_modified_source_y.end()); 
		
				rank_source_y = distance(t_nt_modified_source_y.begin(), find(t_nt_modified_source_y.begin(),t_nt_modified_source_y.end(), t_proposed) );
		
				current_size_modified.at(source_y) = current_size_modified.at(source_y) + 1;
		
				//----------------------------------------------------------------------------------------------------------------//
				t_past = t_nt_modified_source_y.at(rank_source_y-1);
				nt_past_forward.assign(nt_current_source_y.begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE);
		
				//t_proposed = t_e_subject;
		
				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y (take closer seq as the future seq)
		
								switch(t_nt_current_arg[subject_proposed][1]<t_nt_modified_source_y.at(rank_source_y +1)){
								//switch(t_sample_arg.at(subject_proposed)<t_nt_modified_source_y.at(rank_source_y +1)){
		
									case 1:{// take the one from subject as future seq

										t_future  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
																				
									seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
		
									case 0:{ // take the one from source_y as future seq
										t_future = t_nt_modified_source_y.at(rank_source_y+1);
										nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
								}
		
							break;
							}
		
							case 0:{//  inserted seq will be last seq at rank_source_y  in source_y
		
										t_future  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
	
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y
								t_future = t_nt_modified_source_y.at(rank_source_y+1);
								nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
								seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);						
							break;
							}
		
							case 0:{// inserted seq will be last seq at rank_source_y  in source_y
								t_future = t_proposed;
								seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				//nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_subject_seq.begin(), nt_subject_seq.end());
				nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
		
				nt_modified.at(source_y) = nt_modified_source_y;
				t_nt_modified.at(source_y) = t_nt_modified_source_y;

			break;
			}
		
			case 1:{// new bg infection

// 				sample_snull (con_seq_CUPDATE, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE, r_c); //sample a seq for background
// 				log_pr_seq_forward = lh_snull(con_seq_CUPDATE, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);


				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject

						t_past  = t_nt_current_arg[subject_proposed][1]; // yes, t_past!
						nt_past_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
		
						t_future = t_proposed;

						seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);			
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
	
// 						log_pr_seq_forward = n_base_CUPDATE*log(0.25);
// 		
// 						for (int i=0; i<=(n_base_CUPDATE-1);i++){
// 						seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
// 						}
		

						 sample_snull (con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE, r_c); //sample a seq for background
						 log_pr_seq_forward = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);

					break;
					}
				}


			break;
			}
		
		}
		
		//--------------------------------------------
		nt_modified_subject.erase(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE); 
		nt_modified_subject.insert(nt_modified_subject.begin(), seq_proposed.begin(), seq_proposed.end());

		t_nt_modified_subject.erase(t_nt_modified_subject.begin());
		t_nt_modified_subject.push_back(t_proposed); 
		sort( t_nt_modified_subject.begin(),  t_nt_modified_subject.end()); 

		nt_modified.at(subject_proposed) = nt_modified_subject;
		t_nt_modified.at(subject_proposed) = t_nt_modified_subject;
		//---------------------------------------------
		
		switch(source_x==9999){
		
			case 0:{// was 2nd infection

				infecting_size_modified.at(source_x) = infecting_size_modified.at(source_x) - 1;

				int rank_x = distance(infecting_list_current_arg.at(source_x).begin(), find(infecting_list_current_arg.at(source_x).begin(), infecting_list_current_arg.at(source_x).end(), subject_proposed));
			
				infecting_list_modified.at(source_x).erase(infecting_list_modified.at(source_x).begin()+rank_x);

				//----------------------------------------------------------------------------------------------------------------//

				nt_current_source_x = nt_current_arg.at(source_x);
				t_nt_current_source_x = t_nt_current_arg.at(source_x);
		
				rank_source_x = distance( t_nt_current_source_x.begin(), find(t_nt_current_source_x.begin(), t_nt_current_source_x.end(), t_e_arg.at(subject_proposed)) ); 
		
				//----------------------------------------------------------------------------------------------------------------//
		
				t_past_backward = t_nt_current_source_x.at(rank_source_x-1);
				nt_past_backward.assign(nt_current_source_x.begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x)*n_base_CUPDATE);
		
				t_proposed_backward =t_e_arg.at(subject_proposed);
				seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
				
					case 1:{// with2nd seq subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
		
							case 1:{// also NOT last seq at rank_source_x  in source_x (take closer seq as the future seq)
								switch(t_nt_current_arg[subject_proposed][1]<t_nt_current_source_x.at(rank_source_x +1)){
								//switch(t_sample_arg.at(subject_proposed)<t_nt_current_source_x.at(rank_source_x +1)){

									case 1:{// take the one from subject as future seq

										t_future_backward  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
						
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}
									case 0:{// take the one from source_x as future seq
										t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
										nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
						
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
									break;
									}
								}
							break;
							}
		
							case 0:{ //last seq at rank_source_x  in source_x

								t_future_backward  = t_nt_current_arg[subject_proposed][1]; // always takes the 2nd seq in subject as the future sequence (both for forward and backward)
								nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}	
						}
		
		
					break;
					}
					
					case 0:{ // with no 2nd seq in the subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
							case 1:{ // not last seq at rank_source_x  in source_x
								t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
								nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
								
							break;
							}
							case 0:{ // last seq rank_source_x  in source_x
								t_future_backward = t_proposed_backward;
		
								seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				nt_modified_source_x = nt_current_source_x;
				t_nt_modified_source_x = t_nt_current_source_x;
				
				t_nt_modified_source_x.erase(t_nt_modified_source_x.begin() + rank_source_x); // erase the original t_nt entry for source_x
				nt_modified_source_x.erase(nt_modified_source_x.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source_x.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source_x

				nt_modified.at(source_x) = nt_modified_source_x;
				t_nt_modified.at(source_x) = t_nt_modified_source_x;


				current_size_modified.at(source_x) = current_size_modified.at(source_x) - 1;
		
			break;
			}
		
			case 1:{// was bg infection

// 				vector<int> seq(n_base_CUPDATE);
// 				seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE );
// 				log_pr_seq_backward = lh_snull(con_seq_CUPDATE, seq, para_current_arg.p_ber, n_base_CUPDATE);


				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){

					case 1:{// with 2nd seq in subject

						t_past_backward = t_nt_current_arg[subject_proposed][1]; // note: it is "past"!
						nt_past_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
		
						t_proposed_backward = t_e_arg.at(subject_proposed);
						seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
						t_future_backward = t_proposed_backward;

						seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
// 						log_pr_seq_backward = n_base_CUPDATE*log(0.25); // note: this is the proposal density! not refers to the model!
						vector<int> seq(n_base_CUPDATE);
						seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE );
						log_pr_seq_backward = lh_snull(con_seq, seq, para_current_arg.p_ber, n_base_CUPDATE);
					break;
					}
				}


			break;
			}
		
		}
		
		


		//------------- deal with change of likelihood due to change of source and sequences in subject_proposed)---------------//
		
		switch (current_size_arg.at(subject_proposed)>1) {
		
		case 1:{
		
		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below
		
		lh_square_modified.log_f_S.at(subject_proposed) = 0.0;
		
		for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){

			vector<int> seq_1(nt_modified.at(subject_proposed).begin()+j*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE));
			vector<int> seq_2(nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+2)*(n_base_CUPDATE));
			
			lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2,t_nt_modified_subject.at(j), t_nt_modified_subject.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);	
		}
		
		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 
		
		break;
		}
		
		default:{
		break;
		}
		
		}


		//-----------------------------------------------------------------------------------------//

		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
		
		for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again
		
			if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
			
				switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
					case 1:{
						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
					
					break;
					}
					case 0:{
						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
					break;
					}
				}
			lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));
	
			}
		}
		//----------
		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;

		switch(source_y){
		
		case 9999:{ // by background
		lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
		break;
		}
		
		default :{ // not by background
		lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y);
		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		break;
		}
		
		}


		switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes
		
			case 1:{
			
				switch(t_proposed<t_e_arg.at(index_arg.at(0))){ // original indexes would be replace by the chosen subject
				
					case 1:{
				
						index_modified.clear();	
						index_modified.assign(1,subject_proposed);// replace index 
						xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
					
						for ( int i =0; i<= (int) (index_arg.size()-1); i++){
						xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
						}	
					
						log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
						
						lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
						lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
						lh_square_modified.q_E.at(subject_proposed)=0.0;
						lh_square_modified.g_E.at(subject_proposed)=1.0;
						lh_square_modified.h_E.at(subject_proposed)=1.0;
						lh_square_modified.f_E.at(subject_proposed)=1.0;
						
					
						for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
						
						//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
						
						lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // this is not the subject_proposed
						lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
						lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
						lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
						
						log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
						}
						
					break;
					}
					
				
					case 0:{
					
						if (t_proposed==t_e_arg.at(index_arg.at(0))){ // addtion of one more index
						
							log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
							index_modified.push_back(subject_proposed); // add index 
							xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
							lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
							lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
							lh_square_modified.q_E.at(subject_proposed)=0.0;
							lh_square_modified.g_E.at(subject_proposed)=1.0;
							lh_square_modified.h_E.at(subject_proposed)=1.0;
							lh_square_modified.f_E.at(subject_proposed)=1.0;
						
						
						}
						
						if (t_proposed>t_e_arg.at(index_arg.at(0))){ // no shift of cases between xi_E and xi_E_minus

							log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below

							lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
					// 		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
							lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
							lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
						
							log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
					
					
						} // end if t_proposs>t_e_arg.at()
						
					break;
					}
				
				}
			
			break;
			}
		
		
			case 0: { // when chosen subject is one of the indexes

				index_modified.clear();
			
				int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
				double min_t = t_e_modified.at(first_min); // the minimum time of exposure
				
				int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time
			
				switch (num_min>1) {
					case 1: {
						index_modified.reserve(n_CUPDATE);	
						for (int i=0; i<=(n_CUPDATE-1);i++){
						if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
						}
					break;
					}
					case 0:{
						index_modified.assign(1,first_min);
					break;
					}
				}
			
				xi_E_minus_modified = xi_E_arg;
			
			
				for (int i=0;i<= (int) (index_modified.size()-1); i++){
			
				xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));
			
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i))); // this subject would have to be removed from likelihood function ( new index might be orginally an index, but the the log(lh_square_modified.f_E.at(index_modified.at(i)) will be zero in this case)
			
				lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
				lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
				lh_square_modified.q_E.at(index_modified.at(i))=0.0;
				lh_square_modified.g_E.at(index_modified.at(i))=1.0;
				lh_square_modified.h_E.at(index_modified.at(i))=1.0;
				lh_square_modified.f_E.at(index_modified.at(i))=1.0;			
				}
			
				switch(find(index_modified.begin(),index_modified.end(),subject_proposed) ==index_modified.end() ){ //return 1 when the chosen  subject is NO longer an index
					case 1:{
				
						log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 	
					
						lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
						lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
					
						log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
						
					break;
					}
					case 0:{
					break;
					}
				}
			
			break;
			}
		
		}
		
		//--------------------//


		switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
			case 1:{
		
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
				lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
				log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 
				
			break;
			}

			case 0:{
			break;
			}
		}
		
		//----------
		
		switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
			case 1:{
			
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
				lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
				log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 
				
			break;
			}
			case 0:{
			break;
			}
		}


		//-----------------------------------------------------------------------------------------//

		switch(source_x){
		
			case 9999:{ // was background infection
		
				switch(source_y){
					case 9999:{ // new  bg infection 		
					break;
					}
		
					default:{ // new secondary infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) = 0.0;
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); // in fact redudancy, but here for clarity
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
			break;
			}
			
			default :{ // was secondary infection
		
				switch(source_y){
					case 9999:{ // new  bg infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); // redudant as second term must be zero
		
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
// 						lh_square_modified.log_f_Snull.at(subject_proposed) =  n_base_CUPDATE*log(0.25);
						lh_square_modified.log_f_Snull.at(subject_proposed) = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);

		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
			
		
				
					break;
					}
		
					default:{ // new secondary infection
			
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
		
						//--
		
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
		
						//---
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
		
			break;
			}
		
		}
		
		
		
	//------------- end of with change of likelihood (due to change of source and sequences in subject_proposed)---------------//
	
		//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward)*exp(log_pr_seq_backward-log_pr_seq_forward)*exp(log_pr_ds_backward-log_pr_ds_forward));



		acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+(log_pr_seq_backward-log_pr_seq_forward) +(log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_ds_backward-log_pr_ds_forward)) );
		
		
		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		
		lh_square_current_arg = lh_square_modified;
		log_lh_current_arg = log_lh_modified;
		
		//nt_current_arg.at(subject_proposed) = nt_modified_subject;
		nt_current_arg = nt_modified;
		t_nt_current_arg = t_nt_modified;

		delta_mat_current_arg = delta_mat_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;

		current_size_arg = current_size_modified;
		infected_source_current_arg.at(subject_proposed) =  source_y;
		//infected_source_current_arg = infected_source_modified;

		infecting_list_current_arg = infecting_list_modified;
		infecting_size_current_arg = infecting_size_modified;
		
// 			switch (source_x){
// 				
// 				case 9999:{ 
// 				break;
// 				}
// 				
// 				default :{ 
// 				//nt_current_arg.at(source_x) = nt_modified_source_x;
// 				t_nt_current_arg.at(source_x) = t_nt_modified_source_x;	
// 				break;	
// 				}
// 			}
// 		
// 			switch (source_y){
// 				
// 				case 9999:{ 
// 				break;
// 				}
// 				
// 				default :{ 
// 				//nt_current_arg.at(source_y) = nt_modified_source_y;
// 				t_nt_current_arg.at(source_y) = t_nt_modified_source_y;
// 				break;	
// 				}
// 			}

	break;
	}

	case 0: {
	break;
	}
}

break;
}

case 1:{ // source_y==source_x
break;
}
}

//gsl_rng_free(r_c);

//--------------------------------------------------

// ofstream myfile_mcmc_out; 
// 
// myfile_mcmc_out.open((string(path4)+string("subect_proposed_source_te_update.txt")).c_str(),ios::app);
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<subject_proposed << ","<<  source_x <<","<< source_y << "," << infecting_size_current_arg.at(subject_proposed) << ","<< (int) list_update.size() << endl;
// myfile_mcmc_out.close();	
// 
// myfile_mcmc_out.open((string(path4)+string("log_lh_change_source_te_update.txt")).c_str(),ios::app);
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_lh_current_arg <<","<< log_lh_modified << "," << source_x <<","<< source_y <<","<< acp_pr << endl;
// myfile_mcmc_out.close();	
// 
// myfile_mcmc_out.open((string(path4)+string("log_pr_forward_backward_source_te_update.txt")).c_str(),ios::app);
// //if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) )
// //myfile_mcmc_out << exp(log_lh_modified-log_lh_current_arg) <<"," << exp(log_pr_backward-log_pr_forward) << ","<<  exp(log_pr_seq_backward-log_pr_seq_forward) << "," << acp_pr << endl;
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_pr_seq_backward <<","<< log_pr_seq_forward<< "," << log_pr_ds_backward <<","<< log_pr_ds_forward << endl;
// myfile_mcmc_out.close();	


}

/*------------------------------------------------*/

void mcmc_UPDATE::t_e_replct(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, int iter){

int subject_proposed ;
double t_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed




	
subject_proposed = xi_E_arg.at(gsl_rng_uniform_int (r_c, xi_E_arg.size())); // gsl_rng_uniform_int : a int random number in [0,xi_E_arg.size()-1] will be drawn

int subject_source = infected_source_current_arg.at(subject_proposed);


//ofstream myfile_mcmc_out; 
//myfile_mcmc_out.open((string(path4)+string("subject_proposed_replct.txt")).c_str(),ios::app); 		//myfile_mcmc_out << subject_proposed  << endl;
//myfile_mcmc_out.close();		

//--

switch(infected_source_current_arg.at(subject_proposed)){

case 9999:{ // by background
t_proposed= gsl_ran_flat(r_c, 0.0, min( t_sample_arg.at(subject_proposed), min( t_i_arg.at(subject_proposed), t_max_CUPDATE)) );
break;
}

default :{ // not by background
t_proposed= gsl_ran_flat(r_c,  t_i_arg.at(subject_source), min(  t_sample_arg.at(subject_proposed), min( min( t_i_arg.at(subject_proposed), t_r_arg.at(subject_source)),t_max_CUPDATE)) );
break;
}

}

//--

//t_proposed= t_e_arg.at(subject_proposed) + 1.0*gsl_ran_gaussian(r_c,1.0); 
// double ub = min(t_max_CUPDATE,t_i_arg.at(subject_proposed) ) ; // upper bound, note: t_i is an extreme value when it has not gone through class I
// switch((t_proposed>ub)|(t_proposed<0.0)) { 
// case 1: {
// if (t_proposed>ub) t_proposed = ub - ( (t_proposed - ub) - (ub-0.0)*floor((t_proposed-ub)/(ub-0.0)) ); //reflection at ub with period ub-0
// if (t_proposed<0.0)  t_proposed =  0.0 + ( (0.0 - t_proposed) - (ub-0.0)*floor((0.0-t_proposed)/(ub-0.0)) ); //reflection at 0 with period ub-0
// break;
// }
// case 0: {
// break;
// }
// }


//t_proposed = t_e_arg.at(index_arg.at(0)) +10.0;


t_e_modified.at(subject_proposed) = t_proposed;
// ofstream myfile_mcmc; 
// myfile_mcmc.open((string(path4)+string("t_proposed_current_index.txt")).c_str(),ios::app);
// myfile_mcmc <<  t_proposed << ","<<index_arg.at(0) << ","<<t_e_arg.at(index_arg.at(0)) << endl;
// myfile_mcmc.close();

//-----------

lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;

for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again

if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {

switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
case 1:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));

//lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j)) ;

break;
}
case 0:{
delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
break;
}
}//end switch

lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

}
}

//----------
switch(infected_source_current_arg.at(subject_proposed)){

case 9999:{ // by background
//lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
break;
}

default :{ // not by background
lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][infected_source_current_arg.at(subject_proposed)]/norm_const_current_arg.at(infected_source_current_arg.at(subject_proposed)); // update k_sum_E
lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
break;
}

}
//----------
		
switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes

case 1:{

	switch(t_proposed<t_e_arg.at(index_arg.at(0))){ // original indexes would be replace by the chosen subject
	
	case 1:{

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("test_1_1.txt")).c_str(),ios::app);
// myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
// myfile_mcmc_out.close();

	index_modified.clear();	
	index_modified.assign(1,subject_proposed);// replace index 
	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));

	for ( int i =0; i<= (int) (index_arg.size()-1); i++){
	xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
	}	

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
	
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;
	

	
	for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
	
	lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // this is not the subject_proposed
	lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
	lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
	lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
	}
	
	break;
	}
	
	
	case 0:{
	
	if (t_proposed==t_e_arg.at(index_arg.at(0))){ // addtion of one more index

	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_a.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();
	
		//if (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()){ // return 1 if proposed subject not one of the original indexes
	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
		index_modified.push_back(subject_proposed); // add index 
		xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
		lh_square_modified.q_E.at(subject_proposed)=0.0;
		lh_square_modified.g_E.at(subject_proposed)=1.0;
		lh_square_modified.h_E.at(subject_proposed)=1.0;
		lh_square_modified.f_E.at(subject_proposed)=1.0;
	
		//}
	
	}
	
	if (t_proposed>t_e_arg.at(index_arg.at(0))){ // no shift of cases between xi_E and xi_E_minus



	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_1_0_b.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_arg.at(index_arg.at(0)) << endl;
	//myfile_mcmc_out.close();

	
		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below


	
		lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
// 		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
		lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
	
		log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 


	} // end if t_proposs>t_e_arg.at()
	
	break;
	}
	
	}


break;
}


case 0: { // when chosen subject is one of the indexes



	index_modified.clear();

	int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
	double min_t = t_e_modified.at(first_min); // the minimum time of exposure
	
	int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time

// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("test_0.txt")).c_str(),ios::app);
// 	myfile_mcmc_out << first_min << "," << subject_proposed << "," <<  t_proposed << "," << num_min << endl;
// 	myfile_mcmc_out.close();
	
	switch (num_min>1) {
	case 1: {
	index_modified.reserve(n_CUPDATE);	
	for (int i=0; i<=(n_CUPDATE-1);i++){
	if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
	}
	break;
	}
	case 0:{
	index_modified.assign(1,first_min);
	break;
	}
	
	}

	xi_E_minus_modified = xi_E_arg;


	for (int i=0;i<= (int) (index_modified.size()-1); i++){

	xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i))); // this subject would have to be removed from likelihood function ( new index might be orginally an index, but the the log(lh_square_modified.f_E.at(index_modified.at(i)) will be zero in this case)
	//lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	//lh_square_modified.q_E.at(subject_proposed)=0.0;
	//lh_square_modified.g_E.at(subject_proposed)=1.0;
	//lh_square_modified.h_E.at(subject_proposed)=1.0;
	//lh_square_modified.f_E.at(subject_proposed)=1.0;

	lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
	lh_square_modified.q_E.at(index_modified.at(i))=0.0;
	lh_square_modified.g_E.at(index_modified.at(i))=1.0;
	lh_square_modified.h_E.at(index_modified.at(i))=1.0;
	lh_square_modified.f_E.at(index_modified.at(i))=1.0;			
	}

	switch(find(index_modified.begin(),index_modified.end(),subject_proposed) ==index_modified.end() ){ //return 1 when the chosen  subject is NO longer an index
	case 1:{
/*
	ofstream myfile_mcmc_out; 
	myfile_mcmc_out.open((string(path4)+string("test_0_1.txt")).c_str(),ios::app);
	myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(index_modified.at(0)) << endl;
	myfile_mcmc_out.close();*/

	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 	

	lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
	//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
	lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
	lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);

	log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
	
	break;
	}
	case 0:{

	//ofstream myfile_mcmc_out; 
	//myfile_mcmc_out.open((string(path4)+string("test_0_0.txt")).c_str(),ios::app);
	//myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(index_modified.at(0)) << endl;
	//myfile_mcmc_out.close();

	break;
	}
	}

break;
}

}
		
		
//--------------------//

switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_I.txt")).c_str(),ios::app);
// myfile_mcmc_out   <<  subject_proposed << "," << t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_I.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------

switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
case 1:{

log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("f_EnI.txt")).c_str(),ios::app);
// myfile_mcmc_out  <<    subject_proposed << "," <<t_i_arg.at(subject_proposed) << "," <<  t_proposed << "," <<  lh_square_modified.f_EnI.at(subject_proposed) << endl;
// myfile_mcmc_out.close();

break;
}
case 0:{
break;
}
}

//----------
		
acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));

//acp_pr =1.0;

double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
delta_mat_current_arg = delta_mat_modified;
log_lh_current_arg = log_lh_modified;
t_e_arg= t_e_modified;
index_arg = index_modified;
xi_E_minus_arg = xi_E_minus_modified;
break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);

//--------------
		
// 		ofstream myfile_mcmc_out; 
// 		
// 		myfile_mcmc_out.open((string(path4)+string("t_e_proposed_acp.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  subject_proposed  << "," <<  acp_pr  <<   "," << t_proposed <<   ","  << t_e_arg.at(subject_proposed) << endl;
// 		myfile_mcmc_out.close();
		
		// myfile_mcmc_out.open((string(path4)+string("index_modified.txt")).c_str(),ios::app);
		// if (index_modified.empty()==0){
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out <<  index_modified.at(i) << "," <<  t_e_modified.at(index_modified.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg.txt")).c_str(),ios::app);
		// if (index_arg.empty()==0){
		// for (int i=0; i<=(int)(index_arg.size()-1); i++){
		// myfile_mcmc_out << index_arg.at(i) << "," << t_e_arg.at(index_arg.at(i)) <<endl;
		// }
		// }
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("find_xi_E_mibus.txt")).c_str(),ios::app);
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out << 1*(find(xi_E_minus_modified.begin(), xi_E_minus_modified.end(),index_modified.at(i))==xi_E_minus_modified.end()) << endl; // should always equals to 1, i.e., new index has been excluded from xi_E_minus
		// }
		// myfile_mcmc_out.close();
		// 

		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_kt_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.kt_sum_E.at(index_arg.at(0)) << endl; // shouls always be 0
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_k_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.k_sum_E.at(index_arg.at(0)) << endl;  // shouls always be 0
		// myfile_mcmc_out.close();

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_replct.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();

		//---------------

gsl_rng_free(r_c);

}


/*------------------------------------------------*/

void mcmc_UPDATE::t_e_seq_add_del(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg,  vector<int>& xi_EnIS_arg,const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vec2int& infecting_list_current_arg, vector<int>& infecting_size_current_arg, vector<int>& con_seq,int iter){

//double t_back =10.0;

int subject_proposed ;


double t_proposed = 0.0;
double acp_pr = 0.0;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

double log_pr_t_e_forward=0.0; 
double log_pr_t_e_backward=0.0;

double log_pr_seq_forward=0.0; 
double log_pr_seq_backward=0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;

vector <int> index_modified = index_arg;

vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;
vector <int> xi_EnIS_modified = xi_EnIS_arg;

vector<int> infected_source_modified = infected_source_current_arg;

vector<int> current_size_modified = current_size_arg;

vec2int nt_modified = nt_current_arg;
vec2d t_nt_modified = t_nt_current_arg;

vec2int infecting_list_modified= infecting_list_current_arg;
vector<int> infecting_size_modified = infecting_size_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed

//---------------------

double P [2] ={1.0,1.0};
gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (2,P);
int action = gsl_ran_discrete (r_c, g); // a int random number in [0,1] will be drawn according to the weights in P :  0=deletion of an infection into gp of xi_EnI (choose from the gp of xi_U); 1= addition of an infection from gp of xi_EnI (add back  to xi_U)
gsl_ran_discrete_free (g);

if (xi_EnIS_arg.empty()==1){ // return 1 if xi_EnIS is empty
action = 1; // only addtion possible
}

if (xi_U_arg.empty()==1){ // return 1 if xi_U is empty
action = 0; // only deletion is possible
}



switch (action) {

case 0:{ // deletion

	subject_proposed = xi_EnIS_arg.at(gsl_rng_uniform_int (r_c, xi_EnIS_arg.size()));

	int source_x = infected_source_current_arg.at(subject_proposed);

	xi_EnIS_modified.erase(find(xi_EnIS_modified.begin(),xi_EnIS_modified.end(),subject_proposed));
	xi_EnI_modified.erase(find(xi_EnI_modified.begin(),xi_EnI_modified.end(),subject_proposed));	
	xi_E_modified.erase(find(xi_E_modified.begin(),xi_E_modified.end(),subject_proposed));
	xi_U_modified.push_back(subject_proposed);
	
	infected_source_modified.at(subject_proposed) = -99; // this subject become NO infected source	
	t_e_modified.at(subject_proposed) = unassigned_time_CUPDATE;

	nt_modified.at(subject_proposed).clear(); 
	t_nt_modified.at(subject_proposed).clear();
	current_size_modified.at(subject_proposed) = 0;

	vector<int> nt_current_source_x ;
	vector<double> t_nt_current_source_x;

	vector<int> nt_modified_source_x ;
	vector<double> t_nt_modified_source_x ;

	vector<int> source_pool;
	
	source_pool = xi_I_arg;

// 	for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){
// 	
// 		switch(t_i_arg.at(xi_I_arg.at(i))>=(t_max_CUPDATE-10.0)){// 
// 	
// 			case 1:{
// 				source_pool.push_back(xi_I_arg.at(i));	
// 			break;
// 			}
// 			case 0:{
// 			break;
// 			}
// 		}
// 	}

	source_pool.insert(source_pool.begin(),9999);

	int num_infectious = (int) source_pool.size();


	switch(source_x==9999){

		case 0:{// delete a 2nd infection

			current_size_modified.at(source_x) = current_size_arg.at(source_x) - 1;

			infecting_size_modified.at(source_x) = infecting_size_current_arg.at(source_x) - 1;

			int rank_x = distance(infecting_list_current_arg.at(source_x).begin(), find(infecting_list_current_arg.at(source_x).begin(), infecting_list_current_arg.at(source_x).end(), subject_proposed));
		
			infecting_list_modified.at(source_x).erase(infecting_list_modified.at(source_x).begin()+rank_x);

			//---
			nt_current_source_x = nt_current_arg.at(source_x);
			t_nt_current_source_x = t_nt_current_arg.at(source_x);

			int rank_source_x = distance( t_nt_current_source_x.begin(), find(t_nt_current_source_x.begin(), t_nt_current_source_x.end(), t_e_arg.at(subject_proposed)) ); 

			nt_modified_source_x = nt_current_source_x;
			t_nt_modified_source_x = t_nt_current_source_x;
				
			t_nt_modified_source_x.erase(t_nt_modified_source_x.begin() + rank_source_x); // erase the original t_nt entry for source_x
			nt_modified_source_x.erase(nt_modified_source_x.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source_x.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source_x

			nt_modified.at(source_x) = nt_modified_source_x;
			t_nt_modified.at(source_x) = t_nt_modified_source_x;

			//----

			log_pr_backward = log(1.0/((double)(xi_U_arg.size())+1.0)) + log( 1.0/((double)(num_infectious)) );// the 2nd term corresponds to the random choice of a source (including background)
			log_pr_forward = log(1.0/((double)(xi_EnIS_arg.size())));

			double t_up = min(t_r_arg.at(source_x), t_max_CUPDATE);			
//			double t_low = t_i_arg.at(source_x);
			double t_low = max(t_i_arg.at(source_x), t_up- t_back);

		
			log_pr_t_e_backward = log(1.0/(t_up-t_low));			
			log_pr_t_e_forward = 0.0;
			
			//---

			double t_past_backward = t_nt_current_source_x.at(rank_source_x-1);
			vector<int> nt_past_backward;
			nt_past_backward.assign(nt_current_source_x.begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x)*n_base_CUPDATE);
		
			double t_proposed_backward =t_e_arg.at(subject_proposed);
			vector<int> seq_proposed_backward;
			seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);

			switch(current_size_arg.at(source_x)>(rank_source_x+1)){

				case 1:{// NOT last seq at rank_source_x  in source_x 

					double t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
					vector<int> nt_future_backward;
					nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
	
					seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

				break;
				}

				case 0:{ //last seq at rank_source_x  in source_x

					double t_future_backward = t_proposed_backward;

					seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

				break;
				}	
			}
			
			//----
	
		break;
		}

		case 1:{ // delete a background infection

			log_pr_backward = log(1.0/((double)(xi_U_arg.size())+1.0)) + log(1.0/((double)(num_infectious)));// the 2nd term corresponds to the random choice of a source (including background)
			log_pr_forward = log(1.0/((double)(xi_EnIS_arg.size())));

			double t_up = t_max_CUPDATE;
			double t_low = max(0.0,t_up - t_back);				
			//double t_low = max(t_e_arg.at(index_arg.at(0)),t_up - t_back);	

			log_pr_t_e_backward = log(1.0/(t_up-t_low));			
			log_pr_t_e_forward = 0.0;

//			log_pr_seq_backward = n_base_CUPDATE*log(0.25);

			vector<int> seq(n_base_CUPDATE);
			seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE );
			log_pr_seq_backward = lh_snull(con_seq, seq, para_current_arg.p_ber, n_base_CUPDATE);

		break;
		}
	}

	//------------- deal with change of likelihood in f_s and f_null due to deletion of an infection---------------//

	switch(source_x==9999){

		case 0:{// deleted a 2nd infection
			log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x);
		
			lh_square_modified.log_f_S.at(source_x) = 0.0;

			switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
				case 1:{		
					for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
	
					vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
					vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
			
					lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
				
					}
				break;
				}

				case 0:{
				break;
				}
			}			

			log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x			
		break;
		}
		
		case 1:{//deleted a background infection
			log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed);
			lh_square_modified.log_f_Snull.at(subject_proposed) = 0.0;
			log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);
		break;
		}
	}
	//--------------------------------------------------------------------------------------------------------------------------//

	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E 
	
		//if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
		switch (t_r_arg.at(xi_I_arg.at(j))>=t_max_CUPDATE) {
			case 1:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
			//lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)] ;
			break;
			}
			case 0:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
			break;
			}
		}//end switch
		
		lh_square_modified.kt_sum_U.at(subject_proposed) = lh_square_modified.kt_sum_U.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
		//}
	}
	
	
	//-----
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  
	
	lh_square_modified.q_T.at(subject_proposed) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_U.at(subject_proposed);
	lh_square_modified.f_U.at(subject_proposed) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(subject_proposed),1.0);
	
	log_lh_modified = log_lh_modified  + log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes =/=1
	
	//----
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  //f_EnI becomes 1
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0;
	
	//---
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));  //f_E becomes 1
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;

	
	//----- this part needed when do not update index
// 		xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));

	//---------------

		switch(find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) {// returns 1 if it was NOT an index originally
		case 0:{// was an index
	
			index_modified.clear();
		
			int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
			double min_t = t_e_modified.at(first_min); // the minimum time of exposure
			
			int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time
			
				switch (num_min>1) {
				case 1: {
				index_modified.reserve(n_CUPDATE);	
				for (int i=0; i<=(n_CUPDATE-1);i++){
				if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
				}
				break;
				}
				case 0:{
				index_modified.assign(1,first_min);
				break;
				}
				
				}
	
			//xi_E_minus_modified = xi_E_arg;
		
			for (int i=0;i<= (int) (index_modified.size()-1); i++){
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));
			log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i)));
			lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.q_E.at(index_modified.at(i))=0.0;
			lh_square_modified.g_E.at(index_modified.at(i))=1.0;
			lh_square_modified.h_E.at(index_modified.at(i))=1.0;
			lh_square_modified.f_E.at(index_modified.at(i))=1.0;	
			}
					
		break;
		}	
	
		case 1:{// was not an index
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		break;
		}
		}

		//-----------------------------------------//
		acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+ (log_pr_t_e_backward-log_pr_t_e_forward) +(log_pr_seq_backward-log_pr_seq_forward) ) );
	
		
		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
	
		switch(uniform_rv<=acp_pr){
			case 1: {
				lh_square_current_arg = lh_square_modified;
				delta_mat_current_arg = delta_mat_modified;
				log_lh_current_arg = log_lh_modified;
				t_e_arg= t_e_modified;
				index_arg = index_modified;
				xi_E_minus_arg = xi_E_minus_modified;
				xi_U_arg = xi_U_modified;
				xi_E_arg = xi_E_modified;
			
				xi_EnI_arg = xi_EnI_modified;
				xi_EnIS_arg = xi_EnIS_modified;
			
				infected_source_current_arg = infected_source_modified;
				current_size_arg = current_size_modified;
				infecting_list_current_arg = infecting_list_modified;
				infecting_size_current_arg = infecting_size_modified;
			
				nt_current_arg = nt_modified;
				t_nt_current_arg = t_nt_modified;
			
			break;
			}
			
			case 0: {
			break;
			}
		}

// 		ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("11test.txt")).c_str(),ios::app);
// 		if (uniform_rv<=acp_pr) myfile_mcmc_out <<action<< "," << acp_pr << "," << source_x << ","<< subject_proposed << endl;
// 		myfile_mcmc_out.close();	

break;
}// end of deletion

//-----------------------------------------------------------------------------------------------------------------------------------------------//

case 1:{ // addition
	
	subject_proposed = xi_U_arg.at(gsl_rng_uniform_int (r_c, xi_U_arg.size()));
	
	xi_U_modified.erase(find(xi_U_modified.begin(),xi_U_modified.end(),subject_proposed)); // alter the gps xi_U, xi_E, xi_EnI
	xi_EnIS_modified.push_back(subject_proposed);
	xi_EnI_modified.push_back(subject_proposed);
	xi_E_modified.push_back(subject_proposed);

	vector<int> source_pool;

	source_pool = xi_I_arg;

// 	for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){
// 	
// 		switch(t_i_arg.at(xi_I_arg.at(i))>=(t_max_CUPDATE-10.0)){
// 			case 1:{
// 				source_pool.push_back(xi_I_arg.at(i));	
// 			break;
// 			}
// 			case 0:{
// 			break;
// 			}
// 		}
// 	}

	source_pool.insert(source_pool.begin(),9999);

	int num_infectious = (int) source_pool.size();

	int source_y = source_pool.at( gsl_rng_uniform_int (r_c, source_pool.size()) );

	vector<int> seq_proposed(n_base_CUPDATE);

	infected_source_modified.at(subject_proposed) = source_y; 	

	current_size_modified.at(subject_proposed) = 1;

	vector<int> nt_current_source_y ;
	vector<double> t_nt_current_source_y;

	vector<int> nt_modified_source_y ;
	vector<double> t_nt_modified_source_y ;


	//----
	switch(source_y==9999){
		case 0:{// add a 2nd infection

			double t_up = min(t_r_arg.at(source_y), t_max_CUPDATE);			
//			double t_low = t_i_arg.at(source_y);	
			double t_low = max(t_i_arg.at(source_y), t_up -t_back);	

			t_proposed= gsl_ran_flat(r_c, t_low, t_up); 
			t_e_modified.at(subject_proposed) = t_proposed;

			log_pr_backward = log(1.0/((double)(xi_EnIS_arg.size())+1.0));
			log_pr_forward = log(1.0/((double)(xi_U_arg.size()))) + log(1.0/((double)(num_infectious)));// the 2nd term corresponds to the random choice of a source (including background)

			log_pr_t_e_backward = 0.0;
			log_pr_t_e_forward = log(1.0/(t_up-t_low));

			current_size_modified.at(source_y) = current_size_arg.at(source_y) + 1;

			infecting_size_modified.at(source_y) = infecting_size_current_arg.at(source_y) + 1;

			vector<double> t_y(infecting_size_current_arg.at(source_y));
			for (int i=0;i<=(infecting_size_current_arg.at(source_y)-1);i++){
			t_y.at(i) = t_e_arg.at(infecting_list_current_arg[source_y][i]);
			}
			t_y.push_back(t_proposed);
			sort(t_y.begin(), t_y.end());
		
			int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_proposed));
			infecting_list_modified.at(source_y).insert(infecting_list_modified.at(source_y).begin()+rank_y, subject_proposed);

			//----
			nt_current_source_y = nt_current_arg.at(source_y);
			t_nt_current_source_y = t_nt_current_arg.at(source_y);
	
			nt_modified_source_y = nt_current_source_y;
			t_nt_modified_source_y = t_nt_current_source_y;
	
			t_nt_modified_source_y.push_back(t_proposed);
	
			sort( t_nt_modified_source_y.begin(),  t_nt_modified_source_y.end()); 
	
			int rank_source_y = distance(t_nt_modified_source_y.begin(), find(t_nt_modified_source_y.begin(),t_nt_modified_source_y.end(), t_proposed) );

			//----
	
			double t_past = t_nt_modified_source_y.at(rank_source_y-1);

			vector<int> nt_past_forward;
			nt_past_forward.assign(nt_current_source_y.begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE);

			switch(current_size_modified.at(source_y)>(rank_source_y+1)){

				case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y 

					double t_future = t_nt_modified_source_y.at(rank_source_y+1);
					vector<int> nt_future_forward;
					nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);

					seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

				break;
				}

				case 0:{//  inserted seq will be last seq at rank_source_y  in source_y

					double t_future = t_proposed;
					seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

				break;
				}
			}
		
			//----

			nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
	
			nt_modified.at(source_y) = nt_modified_source_y;
			t_nt_modified.at(source_y) = t_nt_modified_source_y;


		break;
		}
		case 1:{// add a background infection

			double t_up = t_max_CUPDATE;
			double t_low = max(0.0,t_up - t_back);				
//			double t_low = max(t_e_arg.at(index_arg.at(0)),t_up - t_back);	
	
			t_proposed= gsl_ran_flat(r_c, t_low, t_up); 
			t_e_modified.at(subject_proposed) = t_proposed;

			log_pr_backward = log(1.0/((double)(xi_EnIS_arg.size())+1.0));
			log_pr_forward = log(1.0/((double)(xi_U_arg.size()))) + log(1.0/((double)(num_infectious)));// the 2nd term corresponds to the random choice of a source (including background)

			log_pr_t_e_backward = 0.0;
			log_pr_t_e_forward = log(1.0/(t_up-t_low));

			//--------

// 			log_pr_seq_forward = n_base_CUPDATE*log(0.25);
// 
// 			for (int i=0; i<=(n_base_CUPDATE-1);i++){
// 			seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
// 			}
			sample_snull (con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE, r_c);//sample a seq for background
			log_pr_seq_forward = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);


		break;
		}
	}

	//-----------//

	vector<int> nt_modified_subject;
	vector<double> t_nt_modified_subject;

	nt_modified_subject = seq_proposed;

	t_nt_modified_subject.push_back(t_proposed); 

	nt_modified.at(subject_proposed) = nt_modified_subject;
	t_nt_modified.at(subject_proposed) = t_nt_modified_subject;

	//------------- deal with change of likelihood in f_s and f_null due to addition of an infection---------------//
	switch(source_y==9999){
		case 0:{// added 2nd infection

			log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y);
		
			lh_square_modified.log_f_S.at(source_y) = 0.0;

			for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){

			vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
			vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));

			lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
			}
	
			log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 


		break;
		}
		case 1:{// added background infection

			log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); // redudant as second term must be zero

//			lh_square_modified.log_f_Snull.at(subject_proposed) =  n_base_CUPDATE*log(0.25);
			lh_square_modified.log_f_Snull.at(subject_proposed) = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE);

			log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);


		break;
		}
	}


	//---------------------------------------------------------------------------------------------------------------------------//	

	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again
	
		if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
			switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
				case 1:{
					delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
				
				break;
				}
				case 0:{
					delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
				break;
				}
			}
		lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));

		}
	}

	//----------

	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;

	switch(source_y){
	
	case 9999:{ // by background
	lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
	break;
	}
	
	default :{ // not by background
	lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y);
	lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
	break;
	}
	
	}

	//--------
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes 1
	lh_square_modified.q_T.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	lh_square_modified.f_U.at(subject_proposed)=1.0;
	
	//--------
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  // f_EnI =/= 1 (since it is from xi_U, it must be not in xi_I before addtion)
	
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 

	//-----this part needed when do not update index
// 	  xi_E_minus_modified.push_back(subject_proposed); 

	  ////log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 
	  
// 	  lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
// 	  //lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
// 	  lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
// 	  lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
// 	  
// 	  log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 		

	//-----------------------
	
	switch (t_proposed<t_e_arg.at(index_arg.at(0))) { // 1 when t_proposed<index infection time
		
		case 1:{ // chosen one is the only index
	
	
			index_modified.clear();	
			index_modified.assign(1,subject_proposed);// replace index 
			//xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		
			for ( int i =0; i<= (int) (index_arg.size()-1); i++){
			xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
			}	
		
	// 		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 		
	// 		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.q_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.g_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.h_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.f_E.at(subject_proposed)=1.0;
			
		
			for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
			
			
			lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // secondary infection plays no role as the infectious times of others are always greater than the infection time of the original index
			lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
			lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
			lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
			}
	
		break;
		}
		
		case 0: {
			if((t_proposed==t_e_arg.at(index_arg.at(0))) ){ // add one more index
	
			index_modified.push_back(subject_proposed); // add index 
		
			}
	
			if((t_proposed>t_e_arg.at(index_arg.at(0))) ){
	
			xi_E_minus_modified.push_back(subject_proposed); 
	
			//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 
			
			lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
			//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
			lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
			lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 			
			}
	
		break;
		}
	
	}

		
	//---------------------------------//

	acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+ (log_pr_t_e_backward-log_pr_t_e_forward) +(log_pr_seq_backward-log_pr_seq_forward) ) );

	double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			lh_square_current_arg = lh_square_modified;
			delta_mat_current_arg = delta_mat_modified;
			log_lh_current_arg = log_lh_modified;
			t_e_arg= t_e_modified;
			index_arg = index_modified;
			xi_E_minus_arg = xi_E_minus_modified;
			xi_U_arg = xi_U_modified;
			xi_E_arg = xi_E_modified;
		
			xi_EnI_arg = xi_EnI_modified;
			xi_EnIS_arg = xi_EnIS_modified;
		
			infected_source_current_arg = infected_source_modified;
			current_size_arg = current_size_modified;
			infecting_list_current_arg = infecting_list_modified;
			infecting_size_current_arg = infecting_size_modified;
		
			nt_current_arg = nt_modified;
			t_nt_current_arg = t_nt_modified;
		
		break;
		}
		
		case 0: {
		break;
		}
	}


// 		ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("00test.txt")).c_str(),ios::app);
// 		if (uniform_rv<=acp_pr) myfile_mcmc_out <<action<< "," << acp_pr << "," << source_y <<","<< subject_proposed << endl;
// 		myfile_mcmc_out.close();		

break;
}// end of addition


} // end of switch on "action"


gsl_rng_free(r_c);

}

/*------------------------------------------------*/

void mcmc_UPDATE::t_e_add_del(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, int iter){

int subject_proposed ;
double t_proposed = 0.0;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;
vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;

//vector<int> infected_source_modified = infected_source_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed

//---------------------


//int action = gsl_rng_uniform_int (r_c, 3);// a int random number in [0,2] will be drawn : 0= replacement of an exsiting infection time; 1=deletion of an infection into gp of xi_EnI (choose from the gp of xi_U); 2= addition of an infection from gp of xi_EnI (add back  to xi_U)

// if (xi_EnI_arg.empty()==1){ // return 1 if xi_EnI is empty
// action = gsl_rng_uniform_int (r_c, 2); // only replacement OR addition is possible
// switch(action){
// case 0:{
// action = 0; //replacement
// break;
// }
// case 1:{
// action = 2; //addition
// break;
// }
// }
// }
// 
// if (xi_U_arg.empty()==1){ // return 1 if xi_U is empty
// action = gsl_rng_uniform_int (r_c, 2);// only replacement OR deletion is possible
// }



double P [2] ={2.0,2.0};
gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (2,P);
int action = gsl_ran_discrete (r_c, g); // a int random number in [0,1] will be drawn according to the weights in P :  0=deletion of an infection into gp of xi_EnI (choose from the gp of xi_U); 1= addition of an infection from gp of xi_EnI (add back  to xi_U)
gsl_ran_discrete_free (g);

if (xi_EnI_arg.empty()==1){ // return 1 if xi_EnI is empty
action = 1; // only addtion possible
}

if (xi_U_arg.empty()==1){ // return 1 if xi_U is empty
action = 0; // only deletion is possible
}


// 	ofstream myfile1_mcmc_out; 
// 	myfile1_mcmc_out.open((string(path4)+string("action.txt")).c_str(),ios::app);
// 	myfile1_mcmc_out << action << endl;
// 	myfile1_mcmc_out.close();


switch (action) {

case 0:{ // deletion

	subject_proposed = xi_EnI_arg.at(gsl_rng_uniform_int (r_c, xi_EnI_arg.size()));
	
// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("subject_proposed_del_before.txt")).c_str(),ios::app);
// 	myfile_mcmc_out <<endl;
// 	myfile_mcmc_out << "sub_del"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 	myfile_mcmc_out << subject_proposed  <<"," <<lh_square_current_arg.kt_sum_U.at(subject_proposed)<<"," <<lh_square_current_arg.f_U.at(subject_proposed)<< "," <<lh_square_current_arg.f_E.at(subject_proposed)<< "," <<lh_square_current_arg.f_EnI.at(subject_proposed)<<endl;
// 	myfile_mcmc_out.close();
	
	xi_EnI_modified.erase(find(xi_EnI_modified.begin(),xi_EnI_modified.end(),subject_proposed));
	xi_E_modified.erase(find(xi_E_modified.begin(),xi_E_modified.end(),subject_proposed));
	xi_U_modified.push_back(subject_proposed);
	
	//infected_source_modified.at(subject_proposed) = -99; // this subject become NO infected source

	//----
	
	t_e_modified.at(subject_proposed) = unassigned_time_CUPDATE;

	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E 
	
		//if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
		switch (t_r_arg.at(xi_I_arg.at(j))>=t_max_CUPDATE) {
			case 1:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_max_CUPDATE - t_i_arg.at(xi_I_arg.at(j));
			//lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)] ;
			break;
			}
			case 0:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
			break;
			}
		}//end switch
		
		lh_square_modified.kt_sum_U.at(subject_proposed) = lh_square_modified.kt_sum_U.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
		//}
	}
	
	
	//-----
	
	//log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  
	
	lh_square_modified.q_T.at(subject_proposed) = para_current_arg.alpha*t_max_CUPDATE + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_U.at(subject_proposed);
	lh_square_modified.f_U.at(subject_proposed) = 1.0 - gsl_cdf_exponential_P(lh_square_modified.q_T.at(subject_proposed),1.0);
	
	log_lh_modified = log_lh_modified  + log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes =/=1
	
	//----
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  //f_EnI becomes 1
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0;
	
	//---
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));  //f_E becomes 1
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.q_E.at(subject_proposed)=0.0;
	lh_square_modified.g_E.at(subject_proposed)=1.0;
	lh_square_modified.h_E.at(subject_proposed)=1.0;
	lh_square_modified.f_E.at(subject_proposed)=1.0;
	
		switch(find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) {// returns 1 if it was NOT an index originally
		case 0:{// was an index
	
			index_modified.clear();
		
			int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
			double min_t = t_e_modified.at(first_min); // the minimum time of exposure
			
			int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time
			
				switch (num_min>1) {
				case 1: {
				index_modified.reserve(n_CUPDATE);	
				for (int i=0; i<=(n_CUPDATE-1);i++){
				if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
				}
				break;
				}
				case 0:{
				index_modified.assign(1,first_min);
				break;
				}
				
				}
	
			//xi_E_minus_modified = xi_E_arg;
		
			for (int i=0;i<= (int) (index_modified.size()-1); i++){
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));
			log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i)));
			lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
			lh_square_modified.q_E.at(index_modified.at(i))=0.0;
			lh_square_modified.g_E.at(index_modified.at(i))=1.0;
			lh_square_modified.h_E.at(index_modified.at(i))=1.0;
			lh_square_modified.f_E.at(index_modified.at(i))=1.0;	
			}
					
		break;
		}	
	
		case 1:{// was not an index
			xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		break;
		}
		}

		//ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_del_after.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<endl;
// 		myfile_mcmc_out << "sub_del"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 		myfile_mcmc_out << subject_proposed  <<"," <<lh_square_modified.kt_sum_U.at(subject_proposed)<<"," <<lh_square_modified.f_U.at(subject_proposed)<< "," <<lh_square_modified.f_E.at(subject_proposed)<< "," <<lh_square_modified.f_EnI.at(subject_proposed)<<endl;
// 		myfile_mcmc_out.close();

		//-----------------------------------------//
				
		acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(2.0/2.0)*(double)(xi_EnI_arg.size())*(1/t_max_CUPDATE)/((double)(xi_U_arg.size())+1.0) );

// 		myfile_mcmc_out.open((string(path4)+string("acp_del.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "acp_pr" << "," <<"log_lh_modified" << "," << "log_lh_unmodified"  << endl;
// 		myfile_mcmc_out <<  acp_pr << "," <<log_lh_modified << "," << log_lh_current_arg  << endl;
// 		myfile_mcmc_out.close();

		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		lh_square_current_arg = lh_square_modified;
		delta_mat_current_arg = delta_mat_modified;
		log_lh_current_arg = log_lh_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;
		xi_U_arg = xi_U_modified;
		xi_E_arg = xi_E_modified;
		xi_EnI_arg = xi_EnI_modified;
		//infected_souce_current_arg = infected_source_modified;
		break;
		}
		
		case 0: {
		break;
		}
		}

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_del.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "log_lh_current_arg"  << "," <<"log_lh_modified" << endl;
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();


break;
}// end of deletion

//------------------------------

case 1:{ // addition
	
	subject_proposed = xi_U_arg.at(gsl_rng_uniform_int (r_c, xi_U_arg.size()));
	
// 	ofstream myfile_mcmc_out; 
// 	myfile_mcmc_out.open((string(path4)+string("subject_proposed_add_before.txt")).c_str(),ios::app);
// 	myfile_mcmc_out <<endl;
// 	myfile_mcmc_out << "sub_add"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 	myfile_mcmc_out << subject_proposed  <<"," <<lh_square_current_arg.kt_sum_U.at(subject_proposed)<<"," <<lh_square_current_arg.f_U.at(subject_proposed)<< "," <<lh_square_current_arg.f_E.at(subject_proposed)<< "," <<lh_square_current_arg.f_EnI.at(subject_proposed)<<endl;
// 	myfile_mcmc_out.close();
	
	xi_U_modified.erase(find(xi_U_modified.begin(),xi_U_modified.end(),subject_proposed)); // alter the gps xi_U, xi_E, xi_EnI
	xi_EnI_modified.push_back(subject_proposed);
	xi_E_modified.push_back(subject_proposed);
	
	//----
	
	t_proposed= gsl_ran_flat(r_c, 0.0, t_max_CUPDATE); // propose an infection time for the added infection
	
	t_e_modified.at(subject_proposed) = t_proposed;
	
	lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	
	for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E 
	
		if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
		
		switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
			case 1:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
			lh_square_modified.k_sum_E.at(subject_proposed) =  lh_square_modified.k_sum_E.at(subject_proposed) + kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));
			break;
			}
			case 0:{
			delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
			break;
			}
		}//end switch
		
		lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j)) ;
		}
	}
	
	//--------
	
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_U.at(subject_proposed));  //f_U becomes 1
	lh_square_modified.q_T.at(subject_proposed)=0.0;
	lh_square_modified.kt_sum_U.at(subject_proposed)=0.0;
	lh_square_modified.f_U.at(subject_proposed)=1.0;
	
	//--------
	log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed));  // f_EnI =/= 1 (since it is from xi_U, it must be not in xi_I before addtion)
	
	lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
	
	log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 
	//-----
	
		switch (t_proposed<t_e_arg.at(index_arg.at(0))) { // 1 when t_proposed<index infection time
		
		case 1:{ // chosen one is the only index
	
	
			index_modified.clear();	
			index_modified.assign(1,subject_proposed);// replace index 
			//xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
		
			for ( int i =0; i<= (int) (index_arg.size()-1); i++){
			xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
			}	
		
	// 		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 		
	// 		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.q_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.g_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.h_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.f_E.at(subject_proposed)=1.0;
			
		
			for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
			
			//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
			
			lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // secondary infection plays no role as the infectious times of others are always greater than the infection time of the original index
			lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
			lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
			lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
			}
	
		break;
		}
		
		case 0: {
			if((t_proposed==t_e_arg.at(index_arg.at(0))) ){ // add one more index
	
			index_modified.push_back(subject_proposed); // add index 
	
	// 		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
	// 		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.q_E.at(subject_proposed)=0.0;
	// 		lh_square_modified.g_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.h_E.at(subject_proposed)=1.0;
	// 		lh_square_modified.f_E.at(subject_proposed)=1.0;		
	
			// //xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
		
			}
	
			if((t_proposed>t_e_arg.at(index_arg.at(0))) ){
	
			xi_E_minus_modified.push_back(subject_proposed); 
	
			//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 
			
			lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
			lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
			lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
			lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
			
			log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 			
			}
	
		break;
		};
	
	}

		//ofstream myfile_mcmc_out; 
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_add_after.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<endl;
// 		myfile_mcmc_out << "sub_add"  <<"," <<"kt_sum_U"<<"," <<"f_U"<< "," <<"f_E"<< ","<<"f_EnI"<<endl;
// 		myfile_mcmc_out << subject_proposed  <<"," <<lh_square_modified.kt_sum_U.at(subject_proposed)<<"," <<lh_square_modified.f_U.at(subject_proposed)<< "," <<lh_square_modified.f_E.at(subject_proposed)<< "," <<lh_square_modified.f_EnI.at(subject_proposed)<<endl;
// 		myfile_mcmc_out.close();
		
		//-----------------------------------------//
		
		acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*(2.0/2.0)*(double)(xi_U_arg.size())*(t_max_CUPDATE)/((double)(xi_EnI_arg.size())+1.0) );

// 		myfile_mcmc_out.open((string(path4)+string("acp_add.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "acp_pr" << "," <<"log_lh_modified" << "," << "log_lh_unmodified_arg"  << endl;
// 		myfile_mcmc_out <<  acp_pr << "," <<log_lh_modified << "," << log_lh_current_arg  << endl;
// 		myfile_mcmc_out.close();

		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		lh_square_current_arg = lh_square_modified;
		delta_mat_current_arg = delta_mat_modified;
		log_lh_current_arg = log_lh_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;
		xi_U_arg = xi_U_modified;
		xi_E_arg = xi_E_modified;
		xi_EnI_arg = xi_EnI_modified;
		break;
		}
		
		case 0: {
		break;
		}
		}

// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_add.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << endl;
// 		myfile_mcmc_out <<  "log_lh_current_arg"  << "," <<"log_lh_modified" << endl;
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << endl;
// 		myfile_mcmc_out.close();

break;
}// end of addition


} // end of switch on "action"


gsl_rng_free(r_c);

}

//-----------------------------------------------
void residual_sub_nucle (vector<double>& u_imputed_sub, const vector<int>& seq_1, const vector<int>& seq_2, const double& t_1, const double& t_2, double mu_1, double mu_2, const int& n_base, base_generator_type& seed){

//base_generator_type seed(time(0));

double dt = t_2 - t_1;

double p_1 = 0.25 + 0.25*exp(-4.0*mu_2*dt) + 0.5*exp(-2.0*(mu_1+mu_2)*dt) ; // pr of a base not changing
double p_2 = 0.25 + 0.25*exp(-4.0*mu_2*dt) - 0.5*exp(-2.0*(mu_1+mu_2)*dt); // pr of a transition of a base
double p_3 = 2.0*(0.25 - 0.25*exp(-4.0*mu_2*dt));  // pr of a transversion (two possible events)
//double p_4 = 1.0*(0.25 - 0.25*exp(-4.0*mu_2*dt));

vector<double> pdf(4);
vector<double> cdf(4);

	for (int i=0;i<=(n_base-1);i++){

		
		pdf.at(0) = 0.0;
		pdf.at(1) = p_1;
		pdf.at(2) = p_2;
		pdf.at(3) = p_3;
		
		
		int type=-1;
		int link =-1;
		double p_link =0.0;


		switch( abs( seq_2.at(i)-seq_1.at(i) ) ){
			case 0:{
				type=1;
			break;
			}
			case 1:{
				switch( ((seq_1.at(i)==2) & (seq_2.at(i)==3)) | ((seq_1.at(i)==3) & (seq_2.at(i)==2)) ){
					case 1:{ // transversion
						type = 3;
					break;
					}
					case 0:{ // transition
						type = 2;
					break;
					}
				}
			break;
			}
			case 2:{
				type = 3;
			break;
			}
			case 3:{
				type = 3;
			break;
			}
		}
			
		p_link= pdf.at(type);
	
		sort(pdf.begin(), pdf.end());
		
		link = distance(  pdf.begin(), find(pdf.begin(),pdf.end(),p_link) );
		
		//---
		int num_dup=1; 
		vector<double> dup_test = pdf;
		dup_test.erase(dup_test.begin()+link);
	
		switch(find(dup_test.begin(),dup_test.end(),p_link)==dup_test.end()) { // return 1 when no duplicate 
		case 0:{ // with dup
		num_dup = count(pdf.begin(),pdf.end(),p_link);
		break;
		}
		case 1:{
		num_dup=1;
		break;
		}
		}
		//--

		
		cdf.at(0) = 0.0;
	
		for (int k=1; k<=(int)(cdf.size()-1);k++){
			cdf.at(k) = pdf.at(k) + cdf.at(k-1);
		}

		switch(link==0){
			case 0:{
//				u_imputed_sub.push_back(gsl_ran_flat(r_c, cdf.at(link-1),  cdf.at(link-1)+num_dup*p_link) );

				boost::uniform_real<> uni_dist(cdf.at(link-1), cdf.at(link-1)+ ((double)num_dup)*p_link);
				boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(seed, uni_dist);
				double u_i = uni();

// 				RNGScope scope;
// 				double u_i= R::runif(cdf.at(link-1), cdf.at(link-1)+num_dup*p_link);

				u_imputed_sub.push_back(u_i);



			break;
			}
			case 1:{
				u_imputed_sub.push_back(0.0);
			break;
			}	
		}

//      ofstream myfile_out; 

// 	myfile_out.open((string(path4)+string("residual_nucle_cdf.txt")).c_str(),ios::app);
// 	myfile_out <<  cdf.at(0) <<"," << cdf.at(1) <<","<< cdf.at(2) << "," << cdf.at(3) << ","<< pdf.at(0) + pdf.at(1) + pdf.at(2) + pdf.at(3)<<  endl;
// 	myfile_out.close();
// 	
// 	myfile_out.open((string(path4)+string("residual_nucle_pdf.txt")).c_str(),ios::app);
// 	myfile_out <<  pdf.at(0) <<"," << pdf.at(1) <<","<< pdf.at(2) << "," << pdf.at(3) << "," <<type<<","<< link << ","<< p_link<< ","<< p_1 << "," << num_dup << endl;
// 	myfile_out.close();


	}




}

//-----------------------------------------------
void residual_nucle (const para_key& para_current_arg, const para_aux& para_other_arg, vector<int>& xi_E_arg,const vector<int>& current_size_arg, const vec2int& nt_current_arg, const vec2d& t_nt_current_arg, int iter){

ofstream myfile_out; 

// const gsl_rng_type* T_c= gsl_rng_ranlux389;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,iter); // set a seed

base_generator_type seed(iter);

vector<long double> u_imputed;

for (int i=0;i<=((int)xi_E_arg.size()-1);i++){

	switch(current_size_arg.at(xi_E_arg.at(i))>1){

		case 1:{

			for (int j=0;j<=(current_size_arg.at(xi_E_arg.at(i))-2);j++){
	
				vector<int> seq_1(nt_current_arg.at(xi_E_arg.at(i)).begin()+j*(para_other_arg.n_base), nt_current_arg.at(xi_E_arg.at(i)).begin()+(j+1)*(para_other_arg.n_base));
				vector<int> seq_2(nt_current_arg.at(xi_E_arg.at(i)).begin()+(j+1)*(para_other_arg.n_base), nt_current_arg.at(xi_E_arg.at(i)).begin()+(j+2)*(para_other_arg.n_base));
				
				double t_1 = t_nt_current_arg[xi_E_arg.at(i)][j];
				double t_2 = t_nt_current_arg[xi_E_arg.at(i)][j+1];

				vector<double> u_imputed_sub;

				residual_sub_nucle (u_imputed_sub, seq_1, seq_2, t_1, t_2, para_current_arg.mu_1, para_current_arg.mu_2, para_other_arg.n_base, seed);

				u_imputed.insert(u_imputed.begin(), u_imputed_sub.begin(), u_imputed_sub.end());
			}

		break;
		}

		case 0:{
			// do nothing
		break;
		}

	}
}

// myfile_out.open((string(path4)+string("residual_nucle.txt")).c_str(),ios::app);
// myfile_out << endl;
// for (int i=0; i<=((int)u_imputed.size()-1);i++){
// if (i!= (int) u_imputed.size()-1) myfile_out << u_imputed.at(i) << ",";
// if (i== (int) u_imputed.size()-1) myfile_out << u_imputed.at(i);
// }
// myfile_out.close();

myfile_out.open((string(path4)+string("residual_nucle.txt")).c_str(),ios::app);
for (int i=0; i<=((int)u_imputed.size()-1);i++){
myfile_out << u_imputed.at(i) << endl;

// 	if( (double) u_imputed.at(i)==1.0){
// 	myfile_out.open((string(path4)+string("00test_2.txt")).c_str(),ios::app);
// 	myfile_out <<  u_imputed.at(i) <<  endl;
// 	myfile_out.close();
// 	}

}
myfile_out.close();

myfile_out.open((string(path4)+string("residual_nucle_num_list.txt")).c_str(),ios::app);
myfile_out <<  u_imputed.size() << endl;
myfile_out.close();


//gsl_rng_free(r_c);

}

//-----------------------------------------------

void residual_kernel_path(const vector< vector<double> >& kernel_mat_current_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg,  const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const para_key& para_current_arg, const para_aux& para_other, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, int iter){

ofstream myfile_out; 

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed


vector <double> u_imputed(xi_E_minus_arg.size()); // the imputed residuals for all infected, exclude index


for (int i=0; i<=((int)xi_E_minus_arg.size()-1);i++){
//for (int i=0; i<=((int)1-1);i++){

	vector<double> ic_link; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/
	//vector<double> pr_link; // probabilities of the infection links between infectious/primary infection and the new infected
	
	vector<double> ic; // instantaneous infective challenge of the infection links between infectious/primary infection and the new infected/susceptibles
	vector<double> pr; // probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	vector<double> cdf; // cumulative probabilities of the infection links between infectious/primary infection and the new infected/susceptibles
	
	double total_ic_link; // the total instantaneous infective challenge to the new infected
	double total_ic; // the total instantaneous infective challenge to the new infected/susceptibles
	
	//int link_imputed; // link imputed
	double ic_imputed; // ic of the imputed link
	int rank_imputed; // the rank of p_imputed among all links; 0 means first element

	ic_link.reserve(xi_I_arg.size()+2); // 2 more positions reserved: one for primary infection and one for faciliating computation

	
	ic_link.push_back(0.0); // 1st element in ic
	ic_link.push_back(para_current_arg.alpha); // 2nd..

	//-----

	total_ic_link = 0.0; 
	for (int j=0; j<=((int)xi_I_arg.size()-1);j++){
	if ((t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_minus_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i)))){
	ic_link.push_back(para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*kernel_mat_current_arg[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j)));
	total_ic_link = total_ic_link + stb_arg.at(xi_E_minus_arg.at(i))*kernel_mat_current_arg[xi_E_minus_arg.at(i)][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
	}
	}
	total_ic_link = para_current_arg.alpha + para_current_arg.beta*total_ic_link;


// 	pr_link.resize(ic_link.size());
// 
// 	pr_link.at(0) = 0.0;
// 	pr_link.at(1) = para_current_arg.alpha/total_ic_link; //pr it was infected by primary infection
// 	
// 	for (int k=2; k<=(int)(pr_link.size()-1);k++){
// 	pr_link.at(k) = ic_link.at(k)/total_ic_link;
// 	}
// 
// 	double *P=&pr_link.at(0);
// 	gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)pr_link.size(),P);
// 	link_imputed= gsl_ran_discrete (r_c, g); // the link imputed
// 	gsl_ran_discrete_free (g);
// 
// 	ic_imputed = ic_link.at(link_imputed);

	//-----
	int source = infected_source_current_arg.at(xi_E_minus_arg.at(i)); // do not need to impute link given transmission path

	switch(source==9999){
		case 0:{
			ic_imputed = para_current_arg.beta*stb_arg.at(xi_E_minus_arg.at(i))*kernel_mat_current_arg[xi_E_minus_arg.at(i)][source]/ norm_const_current_arg.at(source);
		break;
		}
		case 1:{
			ic_imputed = para_current_arg.alpha;
		break;
		}
	}
	//-----
	
	ic.reserve(ic_link.size()+(xi_I_arg.size()+1)*para_other.n);

	total_ic = 0.0; // sum of ic of links from the bunch of infectious to ALL susceptible subject before t_e_arg.at(xi_E_minus_arg.at(i))

	for (int iu=0; iu<=((int)para_other.n-1);iu++){

		double total_ic_iu =0.0; // sum if ic of links from the bunch of infectious to the particular susceptible subject

		switch(t_e_arg.at(iu)>t_e_arg.at(xi_E_minus_arg.at(i))) { // return 1 if the subject is susceptible before t_e_arg.at(xi_E_minus_arg.at(i))note: this excludes the links connected to xi_E_minus.at(i), where the infection actually happened
		case 1:{
		//ic.push_back(0.0); 
		ic.push_back(para_current_arg.alpha); 
	
		for (int j=0; j<=((int)xi_I_arg.size()-1);j++){
		if ((t_i_arg.at(xi_I_arg.at(j))<t_e_arg.at(xi_E_minus_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(j))>=t_e_arg.at(xi_E_minus_arg.at(i)))){ // this ensures we are using same bunch of infectious sources
		ic.push_back(para_current_arg.beta*stb_arg.at(iu)*kernel_mat_current_arg[iu][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j)));
		total_ic_iu = total_ic_iu + stb_arg.at(iu)*kernel_mat_current_arg[iu][xi_I_arg.at(j)]/ norm_const_current_arg.at(xi_I_arg.at(j));
		}
		}
		total_ic_iu = para_current_arg.alpha + para_current_arg.beta*total_ic_iu; 
	
		break;
		}
	
		case 0:{
		break;
		}
		}

		total_ic = total_ic + total_ic_iu;

	}
	
	total_ic = total_ic + total_ic_link; // add the ic from the infectious to infected links

	ic.insert(ic.begin(),ic_link.begin(),ic_link.end()); // combine ic_link and ic

	pr.resize(ic.size());

	for (int k=0; k<=((int)pr.size()-1);k++){
	pr.at(k) = ic.at(k)/total_ic;
	}

	sort(pr.begin(),pr.end());

	double p_imputed = ic_imputed / total_ic;
	
	rank_imputed = distance(  pr.begin(), find(pr.begin(),pr.end(),p_imputed) );


	//---------//

	int num_dup=1; // number of prs same to p_imputed (including p_imputed itself)
	vector<double> dup_test = pr;
	dup_test.erase(dup_test.begin()+rank_imputed);

	switch(find(dup_test.begin(),dup_test.end(),p_imputed)==dup_test.end()) { // return 1 when no duplicate 
	case 0:{ // with dup
	num_dup = count(pr.begin(),pr.end(),p_imputed);
	break;
	}
	case 1:{
	num_dup=1;
	break;
	}
	}

	//---

// 	myfile_out.open((string(path4)+string("t_e_consider.txt")).c_str(),ios::app);
// 	if (find(dup_test.begin(),dup_test.end(),p_imputed)!=dup_test.end()){
// 	myfile_out << t_e_arg.at(xi_E_minus_arg.at(i)) << endl;
// 	}
// 	myfile_out.close();
// 
// 	vector<double> t_e_minus = t_e_arg;
// 	t_e_minus.erase(min_element(t_e_minus.begin(),t_e_minus.end()));
// 	double min_t_e_minus = *min_element(t_e_minus.begin(),t_e_minus.end());
// 	myfile_out.open((string(path4)+string("min_t_e_minus.txt")).c_str(),ios::app);
// 	if (find(dup_test.begin(),dup_test.end(),p_imputed)!=dup_test.end()){
// 	myfile_out << min_t_e_minus << endl;
// 	}
// 	myfile_out.close();

// 	myfile_out.open((string(path4)+string("dup_test.txt")).c_str(),ios::app);
// 	int num_excl = 0; // number that have been infected before the t_e considering
// 	for (int kk=0;kk<=(int)(xi_E_minus_arg.size()-1);kk++){
// 	if (t_e_arg.at(kk)<t_e_arg.at(xi_E_minus_arg.at(i))) num_excl =num_excl + 1;
// 	}
// 	int num_infectious=0; // number of infectious for the t_e considering
// 	for (int ii=0;ii<=(int)(xi_I_arg.size()-1);ii++){
// 	if ( (t_i_arg.at(xi_I_arg.at(ii))<t_e_arg.at(xi_E_minus_arg.at(i))) & (t_r_arg.at(xi_I_arg.at(ii))>=t_e_arg.at(xi_E_minus_arg.at(i))) ) num_infectious =num_infectious+ 1;
// 	}	
// 	vector<double> ic_sort =ic;
// 	sort(ic_sort.begin(),ic_sort.end());
// 	if (find(dup_test.begin(),dup_test.end(),p_imputed)!=dup_test.end()){
// 		if (num_infectious==0 ) myfile_out << p_imputed << "," << i << "," << 1.0/(double) (para_other.n - num_excl)<< endl; // the imputed link is from primary (as no infectious)
// 	if (num_infectious>=1) myfile_out << p_imputed << "," << i << "," << ic_sort.at(rank_imputed)/ total_ic<<"," << para_current_arg.alpha/total_ic <<endl; // if the imputed link is from primary: p_imputed = last entry ;if the imputed link is from 2nd : p_imputed = 2nd last entry;  if 2 last entries are equal, it is from primary
// 	}
// 	myfile_out.close();
// 
// 
// 	myfile_out.open((string(path4)+string("dup_num.txt")).c_str(),ios::app);
// 	if (find(dup_test.begin(),dup_test.end(),p_imputed)!=dup_test.end()){
// 	myfile_out <<  num_dup << "," << p_imputed << "," << i << "," << xi_E_minus_arg.at(i) << endl;
// 	}
// 	myfile_out.close();


	//----

	cdf.resize(pr.size());

	cdf.at(0) = 0.0;
	for (int k=1; k<=(int)(cdf.size()-1);k++){
	cdf.at(k) = pr.at(k) + cdf.at(k-1);
	}

	double u_i=0.0; // the u_i to be imputed
	switch(rank_imputed==0){
	case 1:{
	u_i = 0.0;
	break;
	}
	case 0:{
/*	u_i = gsl_ran_flat(r_c, cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+p_imputed*/
	u_i = gsl_ran_flat(r_c, cdf.at(rank_imputed-1),  cdf.at(rank_imputed-1)+num_dup*p_imputed); // a uniform rv drawn between cdf.at(rank_imputed-1), and cdf.at(rank_imputed-1)+num_dup*p_imputed
	break;
	}
	}

	u_imputed.at(i) =u_i;

	//-------------------
/*
	myfile_out.open((string(path4)+string("size_ic_pr.txt")).c_str(),ios::app);
	myfile_out <<  ic.size() << "," << pr.size() <<endl;
	myfile_out.close();

	myfile_out.open((string(path4)+string("ic.txt")).c_str(),ios::app);
	for (int i=0;i<=(int) (ic.size()-1);i++){
	myfile_out << ic.at(i) << endl; 
	}	
	myfile_out.close();
	
	double total_pr =0.0;
	myfile_out.open((string(path4)+string("pr.txt")).c_str(),ios::app);
	for (int i=0;i<=(int) (pr.size()-1);i++){
	total_pr =total_pr + pr.at(i);
	myfile_out << pr.at(i) << endl; 
	}	
	myfile_out.close();

	myfile_out.open((string(path4)+string("cdf.txt")).c_str(),ios::app);
	for (int i=0;i<=(int) (cdf.size()-1);i++){
	myfile_out << cdf.at(i) << endl; 
	}	
	myfile_out.close();

	myfile_out.open((string(path4)+string("total_ic.txt")).c_str(),ios::app);
	myfile_out << total_ic << endl; 
	myfile_out.close();	


	myfile_out.open((string(path4)+string("link_imputed.txt")).c_str(),ios::app);
	myfile_out <<  link_imputed  << endl; 
	myfile_out.close();

	myfile_out.open((string(path4)+string("rank_imputed.txt")).c_str(),ios::app);
	myfile_out <<  rank_imputed  << endl; 
	myfile_out.close();

	myfile_out.open((string(path4)+string("p_imputed.txt")).c_str(),ios::app);
	myfile_out <<  p_imputed  << endl; 
	myfile_out.close();*/

	//------------------------


	//pr_link.clear();
	ic_link.clear();

	pr.clear();
	ic.clear();
	cdf.clear();

}

	myfile_out.open((string(path4)+string("residual_kernel.txt")).c_str(),ios::app);
	myfile_out << endl;
	for (int i=0; i<=((int)xi_E_minus_arg.size()-1);i++){
	if (i!= (int) xi_E_minus_arg.size()-1) myfile_out << u_imputed.at(i) << ",";
	if (i== (int) xi_E_minus_arg.size()-1) myfile_out << u_imputed.at(i);
	}
	myfile_out.close();

// 	myfile_out.open((string(path4)+string("size_residual.txt")).c_str(),ios::app);
// 	myfile_out <<  u_imputed.size() << endl;
// 	myfile_out.close();

gsl_rng_free(r_c);

}
//-----------------------------------



void residual_lat( const vector<int>& xi_E_arg, const vector<int>& xi_I_arg,  const vector<int>& xi_EnI_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const para_key& para_current_arg, const para_aux& para_other,  int iter){

ofstream myfile_out; 

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed

vector <double> u_imputed; // the imputed residuals for all infected
u_imputed.reserve(xi_E_arg.size());

for (int i=0; i<=((int)xi_I_arg.size()-1);i++){
double u_i = exp(log(1.0-func_latent_cdf( t_i_arg.at(xi_I_arg.at(i))-t_e_arg.at(xi_I_arg.at(i)), para_current_arg.mu_lat, para_current_arg.var_lat )));
u_imputed.push_back(u_i);
}

for (int i=0; i<=((int)xi_EnI_arg.size()-1);i++){
double u_i = gsl_ran_flat( r_c, 0.0, 1.0-func_latent_cdf( para_other.t_max-t_e_arg.at(xi_EnI_arg.at(i)), para_current_arg.mu_lat, para_current_arg.var_lat) );
u_imputed.push_back(u_i);
}

// myfile_out.open((string(path4)+string("S_xi_I_latent.txt")).c_str(),ios::app);
// myfile_out << endl;
// for (int i=0; i<=((int)xi_I_arg.size()-1);i++){
// double s_i = 1.0 - func_latent_cdf( t_i_arg.at(xi_I_arg.at(i))-t_e_arg.at(xi_I_arg.at(i)), para_current_arg.mu_lat, para_current_arg.var_lat );
// myfile_out << s_i  << "," ;
// }
// myfile_out.close();
// 
// myfile_out.open((string(path4)+string("S_xi_EnI_latent.txt")).c_str(),ios::app);
// myfile_out << endl;
// for (int i=0; i<=((int)xi_EnI_arg.size()-1);i++){
// double s_i = gsl_ran_flat( r_c, 0.0, 1.0-func_latent_cdf( para_other.t_max-t_e_arg.at(xi_EnI_arg.at(i)), para_current_arg.mu_lat, para_current_arg.var_lat) );
// myfile_out << s_i  << exp(log(s_i)) << "," ;
// }
// myfile_out.close();


myfile_out.open((string(path4)+string("residual_latent.txt")).c_str(),ios::app);
myfile_out << endl;
for (int i=0; i<=((int)xi_E_arg.size()-1);i++){
if (i!= (int) xi_E_arg.size()-1) myfile_out << u_imputed.at(i) << ",";
if (i== (int) xi_E_arg.size()-1) myfile_out << u_imputed.at(i);
}
myfile_out.close();


gsl_rng_free(r_c);

}

//---------------------

void residual_gnl(const vector< vector<double> >& kernel_mat_current_arg, const vector<int>& xi_E_arg, const vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg,  const vector<double>& t_r_arg, const vector<double>& t_i_arg, const vector<double>& t_e_arg, const para_key& para_current_arg, const para_aux& para_other, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, int iter){

ofstream myfile_out; 

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter*1000); // set a seed

vector <int> xi_E_eff = xi_E_minus_arg; //exclude index
//vector <int> xi_E_eff = xi_E_arg; // include index ( when primary infection assumed presented)

vector <double> q_imputed(xi_E_eff.size() + 1); // the imputed residuals for all infected and last unobserved infection
vector <double> u_imputed(xi_E_eff.size() + 1); // the imputed (transformed) residuals for all infected and last unobserved infection

vector <double> qt_acml(xi_E_eff.size() +1); // the accumulated challenges at each time of infection considered and at t_max

vector<double> t_e_sorted(xi_E_eff.size()); // include (sorted) valid infection times, and would have t_max at the end
for  ( int i=0; i<=((int)t_e_sorted.size()-1);i++){
t_e_sorted.at(i) = t_e_arg.at(xi_E_eff.at(i));
}

sort(t_e_sorted.begin(),t_e_sorted.end());

t_e_sorted.push_back(para_other.t_max);


for (int i=0; i<=((int)t_e_sorted.size()-1);i++){

	vector<double> ic(para_other.n); // the vector contains the infective challange to the individuals in population (would be re-initialized for each infection event being considered

	double t_now = t_e_sorted.at(i);
	
	for (int j=0; j<=(para_other.n-1);j++){ //loop over all individuals and compute the accumulated pressure, finally add them up

		switch(t_e_arg.at(j)>=t_now){ //return 1 if  j corresponds to a susceptible before t_now (this would include the individual whose t_e=t_now)
	
		case 1:{ // for susceptible j before t_now

			ic.at(j) = para_current_arg.alpha*t_now;

			for (int k=0; k<=((int)t_i_arg.size()-1);k++){
			if (t_i_arg.at(k)<t_now){ // if k is once infectious before t_now

				switch(t_r_arg.at(k)<t_now){ // reurn 1 if k has been recovered before t_now
				case 1:{
				ic.at(j) = ic.at(j) + para_current_arg.beta*stb_arg.at(j)*kernel_mat_current_arg[j][k]*(t_r_arg.at(k) - t_i_arg.at(k))/norm_const_current_arg.at(k);
				break;
				}
				case 0:{
				ic.at(j) = ic.at(j) + para_current_arg.beta*stb_arg.at(j)*kernel_mat_current_arg[j][k]*(t_now - t_i_arg.at(k))/norm_const_current_arg.at(k);
				break;
				}			
				}		
			}
			}
		break;
		}
	
		case 0:{ // for infected j before t_now

			ic.at(j) = para_current_arg.alpha*t_e_arg.at(j);

			for (int k=0; k<=((int)t_i_arg.size()-1);k++){
			if (t_i_arg.at(k)<t_e_arg.at(j)){ // if k is once infectious before t_e_arg.at(j)

				switch(t_r_arg.at(k)<t_e_arg.at(j)){ // reurn 1 if k has been recovered before t_e_arg.at(j)
				case 1:{
				ic.at(j) = ic.at(j) + para_current_arg.beta*stb_arg.at(j)*kernel_mat_current_arg[j][k]*(t_r_arg.at(k) - t_i_arg.at(k))/norm_const_current_arg.at(k);
				break;
				}
				case 0:{
				ic.at(j) = ic.at(j) + para_current_arg.beta*stb_arg.at(j)*kernel_mat_current_arg[j][k]*(t_e_arg.at(j)- t_i_arg.at(k))/norm_const_current_arg.at(k);
				break;
				}			
				}		
			}
			}
		
		break;
		}
		}	

	}	

	qt_acml.at(i) = accumulate(ic.begin(),ic.end(),0.0);

	myfile_out.open((string(path4)+string("ic_sum.txt")).c_str(),ios::app);
	myfile_out <<  qt_acml.at(i)  << endl;
	myfile_out.close();

	if (i==0) q_imputed.at(i) = qt_acml.at(i);
	if (i>0){
        q_imputed.at(i) = qt_acml.at(i) - qt_acml.at(i-1) ;
	}

	if (i<((int)(t_e_sorted.size()-1))) u_imputed.at(i) = exp(-q_imputed.at(i));
	if (i==((int)(t_e_sorted.size()-1))) u_imputed.at(i) = gsl_ran_flat( r_c, 0.0,exp(-q_imputed.at(i)));


} // END of loop (i) over the infection times


	myfile_out.open((string(path4)+string("residual_gnl.txt")).c_str(),ios::app);
	myfile_out << endl;
	for (int i=0; i<=((int)u_imputed.size()-1);i++){
	if (i!= (int) (u_imputed.size()-1)) myfile_out << u_imputed.at(i) << ",";
	if (i== (int) (u_imputed.size()-1)) myfile_out << u_imputed.at(i);
	}
	myfile_out.close();


	gsl_rng_free(r_c);

}

//------------------

void count_type_all(nt_struct & nt_data_arg, vector<int>& xi_E_current, int& n_base_arg, int& total_count_1, int& total_count_2, int& total_count_3){ // count number of unchanged, transition, transversion (whole dataset)

for (int i=0;i<= (int)(xi_E_current.size()-1);i++){ // loop over all the infected

	int k_E = xi_E_current.at(i);

	switch (nt_data_arg.current_size.at(k_E)>1) {

	case 1:{


	for (int j=0;j<=(nt_data_arg.current_size.at(k_E)-2);j++){


	vector<int> seq_1(nt_data_arg.nt[k_E].begin()+j*(n_base_arg), nt_data_arg.nt[k_E].begin()+(j+1)*(n_base_arg));
	vector<int> seq_2(nt_data_arg.nt[k_E].begin()+(j+1)*(n_base_arg), nt_data_arg.nt[k_E].begin()+(j+2)*(n_base_arg));

	count_type_seq(seq_1, seq_2, n_base_arg, total_count_1, total_count_2, total_count_3); // count between two sequence

	}

	break;
	}

	default:{
	break;
	}

	}
}

}
//-----------------------

void count_type_seq (vector<int>& seq_1_arg, vector<int> seq_2_arg, int& n_base_arg, int& total_count_1, int& total_count_2, int& total_count_3){

int count_1, count_2, count_3; // count_1=count of unchanged sites, ..transition.., transversion
count_1=count_2=count_3=0; 

for ( int i=0;i<=(n_base_arg-1); i++){

	switch(abs(seq_1_arg.at(i)-seq_2_arg.at(i))){
	case 0:{
	count_1 = count_1 + 1;
	break;
	}
	case 1:{
		switch( ((seq_1_arg.at(i)==2) & (seq_2_arg.at(i)==3)) | ((seq_1_arg.at(i)==3) & (seq_2_arg.at(i)==2)) ){
		case 1:{
		count_3 = count_3 + 1;
		break;
		}
		case 0:{
		count_2 = count_2 + 1;
		break;
		}
		}
	break;
	}
	case 2:{
	count_3 = count_3 + 1;
	break;
	}
	case 3:{
	count_3 = count_3 + 1;
	break;
	}
	}

}

total_count_1 = total_count_1 + count_1;
total_count_2 = total_count_2 + count_2;
total_count_3 = total_count_3 + count_3;

}
//------------------------


/*------------------------------------------------*/

void mcmc_UPDATE::source_t_e_update_tri(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vec2int& infecting_list_current_arg, vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, int& subject_proposed, vector<int>& list_update, int iter){ // consider the "triangle" in proposing a new source and respective sequence

//double t_back =10.0;

double acp_pr = 0.0;
double t_low, t_up;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

double log_pr_t_e_forward=0.0; 
double log_pr_t_e_backward=0.0;

double log_pr_seq_forward=0.0; 
double log_pr_seq_backward=0.0;

double log_pr_ds_forward=0.0;
double log_pr_ds_backward=0.0;


double t_proposed;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
vector<double> t_e_modified = t_e_arg;
vector <int> index_modified = index_arg;
vector <int> xi_E_minus_modified = xi_E_minus_arg;

vector <int> xi_U_modified = xi_U_arg;
vector <int> xi_E_modified = xi_E_arg;
vector <int> xi_EnI_modified = xi_EnI_arg;

vector<int> current_size_modified = current_size_arg;

vec2int nt_modified = nt_current_arg;
vec2d t_nt_modified = t_nt_current_arg;

vec2int infecting_list_modified= infecting_list_current_arg;
vector<int> infecting_size_modified = infecting_size_current_arg;

const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter); // set a seed

//vector<int> infected_source_modified = infected_source_current_arg;

vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

//vector<int> nt_subject_seq; // the orginal first sequence of the subject
//nt_subject_seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);

int source_x = infected_source_current_arg.at(subject_proposed);

int rank_source_x; // the rank of removed sequence in old source (note: it has different meaning in function t_e_seq in which the source remains the same)
vector<int> nt_current_source_x;
vector<double> t_nt_current_source_x;
vector<int> nt_modified_source_x;
vector<double> t_nt_modified_source_x;



//-------------- propose a new source --------------//

int source_y;
int rank_source_y;
vector<int> nt_current_source_y;
vector<double> t_nt_current_source_y;
vector<int> nt_modified_source_y;
vector<double> t_nt_modified_source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject

double t_bound = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

//	switch( (t_i_arg.at(xi_I_arg.at(i))<t_e_subject) & (t_r_arg.at(xi_I_arg.at(i))>=t_e_subject)){
	switch( t_i_arg.at(xi_I_arg.at(i))<t_bound){

		case 1:{
			source_pool.push_back(xi_I_arg.at(i));	
		break;
		}
		case 0:{
		break;
		}
	}
}
	

source_pool.insert(source_pool.begin(),9999);

int num_infectious = (int)source_pool.size();

//-----------------------------propose uniformly-------------------------------------------//

//source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_infectious)); // uniformly choose a new source (including bg)

//-----propose according to infectious challenge--------------------//

vector<double> ic(num_infectious);
ic.at(0) = para_current_arg.alpha;

switch(num_infectious>=2){

	case 1:{ // with 2nd sources from pool

		for (int j=1;j<=(num_infectious-1);j++){
		ic.at(j)= para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_pool.at(j)]/norm_const_current_arg.at(source_pool.at(j)); // a new source will be proposed according to the infectious challenges
		}

		double *P=&ic.at(0); // convert vector to array
		gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
		int link= gsl_ran_discrete (r_c, g);
		gsl_ran_discrete_free (g);

		source_y = source_pool.at(link); // a new source
		log_pr_forward = log(ic.at(link));

		switch(source_x==9999){
			case 0:{
				double ic_source_x =  para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_x]/norm_const_current_arg.at(source_x);
				log_pr_backward =log(ic_source_x);
			break;
			}
			case 1:{
				double ic_source_x =  para_current_arg.alpha;
				log_pr_backward =log(ic_source_x);
			break;
			}
		}


		break;
	}

	case 0:{ // only primary source from pool

		source_y = 9999;
		log_pr_forward = log(para_current_arg.alpha);

		double ic_source_x =  para_current_arg.alpha;
		log_pr_backward =log(ic_source_x);

	break;
	}
}


//--------------//

//----------end of proposing a new source -----//

//----------------------------------------------------------------------------------------------------------------//

switch(source_y==source_x){

	case 0:{


		//-------- propose a new t_e---------------//
		switch(source_y){
		
		case 9999:{ // by background
		
			t_up = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(0.0, t_up-t_back);
		
// 			switch( t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 			
// 					switch(t_temp< t_sample_arg.at(subject_proposed)){
// 						case 1:{
// 							t_low =  max(0.0, t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_sample_arg.at(subject_proposed);
// 							t_low =  max(0.0, t_sample_arg.at(subject_proposed) -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no valid t_s
// 					t_low = max(0.0, t_up-t_back);
// 				break;
// 				}
// 			
// 			
// 			}
			
			t_proposed= gsl_ran_flat(r_c,t_low, t_up );
			
			log_pr_t_e_forward = log(1.0/(t_up-t_low));
		
		break;
		}
		
		default :{ // not by background
		
			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_y)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));			
			t_low = max(t_i_arg.at(source_y), t_up -t_back );

// 			double t_temp_2 = min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_y));
// 			switch( t_temp_2!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s (subject) or t_r(source)
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 					switch(t_temp< t_temp_2){
// 						case 1:{
// 							t_low =  max(t_i_arg.at(source_y), t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_temp_2;
// 							t_low =  max(t_i_arg.at(source_y), t_temp_2 -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no with valid t_s (subject) and t_r(source)
// 					t_low =  max(t_i_arg.at(source_y), t_up - t_back);
// 				break;
// 				}
// 			}

				
			t_proposed= gsl_ran_flat(r_c, t_low, t_up );

			log_pr_t_e_forward = log(1.0/(t_up-t_low));

		break;
		}
		
		}

		//------------//
		t_e_modified.at(subject_proposed) = t_proposed;
		//------------//

		switch(source_x){
		
		case 9999:{ // by background
		
			t_up = min( t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(0.0, t_up-t_back);

// 			switch( t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 			
// 					switch(t_temp< t_sample_arg.at(subject_proposed)){
// 						case 1:{
// 							t_low =  max(0.0, t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_sample_arg.at(subject_proposed);
// 							t_low =  max(0.0, t_sample_arg.at(subject_proposed) -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no valid t_s
// 					t_low = max(0.0, t_up-t_back);
// 				break;
// 				}
// 			
// 			
// 			}

						
			log_pr_t_e_backward = log(1.0/(t_up-t_low));
		
		break;
		}
		
		default :{ // not by background
		
			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_x)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
			t_low = max(t_i_arg.at(source_x), t_up -t_back );

// 			double t_temp_2 = min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_x));	
// 			switch( t_temp_2!=unassigned_time_CUPDATE){
// 			
// 				case 1:{// with valid t_s (subject) or t_r(source)
// 					double t_temp = min( t_i_arg.at(subject_proposed), t_max_CUPDATE) - t_back;
// 					switch(t_temp< t_temp_2){
// 						case 1:{
// 							t_low =  max(t_i_arg.at(source_x), t_temp);
// 						break;
// 						}
// 						case 0:{// should be unlikely if t_back is large enough
// 							double dt = t_temp -  t_temp_2;
// 							t_low =  max(t_i_arg.at(source_x), t_temp_2 -  dt);
// 						break;
// 						}
// 					}
// 				break;
// 				}
// 			
// 				case 0:{ // no with valid t_s (subject) and t_r(source)
// 					t_low =  max(t_i_arg.at(source_x), t_up - t_back);
// 				break;
// 				}
// 			}


			log_pr_t_e_backward = log(1.0/(t_up-t_low));

		break;
		}
		
		}

		//--------------------------------------//

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00test.txt")).c_str(),ios::app);
// myfile_mcmc_out << 0 << endl;
// myfile_mcmc_out.close();		
		
		vector<int> seq_proposed(n_base_CUPDATE); // newly proposed sequence when source changes
		vector<int> nt_past_forward(n_base_CUPDATE);
		vector<int> nt_future_forward(n_base_CUPDATE);
		double t_past, t_future;

		vector<int> nt_tri_1_forward(n_base_CUPDATE);
		vector<int> nt_tri_2_forward(n_base_CUPDATE);
		double t_tri_1_forward, t_tri_2_forward;

		vector<int> nt_tri_1_backward(n_base_CUPDATE);
		vector<int> nt_tri_2_backward(n_base_CUPDATE);
		double t_tri_1_backward, t_tri_2_backward;

		vector<int> seq_proposed_backward(n_base_CUPDATE);
		vector<int> nt_past_backward(n_base_CUPDATE);
		vector<int> nt_future_backward(n_base_CUPDATE);
		double t_proposed_backward, t_past_backward, t_future_backward;
		
		//------------------------------------------------
		
		switch(source_y==9999){
		
			case 0:{// new 2nd infection

				infecting_size_modified.at(source_y) = infecting_size_modified.at(source_y) +1;
				
				vector<double> t_y(infecting_size_current_arg.at(source_y));
				for (int i=0;i<=(infecting_size_current_arg.at(source_y)-1);i++){
				t_y.at(i) = t_e_arg.at(infecting_list_current_arg[source_y][i]);
				}
				t_y.push_back(t_proposed);
				sort(t_y.begin(), t_y.end());
			
				int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_proposed));
				infecting_list_modified.at(source_y).insert(infecting_list_modified.at(source_y).begin()+rank_y, subject_proposed);
				
				//----------------------------------------------------//
	
				nt_current_source_y = nt_current_arg.at(source_y);
				t_nt_current_source_y = t_nt_current_arg.at(source_y);
		
				nt_modified_source_y = nt_current_source_y;
				t_nt_modified_source_y = t_nt_current_source_y;
		
				t_nt_modified_source_y.push_back(t_proposed);
		
				sort( t_nt_modified_source_y.begin(),  t_nt_modified_source_y.end()); 
		
				rank_source_y = distance(t_nt_modified_source_y.begin(), find(t_nt_modified_source_y.begin(),t_nt_modified_source_y.end(), t_proposed) );
		
				current_size_modified.at(source_y) = current_size_modified.at(source_y) + 1;
		
				//----------------------------------------------------------------------------------------------------------------//
				t_past = t_nt_modified_source_y.at(rank_source_y-1);
				nt_past_forward.assign(nt_current_source_y.begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE);
		
				//t_proposed = t_e_subject;
		
				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y (take closer seq as the future seq)

								nt_tri_1_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
								t_tri_1_forward = t_nt_modified_source_y.at(rank_source_y+1);

								nt_tri_2_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
								t_tri_2_forward =  t_nt_current_arg[subject_proposed][1];
						
								seq_propose_tri (seq_proposed,  log_pr_seq_forward, nt_past_forward, nt_tri_1_forward, nt_tri_2_forward, t_past, t_tri_1_forward, t_tri_2_forward, t_proposed,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);


// 								switch(t_nt_current_arg[subject_proposed][1]<t_nt_modified_source_y.at(rank_source_y +1)){
// 		
// 									case 1:{// take the one from subject as future seq
// 
// 										t_future  = t_nt_current_arg[subject_proposed][1]; 
// 										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
// 																				
// 									seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
// 		
// 									break;
// 									}
// 		
// 									case 0:{ // take the one from source_y as future seq
// 										t_future = t_nt_modified_source_y.at(rank_source_y+1);
// 										nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
// 		
// 										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
// 		
// 									break;
// 									}
// 								}
// 		
							break;
							}
		
							case 0:{//  inserted seq will be last seq at rank_source_y  in source_y
		
										t_future  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
	
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y
								t_future = t_nt_modified_source_y.at(rank_source_y+1);
								nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
								seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);						
							break;
							}
		
							case 0:{// inserted seq will be last seq at rank_source_y  in source_y
								t_future = t_proposed;
								seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				//nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_subject_seq.begin(), nt_subject_seq.end());
				nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
		
				nt_modified.at(source_y) = nt_modified_source_y;
				t_nt_modified.at(source_y) = t_nt_modified_source_y;

			break;
			}
		
			case 1:{// new bg infection

				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject

						t_past  = t_nt_current_arg[subject_proposed][1]; // yes, t_past!
						nt_past_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
		
						t_future = t_proposed;

						seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);			
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
		
						log_pr_seq_forward = n_base_CUPDATE*log(0.25);
		
						for (int i=0; i<=(n_base_CUPDATE-1);i++){
						seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
						}
		
					break;
					}
				}
			break;
			}
		
		}
		
		//--------------------------------------------
		nt_modified_subject.erase(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE); 
		nt_modified_subject.insert(nt_modified_subject.begin(), seq_proposed.begin(), seq_proposed.end());

		t_nt_modified_subject.erase(t_nt_modified_subject.begin());
		t_nt_modified_subject.push_back(t_proposed); 
		sort( t_nt_modified_subject.begin(),  t_nt_modified_subject.end()); 

		nt_modified.at(subject_proposed) = nt_modified_subject;
		t_nt_modified.at(subject_proposed) = t_nt_modified_subject;
		//---------------------------------------------
		
		switch(source_x==9999){
		
			case 0:{// was 2nd infection

				infecting_size_modified.at(source_x) = infecting_size_modified.at(source_x) - 1;

				int rank_x = distance(infecting_list_current_arg.at(source_x).begin(), find(infecting_list_current_arg.at(source_x).begin(), infecting_list_current_arg.at(source_x).end(), subject_proposed));
			
				infecting_list_modified.at(source_x).erase(infecting_list_modified.at(source_x).begin()+rank_x);

				//----------------------------------------------------------------------------------------------------------------//

				nt_current_source_x = nt_current_arg.at(source_x);
				t_nt_current_source_x = t_nt_current_arg.at(source_x);
		
				rank_source_x = distance( t_nt_current_source_x.begin(), find(t_nt_current_source_x.begin(), t_nt_current_source_x.end(), t_e_arg.at(subject_proposed)) ); 
		
				//----------------------------------------------------------------------------------------------------------------//
		
				t_past_backward = t_nt_current_source_x.at(rank_source_x-1);
				nt_past_backward.assign(nt_current_source_x.begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x)*n_base_CUPDATE);
		
				t_proposed_backward =t_e_arg.at(subject_proposed);
				seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
				
					case 1:{// with2nd seq subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
		
							case 1:{// also NOT last seq at rank_source_x  in source_x (take closer seq as the future seq)

								nt_tri_1_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
								t_tri_1_backward = t_nt_current_source_x.at(rank_source_x +1);

								nt_tri_2_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
								t_tri_2_backward =  t_nt_current_arg[subject_proposed][1]; 
						
								seq_backward_pr_tri (seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, nt_tri_1_backward, nt_tri_2_backward, t_past_backward, t_tri_1_backward, t_tri_2_backward, t_proposed_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

// 								switch(t_nt_current_arg[subject_proposed][1]<t_nt_current_source_x.at(rank_source_x +1)){
// 								//switch(t_sample_arg.at(subject_proposed)<t_nt_current_source_x.at(rank_source_x +1)){
// 
// 									case 1:{// take the one from subject as future seq
// 
// 										t_future_backward  = t_nt_current_arg[subject_proposed][1]; 
// 										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
// 						
// 										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
// 
// 									break;
// 									}
// 									case 0:{// take the one from source_x as future seq
// 										t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
// 										nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
// 						
// 										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
// 									break;
// 									}
// 								}

							break;
							}
		
							case 0:{ //last seq at rank_source_x  in source_x

								t_future_backward  = t_nt_current_arg[subject_proposed][1]; // always takes the 2nd seq in subject as the future sequence (both for forward and backward)
								nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}	
						}
		
		
					break;
					}
					
					case 0:{ // with no 2nd seq in the subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
							case 1:{ // not last seq at rank_source_x  in source_x
								t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
								nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
								
							break;
							}
							case 0:{ // last seq rank_source_x  in source_x
								t_future_backward = t_proposed_backward;
		
								seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				nt_modified_source_x = nt_current_source_x;
				t_nt_modified_source_x = t_nt_current_source_x;
				
				t_nt_modified_source_x.erase(t_nt_modified_source_x.begin() + rank_source_x); // erase the original t_nt entry for source_x
				nt_modified_source_x.erase(nt_modified_source_x.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source_x.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source_x

				nt_modified.at(source_x) = nt_modified_source_x;
				t_nt_modified.at(source_x) = t_nt_modified_source_x;


				current_size_modified.at(source_x) = current_size_modified.at(source_x) - 1;
		
			break;
			}
		
			case 1:{// was bg infection
		
				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){

					case 1:{// with 2nd seq in subject

						t_past_backward = t_nt_current_arg[subject_proposed][1]; // note: it is "past"!
						nt_past_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
		
						t_proposed_backward = t_e_arg.at(subject_proposed);
						seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
						t_future_backward = t_proposed_backward;

						seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
						log_pr_seq_backward = n_base_CUPDATE*log(0.25); // note: this is the proposal density! not refers to the model!
					break;
					}
				}
		
			break;
			}
		
		}
		
		
		//------------- deal with change of likelihood due to change of source and sequences in subject_proposed)---------------//
		
		switch (current_size_arg.at(subject_proposed)>1) {
		
		case 1:{
		
		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below
		
		lh_square_modified.log_f_S.at(subject_proposed) = 0.0;
		
		for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){

			vector<int> seq_1(nt_modified.at(subject_proposed).begin()+j*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE));
			vector<int> seq_2(nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+2)*(n_base_CUPDATE));
			
			lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2,t_nt_modified_subject.at(j), t_nt_modified_subject.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);	
		}
		
		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 
		
		break;
		}
		
		default:{
		break;
		}
		
		}


		//-----------------------------------------------------------------------------------------//

		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
		
		for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again
		
			if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
			
				switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
					case 1:{
						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
					
					break;
					}
					case 0:{
						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
					break;
					}
				}
			lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));
	
			}
		}
		//----------
		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;

		switch(source_y){
		
		case 9999:{ // by background
		lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
		break;
		}
		
		default :{ // not by background
		lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y);
		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		break;
		}
		
		}
		//----------

		switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes
		
			case 1:{
			
				switch(t_proposed<t_e_arg.at(index_arg.at(0))){ // original indexes would be replace by the chosen subject
				
					case 1:{
				
						index_modified.clear();	
						index_modified.assign(1,subject_proposed);// replace index 
						xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed));
					
						for ( int i =0; i<= (int) (index_arg.size()-1); i++){
						xi_E_minus_modified.push_back(index_arg.at(i)); // the original indexes in xi_E_minus now
						}	
					
						log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
						
						lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
						lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
						lh_square_modified.q_E.at(subject_proposed)=0.0;
						lh_square_modified.g_E.at(subject_proposed)=1.0;
						lh_square_modified.h_E.at(subject_proposed)=1.0;
						lh_square_modified.f_E.at(subject_proposed)=1.0;
						
					
						for (int i=0; i<=(int) (index_arg.size()-1);i++){ // the original indexes have to be acocunted in likelihood now
						
						//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_arg.at(i)));
						
						lh_square_modified.g_E.at(index_arg.at(i)) = para_current_arg.alpha; // this is not the subject_proposed
						lh_square_modified.q_E.at(index_arg.at(i)) =para_current_arg.alpha*t_e_arg.at(index_arg.at(i));
						lh_square_modified.h_E.at(index_arg.at(i)) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(index_arg.at(i)),1.0);
						lh_square_modified.f_E.at(index_arg.at(i)) = lh_square_modified.g_E.at(index_arg.at(i))*lh_square_modified.h_E.at(index_arg.at(i));
						
						log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(index_arg.at(i)));
						}
						
					break;
					}
					
				
					case 0:{
					
						if (t_proposed==t_e_arg.at(index_arg.at(0))){ // addtion of one more index
						
							log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would have to be removed from likelihood function
							index_modified.push_back(subject_proposed); // add index 
							xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),subject_proposed)); // removed from xi_E_minus
							lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
							lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
							lh_square_modified.q_E.at(subject_proposed)=0.0;
							lh_square_modified.g_E.at(subject_proposed)=1.0;
							lh_square_modified.h_E.at(subject_proposed)=1.0;
							lh_square_modified.f_E.at(subject_proposed)=1.0;
						
						
						}
						
						if (t_proposed>t_e_arg.at(index_arg.at(0))){ // no shift of cases between xi_E and xi_E_minus

							log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below

							lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
					// 		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
							lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
							lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
						
							log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
					
					
						} // end if t_proposs>t_e_arg.at()
						
					break;
					}
				
				}
			
			break;
			}
		
		
			case 0: { // when chosen subject is one of the indexes

				index_modified.clear();
			
				int first_min = distance(t_e_modified.begin(), min_element(t_e_modified.begin(), t_e_modified.end()));
				double min_t = t_e_modified.at(first_min); // the minimum time of exposure
				
				int num_min = (int) count(t_e_modified.begin(), t_e_modified.end(), min_t); // numberof subects with the min exposure time
			
				switch (num_min>1) {
					case 1: {
						index_modified.reserve(n_CUPDATE);	
						for (int i=0; i<=(n_CUPDATE-1);i++){
						if (t_e_modified.at(i)==min_t ) index_modified.push_back(i);		
						}
					break;
					}
					case 0:{
						index_modified.assign(1,first_min);
					break;
					}
				}
			
				xi_E_minus_modified = xi_E_arg;
			
			
				for (int i=0;i<= (int) (index_modified.size()-1); i++){
			
				xi_E_minus_modified.erase(find(xi_E_minus_modified.begin(),xi_E_minus_modified.end(),index_modified.at(i)));
			
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(index_modified.at(i))); // this subject would have to be removed from likelihood function ( new index might be orginally an index, but the the log(lh_square_modified.f_E.at(index_modified.at(i)) will be zero in this case)
			
				lh_square_modified.k_sum_E.at(index_modified.at(i))=0.0;
				lh_square_modified.kt_sum_E.at(index_modified.at(i))=0.0;
				lh_square_modified.q_E.at(index_modified.at(i))=0.0;
				lh_square_modified.g_E.at(index_modified.at(i))=1.0;
				lh_square_modified.h_E.at(index_modified.at(i))=1.0;
				lh_square_modified.f_E.at(index_modified.at(i))=1.0;			
				}
			
				switch(find(index_modified.begin(),index_modified.end(),subject_proposed) ==index_modified.end() ){ //return 1 when the chosen  subject is NO longer an index
					case 1:{
				
						log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); 	
					
						lh_square_modified.q_E.at(subject_proposed) = para_current_arg.alpha*t_proposed + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.kt_sum_E.at(subject_proposed);
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.alpha + para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						lh_square_modified.h_E.at(subject_proposed) = gsl_ran_exponential_pdf(lh_square_modified.q_E.at(subject_proposed),1.0);
						lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
					
						log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); //add back the  part of likelihood 
						
					break;
					}
					case 0:{
					break;
					}
				}
			
			break;
			}
		
		}
		
		//--------------------//


		switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
			case 1:{
		
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
				lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
				log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 
				
			break;
			}

			case 0:{
			break;
			}
		}
		
		//----------
		
		switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
			case 1:{
			
				log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
				lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
				log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 
				
			break;
			}
			case 0:{
			break;
			}
		}


		//-----------------------------------------------------------------------------------------//

		switch(source_x){
		
			case 9999:{ // was background infection
		
				switch(source_y){
					case 9999:{ // new  bg infection 		
					break;
					}
		
					default:{ // new secondary infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) = 0.0;
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); // in fact redudancy, but here for clarity
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
			break;
			}
			
			default :{ // was secondary infection
		
				switch(source_y){
					case 9999:{ // new  bg infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); // redudant as second term must be zero
		
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) =  n_base_CUPDATE*log(0.25);
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
			
		
				
					break;
					}
		
					default:{ // new secondary infection
			
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
		
						//--
		
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
		
						//---
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
		
			break;
			}
		
		}
		
		
		
	//------------- end of with change of likelihood (due to change of source and sequences in subject_proposed)---------------//
	
		//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward)*exp(log_pr_seq_backward-log_pr_seq_forward)*exp(log_pr_ds_backward-log_pr_ds_forward));



		acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+(log_pr_seq_backward-log_pr_seq_forward) +(log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_ds_backward-log_pr_ds_forward)) );
		
		
		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		
		lh_square_current_arg = lh_square_modified;
		log_lh_current_arg = log_lh_modified;
		
		//nt_current_arg.at(subject_proposed) = nt_modified_subject;
		nt_current_arg = nt_modified;
		t_nt_current_arg = t_nt_modified;

		delta_mat_current_arg = delta_mat_modified;
		t_e_arg= t_e_modified;
		index_arg = index_modified;
		xi_E_minus_arg = xi_E_minus_modified;

		current_size_arg = current_size_modified;
		infected_source_current_arg.at(subject_proposed) =  source_y;
		//infected_source_current_arg = infected_source_modified;

		infecting_list_current_arg = infecting_list_modified;
		infecting_size_current_arg = infecting_size_modified;
		
// 			switch (source_x){
// 				
// 				case 9999:{ 
// 				break;
// 				}
// 				
// 				default :{ 
// 				//nt_current_arg.at(source_x) = nt_modified_source_x;
// 				t_nt_current_arg.at(source_x) = t_nt_modified_source_x;	
// 				break;	
// 				}
// 			}
// 		
// 			switch (source_y){
// 				
// 				case 9999:{ 
// 				break;
// 				}
// 				
// 				default :{ 
// 				//nt_current_arg.at(source_y) = nt_modified_source_y;
// 				t_nt_current_arg.at(source_y) = t_nt_modified_source_y;
// 				break;	
// 				}
// 			}

	break;
	}

	case 0: {
	break;
	}
}

break;
}

case 1:{ // source_y==source_x
break;
}
}

gsl_rng_free(r_c);

}

//--------------------------------------------------


inline void seq_propose_tri (vector<int>& seq_proposed,  double& log_pr_seq_forward, const vector<int>& nt_past_forward, const vector<int>& nt_tri_1_forward, const vector<int>& nt_tri_2_forward, const double& t_past, const double& t_tri_1_forward, const double& t_tri_2_forward, const double& t_proposed,  const double& mu_1, const double& mu_2, const  int& n_base_CUPDATE, gsl_rng * r_c){

double lambda=0.5;
double total_pr = 0.5;

double dt_1 = fabs(t_proposed - t_past);
double dt_2 = fabs(t_proposed - t_tri_1_forward);
double dt_3 = fabs(t_proposed - t_tri_2_forward);

double total_risk = exp(-lambda*dt_1) + exp(-lambda*dt_2) + exp(-lambda*dt_3);

double P[4] = {total_pr*exp(-lambda*dt_1)/total_risk, total_pr*exp(-lambda*dt_2)/total_risk, total_pr*exp(-lambda*dt_3)/total_risk, 1.0-total_pr};

gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (sizeof(P)/sizeof(P[0]),P);

for (int i=0; i<=(n_base_CUPDATE-1); i++){

	int type = gsl_ran_discrete (r_c, g);

// 
// 		myfile_out.open((string(path4)+string("11.txt")).c_str(),ios::app);
// 		myfile_out <<log_pr_seq_forward << endl;
// 		myfile_out.close();	

	switch(type){
		case 0:{
			seq_proposed.at(i) = nt_past_forward.at(i);
			log_pr_seq_forward =  log_pr_seq_forward + log(P[0]);
		break;
		}
		case 1:{
			seq_proposed.at(i) = nt_tri_1_forward.at(i);
			log_pr_seq_forward =  log_pr_seq_forward + log(P[1]);
		break;
		}
		case 2:{
			seq_proposed.at(i) = nt_tri_2_forward.at(i);
			log_pr_seq_forward =  log_pr_seq_forward + log(P[2]);
		break;
		}
		case 3:{
			seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
			log_pr_seq_forward =  log_pr_seq_forward + log(P[3]) + log(0.25);
		break;
		}
	}

// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// 		myfile_out <<type << "," << nt_past_forward.at(i) << "," <<nt_tri_1_forward.at(i) <<"," << nt_tri_2_forward.at(i) <<","<< seq_proposed.at(i) << endl;
// 		myfile_out.close();	

}


gsl_ran_discrete_free (g);

}

//--------------------------------------------------


inline void seq_backward_pr_tri (const vector<int>& seq_proposed_backward,  double& log_pr_seq_backward, const vector<int>& nt_past_backward, const vector<int>& nt_tri_1_backward, const vector<int>& nt_tri_2_backward, const double& t_past_backward, const double& t_tri_1_backward, const double& t_tri_2_backward, const double& t_proposed_backward, const double& mu_1, const double& mu_2, const int& n_base_CUPDATE){


double lambda=0.5;
double total_pr = 0.5;

double dt_1 = fabs(t_proposed_backward - t_past_backward);
double dt_2 = fabs(t_proposed_backward - t_tri_1_backward);
double dt_3 = fabs(t_proposed_backward - t_tri_2_backward);

double total_risk = exp(-lambda*dt_1) + exp(-lambda*dt_2) + exp(-lambda*dt_3);

double P[4] = {total_pr*exp(-lambda*dt_1)/total_risk, total_pr*exp(-lambda*dt_2)/total_risk, total_pr*exp(-lambda*dt_3)/total_risk, 1.0-total_pr};

//gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (sizeof(P)/sizeof(P[0]),P);


for (int i=0; i<=(n_base_CUPDATE-1); i++){
	
	double pr_sum=0.0; // an "or" probability

	switch( seq_proposed_backward.at(i)==nt_past_backward.at(i)){
		case 1:{
			pr_sum  = pr_sum + P[0];
		break;
		}
		case 0:{
		break;
		}
	}

	switch( seq_proposed_backward.at(i)==nt_tri_1_backward.at(i)){
		case 1:{
			pr_sum  = pr_sum + P[1];
		break;
		}
		case 0:{
		break;
		}
	}

	switch( seq_proposed_backward.at(i)==nt_tri_2_backward.at(i)){
		case 1:{
			pr_sum  = pr_sum + P[2];
		break;
		}
		case 0:{
		break;
		}
	}

	log_pr_seq_backward = log_pr_seq_backward + log(pr_sum+P[3]*0.25);

// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("11.txt")).c_str(),ios::app);
// 		myfile_out <<pr_sum << endl;
// 		myfile_out.close();

}

// 		ofstream myfile_out; 
// 		myfile_out.open((string(path4)+string("33.txt")).c_str(),ios::app);
// 		myfile_out <<log_pr_seq_backward << endl;
// 		myfile_out.close();	

}

//--------------------------------------------------


void mcmc_UPDATE::seq_n_update(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vector<int>&  xi_beta_E_arg, const int& subject_proposed,  gsl_rng* & r_c){

//int subject_proposed ;
double acp_pr = 0.0;

double log_part_x_subject=0.0;
double log_part_y_subject=0.0;

double log_part_x_source=0.0;
double log_part_y_source =0.0;


//int position_proposed, base_proposed;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;


//position_proposed =iter;

int subject_source = infected_source_current_arg.at(subject_proposed);


//vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);

vector<int> seq_proposed (n_base_CUPDATE); 
//seq_proposed.assign(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE );

vector<int> nt_current (n_base_CUPDATE);
nt_current.assign( nt_current_arg.at(subject_proposed).begin() , nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);

vector<int> nt_next_subject (n_base_CUPDATE);

//---

for (int i=0;i<=(n_base_CUPDATE-1);i++){

	int base_current = nt_current_arg[subject_proposed][i]; // always refers to the first sequence
	
	switch(base_current){
		case 1:{
			int type = gsl_rng_uniform_int (r_c, 3);
			switch(type){
				case 0:{
					seq_proposed.at(i)  = 2;
				break;
				}
				case 1:{
					seq_proposed.at(i) = 3;
				break;
				}
				case 2:{
					seq_proposed.at(i) = 4;
				break;
				}
			}
		break;
		}
		case 2:{
			int type = gsl_rng_uniform_int (r_c, 3);
	
			switch(type){
				case 0:{
					seq_proposed.at(i) = 1;
				break;
				}
				case 1:{
					seq_proposed.at(i) = 3;
				break;
				}
				case 2:{
					seq_proposed.at(i) = 4;
				break;
				}
			}	
		break;
		}
		case 3:{
			int type = gsl_rng_uniform_int (r_c, 3);
			switch(type){
				case 0:{
					seq_proposed.at(i) = 1;
				break;
				}
				case 1:{
					seq_proposed.at(i) = 2;
				break;
				}
				case 2:{
					seq_proposed.at(i) = 4;
				break;
				}
			}	
		break;
		}
		case 4:{
			int type = gsl_rng_uniform_int (r_c, 3);
			switch(type){
				case 0:{
					seq_proposed.at(i) = 1;
				break;
				}
				case 1:{
					seq_proposed.at(i) = 2;
				break;
				}
				case 2:{
					seq_proposed.at(i) = 3;
				break;
				}
			}	
		break;
		}
	}

}
//---

//nt base_next_subject =0;

switch (current_size_arg.at(subject_proposed)>1) {

	case 1:{
		//--
		
//		base_next_subject = nt_current_arg[subject_proposed][position_proposed + n_base_CUPDATE];

// 		log_part_x_subject = log_lh_base(base_current, base_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2);
// 
// 		log_part_y_subject = log_lh_base(base_proposed, base_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2);

		nt_next_subject.assign( nt_current_arg.at(subject_proposed).begin() + n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+2*n_base_CUPDATE);

		log_part_x_subject = log_lh_seq (nt_current, nt_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

		log_part_y_subject = log_lh_seq (seq_proposed, nt_next_subject, t_nt_current_arg[subject_proposed][0], t_nt_current_arg[subject_proposed][1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

		//--
		lh_square_modified.log_f_S.at(subject_proposed) =  lh_square_modified.log_f_S.at(subject_proposed) - log_part_x_subject;
		log_lh_modified = log_lh_modified - log_part_x_subject;
		
		lh_square_modified.log_f_S.at(subject_proposed) =  lh_square_modified.log_f_S.at(subject_proposed) + log_part_y_subject;
		log_lh_modified = log_lh_modified + log_part_y_subject;


	break;
	}

	case 0:{
	break;
	}
}
//-------

int rank_source =-1;  //count the rank of the original t_e among t_nt_current_arg.at(subject_source) 

// int base_before_source =0;
// int base_next_source = 0;
vector<int> nt_before_source(n_base_CUPDATE);
vector<int> nt_next_source(n_base_CUPDATE);


switch(subject_source ){

	case 9999:{ // by background
	break;
	}
	
	default :{ // not by background
	
	//nt_modified_source = nt_current_arg.at(subject_source);
	rank_source = distance( t_nt_current_arg.at(subject_source).begin(), find(t_nt_current_arg.at(subject_source).begin(), t_nt_current_arg.at(subject_source).end(), t_e_arg.at(subject_proposed)) );
	//nt_modified_source.erase(nt_modified_source.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source
	
	//base_before_source =  nt_current_arg[subject_source][(rank_source-1)*n_base_CUPDATE + position_proposed];
	nt_before_source.assign(nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*(rank_source-1), nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*rank_source );

		switch(current_size_arg.at(subject_source)>(rank_source+1)){
			case 1:{// there  is a valid base_next_source

// 				base_next_source =  nt_current_arg[subject_source][(rank_source+1)*n_base_CUPDATE + position_proposed];
	
// 				log_part_x_source = log_lh_base(base_before_source, base_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2) + log_lh_base(base_current, base_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2);
// 	
// 				log_part_y_source = log_lh_base(base_before_source, base_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2) + log_lh_base(base_proposed, base_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2);
// 	

				nt_next_source.assign(nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*(rank_source+1), nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*(rank_source+2) );	
	
				log_part_x_source = log_lh_seq (nt_before_source, nt_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE) + log_lh_seq (nt_current, nt_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
				log_part_y_source = log_lh_seq (nt_before_source, seq_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE) +log_lh_seq (seq_proposed, nt_next_source, t_nt_current_arg[subject_source][rank_source], t_nt_current_arg[subject_source][rank_source+1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);				
	
			break;
			}
			
			case 0:{
	
// 				log_part_x_source = log_lh_base(base_before_source, base_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2);
// 	
// 				log_part_y_source = log_lh_base(base_before_source, base_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2);

				log_part_x_source = log_lh_seq (nt_before_source, nt_current, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE) ;
		
				log_part_y_source = log_lh_seq (nt_before_source, seq_proposed, t_nt_current_arg[subject_source][rank_source-1], t_nt_current_arg[subject_source][rank_source], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);		
		
				break;
			}
		}
	
	lh_square_modified.log_f_S.at(subject_source) =  lh_square_modified.log_f_S.at(subject_source) - log_part_x_source;
	log_lh_modified = log_lh_modified - log_part_x_source;
	
	
	lh_square_modified.log_f_S.at(subject_source) =  lh_square_modified.log_f_S.at(subject_source) + log_part_y_source;
	log_lh_modified = log_lh_modified + log_part_y_source;
	
	break;
	}

}

//------

double log_part_y = log_part_y_subject + log_part_y_source;
double log_part_x = log_part_x_subject + log_part_x_source;

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg));
acp_pr = min(1.0,exp(log_part_y-log_part_x));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();
// 
// myfile_mcmc_out.open((string(path4)+string("11.txt")).c_str(),ios::app);
// myfile_mcmc_out << log_part_x_subject <<"," << log_part_x_source << "," << log_part_y_subject <<","<< log_part_y_source<< endl;
// myfile_mcmc_out.close();


switch(uniform_rv<=acp_pr){
case 1: {

lh_square_current_arg = lh_square_modified;
log_lh_current_arg = log_lh_modified;

//nt_current_arg[subject_proposed][position_proposed]= base_proposed;

nt_current_arg.at(subject_proposed).erase(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);  //erase the original nt entry for source

nt_current_arg.at(subject_proposed).insert(nt_current_arg.at(subject_proposed).begin(), seq_proposed.begin(), seq_proposed.end());  //insert  the new nt 

	switch (subject_source){
		
		case 9999:{ // by background
		break;
		}
		
		default :{ // not by background
//			nt_current_arg[subject_source][(rank_source)*n_base_CUPDATE + position_proposed] = base_proposed;

			nt_current_arg.at(subject_source).erase(nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*rank_source, nt_current_arg.at(subject_source).begin()+n_base_CUPDATE*(rank_source+1) );  //erase the original nt entry for source
			
			nt_current_arg.at(subject_source).insert(nt_current_arg.at(subject_source).begin()+(rank_source)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());  //insert  the new nt 
		break;
		}
	}
			

break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);


}


/*------------------------------------------------*/

void delete_seq_samples(para_key& para_current, para_aux& para_other,  vector<int>& xi_I_current, vector<int>& xi_U_current, vector<int>& xi_E_current, vector<int>& xi_E_minus_current, vector<int>& xi_R_current,  vector<int>& xi_EnI_current,  vector<int>& xi_EnIS_current,  vector<int>& xi_InR_current,  vector<double>& t_e_current, vector<double>& t_i_current, vector<double>& t_r_current, vector<int>& index_current, vector<double>& stb_current, vector<int>& gp_stb_current, vector<int>& infected_source_current, vector < vector<double> >& kernel_mat_current, vector <double>& norm_const_current, vec2int& sample_data, vector<double>& t_onset, nt_struct& nt_data_current){ // keep the true config, only delete the sampled sequences


const gsl_rng_type* T= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r = gsl_rng_alloc (T); // r is pointer points to an object with Type T

int seed_intial_mcmc = 1;  //1,999,-1000,-10000,123456

gsl_rng_set (r, seed_intial_mcmc); // set a seed

ofstream myfile_out; 




para_current.alpha = para_current.alpha*1; //initialization of parameter to be estimated
para_current.beta = para_current.beta *1 ; //initialization of parameter to be estimated

// para_current.mu_lat = para_current.a*para_current.b*1; //initialization of parameter to be estimated
// para_current.var_lat = para_current.a*para_current.b*para_current.b*1; //initialization of parameter to be estimated
para_current.mu_lat = para_current.a; //initialization of parameter to be estimated
para_current.var_lat = para_current.b; //initialization of parameter to be estimated

para_current.c = para_current.c*1; //initialization of parameter to be estimated
para_current.d = para_current.d*1; //initialization of parameter to be estimated

para_current.k_1 = para_current.k_1*1; //initialization of parameter to be estimated
para_current.k_2 =  para_current.k_2*1; //initialization of parameter to be estimated

para_current.mu_1 = 1.0*para_current.mu_1; //initialization of parameter to be estimated
para_current.mu_2 = 1.0*para_current.mu_2; //initialization of parameter to be estimated

para_current.p_ber = 1.0*para_current.p_ber; //initialization of parameter to be estimated

para_current.stb_1 = 1.0; //initialization of parameter to be estimated
para_current.stb_2 = 1.0; //initialization of parameter to be estimated

//-------- for partial deletion----///

double p_sample = 0.12;   // pr a sample included for an exposure

//--- this part is to match the usage of seed and xi_E_current.size() in initializemcmc such that the unsampled individuals match
xi_E_current =xi_I_current; // individuals gone through I (assumed known here) would be initialized as infected; sampled would also be infected (see below); assume index has gone through I

xi_EnI_current.clear() ;

for (int i=0; i<= (int)(para_other.n-1);i++){

	switch((nt_data_current.t_sample.at(i)!=para_other.unassigned_time)&(find(xi_I_current.begin(),xi_I_current.end(),i)==xi_I_current.end())){

		case 1:{ // sampled(infected) but not in xi_I_current
			xi_E_current.push_back(i);
			xi_EnI_current.push_back(i);
		break;
		}
		case 0:{
		break;
		}

	}

}

xi_EnIS_current.clear() ; // this would be empty as initial value

xi_U_current.clear();

for (int i=0; i<= (int)(para_other.n-1);i++){
	if(find(xi_E_current.begin(),xi_E_current.end(),i)==xi_E_current.end()){
	xi_U_current.push_back(i);
	t_e_current.at(i)=para_other.unassigned_time;
 	infected_source_current.at(i) = -99;
	}
}
for (int i=0; i<= (int)(xi_I_current.size()-1);i++){
	gsl_ran_flat(r, 0,1);
}

//---
	
for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections

	int subject= xi_E_current.at(i);

	double PS[2] = {1.0 - p_sample, p_sample};
	gsl_ran_discrete_t * gs = gsl_ran_discrete_preproc (sizeof(PS)/sizeof(PS[0]),PS);
	int sample_ind = gsl_ran_discrete (r, gs); // 1 = would include the sample

	
	gsl_ran_discrete_free(gs);

	switch(sample_ind){
		case 0:{// exclude
			
			if ((find( index_current.begin(), index_current.end(), xi_E_current.at(i) )==index_current.end())&(nt_data_current.t_sample.at(xi_E_current.at(i)) !=para_other.unassigned_time)){ // if it is an non-index & it has sample
		
				if((find(xi_EnI_current.begin(),xi_EnI_current.end(),xi_E_current.at(i))!=xi_EnI_current.end())&(nt_data_current.t_sample.at(xi_E_current.at(i)) !=para_other.unassigned_time)){// this infection is in xi_EnI& had a sample(which the sample gotta be deleted)
				xi_EnIS_current.push_back(xi_E_current.at(i));
				}

				int rank = distance( nt_data_current.t_nt.at(subject).begin(), find(nt_data_current.t_nt.at(subject).begin(),nt_data_current.t_nt.at(subject).end(),nt_data_current.t_sample.at(subject) ) ); 

				nt_data_current.t_nt.at(subject).erase(nt_data_current.t_nt.at(subject).begin() + rank);
				nt_data_current.nt.at(subject).erase(nt_data_current.nt.at(subject).begin()+para_other.n_base*rank , nt_data_current.nt.at(subject).begin()+para_other.n_base*(rank+1) );

				nt_data_current.t_sample.at(xi_E_current.at(i)) =para_other.unassigned_time; // exclude the sample

				nt_data_current.current_size.at(xi_E_current.at(i)) = nt_data_current.current_size.at(xi_E_current.at(i)) - 1;
			}

			//----
		break;
		}
		case 1:{ // keep the original t_sample (could be unassigned_time)
			//nt_data_current.t_sample.at(xi_E_current.at(i)) = nt_data_current.t_sample.at(xi_E_current.at(i));
		break;
		}
	}
}

	myfile_out.open((string(path4)+string("seed_delete_seq.txt")).c_str(),ios::out);
	myfile_out << seed_intial_mcmc;
	myfile_out.close();
	

	myfile_out.open((string(path4)+string("p_sample.txt")).c_str(),ios::out);
	myfile_out << p_sample;
	myfile_out.close();

	int total_sample =0;
	for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections
		if (nt_data_current.t_sample.at(xi_E_current.at(i))!=para_other.unassigned_time) total_sample = total_sample+1;
	}
	double p_sample_actual = ((double)  total_sample)/ ((double) xi_E_current.size() );
	myfile_out.open((string(path4)+string("p_sample_actual.txt")).c_str(),ios::out);
	myfile_out << p_sample_actual;
	myfile_out.close();

	myfile_out.open((string(path4)+string("xi_EnIS_initial.txt")).c_str(),ios::out);
	myfile_out << "k" << endl;
	if (xi_EnIS_current.empty()!=1){
	for (int i=0; i<=((int)xi_EnIS_current.size()-1);i++){
	myfile_out << xi_EnIS_current.at(i) << endl;
	}
	}
	myfile_out << "size" << endl;
	myfile_out << xi_EnIS_current.size();
	myfile_out.close();

	myfile_out.open((string(path4)+string("t_sample_initial.txt")).c_str(),ios::app);
	for (int i=0; i<=(para_other.n-1);i++){
	myfile_out << nt_data_current.t_sample.at(i) << endl;
	}
	myfile_out.close();



//-----------------------

//----- this part works for deleting all samples-----//

/*
for (int i=0; i<=(int) xi_E_current.size() -1 ; i++){

	int subject= xi_E_current.at(i);

	switch(nt_data_current.t_sample.at(subject)==para_other.unassigned_time){
		case 0:{// with sample, delete the sample and udpate current_size s correspondingly
				int rank = distance( nt_data_current.t_nt.at(subject).begin(), find(nt_data_current.t_nt.at(subject).begin(),nt_data_current.t_nt.at(subject).end(),nt_data_current.t_sample.at(subject) ) ); 

				nt_data_current.t_nt.at(subject).erase(nt_data_current.t_nt.at(subject).begin() + rank);
				nt_data_current.nt.at(subject).erase(nt_data_current.nt.at(subject).begin()+para_other.n_base*rank , nt_data_current.nt.at(subject).begin()+para_other.n_base*(rank+1) );

				nt_data_current.current_size.at(subject) = nt_data_current.current_size.at(subject) - 1;

				nt_data_current.t_sample.at(subject)=para_other.unassigned_time;

		break;
		}
		case 1:{// without sample, keep it
		break;
		}
	}
}

xi_EnIS_current.clear() ; // this would be empty as initial value

xi_EnIS_current = xi_EnI_current;

//--------------

	myfile_out.open((string(path4)+string("t_sample_initial.txt")).c_str(),ios::app);
	for (int i=0; i<=(para_other.n-1);i++){
	myfile_out << nt_data_current.t_sample.at(i) << endl;
	}
	myfile_out.close();

	myfile_out.open((string(path4)+string("xi_EnIS_initial.txt")).c_str(),ios::out);
	myfile_out << "k" << endl;
	if (xi_EnIS_current.empty()!=1){
	for (int i=0; i<=((int)xi_EnIS_current.size()-1);i++){
	myfile_out << xi_EnIS_current.at(i) << endl;
	}
	}
	myfile_out.close();

	myfile_out.open((string(path4)+string("current_size_initial.txt")).c_str(),ios::app);
	for (int i=0; i<=(para_other.n-1);i++){
	myfile_out << nt_data_current.current_size.at(i) << endl;
	}
	myfile_out.close();
*/

}

//------------------------------------------------------------------------//

void mcmc_UPDATE::source_t_e_update_V2(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg, vec2int& infecting_list_current_arg, vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, int& subject_proposed, vector<int>& list_update, int iter, gsl_rng* & r_c){ // do not update t_e compare to source_t_e_update 

//double t_back =10.0;

double acp_pr = 0.0;
//double t_low, t_up;

double log_pr_forward=0.0; 
double log_pr_backward=0.0;

// double log_pr_t_e_forward=0.0;  
// double log_pr_t_e_backward=0.0;

double log_pr_seq_forward=0.0; 
double log_pr_seq_backward=0.0;

double log_pr_ds_forward=0.0;
double log_pr_ds_backward=0.0;


double t_proposed = t_e_arg.at(subject_proposed); // not changing t_e

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

//vector< vector<double> > delta_mat_modified = delta_mat_current_arg;

//vector<double> t_e_modified = t_e_arg;
//vector <int> index_modified = index_arg;
//vector <int> xi_E_minus_modified = xi_E_minus_arg;

// vector <int> xi_U_modified = xi_U_arg;
// vector <int> xi_E_modified = xi_E_arg;
// vector <int> xi_EnI_modified = xi_EnI_arg;

vector<int> current_size_modified = current_size_arg;

vec2int nt_modified = nt_current_arg;
vec2d t_nt_modified = t_nt_current_arg;

vec2int infecting_list_modified= infecting_list_current_arg;
vector<int> infecting_size_modified = infecting_size_current_arg;

// const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
// gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
// gsl_rng_set (r_c,iter); // set a seed

//vector<int> infected_source_modified = infected_source_current_arg;

vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

//vector<int> nt_subject_seq; // the orginal first sequence of the subject
//nt_subject_seq.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);

int source_x = infected_source_current_arg.at(subject_proposed);

int rank_source_x; // the rank of removed sequence in old source (note: it has different meaning in function t_e_seq in which the source remains the same)
vector<int> nt_current_source_x;
vector<double> t_nt_current_source_x;
vector<int> nt_modified_source_x;
vector<double> t_nt_modified_source_x;



//-------------- propose a new source --------------//

int source_y;
int rank_source_y;
vector<int> nt_current_source_y;
vector<double> t_nt_current_source_y;
vector<int> nt_modified_source_y;
vector<double> t_nt_modified_source_y;

vector<int> source_pool; // vector contains the indices of possible source for the subject

//double t_bound = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
double t_bound = t_e_arg.at(subject_proposed);

for (int i=0;i<=(int)(xi_I_arg.size()-1);i++){

//	switch( t_i_arg.at(xi_I_arg.at(i))<t_bound){
//	switch( (t_i_arg.at(xi_I_arg.at(i))<t_bound) & ((t_bound  - t_r_arg.at(xi_I_arg.at(i)))<=t_back) ){
	switch( (t_i_arg.at(xi_I_arg.at(i))<=t_bound) &  (t_r_arg.at(xi_I_arg.at(i))>t_bound)){
		case 1:{
			source_pool.push_back(xi_I_arg.at(i));	
		break;
		}
		case 0:{
		break;
		}
	}
}
	

source_pool.insert(source_pool.begin(),9999);

int num_infectious = (int)source_pool.size();

//-----------------------------propose uniformly-------------------------------------------//

//source_y = source_pool.at(gsl_rng_uniform_int(r_c, num_infectious)); // uniformly choose a new source (including bg)

//-----propose according to infectious challenge--------------------//

vector<double> ic(num_infectious);
ic.at(0) = para_current_arg.alpha;

switch(num_infectious>=2){

	case 1:{ // with 2nd sources from pool

		for (int j=1;j<=(num_infectious-1);j++){
		ic.at(j)= para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_pool.at(j)]/norm_const_current_arg.at(source_pool.at(j)); // a new source will be proposed according to the infectious challenges
		}

		double *P=&ic.at(0); // convert vector to array
		gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
		int link= gsl_ran_discrete (r_c, g);
		gsl_ran_discrete_free (g);

		source_y = source_pool.at(link); // a new source
		log_pr_forward = log(ic.at(link));

		switch(source_x==9999){
			case 0:{
				double ic_source_x =  para_current_arg.beta*stb_arg.at(subject_proposed)*kernel_mat_current_arg[subject_proposed][source_x]/norm_const_current_arg.at(source_x);
				log_pr_backward =log(ic_source_x);
			break;
			}
			case 1:{
				double ic_source_x =  para_current_arg.alpha;
				log_pr_backward =log(ic_source_x);
			break;
			}
		}


		break;
	}

	case 0:{ // only primary source from pool

		source_y = 9999;
		log_pr_forward = log(para_current_arg.alpha);

		double ic_source_x =  para_current_arg.alpha;
		log_pr_backward =log(ic_source_x);

	break;
	}
}


//--------------//

//----------end of proposing a new source -----//

//----------------------------------------------------------------------------------------------------------------//

switch(source_y==source_x){

	case 0:{


		//-------- propose a new t_e---------------//
// 		switch(source_y){
// 		
// 		case 9999:{ // by background
// 		
// 			t_up = min(t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
// 			t_low = max(0.0, t_up-t_back);
// 		
// 			
// 			t_proposed= gsl_ran_flat(r_c,t_low, t_up );
// 			
// 			log_pr_t_e_forward = log(1.0/(t_up-t_low));
// 		
// 		break;
// 		}
// 		
// 		default :{ // not by background
// 		
// 			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_y)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));			
// 			t_low = max(t_i_arg.at(source_y), t_up -t_back );
// 			
// 			t_proposed= gsl_ran_flat(r_c, t_low, t_up );
// 
// 			log_pr_t_e_forward = log(1.0/(t_up-t_low));
// 
// 		break;
// 		}
// 		
// 		}

		//------------//
//		t_e_modified.at(subject_proposed) = t_proposed;
		//------------//

// 		switch(source_x){
// 		
// 		case 9999:{ // by background
// 		
// 			t_up = min( t_sample_arg.at(subject_proposed), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
// 			t_low = max(0.0, t_up-t_back);
// 						
// 			log_pr_t_e_backward = log(1.0/(t_up-t_low));
// 		
// 		break;
// 		}
// 		
// 		default :{ // not by background
// 		
// 			t_up = min(min(t_sample_arg.at(subject_proposed), t_r_arg.at(source_x)), min(t_i_arg.at(subject_proposed), t_max_CUPDATE));
// 			t_low = max(t_i_arg.at(source_x), t_up -t_back );
// 
// 			log_pr_t_e_backward = log(1.0/(t_up-t_low));
// 
// 		break;
// 		}
// 		
// 		}


		//--------------------------------------//

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("00test.txt")).c_str(),ios::app);
// myfile_mcmc_out << 0 << endl;
// myfile_mcmc_out.close();		
		
		
		vector<int> seq_proposed(n_base_CUPDATE); // newly proposed sequence when source changes
		vector<int> nt_past_forward(n_base_CUPDATE);
		vector<int> nt_future_forward(n_base_CUPDATE);
		double t_past, t_future;
		
		vector<int> seq_proposed_backward(n_base_CUPDATE);
		vector<int> nt_past_backward(n_base_CUPDATE);
		vector<int> nt_future_backward(n_base_CUPDATE);
		double t_proposed_backward, t_past_backward, t_future_backward;
		
		//------------------------------------------------
		
		switch(source_y==9999){
		
			case 0:{// new 2nd infection

				infecting_size_modified.at(source_y) = infecting_size_modified.at(source_y) +1;
				
				vector<double> t_y(infecting_size_current_arg.at(source_y));
				for (int i=0;i<=(infecting_size_current_arg.at(source_y)-1);i++){
				t_y.at(i) = t_e_arg.at(infecting_list_current_arg[source_y][i]);
				}
				t_y.push_back(t_proposed);
				sort(t_y.begin(), t_y.end());
			
				int rank_y = distance(t_y.begin(), find(t_y.begin(), t_y.end(), t_proposed));
				infecting_list_modified.at(source_y).insert(infecting_list_modified.at(source_y).begin()+rank_y, subject_proposed);
				
				//----------------------------------------------------//
	
				nt_current_source_y = nt_current_arg.at(source_y);
				t_nt_current_source_y = t_nt_current_arg.at(source_y);
		
				nt_modified_source_y = nt_current_source_y;
				t_nt_modified_source_y = t_nt_current_source_y;
		
				t_nt_modified_source_y.push_back(t_proposed);
		
				sort( t_nt_modified_source_y.begin(),  t_nt_modified_source_y.end()); 
		
				rank_source_y = distance(t_nt_modified_source_y.begin(), find(t_nt_modified_source_y.begin(),t_nt_modified_source_y.end(), t_proposed) );
		
				current_size_modified.at(source_y) = current_size_modified.at(source_y) + 1;
		
				//----------------------------------------------------------------------------------------------------------------//
				t_past = t_nt_modified_source_y.at(rank_source_y-1);
				nt_past_forward.assign(nt_current_source_y.begin()+(rank_source_y-1)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE);
		
				//t_proposed = t_e_subject;
		
				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y (take closer seq as the future seq)
		
								switch(t_nt_current_arg[subject_proposed][1]<t_nt_modified_source_y.at(rank_source_y +1)){
								//switch(t_sample_arg.at(subject_proposed)<t_nt_modified_source_y.at(rank_source_y +1)){
		
									case 1:{// take the one from subject as future seq

										t_future  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
																				
									seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
		
									case 0:{ // take the one from source_y as future seq
										t_future = t_nt_modified_source_y.at(rank_source_y+1);
										nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
		
									break;
									}
								}
		
							break;
							}
		
							case 0:{//  inserted seq will be last seq at rank_source_y  in source_y
		
										t_future  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
	
										seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
		
						switch(current_size_modified.at(source_y)>(rank_source_y+1)){
		
							case 1:{// inserted seq will NOT be last seq at rank_source_y  in source_y
								t_future = t_nt_modified_source_y.at(rank_source_y+1);
								nt_future_forward.assign(nt_current_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_current_source_y.begin()+(rank_source_y+1)*n_base_CUPDATE);
		
								seq_propose_cond(seq_proposed, log_pr_seq_forward, nt_past_forward,  nt_future_forward,  t_proposed, t_past, t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);						
							break;
							}
		
							case 0:{// inserted seq will be last seq at rank_source_y  in source_y
								t_future = t_proposed;
								seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				//nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, nt_subject_seq.begin(), nt_subject_seq.end());
				nt_modified_source_y.insert(nt_modified_source_y.begin()+(rank_source_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());
		
				nt_modified.at(source_y) = nt_modified_source_y;
				t_nt_modified.at(source_y) = t_nt_modified_source_y;

			break;
			}
		
			case 1:{// new bg infection

				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
		
					case 1:{// with 2nd seq in subject

						t_past  = t_nt_current_arg[subject_proposed][1]; // yes, t_past!
						nt_past_forward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE);
		
						t_future = t_proposed;

						seq_propose_uncond( seq_proposed, log_pr_seq_forward,  nt_past_forward, t_proposed, t_past,  t_future, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);			
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
		
						log_pr_seq_forward = n_base_CUPDATE*log(0.25);
		
						for (int i=0; i<=(n_base_CUPDATE-1);i++){
						seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
						}
		
					break;
					}
				}
			break;
			}
		
		}
		
		//--------------------------------------------
		nt_modified_subject.erase(nt_modified_subject.begin() , nt_modified_subject.begin()+n_base_CUPDATE); 
		nt_modified_subject.insert(nt_modified_subject.begin(), seq_proposed.begin(), seq_proposed.end());

		t_nt_modified_subject.erase(t_nt_modified_subject.begin());
		t_nt_modified_subject.push_back(t_proposed); 
		sort( t_nt_modified_subject.begin(),  t_nt_modified_subject.end()); 

		nt_modified.at(subject_proposed) = nt_modified_subject;
		t_nt_modified.at(subject_proposed) = t_nt_modified_subject;
		//---------------------------------------------
		
		switch(source_x==9999){
		
			case 0:{// was 2nd infection

				infecting_size_modified.at(source_x) = infecting_size_modified.at(source_x) - 1;

				int rank_x = distance(infecting_list_current_arg.at(source_x).begin(), find(infecting_list_current_arg.at(source_x).begin(), infecting_list_current_arg.at(source_x).end(), subject_proposed));
			
				infecting_list_modified.at(source_x).erase(infecting_list_modified.at(source_x).begin()+rank_x);

				//----------------------------------------------------------------------------------------------------------------//

				nt_current_source_x = nt_current_arg.at(source_x);
				t_nt_current_source_x = t_nt_current_arg.at(source_x);
		
				rank_source_x = distance( t_nt_current_source_x.begin(), find(t_nt_current_source_x.begin(), t_nt_current_source_x.end(), t_e_arg.at(subject_proposed)) ); 
		
				//----------------------------------------------------------------------------------------------------------------//
		
				t_past_backward = t_nt_current_source_x.at(rank_source_x-1);
				nt_past_backward.assign(nt_current_source_x.begin()+(rank_source_x-1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x)*n_base_CUPDATE);
		
				t_proposed_backward =t_e_arg.at(subject_proposed);
				seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
				switch(current_size_arg.at(subject_proposed)>1){ 
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){
				
					case 1:{// with2nd seq subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
		
							case 1:{// also NOT last seq at rank_source_x  in source_x (take closer seq as the future seq)
								switch(t_nt_current_arg[subject_proposed][1]<t_nt_current_source_x.at(rank_source_x +1)){
								//switch(t_sample_arg.at(subject_proposed)<t_nt_current_source_x.at(rank_source_x +1)){

									case 1:{// take the one from subject as future seq

										t_future_backward  = t_nt_current_arg[subject_proposed][1]; 
										nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
						
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

									break;
									}
									case 0:{// take the one from source_x as future seq
										t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
										nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
						
										seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
									break;
									}
								}
							break;
							}
		
							case 0:{ //last seq at rank_source_x  in source_x

								t_future_backward  = t_nt_current_arg[subject_proposed][1]; // always takes the 2nd seq in subject as the future sequence (both for forward and backward)
								nt_future_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}	
						}
		
		
					break;
					}
					
					case 0:{ // with no 2nd seq in the subject
		
						switch(current_size_arg.at(source_x)>(rank_source_x+1)){
							case 1:{ // not last seq at rank_source_x  in source_x
								t_future_backward  = t_nt_current_source_x.at(rank_source_x +1);
								nt_future_backward.assign(nt_current_source_x.begin()+(rank_source_x+1)*n_base_CUPDATE, nt_current_source_x.begin()+(rank_source_x+2)*n_base_CUPDATE);
				
								seq_backward_pr_cond(seq_proposed_backward, log_pr_seq_backward, nt_past_backward, nt_future_backward, t_proposed_backward, t_past_backward,  t_future_backward, para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
								
							break;
							}
							case 0:{ // last seq rank_source_x  in source_x
								t_future_backward = t_proposed_backward;
		
								seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							break;
							}
						}
		
					break;
					}
				}
		
				//----------------------------------------------------------------------------------------------------------------//
		
				nt_modified_source_x = nt_current_source_x;
				t_nt_modified_source_x = t_nt_current_source_x;
				
				t_nt_modified_source_x.erase(t_nt_modified_source_x.begin() + rank_source_x); // erase the original t_nt entry for source_x
				nt_modified_source_x.erase(nt_modified_source_x.begin()+n_base_CUPDATE*rank_source_x , nt_modified_source_x.begin()+n_base_CUPDATE*(rank_source_x+1) );  //erase the original nt entry for source_x

				nt_modified.at(source_x) = nt_modified_source_x;
				t_nt_modified.at(source_x) = t_nt_modified_source_x;


				current_size_modified.at(source_x) = current_size_modified.at(source_x) - 1;
		
			break;
			}
		
			case 1:{// was bg infection
		
				switch(current_size_arg.at(subject_proposed)>1){
				//switch(t_sample_arg.at(subject_proposed)!=unassigned_time_CUPDATE){

					case 1:{// with 2nd seq in subject

						t_past_backward = t_nt_current_arg[subject_proposed][1]; // note: it is "past"!
						nt_past_backward.assign(nt_current_arg.at(subject_proposed).begin()+ n_base_CUPDATE, nt_current_arg.at(subject_proposed).begin()+ 2*n_base_CUPDATE );
		
						t_proposed_backward = t_e_arg.at(subject_proposed);
						seq_proposed_backward.assign(nt_current_arg.at(subject_proposed).begin(), nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE);
		
						t_future_backward = t_proposed_backward;

						seq_backward_pr_uncond(seq_proposed_backward,  log_pr_seq_backward, nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
		
		
					break;
					}
		
					case 0:{// with no 2nd seq in subject
						log_pr_seq_backward = n_base_CUPDATE*log(0.25); // note: this is the proposal density! not refers to the model!
					break;
					}
				}
		
			break;
			}
		
		}
		
		

		//------------- deal with change of likelihood due to change of source and sequences in subject_proposed)---------------//
		
		switch (current_size_arg.at(subject_proposed)>1) {
		
		case 1:{
		
		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below
		
		lh_square_modified.log_f_S.at(subject_proposed) = 0.0;
		
		for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){

			vector<int> seq_1(nt_modified.at(subject_proposed).begin()+j*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE));
			vector<int> seq_2(nt_modified.at(subject_proposed).begin()+(j+1)*(n_base_CUPDATE), nt_modified.at(subject_proposed).begin()+(j+2)*(n_base_CUPDATE));
			
			lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2,t_nt_modified_subject.at(j), t_nt_modified_subject.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);	
		}
		
		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 
		
		break;
		}
		
		default:{
		break;
		}
		
		}


		//-----------------------------------------------------------------------------------------//

// 		lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
// 		
// 		for (int j=0;j<=(int) (xi_I_arg.size()-1);j++){ // update delta_mat; and pre-assign k_sum_E & kt_sum_E which might be changed later again
// 		
// 			if (t_i_arg.at(xi_I_arg.at(j))<t_proposed) {
// 			
// 				switch (t_r_arg.at(xi_I_arg.at(j))>=t_proposed) {
// 					case 1:{
// 						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_proposed - t_i_arg.at(xi_I_arg.at(j));
// 					
// 					break;
// 					}
// 					case 0:{
// 						delta_mat_modified[subject_proposed][xi_I_arg.at(j)] = t_r_arg.at(xi_I_arg.at(j)) - t_i_arg.at(xi_I_arg.at(j));
// 					break;
// 					}
// 				}
// 			lh_square_modified.kt_sum_E.at(subject_proposed) = lh_square_modified.kt_sum_E.at(subject_proposed) + delta_mat_modified[subject_proposed][xi_I_arg.at(j)]*kernel_mat_current_arg[subject_proposed][xi_I_arg.at(j)]/norm_const_current_arg.at(xi_I_arg.at(j));
// 	
// 			}
// 		}

		//----------

		lh_square_modified.k_sum_E.at(subject_proposed)=0.0;

		switch(source_y){
		
		case 9999:{ // by background
		lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
		break;
		}
		
		default :{ // not by background
		lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y);
		lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
		break;
		}
		
		}
		//----------

		log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below

		lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
	
		log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 

// 		switch (find(index_arg.begin(),index_arg.end(),subject_proposed)==index_arg.end()) { // return 1 if proposed subject not one of the original indexes
// 		
// 			case 1:{
// 
// 				log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); //subtract part of likelihood that would be updated below
// 
// 				lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
// 			
// 				log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); 
// 
// 			break;
// 			}
// 		
// 		
// 			case 0: { // when chosen subject is one of the indexes
// 
// 				//log_lh_modified = log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed)); // this subject would 
// 				lh_square_modified.k_sum_E.at(subject_proposed)=0.0;
// 				lh_square_modified.kt_sum_E.at(subject_proposed)=0.0;
// 				lh_square_modified.q_E.at(subject_proposed)=0.0;
// 				lh_square_modified.g_E.at(subject_proposed)=1.0;
// 				lh_square_modified.h_E.at(subject_proposed)=1.0;
// 				lh_square_modified.f_E.at(subject_proposed)=1.0;
// 
// 				log_lh_modified = log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed)); // this subject would 
// 			
// 			break;
// 			}
// 		
// 		}
// // 		
		//--------------------//


// 		switch ( find(xi_I_arg.begin(), xi_I_arg.end(),subject_proposed) != (xi_I_arg.end()) ) { //return 1 when the subject is also in xi_I
// 			case 1:{
// 		
// 				log_lh_modified = log_lh_modified - log(lh_square_modified.f_I.at(subject_proposed)); //subtract part of likelihood that would be updated below
// 				lh_square_modified.f_I.at(subject_proposed) = func_latent_pdf(t_i_arg.at(subject_proposed) - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
// 				log_lh_modified = log_lh_modified + log(lh_square_modified.f_I.at(subject_proposed)); 
// 				
// 			break;
// 			}
// 
// 			case 0:{
// 			break;
// 			}
// 		}
// 		
// 		//----------
// 		
// 		switch ( find(xi_EnI_arg.begin(), xi_EnI_arg.end(),subject_proposed) != (xi_EnI_arg.end()) ) { //return 1 when the subject is also in xi_EnI
// 			case 1:{
// 			
// 				log_lh_modified = log_lh_modified - log(lh_square_modified.f_EnI.at(subject_proposed)); //subtract part of likelihood that would be updated below
// 				lh_square_modified.f_EnI.at(subject_proposed) = 1.0 - func_latent_cdf( t_max_CUPDATE - t_proposed, para_current_arg.mu_lat, para_current_arg.var_lat);
// 				log_lh_modified = log_lh_modified + log(lh_square_modified.f_EnI.at(subject_proposed)); 
// 				
// 			break;
// 			}
// 			case 0:{
// 			break;
// 			}
// 		}


		//-----------------------------------------------------------------------------------------//

		switch(source_x){
		
			case 9999:{ // was background infection
		
				switch(source_y){
					case 9999:{ // new  bg infection 		
					break;
					}
		
					default:{ // new secondary infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); 
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) = 0.0;
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); // in fact redudancy, but here for clarity
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
			break;
			}
			
			default :{ // was secondary infection
		
				switch(source_y){
					case 9999:{ // new  bg infection
		
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); // redudant as second term must be zero
		
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = 0.0; // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) =  para_current_arg.alpha;
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						lh_square_modified.log_f_Snull.at(subject_proposed) =  n_base_CUPDATE*log(0.25);
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));

						log_lh_modified =  log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed);
		
						//--
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
			
		
				
					break;
					}
		
					default:{ // new secondary infection
			
						//log_lh_modified =  log_lh_modified - log(lh_square_modified.f_E.at(subject_proposed));
		
						//lh_square_modified.k_sum_E.at(subject_proposed) = kernel_mat_current_arg[subject_proposed][source_y]/norm_const_current_arg.at(source_y); // update k_sum_E
						//lh_square_modified.g_E.at(subject_proposed) = para_current_arg.beta*stb_arg.at(subject_proposed)*lh_square_modified.k_sum_E.at(subject_proposed);
						//lh_square_modified.f_E.at(subject_proposed) = lh_square_modified.g_E.at(subject_proposed)*lh_square_modified.h_E.at(subject_proposed);
		
						//log_lh_modified =  log_lh_modified + log(lh_square_modified.f_E.at(subject_proposed));
		
						//--
		
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_x); // source_x must had more than or equal to 2 sequences
					
						lh_square_modified.log_f_S.at(source_x) = 0.0;
		
						switch(current_size_modified.at(source_x)>1){// only have to count log_f_S if there are more than or equal to 2 seq left in source_x
							case 1:{		
								for (int j=0;j<=(current_size_modified.at(source_x)-2);j++){
				
								vector<int> seq_1(nt_modified_source_x.begin()+j*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE));
								vector<int> seq_2(nt_modified_source_x.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_x.begin()+(j+2)*(n_base_CUPDATE));
						
								lh_square_modified.log_f_S.at(source_x) =lh_square_modified.log_f_S.at(source_x)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_x.at(j), t_nt_modified_source_x.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
							
								}
							break;
							}
		
							case 0:{
							break;
							}
						}			
		
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_x); // 2nd term migt be zero depends on if there are more than or equal to 2 seq left in source_x
		
						//---
		
						log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(source_y); // if source_y had one seq only, log_f_s would be zero anyway
					
						lh_square_modified.log_f_S.at(source_y) = 0.0;
					
						for (int j=0;j<=(current_size_modified.at(source_y)-2);j++){
		
						vector<int> seq_1(nt_modified_source_y.begin()+j*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE));
						vector<int> seq_2(nt_modified_source_y.begin()+(j+1)*(n_base_CUPDATE), nt_modified_source_y.begin()+(j+2)*(n_base_CUPDATE));
				
						lh_square_modified.log_f_S.at(source_y) =lh_square_modified.log_f_S.at(source_y)+log_lh_seq(seq_1, seq_2, t_nt_modified_source_y.at(j), t_nt_modified_source_y.at(j+1), para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);
					
						}
					
						log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(source_y); 
			
		
					break;
					}
				}
		
		
			break;
			}
		
		}
		
		
		
	//------------- end of with change of likelihood (due to change of source and sequences in subject_proposed)---------------//
	
		//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward)*exp(log_pr_seq_backward-log_pr_seq_forward)*exp(log_pr_ds_backward-log_pr_ds_forward));



		//acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+(log_pr_seq_backward-log_pr_seq_forward) +(log_pr_t_e_backward-log_pr_t_e_forward)+ (log_pr_ds_backward-log_pr_ds_forward)) );
		
		acp_pr = min(1.0,exp( (log_lh_modified-log_lh_current_arg) + (log_pr_backward-log_pr_forward)+(log_pr_seq_backward-log_pr_seq_forward) + (log_pr_ds_backward-log_pr_ds_forward)) );
				
		double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);
		
		switch(uniform_rv<=acp_pr){
		case 1: {
		
		lh_square_current_arg = lh_square_modified;
		log_lh_current_arg = log_lh_modified;
		
		nt_current_arg = nt_modified;
		t_nt_current_arg = t_nt_modified;

// 		delta_mat_current_arg = delta_mat_modified;
// 		t_e_arg= t_e_modified;
// 		index_arg = index_modified;
// 		xi_E_minus_arg = xi_E_minus_modified;

		current_size_arg = current_size_modified;
		infected_source_current_arg.at(subject_proposed) =  source_y;
		//infected_source_current_arg = infected_source_modified;

		infecting_list_current_arg = infecting_list_modified;
		infecting_size_current_arg = infecting_size_modified;

	break;
	}

	case 0: {
	break;
	}
}

break;
}

case 1:{ // source_y==source_x
break;
}
}

//gsl_rng_free(r_c);

//--------------------------------------------------

// ofstream myfile_mcmc_out; 
// 
// myfile_mcmc_out.open((string(path4)+string("subect_proposed_source_te_update.txt")).c_str(),ios::app);
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<subject_proposed << ","<<  source_x <<","<< source_y << "," << infecting_size_current_arg.at(subject_proposed) << ","<< (int) list_update.size() << endl;
// myfile_mcmc_out.close();	
// 
// myfile_mcmc_out.open((string(path4)+string("log_lh_change_source_te_update.txt")).c_str(),ios::app);
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_lh_current_arg <<","<< log_lh_modified << "," << source_x <<","<< source_y <<","<< acp_pr << endl;
// myfile_mcmc_out.close();	
// 
// myfile_mcmc_out.open((string(path4)+string("log_pr_forward_backward_source_te_update.txt")).c_str(),ios::app);
// //if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) )
// //myfile_mcmc_out << exp(log_lh_modified-log_lh_current_arg) <<"," << exp(log_pr_backward-log_pr_forward) << ","<<  exp(log_pr_seq_backward-log_pr_seq_forward) << "," << acp_pr << endl;
// if ( (source_x!=source_y) & (log_lh_current_arg == log_lh_modified) ) myfile_mcmc_out <<log_pr_seq_backward <<","<< log_pr_seq_forward<< "," << log_pr_ds_backward <<","<< log_pr_ds_forward << endl;
// myfile_mcmc_out.close();	


}

/*------------------------------------------------*/


void mcmc_UPDATE::index_first_seq(lh_SQUARE& lh_square_current_arg, double& log_lh_current_arg, const vector< vector<double> >& kernel_mat_current_arg, vector< vector<double> >& delta_mat_current_arg, vector<int>& xi_U_arg, vector<int>& xi_E_arg, vector<int>& xi_E_minus_arg, const vector<int>& xi_I_arg, vector<int>& xi_EnI_arg, const vector<double>& t_r_arg, const vector<double>& t_i_arg, vector<double>& t_e_arg, vector<int>& index_arg, const para_key& para_current_arg, const vector<double>& stb_arg, const vector<double>& norm_const_current_arg, const vector<int>& infected_source_current_arg, const vector<double>& t_sample_arg, const vector<int>& current_size_arg, vec2int& nt_current_arg , vec2d& t_nt_current_arg,  vec2int& infecting_list_current_arg, const vector<int>& infecting_size_current_arg, vector<int>&  xi_beta_E_arg, vector<int>& con_seq,int& subject_proposed,int iter,gsl_rng* & r_c){

//double t_back =10.0;

//int subject_proposed ;

//double t_low, t_up;
double acp_pr = 0.0;

lh_SQUARE lh_square_modified = lh_square_current_arg;
double log_lh_modified =  log_lh_current_arg;

//vector< vector<double> > delta_mat_modified = delta_mat_current_arg;
//vector<double> t_e_modified = t_e_arg;
//vector <int> index_modified = index_arg;
//vector <int> xi_E_minus_modified = xi_E_minus_arg;

// vector <int> xi_U_modified = xi_U_arg;
// vector <int> xi_E_modified = xi_E_arg;
// vector <int> xi_EnI_modified = xi_EnI_arg;

/*const gsl_rng_type* T_c= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_rng_set (r_c,iter); // set a see*/



vector<int> nt_modified_subject = nt_current_arg.at(subject_proposed);
//vector<double> t_nt_modified_subject = t_nt_current_arg.at(subject_proposed);

vector<int> nt_current_seq; // the orginal first sequence of the subject which would be updated
vector<int> nt_second_seq; // the second seq of subject (if available)

//int subject_source = infected_source_current_arg.at(subject_proposed);


//int rank_subject_x =distance( t_nt_current_arg.at(subject_proposed).begin(), find(t_nt_current_arg.at(subject_proposed).begin(), t_nt_current_arg.at(subject_proposed).end(), t_e_arg.at(subject_proposed)) ); //count the rank (distance from the first element) of the original t_e among t_nt_current_arg.at(subject_proposed) 
int rank_subject_x =0; // it is always zero as we are updating the sequence at its own infection

//t_nt_modified_subject.erase(t_nt_modified_subject.begin() + rank_subject_x); // erase the original t_nt entry for subject_proposed

//---

nt_current_seq.assign(nt_modified_subject.begin()+n_base_CUPDATE*rank_subject_x , nt_modified_subject.begin()+n_base_CUPDATE*(rank_subject_x+1) ); //copy the original nt before erasing

nt_modified_subject.erase(nt_modified_subject.begin()+n_base_CUPDATE*rank_subject_x , nt_modified_subject.begin()+n_base_CUPDATE*(rank_subject_x+1) );  //erase the original nt entry for subject_proposed

//---------------------------------------- proposing a new sequence & the proposal probability ----------------------------------------------------------//

vector<int> seq_proposed(n_base_CUPDATE);

//double dt;

double t_proposed = 0.0; // not really gonna be used

double t_past, t_future;

double log_pr_forward=0.0; // the log of proposal probability 

vector<int> nt_past_forward(n_base_CUPDATE); // the sequence at the nearest past (in the direction of time change) compared to the time of the proposed sequence; this might be or might not be the original sequence which gotta be replaced
vector<int> nt_future_forward(n_base_CUPDATE); // the sequence at the nearest future(in the direction of time change) compared to the time of the proposed sequence

//dt = t_proposed - t_e_arg.at(subject_proposed); // the dimension of dt tells the direction of time change

//

double t_proposed_backward = 0.0;

vector<int> seq_proposed_backward(n_base_CUPDATE);

double t_past_backward, t_future_backward;

double log_pr_backward=0.0; // the log of proposal probability 

vector<int> nt_past_backward(n_base_CUPDATE); 
vector<int> nt_future_backward(n_base_CUPDATE); 

//-----------------

switch(current_size_arg.at(subject_proposed)>1){ // return 1 when the subject has more than one sequence available

	case 0:{ //  ONLY one sequence available for the subject; XX not relevant as we assume #seq on index>=2 XX

// 		for (int i=0; i<=(n_base_CUPDATE-1);i++){
// 		seq_proposed.at(i) = gsl_rng_uniform_int(r_c, 4) +1;
// 		}
		
	break;
	}

	case 1:{ // MORE than one sequence available for the subject

		nt_second_seq.assign(nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE , nt_current_arg.at(subject_proposed).begin()+n_base_CUPDATE*2 );

		t_past =   t_nt_current_arg[subject_proposed][1];
		nt_past_forward =nt_second_seq;
		t_future = t_e_arg.at(subject_proposed);//=0.0
		
		seq_propose_uncond(seq_proposed,  log_pr_forward, nt_past_forward, t_proposed, t_past, t_future,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE, r_c);

		//---
		t_past_backward =  t_nt_current_arg[subject_proposed][1];
		nt_past_backward =  nt_second_seq;
		t_future_backward = t_e_arg.at(subject_proposed);//=0.0
		
		seq_proposed_backward = nt_current_seq;

		seq_backward_pr_uncond(seq_proposed_backward,  log_pr_backward,nt_past_backward, t_proposed_backward, t_past_backward, t_future_backward,  para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE); 
				
	break;
	}

}


//-------------------------------------------end of proposing a new sequence ---------------------------------------------------------------------------------//

int rank_subject_y = 0;

nt_modified_subject.insert(nt_modified_subject.begin()+(rank_subject_y)*n_base_CUPDATE, seq_proposed.begin(), seq_proposed.end());  //insert  the new nt 

//----------------------------------------------------------------------------------//

switch (current_size_arg.at(subject_proposed)>1) {

	case 1:{

		log_lh_modified = log_lh_modified - lh_square_modified.log_f_S.at(subject_proposed); //subtract part of likelihood that would be updated below

		lh_square_modified.log_f_S.at(subject_proposed) = 0.0;

		for (int j=0;j<=(current_size_arg.at(subject_proposed)-2);j++){


		vector<int> seq_1(nt_modified_subject.begin()+j*(n_base_CUPDATE), nt_modified_subject.begin()+(j+1)*(n_base_CUPDATE));
		vector<int> seq_2(nt_modified_subject.begin()+(j+1)*(n_base_CUPDATE), nt_modified_subject.begin()+(j+2)*(n_base_CUPDATE));

		lh_square_modified.log_f_S.at(subject_proposed) =lh_square_modified.log_f_S.at(subject_proposed) + log_lh_seq(seq_1, seq_2, t_nt_current_arg[subject_proposed][j],t_nt_current_arg[subject_proposed][j+1], para_current_arg.mu_1, para_current_arg.mu_2, n_base_CUPDATE);

		}

		log_lh_modified = log_lh_modified + lh_square_modified.log_f_S.at(subject_proposed); 

	break;
	}

	default:{
	break;
	}

}


//----------------------------------------------

log_lh_modified = log_lh_modified - lh_square_modified.log_f_Snull.at(subject_proposed); //subtract part of likelihood that would be updated below

lh_square_modified.log_f_Snull.at(subject_proposed) = lh_snull(con_seq, seq_proposed, para_current_arg.p_ber, n_base_CUPDATE); // compute the log pr a seq for background

log_lh_modified = log_lh_modified + lh_square_modified.log_f_Snull.at(subject_proposed); 
	
//----------------------------------------------------------------------------------//

//acp_pr = min(1.0,exp(log_lh_modified-log_lh_current_arg)*exp(log_pr_backward-log_pr_forward));
acp_pr = min(1.0,exp((log_lh_modified-log_lh_current_arg)+(log_pr_backward-log_pr_forward)));


double uniform_rv = gsl_ran_flat(r_c, 0.0, 1.0);

// ofstream myfile_mcmc_out; 
// myfile_mcmc_out.open((string(path4)+string("22.txt")).c_str(),ios::app);
// myfile_mcmc_out << uniform_rv << "," << acp_pr << endl;
// myfile_mcmc_out.close();

switch(uniform_rv<=acp_pr){
case 1: {
lh_square_current_arg = lh_square_modified;
// delta_mat_current_arg = delta_mat_modified;
log_lh_current_arg = log_lh_modified;
// t_e_arg= t_e_modified;
//index_arg = index_modified;
//xi_E_minus_arg = xi_E_minus_modified;

// nt_current_arg = nt_modified;
// t_nt_current_arg = t_nt_modified;

nt_current_arg.at(subject_proposed) = nt_modified_subject;
//t_nt_current_arg.at(subject_proposed) = t_nt_modified_subject;


/*	switch (subject_source){
		
		case 9999:{ // by background
		break;
		}
		
		default :{ // not by background
		nt_current_arg.at(subject_source) = nt_modified_source;
		t_nt_current_arg.at(subject_source) = t_nt_modified_source;
		infecting_list_current_arg.at(subject_source) = infecting_list_modified_source;		
		}
	}
		*/	

break;
}

case 0: {
break;
}
}

//gsl_rng_free(r_c);

//--------------
		



		// myfile_mcmc_out.open((string(path4)+string("find_xi_E_mibus.txt")).c_str(),ios::app);
		// for (int i=0; i<=(int)(index_modified.size()-1); i++){
		// myfile_mcmc_out << 1*(find(xi_E_minus_modified.begin(), xi_E_minus_modified.end(),index_modified.at(i))==xi_E_minus_modified.end()) << endl; // should always equals to 1, i.e., new index has been excluded from xi_E_minus
		// }
		// myfile_mcmc_out.close();
		// 

		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_kt_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.kt_sum_E.at(index_arg.at(0)) << endl; // shouls always be 0
		// myfile_mcmc_out.close();
		// 
		// myfile_mcmc_out.open((string(path4)+string("index_arg_k_sum_E.txt")).c_str(),ios::app);
		// myfile_mcmc_out << lh_square_modified.k_sum_E.at(index_arg.at(0)) << endl;  // shouls always be 0
		// myfile_mcmc_out.close();

//  		ofstream myfile_mcmc_out; 
// 
// 		myfile_mcmc_out.open((string(path4)+string("index_t_e_seq.txt")).c_str(),ios::app);
// 		if (index_arg.empty()==0){
// 		for (int i=0; i<=(int)(index_arg.size()-1); i++){
// 		myfile_mcmc_out << index_arg.at(i) << "," << infected_source_current_arg.at(index_arg.at(i))  << endl;
// 		}
// 		}
// 		myfile_mcmc_out.close();

// 		myfile_mcmc_out.open((string(path4)+string("log_pr_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_pr_backward << "," << log_pr_forward << "," << exp(log_pr_backward-log_pr_forward) <<endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("log_lh_change_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  log_lh_current_arg  << "," <<log_lh_modified << "," << acp_pr << endl;
// 		myfile_mcmc_out.close();
// 
// 		
// 		myfile_mcmc_out.open((string(path4)+string("subject_proposed_t_e_seq.txt")).c_str(),ios::app); 		
// 		myfile_mcmc_out << subject_proposed  << ","<< current_size_arg.at(subject_proposed) << "," << infected_source_current_arg.at(subject_proposed) <<endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("t_proposed_t_e_seq.txt")).c_str(),ios::app);
// 		myfile_mcmc_out <<  t_proposed << "," << t_e_modified.at(subject_proposed) << "," << t_e_arg.at(subject_proposed) << endl;
// 		myfile_mcmc_out.close();


// 		myfile_mcmc_out.open((string(path4)+string("p_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << p_1  << "," << p_2 << "," << p_3 << endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("m_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << m_1  << "," << m_2 << "," << m_3 << "," << (long double) pow( p_1,m_1)*(long double) pow( p_2, m_2)*(long double) pow( p_3,m_3)<< endl;
// 		myfile_mcmc_out.close();
// 
// 		myfile_mcmc_out.open((string(path4)+string("n_i.txt")).c_str(),ios::app);
// 		myfile_mcmc_out << n_1  << "," << n_2 << "," << n_3 << "," <<(long double)  pow( p_1,n_1)*(long double) pow( p_2, n_2)*(long double) pow( p_3,n_3)<< endl;
// 		myfile_mcmc_out.close();





		//---------------

}

/*------------------------------------------------*/

//------------------------ initialization of parameters/unobserved data for McMC ---------------------------------//

/* Assumptions: t_r known; t_i is known within a range */

void initialize_mcmc(para_key& para_current, para_aux& para_other,  vector<int>& xi_I_current, vector<int>& xi_U_current, vector<int>& xi_E_current, vector<int>& xi_E_minus_current, vector<int>& xi_R_current,  vector<int>& xi_EnI_current,  vector<int>& xi_EnIS_current,  vector<int>& xi_InR_current,  vector<double>& t_e_current, vector<double>& t_i_current, vector<double>& t_r_current, vector<int>& index_current, vector<double>& stb_current, vector<int>& gp_stb_current, vector<int>& infected_source_current, vector < vector<double> >& kernel_mat_current, vector <double>& norm_const_current, vec2int& sample_data, vector<double>& t_onset, nt_struct& nt_data_current, vector<int>& con_seq){


const gsl_rng_type* T= gsl_rng_default;  // T is pointer points to the type of generator
gsl_rng *r = gsl_rng_alloc (T); // r is pointer points to an object with Type T

int seed_intial_mcmc = 1;  //1,999,-1000,-10000,123456

gsl_rng_set (r, seed_intial_mcmc); // set a seed

ofstream myfile_out; 

myfile_out.open((string(path4)+string("seed_initial_mcmc.txt")).c_str(),ios::out);
myfile_out << seed_intial_mcmc;
myfile_out.close();


para_current.alpha = 0.0002; //initialization of parameter to be estimated
para_current.beta = 3.0 ; //initialization of parameter to be estimated

// para_current.mu_lat = para_current.a*para_current.b*1; //initialization of parameter to be estimated
// para_current.var_lat = para_current.a*para_current.b*para_current.b*1; //initialization of parameter to be estimated
para_current.mu_lat = 2.0; //initialization of parameter to be estimated
para_current.var_lat = 2.0; //initialization of parameter to be estimated


para_current.c = 4.0; //initialization of parameter to be estimated
para_current.d = 2.0; //initialization of parameter to be estimated

para_current.k_1 =0.1; //initialization of parameter to be estimated
para_current.k_2 = 10.0; //initialization of parameter to be estimated

para_current.mu_1 = 1e-04; //initialization of parameter to be estimated
para_current.mu_2 = 5e-05; //initialization of parameter to be estimated

para_current.p_ber =  0.05;
//para_current.p_ber = p_ber_estm;

para_current.stb_1 = 1.0; //initialization of parameter to be estimated
para_current.stb_2 = 1.0; //initialization of parameter to be estimated


//--------

nt_data_current.t_nt.resize(para_other.n);
nt_data_current.nt.resize(para_other.n);

nt_data_current.current_size.resize(para_other.n);
nt_data_current.infecting_size.resize(para_other.n);
nt_data_current.infecting_list.resize(para_other.n);

for (int i=0; i<= (int)(para_other.n-1);i++){
	nt_data_current.t_nt.at(i).clear();
	nt_data_current.nt.at(i).clear();

	nt_data_current.current_size.at(i)  = 0;

	nt_data_current.infecting_size.at(i) = 0; // initially infecting no one
	nt_data_current.infecting_list.at(i).clear();

}

////------------------------initalization of t_i------------------//

// double t_range_start = 0.1;
// 
// for (int i=0; i<= (int)(xi_I_current.size()-1);i++){// loop over infections
// 
// 	double t_low = max(0.0,t_onset.at(xi_I_current.at(i)) - t_range_start);
// 	double t_up = min(t_r_current.at(xi_I_current.at(i)), min(para_other.t_max, t_onset.at(xi_I_current.at(i)) + t_range_start));
// 
// 	t_i_current.at(xi_I_current.at(i)) =  gsl_ran_flat(r, t_low, t_up);
// }
//----------------------------------------------------------------------//

/*
////---- intialization of xi_E and xi_EnI and xi_EnIS & insert the sampled seq and sampling times for infected----////

//-- (below)if assume xi_E and xi_EnI and  xi_EnIS are known --//
// xi_E_current = xi_E_current;
// xi_EnI_current = xi_EnI_current;
// xi_EnIS_current = xi_EnIS_current;


//-- (below) if assume xi_E and xi_EnI are known  based on whether or not sampled --//

xi_E_current =xi_I_current; // individuals gone through I (assumed known here) would be initialized as infected; sampled would also be infected (see below); assume index has gone through I

xi_EnI_current.clear() ;

for (int i=0; i<= (int)(para_other.n-1);i++){

	switch((nt_data_current.t_sample.at(i)!=para_other.unassigned_time)&(find(xi_I_current.begin(),xi_I_current.end(),i)==xi_I_current.end())){

		case 1:{ // sampled(infected) but not in xi_I_current

// 			nt_data_current.current_size.at(i) = 2;//  (initial size) infection + the sample
// 			nt_data_current.t_nt.at(i).push_back(nt_data_current.t_sample.at(i));
// 			nt_data_current.nt.at(i).insert(nt_data_current.nt.at(i).begin(), sample_data.at(i).begin(), sample_data.at(i).end());

			xi_E_current.push_back(i);
			xi_EnI_current.push_back(i);

		break;
		}

		case 0:{
		break;
		}

	}

}

xi_EnIS_current.clear() ; // this would be empty as initial value

//sort(xi_E_current.begin(), xi_E_current.end());
//sort(xi_EnI_current.begin(), xi_EnI_current.end());
//sort(xi_EnIS_current.begin(), xi_EnIS_current.end());

////---- intialization of xi_U  ----////

xi_U_current.clear();

for (int i=0; i<= (int)(para_other.n-1);i++){
	if(find(xi_E_current.begin(),xi_E_current.end(),i)==xi_E_current.end()){
	xi_U_current.push_back(i);
	t_e_current.at(i)=para_other.unassigned_time;
 	infected_source_current.at(i) = -99;
	}
}

//sort(xi_U_current.begin(), xi_U_current.end());

//---------------------------------------------///


myfile_out.open((string(path4)+string("xi_E_initial.txt")).c_str(),ios::out);
myfile_out << "k"  << endl;
if (xi_E_current.empty()!=1){
for (int i=0; i<=((int)xi_E_current.size()-1);i++){
myfile_out << xi_E_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_E_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_I_initial.txt")).c_str(),ios::out);
myfile_out << "k"  << endl;
if (xi_I_current.empty()!=1){
for (int i=0; i<=((int)xi_I_current.size()-1);i++){
myfile_out << xi_I_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_I_current.size();
myfile_out.close();

myfile_out.open((string(path4)+string("xi_EnI_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_EnI_current.empty()!=1){
for (int i=0; i<=((int)xi_EnI_current.size()-1);i++){
myfile_out << xi_EnI_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_EnI_current.size();
myfile_out.close();



myfile_out.open((string(path4)+string("xi_U_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_U_current.empty()!=1){
for (int i=0; i<=((int)xi_U_current.size()-1);i++){
myfile_out << xi_U_current.at(i)  << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_U_current.size();
myfile_out.close();


myfile_out.open((string(path4)+string("xi_E_minus_initial.txt")).c_str(),ios::out);
myfile_out << "k" << endl;
if (xi_E_minus_current.empty()!=1){
for (int i=0; i<=((int)xi_E_minus_current.size()-1);i++){
myfile_out << xi_E_minus_current.at(i) << endl;
}
}
myfile_out << "size" << endl;
myfile_out << xi_E_minus_current.size();
myfile_out.close();

//---------------------------//*/

////---- intialization of infected_source, infecting_list, infecting_size & t_e & t_nt for infected ----////

/*
	//--- when RANDOMLY exclude a part of sampled sequences --//
	
	double p_sample = 0.05;// pr a sample included for an exposure
	
	for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections
	
		double PS[2] = {1.0 - p_sample, p_sample};
		gsl_ran_discrete_t * gs = gsl_ran_discrete_preproc (sizeof(PS)/sizeof(PS[0]),PS);
		int sample_ind = gsl_ran_discrete (r, gs); // 1 = would include the sample

		gsl_ran_discrete_free(gs);
	
		switch(sample_ind){
			case 0:{// exclude

			if ((find( index_current.begin(), index_current.end(), xi_E_current.at(i) )==index_current.end())&(nt_data_current.t_sample.at(xi_E_current.at(i)) !=para_other.unassigned_time)){ // if it is an non-index & it has sample
		
				if((find(xi_EnI_current.begin(),xi_EnI_current.end(),xi_E_current.at(i))!=xi_EnI_current.end())&(nt_data_current.t_sample.at(xi_E_current.at(i)) !=para_other.unassigned_time)){// this infection is in xi_EnI& had a sample(which the sample gotta be deleted)
				xi_EnIS_current.push_back(xi_E_current.at(i));
				}

				nt_data_current.t_sample.at(xi_E_current.at(i)) =para_other.unassigned_time; // exclude the sample

				nt_data_current.current_size.at(xi_E_current.at(i)) = nt_data_current.current_size.at(xi_E_current.at(i)) - 1;
			}

			break;
			}
			case 1:{ // keep the original t_sample (could be unassigned_time)
				nt_data_current.t_sample.at(xi_E_current.at(i)) = nt_data_current.t_sample.at(xi_E_current.at(i));
			break;
			}
		}
	}

	myfile_out.open((string(path4)+string("p_sample.txt")).c_str(),ios::out);
	myfile_out << p_sample;
	myfile_out.close();

	myfile_out.open((string(path4)+string("t_sample_initial.txt")).c_str(),ios::app);
	for (int i=0; i<=(para_other.n-1);i++){
	myfile_out << nt_data_current.t_sample.at(i) << endl;
	}
	myfile_out.close();

	int total_sample =0;
	for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections
		if (nt_data_current.t_sample.at(xi_E_current.at(i))!=para_other.unassigned_time) total_sample = total_sample+1;
	}
	double p_sample_actual = ((double)  total_sample)/ ((double) xi_E_current.size() );
	myfile_out.open((string(path4)+string("p_sample_actual.txt")).c_str(),ios::out);
	myfile_out << p_sample_actual;
	myfile_out.close();

	myfile_out.open((string(path4)+string("xi_EnIS_initial.txt")).c_str(),ios::out);
	myfile_out << "k" << endl;
	if (xi_EnIS_current.empty()!=1){
	for (int i=0; i<=((int)xi_EnIS_current.size()-1);i++){
	myfile_out << xi_EnIS_current.at(i) << endl;
	}
	}
	myfile_out << "size" << endl;
	myfile_out << xi_EnIS_current.size();
	myfile_out.close();


	//---------------------------------------------------------------------------------------//

*/


	//--- when exclude INITIAL part of sampled sequences according to order of infectious time  --//
	
// 	double p_sample = 0.56;// pr a sample included for an exposure
// 	int n_sample = (int) trunc( xi_I_current.size()*p_sample);
// 
// 	vector<double> t_i_current_sort = t_i_current; 
// 	sort(t_i_current_sort.begin(), t_i_current_sort.end());
// 
// 	vector<int> xi_I_current_sort ((int)xi_I_current.size());
// 
// 	for (int i=0; i<= (int)(xi_I_current.size()-1);i++){
// 		int rank = distance(t_i_current_sort.begin(), find(t_i_current_sort.begin(), t_i_current_sort.end(), t_i_current.at(xi_I_current.at(i))) );
// 
// 		xi_I_current_sort.at(rank) = xi_I_current.at(i); // needa make sure no same t_i
// 	}
// 	
// 	for (int i=(n_sample); i<= (int)( xi_I_current_sort.size()-1);i++){// loop over infections
// 		nt_data_current.t_sample.at(xi_I_current_sort.at(i)) =para_other.unassigned_time; // exclude the sample
// 	}
// 
// 	int total_sample =0;
// 	for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections
// 		if (nt_data_current.t_sample.at(xi_E_current.at(i))!=para_other.unassigned_time) total_sample = total_sample+1;
// 	}
// 	double p_sample_actual = ((double)  total_sample)/ ((double) xi_E_current.size() );
// 	myfile_out.open((string(path4)+string("p_sample_actual.txt")).c_str(),ios::out);
// 	myfile_out << p_sample_actual;
// 	myfile_out.close();
// 
// 	myfile_out.open((string(path4)+string("xi_EnIS_initial.txt")).c_str(),ios::out);
// 	myfile_out << "k" << endl;
// 	if (xi_EnIS_current.empty()!=1){
// 	for (int i=0; i<=((int)xi_EnIS_current.size()-1);i++){
// 	myfile_out << xi_EnIS_current.at(i) << endl;
// 	}
// 	}
// 	myfile_out << "size" << endl;
// 	myfile_out << xi_EnIS_current.size();
// 	myfile_out.close();

	//---------------------------------------------------------------------------------------------------------------------//
/*
	//-----  import source from the posterior sample using only epidemic data----//

	vector<int> infected_source_start(para_other.n);

	ifstream myfile_in;
	string line, field;

	myfile_in.open((string(path4)+string("infected_source_start.txt")).c_str(),ios::in);
	int line_count=0;
	
	while (getline(myfile_in,line)) {
	
	stringstream ss(line);
	int field_count=0;
	
	while (getline(ss, field, ',' )) {
	stringstream fs (field);

// 	if ((line_count>=1) & (field_count==0)) fs >> epi_final_arg.k.at(line_count-1);
// 	if ((line_count>=1) & (field_count==1)) fs >> epi_final_arg.q.at(line_count-1);
// 	if ((line_count>=1) & (field_count==2)) fs >> epi_final_arg.coor_x.at(line_count-1); //note: k_1 and k_2 were named as par_kernel_1 and par_kernel_2 in simulation
// 	if ((line_count>=1) & (field_count==3)) fs >> epi_final_arg.coor_y.at(line_count-1);
// 	if ((line_count>=1) & (field_count==4)) fs >> epi_final_arg.t_e.at(line_count-1);
// 	if ((line_count>=1) & (field_count==5)) fs >> epi_final_arg.t_i.at(line_count-1);
// 	if ((line_count>=1) & (field_count==6)) fs >> epi_final_arg.t_r.at(line_count-1);
// 	if ((line_count>=1) & (field_count==7)) fs >> epi_final_arg.status.at(line_count-1);
// 	if ((line_count>=1) & (field_count==8)) fs >> epi_final_arg.gp_stb.at(line_count-1);
// 	if ((line_count>=1) & (field_count==9)) fs >> epi_final_arg.stb.at(line_count-1);
// 	if ((line_count>=1) & (field_count==10)) fs >> epi_final_arg.infected_source.at(line_count-1);
	
	fs >> infected_source_start.at(field_count);
	
	field_count = field_count + 1;
	}
	
	line_count = line_count + 1 ;
	}

	myfile_in.close();

	//---------------------------------------------------------------------------------------//*/

for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections

	switch(nt_data_current.t_sample.at(xi_E_current.at(i))!=para_other.unassigned_time){

		case 1:{ // with sample, must be infected

			nt_data_current.current_size.at(xi_E_current.at(i)) = 2;//  (initial size) infection + the sample

			nt_data_current.t_nt.at(xi_E_current.at(i)).push_back(nt_data_current.t_sample.at(xi_E_current.at(i)));
			//nt_data_current.nt.at(xi_E_current.at(i)).insert(nt_data_current.nt.at(xi_E_current.at(i)).begin(), sample_data.at(xi_E_current.at(i)).begin(), sample_data.at(xi_E_current.at(i)).end());
		break;
		}

		case 0:{
			nt_data_current.current_size.at(xi_E_current.at(i)) = 1;//  (initial size) infection
		break;
		}

	}
}

//--

for (int i=0; i<= (int)(xi_E_current.size()-1);i++){// loop over infections
	
	int source;

	//--- when propose t_e before deciding source_pool---//
// 	double t_up =   min( nt_data_current.t_sample.at(xi_E_current.at(i)), min( t_i_current.at(xi_E_current.at(i)),para_other.t_max));
// 	double t_low = max(0.0, t_up - 2.0 );
// 
// 	//t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,t_low, t_up);
// 	t_e_current.at(xi_E_current.at(i)) =  max(0.0, t_up - gsl_ran_gamma(r, 10.0, 0.5) );
// 
// 	//t_e_current.at(xi_E_current.at(i)) =t_e_current.at(xi_E_current.at(i));
	//----------------------------------------------------------//

	vector<int> source_pool; 

	for (int j=0;j<=(int)(xi_I_current.size()-1);j++){

		//double t_low = min(nt_data_current.t_sample.at(xi_E_current.at(i)),t_i_current.at(xi_E_current.at(i)));
		double t_bound = min(nt_data_current.t_sample.at(xi_E_current.at(i)), min(t_i_current.at(xi_E_current.at(i)), para_other.t_max)); //
		
		switch( (t_i_current.at(xi_I_current.at(j))<t_bound) & (t_i_current.at(xi_I_current.at(j))>=(t_bound-10.0)) ){ // when propose t_e after deciding source_pool; 2nd condition only allows a subset of xi_I to be possible sources (avoid f_EnI=0)
//		switch( (t_i_current.at(xi_I_current.at(j))<t_e_current.at(xi_E_current.at(i))) &  (t_r_current.at(xi_I_current.at(j))>=t_e_current.at(xi_E_current.at(i)))  ){ //when propose t_e before deciding source_pool


			case 1:{
				source_pool.push_back(xi_I_current.at(j));	
			break;
			}
			case 0:{
			break;
			}
		}
	}
	
	source_pool.insert(source_pool.begin(),9999);
	
	int num_infectious = (int)source_pool.size();

	//-----------propose uniformly ----------------//

	source = source_pool.at(gsl_rng_uniform_int(r, num_infectious)); // uniformly choose a new source (including bg)

	
	//-----propose according to infectious challenge--------------------//
	
// 	vector<double> ic(num_infectious);
// 	ic.at(0) = para_current.alpha;
// 	
// 	switch(num_infectious>=2){
// 	
// 		case 1:{ // with 2nd sources from pool
// 	
// 			for (int j=1;j<=(num_infectious-1);j++){
// 			ic.at(j)= para_current.beta*stb_current.at(xi_E_current.at(i))*kernel_mat_current[xi_E_current.at(i)][source_pool.at(j)]/norm_const_current.at(source_pool.at(j)); // a new source will be proposed according to the infectious challenges
// 			}
// 	
// 			double *P=&ic.at(0); // convert vector to array
// 			gsl_ran_discrete_t * g = gsl_ran_discrete_preproc ((int)ic.size(),P);
// 			int link= gsl_ran_discrete (r, g);
// 			gsl_ran_discrete_free (g);
// 	
// 			source = source_pool.at(link); // a new source
// 			
// 			break;
// 		}
// 	
// 		case 0:{ // only primary source from pool
// 	
// 			source = 9999;
// 	
// 		break;
// 		}
// 	}
// 	
// 	myfile_out.open((string(path4)+string("initial_ic.txt")).c_str(),ios::app);
// 	for (int j=0; j<= (int) (ic.size()-1);j++){
// 	myfile_out <<ic.at(j) << endl;
// 	}
// 	myfile_out.close();

// 	//----set source from posterior sample  with only epidemic data---//
// 	if( find(source_pool.begin(), source_pool.end(), infected_source_start.at(xi_E_current.at(i)))!=source_pool.end()){
// 	source = infected_source_start.at(xi_E_current.at(i)); 
// 	}
// 
// 	//-------------------------------------------------------------------------------------------------//

	if ( find( index_current.begin(), index_current.end(), xi_E_current.at(i) )!=index_current.end()) source = 9999; // index is known, source = 9999

	//-------------------------------------------------------------------------------------------------//

	infected_source_current.at(xi_E_current.at(i)) = source;

	switch(source==9999){
		case 0:{// 2nd infection

			//--------when propose t_e after deciding source_pool//

			double t_up = min(min(nt_data_current.t_sample.at(xi_E_current.at(i)), t_r_current.at(source)), min(t_i_current.at(xi_E_current.at(i)), para_other.t_max));			
			double t_low = max(t_i_current.at(source), t_up -10.0 );
			t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,t_low, t_up);

			//-------------------------------------------------------------------------//

			nt_data_current.t_nt.at(xi_E_current.at(i)).push_back(t_e_current.at(xi_E_current.at(i)));
			sort(nt_data_current.t_nt.at(xi_E_current.at(i)).begin(), nt_data_current.t_nt.at(xi_E_current.at(i)).end() );

			nt_data_current.t_nt.at(source).push_back(t_e_current.at(xi_E_current.at(i)));
			sort(nt_data_current.t_nt.at(source).begin(), nt_data_current.t_nt.at(source).end() );


			switch(nt_data_current.infecting_size.at(source)>=1){
				case 1:{
					vector<double> t_y(nt_data_current.infecting_size.at(source));
					for (int k=0;k<=(nt_data_current.infecting_size.at(source)-1);k++){
					t_y.at(k) = t_e_current.at(nt_data_current.infecting_list[source][k]);
					}
					t_y.push_back(t_e_current.at(xi_E_current.at(i)));
					sort(t_y.begin(), t_y.end());
				
					int rank_source = distance(t_y.begin(), find(t_y.begin(), t_y.end(),t_e_current.at(xi_E_current.at(i))));
		
					nt_data_current.infecting_list.at(source).insert( nt_data_current.infecting_list.at(source).begin() + rank_source, xi_E_current.at(i));

				break;
				}
		
				case 0:{
					nt_data_current.infecting_list.at(source).push_back(xi_E_current.at(i));
				break;
				}
			}


			nt_data_current.infecting_size.at(source) = nt_data_current.infecting_size.at(source) + 1;
			nt_data_current.current_size.at(source) = nt_data_current.current_size.at(source) + 1;

		break;
		}
		case 1:{// bg infection

			//--------when propose t_e after deciding source_pool//
// 			double t_up = min(nt_data_current.t_sample.at(xi_E_current.at(i)), min(t_i_current.at(xi_E_current.at(i)), para_other.t_max));			
// 			double t_low = max(0.0, t_up -10.0 );
// 			t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,t_low, t_up);
// 			
			switch( find( index_current.begin(), index_current.end(), xi_E_current.at(i) )!=index_current.end()){// return 1 when it is an index
				case 0:{
					double t_up = min(nt_data_current.t_sample.at(xi_E_current.at(i)), min(t_i_current.at(xi_E_current.at(i)), para_other.t_max));			
					double t_low = max(0.0, t_up -10.0 );
					t_e_current.at(xi_E_current.at(i)) = gsl_ran_flat(r,t_low, t_up);
				break;
				}
				case 1:{
					t_e_current.at(xi_E_current.at(i)) = 0.0;
				break;
				}
			 }
			//-------------------------------------------------------------------------//

			nt_data_current.t_nt.at(xi_E_current.at(i)).push_back(t_e_current.at(xi_E_current.at(i)));
			sort(nt_data_current.t_nt.at(xi_E_current.at(i)).begin(), nt_data_current.t_nt.at(xi_E_current.at(i)).end());

		break;
		}		
	}

}// end of loop over infections


////----initialization of index_current and xi_E_minus --///

index_current.clear();
xi_E_minus_current = xi_E_current;

double min_t = *min_element(t_e_current.begin(),t_e_current.end());
for (int i=0; i<= (int)(xi_E_current.size()-1);i++){
if (t_e_current.at(xi_E_current.at(i))==min_t) {
index_current.push_back(xi_E_current.at(i));
}
}
for (int i=0; i<= (int)(index_current.size()-1);i++){
xi_E_minus_current.erase(find(xi_E_minus_current.begin(),xi_E_minus_current.end(),index_current.at(i)));
}


myfile_out.open((string(path4)+string("initial_index.txt")).c_str(),ios::app);
for (int i=0; i<= (int) (index_current.size()-1);i++){
myfile_out << index_current.at(i) << "," << t_e_current.at(index_current.at(i)) << "," << infected_source_current.at(index_current.at(i)) << endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_current_size.txt")).c_str(),ios::app);
for (int i=0; i<= (int) (para_other.n-1);i++){
myfile_out << nt_data_current.current_size.at(i) << endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_t_e.txt")).c_str(),ios::app);
for (int i=0; i<= (int) (para_other.n-1);i++){
myfile_out << t_e_current.at(i) << endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_t_i.txt")).c_str(),ios::app);
for (int i=0; i<= (int) (para_other.n-1);i++){
myfile_out << t_i_current.at(i) << endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_infected_source.txt")).c_str(),ios::app);
for (int i=0; i<= (int) (para_other.n-1);i++){
myfile_out << infected_source_current.at(i) << endl;
}
myfile_out.close();

myfile_out.open((string(path4)+string("initial_infecting_size.txt")).c_str(),ios::app);
int sum_infected=0;
for (int i=0; i<= (int) (para_other.n-1);i++){
sum_infected = sum_infected +  nt_data_current.infecting_size.at(i) ;
myfile_out << nt_data_current.infecting_size.at(i) << endl;
}
myfile_out <<endl;
myfile_out <<"totoal" <<endl;
myfile_out << sum_infected;
myfile_out.close();


////---- intialization of nt_data_current.nt ////

vector<double> t_e_sort; // sorted t_e excluding susceptible

for (int i=0; i<=(int) xi_E_current.size() -1 ; i++){
	if (t_e_current.at(xi_E_current.at(i))!=para_other.unassigned_time)  t_e_sort.push_back(t_e_current.at(xi_E_current.at(i)));
}

sort( t_e_sort.begin(),  t_e_sort.end()); 

vector<int> xi_E_sort((int) xi_E_current.size());

for (int i=0; i<=(int) xi_E_current.size() -1 ; i++){

	int rank_t = distance(t_e_sort.begin(), find(t_e_sort.begin(),t_e_sort.end(), t_e_current.at(xi_E_current.at(i))) );

	xi_E_sort.at(rank_t) = xi_E_current.at(i);
}


	myfile_out.open((string(path4)+string("xi_E_sort_initial.txt")).c_str(),ios::out);
	for (int i=0; i<=((int)xi_E_sort.size()-1);i++){
	myfile_out << xi_E_sort.at(i) << endl;
	}
	myfile_out.close();


	myfile_out.open((string(path4)+string("t_e_sort_initial.txt")).c_str(),ios::out);
	for (int i=0; i<=((int)xi_E_sort.size()-1);i++){
	myfile_out << t_e_current.at(xi_E_sort.at(i)) << endl;
	}
	myfile_out.close();


vector<int> ind_sample (para_other.n, 0); // indicate if the sample a sampled case  has been included


for (int i=0; i<=(int) xi_E_sort.size() -1 ; i++){

	int subject = xi_E_sort.at(i);
	int source = infected_source_current.at(subject);

// 	nt_data_current.nt.at(subject).erase(nt_data_current.nt.at(subject).begin(), nt_data_current.nt.at(subject).end()); // erase the whole sequences record

	nt_data_current.nt.at(subject).clear(); // erase the whole sequences record

// 	nt_data_current.nt.at(subject).resize(nt_data_current.current_size.at(subject) *para_other.n_base); 


	switch(source==9999) {
		case 0:{// secondary
			switch(nt_data_current.t_sample.at(source)==para_other.unassigned_time){
				case 1:{// source has no sample
					seq_initialize_pair(nt_data_current, source, subject, t_e_current.at(subject), para_other.n_base, para_current);


				break;
				}
				case 0:{// source has sample
					switch((nt_data_current.t_sample.at(source)<= t_e_current.at(subject)) & (ind_sample.at(source)==0)){
						case 1:{ // has to insert sample first
							nt_data_current.nt.at(source).insert(nt_data_current.nt.at(source).end(),sample_data.at(source).begin(), sample_data.at(source).end());
							ind_sample.at(source) = 1;
							seq_initialize_pair(nt_data_current, source, subject, t_e_current.at(subject), para_other.n_base, para_current);	


			
						break;
						}
						case 0:{
							seq_initialize_pair(nt_data_current, source, subject, t_e_current.at(subject), para_other.n_base, para_current);

					break;
						}
					}
				break;
				}
			}
		break;
		}
		case 1:{// bg

			vector<int> seq_new(para_other.n_base);
			sample_snull(con_seq, seq_new, para_current.p_ber, para_other.n_base, r);
			nt_data_current.nt.at(subject).insert(nt_data_current.nt.at(subject).begin(),seq_new.begin(), seq_new.begin()+para_other.n_base);


		break;
		}
	}
}


for (int i=0; i<=(int) xi_E_sort.size() -1 ; i++){

	int subject = xi_E_sort.at(i);

	switch((nt_data_current.t_sample.at(subject)!=para_other.unassigned_time) &(ind_sample.at(subject)==0)){
		case 1:{ // has sample but not yet inclued
			nt_data_current.nt.at(subject).insert(nt_data_current.nt.at(subject).end(),sample_data.at(subject).begin(), sample_data.at(subject).end());
			ind_sample.at(subject)=1;	


		break;
		}
		case 0:{
		break;
		}

	}

}


//---------------------end of  intialization of nt_data_current.nt-------------------------//


}


//---------------------------------//

void seq_initialize_pair(nt_struct& nt_data_arg, int  k_source_arg, int k_arg, double t_now_arg, int n_base, para_key& para_current){ // update sequences of infectious-infected pair

// int rank_source= nt_data_arg.current_size.at(k_source_arg);
// double t_null = nt_data_arg.t_nt[k_source_arg][rank_source-1];

vector<int> seq_new(n_base);

int rank_source = distance(nt_data_arg.t_nt.at(k_source_arg).begin(), find(nt_data_arg.t_nt.at(k_source_arg).begin(), nt_data_arg.t_nt.at(k_source_arg).end(), t_now_arg));
double t_null = nt_data_arg.t_nt[k_source_arg][rank_source-1];

double dt = t_now_arg - t_null;

double p_1 = 0.25 + 0.25*exp(-4.0*para_current.mu_2*dt) + 0.5*exp(-2.0*(para_current.mu_1+para_current.mu_2)*dt) ; // pr of a base not changing
double p_2 = 0.25 + 0.25*exp(-4.0*para_current.mu_2*dt) - 0.5*exp(-2.0*(para_current.mu_1+para_current.mu_2)*dt); // pr of a transition of a base
double p_3 = 1.0*(0.25 - 0.25*exp(-4.0*para_current.mu_2*dt));  // pr of a transversion (two possible events)
double p_4 = p_3;
//double p_1 = 1.0- p_2 - 2.0*p_3;

double P[4] = {p_1, p_2, p_3,p_4};

const gsl_rng_type* T_c= gsl_rng_ranlux;  // T is pointer points to the type of generator
gsl_rng *r_c = gsl_rng_alloc (T_c); // r is pointer points to an object with Type T
gsl_ran_discrete_t * g = gsl_ran_discrete_preproc (sizeof(P)/sizeof(P[0]),P);

int count_1, count_2, count_3, count_4;

count_1 = count_2 = count_3 =count_4 =0;

gsl_rng_set (r_c,-1000*k_arg*k_source_arg*(rank_source+1)); // set a seed


for (int j=0;j<=(n_base-1);j++){

	int type= gsl_ran_discrete (r_c, g) + 1; 


	switch(nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j]){

	case 1:{ // an A
		switch(type){
		case 1:{
//		nt_data_arg.nt[k_source_arg].push_back(nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j]);
		seq_new.at(j) = nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j];
		break;
		}
		case 2:{
// 		nt_data_arg.nt[k_source_arg].push_back(2);
		seq_new.at(j) = 2;
		break;
		}
		case 3:{
// 		nt_data_arg.nt[k_source_arg].push_back(3);
		seq_new.at(j) = 3;
		break;
		}					
		case 4:{
// 		nt_data_arg.nt[k_source_arg].push_back(4);
		seq_new.at(j) = 4;
		break;
		}	
		}
	break;
	}
	
	case 2:{ // a G
		switch(type){
		case 1:{
// 		nt_data_arg.nt[k_source_arg].push_back(nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j]);
		seq_new.at(j) = nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j];
		break;
		}
		case 2:{
// 		nt_data_arg.nt[k_source_arg].push_back(1);
		seq_new.at(j) = 1;
		break;
		}		
		case 3:{
// 		nt_data_arg.nt[k_source_arg].push_back(3);
		seq_new.at(j) = 3;
		break;
		}			
		case 4:{
// 		nt_data_arg.nt[k_source_arg].push_back(4);
		seq_new.at(j) = 4;
		break;
		}	
		}
	break;
	}
	
	case 3:{ // a T
		switch(type){
		case 1:{
// 		nt_data_arg.nt[k_source_arg].push_back(nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j]);
		seq_new.at(j) = nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j];
		break;
		}
		case 2:{
// 		nt_data_arg.nt[k_source_arg].push_back(4);
		seq_new.at(j) = 4;
		break;
		}		
		case 3:{
// 		nt_data_arg.nt[k_source_arg].push_back(1);
		seq_new.at(j) = 1;
		break;
		}		
		case 4:{
// 		nt_data_arg.nt[k_source_arg].push_back(2);
		seq_new.at(j) = 2;
		break;
		}	
		}
	break;
	}

	case 4:{ // a C
		switch(type){
		case 1:{
// 		nt_data_arg.nt[k_source_arg].push_back(nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j]);
		seq_new.at(j) = nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+j];
		break;
		}
		case 2:{
// 		nt_data_arg.nt[k_source_arg].push_back(3);
		seq_new.at(j) = 3;
		break;
		}		
		case 3:{
// 		nt_data_arg.nt[k_source_arg].push_back(1); 
		seq_new.at(j) = 1;
		break;
		}		
		case 4:{
// 		nt_data_arg.nt[k_source_arg].push_back(2);
		seq_new.at(j) = 2;
		break;
		}	
		}
	break;
	}

	}
}

gsl_ran_discrete_free (g);

nt_data_arg.nt.at(k_source_arg).insert(nt_data_arg.nt.at(k_source_arg).end(), seq_new.begin(), seq_new.end());
nt_data_arg.nt.at(k_arg).insert(nt_data_arg.nt.at(k_arg).end(), seq_new.begin(), seq_new.end());

// ofstream myfile; // this object can be recycled after myfile.close()
// myfile.open((string(path4)+string("00.txt")).c_str(),ios::app);
// myfile <<  k_source_arg << "," << k_arg <<"," << nt_data_arg.nt[k_source_arg][(rank_source-1)*n_base+0] << "," << nt_data_arg.nt[k_source_arg][(rank_source)*n_base] << "," <<  nt_data_arg.nt[k_arg][0] <<"," << seq_new.at(0) << endl;
// myfile.close();


}



