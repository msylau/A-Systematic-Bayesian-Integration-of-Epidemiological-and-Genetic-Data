
//This source file is to perform some tedious and lengthy operations such as read in simulated data //

#include "phylo_seir_fmd_inference_header.h"


void IO_simpara (para_key& para_true_arg, para_aux& para_other_arg){ // this function is called to input/output parameters from simulation


ifstream myfile_in_simpara;
ofstream myfile_out_simpara; 

string line, field;
int line_count=0, field_count=0;


myfile_in_simpara.open((string(path2)+string("parameters_other.txt")).c_str(),ios::in);
line_count =0;

while (getline(myfile_in_simpara,line)) {

//getline(myfile_in_simdata,line);
stringstream ss(line);
//string field;
field_count=0;

while (getline(ss, field, ',' )) {
stringstream fs (field);
if ((line_count==1) & (field_count==0)) fs >> para_other_arg.n;
if ((line_count==1) & (field_count==1)) fs >> para_other_arg.kernel_type;
if ((line_count==1) & (field_count==2)) fs >> para_other_arg.dimen_x;
if ((line_count==1) & (field_count==3)) fs >> para_other_arg.dimen_y;
if ((line_count==1) & (field_count==4)) fs >> para_other_arg.t_max;
if ((line_count==1) & (field_count==5)) fs >> para_other_arg.unassigned_time;
if ((line_count==1) & (field_count==6)) fs >> para_other_arg.seed;
if ((line_count==1) & (field_count==7)) fs >> para_other_arg.n_seq;

//if ((line_count==1) & (field_count==8)) fs >> para_other_arg.n_base;

if ((line_count==1) & (field_count==8)){
	switch(ind_n_base_part==1){
		case 0:{
		fs >> para_other_arg.n_base;
		break;
		}
		case 1:{
		para_other_arg.n_base = n_base_part; // this change the para_other.n_base to be n_base_part
		break;
		}
	}
}

field_count = field_count + 1;
}

line_count =line_count + 1;

}

myfile_in_simpara.close();



myfile_out_simpara.open((string(path4)+string("parameters_other.txt")).c_str(),ios::out);
myfile_out_simpara << "n" << "," << "kernel_type" << "," << "dimen_x" << "," << "dimen_y" << "," << "t_max" << "," << "unassigned_time" << "," << "seed" << "," << "n_seq" << "," << "n_base" << endl;
myfile_out_simpara << para_other_arg.n << "," << para_other_arg.kernel_type << "," << para_other_arg.dimen_x << "," << para_other_arg.dimen_y << "," << para_other_arg.t_max << "," << para_other_arg.unassigned_time << "," << para_other_arg.seed << "," << para_other_arg.n_seq<< "," << para_other_arg.n_base << endl;
myfile_out_simpara.close();


}


/*--------------------------------------------*/
//-------------
void IO_simdata (para_key para_true_arg, para_aux para_other_arg, vector < vector<double> >& coordinate_arg, epi_struct& epi_final_arg, nt_struct& nt_data_arg, vector<int>& index_arg, vector<int>& con_seq){ // this function is called to input/output data from simulation

ifstream myfile_in_simdata;
ofstream myfile_out_simdata; 

string line, field;
int line_count=0, field_count=0;


myfile_in_simdata.open((string(path2)+string("coordinate.txt")).c_str(),ios::in);
//string line;
line_count=0;

//coordinate_arg.resize(para_other_arg.n);

while (getline(myfile_in_simdata,line)) {

 //getline(myfile_in_simdata,line);
 stringstream ss(line);
 //string field;
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);
 fs >> coordinate_arg[line_count][field_count];
 field_count = field_count + 1;
 }

line_count = line_count + 1;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("coordinate.txt")).c_str(),ios::out);
for (int i=0;i<=(para_other_arg.n-1);i++){
myfile_out_simdata << coordinate_arg[i][0] << "," << coordinate_arg[i][1]<< endl;
}
myfile_out_simdata.close();

/*--------------------------------------------*/


epi_final_arg.k.resize(para_other_arg.n);
epi_final_arg.q.resize(para_other_arg.n);
epi_final_arg.t_e.resize(para_other_arg.n);
epi_final_arg.t_i.resize(para_other_arg.n);
epi_final_arg.t_r.resize(para_other_arg.n);
epi_final_arg.status.resize(para_other_arg.n);
epi_final_arg.coor_x.resize(para_other_arg.n);
epi_final_arg.coor_y.resize(para_other_arg.n);

epi_final_arg.gp_stb.resize(para_other_arg.n);
epi_final_arg.stb.resize(para_other_arg.n);

epi_final_arg.infected_source.resize(para_other_arg.n);

myfile_in_simdata.open((string(path2)+string("epi_fmd.txt")).c_str(),ios::in);
line_count=0;

while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);
 field_count=0;

 while (getline(ss, field, ',' )) {
 stringstream fs (field);
 if ((line_count>=1) & (field_count==0)) fs >> epi_final_arg.k.at(line_count-1);
 if ((line_count>=1) & (field_count==1)) fs >> epi_final_arg.coor_x.at(line_count-1); //note: k_1 and k_2 were named as par_kernel_1 and par_kernel_2 in simulation
 if ((line_count>=1) & (field_count==2)) fs >> epi_final_arg.coor_y.at(line_count-1);
 if ((line_count>=1) & (field_count==3)) fs >> epi_final_arg.t_e.at(line_count-1);
 if ((line_count>=1) & (field_count==4)) fs >> epi_final_arg.t_i.at(line_count-1);
 if ((line_count>=1) & (field_count==5)) fs >> epi_final_arg.t_r.at(line_count-1);
 if ((line_count>=1) & (field_count==6)) fs >> epi_final_arg.stb.at(line_count-1);
 if ((line_count>=1) & (field_count==7)) fs >> epi_final_arg.gp_stb.at(line_count-1);

if ((line_count>=1) & (field_count==8)) fs >> epi_final_arg.infected_source.at(line_count-1);


 field_count = field_count + 1;
 }

line_count = line_count + 1 ;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("epi_final.txt")).c_str(),ios::app);
myfile_out_simdata << "k" << "," << "coor_x" << "," << "coor_y" << "," << "t_e" << "," << "t_i"<< "," << "t_r" <<"," << "stb" << "," << "gp_stb"<<  "," <<  "infected_source"<< endl;
for (int i=0; i<=(para_other_arg.n-1);i++){
myfile_out_simdata << epi_final_arg.k.at(i) << "," << epi_final_arg.coor_x.at(i) << "," << epi_final_arg.coor_y.at(i)<< "," << epi_final_arg.t_e.at(i) << "," << epi_final_arg.t_i.at(i)<< "," << epi_final_arg.t_r.at(i)  <<  "," << epi_final_arg.stb.at(i) << "," << epi_final_arg.gp_stb.at(i) << ","  << epi_final_arg.infected_source.at(i) <<endl;
}
myfile_out_simdata.close();




/*--------------------------------------------*/

index_arg.reserve(para_other_arg.n);

myfile_in_simdata.open((string(path2)+string("index.txt")).c_str(),ios::in);
line_count=0;


while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);

 while (getline(ss, field)) {
 stringstream fs (field);
 if (line_count>=1) {
 int ind;
 fs >> ind;
 index_arg.push_back(ind); 

 }
 }

line_count = line_count + 1 ;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("index.txt")).c_str(),ios::app);
myfile_out_simdata << "k" << endl;
for (int i=0; i<=((int)index_arg.size()-1);i++){
myfile_out_simdata << index_arg.at(i) << endl;
}
myfile_out_simdata.close();


/*--------------------------------------------*/



nt_data_arg.nt.resize(para_other_arg.n);


for (int i=0; i<=(para_other_arg.n-1); i++){
nt_data_arg.nt[i].reserve(para_other_arg.n_seq*para_other_arg.n_base); 
}

nt_data_arg.t_nt.resize(para_other_arg.n);
for (int i=0; i<=(para_other_arg.n-1); i++){
nt_data_arg.t_nt[i].reserve(para_other_arg.n_seq); 
}

nt_data_arg.current_size.resize(para_other_arg.n);
nt_data_arg.t_sample.resize(para_other_arg.n); 



myfile_in_simdata.open((string(path2)+string("t_sample.txt")).c_str(),ios::in);

line_count=0;

while (getline(myfile_in_simdata,line)) {

 stringstream ss(line);

 while (getline(ss, field)) {
 stringstream fs (field);
double t;
 fs >>t;
 nt_data_arg.t_sample.at(line_count) = t; 
 }
 
line_count = line_count + 1 ;
}

myfile_in_simdata.close();


myfile_out_simdata.open((string(path4)+string("t_sample.txt")).c_str(),ios::app);
for (int i=0; i<=(para_other_arg.n-1);i++){
myfile_out_simdata << nt_data_arg.t_sample.at(i) << endl;
}
myfile_out_simdata.close();

/*--------------------------------------------*/
	
myfile_in_simdata.open((string(path2)+ string("con_seq_estm.txt")).c_str(),ios::in);

line_count = 0;

while (getline(myfile_in_simdata,line)) {

stringstream ss(line);
field_count=0;

while (getline(ss, field, ',')) {
stringstream fs (field);
int nt;
fs >> nt;
con_seq.push_back(nt);

field_count = field_count + 1;
}

line_count = line_count + 1 ;
}

myfile_in_simdata.close();

myfile_out_simdata.open((string(path4)+ string("con_seq_estm.txt")).c_str(),ios::app);
for (int j=0;j<=(para_other_arg.n_base-1);j++){
int rem = (j+1)%para_other_arg.n_base;
if ((rem!=0) | (j==0)) myfile_out_simdata << con_seq.at(j) << ",";
if ((rem==0) & (j!=0)) myfile_out_simdata <<  con_seq.at(j)  << " " << endl;
}
myfile_out_simdata.close();

}
