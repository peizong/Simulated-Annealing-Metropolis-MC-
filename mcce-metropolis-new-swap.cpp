//#include <mpi.h>
#include <string>
#include <iostream>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include "parse.h"
#include "clus_str.h"
#include "getvalue.h"
#include "ctype.h"
#include "version.h"
#include "plugin.h"
#include <random>
#include <vector>
#include <omp.h>
using namespace std;

//define constants
const int binNum=30; //66; //int(argv[2])-int(argv[1]); //46;
const int NumEle=2; //3;
//define new functions
void show_neighbors(vector<vector<vector<int>>> neighbor_list);
int genRandom(int range_from, int range_to);
void readLog(double gamma,array<double,binNum> dos, array<int,binNum> hist, array<int,binNum> count, double sro[binNum][NumEle], const char *logfile);
void writeLog(double gamma,array<double,binNum> dos, array<int,binNum>hist, array<int,binNum> count, double sro[binNum][NumEle], const char *logfile);
float maxArray(array<float,binNum>rest, float maskedValue);
float etot(const char* cmd);
//void swap_write_struct(int atom_a,int atom_b, Structure &str, Array<AutoString> &atom_label, rMatrix3d &axes, char *file);
Structure swapped_str(int atom_a,int atom_b, Structure &str);
float cal_etot(Structure ideal_str,LinkedList<MultiCluster> clusterlist,SpaceGroup spacegroupi, Array<Arrayint> labellookup, Array<Real> eci);
vector<double> short_range_order_bcc(Structure str, double latt, int ith_NN, int NumEle);
//vector<double> short_range_order_fcc(Structure str, double latt, int ith_NN, int NumEle);
vector<vector<double>> short_range_order(Structure str, double latt, int ith_NN, vector<int>component, int cryStru);
//new functions
vector<vector<vector<int>>> find_neighbor_list(Structure str, double latt, int ith_NN, vector<int> component, int cryStru);
vector<vector<vector<int>>> find_neighbor_list_swapped(vector<vector<vector<int>>> neighbor_list, int a, int b);
vector<vector<double>> short_range_order2(vector<vector<vector<int>>> neighbor_list,vector<int> component);
Vector3d<double> NNNbcc(double latt,int ithNN,int jth_ithNN);
Vector3d<double> NNNfcc(double latt,int ithNN,int jth_ithNN);
void moveInBox(Vector3d<double> &Atom, double Box[3]);
int matchVector(Vector3d<double> atom1, Vector3d<double> atom2, double Box[3]);
void hello(int i);
vector<int> random_pick_two_numbers(vector<int> component);
int countSameElement(Arrayint labels,int elementLabel);
//double Temp_scheme(int Ni, int N_final, int NT_sample, double T);
double Temp_scheme(int Ni,int N_final, int NT_sample, double T_init, double T_final,int Scheme);
vector<double> cal_specific_heat(vector<double>temperatures,vector<double> energies, int NT_sample, int size);
//double cal_energy_nn(vector<vector<double>> sro, vector<double>eci_nn, vector<int> component, int cryStru);
double cal_energy_nn(vector<vector<double>> sro, vector <vector<double>> eci_matrix, vector<int> component, int cryStru);
vector<string> split_fun(string str, char delimiter);
vector<double> long_range_order(Structure str, vector<int> component, int cryStru);
vector<double> long_range_order2(Structure str,vector<vector<vector<int>>>neighbor_list, vector<int> component);
vector<vector<double>> cal_LRO(int NT_sample, vector<vector<double>>lros);
vector<double> cal_internal_energy(vector<double>temperatures,vector<double> energies, int NT_sample, int size);
vector<double> cal_autocorrelation(vector<double> energies);
double cal_dE_nn(vector<vector<vector<int>>>neighbor_list, int a, int b, vector<vector<double>> eci_matrix);

int main(int argc, char **argv) {
int N_final,T_init,T_final,NT_sample,Scheme;
int cryStruct, numEle;
//int show_detail_step, cal_short_range_order, restart;
const char *log="log.wlce";
const char *restartLog="restart.log.wlce";
const char * split = ",";
char * p;
vector<char> pp;
pp.resize(2);
char ch[101];
int *ptr_NumEle;
int NNeighbor,test_label;
ptr_NumEle = (int*)( &NumEle );
//-------------read input data-----------------
fstream fin("incar", ios::in);
int row=0;
if (!fin) {cout<<"Sorry bad file"; exit(0);}
while (!fin.eof())
  {
    row=row+1;
    fin.getline(ch,100);
    if (row==2) { p=strtok(ch,split);  cryStruct=atoi(p);}
    if (row==3) { p=strtok(ch,split);  numEle=atoi(p);}//*ptr_NumEle=atoi(p);}
    if (row==4) { p=strtok(ch,split);  N_final=atoi(p);}
    if (row==5) { p=strtok(ch,split);  T_init=atoi(p);}
    if (row==6) { p=strtok(ch,split);  T_final=atoi(p);}
    if (row==7) { p=strtok(ch,split);  NT_sample=atoi(p);}
    if (row==8) { p=strtok(ch,split);  NNeighbor=atoi(p);}
    if (row==9) { p=strtok(ch,split);  Scheme=atoi(p);}
    if (row==10) {p=strtok(ch,split);  test_label=atoi(p);}
cout <<"I finished reading (Row, Value) "<<row<<" "<<ch<<'\n';
//cout<<lb<<" "<<ub<<" "<<gamma_start<<" "<<gamma_end<<" "<<show_detail_step<<'\n';
//      for (int i=0;i<D;i++){
//          p_range[0][i]=atof(p);
//          p = strtok(NULL,split); }}
//     if (row==7) { p=strtok(ch,split);
//      for (int i=0;i<D;i++){
//          p_range[1][i]=atof(p);
//          p = strtok(NULL,split); }}
  }
  fin.close();
vector<double> eci_nn;
eci_nn.resize(numEle*(numEle-1)/2);
if (NNeighbor){
// vector<double> eci_nn;
// eci_nn.resize(numEle*(numEle-1)/2);
 fstream ECINNf("eci_nn",ios::in);
 row=0;
 while (!ECINNf.eof()){
   row +=1;
   ECINNf.getline(ch,100);
   if ((row>1) and (row<eci_nn.size()+2)){ p=strtok(ch,split);  eci_nn[row-2]=atof(p);}
 }
 ECINNf.close();
}
// prepare the eci matrix
vector <vector<double>> eci_matrix;
int eci_id=0;
eci_matrix.resize(numEle,vector<double>(numEle,0));
for (int i=0;i<numEle;i++){                                                                                        
  for (int j=i+1;j<numEle;j++){
    eci_matrix[i][j]=eci_nn[eci_id];
    eci_id +=1;}
  for (int j=0;j<i;j++){
    eci_matrix[i][j]=eci_matrix[j][i];}  
}
/*
for (int i=0;i<numEle;i++){
  for (int j=0;j<numEle;j++){
    cout <<eci_matrix[i][j]<<"	";
  }
  cout << endl;
} */

//---------------------------------------------
//int inNum=atoi(argv[2])-atoi(argv[1])+1;
//cout <<argv[1]<<" "<<argv[2]<<" "<<inNum<<'\n';
vector<int> component;
component.resize(numEle);
float E0=0;
float E1=0;
const char *strfilename="init.struct";
const char *restart_str_file="restart.struct";
/*if (restart)
  {strfilename=restart_str_file;}
else
  {strfilename="init.struct";}
*/
const char *latticefilename="lat.in";
const char *ecifile="eci.out";
//-----------------read lattice parameter---------------------
double latt=1.0;
/*
fstream flatt(latticefilename, ios::in);
int line=0;
if (!flatt) {cout<<"Sorry bad file"; exit(0);}
while (!flatt.eof()){
  line=line+1;
  flatt.getline(ch,100);
  if (line==1) { p=strtok(ch,split);  latt=atof(p);}
  cout <<"lattice para"<<line<<latt<<'\n';
  break;
}
flatt.close(); */
//-----------------read lattice file and get atom labels and space groups-------------------
int readocc=0;
  Structure lattice;                                                             
  Array<Arrayint> labellookup;                                                   
  Array<AutoString> label;                                                       
  rMatrix3d axes;                                                                
  {                                                                              
    ifstream latticefile(latticefilename);                                       
    if (!latticefile) ERRORQUIT("Unable to open lattice file");                  
    if (readocc) {                                                               
      Array<Array<Real> > occ;                                                   
      parse_rndstr_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &occ, &labellookup, &label, latticefile, &axes);
    }                                                                            
    else {                                                                       
      parse_lattice_file(&lattice.cell, &lattice.atom_pos, &lattice.atom_type, &labellookup, &label, latticefile, &axes);                                      
    }                                                                            
    wrap_inside_cell(&lattice.atom_pos,lattice.atom_pos,lattice.cell);           
  } 
  SpaceGroup spacegroup;                                                                 
  spacegroup.cell=lattice.cell;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lattice.cell,lattice.atom_pos,lattice.atom_type);
//--------------------read initial structure files------------------------
  Structure str;
  ifstream strfile(strfilename);
  if (!strfile) ERRORQUIT("Unable to open structure file");
  int strnum=1;
  while (!strfile.eof()) {
    parse_structure_file(&str.cell,&str.atom_pos,&str.atom_type,label,strfile);
    skip_to_next_structure(strfile);
    wrap_inside_cell(&str.atom_pos,str.atom_pos,str.cell);
    rMatrix3d supercell=(!lattice.cell)*str.cell;
    rMatrix3d rounded_supercell=to_real(to_int(supercell));
    rMatrix3d transfo=lattice.cell*rounded_supercell*(!str.cell);
    str.cell=transfo*str.cell;
   // cout<<str.cell(0,0)<<endl;
//    cout<<str.atom_type<<endl;
//    cout<<str.atom_pos<<endl;
   // for (int i=0; i<str.atom_pos.get_size(); i++) {
   //   cout <<str.atom_type<<endl; //str.atom_pos.get_size()<<'\n';
   // }
  }
  if (test_label==1){
    for (int i=0; i<str.atom_pos.get_size(); i++) {
      cout <<str.atom_type[i]<<endl; //str.atom_pos.get_size()<<'\n';
  }}
  for (int i=0;i<numEle;i++){
  component[i]=countSameElement(str.atom_type,i);                                                                                                                           
  }
//--------------------read clusterlist from cluster.out and eci from eci.out-------------------------------
  LinkedList<MultiCluster> clusterlist;
  ifstream clusterfile("clusters.out"); 
  if (!clusterfile) ERRORQUIT("Unable to open clusters.out file");
  read_clusters_and_eci(&clusterlist,NULL,clusterfile,clusterfile,axes);
 // cout<<clusterlist<<endl;
  LinkedListIterator<MultiCluster> icluster(clusterlist);
  Array<Real> eci;
  if (strlen(ecifile)>0) {
    ifstream file(ecifile);
    if (!file) ERRORQUIT("Unable to open ECI file.");
    LinkedList<Real> ecilist;
    while (skip_delim(file)) {
      Real e; 
      file >> e;
      ecilist << new Real(e);
    }         
    LinkedList_to_Array(&eci,ecilist);
    if (clusterlist.get_size()!=eci.get_size()) ERRORQUIT("Number of ECI does not match number of clusters."); 
  }
/*//---------------------test---------------------------
//  LinkedList<MultiCluster> clusterlist;               
LinkedList<Cluster> *clusterlist;
  ifstream clusterfile("clusters.out");
// LinkedList<Real> ecilist; 
LinkedList<Real> *ecilist;
 Array<Real> eci; 
//Real read_clusters_and_eci(LinkedList<Cluster> *clusterlist, LinkedList<Real> *ecilist,                                                                            
//                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes)
read_clusters_and_eci(&clusterlist, &ecilist, clusterfile, ecifile, axes);
LinkedList_to_Array(&eci,ecilist);
if (clusterlist.get_size()!=eci.get_size()) ERRORQUIT("Number of ECI does not match number of clusters.");
//----------------------------------------------------
*/
//---------------count number of bonds in NN case-------------------------
//cout <<"I finished preparing the calculations, but do not start it yet!"<<'\n'; 

vector<vector<double>> sro;
vector<double> lro;
lro.resize(numEle);
vector<vector<vector<int>>> neighbor_list, neighbor_list_swapped;
if (NNeighbor==1){
//  vector<vector<double>> sro;
  neighbor_list=find_neighbor_list(str,latt,1,component,cryStruct);
  show_neighbors(neighbor_list);
  sro=short_range_order2(neighbor_list,component);
  //sro= short_range_order(str,latt,1,component,cryStruct);
  E0=cal_energy_nn(sro,eci_matrix, component, cryStruct);
  //lro= long_range_order(str,component,cryStruct);
  lro=long_range_order2(str,neighbor_list,component);
/*  cout<<"The initial structure has long range order parameters:"<<endl;
  for (int i=0;i<lro.size();i++){
    cout<<lro[i]<<"	"; 
  }
  cout<<"\n"; */
  //cout<<lro[0]<<" "<<lro[1]<<" "<<lro[2]<<" "<<lro[3]<<" "<<lro[4]<<endl;
}
else if(NNeighbor==0){
E0=cal_etot(str, clusterlist,spacegroup, labellookup, eci);
E0=E0*str.atom_pos.get_size(); // the energy by CE is per atom, has to change it;
}
Structure initial_struct;
//int bin_now;
int Ni=0; //record iteration step
vector<int> two_atom_ids(2);
int atom_a, atom_b;
int offset_a, offset_b; // for the choice of atom swap
int num_species=str.atom_pos.get_size()/NumEle;
Structure ideal_str;
const double k=8.617e-5; //eV/K
//double Ti;
vector<double> energies, temperatures;
vector<vector<double>>lros;
//lros.resize(N_final, vector<double> (numEle,0.0));
energies.resize(N_final);
temperatures.resize(N_final);
double r;
vector<vector<vector<double>>> sros;
//vector<vector<vector<int>>> neighbor_list, neighbor_list_swapped;
//neighbor_list=find_neighbor_list(str,latt,1,component,cryStruct);
while (Ni<N_final){
  temperatures[Ni]=Temp_scheme(Ni, N_final, NT_sample, T_init, T_final,Scheme);
  two_atom_ids=random_pick_two_numbers(component);
  atom_a=two_atom_ids[0];
  atom_b=two_atom_ids[1];
  r = ((double) rand() / (RAND_MAX));
  //ideal_str=swapped_str(atom_a,atom_b,str);
  if (NNeighbor){
    //sro= short_range_order(ideal_str,latt,1,component,cryStruct);
    neighbor_list_swapped=find_neighbor_list_swapped(neighbor_list,atom_a,atom_b);
    //show_neighbors(neighbor_list_swapped);
    sro=short_range_order2(neighbor_list_swapped,component);
    E1=cal_energy_nn(sro,eci_matrix, component, cryStruct);
    //cout<<"The two numbers should be the same: "<<E1-E0<<" "<<cal_dE_nn(neighbor_list,atom_a,atom_b,eci_matrix)<<endl;
    if (exp(-(E1-E0)/(k*temperatures[Ni]))>r){
      neighbor_list=neighbor_list_swapped;
      E0=E1;}
    energies[Ni]=E0;
    sro=short_range_order2(neighbor_list,component);
    sros.push_back(sro);
    lro= long_range_order2(str,neighbor_list,component);
    lros.push_back(lro);
  }
  else {
    ideal_str=swapped_str(atom_a,atom_b,str);
    E1=cal_etot(ideal_str, clusterlist,spacegroup, labellookup, eci)*str.atom_pos.get_size(); 
    if (exp(-(E1-E0)/(k*temperatures[Ni]))>r){
      str=ideal_str;
      E0=E1;
    }
    energies[Ni]=E0;
    sro=short_range_order(str, latt, 1, component,cryStruct);
    sros.push_back(sro);
    lro= long_range_order(str, component,cryStruct);
    lros.push_back(lro);
  }
  //srand(time(NULL));
  //r = ((double) rand() / (RAND_MAX));
  //cout <<"T, E1, E0, kT, rate, r:" <<temperatures[Ni]<<" "<<E1<<" "<<E0<<" "<<k*temperatures[Ni]<<" "<<exp(-(E1-E0)/(k*temperatures[Ni]))<<" "<<r<<endl;
/*
  if (exp(-(E1-E0)/(k*temperatures[Ni]))>r){ 
    //str=ideal_str;
    neighbor_list=neighbor_list_swapped;
    E0=E1;}
  energies[Ni]=E0;
//  cout<<"T vs energy: "<<temperatures[Ni]<<" "<<E0<<endl;
  sro=short_range_order2(neighbor_list,component);
  sros.push_back(sro);
  //cout<<"current_step:"<<" "<<Ni<<" "<<"swap a,b:"<<" "<<atom_a<<" "<<atom_b<<endl;
  lro= long_range_order2(str,neighbor_list,component);
  //lros[Ni]=lro;
  lros.push_back(lro);
  //cout<<"Steps:"<<"	"<<Ni<<"	"<<"Temperature:"<<"	"<<temperatures[Ni]<<"	"<<"energy:"<<"	"<<energies[Ni]<<'\n';
*/
  Ni +=1;
}
vector<double> specific_heat=cal_specific_heat(temperatures, energies, NT_sample,str.atom_pos.get_size());
vector<double> U=cal_internal_energy(temperatures, energies, NT_sample,str.atom_pos.get_size());
/* cout <<"Specific heat:"<<endl;
 for (int j=0; j<specific_heat.size();j++){
   cout<<temperatures[N_final/NT_sample*(j+1)-1]<<"	"<<specific_heat[j]<<endl;
 }
*/
const char *logfile;
//write energies
logfile="energies.dat";
ofstream logf(logfile);
for (int j=0;j<energies.size();j++){
  energies[j]=energies[j]/str.atom_pos.get_size();
  logf<<temperatures[j]<<"	"<<energies[j]<<endl;
}
logf.close();

// write specific heat
logfile="Cv.dat";
int step;
double Tj;
ofstream logf1(logfile);
step=temperatures.size()/specific_heat.size();
for (int j=0; j<specific_heat.size();j++){
  Tj=Temp_scheme(j*step, N_final, NT_sample,T_init,T_final,Scheme);
  logf1<<Tj<<"	"<<specific_heat[j]<<"	"<<U[j]<<endl;
} 
logf1.close();

//write long range order paramters
logfile="lro.dat";
ofstream logf2(logfile);  
int T_step=N_final/NT_sample;
vector<vector<double>> lroT;
lroT.resize(NT_sample, vector<double> (numEle,0.0));
//cout <<"T_step: "<<T_step<<"	"<<endl;
for (int j=0;j<lros.size();j++){
//  cout<<"j/T_step: "<<j/T_step<<" "<<"temperature: "<<temperatures[j]<<endl;
  for (int k=0;k<numEle;k++){
    if (j%T_step>T_step/2.0){
      //cout <<j<<" "<<NT_sample<<" "<< j%NT_sample <<"  "<<T_step/2.0<<endl;
      lroT[j/T_step][k] += 1.0*lros[j][k]/(T_step/2.0);}
  }
}
for (int j=0;j<NT_sample;j++){
  logf2 << temperatures[j*T_step]<<"   ";
  for (int k=0;k<numEle;k++){
   logf2 << lroT[j][k]<<"       ";}
  logf2 << '\n';}
logf2.close();

logfile="autocorr.dat";
vector<double> autocorr= cal_autocorrelation(energies);
ofstream logf3(logfile);   
for (int j=0;j<autocorr.size();j++){                    
  logf3 <<j<<"	"<< autocorr[j]<<endl; }                 
logf3.close();  

logfile="sro.dat";
ofstream logf4(logfile);
vector<vector<vector<double>>> sroT;
sroT.resize(NT_sample,vector<vector<double>>(numEle,vector<double>(numEle,0.0)));
for (int i=0;i<sros.size();i++){
  if (i%T_step>T_step/2.0-1) {
    for (int j=0;j<sros[0].size();j++){
      for (int k=0;k<sros[0][0].size();k++){
        sroT[i/T_step][j][k] += 1.0*sros[i][j][k]/(T_step/2.0);}}
  }
  if ((i+1)%T_step == 0){
    logf4<<temperatures[i]<<"	";
    for (int j=0;j<sros[0].size();j++){
      for (int k=0;k<sros[0][0].size();k++){
        logf4<<sroT[i/T_step][j][k]<<"	";}}
    logf4<<"\n";
  }
}
logf4.close();
/*
cout <<"Long range order:"<<endl;
vector<vector<double>>lroT;
lroT.resize(NT_sample, vector<double> (numEle,0.0));
lroT=cal_LRO(NT_sample,lros);
 for (int j=0; j<NT_sample;j++){
  for (int k=0;k<numEle;k++){
   cout <<lroT[j][k]<<"	";}
  cout<<endl;
}*/

if (NNeighbor){
//prepare str
double Box[3]={str.cell(0,0),str.cell(1,1),str.cell(2,2)};
Structure wstr=str;
for (int i=0;i<neighbor_list.size();i++){
  for (int j=0;j<neighbor_list[0][0].size()-1;j++){
    //cout<<i<<" "<<neighbor_list[i][neighbor_list[0].size()-1][j]<<endl;
    if (neighbor_list[i][neighbor_list[0].size()-1][j]==1){
      //cout<<i<<" "<<j<<endl;
      wstr.atom_type[i]=j;}
  }
}
Vector3d<double> mid;
int label_mid;
for (int i=0;i<wstr.atom_pos.get_size();i++){
  moveInBox(wstr.atom_pos(i),Box);
  while (wstr.atom_type[i] !=str.atom_type[i]){
    for (int j=i+1;j<wstr.atom_pos.get_size();j++){
      mid=wstr.atom_pos(i);
      wstr.atom_pos(i)=wstr.atom_pos(j);
      wstr.atom_pos(j)=mid;
      label_mid=wstr.atom_type[i];
      wstr.atom_type[i]=wstr.atom_type[j];
      wstr.atom_type[j]=label_mid;
    }
  }
}
//write restart.struct
ofstream wfile(restart_str_file);       
write_structure(wstr, label, axes, wfile);
wfile.close();
}
else{
ofstream wfile(restart_str_file);
write_structure(str, label, axes, wfile);
wfile.close();
}
return 0;
}

//----------------define functions here-----------------------
void show_neighbors(vector<vector<vector<int>>> neighbor_list){
for (int i=0;i<neighbor_list.size();i++){
  cout <<"Neighbor of atom "<<i<<endl;
  for (int j=0;j<neighbor_list[0].size();j++){
    for (int k=0;k<neighbor_list[0][0].size();k++){
      cout <<neighbor_list[i][j][k]<<"	";
    }
  cout<<endl;
  }
  cout << endl;
}
}

double cal_dE_nn(vector<vector<vector<int>>>neighbor_list, int a, int b, vector<vector<double>> eci_matrix){
double dE0,dE1,dE;
dE0=0.0;dE1=0.0;dE=0.0;
int ai,aj,bi,bj;
int num_neighbor=neighbor_list[0].size()-1;
int numEle=neighbor_list[0][0].size()-1;
for (int j=0;j<num_neighbor;j++){
  for (int k1=0;k1<numEle;k1++){
    for (int k2=0;k2<numEle;k2++){
      ai=neighbor_list[a][num_neighbor][k1];
      aj=neighbor_list[a][j][k2];
      dE0 += ai*aj*eci_matrix[k1][k2];
      bi=neighbor_list[b][num_neighbor][k1];
      bj=neighbor_list[b][j][k2];
      dE0 += bi*bj*eci_matrix[k1][k2];
      //energy after a, b exchanged
      dE1 += ai*bj*eci_matrix[k1][k2];
      dE1 += bi*aj*eci_matrix[k1][k2];
    }
  }
}
dE=dE1-dE0;
return dE;
}

double cal_energy_nn(vector<vector<double>> sro, vector <vector<double>> eci_matrix, vector<int> component, int cryStru){
  int numEle=component.size();
  int N_neighbor=0;
  if (cryStru==1){N_neighbor=8;}   
  else if (cryStru==2){N_neighbor=12;}                  
  else {cout<<"Crystal Structure not supported now!"<<endl;}
  //cal total atom number
  int sum=0;
  for (int i=0;i<numEle;i++){sum += component[i];}
/*
  // prepare the eci matrix
  vector <vector<double>> eci_matrix;
  int eci_id=0;
  eci_matrix.resize(numEle,vector<double>(numEle,0));
  for (int i=0;i<numEle;i++){
    for (int j=i+1;j<numEle;j++){
      eci_matrix[i][j]=eci_nn[eci_id];
      eci_id +=1;}
    for (int j=0;j<i;j++){
      eci_matrix[i][j]=eci_matrix[j][i];}  
  } */
  // cal energy
  double E0=0;
  for (int i=0;i<numEle;i++){
  for (int j=0;j<numEle;j++){          
      E0 += 0.5*sro[i][j]*eci_matrix[i][j]*component[i]*N_neighbor/sum;
  }}  
  return E0*sum; // because of the double counting of the bonds
}
double Temp_scheme(int Ni,int N_final, int NT_sample, double T_init, double T_final,int Scheme){
double dT;
if (NT_sample>1){
  dT=1.0*(T_final-T_init)/(NT_sample-1);
}
else{
  dT=T_final-T_init;
}
int dStep=N_final/NT_sample;
int i=Ni/dStep;
if (Scheme) {
  return T_final - i*dT;
}
else {
  return T_init + i*dT;
}
/*  double u_T=(T_final-T_init)/NT_sample;
  double u_T_step=(T_final-T_init)/N_final*2;
  int u_STEP=N_final/NT_sample;
 // if (Scheme){Ni=N_final-1-Ni;}
  int Ni_mod=Ni%u_STEP;
  int Ni_whole=Ni/u_STEP;
 // if (Scheme){
 //   if (Ni_mod<u_STEP/2){Ni_mod=0;}}
 // else{
  if (Ni_mod>u_STEP/2){Ni_mod=u_STEP/2;}//}
  if (Scheme){
    return T_final-(Ni_whole*u_T+Ni_mod*u_T_step);}
  else{
    return T_init+Ni_whole*u_T+Ni_mod*u_T_step;}
*/
}

vector<double> cal_autocorrelation(vector<double> energies){
int sampling_points=10;
int sample_size=energies.size()/sampling_points;
//double gamma0, gammat; 
int N=energies.size();
double Q_k_kt, Q_k, Q_2k,Q_k2;
vector<double> autocorr;
//double avg=0;
Q_k=0; Q_2k=0; Q_k2=0;
for (int i=0;i<energies.size();i++){
  Q_k +=energies[i]/N;
  Q_2k += pow(energies[i],2)/N;
}
Q_k2=pow(Q_k,2);

//Q_k=avg; ///energies.size();
/*
gamma0=0;
for (int i=0;i<energies.size();i++){
  gamma0 +=1.0/energies.size()*pow((energies[i]-avg),2);
} */
for (int i=0;i<sampling_points;i++){
  Q_k_kt=0; N=0; // Q_2k=0; Q_k2=0; //Q_k=0; N=0;
  //gammat=0.0;
  for (int j=0;j<(sampling_points-i)*sample_size;j++){     
    N +=1;
    Q_k_kt +=energies[j]*energies[j+i*sample_size];
   // Q_2k+=pow(energies[j],2);
    //Q_k+=energies[j];
    //cout<<Q_k_kt<<endl;
   // gamma0 +=1.0/energies.size()*pow((energies[i]-avg),2)
   /* if (i==sample_size-1){
      Q22_bar +=1.0/sample_size*pow(energies[0*sample_size+j],2);
      QQt_bar +=1.0/sample_size*energies[i*sample_size+j]*energies[0*sample_size+j];}
    else {*/
    //  gammat += 1.0/energies.size()*(energies[j]-avg)*(energies[j+i*sample_size]-avg);
     // Q22_bar +=1.0/sample_size*pow(energies[(i+1)*sample_size+j],2);
     // QQt_bar +=1.0/sample_size*energies[i*sample_size+j]*energies[(i+1)*sample_size+j];}
  }
  Q_k_kt=Q_k_kt/N; //Q_2k=Q_2k/N; Q_k2=pow(Q_k/N,2);
  cout<<Q_k_kt<<" "<<Q_2k<<" "<<Q_k<<"	"<<Q_k2<<"	"<<N<<endl;
 // autocorr.push_back(QQt_bar/sqrt(Q12_bar*Q22_bar));
  autocorr.push_back((Q_k_kt-Q_k2)/(Q_2k-Q_k2)); //(gammat/gamma0);
}
return autocorr;
}
vector<double> cal_specific_heat(vector<double>temperatures,vector<double> energies, int NT_sample, int size){
vector<double> specific_heat;
double E2=0;
double E =0;
double Z=0;
double k=8.6173303e-5;
const int steps=energies.size()/NT_sample/2;
double E_ref=energies[0];
if (energies[energies.size()-1]<E_ref){
  E_ref=energies[energies.size()-1];}
for (int i=0;i<energies.size();i++){
  if ((i%(2*steps))>(steps-1)){
    E2 += pow(energies[i],2)*exp((energies[i]-E_ref)/k/temperatures[i]);
    E += (energies[i])*exp((energies[i]-E_ref)/k/temperatures[i]);
    Z += exp((energies[i]-E_ref)/k/temperatures[i]);
  }
  if ((i+1)%(2*steps)==0) {
    //specific_heat.push_back((E2-pow(E,2))/k/pow(temperatures[i],2));
    specific_heat.push_back((E2/(Z)-pow(E/Z,2))/k/pow(temperatures[i],2));
    E2=0; E=0; Z=0;}
}
return specific_heat;
}

vector<double> cal_internal_energy(vector<double>temperatures,vector<double> energies, int NT_sample, int size){
vector<double> U;
double E =0;       
double Z=0;        
double k=8.6173303e-5;
const int steps=energies.size()/NT_sample/2;
double E_ref=energies[0];                   
if (energies[energies.size()-1]<E_ref){     
  E_ref=energies[energies.size()-1];} 
for (int i=0;i<energies.size();i++){
  if ((i%(2*steps))>(steps-1)){
    E += (energies[i])*exp((energies[i]-E_ref)/k/temperatures[i]);
    Z += exp((energies[i]-E_ref)/k/temperatures[i]);
  } 
  if ((i+1)%(2*steps)==0) {
    U.push_back(1.0*E/Z/size);
    E=0; Z=0;}
}
return U;
}

vector<vector<double>> cal_LRO(int NT_sample, vector<vector<double>>lros){
vector<vector<double>> lroT;                     
vector<double> avg(lros[0].size(),0.0);
//lroT.resize(NT_sample,vector<double>(lros[0].size(),0.0));
const int steps=lros.size()/NT_sample/2;     
for (int i=0;i<lros.size();i++){             
  if ((i%(2*steps))>(steps-1)){                  
    for (int k=0;k<lros[0].size();k++){          
      avg[k]+=lros[i][k]*1.0/steps; }
  }
  if ((i+1)%(2*steps)==0) {
    lroT.push_back(avg);
    for (int j=0;j<avg.size();j++){
      avg[j]=0;}}
    //avg.resize(lros[0].size(),0.0);}            
}                    
return lroT;
}
Vector3d<double> NNNbcc(double latt,int ithNN,int jth_ithNN){         
Vector3d<double> offset;                                              
switch(ithNN){ 
  case 1:
    switch(jth_ithNN){
      case 0: //offset(0.5*latt,0.5*latt,0.5*latt);      
        {offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;} 
      case 1:         
        {offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;} 
      case 2:         
        {offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;} 
      case 3:         
        {offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;} 
      case 4:         
        {offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;} 
      case 5:         
        {offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;} 
      case 6:         
        {offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;} 
      case 7:         
        {offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;} 
    }   
  case 2:
    exit;
}
//cout<<"I am offset: "<<jth_ithNN<<" "<<offset[0]<<" "<<offset[1]<<" "<<offset[2]<<endl;
return offset;
}

Vector3d<double> NNNfcc(double latt,int ithNN,int jth_ithNN){
Vector3d<double> offset;
switch(ithNN){                                          
  case 1:                  
    switch (jth_ithNN){                                   
     case 0: 
       offset[0]= 0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.0*latt; break;
     case 1:                                       
       offset[0]= 0.5*latt; offset[1]= 0.0*latt; offset[2]= 0.5*latt; break;
     case 2:                                       
       offset[0]= 0.0*latt; offset[1]= 0.5*latt; offset[2]= 0.5*latt; break;
     case 3:                                       
       offset[0]= 0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.0*latt; break;
     case 4:                                       
       offset[0]= 0.5*latt; offset[1]= 0.0*latt; offset[2]=-0.5*latt; break;
     case 5:                                       
       offset[0]= 0.0*latt; offset[1]= 0.5*latt; offset[2]=-0.5*latt; break;
     case 6:                                       
       offset[0]=-0.5*latt; offset[1]= 0.5*latt; offset[2]= 0.0*latt; break;
     case 7:                                       
       offset[0]=-0.5*latt; offset[1]= 0.0*latt; offset[2]= 0.5*latt; break;
     case 8:                                       
       offset[0]= 0.0*latt; offset[1]=-0.5*latt; offset[2]= 0.5*latt; break;
     case 9:                                       
       offset[0]=-0.5*latt; offset[1]=-0.5*latt; offset[2]= 0.0*latt; break;
     case 10:                                      
       offset[0]=-0.5*latt; offset[1]= 0.0*latt; offset[2]=-0.5*latt; break;
     case 11:                                      
       offset[0]= 0.0*latt; offset[1]=-0.5*latt; offset[2]=-0.5*latt; break;
    }
  case 2:
    exit;
}
return offset;
}

//vector<int>site_types(Structure str){

//}
vector<double> long_range_order2(Structure str,vector<vector<vector<int>>>neighbor_list, vector<int> component){
int numEle=component.size();
vector<double> lro, lro_alpha, lro_beta;
int current_label;
double atomX,rounded;
lro.resize(numEle,0.0);
lro_alpha.resize(numEle,0.0);
lro_beta.resize(numEle,0.0);
//#pragma omp parallel for private(current_label)
for (int i=0; i<str.atom_pos.get_size(); i++) {
  //current_label=str.atom_type[i];
  //atomX=str.atom_pos[i][0]-(int)str.atom_pos[i][0];
  //rounded = ((int)(atomX * 100 + .5) / 100.0);
  atomX=str.atom_pos[i][0]-round(str.atom_pos[i][0]);
  rounded=round(atomX*10+5)/10;
  //cout << rounded<<"	"<<str.atom_pos[i][0]<<endl;
  // if (atomX-floor(atomX) ==0.0) {
  if (rounded ==0.0) {
    //cout <<"check rounded==0.0"<<rounded<<" "<<str.atom_pos[i][0]<<endl;
    for (int j=0;j<neighbor_list[0][0].size()-1;j++){
      //cout <<"ID:"<<i<<" "<<neighbor_list[i][neighbor_list[0].size()-1][j]<<"	";
      lro_alpha[j] += neighbor_list[i][neighbor_list[0].size()-1][j];}
    //cout<<endl;
  }
  else {
    for (int j=0;j<neighbor_list[0][0].size()-1;j++){
      lro_beta[j] += neighbor_list[i][neighbor_list[0].size()-1][j];}
  }
}
//cout<<neighbor_list.size()<<" "<<neighbor_list[0].size()<<" "<<neighbor_list[0][0].size()<<endl;
for (int i=0; i<numEle;i++){
  //lro[i] =2.0*lro[i]/component[i]-1.0;}
  lro[i] = 1.0*(lro_alpha[i]-lro_beta[i])/(lro_alpha[i]+lro_beta[i]);
  //cout<<lro_alpha[i]<<" "<<lro_beta[i]<<endl;}
  }
return lro;
}

vector<double> long_range_order(Structure str, vector<int> component, int cryStru){
int numEle=component.size();
vector<double> lro;
int current_label;
double atomX,rounded;
lro.resize(numEle,0.0);
//#pragma omp parallel for private(current_label)
for (int i=0; i<str.atom_pos.get_size(); i++) {
 current_label=str.atom_type[i];
 atomX=str.atom_pos[i][0]-(int)str.atom_pos[i][0];
 rounded = ((int)(atomX * 100 + .5) / 100.0);
// if (atomX-floor(atomX) ==0.0) {
 if (rounded ==0.00) {
  lro[current_label] +=1;}
}
for (int i=0; i<numEle;i++){
  lro[i] =1.0*lro[i]/component[i];}
return lro;
}

vector<vector<vector<int>>> find_neighbor_list(Structure str, double latt, int ith_NN, vector<int> component, int cryStru){
vector<vector<vector<int>>> neighbor_list;
int numEle=component.size();
int size=str.atom_pos.get_size();
//vector<vector<double>> glo_avg_i_j;
//vector<vector<vector<double>>> loc_avg_i_j;
//loc_avg_i_j.resize(size,vector<vector<double>>(numEle,vector<double>(numEle,0)));
//glo_avg_i_j.resize(numEle,vector<double>(numEle,0));
//vector<double>glo_avg;
//glo_avg.resize(numEle*numEle,0.0);
//vector<vector<double>>loc_avg;
//loc_avg.resize(size,vector<double>(numEle*numEle,0.0));
//Vector3d<double> offset,neighbour;
double e1,e2;
double Box[3]={str.cell(0,0),str.cell(1,1),str.cell(2,2)};
int N_neighbor = 12;
if (cryStru==1){N_neighbor=8;}
else if (cryStru==2){N_neighbor=12;}
neighbor_list.resize(size,vector<vector<int>>((N_neighbor+1),vector<int>(numEle)));
#pragma omp parallel for private(e1,e2)
for (int i=0; i<str.atom_pos.get_size(); i++) {
  vector<int> occ_i(numEle+1);
  vector<vector<int>> occ_j;
  vector<vector<int>> neighbor_list_i;
  Vector3d<double> offset,neighbour;
  occ_j.resize(N_neighbor,vector<int>(numEle+1,0));
  e1=int(str.atom_type[i]);
  // initialize occ_i, occ_j
  for (int i1=0;i1<numEle;i1++){occ_i[i1]=0;}
  for (int j1=0;j1<N_neighbor;j1++){
    for (int i1=0;i1<numEle;i1++){occ_j[j1][i1]=0;}}
  occ_i[e1]=1; occ_i[numEle]=i;
  for (int j=0; j<N_neighbor; j++){
    // get the nearest neighbor
    if (cryStru==1) {offset=NNNbcc(latt,ith_NN,j);}
    else if (cryStru==2){offset=NNNfcc(latt,ith_NN,j);}
    else {cout<<"Please put a smaller number for crystal structure!"<<endl;
          exit;}
    for (int i2=0;i2<3;i2++){neighbour[i2] =str.atom_pos[i][i2]+offset[i2];}
    //cout <<"neighbor "<<j<<" "<<offset[0]<<" "<<offset[1]<<" "<<offset[2]<<endl;
    moveInBox(neighbour,Box);
    //cout <<"neighbor "<<j<<" "<<neighbour[0]<<" "<<neighbour[1]<<" "<<neighbour[2]<<endl;
    // find the label of the element occupying the NN
    int k=0;
    while(k<str.atom_pos.get_size() and !matchVector(str.atom_pos[k],neighbour,Box)){ 
      k+=1;}
    if (k==str.atom_pos.get_size()){cout<<"Warning: no match of neighbors found!"; 
       exit;}
    e2=int(str.atom_type[k]);
    occ_j[j][e2]=1; occ_j[j][numEle]=k;
    //cout<<"neighbor "<<j<<" "<<e2<<" "<<k<<" "<<str.atom_type[k]<<endl;
  }
  neighbor_list_i=occ_j;
  neighbor_list_i.push_back(occ_i); // remember the last element is for the site itself;
  neighbor_list[i]=neighbor_list_i;
  //cout<<neighbor_list_i[0][numEle]<<endl;
  //neighbor_list_i.resize(N_neighbor,vector<int>(numEle+1,0));
}
return neighbor_list;
}

vector<vector<vector<int>>> find_neighbor_list_swapped(vector<vector<vector<int>>> neighbor_list, int a, int b){
vector<vector<vector<int>>> neighbor_list_swapped=neighbor_list;
/* data structure of neighbor_list:
 * [atom_ID][neighbor_list_i]
 * neighbor_list_i:
 * [[0-num_of_neighbors-1][0-numEle-1 for labels, the last element is neighbor_ID]
 *  [num_of_neighbors for current central atom][labels and current atom_ID]]
 * removing the last column and last row gives occ_j for all neighbors' labels!!!
 */
//exchange the neighbor list of a and b
int num_neighbor=neighbor_list[0].size()-1;
int numEle=neighbor_list[0][0].size()-1;
int ID;
for (int i=0;i<numEle;i++){
  neighbor_list_swapped[a][num_neighbor][i]=neighbor_list[b][num_neighbor][i];
  neighbor_list_swapped[b][num_neighbor][i]=neighbor_list[a][num_neighbor][i];
}
for (int i=0;i<num_neighbor;i++){
  // update the matrix for neighbors of a
  ID=neighbor_list[a][i][numEle];
  for (int j=0;j<num_neighbor;j++){
    if (neighbor_list[ID][j][numEle]==a){
      for (int k=0;k<numEle;k++){
        neighbor_list_swapped[ID][j][k]=neighbor_list[b][num_neighbor][k];}
      //break;
    }
  }
  // update the matrix for neighbors of b
  ID=neighbor_list[b][i][numEle];
  for (int j=0;j<num_neighbor;j++){
    if (neighbor_list[ID][j][numEle]==b){
      for (int k=0;k<numEle;k++){
        neighbor_list_swapped[ID][j][k]=neighbor_list[a][num_neighbor][k];}
      //break;
    }
  }
}
return neighbor_list_swapped;
vector<vector<vector<int>>>().swap(neighbor_list_swapped);
}

vector<vector<double>> short_range_order2(vector<vector<vector<int>>> neighbor_list,vector<int> component){
int numEle=neighbor_list[0][0].size()-1;
vector<vector<double>> glo_avg_i_j;
glo_avg_i_j.resize(numEle,vector<double>(numEle,0));
vector<double>glo_avg;
glo_avg.resize(numEle*numEle,0.0);
vector<vector<double>>loc_avg;
int size=neighbor_list.size();
loc_avg.resize(size,vector<double>(numEle*numEle,0.0));
int N_neighbor=neighbor_list[0].size()-1;
/*
vector<int> occ_i(numEle);
vector<vector<int>> occ_j;
occ_j.resize(N_neighbor,vector<int>(numEle+1,0));
*/
#pragma omp parallel for 
for (int i=0; i<size; i++) {
  vector<int> occ_i(numEle);
  vector<vector<int>> occ_j;
  occ_j.resize(N_neighbor,vector<int>(numEle+1,0));
  occ_i=neighbor_list[i][N_neighbor];
  for (int j=0;j<N_neighbor;j++){
    occ_j[j]=neighbor_list[i][j];
    for (int elem1=0;elem1<numEle;elem1++){
    for (int elem2=0;elem2<numEle;elem2++){
      loc_avg[i][elem1*numEle+elem2] += occ_i[elem1]*occ_j[j][elem2];}}
  }
}
#pragma omp parallel for
for (int j=0;j<numEle*numEle;j++){
for (int i=0;i<size;i++){
  glo_avg[j] += loc_avg[i][j];
}
}

for (int elem1=0;elem1<numEle;elem1++){
for (int elem2=0;elem2<numEle;elem2++){
  glo_avg_i_j[elem1][elem2]=1.0*glo_avg[elem1*numEle+elem2]/(component[elem1]*N_neighbor);
}}
return glo_avg_i_j;
}

vector<vector<double>> short_range_order(Structure str, double latt, int ith_NN, vector<int> component, int cryStru){
int numEle=component.size();
int size=str.atom_pos.get_size();
vector<vector<double>> glo_avg_i_j;
glo_avg_i_j.resize(numEle,vector<double>(numEle,0));
vector<double>glo_avg;
glo_avg.resize(numEle*numEle,0.0);
vector<vector<double>>loc_avg;
loc_avg.resize(size,vector<double>(numEle*numEle,0.0));
double e1,e2;
double Box[3]={str.cell(0,0),str.cell(1,1),str.cell(2,2)};
int N_neighbor = 12;
if (cryStru==1){N_neighbor=8;}
else if (cryStru==2){N_neighbor=12;}
#pragma omp parallel for private(e1,e2) //,occ_i,occ_j,offset,neighbour)
for (int i=0; i<str.atom_pos.get_size(); i++) {
  vector<int> occ_i(numEle);                   
  vector<vector<int>> occ_j;
  Vector3d<double> offset,neighbour;
  occ_j.resize(N_neighbor,vector<int>(numEle,0));
  e1=int(str.atom_type[i]); //(str.atom_type[i]+1)%numEle;
  // initialize occ_i, occ_j
  for (int i1=0;i1<numEle;i1++){occ_i[i1]=0;}
  for (int j1=0;j1<N_neighbor;j1++){
  for (int i1=0;i1<numEle;i1++){occ_j[j1][i1]=0;}}
  occ_i[e1]=1;
  for (int j=0; j<N_neighbor; j++){     
    // get the nearest neighbor
    if (cryStru==1) {offset=NNNbcc(latt,ith_NN,j);}
    else if (cryStru==2){offset=NNNfcc(latt,ith_NN,j);}
    else {cout<<"Please put a smaller number for crystal structure!"<<endl;
          exit;}
    for (int i2=0;i2<3;i2++){neighbour[i2] =str.atom_pos[i][i2]+offset[i2];}
    moveInBox(neighbour,Box);    
    // find the label of the element occupying the NN
    int k=0;     
    while(k<str.atom_pos.get_size() and !matchVector(str.atom_pos[k],neighbour,Box)){ k+=1;}
    if (k==str.atom_pos.get_size()){cout<<"Warning: no match of neighbors found!"; 
       exit;}
    e2=int(str.atom_type[k]); //(str.atom_type[k]+1)%numEle;
    occ_j[j][e2]=1;
  }             
  for (int j=0;j<N_neighbor;j++){  
    for (int elem1=0;elem1<numEle;elem1++){        
    for (int elem2=0;elem2<numEle;elem2++){
      loc_avg[i][elem1*numEle+elem2] += occ_i[elem1]*occ_j[j][elem2];}}
      //cout <<occ_j[j][elem2]<<"	";}} //<<float(occ_j[j][elem2])<<"	";}}
      //cout<<elem1<<"	"<<elem2<<"	"<<occ_i[elem1]*occ_j[j][elem2]<<"	";}} 
    //cout<<"\n";
  }               
}     
#pragma omp parallel for
for (int j=0;j<numEle*numEle;j++){
for (int i=0;i<str.atom_pos.get_size();i++){
  glo_avg[j] += loc_avg[i][j];
}
}

for (int elem1=0;elem1<numEle;elem1++){
for (int elem2=0;elem2<numEle;elem2++){
  glo_avg_i_j[elem1][elem2]=1.0*glo_avg[elem1*numEle+elem2]/(component[elem1]*N_neighbor);
}}
//int total_bonds=0; 
/*
for (int i=0;i<str.atom_pos.get_size();i++){
for (int elem1=0;elem1<numEle;elem1++){        
for (int elem2=0;elem2<numEle;elem2++){
  glo_avg_i_j[elem1][elem2] += 1.0*loc_avg_i_j[i][elem1][elem2]/(component[elem1]*N_neighbor);}}
// total_bonds+=loc_avg_i_j[i][elem1][elem2]; }}
}*/
  //cout <<loc_avg_i_j[elem1][elem2]<<"	"; total_bonds+=loc_avg_i_j[elem1][elem2];} cout<<total_bonds<<endl;}
  //cout <<"TOTAL BONDS: "<<total_bonds/2.0<<"Should have: "<<str.atom_pos.get_size()*N_neighbor/2.0<<endl; 
return glo_avg_i_j; 
}


int matchVector(Vector3d<double> atom1, Vector3d<double> atom2, double Box[3]){
 int boole[3];
 for (int i=0;i<3;i++){
   if (abs(atom1[i] -atom2[i])<1e-4 or abs(abs(atom1[i] -atom2[i])-Box[i])<1e-4) 
     {boole[i]=1;} 
   else {boole[i]=0;}}
 return boole[0]*boole[1]*boole[2];
}

void moveInBox(Vector3d<double> &Atom, double Box[3]){
for (int i=0;i<3;i++){
 while (Atom[i]<0) {Atom[i] += Box[i];}
 while (Atom[i]>=Box[i]) {Atom[i] -= Box[i];}
}}

float cal_etot(Structure ideal_str,LinkedList<MultiCluster> clusterlist,SpaceGroup spacegroup, Array<Arrayint> labellookup, Array<Real> eci){
Real pred=0.;
  LinkedList<Array<MultiCluster> > eq_clusterlist;
  LinkedListIterator<MultiCluster> icluster(clusterlist);
  for ( ; icluster; icluster++) {  
    Array<MultiCluster> *pmulticlus=new Array<MultiCluster>;
    find_equivalent_clusters(pmulticlus, *icluster, spacegroup.cell, spacegroup.point_op, spacegroup.trans);
    eq_clusterlist << pmulticlus;                                                                                                                   
  } 
const char *corrfunc_label="trigo";
int multincl=0;
CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
pcorrfunc->init_from_site_type_list(labellookup);
LinkedListIterator<Array<MultiCluster> > icluster_eq(eq_clusterlist);
int ieci=0;
for ( ; icluster_eq; icluster_eq++, ieci++) {
  Real rho=calc_correlation(ideal_str, *icluster_eq, spacegroup.cell, *pcorrfunc);
//  if (strlen(ecifile)>0) {
    pred+=eci(ieci)*(multincl ? 1 : icluster_eq->get_size())*rho;
//  }        
//  else {   
//    if (!doconc) {cout << rho << delim;}
  }        
//if (strlen(ecifile)>0) {
  return float(pred);
//}
}
int genRandom(int range_from, int range_to){
std::random_device                  rand_dev;
std::mt19937                        generator(rand_dev());
std::uniform_int_distribution<int>  distr(range_from, range_to);
return distr(generator);
}
void writeLog(double gamma,array<double,binNum> dos, array<int,binNum>hist, array<int,binNum>count, double sro[binNum][NumEle], char *logfile){
ofstream logf(logfile);  
logf <<gamma<<'\n';                                                
logf <<"DOS:" <<'\n';                                                
for (int i = 0; i < binNum; i++){                                          
   logf<< dos[i] <<'\n';}                                                     
logf <<"HISTOGRAM:" <<'\n';                                          
for (int i = 0; i < binNum; i++){                                          
   logf<< hist[i] <<'\n';} 
logf <<"Short Range Parameter (sro):" <<'\n'; 
for (int i = 0; i < binNum; i++){
   for (int j=0; j < NumEle; j++){                                          
     logf<< sro[i][j]/count[i]<<",";}
     logf<<count[i]<<'\n';}
}
void readLog(double gamma,array<double,binNum> dos, array<int,binNum> hist, array<int,binNum> count, double sro[binNum][NumEle], char *logfile){
const char * split = ",";
char * p;    
char ch[101];
fstream logf(logfile,ios::in);   
int row=0;
  if (!logf) {cout<<"Sorry bad file"; exit(0);}
  while (!logf.eof())    
   {                    
     row=row+1;         
     logf.getline(ch,100);
     //cout <<row<<" "<<ch<<endl;
     if (row==1) { p=strtok(ch,split);  gamma=atof(p); }
     if (row>2 and row<binNum+3) { p=strtok(ch,split);  dos[row-3]=atof(p); }
     if (row>binNum+3 and row<2*binNum+4) { p=strtok(ch,split);  hist[row-binNum-4]=atof(p); }
     if (row>2*binNum+4 and row<3*binNum+5) { 
        p=strtok(ch,split); int i=0;           
        while (i<NumEle){                      
          sro[row-2*binNum-5][i]=atof(p);i++;} 
        count[row-2*binNum-5]=atoi(p);}
   }
}

//write structure
/*
void swap_write_struct(int atom_a, int atom_b, Structure &str, Array<AutoString> &atom_label, rMatrix3d &axes, char *file){
    Vector3d<double> mid;
    ofstream wfile(file); 
    mid=str.atom_pos(atom_a);
    str.atom_pos(atom_a)=str.atom_pos(atom_b);
    str.atom_pos(atom_b)=mid;
    write_structure(str,atom_label, axes, wfile);
}
*/
Structure swapped_str(int atom_a,int atom_b, Structure &str){
    Structure new_str=str;
    Vector3d<double> mid;
    mid=new_str.atom_pos(atom_a);
    new_str.atom_pos(atom_a)=new_str.atom_pos(atom_b);
    new_str.atom_pos(atom_b)=mid;
    return new_str;
}

void hello(int i){
cout <<"Hello!"<<"	"<<i<<'\n';
}
//total energy
float etot(const char* cmd) {
    array<char, 128> buffer;
    string result;
    shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != NULL)
            result += buffer.data();
    }
    return stof(result);
}
/*
vector<int> random_pick_two_numbers(vector<int> component){
  int atom_a, atom_b;
  int component_a, component_b;
  int id_lb_a, id_lb_b, id_ub_a, id_ub_b;
  id_lb_a=0; id_lb_b=0; 
  int numEle = component.size();
  vector<int> two_int;
  component_a=genRandom(0,numEle-1);    
  component_b=genRandom(0,numEle-1);    
  while (component_a==component_b){
    component_b=genRandom(0,numEle-1);
  } 
  if (component_a==0) {id_ub_a=component[component_a]-1;}
  else {  
    for (int j=0; j<component_a; j++){    
      id_lb_a += component[j];}  
    id_ub_a =id_lb_a +component[component_a]-1;    
  }  
  if (component_b==0) {id_ub_b=component[component_b]-1;}
  else {                         
    for (int j=0; j<component_b; j++){    
      id_lb_b += component[j];}  
    id_ub_b =id_lb_b+component[component_b]-1;
  }  
  atom_a= genRandom(id_lb_a,id_ub_a);     
  atom_b= genRandom(id_lb_b,id_ub_b);     
  two_int.push_back(atom_a);
  two_int.push_back(atom_b);
  return two_int;
} */
// new swap cheme
vector<int> random_pick_two_numbers(vector<int> component){
  int atom_a, atom_b;
  int component_a, component_b;
  int id_lb_a, id_lb_b, id_ub_a, id_ub_b;
  id_lb_a=0; id_lb_b=0;
  int numEle = component.size();
  int size=0;
  for (int i=0;i<numEle;i++){
   size += component[i];
  }
  vector<int> two_int;
  atom_a=genRandom(0,size-1);
  atom_b=genRandom(0,size-1);
  while (atom_a==atom_b){
    atom_b=genRandom(0,size-1);
  }
/*
  if (component_a==0) {id_ub_a=component[component_a]-1;}
  else {
    for (int j=0; j<component_a; j++){
      id_lb_a += component[j];}
    id_ub_a =id_lb_a +component[component_a]-1;
  }
  if (component_b==0) {id_ub_b=component[component_b]-1;}
  else {
    for (int j=0; j<component_b; j++){
      id_lb_b += component[j];}
    id_ub_b =id_lb_b+component[component_b]-1;
  }
  atom_a= genRandom(id_lb_a,id_ub_a);
  atom_b= genRandom(id_lb_b,id_ub_b);
*/
  two_int.push_back(atom_a);
  two_int.push_back(atom_b);
  return two_int;
}
int countSameElement(Arrayint labels,int elementLabel){                                                                                                                     
  int count=0;
  for (int i=0;i<labels.get_size();i++){
    //cout<<labels(i)<<endl;
    if (labels(i) == elementLabel){count +=1;}
  }  
  return count;
}
vector<string> split_fun(string str, char delimiter) {
  vector<string> internal;
  stringstream ss(str); // Turn the string into a stream.
  string tok; 
  while(getline(ss, tok, delimiter)) {
    internal.push_back(tok);
  }}
