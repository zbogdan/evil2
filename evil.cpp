////////////////////////////////////////////////
//
// (EVIL) Extremely-hard and Versatile Instance
// Library for Clique Search Benchmarks
//
// version: 1.0
// by: Bogdan Zavalnij
//
// Academic Free License 3.0
//
////////////////////////////////////////////////
//
// Given several input graphs and their multiplicity
// the program constructs a hard benchmark instance.
// The nodes of the different input graphs are connected with
// $p$ probability, then one maximum clique from each input
// graph connected all-together.
//
// The maximum clique size of the resulting graph is the
// sum of the clique sizes of the separate input graphs.
// If $p$ is close to $1$, then the chromatic number of the
// resulting graph is close to the sum of the chromatic
// numbers of the separate graphs.
//
// Thus, if graph instances where the clique size far
// from the chromatic number given, then the resulting
// graph will represent a hard maximum clique search
// benchmark instance.
//
//////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <vector>

using namespace std;

bool debug=false;

double p=.98;   // default edge probability
int ip=p*100;   // for file naming
vector<vector<int> > clx; //nodes of the embeded graphs' maximum clique

string comment;
string comment_pbm;
bool read_comments=true;

vector<vector<bool> > inc;        // the incidence matrix of input graphs
vector<vector<bool> > inc_out;    // the incidence matrix of output graph
int N_out;
int n0=0;
vector<int> solution;
vector<int> perm;

// search for maximum clique in the input graph
// simple Carraghan-Pardalos algorithm
void clique(int &N){
  vector<vector<bool> > stack_in;  // neigburs of the clique
  vector<vector<bool> > stack_out; // points in the clique
  vector<int> st_i;                // point examined
  int level;
  int i,j, max_level=0;

  stack_in.resize(N+1);
  stack_out.resize(N+1);
  st_i.resize(N+1);
  for(int i=0;i<N+1;++i){
    st_i[i]=0;
    stack_in[i].resize(N);
    stack_out[i].resize(N);
    for(j=0;j<N;++j){
      stack_in[i][j]=0;
      stack_out[i][j]=0;
    }
  }
  for(int i=0;i<N;++i)
      stack_in[0][i]=1;
  int level_count=N;

  //back tracking
  st_i[0]=-1;
  level=0;
  while(level>=0){
    if(level_count+level<max_level){
      //back
      level--;
    }else{
      //get next point
      ++(st_i[level]);
      for(; st_i[level]<N
            && (stack_in[level][st_i[level]]==0
                || stack_out[level][st_i[level]]==1)
            ; ++(st_i[level]));
      if(st_i[level]<N){
        //forward
	level_count=0;
	for(i=0;i<N;++i){
	  stack_in[level+1][i] = stack_in[level][i]&&inc[st_i[level]][i];
	  if(stack_in[level+1][i]) ++level_count;
	  stack_out[level+1][i]=stack_out[level][i];
	}
        stack_out[level+1][st_i[level]]=1;
        st_i[level+1]=st_i[level];
        ++level;
        if(level>max_level){
          max_level=level;
	  solution.clear();
	  for(i=0;i<N;++i)
	    if(stack_out[level][i]){
	      solution.push_back(i);
	    }
       }
      }
      else{
	level--;
      }
    }
  }
}

// connect the maximum cliques of the input graphs
// clx[][] contains the nodes of the maximum clique
void fillinc(){
  int i, j, xi, xj;

  if(debug){
    for(i=0;i<clx.size();++i){
      for(xi=0;xi<clx[i].size();++xi)
	cout<<clx[i][xi]<<", ";
      cout<<endl;
    }
  }

  for(i=0;i<clx.size();++i)
    for(j=i+1;j<clx.size();++j)

      for(xi=0;xi<clx[i].size();++xi)
	for(xj=0;xj<clx[j].size();++xj){
	  inc_out[clx[i][xi]][clx[j][xj]]=1;
	  inc_out[clx[j][xj]][clx[i][xi]]=1;
	}
}

// fill the incidence graph with
// k number of N sized input graphs
// inc[][] is the input graph
void fillemb(int N, int k){
  int i,j,z;
  for(z=0;z<k;++z)
    for(i=0;i<N;++i)
      for(j=0;j<N;++j)
        inc_out[n0+i+N*z][n0+j+N*z]=inc[i][j];
}

void read_clq(string file_name, int &N){
  int i,j,edge_number,x,y;
  string type,line;
  ifstream fin(file_name.c_str());
  if(!fin.is_open()){
    cout<<"ERROR opening file: "<<file_name<<endl;
    exit(1);
  }

  //eliminate comment lines
  while(fin.peek()=='c'){
    getline(fin,line);
    if(read_comments){
      comment=comment+line+'\n';
      comment_pbm=comment_pbm+"#"+line+'\n';
    }
  }
  fin>>type; // p
  fin>>type; // edge
  fin>>N; // vertexes
  fin>>edge_number; // edges

  inc.resize(N);
  for(i=0;i<N;i++)inc[i].resize(N);

  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      inc[i][j]=0;

  for(i=0;i<edge_number;i++){
    fin>>type; // e
    if(type=="e"){
      fin>>x;fin>>y;
      //cout<<"x: "<<x<<"; y: "<<y<<endl;
      inc[x-1][y-1]=1;
      inc[y-1][x-1]=1;
    }
  }
  fin.close();
}

void write_pbm(string file_name){
  int i,j;
  file_name=file_name+".pbm";
  ofstream fout(file_name.c_str());
  if(!fout.is_open()){
    cout<<"ERROR opening file"<<endl;
    exit(1);
  }
  fout<<"P1"<<endl;
  fout<<comment_pbm;
  fout<<N_out<<" "<<N_out<<endl;

  for(i=0; i<N_out; i++){
    for(j=0; j<N_out; j++){
      fout<<inc_out[i][j]<<" ";
    }
    fout<<endl;
  }  
  fout.close();
}
void write_clq(string file_name){
  int i,j;
  long long edge_number=0,tmp=0;

  ofstream fout(file_name.c_str());
  if(!fout.is_open()){
    cout<<"ERROR opening file"<<endl;
    exit(1);
  }

  for(i=0;i<N_out;i++)
    for(j=i+1;j<N_out;j++)
      edge_number += inc_out[i][j];

  fout<<comment;
  fout<<"p edge "<<N_out<<" "<<edge_number<<endl;

  for(i=0;i<N_out;i++)
    for(j=i+1;j<N_out;j++)
      if(inc_out[i][j]){
	fout<<"e "<<i+1<<" "<<j+1<<endl;
	++tmp;  
      }
  fout.close();
}

void permutation(){
  int k,r;
  srand(time(NULL));
  vector<int> nums;
  
  for(k=0;k<N_out;k++) nums.push_back(k);
  for(;k>0;k--){
    r=rand()%k;
    perm.push_back(nums[r]);
    nums.erase(nums.begin()+r);
  }
}

void permto(){
  vector<vector<bool> > tmp;
  int i,j;
  tmp.resize(N_out);
  for(i=0;i<N_out;i++)
    tmp[i].resize(N_out);

  for(i=0;i<N_out;i++)
    for(j=0;j<N_out;j++)
      tmp[i][j]=inc_out[perm[i]][perm[j]];
  for(i=0;i<N_out;i++)
    for(j=0;j<N_out;j++)
      inc_out[i][j]=tmp[i][j];
}

int main(int argc, char **argv){
  int i,j, l, N;
  bool good_format=true;
  vector<string> file_name;
  vector<int> k;
  vector<int> sol_tmp;
  vector<int> outN;
  int clique_sum=0;

  comment = comment + "c\nc (EVIL) Extremely-hard and Versatile Instance\n" +
    "c Library for Clique Search Benchmarks\nc\n" +
    "c version: 0.9\n" +
    "c by: Bogdan Zavalnij\nc\n" +
    "c Academic Free License 3.0\nc\n";
  comment_pbm = comment_pbm + "#\n# (EVIL) Extremely-hard" +
    " and Versatile Instance\n" +
    "# Library for Clique Search Benchmarks\n#\n" +
    "# version: 0.9\n" +
    "# by: Bogdan Zavalnij\n#\n" +
    "# Academic Free License 3.0\n#\n";
 
  for(i=1; i<argc; ++i)
    if(string(argv[i])==string("-p")){
      if(i+1>=argc){
	good_format=false;
	break;
      }
      ip=atoi(argv[++i]);
      p = (double)ip/100;
    }else if(string(argv[i])==string("-g")){
      if(i+3>=argc){
	good_format=false;
	break;
      }
      file_name.push_back(argv[++i]);
      if(string(argv[++i])==string("-k"))
	k.push_back(atoi(argv[++i]));
      else{
	good_format=false;
	break;
      }
    }else{
      good_format=false;
      break;
    }


  if(!good_format || k.size()==0){
    cerr<<"Usage: "<<argv[0];
    cerr<<" [-p PP] -g file.clq -k K [-g file.clq -k K]"<<endl;
    cerr<<endl<<"PP means edge probability in percent. (default: 98%)"<<endl;
    cerr<<"For each input graph a K multiplicity must be given!"<<endl;
    cerr<<endl<<"At least one input graph is needed!"<<endl<<endl;
    exit (1);
  }

  N_out=0;
  for(i=0;i<k.size();++i){
    read_clq(file_name[i], N);
    N_out += k[i]*N;
    outN.push_back(N);
  }
  read_comments=false;

  inc_out.resize(N_out);
  for(i=0;i<N_out;i++)
    inc_out[i].resize(N_out);

  // fill the adjacency matrix with edges of $p$ probability
  srand(time(NULL));
  for(i=0;i<N_out;++i)
    for(j=0;j<N_out;++j)
      if(i==j)
	inc_out[i][j]=0;
      else if(rand()>(RAND_MAX*p))
	{inc_out[i][j]=0; inc_out[j][i]=0;}
      else
	{inc_out[i][j]=1; inc_out[j][i]=1;}

  clx.clear();
  for(i=0;i<k.size();++i){
    read_clq(file_name[i], N);
    clique(N);
    clique_sum += solution.size()*k[i];
    for(j=0;j<k[i];++j){
      sol_tmp.clear();
      for(l=0;l<solution.size();++l)
	sol_tmp.push_back(solution[l] + N*j + n0);
      clx.push_back(sol_tmp);
    }
    // fill in by the sequence of input grtaphs
    fillemb(N, k[i]);
    n0 += k[i]*N;
  }

  long long edge_number=0;
  for(i=0;i<N_out;i++)
    for(j=i+1;j<N_out;j++)
      edge_number += inc_out[i][j];

  cout<<endl<<"The resulting graph has "<<N_out<<" nodes, ";
  cout<<edge_number<<" edges."<<endl<<"Its density is ";
  cout<<(double)edge_number/(N_out*(N_out-1)/2)<<"%"<<endl<<endl;
  cout<<"The size of its maximum clique is "<<clique_sum<<endl<<endl;


  ostringstream comm_ost;
  comm_ost << "|V|=" << N_out<< ", density=";
  comm_ost <<(double)edge_number/(N_out*(N_out-1)/2);
  comm_ost << "%, omega(E)="<<clique_sum;

  comment = comment + "c The evil graph:\nc " +
    comm_ost.str() +"\nc\n";
  comment_pbm = comment_pbm  + "# The evil graph:\n# " +
    comm_ost.str() +"\n#\n";

  if(n0!=N_out){
    cerr<<"Something wrong!"<<endl;
    exit (1);
  }

  // connect the maximum cliques of the sequence graphs
  fillinc();

  ostringstream ost;
  ost << "evil-N" << N_out<<"-p"<<ip;
  for(i=0;i<k.size();++i){
    ost<<"-"<<file_name[i].substr(0,3)<<outN[i]<<"x"<<k[i];
  }
  // write_pbm(ost.str()+"-noperm");
  // write_clq(ost.str()+"-noperm.clq");

  // permute the resulting graph
  permutation();
  permto();
  write_clq(ost.str()+".clq");
  cout<<"The name of the resulting graph is "<<ost.str()<<".clq"<<endl<<endl;
}
