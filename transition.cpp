#include "stdafx.h"
#include <queue>
#include <iterator>
#include <unordered_map>
#include <map>
#include <omp.h>
#include <complex>
#include <stdlib.h>


bool columnSort(const std::vector<int>& lhs, const std::vector<int>& rhs);

std::vector< std::complex<double> > A;
std::vector< int > IA;
std::vector< int > JA;

std::vector< std::complex<double> > A_short;
std::vector< int > IA_short;
std::vector< int > JA_short;

std::unordered_map<int,int> mapper;
std::unordered_map<int,int> mapper2;

std::vector<int> indexer;
std::vector<int> actual_indices;
std::vector<int> diagonal;
std::vector<int> trace_indices;
std::vector<int> critical_indices;
std::vector<int> reverse_indices;

std::vector<int> degrees;

std::complex<double> summing( 0.0, 0.0 );
double clock_start=omp_get_wtime();

int N;
int M;
int Size;
int problem_size;
int actual_size;
int photon_size;
int size[5];
int nnz;

std::queue<int> Q;
std::vector<int> Order;
std::vector<int> S;

double START;

void bfs(){
	
	nnz=0;
	Q.push(0);
	Order.push_back(0);
	S.resize(Size,-1);
	
	S[0]=0;
	
	nnz=1;	
	
	while(!Q.empty()){
		int current=Q.front();
		Q.pop();
		std::vector<std::vector<int> > temp;
		for(int i=IA[current];i<IA[current+1];i++){
			std::vector<int> tempRow={JA[i],degrees[JA[i]]};
			temp.push_back(tempRow);			
			//std::cout<<tempRow.size()<<std::endl;

		}
		
		std::sort(temp.begin(),temp.end(),columnSort);
		
		for(int i=0;i<temp.size();i++){
			int j=temp[i][0];
			if(S[j]==-1){
				S[j]=nnz;
				nnz++;
				Q.push(j);
				Order.push_back(j);	
			}
		}		
	}
	//std::cout<<critical_indices.size()<<std::endl;
	critical_indices.reserve(nnz);
	for(int i=0; i < S.size();i++){
		if(S[i]>=0)	critical_indices.push_back(i);
	}
	//std::cout<<"dione"<<std::endl;
	/*
	for(int i=0;i<Order.size();i++){
		//std::cout<<i<<" "<<critical_indices[i]<<" "<<reverse_indices[critical_indices[i]]<<std::endl;
		std::vector<int> temp=reverse(reverse_indices[Order[i]]);		
		std::cout<< temp[0]<<"  "<< temp[1]<<"  "<< temp[2]<<"  "<< temp[3]<<"  "<< temp[4]<<std::endl;
	}
	*/	
	std::cout<<"Size: "<<nnz<<std::endl;

	//for(int i=0;i<Order.size();i++)	std::cout<<Order[i]<<std::endl;	
}
	
std::vector<int> reverse(int index){
	std::vector<int> indices(5);
	
	std::vector<int> s={(N+1)*(N+1)*(M+1)*(M+1),(N+1)*(M+1)*(M+1), (M+1)*(M+1) ,(M+1) ,1};
	for(int i=0;i<5;i++){
    		int d=floor(index/s[i]);
    		indices[i]=d;
    		index=index-d*s[i];
	}   
	return indices;
} 


bool columnSort(const std::vector<int>& lhs, const std::vector<int>& rhs)
{
       	return lhs[1] < rhs[1];
}


void transition(const state_type &x, state_type &dxdt, double t)
{	
	double start = omp_get_wtime();
		
	omp_set_dynamic(0);
	omp_set_num_threads(24);
	#pragma omp parallel for schedule(static,4000)	
	for (int pos = 0; pos < x.size(); pos = pos + 1){
		std::complex<double> sum (0.0,0.0);
		for (int k = IA_short[pos]; k < IA_short[pos + 1]; k = k + 1){
			sum += A_short[k] * x[JA_short[k]];
		}
		dxdt[pos]=sum;
		int thread=omp_get_thread_num();
		//std::cout<<"This is thread number: "<<thread<<std::endl;
	}
	
	
	#ifdef derivative
    	std::cout << "derivative: "<<'\t' << dxdt[0] << std::endl;
	#endif

	#ifdef time
	double end = omp_get_wtime();
	std::cout << "transition done in " << double(end - start) << std::endl;
	#endif
}

void init_vars(int passed1,int passed2){
	N = passed1;
	M = passed2;
	problem_size = (N + 1)*(N + 1)*(N + 1)*(M + 1);
	actual_size = (N*N*N + 6 * N*N + 11 * N + 6)*(M + 1)/ 6;
	photon_size = (M + 1);
	size[0]= (N + 1)*(N + 1)*(M + 1);
	size[1]= (N + 1)*(M + 1);
	size[2]= (M + 1);
	size[3]= 1;
	size[4]=0;	
	//std::cout<<"done"<<std::endl;
}

void initialize(int passed1, int passed2){
	START=omp_get_wtime();
	//std::cout<<"Checkpoint 1"<<std::endl;
	init_vars(passed1,passed2);
	//std::cout<<"Initialized!"<<std::endl;
	int index = 0;
	const std::complex< double > I(0.0, 1.0);
	std::complex< double > sum(0.0, 0.0);
	int p, q, r = 0;
	int temp = 0;
	std::cout<<actual_size<<std::endl;
	IA.reserve(actual_size+1);
	std::cout<<"Checkpoint 2; "<<" Max Size:"<<A.max_size()<<std::endl;
	long long reservation=13*(long long)actual_size;
	std::cout<<reservation<<std::endl;
	A.reserve(reservation);
	std::cout<<"OK"<<std::endl;
	JA.reserve(reservation);
	
	degrees.reserve(actual_size);
	
	indexer.resize(problem_size, 0);
	actual_indices.resize(problem_size,0);
	trace_indices.reserve((N+1)*(M+1));
	diagonal.reserve((N + 1)*(M + 1));
	
	for(int a = 0; a < N+1 ; a++ ){
		for(int b = a; b < N+1 ; b++ ){
			for(int c = b; c < N+1 ; c++ ){
				for(int d = 0; d < M+1 ; d++ ){
					p = a;
					q = b - a;
					r = c - b;
					
					int e;
					if(q>r)	e=d+(q-r);
					else	e=d-(r-q);
					if(e>=0&&e<M+1){
						index = p*size[0] + q*size[1] + r*size[2] + d*size[3];							
						if (q == 0 && r == 0 && d == e){
							diagonal.push_back(temp);
						}
						mapper[index]=temp;
						temp++;
					}									
				}
			}
		}
	}
	
	
	
	IA.push_back(0);
	Size=temp;
	int ind_ia=0;
	std::cout<<"Size:"<<temp<<std::endl;
	clock_t start = clock();
	for (int a = 0; a < N + 1; a = a + 1){
		for (int b = a; b < N + 1; b = b + 1){
			for (int c = b; c < N + 1; c = c + 1){
				for (int d = 0; d < M + 1; d = d + 1){
					//for (int e = 0; e < M + 1; e = e + 1){
					p = a;
					q = b - a;
					r = c - b;
						
					int e;
					if(q>r)	e=d+(q-r);
					else	e=d-(r-q);
						
					if( e>=0 && e<M+1){
						int degree=0;
						index = p*size[0] + q*size[1] + r*size[2] + d*size[3];							
						//if(index==112211)	std::cout<<"DAng "<<p<<q<<r<<d<<e<<std::endl;
						
						
						//D10

						if (a > 0){
							temp = mapper[index - size[0]];
							sum = (N - c + 1)*Pumping;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}				


						if (p > 0 && e < M){
							temp = mapper[index - size[0] + size[1]];
							sum = -(q + 1)*sqrt(e + 1)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
							
						}

						//HI

						if (p > 0 && d < M){
							temp = mapper[index - size[0] + size[2] + size[3]];
							sum = (r + 1)*sqrt(d + 1)*g * I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						if (q > 0 && d < M){
							temp = mapper[index - size[1] + size[3]];
							sum = (N - c + 1)*sqrt(d + 1)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						if (r > 0 && e < M){
							temp = mapper[index - size[2]];
							sum = -(N - c + 1)*sqrt(e + 1)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}


						
						//Db*						
						if (d > 0 && e>0){
							temp = mapper[index - size[3] - size[4]];
							sum = sqrt(d*e)*gamma_sp2;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						
						//Middle
						temp = mapper[index];
						sum = (((e - d)*omega_sp + (r - q)*omega_x) * I - (p + c)*gamma_x / 2 - (2 * N - (p + c))*Pumping / 2 - (c - p)*gamma_pd - (d + e)*gamma_sp1 / 2 - (d + e + 2)*gamma_sp2 / 2);
						A.push_back(sum);
						JA.push_back(temp);
						ind_ia++;

						//Db

						if (d < M && e < M){
							temp = mapper[index + size[3] + size[4]];
							sum = sqrt((d + 1)*(e + 1))*gamma_sp1;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						if (c < N && e>0){
							temp = mapper[index + size[2] - size[4]];
							sum = -(r + 1)*sqrt(e)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						if (c < N && d>0){
							temp = mapper[index + size[1] - size[3]];
							sum = (q + 1)*sqrt(d)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						if (r > 0 && d > 0){
							temp = mapper[index + size[0] - size[2] - size[3]];
							sum = (p + 1)*sqrt(d)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						//HI
						if (q > 0 && e > 0){
							temp = mapper[index + size[0] - size[1] - size[4]];
							sum = -(p + 1)*sqrt(e)*g *I;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}

						//D01						
						if (c < N){
							temp = mapper[index + size[0]];
							sum = (p + 1)*gamma_x;
							A.push_back(sum);
							JA.push_back(temp);
							ind_ia++;
						}
						degrees.push_back(ind_ia-IA.back());
						IA.push_back(ind_ia);
						//std::cout << index << std::endl;
					}
				}
			}
		}
	}
	/*
	FILE *fp;
	fp=fopen("A.txt", "w");
	int j=0;
	for(j=0;j<A.size();j++)
	{
		fprintf(fp,"%g  %g\n",real(A[j]),imag(A[j]));
	}
	fclose(fp);
	
	fp=fopen("IA.txt", "w");
	for(j=0;j<IA.size();j++)
	{
		fprintf(fp,"%d\n",IA[j]);
	}
	fclose(fp);
	
	fp=fopen("JA.txt", "w");
	for(j=0;j<JA.size();j++)
	{
		fprintf(fp,"%d\n",JA[j]);
	}
	fclose(fp);
	
	std::cout<<"Done Printing"<<std::endl;
	*/	
	bfs();
	
	//std::cout<<"done"<<std::endl;
	long long reservation1=13*(long long)Order.size();
	A_short.reserve(reservation1);
	IA_short.reserve(Order.size()+1);
	JA_short.reserve(reservation1);
	
	IA_short.push_back(0);
	
	for(int a = 0; a < N+1 ; a++ ){
		for(int b = a; b < N+1 ; b++ ){
			for(int c = b; c < N+1 ; c++ ){
				for(int d = 0; d < M+1 ; d++ ){
					p = a;
					q = b - a;
					r = c - b;
					
					int e;
					if(q>r)	e=d+(q-r);
					else	e=d-(r-q);
					if(e>=0&&e<M+1){
						index = p*size[0] + q*size[1] + r*size[2] + d*size[3];							
						mapper2[index]=S[mapper[index]];						
					}									
				}
			}
		}
	}
	
	for(int i=0; i< JA.size(); i++){
		JA[i]=S[JA[i]];				
	}
	
	for(int i=0 ; i< Order.size(); i++){
		std::map<int , std::complex<double> > tempvec;
		for(int j=IA[Order[i]];j<IA[Order[i]+1];j++){
			tempvec[JA[j]]= A[j];						
		}
		std::map<int,std::complex<double> >::iterator it;
		for(it=tempvec.begin();it!=tempvec.end();it++){
			A_short.push_back(it->second);
			JA_short.push_back(it->first);
			//std::cout<< it->first <<std::endl;
		}
		//std::cout<<"ENDL" <<std::endl;
		IA_short.push_back(IA_short.back()+tempvec.size());	
	}
	
	for(int i=0; i<N+1 ; i++){
		for(int j=0; j< M+1 ; j++){
			trace_indices.push_back(mapper2[size[0]*i + size[3]*j]);
		}
	}	

	std::vector<std::complex<double> >().swap(A);
	std::vector<int>().swap(IA);
	std::vector<int>().swap(JA);
	
	
	matrix.resize(Order.size());
	double end = omp_get_wtime();
	std::cout << "house keeping done in " << double(end - start)  << std::endl;
}


void write_output(const state_type &x, const double t)
{
	summing = 0;
	for_each(trace_indices.begin(), trace_indices.end(), summer);
	#ifdef time
	std::cout << "Total Time " << double(omp_get_wtime() - START) << std::endl;
	std::cout << "cycle done in " << double(omp_get_wtime() - clock_start) << std::endl;
	#endif

	#ifdef values
	
	std::cout << t <<'\t'<< "sum is : "<<summing <<std::endl;
	/*
	std::cout << matrix[actual_indices[size[0]*0]] <<'\t'<<matrix[actual_indices[size[0]*1]] <<'\t'<< matrix[actual_indices[size[0]*2]] <<'\t'<<matrix[actual_indices[size[0]*3]] <<'\t'<<matrix[actual_indices[size[0]*4]] <<'\t'<<matrix[actual_indices[size[0]*5]] <<'\t'<<matrix[actual_indices[size[0]*6]] <<'\t'<<matrix[actual_indices[size[0]*7]] <<'\t'<<matrix[actual_indices[size[0]*8]] <<'\t'<<matrix[actual_indices[size[0]*9]] <<'\t'<<matrix[actual_indices[size[0]*10]] <<'\t'<<std::endl<<std::endl;
	std::cout << matrix[actual_indices[size[3]+size[4]]]<<std::endl<<std::endl;
	*/
	std::complex<double> ex[N+1];
	std::complex<double> excitons(0.0,0.0);
	std::complex<double> g_ex(0.0,0.0);
	get_excitons(x,ex);
	for(int i=0 ; i<N+1 ; i++){
		std::cout<<ex[i]<<'\t';
		excitons+=ex[i]*(double)i;
		g_ex+=ex[i]*double(i)*double(i-1);
	}	
	g_ex=g_ex/(excitons*excitons);
	std::cout<<std::endl;
	std::cout<<"EXCITONS: "<<excitons<< std::endl;
	std::cout<<"EXCITON G: "<<g_ex<<std::endl;
	
	std::complex<double> pl[M+1];
	std::complex<double> plasmons(0.0,0.0);
	std::complex<double> g_pl(0.0,0.0);
	get_plasmons(x,pl);
	for(int i=0 ; i<M+1 ; i++){
		std::cout<<pl[i]<<'\t';
		plasmons+=pl[i]*(double)i;
		g_pl+=pl[i]*double(i)*double(i-1);
	}
	g_pl=g_pl/(plasmons*plasmons);	
	std::cout<<std::endl;
	std::cout<<"PLASMONS: "<<plasmons<< std::endl;
	std::cout<<"PLASMON G: "<<g_pl<<std::endl;
	
	std::cout<<std::endl<<std::endl;
	#endif
	clock_start = omp_get_wtime();
}

void summer(int value){
	summing += matrix[value];
}

int get_index(int a11,int a01,int a10, int p, int q){
	return actual_indices[a11*size[0]+a01*size[1]+a10*size[2]+p*size[3]+q*size[4]];	
}

void get_excitons(const state_type &x, std::complex<double> ex[]){
	for (int a = 0; a < N + 1; a = a + 1){
		ex[a]=std::complex<double>(0.0,0.0);
		for (int d = 0; d < M + 1; d = d + 1){
			int index = a*size[0] + d*size[3];							
			ex[a]=ex[a]+x[mapper2[index]];			
		}
	}	
}
void get_plasmons(const state_type &x, std::complex<double> pl[]){
	for (int d = 0; d < M + 1; d = d + 1){
		pl[d]=std::complex<double>(0.0,0.0);
		for (int a = 0;a < N + 1; a = a + 1){
			int index = a*size[0] + d*size[3];							
			pl[d]=pl[d]+x[mapper2[index]];			
		}
	}	
}

