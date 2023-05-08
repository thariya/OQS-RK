#include "stdafx.h"
#include <fstream>
#include <iterator>
#include <sys/time.h>
#include <sys/resource.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>

using namespace std;
using namespace boost::numeric::odeint;

//typedef runge_kutta_cash_karp54<state_type> error_stepper;
typedef runge_kutta4<state_type> general_stepper;
typedef runge_kutta_cash_karp54<state_type , double ,state_type , double ,openmp_range_algebra > error_stepper;
typedef controlled_runge_kutta<error_stepper> controlled_stepper;

double Pumping;

const struct rlimit limit={100000000,100000000};

state_type matrix;

int main(int argc, const char *argv[])
{
        setrlimit(RLIMIT_STACK , &limit); 	
	//int n=10;
	//int m=100;
	//double end_time=0.4;
	//cout<<"Initialized"<<endl;
	int n=atoi(argv[1]);
	int m=atoi(argv[2]);
	Pumping=atof(argv[3]);
	double end_time=atof(argv[4]);
	//int s=(n*n*n + 6 * n*n + 11 * n + 6)*(m + 1)*(m + 1) / 6;
	//state_type matrix(s);
	initialize(n,m);
	//cout<<"Initialized"<<endl;	
	matrix[0] = { 1.0, 0.0 };
	const double dt =1.0e-12;
	
	//const complex<double> abs_err(1.0e-6,0.0);
	
	//general_stepper stepper=general_stepper();
	error_stepper estepper = error_stepper();
	//controlled_stepper cstepper(default_error_checker<double,range_algebra,default_operations >(1.0e-6, 1.0e-6,1.0,0.0));
	//controlled_stepper cstepper(default_error_checker<complex<double>,range_algebra,default_operations>(complex<double>( 1.0e-6,0.0), complex<double>( 1.0e-6,0.0),complex<double>( 0.0,0.0),complex<double>( 1.0,0.0)));
	
	/*for (int i = 1; i < 1000; i = i + 1){
	stepper.do_step(transition, x, 0.0, dt);
	std::cout << x[0];
	std::cin.ignore();
	}*/

	//integrate_n_steps(stepper, transition, x, 0.0, dt, 100);
	//integrate_adaptive(cstepper, transition, matrix, 0.0, 1.0, dt, write_output);
	integrate_adaptive(make_controlled(1.0e-8, 1.0e-12,estepper), transition, matrix, 0.0, end_time, dt, write_output);
	//integrate_const(stepper,transition, matrix, 0.0, 1.0, dt, write_output);
	//integrate_const(bs_stepper, transition, matrix, 0.0, 1.0e-12, dt, write_output);
		
	for(int i=0;i<matrix.size();i++){
		std::cout<<i<<'\t'<<matrix[i]<<std::endl;
	}
	
        /*
	std::ofstream output_file("./data.txt");
	std::ostream_iterator<std::complex<double> > output_iterator(output_file, "\n");
	output_file.precision(10);
	std::copy(matrix.begin(), matrix.end(), output_iterator);

	std::ofstream output_file_error("./error.txt");
	std::ostream_iterator<std::complex<double> > output_iterator_error(output_file_error, "\n");
	output_file_error.precision(10);
	std::copy(err_matrix.begin(), err_matrix.end(), output_iterator_error);
	*/	
	return 0;
}

