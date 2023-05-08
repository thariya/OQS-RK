// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <omp.h>

#define values
#define time

typedef std::vector<std::complex<double> > state_type;

void initialize(int passed1, int passed2);
void transition(const state_type &x, state_type &dxdt, double t);
void write_output(const state_type &x, const double t);
void summer(int value);
int get_index(int a11, int a01, int a10, int p, int q);
void get_excitons(const state_type &x, std::complex<double> ex[]);
void get_plasmons(const state_type &x, std::complex<double> pl[]);
void init_vars(int passed1, int passed2);
void compute(int number, const state_type &x, state_type &dxdt);
void addElements(int start_row, int start_column, int end_row, int end_column);

const double hbar = 6.582119e-16;
const double kT = 25.7;  //kT at room temperature 25C in meV

const double omega_sp = 2.88773*1e-12 / hbar;
const double omega_x = omega_sp;
const double population = 1 / (exp(2887.73 / kT) - 1);

const double gamma_sp = 80;
const double gamma_x = 0.003;
const double gamma_pd = 3.0;
const double gamma_sp1 = (population + 1)*gamma_sp;
const double gamma_sp2 = population*gamma_sp;

extern double Pumping;

const double g = 19.7*1.0e-15/ hbar;

const int numcores=4;

extern state_type matrix;

// TODO: reference additional headers your program requires here
