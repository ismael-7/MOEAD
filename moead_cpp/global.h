#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>


using namespace std;

#include "random.h"
#include "SNP.h"


//-----datos de la instancia---
SNP instancia;

//----Paralelizacion---
int numHilos;

//----MersenneTwister---
MTRand r;

//----variables de operadores---
int probCross;
int lmut;
int fmut;
int probMutation;
int rangoMut;

//prueba///medir tiempo (borrar result luego.)
int contador=0;
double result=0;
clock_t t_ini,t_fin;

//------------- Parameters in test instance ------------------
int     nvar,  nobj=2, max_eval;                    //  the number of variables and objectives
char    strTestInstance[256];

//------------- Parameters in random number ------------------
long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------
vector <double> idealpoint; 

#endif
