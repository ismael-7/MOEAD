#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"
#include "SNP.h"



class CIndividual{
public:
	CIndividual();
	virtual ~CIndividual();

	solution_t solucion; //almacena los datos de una solucion
	int    count;

	void   rnd_init();
	void   obj_eval(double *time,std::fstream &fout);
	void   show_objective();
	void   show_variable();

    bool   operator<(const CIndividual &ind2);
	bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};

CIndividual::CIndividual()
{

}

CIndividual::~CIndividual()
{

}

void CIndividual::rnd_init()
{
      int i;
      for(i=0;i<DIM_EPI;i++)
        solucion.tabu[i]=rand() % instancia.locisize;   
      
      std::qsort(solucion.tabu, DIM_EPI, sizeof(int), compare);
      solucion.Validar(instancia.locisize);
}

void CIndividual::obj_eval(double *time,std::fstream &fout)
{
    //red bayesiana y verosimilitud
    solucion.score[0]=Bayesian_score(solucion.tabu,DIM_EPI,instancia);
    solucion.score[1]=logistic_score(solucion.tabu,DIM_EPI,instancia, time,fout);
}


void CIndividual::show_objective()
{
    solucion.print_scores();
}

void CIndividual::show_variable()
{
   solucion.print_ind();
}

 //Sobrecarga del operador '=' (asignación)
void CIndividual::operator=(const CIndividual &ind2)
{
   solucion.operator=(ind2.solucion);
}

// Sobrecarga del operador '<' (dominancia)
bool CIndividual::operator<(const CIndividual &ind2)
{
   solucion.operator<(ind2.solucion);
}

// defining subproblem 

class CSubproblem  
{
public:
	CSubproblem();
	virtual ~CSubproblem();

	void show();

	CIndividual     indiv;     // best solution
	CIndividual     saved;     // last solution
	vector <double> namda;     // weight vector
	vector <int>    table;     // neighbourhood table

	double          fitness;

    void  operator=(const CSubproblem &sub2);
};

CSubproblem::CSubproblem()
{
    namda = vector<double>(nobj, 0);
}

CSubproblem::~CSubproblem(){
}

void CSubproblem::show()
{
   for(int n=0; n<namda.size(); n++)
   {
       printf("%f ",namda[n]);
   }
   printf("\n");
}

void CSubproblem::operator=(const CSubproblem &sub2)
{
    indiv  = sub2.indiv;
	table  = sub2.table;
	namda  = sub2.namda;
}


#endif

