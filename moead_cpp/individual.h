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
	void   obj_eval();
	void   show_objective();
	void   show_variable();

    bool   operator<(const CIndividual &ind2);
	bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};

CIndividual::CIndividual()
{
/*	solucion.tabu = vector<int>(DIM_EPI_MAX, 0);
    solucion.score = vector<int>(NUM_OBJETIVES, 0);
	solucion.rank = 0;*/
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

void CIndividual::obj_eval()
{
    //red bayesiana y verosimilitud
    //int i=0;

    solucion.score[0]=Bayesian_score(solucion.tabu,DIM_EPI,instancia);
//t_ini=clock();
    solucion.score[1]=logistic_score(solucion.tabu,DIM_EPI,instancia);
//t_fin=clock();
//result=result+(double)(t_fin-t_ini)/CLOCKS_PER_SEC;
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

/*
bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}*/

/*
bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}
*/

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

