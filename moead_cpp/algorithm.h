#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <iomanip>
#include "global.h"
#include "common.h"
#include "individual.h"

class CMOEAD
{
public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_neighbourhood();               // calculate the neighbourhood of each subproblem
	void init_population(double *time,std::fstream &fout);                  // initialize the population

	//void load_parameter();

	void update_reference(CIndividual &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	void update_problem(CIndividual &ind, int &id, int &type); // update current solutions in the neighbourhood

	void evol_population(double *time,std::fstream &fout);                                      // DE-based recombination
	void mate_selection(vector<int> &list, int cid, int size, int type);  // select mating parents

	// execute MOEAD
	void exec_emo(string nombre_fichero);

	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_pos(char savefilename[1024]);

	void tour_selection(int depth, vector<int> &selected);
	void comp_utility();

    vector <CSubproblem> population;
	vector <double>      utility;

	void operator=(const CMOEAD &moea);

public:

	// algorithm parameters
    int		max_gen;       //  the maximal number of generations
	int     pops;          //  the population size
    int	    niche;         //  the neighborhood size
	int     limit;         //  the maximal number of solutions replaced
	double  prob;          //  the neighboring selection probability
//	double  rate;          //  the differential rate between solutions

	int     nfes;          //  the number of function evluations

};

CMOEAD::CMOEAD()
{

}

CMOEAD::~CMOEAD()
{

}

void CMOEAD::init_population(double *time,std::fstream &fout)
{
	idealpoint = vector<double>(nobj, 1.0e+30);
	utility    = vector<double>(pops, 1.0);

	char filename[1024];
	// Read weight vectors from a data file
	sprintf(filename,"Weight/W%dD_%d.dat", nobj, pops);
	std::ifstream readf(filename);

    for(int i=0; i<pops; i++)
	{

		CSubproblem sub;

		// Randomize and evaluate solution
		sub.indiv.rnd_init();
		sub.indiv.obj_eval(time,fout);

		sub.saved = sub.indiv;

		// Initialize the reference point
		update_reference(sub.indiv);

		// Load weight vectors
		for(int j=0; j<nobj; j++)
		{
		    readf>>sub.namda[j];
			//printf("%f ", sub.namda[j]);
		}
		//printf("\n"); getchar();

		// Save in the population
		population.push_back(sub);
		nfes++;
	}

	readf.close( );
}

void CMOEAD::operator=(const CMOEAD &alg)
{
	//population = alg.population;
}

void CMOEAD::init_neighbourhood()
{
    vector<double> dist   = vector<double>(pops, 0);
	vector<int>    indx   = vector<int>(pops, 0);

	for(int i=0; i<pops; i++)
	{
		// calculate the distances based on weight vectors
		for(int j=0; j<pops; j++)
		{
		    dist[j]    = dist_vector(population[i].namda,population[j].namda);
			indx[j]  = j;
		}

		// find 'niche' nearest neighboring subproblems
		minfastsort(dist,indx,population.size(),niche);

		// save the indexes of the nearest 'niche' neighboring weight vectors
		for(int k=0; k<niche; k++)
		{
			population[i].table.push_back(indx[k]);
		}

	}
    dist.clear();
	indx.clear();
}

void CMOEAD::tour_selection(int depth, vector<int> &selected)
{
	// selection based on utility
	vector<int> candidate;
	for(int k=0;    k<nobj; k++)    selected.push_back(k);   // select first m weights
	for(int n=nobj; n<pops; n++)    candidate.push_back(n);  // set of unselected weights

	while(selected.size()<int(pops/5.0))
	{
	    int best_idd = int(rnd_uni(&rnd_uni_init)*candidate.size()), i2;
		int best_sub = candidate[best_idd], s2;
		for(int i=1; i<depth; i++)
		{
		    i2  = int(rnd_uni(&rnd_uni_init)*candidate.size());
			s2  = candidate[i2];
			if(utility[s2]>utility[best_sub])
			{
				best_idd = i2;
			    best_sub = s2;
			}
		}
		selected.push_back(best_sub);
		candidate.erase(candidate.begin()+best_idd);
	}
}

void CMOEAD::comp_utility()
{
	double f1, f2, uti;
        double delta;
    for(int n=0; n<pops; n++)
	{
		f1 = fitnessfunction(population[n].indiv.solucion, population[n].namda);
		f2 = fitnessfunction(population[n].saved.solucion, population[n].namda);
		delta = (f2 - f1)/f2;
		//delta = (f2 - f1);
        if(delta>0.001)
		{
			utility[n] = 1.0;
		}
		else
		{
            //uti        = 0.95*(1.0+delta/0.001)*utility[n];
			//utility[n] = uti<1.0?uti:1.0;
			utility[n] = (0.95 + 0.05*delta/0.001)*utility[n];
		}
        population[n].saved = population[n].indiv;
	}
}

void CMOEAD::update_problem(CIndividual &indiv, int &id, int &type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if(type==1)	size = population[id].table.size();  // from neighborhood
	else        size = population.size();            // from whole population

	// a random order to update
	std::vector<int> perm(std::vector<int>(size, 0));
	for(int k=0; k<size; k++) perm[k] = k;

	//random_shuffle(perm.begin(), perm.end());
	// replaced by Aimin 2011.04.29
	permutation(perm);


    for(int i=0; i<size; i++)
	{
		// Pick a subproblem to update
		int k;
		if(type==1) k = population[id].table[perm[i]];
		else        k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
                ////fitnessfunction trabaja con double, pasar todo a entero??
		f1 = fitnessfunction(population[k].indiv.solucion, population[k].namda);
		f2 = fitnessfunction(indiv.solucion, population[k].namda);
		if(f2<f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if(time>=limit)
		{
			return;
		}
	}
	perm.clear();
}

void CMOEAD::update_reference(CIndividual &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.solucion.score[n]<idealpoint[n])
		{
			idealpoint[n] = ind.solucion.score[n];
		}
	}
}

void CMOEAD::mate_selection(vector<int> &list, int cid, int size, int type){
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss   = population[cid].table.size(), id, parent;
        int contador=0;
    while(list.size()<size)
	{
		if(type==1 && contador < ss){
		    id      = int(ss*rnd_uni(&rnd_uni_init));
			parent  = population[cid].table[id];
                      
                       // printf("busco en los vecinos, contador= %d,ss=%d\n",contador,ss);
		}
		else{
			parent  = int(population.size()*rnd_uni(&rnd_uni_init));
                        // printf("busco en toda la poblacion\n");
                }

                  contador++;
                  
                
                // avoid the repeated selection
		bool flag = true;
		for(int i=0; i<list.size(); i++)
		{
                    if(list[i]==parent || population[list[i]].indiv.solucion == population[parent].indiv.solucion) // parent is in the list
                    {
                        flag = false;
                        break;
                    }
		}

		if(flag || contador == pops ) list.push_back(parent);
	}
        
}

void CMOEAD::evol_population(double *time,std::fstream &fout)
{

	// random order of subproblems at each generation
	//vector<int> order(vector<int>(pops,0));
	//for(int i=0; i<pops; i++)  order[i] = i;
	//random_shuffle(order.begin(), order.end());

	vector<int> order;	this->tour_selection(10, order);

    for(int sub=0; sub<order.size(); sub++)
	{

		int c_sub = order[sub];    // random order
		int type;
        double rnd = rnd_uni(&rnd_uni_init);

		// mating selection based on probability
		if(rnd<prob)     type = 1;   // from neighborhood
		else             type = 2;   // from population
               
		// select the indexes of mating parents
		vector<int> plist;
               
		mate_selection(plist, c_sub, 2, type);  // neighborhood selection
		// produce a child solution
		CIndividual child,child1;
		//double rate2 = 0.5; //rate + 0.25*(rnd_uni(&rnd_uni_init) - 0.5);
                bool c = crossover(&population[plist[0]].indiv.solucion, &population[plist[1]].indiv.solucion, &child.solucion, &child1.solucion);
		plist.clear();
                
                
		// apply polynomial mutation
                if((child.solucion == population[plist[0]].indiv.solucion || child.solucion == population[plist[1]].indiv.solucion) && c)
                     child.solucion.tabu[1]=child1.solucion.tabu[0];
                
                //realizacion de la mutacion
		bool m = mutation(&child.solucion);
                
                if(c && !m)///con prob de cruce =0 c es siempre falso.
                    child.solucion.Validar(instancia.locisize);
                
		// evaluate the child solution
		child.obj_eval(time,fout);

		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		update_problem(child, c_sub, type);

		nfes++; //aumenta dentro del while de exec_emo
	}
}


void CMOEAD::exec_emo(string nombre_fichero)
{
	///--
	char filename_aic[1024];
	sprintf(filename_aic,"valores_aic/result.txt");
	std::fstream fout;
	fout.open(filename_aic,std::ios::out);
	///

double time = 0.0;
    char filename[1024];
	// initialization
	nfes      = 0;
        init_population(&time,fout);
        init_neighbourhood();

	int gen = 0;
	
	while(nfes < max_eval) //
	{
		evol_population(&time,fout);

		gen++;

		if(gen%50==0)
		{
		    comp_utility();
		}
		//printf("%d \n",nfes);
	}

	fout.close();

	sprintf(filename,"POFS/%s_%d_%d_%d.dat",nombre_fichero.c_str(),probCross,lmut,fmut);
	save_pos(filename);
	//sprintf(filename,"POF/%s_%d_%d_%d.dat",nombre_fichero.c_str() ,probCross,lmut,fmut);
	sprintf(filename,"POF/res_func_obj_%d.dat",numHilos);
	save_front(filename);

	population.clear();
	idealpoint.clear();
printf("Time in parallel section %lf\n", time);
}
void CMOEAD::save_front(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	//fout<<setprecision(3) << scientific;
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.solucion.score[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEAD::save_pos(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	//fout<<setprecision(3) << scientific;
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.solucion.tabu[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

#endif
