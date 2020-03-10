#ifndef SNP_H
#define SNP_H

#include <cstdio>	//printf
#include <cmath>	//log
#include <vector>  	//std::vector
#include <limits>  	//to set infinity value
#include "MersenneTwister.h"
#include <omp.h>
#include <fstream>
#include <cstring>
#include <sstream>
#include <eigen3/Eigen/Dense>

#define NUM_OBJETIVES 2
#define DIM_EPI_MAX 10
extern int DIM_EPI;

using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;

class SNP { // Contiene los datos de la instancia
public:
	int samplesize; //Número de individuos de la instancia
	int locisize; // Número de SNPs de la instancia
	int data_col; // Número de columnas de la instancia (SNPs+Clase[0-control, 1-caso])
	int** data; // Datos de la instancia
	char** SNPnames; // Nombres de los SNPs
	//int classvalues[2];
	void input_data(const char* path); // Lee los datos del fichero path y los almacena en los atributos correspondientes
	void destroy();
};

int compareAIC(const void *a, const void *b); // Compara dos soluciones usando el objetivo de la verosimilitud; necesaria para hacer el qsort
int compareBayesian(const void *a, const void *b); // Compara dos soluciones usando el objetivo de la red bayesiana; necesaria para hacer el qsort
int compare (const void *a, const void *b); // Compara dos enteros; necesaria para hacer el qsort
double Bayesian_score(int* selectedSNPSet, int k, SNP SNPdata); // Calcula el valor de la función objetivo de la red bayesiana de un conjunto de SNPs dado (selectedSNPSet). k es el tamaño de la epistasia y SNPData los datos de la instancia
double logistic_score(int* selectedSNPSet, int k, SNP SNPdata); // Calcula el valor de la función objetivo de la verosimilitud de un conjunto de SNPs dado (selectedSNPSet). k es el tamaño de la epistasia y SNPData los datos de la instancia



struct solution_t { //Contiene los datos de una solución
	int tabu[DIM_EPI_MAX]; // Conjunto de SNPs de la solución. Tendrá valores entre 0 y X-1, siendo X el tamaño de la epistasia. Como máximo, podrá tener DIM_EPI_MAX
	int id; // Identificador de la solución
	double score[NUM_OBJETIVES]; // Valores de las funciones objetivos de la solución
	int rank;	// Ranking no dominado de la solución
	double distance; // Distancia de crowding de la solución

    //Sobrecarga del operador '=' (asignación)
    solution_t& operator=(const solution_t& a) {
		for(int sc = 0; sc < NUM_OBJETIVES; ++sc)
			score[sc] = a.score[sc];
		for(int tab = 0; tab < DIM_EPI; ++tab)
			tabu[tab] = a.tabu[tab];
		rank = a.rank;
		id =a.id;
		distance = a.distance;
		return *this;
    }

	/*
	*	Se dice que la solution A1 domina a la solución A2 si satisface las siguientes condiciones (funciones objetivo a minimizar):
	*	1. fx(A1)<=fx(A2), donde fx son las funciones objetivo del problema.
	*	2. fx(A1)<fx(A2) para al menos un objetivo del problema.
	*/
	// Sobrecarga del operador '<' (dominancia)
    bool operator<(const solution_t& a) const {
		int dominate = 0;
		for(int sc = 0; sc < NUM_OBJETIVES; ++sc)
		{
			if(score[sc] < a.score[sc])
			{
				dominate += 2;
			}
			else if(score[sc] == a.score[sc])
			{
				dominate += 1;
			}
		}
        return (dominate > 2);
    }
    
    bool operator==(const solution_t& a) {
        for(int tab = 0; tab < DIM_EPI; tab++)
            if (tabu[tab] != a.tabu[tab])
                return  false;
        return true;
     } 

     // Sobrecarga del operador '<<' (comparación de la distancia de crowding)
	bool operator<<(const solution_t& a) const {
		return (rank < a.rank or (rank == a.rank and distance > a.distance));
	}

	// Imprime los valores de las funciones objetivo de la solución
	void print_scores() const {
		printf("%.2f %.2f\n",score[0], score[1]);
	}
	
	// Imprime los SNPs de la solución
	void print_ind() const {
		for(int x=0; x<DIM_EPI;++x)
			printf("%d ",tabu[x]);
		printf("\n");
	}
	
	// Imprime los valores de ranking no dominado y distancia de crowding de la solución
	void print_rank_crow() const {
 		printf("%d %f\n",rank, distance);
	}
	
	// Indica si el snp pasado por parámetro está incluido en la solución excluyendo la posición pos
	bool buscar_snp(int snp, int pos) {
		int i;
		for (i=0; i<DIM_EPI; i++)
			if (snp==tabu[i] && pos!=i)
				return true;
		return false;
	}
	
	// Comprueba si la solución es válida o no (SNPs repetidos). En caso de no ser válida, la arregla. locisize es el número de SNPs de la instancia
	void Validar(int locisize) {
		int i;
		bool valido=true;
		for (i=0; i<DIM_EPI; i++)
			while (buscar_snp(tabu[i], i)) {
				valido=false;
				tabu[i]=(tabu[i]+1)%locisize;
			}
		if (!valido)
			std::qsort(tabu, DIM_EPI, sizeof(int), compare);
	}
};

//-------operadores------
bool crossover(solution_t *p1, solution_t *p2, solution_t *q1, solution_t *q2);

bool mutation(solution_t *q);

#endif // SNP_H
