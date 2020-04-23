/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/
#include "time.h"
#include "global.h"
#include "algorithm.h"
#include "SNP.h"
#include <string>

extern int DIM_EPI;

int main(int argc,char *argv[])
{
std::chrono::time_point<std::chrono::system_clock> start, end;
start = std::chrono::system_clock::now();
std::chrono::duration<double> elapsed_segment;
    
    if(argc < 9){
      printf("error al lanzar la ejecución. La forma correcta es %s nombre_fichero proCruce lambamut factorMut DIM_EPI max_eval tamaño_poblacion numero_de_hilos\n",argv[0]);
      
      return 0;
    }
	///Paralelización
	numHilos=atoi(argv[8]);
        
        string nombre_fichero =argv[1];
        instancia.input_data(argv[1]);//fichero por parametro
        
//Parametros necesarios
        probCross = atoi(argv[2]);
        lmut = atoi(argv[3]);
        fmut = atoi(argv[4]);
        DIM_EPI=atoi(argv[5]);
        max_eval=atoi(argv[6]);
        probMutation = 100/DIM_EPI*(100+lmut)/100;
        rangoMut = (instancia.locisize*fmut)/100;
	
	// random number init
        rnd_uni_init = 1000000;//time(NULL);
	r.seed(rnd_uni_init);
        // the parameter setting of test instance
        nvar=DIM_EPI;
        CMOEAD MOEAD;
        MOEAD.pops=atoi(argv[7]);
        MOEAD.niche=MOEAD.pops/10;
        MOEAD.limit=MOEAD.pops/100;
        if(MOEAD.limit == 0)
            MOEAD.limit=1;
        MOEAD.prob=0.7;
       // MOEAD.load_parameter();
        MOEAD.exec_emo(nombre_fichero);//run por parametro
end = std::chrono::system_clock::now();
elapsed_segment = end - start;
printf("Total time: %lf\n", elapsed_segment.count());
	return 0;
}
