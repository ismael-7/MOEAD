#include "SNP.h"

extern MTRand r;
extern int probMutation;
extern int probCross;
extern int rangoMut;
extern SNP instancia;
////
extern int numHilos;

using namespace std;
int DIM_EPI;

void SNP::destroy()
{
	printf("echo destroy\n");
	int i;
	for(i = 0; i < locisize; i++)
		delete []SNPnames[i];
	delete []SNPnames;
	for(i = 0; i < samplesize; i++)
		delete []data[i];
	delete []data;
}

void SNP::input_data(const char* path)
{
	
	int i,j,temp;
	string line;
	ifstream in(path);
	getline(in,line);
	istringstream test(line);
	i=0;
	string word;
	while(!test.eof())
	{
		getline(test,word,',');
		i++;
	}
	data_col = i;
	locisize = i-1;
	j=0;
	while(!in.eof())
	{   
		getline(in,line);
		j++;
	}
	
	samplesize = j-1;
	in.close();
	
	SNPnames=new char*[data_col];
	for(i=0;i<data_col;i++)
		SNPnames[i]=new char[50];
	data=new int* [samplesize];
	for(i=0;i<samplesize;i++)
		data[i]=new int[data_col];
	
	//classvalues[0]=0;
	//classvalues[1]=0;
	ifstream in1(path);
	getline(in1,line);
	istringstream test1(line);
	i=0;
	while(!test1.eof())
	{
		if(i==data_col)
			break;
		getline(test1,word,',');
		strcpy(SNPnames[i],word.c_str());
		i++;
	}
	
	
	i=0;
	while(!in1.eof())
	{
		if(i==samplesize)
		{
			break;
		}
		getline(in1,line);
		istringstream values(line);
		j=0;
		while(!values.eof())
		{
			getline(values,word,',');
			istringstream int_iss(word);
			int_iss>>temp;
			data[i][j++] = temp;
			/*if(j==data_col-1)
			{
				classvalues[temp]++;
			}*/

		}
		i++;
	}
	
	in1.close();
}

int compare (const void * a, const void * b) {
	return ( *(int*)a - *(int*)b );
}

int compareBayesian(const void * a, const void * b) {
	if ( (*(solution_t*)a).score[0] <  (*(solution_t*)b).score[0] ) return -1;
	if ( (*(solution_t*)a).score[0] == (*(solution_t*)b).score[0] ) return 0;
	if ( (*(solution_t*)a).score[0] >  (*(solution_t*)b).score[0] ) return 1;
	return 0;
}

int compareAIC(const void * a, const void * b) { 
	if ( (*(solution_t*)a).score[1] <  (*(solution_t*)b).score[1] ) return -1;
	if ( (*(solution_t*)a).score[1] == (*(solution_t*)b).score[1] ) return 0;
	if ( (*(solution_t*)a).score[1] >  (*(solution_t*)b).score[1] ) return 1;
	return 0;
}

//**************************************************************************
// K2 score
//**************************************************************************
double My_factorial(double n) {
	double i;
	double z = 0;
	if( n < 0 ) {
		printf("Illegal n, n should be a non-negative number[/*omp_get_thread_num()*/ 0].\n");
		return 0;
	}
	if( n == 0 )
		return 0;
	for( i = 1.0 ; i < ( n + 1.0 ) ; i ++ )
		z += log( i );
	return z;
}

double Bayesian_score(int* selectedSNPSet, int k, SNP SNPdata) {
	int comb = (int)pow(3.0, k);
	double observedValues[2][comb];
	double colSumTable[comb];
	int i,j,index;
	double score = 0;
	for(i=0;i<comb;i++) {
		observedValues[0][i] = 0;   
		observedValues[1][i] = 0;  
		colSumTable[i] = 0;         
	}
	/*constructing observed freq table*/
	bool cont;
	for(i=0;i<SNPdata.samplesize;i++) {
		index = 0;
		cont = 1;
		for(j=0;j<k;j++) {
			if(SNPdata.data[i][selectedSNPSet[j]] == 3) {
				cont = 0;
				break;
			}
			else
				index = index + SNPdata.data[i][selectedSNPSet[j]]*(int)pow(3.0,(k-1-j));
		}
		if(cont) {
			observedValues[SNPdata.data[i][SNPdata.data_col-1]][index]++;
			colSumTable[index]++;
		}
	}
	
	for (i=0; i<comb; i++)
		score=score+My_factorial(observedValues[0][i])+My_factorial(observedValues[1][i])-My_factorial(colSumTable[i]+1);
	score = fabs(score);
	return score;
}

//**************************************************************************
// AIC score
//**************************************************************************
double logistic_score(int* selectedSNPSet, int k, SNP SNPdata) {
	double delta = 0.001;
	int maxiter = 30;
	double aic = 0;
	double lossold,lossnew,loss;
	int i, j, s, theta_size, iter;
	int testsample = SNPdata.samplesize;
	theta_size = k+2;
	int newdata[testsample][theta_size];
	MatrixXd xtwx(theta_size,theta_size);
	// allocation
	double theta[theta_size];
	double ypre[testsample];
	double pi[testsample]; 
	double w[testsample];
	double wz[testsample];
	double xwz[theta_size];
	double pre[testsample];
	// initialization
	// par

//printf("tamaño de testsample: %d\n", testsample);
//printf("tamaño de theta_size: %d\n", theta_size);
#pragma omp parallel num_threads(numHilos)
{
	//int id= omp_get_thread_num();
	#pragma omp for private(i) //shared(newdata)
	for (i = 0; i < testsample; i++) {
		//printf("Soy el hilo %d, bucle 1.\n",id);
		newdata[i][0] = 1;
		newdata[i][theta_size-1] = 1;
		for(int h=1;h<theta_size-1;h++) {
//printf("Estoy en el bucle 1.2 par\n");
			newdata[i][h] = SNPdata.data[i][selectedSNPSet[h-1]];
			#pragma omp barrier
			newdata[i][theta_size-1] *=newdata[i][h];
		}
	}
//}
	for (i = 0; i < theta_size; i++)
		theta[i] = 0;
	// iteration
	iter = 0;
	loss = 1;
	lossnew = 0;
	while(loss > delta && iter < maxiter) {
		lossold = lossnew;
		for (s = 0; s < testsample; s++) {
			ypre[s] = 0;
			for (j = 0; j < theta_size; j++)
				ypre[s] += newdata[s][j] * theta[j];
			if (ypre[s]>= -708 && ypre[s] <= 709)
				pi[s] = exp(ypre[s])/(1+exp(ypre[s]));
			else if (ypre[s]< -708)
				pi[s]=0;
			else
				pi[s]=1;
			w[s] = pi[s]*(1-pi[s]);
			wz[s] = w[s] * ypre[s] + SNPdata.data[s][SNPdata.data_col-1]-pi[s];
		}
		xtwx.setZero();
		double red=0.0;
		for (i = 0; i < theta_size; i++)
			for (j = 0; j < theta_size; j++) {
				red=0.0;
				for (s = 0; s < testsample; s++)
					red += w[s] * newdata[s][i] * newdata[s][j];
				xtwx(i,j) = red;
			}
		xtwx = xtwx.inverse();
		for (i = 0; i < theta_size; i++) {
			xwz[i] = 0;
			for (s = 0; s < testsample; s++)
				xwz[i] += newdata[s][i] * wz[s];
		}
		
		for (i = 0; i < theta_size; i++) {
			theta[i] = 0;
			for (j = 0; j < theta_size; j++)
				theta[i] += xtwx(i,j)*xwz[j];
		}
		lossnew = 0;
		for (i = 0; i < theta_size; i++)
			lossnew += std::abs(theta[i]);
		loss = std::abs(lossnew - lossold);
		iter++;
	}
#pragma omp parallel num_threads(numHilos)
{
int id= omp_get_thread_num();
	#pragma omp for private(s) reduction(+:aic) 
	for (s = 0; s < testsample; s++) {
//printf("Soy el hilo %d, bucle 2.\n",id);
		//printf("Estoy en el bucle 2 par\n");
		pre[s] = std::abs(1 - SNPdata.data[s][SNPdata.data_col-1] - pi[s]);
		aic += -2*log(pre[s]);
	}
}
	//#pragma omp atomic
	aic = aic + 2*theta_size;
//}	
	return aic;
}


bool crossover(solution_t *p1, solution_t *p2, solution_t *q1, solution_t *q2) {
	int pto1, pto2, i;
	bool cruzado=false;
	if (r.randInt(100) < probCross) {
		if (DIM_EPI>2) {
			pto1=r.randInt(DIM_EPI-2);
			pto2=r.randInt(DIM_EPI-pto1-1)+pto1;
		}
		else {
			pto1=0;
			pto2=0;
		}
		for (i=0; i<DIM_EPI; i++)
			if (i<pto1 || i>pto2) {
				q1->tabu[i]=p1->tabu[i];
				q2->tabu[i]=p2->tabu[i];
			}
			else {
				q1->tabu[i]=p2->tabu[i];
				q2->tabu[i]=p1->tabu[i];
			}
		std::qsort(q1->tabu, DIM_EPI, sizeof(int), compare);
		std::qsort(q2->tabu, DIM_EPI, sizeof(int), compare);
		cruzado=true;
	}
	else {
		(*q1)=(*p1);
		(*q2)=(*p2);
		cruzado=false;
	}
	return cruzado;
}


bool mutation(solution_t *q) {
	int i, mut;
	bool mutado=false;
	for (i=0; i<DIM_EPI; i++)
		if (r.randInt(100) < probMutation) {
			mutado=true;
			do {
				mut=0;
				while (mut==0)
					mut=rangoMut-r.randInt(rangoMut*2);
				q->tabu[i]+=mut;
				if (q->tabu[i]<0)
					q->tabu[i]=0;
				else if (q->tabu[i]>instancia.locisize-1)
					q->tabu[i]=instancia.locisize-1;
			} while (q->buscar_snp(q->tabu[i], i)); //Se hace mientras que el snp generado no estuviera ya en el cromosoma
			std::qsort(q->tabu, DIM_EPI, sizeof(int), compare);
		}
	return mutado;
}

