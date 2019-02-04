#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "borgms.h"
#include "fex.h"
typedef struct weightEntry *node;

struct weightEntry{
	double *weights;
	double *objectives;
	struct weightEntry *next;
	struct weightEntry *previous;
	struct weightEntry *head;
};

node removeEntry(node p)
{
	node next = NULL;
	if (p != NULL)
	{
		if (p->next != NULL)
		{
			if (p->previous != NULL)
			{
				p->next->previous = p->previous;
				p->previous->next = p->next;
				next = p->next;
			}
			else
			{
				p->next->previous = NULL;
				node p_ = p->next;
				next = p->next;
				while (p_ != NULL)
				{
					p_->head = p->next;
					p_ = p_->next;
				}
			}
		}
		else
		{
			p->previous->next = NULL;
		}
		
	}
	
	free(p->weights);
	free(p->objectives);
	free(p);
	return next;
}
node createEntry(){
	node temp;
	temp = (node)malloc(sizeof(struct weightEntry));
	temp->next = NULL;
	temp->previous = NULL;
	temp->head = NULL;
	return temp;
}

node addEntry(node head, double *weights, double *objectives){
	node temp, p;
	temp = createEntry();
	temp->weights = weights;
	temp->objectives = objectives;
	if (head == NULL){
		head = temp;
	}
	else
	{
		p = head;
		while (p->next != NULL){
			p = p->next;
		}
		p->next = temp;
		temp->previous = p;
	}
	temp->head = head;
	return head;
}

node resetNode(node head){
	node current = head;
	node next;
	while (current != NULL) {
		next = current->next;
		free(current->weights);
		free(current->objectives);
		free(current);
		current = next;
	}
	return current;
}
int main(int argc, char* argv[]) {
	int i, j , k;
	int simType = atoi(argv[2]), framing = atoi(argv[1]);
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize); //worldSize and worldRank are globals declared by borgms.c or h .. I'm hijacking them so I can use them outside of Borg as well as inside
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	// int numWeights, simType = 0, framing = 6, numFramings = 6;
	int numWeights, numFramings = 8;
	
	//-------------------------------------------------------------------//
	// Setup 
	//-------------------------------------------------------------------//
	
	struct fexEnsemble myEnsemble;
	myEnsemble.ensembleSize = 1;
	myEnsemble.ensembleT = 1000;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 50;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = C1 only, 4 = un-cooperative
	createEnsemble(&myEnsemble, framing);	
	
	double *epsilons = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
	double *epsilonUB = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
	double *epsilonLB = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
	char epstrIni[100];
	
	int betaInd, numNetHidden;
	switch (simType)
	{
		case 0:
			numWeights = myEnsemble.ensemble[0]->cooperativeNet->numWeights;
			betaInd = myEnsemble.ensemble[0]->cooperativeNet->numNetInputs*myEnsemble.ensemble[0]->cooperativeNet->numNetHidden;
			numNetHidden = myEnsemble.ensemble[0]->cooperativeNet->numNetHidden;
			
			epsilons[0] = 0.11; // D2 def
			epsilons[1] = 0.13; // C1 def
			epsilons[2] = 1.3; // c2£
			epsilons[3] = 9; // C1£
			
			epsilonUB[0] = epsilons[0]*30;
			epsilonUB[1] = epsilons[1]*30;
			epsilonUB[2] = epsilons[2]*30;
			epsilonUB[3] = epsilons[3]*30;
			
			epsilonLB[0] = epsilons[0]/10;
			epsilonLB[1] = epsilons[1]/10;
			epsilonLB[2] = epsilons[2]/10;
			epsilonLB[3] = epsilons[3]/10;
			
			sprintf(epstrIni,"%f_%f_%f_%f",epsilons[0],epsilons[1],epsilons[2],epsilons[3]);
			
			if (framing < 5)
			{
				printf("Error: cooperative simulation is not intended for uncooperative framings\n");
				return EXIT_FAILURE;
			}
			break;
		case 1:
			numWeights = myEnsemble.ensemble[0]->c2Net->numWeights;
			betaInd = myEnsemble.ensemble[0]->c2Net->numNetInputs*myEnsemble.ensemble[0]->c2Net->numNetHidden;
			numNetHidden = myEnsemble.ensemble[0]->c2Net->numNetHidden;
			
			epsilons[0] = 0.11; // D2 def
			epsilons[1] = 1.3; // c2£
			
			epsilonUB[0] = epsilons[0]*30;
			epsilonUB[1] = epsilons[1]*30;
			
			epsilonLB[0] = epsilons[0]/10;
			epsilonLB[1] = epsilons[1]/10;
			
			sprintf(epstrIni,"%f_%f",epsilons[0],epsilons[1]);
			
			if (framing > 4)
			{
				printf("Error: uncooperative simulation is not intended for cooperative framings\n");
				return EXIT_FAILURE;
			}
			break;
		case 2:
			numWeights = myEnsemble.ensemble[0]->c1Net->numWeights;
			betaInd = myEnsemble.ensemble[0]->c1Net->numNetInputs*myEnsemble.ensemble[0]->c1Net->numNetHidden;
			numNetHidden = myEnsemble.ensemble[0]->c1Net->numNetHidden;
			
			epsilons[0] = 0.17; // C1 def
			epsilons[1] = 9; // C1£
			
			
			epsilonUB[0] = epsilons[0]*30;
			epsilonUB[1] = epsilons[1]*30;
			
			epsilonLB[0] = epsilons[0]/10;
			epsilonLB[1] = epsilons[1]/10;
			
			sprintf(epstrIni,"%f_%f",epsilons[0],epsilons[1]);
			
			if (framing > 4)
			{
				printf("Error: uncooperative simulation is not intended for cooperative framings\n");
				return EXIT_FAILURE;
			}
			break;
		case 3:
			printf("Error: uncooperative simulation is not intended for optimization\n");
			return EXIT_FAILURE;
	}
	destroyEnsemble(&myEnsemble);
	int NFE, NRandom;
	int nseeds;
	BORG_Problem problem;
	BORG_Archive result;
	struct timeval tv1, tv2;
	struct rusage ru_before, ru_after;
	//-------------------------------------------------------------------//
	
	
	//-------------------------------------------------------------------//
	// Run actual optimization
	//-------------------------------------------------------------------//
	myEnsemble.ensembleSize = 400;
	myEnsemble.ensembleT = 300;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 50;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = C1 only, 4 = un-cooperative
	// myEnsemble.includeSimulationData = 0; //This can be = 0 or 24 - you must not include simulation data if running the ensemble more than once
	
	createEnsemble(&myEnsemble, framing);
	NFE = 1000;
	// Create the rbf_sim problem, defining the number of decision variables,
	// objective and constraints.  The last argument, rbf_sim, references
	// the function that evaluates the rbf_sim problem.
	problem = BORG_Problem_create(numWeights, myEnsemble.ensemble[0]->numObjectives, 0, evaluateEnsemble, &myEnsemble);
	// Set the lower and upper bounds for each decision variable.
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -1, 1);
	}
	for (i=betaInd; i< betaInd + numNetHidden ; i++){
		BORG_Problem_set_bounds(problem, i, 0, 1);
	}
	// Set the epsilon values used by the Borg MOEA.  Epsilons define the
	// problem resolution, which controls how many Pareto optimal solutions
	// are generated and how far apart they are spaced. Set num function evals
	for (i = 0; i < myEnsemble.ensemble[0]->numObjectives; i++)
	{
		BORG_Problem_set_epsilon(problem, i, epsilons[i]);
	}
	
	
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);
	
	// Run the Borg MOEA 
	char runtimeFinal[512];
	sprintf(runtimeFinal,"optimization_framing%d_simType%d_ensembleSize%d_T%d_NFE%d_epsilons%s.runtime",framing,myEnsemble.simulationType,myEnsemble.ensembleSize,myEnsemble.ensembleT,NFE,epstrIni);
	BORG_Algorithm_output_runtime(runtimeFinal); //create runtime (i.e. evolution of performance and weights over time
	result = BORG_Algorithm_ms_run(problem);
	if (result != NULL)
	{
		BORG_Archive_destroy(result);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	BORG_Problem_destroy(problem);
	
	MPI_Finalize();
	return EXIT_SUCCESS;
	
	//-------------------------------------------------------------------//
	// Run fex with random weights (to show benefits of optimization) - maybe print only the paretofront please (to epsilon?)
	//-------------------------------------------------------------------//
	if (framing > 4)
	{
		NRandom = 100000;
	}
	else
	{
		NRandom = 100;
	}
	char randomRuntime[256];
	FILE *fidRandomRuntime;

	getrusage (RUSAGE_SELF, &ru_before);

	node randSet = NULL;
	for (i = 0;i < NRandom; i++)
	{
		double *weights = calloc(numWeights,sizeof(double));
		double *objectives = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
		srand(i);
		//e.g. if 10 processors are available, each processor evaluates the performance of every 10th random weight
		for (j = 0; j < numWeights; j++)
		{
			
			//You can move this into if (j%worldSize..) statement if you are not certain that random numbers generated on different processors with the same seed will be identical (tests show this is the case for blue crystal)
			weights[j] = (2*(double)rand()/(double)RAND_MAX)-1; //random weight between -1 and 1
			if (j >= betaInd)
			{
				if (j< betaInd + numNetHidden)
				{
					weights[j] = (double)rand()/(double)RAND_MAX; //random weight between 0 and 1
				}
			}
		}
		if (i%worldSize == worldRank) 
		{
			
			evaluateEnsemble(weights,objectives,0,&myEnsemble);
			if (worldRank != 0)
			{
				MPI_Send(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
			}
		}
		if (worldRank == 0)
		{
			if (i%worldSize != worldRank)
			{
				MPI_Recv(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,i%worldSize,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			randSet = addEntry(randSet,weights,objectives);
		}
		// if (worldRank == 0)
		// {
			// fprintf(fidRandomRuntime,"\n");
		// }
		if (worldRank != 0)
		{
			free(weights);
			free(objectives);
		}
	}
	
	getrusage (RUSAGE_SELF, &ru_after);
	tv1 = ru_before.ru_utime;
	tv2 = ru_after.ru_utime;
	if (worldRank == 0)
	{
		
		
		// int *idxs = calloc(NRandom,sizeof(int));
		// for (i=0;i<NRandom;i++)
		// {
			// idxs[i] = i;
		// }
		node p1 = randSet;
		while (p1 != NULL)
		{
			node p2 = p1->head;
			while (p2 != NULL)
			{
				int index = 0;
				for (i = 0; i < myEnsemble.ensemble[0]->numObjectives; i++)
				{
					if (p1->objectives[i] <= p2->objectives[i])
					{
						index++;
					}
				}
				if (index == myEnsemble.ensemble[0]->numObjectives && p1 != p2)
				{
					p2 = removeEntry(p2);
					// printf("I removed %d \n",l);
					// l++;
					// fflush(stdout);
				}
				else
				{
					p2 = p2->next;
				}
			}
			randSet = p1->head;
			p1 = p1->next;
		}
		sprintf(randomRuntime,"random_weights_framing%d_simType%d_ensembleSize%d_T%d_NRandom%d.runtime",framing,myEnsemble.simulationType,myEnsemble.ensembleSize,myEnsemble.ensembleT,NRandom);
		fidRandomRuntime = fopen(randomRuntime,"w");
		fprintf(fidRandomRuntime,"//Have a nice day, today we are using %d processors\n",worldSize);
		fprintf(fidRandomRuntime,"// total time elapsed = %f (s)\n",tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
		fflush(fidRandomRuntime);
		p1 = randSet;
		while (randSet != NULL)
		{
			for (i = 0; i < numWeights ; i++)
			{
				fprintf(fidRandomRuntime, "%f ", randSet->weights[i]);
			}
			for (i = 0; i < myEnsemble.ensemble[0]->numObjectives ; i++)
			{
				fprintf(fidRandomRuntime, "%f", randSet->objectives[i]);
				if (i != myEnsemble.ensemble[0]->numObjectives - 1)
				{
					fprintf(fidRandomRuntime, " ");
				}
				else
				{
					fprintf(fidRandomRuntime, "\n");
				}
			}
			randSet = randSet->next;
		}
		fprintf(fidRandomRuntime, "#");
		fclose(fidRandomRuntime);
		p1 = resetNode(p1);
		free(p1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	destroyEnsemble(&myEnsemble);
	//-------------------------------------------------------------------//
	// Check convergence of objectives 
	// IF YOU WANT TO USE A LARGER ENSEMBLESIZE THAN 1 YOU HAVE TWO OPTIONS:
	// 1) CREATE AN ARRAY OF ENSEMBLES - THIS WILL BE QUICK BUT MEMORY INTENSIVE
	// 2) CREATE AND DESTROY A NEW ENSEMBLE ON EACH ITERATION INSIDE THE J LOOP - THIS WILL BE SLOWER BUT MORE MEMORY EFFICIENT
	//-------------------------------------------------------------------//
	NRandom = worldSize;
	int T0 = 1000, TMax = 400000, Nsteps = 5;
	
	myEnsemble.ensembleSize = 1000;
	myEnsemble.ensembleT = T0;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 50;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = C1 only, 4 = un-cooperative
	createEnsemble(&myEnsemble, framing);
	char convergence[256];
	FILE *fidConvergence;
	
	if (worldRank==0)
	{
		sprintf(convergence,"convergence_framing%d_simType%d_ensembleSize%d_TMax%d_T0%d_Nsteps%d_NRandom%d.runtime",framing,myEnsemble.simulationType,myEnsemble.ensembleSize,TMax,T0,Nsteps,NRandom);
		fidConvergence = fopen(convergence,"w");
		fprintf(fidConvergence,"//Note, this isn't your average .runtime file!\n");
		fprintf(fidConvergence,"//Each new line is the weights and objectives associated with a simulation of the same weights over simulation length increasing with logspace(T0,TMax,Nsteps) file\n");
		fprintf(fidConvergence,"//The # separates each iteration of the experiment with different weights, today we are using %d processors\n",worldSize);
	}
	
	
	double totaltime = 0;
	for (i = 0; i < Nsteps ;i++)
	{
		getrusage (RUSAGE_SELF, &ru_before);
		
		double *objectives = calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
		
		for (j = 0; j < NRandom; j++)
		{
			
			
			myEnsemble.ensembleT = pow(10,log10(T0) + (log10(TMax) - log10(T0))*i/(Nsteps-1)); //logspace
			//THIS TAKES A CHUNK OF TIME AND SHOULD BE MOVED INSIDE if (j%worldSize == worldRank)
			//#######
			// for (k = 0; k < myEnsemble.ensembleSize; k++)
			// {
				// repopulateFex(myEnsemble.ensemble[k], &myEnsemble,(Nsteps*j+i)*myEnsemble.ensembleSize + k +1);
			// }
			//#######
			srand(j+1);
			double *weights = calloc(numWeights,sizeof(double));
			for (k = 0; k < numWeights; k++)
			{
				//You can move this into if (j%worldSize..) statement if you are not certain that random numbers generated on different processors with the same seed will be identical (tests show this is the case for blue crystal)
				weights[k] = (2*(double)rand()/(double)RAND_MAX)-1; //random weight between -1 and 1
				if (k >= betaInd)
				{
					if (k< betaInd + numNetHidden)
					{
						weights[k] = (double)rand()/(double)RAND_MAX; //random weight between 0 and 1
					}
				}
			}
			// myEnsemble.ensemble[0]->fexT = T0 + (TMax - T0)*i/(Nsteps-1); //linspace
			
			if (j%worldSize == worldRank)
			{
				for (k = 0; k < myEnsemble.ensembleSize; k++)
				{
					repopulateFex(myEnsemble.ensemble[k], &myEnsemble,(Nsteps*j+i)*myEnsemble.ensembleSize + k +1);
				}
				evaluateEnsemble(weights,objectives,0,&myEnsemble);
				if (worldRank != 0)
				{
					MPI_Send(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
				}
			}
			if (worldRank == 0)
			{
				if (j%worldSize != 0)
				{
					MPI_Recv(objectives,myEnsemble.ensemble[0]->numObjectives,MPI_DOUBLE,j%worldSize,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				for (k = 0; k < numWeights; k++)
				{
					fprintf(fidConvergence,"%f ", weights[k]);
				}
				for (k = 0; k < myEnsemble.ensemble[0]->numObjectives; k++)
				{
					fprintf(fidConvergence,"%f", objectives[k]);
					if (k == (myEnsemble.ensemble[0]->numObjectives - 1))
					{
						fprintf(fidConvergence,"\n");
					}
					else
					{
						fprintf(fidConvergence," ");
					}
				}
				fflush(fidConvergence);
			}
			free(weights);
		}
		if (worldRank == 0)
		{
			getrusage (RUSAGE_SELF, &ru_after);
			tv1 = ru_before.ru_utime;
			tv2 = ru_after.ru_utime;
			totaltime += tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6;
			fprintf(fidConvergence,"# //time elapsed = %f (s)\n",totaltime);
		}
		free(objectives);
		
	}
		
	if (worldRank == 0)
	{
		// fprintf(fidConvergence,"#");
		fclose(fidConvergence);
	}	
	MPI_Barrier(MPI_COMM_WORLD);
	destroyEnsemble(&myEnsemble);
	
	
	
	
	
	MPI_Finalize();
	return EXIT_SUCCESS;
	
	//-------------------------------------------------------------------//
	// Run multiseed optimization
	//-------------------------------------------------------------------//
	myEnsemble.ensembleSize = 1000;
	myEnsemble.ensembleT = 20000;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 50;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = C1 only, 3 = un-cooperative
	
	createEnsemble(&myEnsemble, framing);	
	NFE = 10000;
	// Create the rbf_sim problem, defining the number of decision variables,
	// objective and constraints.  The last argument, rbf_sim, references
	// the function that evaluates the rbf_sim problem.
	problem = BORG_Problem_create(numWeights, myEnsemble.ensemble[0]->numObjectives, 0, evaluateEnsemble, &myEnsemble);
	// Set the lower and upper bounds for each decision variable.
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -1, 1);
	}
	for (i=betaInd; i< betaInd + numNetHidden ; i++){
		BORG_Problem_set_bounds(problem, i, 0, 1);
	}
	// Set the epsilon values used by the Borg MOEA.  Epsilons define the
	// problem resolution, which controls how many Pareto optimal solutions
	// are generated and how far apart they are spaced. Set num function evals
	for (i = 0; i < myEnsemble.ensemble[0]->numObjectives; i++)
	{
		BORG_Problem_set_epsilon(problem, i, epsilons[i]);
	}
	
	
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/10000);
	
	// Run the Borg MOEA across multiple seeds
	nseeds = 20;
	for (i = 0;i < nseeds;i++)
	{
		char runtime[512];
		sprintf(runtime,"borg_multiseed_framing%d_simType%d_ensembleSize%d_T%d_NFE%d_epsilons%s_seed%d.runtime",framing,myEnsemble.simulationType,myEnsemble.ensembleSize,myEnsemble.ensembleT,NFE,epstrIni,i);
		BORG_Algorithm_output_runtime(runtime); //create runtime (i.e. evolution of performance and weights over time
		BORG_Random_seed(worldSize*i*(worldRank+1));
		result = BORG_Algorithm_ms_run(problem);
		if (result != NULL)
		{
			BORG_Archive_destroy(result);	
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	// Free any allocated memory and ensure processors are caught up
	BORG_Problem_destroy(problem);
	destroyEnsemble(&myEnsemble);
	//-------------------------------------------------------------------//
	
	
	
	//-------------------------------------------------------------------//
	// Check multiple epsilons
	//-------------------------------------------------------------------//
	int nepsilon = 20;
	// double epsilonUB[5] = {1000, 5000, 100000, 100000, 10};//Base epsilons to multiply up from
	// double epsilonLB[5] = {0.1, 0.5, 10, 10, 0.1};//Base epsilons to multiply up from
	myEnsemble.ensembleSize = 1000;
	myEnsemble.ensembleT = 20000;
	myEnsemble.baseSeed = 1;
	myEnsemble.percentile = 50;
	myEnsemble.simulationType = simType; //0 = cooperative, 1 = c2 only, 2 = C1 only, 4 = un-cooperative
	// myEnsemble.includeSimulationData = 0; //This can be = 0 or 24 - you must not include simulation data if running the ensemble more than once
	
	createEnsemble(&myEnsemble, framing);
		
	NFE = 10000;
	// Create the rbf_sim problem, defining the number of decision variables,
	// objective and constraints.  The last argument, rbf_sim, references
	// the function that evaluates the rbf_sim problem.
	problem = BORG_Problem_create(numWeights, myEnsemble.ensemble[0]->numObjectives, 0, evaluateEnsemble, &myEnsemble);
	// Set the lower and upper bounds for each decision variable.
	for (i=0; i<numWeights; i++) {
		BORG_Problem_set_bounds(problem, i, -1, 1);
	}
	for (i=betaInd; i< betaInd + numNetHidden ; i++){
		BORG_Problem_set_bounds(problem, i, 0, 1);
	}
	BORG_Algorithm_ms_max_evaluations(NFE);
	BORG_Algorithm_output_frequency(NFE/100);
	
	// Run the Borg MOEA across multiple epsilons
	
	for (i = 0;i < nepsilon;i++)
	{
		char epstr[100];
		sprintf(epstr,"epsilons");
		double *epsilonsRand =calloc(myEnsemble.ensemble[0]->numObjectives,sizeof(double));
		srand(i);
		for (j = 0; j < myEnsemble.ensemble[0]->numObjectives; j++)
		{
			epsilonsRand[j] = exp(log(epsilonLB[j])+log(epsilonUB[j]/epsilonLB[j])*(double)rand()/(double)RAND_MAX);
			sprintf(epstr + strlen(epstr),"_%0.2f",epsilonsRand[j]);
			BORG_Problem_set_epsilon(problem, j, epsilonsRand[j]);
		}
		char runtime[512];
		//no need to add seed here - no seed is used - will also need to change MATLAB plotting code
		sprintf(runtime,"epislon_test_framing%d_simType%d_ensembleSize%d_T%d_NFE%d_%s.runtime",framing,myEnsemble.simulationType,myEnsemble.ensembleSize,myEnsemble.ensembleT,NFE,epstr);
		BORG_Algorithm_output_runtime(runtime); //create runtime (i.e. evolution of performance and weights over time
		result = BORG_Algorithm_ms_run(problem);
		if (result != NULL)
		{
			BORG_Archive_destroy(result);	
		}
		MPI_Barrier(MPI_COMM_WORLD);
		free(epsilonsRand);
	}
	// Free any allocated memory
	
	BORG_Problem_destroy(problem);
	destroyEnsemble(&myEnsemble);
	
	
}