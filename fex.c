#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
// #include "qr_solve.h"
#include "fex.h"


double doumin(double *vec,int len)
{
	//Minimum value in a vector (vec) of doubles (length = len)
	double minno = DBL_MAX;
	int i;
	for (i=0;i<len;i++)
	{
		if (minno > vec[i])
		{
			minno = vec[i];
		}
	}
	return minno;
}
double doumax(double *vec,int len)
{
	//Maximum value in a vector (vec) of doubles (length = len)
	int minno = DBL_MIN;
	int i;
	for (i=0;i<len;i++)
	{
		if (minno < vec[i])
		{
			minno = vec[i];
		}
	}
	return minno;
}
void changeFexSimType(struct fexProblem *myFex, int simulationType)
{
	//Changes a fexProblem's simulation type (0 = cooperative, 1 = c2 only, 2 = c1 only, 3 = uncooperative (weights for the two policies are joined into one weights vector with c2 policy weights first then c1 policy weights)
	myFex->simulationType = simulationType;
	if (simulationType == 0 || simulationType == 3)
	{
		myFex->numObjectives=4;//do better
		
		myFex->objectiveSpecifier[0] = 0; //C2 reliability
		myFex->objectiveSpecifier[1] = 0; //C1 reliability
		myFex->objectiveSpecifier[2] = 1; //C2 linear def
		myFex->objectiveSpecifier[3] = 0; //C2 sqr def
		myFex->objectiveSpecifier[4] = 1; //C1 linear def
		myFex->objectiveSpecifier[5] = 0; //C1 sqr def
		myFex->objectiveSpecifier[6] = 0; //net£NPV
		myFex->objectiveSpecifier[7] = 0; //C2£NPV
		myFex->objectiveSpecifier[8] = 0; //C1£NPV
		myFex->objectiveSpecifier[9] = 0; //C1 lin gs
		myFex->objectiveSpecifier[10] = 0; //C1 sqr gs
		myFex->objectiveSpecifier[11] = 0; //net£
		myFex->objectiveSpecifier[12] = 1; //C2£
		myFex->objectiveSpecifier[13] = 1; //C1£

	}
	else if (simulationType == 1)
	{
		myFex->numObjectives=2;
		
		myFex->objectiveSpecifier[0] = 0; 
		myFex->objectiveSpecifier[1] = 0; 
		myFex->objectiveSpecifier[2] = 1; 
		myFex->objectiveSpecifier[3] = 0; 
		myFex->objectiveSpecifier[4] = 0; 
		myFex->objectiveSpecifier[5] = 0; 
		myFex->objectiveSpecifier[6] = 0; 
		myFex->objectiveSpecifier[7] = 0; 
		myFex->objectiveSpecifier[8] = 0; 
		myFex->objectiveSpecifier[9] = 0; 
		myFex->objectiveSpecifier[10] = 0;
		myFex->objectiveSpecifier[11] = 0;
		myFex->objectiveSpecifier[12] = 1;
		myFex->objectiveSpecifier[13] = 0;
		

	}
	else if (simulationType == 2)
	{
		myFex->numObjectives=2;
		
		myFex->objectiveSpecifier[0] = 0; 
		myFex->objectiveSpecifier[1] = 0;
		myFex->objectiveSpecifier[2] = 0;
		myFex->objectiveSpecifier[3] = 0;
		myFex->objectiveSpecifier[4] = 1;
		myFex->objectiveSpecifier[5] = 0;
		myFex->objectiveSpecifier[6] = 0;
		myFex->objectiveSpecifier[7] = 0;
		myFex->objectiveSpecifier[8] = 0;
		myFex->objectiveSpecifier[9] = 0;
		myFex->objectiveSpecifier[10] = 0;
		myFex->objectiveSpecifier[11] = 0;
		myFex->objectiveSpecifier[12] = 0;
		myFex->objectiveSpecifier[13] = 1;
	
	}
	
}
int sumInt(int *vec, int T)
{
	//Sum of vector (vec) of ints (length = T)
	int i;
	int sum = 0;
	for (i=0;i<T;i++)
	{
		sum = sum + vec[i];
	}
	return sum;
}

void createFex(struct fexProblem *myFex, struct fexEnsemble *myEnsemble, int T, int baseSeed, int simulationType)
{
	//Create a fex problem - in this function you can alter the parameters for simulation, choose what variables are used in a policy and specify the architecture of the policy network.
	//Memory is allocated here for forcing variables (although these forcing variables are populated in fex.c/populateSeededVariables)
	//fexEnsemble is only required here since it contains historical data - nothing else from it is used
	//NO NEED FOR BASESEED HERE
	int i;
	/// -------------- SYSTEM PARAMETERS --------------- ///
	myFex->linkCapS2D2=36.0; // daily link cap in Ml
	myFex->linkCapS1D2=45.0; // daily link cap in Ml
	myFex->linkCapR1S1=150; // daily link cap in Ml
	myFex->linkCapS1R1=600;// daily link cap in Ml
	myFex->linkCapBorehole=50;// daily link cap in Ml
	myFex->linkAnnualCapS1D2 = 14917.26;// annual link cap in Ml
	myFex->linkAnnualCapR1S1 = 13633; // annual link cap in Ml
	myFex->linkAnnualCapS2D2 = 10000;// annual link cap in Ml
	myFex->linkMonthlyCapS2D2 = 1120;// monthly link cap in Ml
	myFex->linkAnnualCapR1S1 = 12585;// annual link cap in Ml - note that this limit is ONLY on water that is later consumed by c1
	
	myFex->storCapS2=5662.21;
	myFex->storCapS1=21365;
	
	myFex->discountRate = 0.03/365;
	
	myFex->deadStorS2=0;
	myFex->deadStorS1=0;
	
	myFex->startDoyNoR1S1 = 91; //1st april -- no R1S1 between these dates
	myFex->endDoyNoR1S1 = 304; //31st october
	
	myFex->startDoyMaxR1D2 = 305; //In simType 2, this is when uS1D2 = max
	myFex->endDoyMaxR1D2 = 90; //the rest of the annual license is evenly spread outside these dates
	
	myFex->storElevationS2[0] = 1.564; //Elev = param[0]*stor^param[1]
	myFex->storElevationS2[1] = 0.3418;
	myFex->elevationAreaS2[0] = 3200;
	myFex->elevationAreaS2[1] = 1.5;
	
	myFex->storElevationS1[0] = 3.849;
	myFex->storElevationS1[1] = 0.2569; 
	myFex->elevationAreaS1[0] = 4250;
	myFex->elevationAreaS1[1] = 1.5;
	
	
	myFex->s1CompensationFixed = 1;
	myFex->s2CompensationFixed = 5;
	
	myFex->R1Thresh = 1.16*86.4;
	myFex->R2Thresh = 3.158*86.4;
		
	myFex->fexDemandD1Fixed = 200;
	
	myFex->r2r1ConversionFactor = 1/0.746; //(Flow at R2) = r2r1ConversionFactor * (Flow at R1), C2 use this conversion factor
	/// ------------------------------------------------ ///
	
	
	/// ------------- SIMULATION PARAMETERS ------------ ///
	myFex->fexInitialStorS2=5000; //Note - this may be changed in fex.c/populateSeededVariables
	myFex->fexInitialStorS1=18000; //Note - this may be changed in fex.c/populateSeededVariables
	myFex->simulationType = simulationType; //don't change this here
	myFex->printDynamics = 0; //don't change this here
	// myFex->includetimeObjectives = 0;
	myFex->fexT = T;
	
	/// ------------------------------------------------ ///
	
	
	/// ------------- EVALUATION PARAMETERS ------------ ///
	myFex->storageReliabilityThresholdS2=2000;
	myFex->storageReliabilityThresholdS1=15000;
	
	myFex->numObjectiveOptions = 14;
	myFex->objectiveSpecifier = calloc(myFex->numObjectiveOptions,sizeof(int));
	// myFex->objectiveSpecifierIndex = calloc(myFex->numObjectiveOptions,sizeof(int));
	changeFexSimType(myFex,simulationType);
	myFex->linkCostS2D2 = 0;
	myFex->linkCostS1D2 = 18;
	myFex->linkCostR1S1 = 50;
	myFex->linkCostS1R1 = 0;
	myFex->linkCostGroundBorehole = 40;
	
	/// ------------------------------------------------ ///
	
	/// --------------- POLICY PARAMETERS -------------- ///
	
	//Describe what to normalize inputs between
	myFex->inflowNormMaxMinS2[0] = 200;
	myFex->inflowNormMaxMinS2[1] = 0;
	myFex->inflowNormMaxMinS1[0] = 500;
	myFex->inflowNormMaxMinS1[1] = 0;
	myFex->inflowNormMaxMinR1[0] = 4000;
	myFex->inflowNormMaxMinR1[1] = 0;
	myFex->inflowNormMaxMinR2[0] = 4000;
	myFex->inflowNormMaxMinR2[1] = 0;
	
	myFex->minS2D2 = 12;
	myFex->minS1D2 = 15;
	
	myFex->petNormMaxMinS1[0] = 4.5;
	myFex->petNormMaxMinS1[1] = 0;

	myFex->demandNormMaxMinD2[0] = 73;
	myFex->demandNormMaxMinD2[1] = 36;
	
	//Specify what policy inputs are...
	myFex->numInputOptions = 16;
	//...for cooperative net
	myFex->cooperativeNet = calloc(1,sizeof(struct netOptions));
	myFex->cooperativeNet->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->cooperativeNet->policyInputSpecifier[0] = 1; //S2 storage
	myFex->cooperativeNet->policyInputSpecifier[1] = 1; //S2 inflow
	myFex->cooperativeNet->policyInputSpecifier[2] = 1; //s1 storage
	myFex->cooperativeNet->policyInputSpecifier[3] = 1; //s1 inflow
	myFex->cooperativeNet->policyInputSpecifier[4] = 1; //d2 demand
	myFex->cooperativeNet->policyInputSpecifier[5] = 0; //pet
	myFex->cooperativeNet->policyInputSpecifier[6] = 0; //CM link cap
	myFex->cooperativeNet->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->cooperativeNet->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//
	myFex->cooperativeNet->policyInputSpecifier[9] = 0; //EW link cap
	myFex->cooperativeNet->policyInputSpecifier[10] = 1; //R1 flow
	myFex->cooperativeNet->policyInputSpecifier[11] = 0; //WE link cap
	myFex->cooperativeNet->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->cooperativeNet->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->cooperativeNet->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->cooperativeNet->policyInputSpecifier[15] = 0; //R2 flow
	myFex->cooperativeNet->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	int l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->cooperativeNet->policyInputSpecifier[i] == 1)
		{
			myFex->cooperativeNet->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}

	//...for c2 net
	myFex->c2Net = calloc(1,sizeof(struct netOptions));
	myFex->c2Net->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->c2Net->policyInputSpecifier[0] = 1; //s2 storage
	myFex->c2Net->policyInputSpecifier[1] = 1; //s2 inflow
	myFex->c2Net->policyInputSpecifier[2] = 0; //s1 storage
	myFex->c2Net->policyInputSpecifier[3] = 0; //s1 inflow
	myFex->c2Net->policyInputSpecifier[4] = 1; //d2 demand
	myFex->c2Net->policyInputSpecifier[5] = 0; //pet
	myFex->c2Net->policyInputSpecifier[6] = 0; //CM link cap
	myFex->c2Net->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->c2Net->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//
	myFex->c2Net->policyInputSpecifier[9] = 0; //EW link cap
	myFex->c2Net->policyInputSpecifier[10] = 0; //R1 flow
	myFex->c2Net->policyInputSpecifier[11] = 0; //WE link cap
	myFex->c2Net->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->c2Net->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->c2Net->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->c2Net->policyInputSpecifier[15] = 0; //R2 flow
	myFex->c2Net->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->c2Net->policyInputSpecifier[i] == 1)
		{
			myFex->c2Net->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}

	//...for c1 net
	myFex->c1Net = calloc(1,sizeof(struct netOptions));
	myFex->c1Net->policyInputSpecifier = calloc(myFex->numInputOptions,sizeof(int));
	myFex->c1Net->policyInputSpecifier[0] = 0; //s2 storage
	myFex->c1Net->policyInputSpecifier[1] = 0; //s2 inflow
	myFex->c1Net->policyInputSpecifier[2] = 1; //s1 storage
	myFex->c1Net->policyInputSpecifier[3] = 1; //s1 inflow
	myFex->c1Net->policyInputSpecifier[4] = 0; //d2 demand
	myFex->c1Net->policyInputSpecifier[5] = 0; //pet
	myFex->c1Net->policyInputSpecifier[6] = 0; //CM link cap
	myFex->c1Net->policyInputSpecifier[7] = 0; //WM link cap
	//------NOTE:THIS IS CURRENTLY NON-FUNCTIONAL SINCE TIME VARYING c1 DEMANDS ARE NOT GIVEN, DO NOT SET THIS TO 1-----------//
	myFex->c1Net->policyInputSpecifier[8] = 0; //c1 demand
	//-----------------//	myFex->c1Net->policyInputSpecifier[9] = 0; //EW link cap
	myFex->c1Net->policyInputSpecifier[10] = 1; //R1 flow
	myFex->c1Net->policyInputSpecifier[11] = 0; //WE link cap
	myFex->c1Net->policyInputSpecifier[12] = 0; //GS supplement cap
	myFex->c1Net->policyInputSpecifier[13] = 0; //sin(doy*pi*2/365)
	myFex->c1Net->policyInputSpecifier[14] = 0; //cos(doy*pi*2/365)
	myFex->c1Net->policyInputSpecifier[15] = 0; //R2 flow
	myFex->c1Net->policyInputSpecifierIndex = calloc(myFex->numInputOptions,sizeof(int));
	l = 0;
	for (i=0;i<myFex->numInputOptions;i++)
	{
		
		if (myFex->c1Net->policyInputSpecifier[i] == 1)
		{
			myFex->c1Net->policyInputSpecifierIndex[i] = l;
			l = l + 1;
		}
	}
	

	//Handy numbers for...
	//...cooperative net
	myFex->cooperativeNet->numNetOutputs = 4;
	myFex->cooperativeNet->numNetHidden = 2;
	myFex->cooperativeNet->numNetInputs = sumInt(myFex->cooperativeNet->policyInputSpecifier,myFex->numInputOptions);
	myFex->cooperativeNet->numWeights = (myFex->cooperativeNet->numNetInputs + 1)*(myFex->cooperativeNet->numNetHidden) + (1 + myFex->cooperativeNet->numNetHidden)*(myFex->cooperativeNet->numNetOutputs);
	//...c2 net
	myFex->c2Net->numNetOutputs = 2;
	myFex->c2Net->numNetHidden = 2;
	myFex->c2Net->numNetInputs = sumInt(myFex->c2Net->policyInputSpecifier,myFex->numInputOptions);
	myFex->c2Net->numWeights = (myFex->c2Net->numNetInputs + 1)*(myFex->c2Net->numNetHidden) + (1 + myFex->c2Net->numNetHidden)*(myFex->c2Net->numNetOutputs);
	//...c1 net
	myFex->c1Net->numNetOutputs = 2;
	myFex->c1Net->numNetHidden = 2;
	myFex->c1Net->numNetInputs = sumInt(myFex->c1Net->policyInputSpecifier,myFex->numInputOptions);
	myFex->c1Net->numWeights = (myFex->c1Net->numNetInputs + 1)*(myFex->c1Net->numNetHidden) + (1 + myFex->c1Net->numNetHidden)*(myFex->c1Net->numNetOutputs);
	
	/// ------------------------------------------------ ///
	
	/// ---------------- ALLOCATE MEMORY --------------- /// -- would be good to also have a reallocate function (i.e. frees memory, then allocates to a new fexT
	myFex->warmup = 10000;
	myFex->fexFlowsS2 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsS1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsR1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexFlowsR2 = calloc(myFex->fexT,sizeof(double));
	myFex->s1FishReleases = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS2D2TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS1D2TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapR1S1TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapS1R1TimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->linkCapBoreholeTimeVarying = calloc(myFex->fexT,sizeof(double));
	myFex->fexPets= calloc(myFex->fexT,sizeof(double));
	// myFex->fexDemandsc1 = calloc(myFex->fexT,sizeof(double));
	myFex->fexDemandsD2 = calloc(myFex->fexT,sizeof(double));
	myFex->fexDoys = calloc(myFex->fexT,sizeof(int));
	/// ------------------------------------------------ ///
}

void matXvec(double *mat, double *vec, int M, int N, double *product)
{
	//Performs matrix multipilcation of a [MxN] matrix (flattened such that mat[i][j] is equivalent to mat[j*M + i]) with a [N] vector (vec) and stores result in pointer product (a [M] vector)
	int i,j;
	for (i=0;i<M;i++)
	{
		product[i] = 0;
		for (j=0;j<N;j++)
		{
			product[i] = product[i] + mat[j*M + i]*vec[j];
		}
	}
}

int poisRnd(int lambda)
{
	int n = 0;
	double limit = exp(-lambda);
	double x = (double)rand()/(double)RAND_MAX;
	while (x > limit)
	{
		n++;
		x *= (double)rand()/(double)RAND_MAX;
	}
	return n;
}


void generateBreaks(double *timeVaryingCapacity, double breakLambdaProbability[2], int syntheticT, double capacityHard,int seed)
{
	//Breaks re-occur with poisson distribution (lambda = breakLambdaProbability[1])
	//Breaks last with duration poisson(lambda = breakLambdaProbability[0])
	//Currently a break is equivalent to setting link capacity = 0
	int condition = 0;
	srand(seed);
	int i;
	int timeUntilNextBreak = poisRnd(breakLambdaProbability[1]);
	for (i=0;i < syntheticT;i++)
	{
		if (condition > 0)
		{
			condition = condition - 1;
			if (condition == 0)
			{
				timeUntilNextBreak = i + poisRnd(breakLambdaProbability[1]);
			}
		}
		if ((i - timeUntilNextBreak) == 0)
		{
			condition = poisRnd(breakLambdaProbability[0]);
		}
		
		
		if (condition == 0)
		{
			timeVaryingCapacity[i] = capacityHard;
		}
		else
		{
			timeVaryingCapacity[i] = 0;
		}
	}
}

double normRnd (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

void generateParma(double *synthetic, struct parma *myParma, double *innovation, int *doys, int T, int logOpt,int offset)
{
	int i, j;
	double *residual = malloc((T + offset) * sizeof(double));
	for (i = 0; i < offset; i++)
	{	
		residual[i] = 0;
	}
	double sigma = sqrt(myParma->varc);
	for (i = offset; i < (T + offset); i++)
	{
		double ar = 0;
		for (j = i - myParma->numAR; j < i; j++)
		{
			ar = ar + residual[j]*myParma->AR[ j - i + myParma->numAR ];
		}
		
		double ma = 0;
		for (j = i - myParma->numMA; j < i; j++)
		{
			ma = ma + innovation[j]*myParma->MA[ j - i + myParma->numMA ]*sigma;
		}
		residual[i] = myParma->cons + ar + ma + innovation[i]*sigma;
	}
	for (i = 0; i < T; i++)
	{
		double periodic = myParma->thetas[0];
		for (j = 0; j < myParma->numHarmonics; j++)
		{
			periodic = periodic + myParma->thetas[(j+1)*2-1]*sin(myParma->harmonics[j]*M_PI*doys[i]/365);
			periodic = periodic + myParma->thetas[(j+1)*2]*cos(myParma->harmonics[j]*M_PI*doys[i]/365);
		}
		if (logOpt == 0)
		{
			synthetic[i] = fmax(periodic + residual[i + offset],0);
		}
		else
		{
			synthetic[i] = fmax(periodic*exp(residual[i + offset]),0);
		}
	}
	free(residual);
}
void populateSeededVariables(struct fexProblem *myFex, struct fexEnsemble *myEnsemble, int seed)
{
	//Populate a fexProblem's seeded variables (i.e. forcing).
	//myEnsemble is used since it contains historic data, synthetic parameters, break parameters and fisheries parameters
	int i,j,k,doy;
	
	srand(seed);
	int startDoy = (rand()+1)%365;
	if (startDoy == 0)
	{
		startDoy = 365;
	}
	for (i=0;i<(myFex->fexT);i++)
	{
		doy = (i + startDoy)%365;
		if (doy == 0)
		{
			doy = 365;
		}
		myFex->fexDoys[i] = doy; //set days of year, leap days don't exist
	}
	int numSeededVariables = 11;
	int *seeds = calloc(numSeededVariables,sizeof(int));
	// srand(seed);
	//endedit
	
	
	for (i = 0; i < numSeededVariables; i++)
	{
		seeds[i] = rand()+i; //create seeds
	}
	
	/// ------- CREATE AUTOCORRELATED VARIABLES --------- ///
	//Generate cross correlated innovation
	srand(seeds[0]);
	double **randomNormal = malloc(myEnsemble->NUM_AUTOCORRELATED_PROCESSES * sizeof(double*));
	double **randomCorrNormal = malloc(myEnsemble->NUM_AUTOCORRELATED_PROCESSES * sizeof(double*));
	for (i = 0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		randomNormal[i] = malloc( (myFex->fexT + myEnsemble->offset) * sizeof(double));
		randomCorrNormal[i] = malloc( (myFex->fexT + myEnsemble->offset) * sizeof(double));
		for (j = 0; j < (myFex->fexT + myEnsemble->offset); j++)
		{
			randomNormal[i][j] = normRnd(0,1);
		}
	}
	//matrix[cols][rows]
	//randomNormal[numAuto][T+offset]
	//chol[numAuto][numAuto]
	for (j = 0; j < (myFex->fexT + myEnsemble->offset); j++)
	{
		for (i = 0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
		{
			randomCorrNormal[i][j] = 0;
			for (k = 0; k < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; k++)
			{
				randomCorrNormal[i][j] = randomCorrNormal[i][j] + randomNormal[k][j]*myEnsemble->cholesky[k][i];
			}
		}
	}
	
	//generate synthetics
	generateParma(myFex->fexFlowsS2, &myEnsemble->S2Parma, randomCorrNormal[0], myFex->fexDoys, myFex->fexT, 1, myEnsemble->offset);
	generateParma(myFex->fexFlowsS1, &myEnsemble->S1Parma, randomCorrNormal[1], myFex->fexDoys, myFex->fexT, 1, myEnsemble->offset);
	generateParma(myFex->fexFlowsR1, &myEnsemble->R1Parma, randomCorrNormal[2], myFex->fexDoys, myFex->fexT, 1, myEnsemble->offset);
	generateParma(myFex->fexFlowsR2, &myEnsemble->R2Parma, randomCorrNormal[3], myFex->fexDoys, myFex->fexT, 1, myEnsemble->offset);
	generateParma(myFex->fexPets, &myEnsemble->petParma, randomCorrNormal[4], myFex->fexDoys, myFex->fexT, 0, myEnsemble->offset);
	generateParma(myFex->fexDemandsD2, &myEnsemble->D2Parma, randomCorrNormal[5], myFex->fexDoys, myFex->fexT, 0, myEnsemble->offset);
	
	for (i = 0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		free(randomNormal[i]);
		free(randomCorrNormal[i]);
	}
	free(randomNormal);
	free(randomCorrNormal);
	
	/// ---------------- GENERATE BREAKS --------------- ///
	//breaks are applied as a temporary drop in link capacity
	generateBreaks(myFex->linkCapS2D2TimeVarying,myEnsemble->breakLambdaProbabilityS2D2, myFex->fexT,myFex->linkCapS2D2,seeds[1]);
	
	generateBreaks(myFex->linkCapS1D2TimeVarying,myEnsemble->breakLambdaProbabilityS1D2, myFex->fexT,myFex->linkCapS1D2,seeds[2]);
	
	generateBreaks(myFex->linkCapR1S1TimeVarying,myEnsemble->breakLambdaProbabilityR1S1, myFex->fexT,myFex->linkCapR1S1,seeds[3]);
	
	generateBreaks(myFex->linkCapS1R1TimeVarying,myEnsemble->breakLambdaProbabilityS1R1, myFex->fexT,myFex->linkCapS1R1,seeds[4]);
	
	generateBreaks(myFex->linkCapBoreholeTimeVarying,myEnsemble->breakLambdaProbabilityBorehole, myFex->fexT,myFex->linkCapBorehole,seeds[5]);

	/// ------------------------------------------------ ///
	
	
	/// --------------- INITIAL STORAGES --------------- ///
	
	srand(seeds[9]);
	myFex->fexInitialStorS2=((double)rand()/(double)RAND_MAX)*(myFex->storCapS2 - myFex->deadStorS2) + myFex->deadStorS2;
	myFex->fexInitialStorS1=((double)rand()/(double)RAND_MAX)*(myFex->storCapS1 - myFex->deadStorS1) + myFex->deadStorS1;
	/// ------------------------------------------------ ///
	
	/// ------------- S1 FISHERIES ------------- ///
	
	
	srand(seeds[10]);
	i=0;
	while (i <myFex->fexT)
	{
		myFex->s1FishReleases[i] = 0;
		if (myFex->fexDoys[i] == myEnsemble->startAllowableFishReleaseDoy)
		{
			int lengthOfReleases = rand()%(myEnsemble->maxFishReleaseLength - myEnsemble->minFishReleaseLength + 1);
			lengthOfReleases = lengthOfReleases + myEnsemble->minFishReleaseLength;
			int startOfReleases = rand()%(myEnsemble->endAllowableFishReleaseDoy - myEnsemble->startAllowableFishReleaseDoy - lengthOfReleases + 2);
			j = i;
			while (i < fmin(j+startOfReleases,myFex->fexT))
			{
				myFex->s1FishReleases[i] = 0;
				i++;
			}
			j = i;
			while (i < fmin(j+lengthOfReleases,myFex->fexT))
			{
				myFex->s1FishReleases[i] = myEnsemble->s1FisheriesAnnual/lengthOfReleases;
				i++;
			}
			if (i < myFex->fexT)
			{
				myFex->s1FishReleases[i] = 0;
			}
		}
		i++;
	}
	
	/// ------------------------------------------------ ///
	
	free(seeds);
}
FILE * readParm(FILE *f, struct parma *myParma, int *offset)
{
	int i;
	char line[256];
	fgets(line,sizeof(line),f);
	myParma->numHarmonics = atoi(line);
	
	
	myParma->harmonics = malloc( myParma->numHarmonics * sizeof(double));
	for (i = 0; i < myParma->numHarmonics; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->harmonics[i] = atof(line);
		
	}
	
	myParma->thetas = malloc( (2*myParma->numHarmonics + 1) * sizeof(double));
	for (i = 0; i < (2*myParma->numHarmonics + 1); i++)
	{
		fgets(line,sizeof(line),f);
		myParma->thetas[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->numAR = atoi(line);
	*offset = fmax(myParma->numAR,*offset);
	
	myParma->AR = malloc( myParma->numAR * sizeof(double));
	for (i = 0; i < myParma->numAR ; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->AR[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->numMA = atoi(line);
	*offset = fmax(myParma->numMA,*offset);
	
	myParma->MA = malloc( myParma->numMA * sizeof(double));
	for (i = 0; i < myParma->numMA ; i++)
	{
		fgets(line,sizeof(line),f);
		myParma->MA[i] = atof(line);
	}
	
	fgets(line,sizeof(line),f);
	myParma->cons = atof(line);
	
	fgets(line,sizeof(line),f);
	myParma->varc = atof(line);
	fgets(line,sizeof(line),f);
	if (line[0] != '#')
	{
		printf("warning: formatting problem in parma file\n");
	}
	
	return f;
}
void createEnsemble(struct fexEnsemble *myEnsemble, int framing)
{
	//AN ENSEMBLE (fexEnsemble) CONTAINS EVERYTHING NOT USED IN THE SIMULATION AND POINTERS (fexEnsemble.ensemble) TO THE PROBLEMS
	//A PROBLEM (fexProblem) CONTAINS ONLY WHAT IS NEEDED FOR THE SIMULATION
	
	int i, j;
	myEnsemble->NUM_AUTOCORRELATED_PROCESSES = 6;
	

	char parmaFid[256];
	char choleskyFid[256];
	if (framing == 1 || framing == 2 || framing == 5 || framing == 6)
	{
		sprintf(parmaFid,"parma10.arma");
		sprintf(choleskyFid,"cholesky10.arma");
	}
	else
	{
		sprintf(parmaFid,"parma20.arma");
		sprintf(choleskyFid,"cholesky20.arma");
	}
	
	/// --------------- READ CHOLESKY ------------------ ///
	myEnsemble->cholesky = malloc( myEnsemble->NUM_AUTOCORRELATED_PROCESSES * (sizeof(double*)) );
	FILE *fChol = fopen(choleskyFid,"r");
	for (i=0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		myEnsemble->cholesky[i] = malloc( myEnsemble->NUM_AUTOCORRELATED_PROCESSES * (sizeof(double)) );
		for (j=0; j < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; j++)
		{
			if (!fscanf(fChol, "%lf", &myEnsemble->cholesky[i][j]))
			break;
		}
	}
	fclose(fChol);
	/// ----------- READ AND ALLOCATE PARMA ------------ ///
	
	FILE *fParm = fopen(parmaFid,"r");
	myEnsemble->offset = 0;
	fParm = readParm(fParm,&myEnsemble->S2Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->S1Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->R1Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->R2Parma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->petParma,&myEnsemble->offset);
	fParm = readParm(fParm,&myEnsemble->D2Parma,&myEnsemble->offset);
	fclose(fParm);
	/// --------------- BREAK PARAMETERS --------------- ///
	if (framing == 1 || framing == 5)
	{
		myEnsemble->breakLambdaProbabilityS2D2[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS2D2[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1D2[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1D2[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityR1S1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityR1S1[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1R1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1R1[1] = 1000000000000; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityBorehole[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityBorehole[1] = 1000000000000; //lambda for poisson duration between break
	}
	else
	{
		myEnsemble->breakLambdaProbabilityS2D2[0] = 5; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS2D2[1] = 800; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1D2[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1D2[1] = 300; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityR1S1[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityR1S1[1] = 300; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityS1R1[0] = 0; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityS1R1[1] = 0; //lambda for poisson duration between break
		myEnsemble->breakLambdaProbabilityBorehole[0] = 3; //lambda for poisson duration of break
		myEnsemble->breakLambdaProbabilityBorehole[1] = 500; //lambda for poisson duration between break
	}
	/// ------------------------------------------------ ///

	
	/// ------------- S1 FISHERIES ------------- ///
	
	//s1FisheriesAnnual will be uniformly released over a random length between min/maxDailyFishRelease during start/endAllowableFishReleaseDoy
	myEnsemble->s1FisheriesAnnual = 900; //Ml/y
	if (framing == 4 || framing == 8)
	{
		myEnsemble->startAllowableFishReleaseDoy = 60; //i.e. March 1
		myEnsemble->endAllowableFishReleaseDoy = 365; //i.e. Dec 31
	}
	else
	{
		myEnsemble->startAllowableFishReleaseDoy = 244; //i.e. Sept 1
		myEnsemble->endAllowableFishReleaseDoy = 273; //i.e. Sept 31	
	}
	myEnsemble->maxFishReleaseLength = 14; //Entire annual allowance must be released over a maximum of 14 consecutive days
	myEnsemble->maxDailyFishRelease = 360; //Ml/d - note- that since this isn't divisible by 900, the actual max is 300
	myEnsemble->minFishReleaseLength = (int)ceil(myEnsemble->s1FisheriesAnnual/myEnsemble->maxDailyFishRelease);
	
	/// ------------------------------------------------ ///
	
	/// ---------------- CREATE ENSEMBLE --------------- ///
	myEnsemble->ensemble = calloc(myEnsemble->ensembleSize,sizeof(struct fexProblem*));
	srand(myEnsemble->baseSeed);
	for (i=0;i<myEnsemble->ensembleSize;i++)
	{
		myEnsemble->ensemble[i] = calloc(1,sizeof(struct fexProblem));
		createFex(myEnsemble->ensemble[i], myEnsemble, myEnsemble->ensembleT, myEnsemble->baseSeed,myEnsemble->simulationType);
		int seed = rand()+i;
		populateSeededVariables(myEnsemble->ensemble[i],myEnsemble,seed); //what it says on the tin
	}
	/// ------------------------------------------------ ///

	
	/// ---------------- CHANGE FRAMING --------------- ///
	myEnsemble->framing = framing;

}

void destroyFex(struct fexProblem *myFex)
{
	
	free(myFex->fexFlowsS2);
	free(myFex->fexFlowsS1);
	free(myFex->fexFlowsR1);
	free(myFex->fexFlowsR2);
	
	free(myFex->s1FishReleases);
	
	free(myFex->fexDoys);
	
	free(myFex->linkCapS2D2TimeVarying);
	free(myFex->linkCapS1D2TimeVarying);
	free(myFex->linkCapR1S1TimeVarying);
	free(myFex->linkCapS1R1TimeVarying);
	free(myFex->linkCapBoreholeTimeVarying);
	
	free(myFex->fexPets);	
	free(myFex->fexDemandsD2);
	// free(myFex->fexDemandsc1);
	
	free(myFex->objectiveSpecifier);
	
	free(myFex->cooperativeNet->policyInputSpecifier);
	free(myFex->cooperativeNet->policyInputSpecifierIndex);
	free(myFex->cooperativeNet);	
	free(myFex->c2Net->policyInputSpecifier);
	free(myFex->c2Net->policyInputSpecifierIndex);
	free(myFex->c2Net);
	free(myFex->c1Net->policyInputSpecifier);
	free(myFex->c1Net->policyInputSpecifierIndex);
	free(myFex->c1Net);
}

void destroyParma(struct parma *myParma)
{
	free(myParma->harmonics);
	free(myParma->thetas);
	free(myParma->AR);
	free(myParma->MA);
}
void destroyEnsemble(struct fexEnsemble *myEnsemble)
{
	//Free fexEnsemble (and fex problems within fexEnsemble.ensemble)
	int i;
	for (i = 0;i<myEnsemble->ensembleSize;i++)
	{
		destroyFex(myEnsemble->ensemble[i]);
		free(myEnsemble->ensemble[i]);
	}
	
	free(myEnsemble->ensemble);
	
	destroyParma(&myEnsemble->S2Parma);
	destroyParma(&myEnsemble->S1Parma);
	destroyParma(&myEnsemble->R2Parma);
	destroyParma(&myEnsemble->R1Parma);
	destroyParma(&myEnsemble->D2Parma);
	destroyParma(&myEnsemble->petParma);
	
	for (i = 0; i < myEnsemble->NUM_AUTOCORRELATED_PROCESSES; i++)
	{
		free(myEnsemble->cholesky[i]);
	}
	free(myEnsemble->cholesky);
}
void evaluateNet(double *outputs, double *inputs, int numNetInputs, int numNetOutputs, int numNetHidden, double *weights)
{
	//evaluate an RBF net with described architecture and inputs, storing result in outputs
	//TOTAL NUMBER OF WEIGHTS = (numNetInputs+1)*numNetHidden + numNetOutputs*(numNetHidden+1)
	//Implements rbfn according to http://mccormickml.com/2013/08/15/radial-basis-function-network-rbfn-tutorial/

	int i,j;
	
	double *q = calloc(numNetHidden,sizeof(double));
	
	for (i=0; i < numNetHidden; i++)
	{
		double sqrDiff = 0;
		for (j = 0; j < numNetInputs; j++)
		{
			sqrDiff = sqrDiff + (weights[i+j*numNetHidden] - inputs[j])*(weights[i+j*numNetHidden] - inputs[j]);
		}
		q[i] = exp( -weights[i+numNetHidden*numNetInputs]*sqrDiff );
	}
	
	// // May wish to normalize here
	// double sumQ = sumDouble(q,numNetHidden);
	// for (i = 0; i < numNetHidden; i++)
	// {
		// q[i] = q[i]/sumQ;
	// }
	
	for (i=0; i < numNetOutputs; i++)
	{
		double Y = weights[ numNetHidden*(numNetInputs +1) + i ];
		for (j=0; j < numNetHidden; j++)
		{
			Y = Y + q[j]*weights[ numNetHidden*(numNetInputs +1) + numNetOutputs*(j+1) + i ];
		}
		if (Y < -1)
		{
			outputs[i] = 0;
		}
		else if (Y < 1)
		{
			outputs[i] = (Y + 1)/2;
		}
		else
		{
			outputs[i] = 1;
		}
	}
	free(q);
}
double sumDouble(double *vec, int T)
{
	//Sum of vector (vec) of doubles (length = T)
	int i;
	double sum = 0;
	for (i=0;i<T;i++)
	{
		sum = sum + vec[i];
	}
	return sum;
}
void simFex(double* weights, double* objs, double* consts,struct fexProblem *myFex)
{
	// simulate a fexProblem for weights, storing objectives in objs
	// consts is currently unused
	
	/// ---------------- ALLOCATE MEMORY --------------- ///
	
	int i, j, t, l;
	double *storageS2 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //S2 stor
	double *spillS2 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //S2 spill
	double *storageS1 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //s1 stor
	double *spillS1 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //s1 spill
	double *supplyToD2 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //How much to supply to D2
	double *supplyToC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //How much to supply to D1
	double *volumeS2D2 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //s2 release
	double *volumeS1D2 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //s1 release
	double *volumeR1S1 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //river R1 pumped abstraction to s1
	double *volumeS1R1 = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //s1all release to R1 for DE2
	double *volumeBorhole = calloc(myFex->warmup + myFex->fexT,sizeof(double)); //c1 groundwater supplement
	
	double elevation; double area; double evap;
	// double pumpCostTotal = 0;
	// double pumpCostC2 = 0;
	// double pumpCostC1 = 0;
	// double NPVc2 = 0;
	// double NPVC1 = 0;
	// double deficitSqrD2 = 0;
	// double deficitLinD2 = 0;
	// double deficitSqrC1 = 0;
	// double deficitLinC1 = 0;
	// double gsSqrC1 = 0;
	// double gsLinC1 = 0;
	// double storageReliabilityS2 = 0;
	// double storageReliabilityS1 = 0;
	double S2D2AnnualTotal = 0;
	double S1D2AnnualTotal = 0;
	double R1S1AnnualTotal = 0;
	double S2D2MonthlyTotal = 0;
	double S1R1AnnualTotal = 0;
	
	double *pumpCostTotal = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *pumpCostC2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *pumpCostC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVc2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *NPVtotal = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitSqrD2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitLinD2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitSqrC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *deficitLinC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *gsSqrC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *gsLinC1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *storageReliabilityS2 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	double *storageReliabilityS1 = calloc(myFex->warmup + myFex->fexT,sizeof(double));
	
	const double DOUBLE_PRECISION = pow(2,-52);
	const int DAYS_IN_YEAR = 365;
	const int DAYS_IN_MONTH = 30;
	const double MM_M2_TO_ML = pow(10,-6);
	// NOTE: outputting objectives over time through use of 'timeObjectives' and 'includetimeObjectives' is currently untested and commented out
	// timeObjectives = calloc(myFex->warmup + myFex->numObjectives,sizeof(double)); // NOTE: if includetimeObjectives == 1, then this data and the associated objectives will not be freed (not even in destroyFex)
	
	/// ------------------------------------------------ ///
	
	
	/// ---------------- SETUP PRINTING ---------------- ///
	FILE *fid;
	if (myFex->printDynamics == 1)
	{
		fid = fopen(myFex->dynamicFid,"w");
		fprintf(fid,"//Sc Sw Ic Iw Iex Ith Dm Ds ucm uwm uew uwe gs wc ww doy evapC evapW s2Env s1Env");
		// for (i = 0; i < myFex->numObjectives; i++)
		// {
			// fprintf(fid, " obj%d", i);
		// }
		fprintf(fid,"\n");
	}
	/// ------------------------------------------------ ///
	
	
	/// ---------------- SETUP NETWORKS ---------------- ///
	
	double *c2Weights;
	double *C1Weights;
	int numNetInputs;
	int numNetOutputs;
	int numNetHidden;
	int numWeights;
	int *policyInputSpecifier;
	int *policyInputSpecifierIndex;
	switch (myFex->simulationType)
	{
		case 0:
			numNetInputs = myFex->cooperativeNet->numNetInputs;
			numNetOutputs = myFex->cooperativeNet->numNetOutputs;
			numNetHidden = myFex->cooperativeNet->numNetHidden;
			numWeights = myFex->cooperativeNet->numWeights;
			policyInputSpecifier = myFex->cooperativeNet->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->cooperativeNet->policyInputSpecifierIndex;
			break;
		case 1:
			numNetInputs = myFex->c2Net->numNetInputs;
			numNetOutputs = myFex->c2Net->numNetOutputs;
			numNetHidden = myFex->c2Net->numNetHidden;
			numWeights = myFex->c2Net->numWeights;
			policyInputSpecifier = myFex->c2Net->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->c2Net->policyInputSpecifierIndex;
			break;
		case 2:
			numNetInputs = myFex->c1Net->numNetInputs;
			numNetOutputs = myFex->c1Net->numNetOutputs;
			numNetHidden = myFex->c1Net->numNetHidden;
			numWeights = myFex->c1Net->numWeights;
			policyInputSpecifier = myFex->c1Net->policyInputSpecifier;
			policyInputSpecifierIndex = myFex->c1Net->policyInputSpecifierIndex;
			break;
		case 3:
			c2Weights = calloc(myFex->c2Net->numWeights,sizeof(double));
			for (i = 0; i < myFex->c2Net->numWeights;i++)
			{
				c2Weights[i] = weights[i];
			}
			C1Weights = calloc(myFex->c1Net->numWeights,sizeof(double));
			for (j = 0; j< myFex->c1Net->numWeights;j++)
			{
				C1Weights[j] = weights[i+j];
			}	
			numNetOutputs = myFex->cooperativeNet->numNetOutputs;
			break;
	}		
	
	/// ------------------------------------------------ ///
	
	//Begin simulation
	for (t=0;t<(myFex->warmup+ myFex->fexT);t++)
	{
		/// ---------- APPLY PHYSICAL FORCING ---------- ///
		double s2Env;
		double s1Env;
		if (t == 0)
		{
			storageS2[t] = myFex->fexInitialStorS2 + myFex->fexFlowsS2[t] - myFex->s2CompensationFixed;
			storageS1[t] = myFex->fexInitialStorS1 + myFex->fexFlowsS1[t] - myFex->s1CompensationFixed - myFex->s1FishReleases[t];
		}
		else
		{
			storageS2[t] = storageS2[t-1] + myFex->fexFlowsS2[t] - myFex->s2CompensationFixed;
			storageS1[t] = storageS1[t-1] + myFex->fexFlowsS1[t] - myFex->s1CompensationFixed - myFex->s1FishReleases[t];
		}
		if (storageS2[t] < 0)
		{
			s2Env = fmax(myFex->s2CompensationFixed + storageS2[t],0);
			storageS2[t] = 0;
		}
		else
		{
			s2Env = myFex->s2CompensationFixed;
		}
		if (storageS1[t] < 0)
		{
			s1Env = fmax(myFex->s1CompensationFixed + myFex->s1FishReleases[t] + storageS1[t],0);
			storageS1[t] = 0;
		}
		else
		{
			s1Env = myFex->s1CompensationFixed + myFex->s1FishReleases[t];
		}
		
		elevation = myFex->storElevationS2[0]*(pow(storageS2[t],myFex->storElevationS2[1]));
		area = myFex->elevationAreaS2[0]*(pow(elevation,myFex->elevationAreaS2[1]));
		
		evap = myFex->fexPets[t];

		double evapS2 = (evap*area)*MM_M2_TO_ML;
		storageS2[t] = fmax((storageS2[t] - evapS2),0);
		
		elevation = myFex->storElevationS1[0]*(pow(storageS1[t],myFex->storElevationS1[1]));
		area = myFex->elevationAreaS1[0]*(pow(elevation,myFex->elevationAreaS1[1]));
		evap = myFex->fexPets[t];
		double evapS1 = (evap*area)*MM_M2_TO_ML;
		storageS1[t] = fmax((storageS1[t] - evapS1),0);
	
		/// -------------------------------------------- ///
		
		
		/// ------- CREATE AND EVALUATE POLICIES ------- ///
		
		double *targetReleases = calloc(numNetOutputs,sizeof(double));
		if (myFex->simulationType != 3)
		{	
			double *netInputs = calloc(numNetInputs,sizeof(double));
			if (policyInputSpecifier[0] == 1)
			{			
				//Normalize storage
				netInputs[policyInputSpecifierIndex[0]] = (storageS2[t] - myFex->deadStorS2)/(myFex->storCapS2 - myFex->deadStorS2 + myFex->inflowNormMaxMinS2[0]);
			}
			if (policyInputSpecifier[1] == 1)
			{
				//Normalize s2 inflows
				netInputs[policyInputSpecifierIndex[1]] = (myFex->fexFlowsS2[t] - myFex->inflowNormMaxMinS2[1])/(myFex->inflowNormMaxMinS2[0] - myFex->inflowNormMaxMinS2[1]);
			}
			
			if (policyInputSpecifier[2] == 1)
			{			
				//Normalize s1 storage
				netInputs[policyInputSpecifierIndex[2]] = (storageS1[t] - myFex->deadStorS1)/(myFex->storCapS1 - myFex->deadStorS1 + myFex->inflowNormMaxMinS1[0]);
			}
			if (policyInputSpecifier[3] == 1)
			{
				//Normalize s1 inflows
				netInputs[policyInputSpecifierIndex[3]] = (myFex->fexFlowsS1[t] - myFex->inflowNormMaxMinS1[1])/(myFex->inflowNormMaxMinS1[0] - myFex->inflowNormMaxMinS1[1]);
			}
			
			if (policyInputSpecifier[4] == 1)
			{
				//Normalize d2down demand
				netInputs[policyInputSpecifierIndex[4]] = (myFex->fexDemandsD2[t] - myFex->demandNormMaxMinD2[1])/(myFex->demandNormMaxMinD2[0] - myFex->demandNormMaxMinD2[1]);
			}
			if (policyInputSpecifier[5] == 1)
			{
				//Normalize pet
				netInputs[policyInputSpecifierIndex[5]] = (myFex->fexPets[t] - myFex->petNormMaxMinS1[1])/(myFex->petNormMaxMinS1[0] - myFex->petNormMaxMinS1[1]);
			}
			if (policyInputSpecifier[6] == 1)
			{
				//Normalize s2 release capacity
				netInputs[policyInputSpecifierIndex[6]] = myFex->linkCapS2D2TimeVarying[t]/myFex->linkCapS2D2;
			}
			if (policyInputSpecifier[7] == 1)
			{
				//Normalize s1 release capacity
				netInputs[policyInputSpecifierIndex[7]] = myFex->linkCapS1D2TimeVarying[t]/myFex->linkCapS1D2;
			}
			if (policyInputSpecifier[8] == 1)
			{
				//Normalize C1 demand
				// netInputs[policyInputSpecifierIndex[8]] = (myFex->fexDemandsC1[t] - myFex->demandNormMaxMinC1[1])/(myFex->demandNormMaxMinC1[0] - myFex->demandNormMaxMinC1[1]);
			}
			if (policyInputSpecifier[9] == 1)
			{
				//Normalize R1 release capacity
				netInputs[policyInputSpecifierIndex[9]] = myFex->linkCapR1S1TimeVarying[t]/myFex->linkCapR1S1;
			}
			if (policyInputSpecifier[10] == 1)
			{
				//Normalize R1 inflows
				netInputs[policyInputSpecifierIndex[10]] = (myFex->fexFlowsR1[t] - myFex->inflowNormMaxMinR1[1])/(myFex->inflowNormMaxMinR1[0] - myFex->inflowNormMaxMinR1[1]);
			}
			if (policyInputSpecifier[11] == 1)
			{
				//Normalize s1 release to R1 capacity
				netInputs[policyInputSpecifierIndex[11]] = myFex->linkCapS1R1TimeVarying[t]/myFex->linkCapS1R1;
			}
			if (policyInputSpecifier[12] == 1)
			{
				//Normalize C1 supplement capacity
				netInputs[policyInputSpecifierIndex[12]] = myFex->linkCapBoreholeTimeVarying[t]/myFex->linkCapBorehole;
			}
			if (policyInputSpecifier[13] == 1)
			{
				//time component sin
				netInputs[policyInputSpecifierIndex[13]] = sin(2*M_PI*myFex->fexDoys[t]/365);
			}
			if (policyInputSpecifier[14] == 1)
			{
				//time component cos
				netInputs[policyInputSpecifierIndex[14]] = cos(2*M_PI*myFex->fexDoys[t]/365);
			}
			if (policyInputSpecifier[15] == 1)
			{
				//normalize R2 flow
				netInputs[policyInputSpecifierIndex[15]] = (myFex->fexFlowsR2[t] - myFex->inflowNormMaxMinR2[1])/(myFex->inflowNormMaxMinR2[0] - myFex->inflowNormMaxMinR2[1]);
			}
			for (i = 0; i < numNetInputs;i++)
			{
				netInputs[i] = fmax(0, fmin(netInputs[i],1)); //In case variables have exceeded their normalization limits
			}
			evaluateNet(targetReleases,netInputs,numNetInputs,numNetOutputs,numNetHidden,weights);
			
			free(netInputs);
		}
		else
		{
			double *netInputsC2 = calloc(myFex->c2Net->numNetInputs,sizeof(double));
			if (myFex->c2Net->policyInputSpecifier[0] == 1)
			{			
				//Normalize storage
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[0]] = (storageS2[t] - myFex->deadStorS2)/(myFex->storCapS2 - myFex->deadStorS2 + myFex->inflowNormMaxMinS2[0]);
			}
			if (myFex->c2Net->policyInputSpecifier[1] == 1)
			{
				//Normalize s2 inflows
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[1]] = (myFex->fexFlowsS2[t] - myFex->inflowNormMaxMinS2[1])/(myFex->inflowNormMaxMinS2[0] - myFex->inflowNormMaxMinS2[1]);
			}
			if (myFex->c2Net->policyInputSpecifier[4] == 1)
			{
				//Normalize d2down demand
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[4]] = (myFex->fexDemandsD2[t] - myFex->demandNormMaxMinD2[1])/(myFex->demandNormMaxMinD2[0] - myFex->demandNormMaxMinD2[1]);
			}
			if (myFex->c2Net->policyInputSpecifier[5] == 1)
			{
				//Normalize pet
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[5]] = (myFex->fexPets[t] - myFex->petNormMaxMinS1[1])/(myFex->petNormMaxMinS1[0] - myFex->petNormMaxMinS1[1]);
			}
			if (myFex->c2Net->policyInputSpecifier[6] == 1)
			{
				//Normalize s2 release capacity
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[6]] = myFex->linkCapS2D2TimeVarying[t]/myFex->linkCapS2D2;
			}
			if (myFex->c2Net->policyInputSpecifier[7] == 1)
			{
				//Normalize s1 release capacity
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[7]] = myFex->linkCapS1D2TimeVarying[t]/myFex->linkCapS1D2;
			}
			
			if (myFex->c2Net->policyInputSpecifier[13] == 1)
			{
				//time component sin
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[13]] = sin(2*M_PI*myFex->fexDoys[t]/365);
			}
			if (myFex->c2Net->policyInputSpecifier[14] == 1)
			{
				//time component cos
				netInputsC2[myFex->c2Net->policyInputSpecifierIndex[14]] = cos(2*M_PI*myFex->fexDoys[t]/365);
			}
			double *targetReleasesC2 = calloc(myFex->c2Net->numNetOutputs, sizeof(double));
			for (i = 0; i < myFex->c2Net->numNetInputs;i++)
			{
				netInputsC2[i] = fmax(0, fmin(netInputsC2[i],1)); //In case variables have exceeded their normalization limits
			}
			evaluateNet(targetReleasesC2,netInputsC2,myFex->c2Net->numNetInputs,myFex->c2Net->numNetOutputs,myFex->c2Net->numNetHidden,c2Weights);
			free(netInputsC2);
			targetReleases[0] = targetReleasesC2[0];
			targetReleases[1] = targetReleasesC2[1];
			free(targetReleasesC2);
			double *netInputsC1 = calloc(myFex->c1Net->numNetInputs,sizeof(double));
			if (myFex->c1Net->policyInputSpecifier[2] == 1)
			{			
				//Normalize s1 storage
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[2]] = (storageS1[t] - myFex->deadStorS1)/(myFex->storCapS1 - myFex->deadStorS1 + myFex->inflowNormMaxMinS1[0]);
			}
			if (myFex->c1Net->policyInputSpecifier[3] == 1)
			{
				//Normalize s1 inflows
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[3]] = (myFex->fexFlowsS1[t] - myFex->inflowNormMaxMinS1[1])/(myFex->inflowNormMaxMinS1[0] - myFex->inflowNormMaxMinS1[1]);
			}

			if (myFex->c1Net->policyInputSpecifier[5] == 1)
			{
				//Normalize pet
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[5]] = (myFex->fexPets[t] - myFex->petNormMaxMinS1[1])/(myFex->petNormMaxMinS1[0] - myFex->petNormMaxMinS1[1]);
			}
		
			if (myFex->c1Net->policyInputSpecifier[7] == 1)
			{
				//Normalize s1 release capacity
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[7]] = myFex->linkCapS1D2TimeVarying[t]/myFex->linkCapS1D2;
			}
			if (myFex->c1Net->policyInputSpecifier[8] == 1)
			{
				//Normalize C1 demand
				// netInputsC1[myFex->c1Net->policyInputSpecifierIndex[8]] = (myFex->fexDemandsC1[t] - myFex->demandNormMaxMinC1[1])/(myFex->demandNormMaxMinC1[0] - myFex->demandNormMaxMinC1[1]);
			}
			if (myFex->c1Net->policyInputSpecifier[9] == 1)
			{
				//Normalize R1 release capacity
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[9]] = myFex->linkCapR1S1TimeVarying[t]/myFex->linkCapR1S1;
			}
			if (myFex->c1Net->policyInputSpecifier[10] == 1)
			{
				//Normalize R1 inflows
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[10]] =  (myFex->fexFlowsR1[t] - myFex->inflowNormMaxMinR1[1])/(myFex->inflowNormMaxMinR1[0] - myFex->inflowNormMaxMinR1[1]);
			}
			if (myFex->c1Net->policyInputSpecifier[11] == 1)
			{
				//Normalize s1 release to R1 capacity
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[11]] = myFex->linkCapS1R1TimeVarying[t]/myFex->linkCapS1R1;
			}
			if (myFex->c1Net->policyInputSpecifier[12] == 1)
			{
				//Normalize C1 supplement capacity
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[12]] = myFex->linkCapBoreholeTimeVarying[t]/myFex->linkCapBorehole;
			}
			if (myFex->c1Net->policyInputSpecifier[13] == 1)
			{
				//time component sin
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[13]] = sin(2*M_PI*myFex->fexDoys[t]/365);
			}
			if (myFex->c1Net->policyInputSpecifier[14] == 1)
			{
				//time component cos
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[14]] = cos(2*M_PI*myFex->fexDoys[t]/365);
			}
			if (myFex->c1Net->policyInputSpecifier[15] == 1)
			{
				//normalize R2 flow
				netInputsC1[myFex->c1Net->policyInputSpecifierIndex[15]] = (myFex->fexFlowsR2[t] - myFex->inflowNormMaxMinR2[1])/(myFex->inflowNormMaxMinR2[0] - myFex->inflowNormMaxMinR2[1]);
			}
			double *targetReleasesC1 = calloc(myFex->c1Net->numNetOutputs, sizeof(double));
			for (i = 0; i < myFex->c1Net->numNetInputs;i++)
			{
				netInputsC1[i] = fmax(0, fmin(netInputsC1[i],1)); //In case variables have exceeded their normalization limits
			}
			evaluateNet(targetReleasesC1,netInputsC1,myFex->c1Net->numNetInputs,myFex->c1Net->numNetOutputs,myFex->c1Net->numNetHidden,C1Weights);
			free(netInputsC1);
			targetReleases[2] = targetReleasesC1[0];
			targetReleases[3] = targetReleasesC1[1];
			free(targetReleasesC1);
		}
		double targS2D2, targS1D2, targR1S1, targS1R1;
		if (myFex->simulationType != 2)
		{ 
			targS2D2 = targetReleases[0];
			targS1D2 = targetReleases[1];
			if (myFex->simulationType != 1)
			{
				targR1S1 = targetReleases[2];
				targS1R1 = targetReleases[3];
			}
			else 
			{
				targR1S1 = 0;
				targS1R1 = 0;
			}
		}
		else
		{
			targS2D2 = 0;
			//NOTE this time gap loops overyear - i.e. Nov->Mar
			//We formulate targS1D2 here from C1 drought plan (publicly available)
			if (myFex->fexDoys[t] >= myFex->startDoyMaxR1D2 || myFex->fexDoys[t] <= myFex->endDoyMaxR1D2)
			{
				targS1D2 = 1;
			}
			else
			{
				int durationNormRelease = myFex->startDoyMaxR1D2-myFex->endDoyMaxR1D2-1;
				int durationMaxRelease = DAYS_IN_YEAR - durationNormRelease;
				double normRelease = (myFex->linkAnnualCapS1D2 - myFex->linkCapS1D2*durationMaxRelease)/durationNormRelease;
				targS1D2 = normRelease/myFex->linkCapS1D2;
			}
			targR1S1 = targetReleases[0];
			targS1R1 = targetReleases[1];
		}
		free(targetReleases);
		
		/// -------------------------------------------- ///
		
		
		/// ---- CREATE FLOWS FROM TARGET DECISIONS ---- ///
				
		volumeS2D2[t] = fmax(targS2D2*myFex->linkCapS2D2,myFex->minS2D2);
		volumeS2D2[t] = fmin(volumeS2D2[t],myFex->linkCapS2D2TimeVarying[t]);
		volumeS2D2[t] = fmin(volumeS2D2[t],storageS2[t]);

		volumeS1D2[t] = fmax(targS1D2*myFex->linkCapS1D2,myFex->minS1D2);
		volumeS1D2[t] = fmin(volumeS1D2[t],myFex->linkCapS1D2TimeVarying[t]);

		//NOTE this time gap loops within year - i.e. Apr->Nov
		if (myFex->fexDoys[t] >= myFex->startDoyNoR1S1 && myFex->fexDoys[t] <= myFex->endDoyNoR1S1)
		{
			volumeR1S1[t] = 0;
		}
		else
		{
			volumeR1S1[t] = targR1S1*myFex->linkCapR1S1;
			volumeR1S1[t] = fmin(volumeR1S1[t],myFex->linkCapR1S1TimeVarying[t]);
			
			volumeR1S1[t] = fmin(volumeR1S1[t],fmax((myFex->fexFlowsR1[t] - myFex->R1Thresh)/2,0)); //where myFex->R1Thresh = (1.16*86.4)Ml/d
			volumeR1S1[t] = fmin(volumeR1S1[t],fmax((myFex->fexFlowsR2[t] + myFex->fexFlowsR1[t] - myFex->R2Thresh),0)); //where myFex->R2Thresh = (3.158*86.4)Ml/d
		}

		double waterAvailableAtR2 = fmax((myFex->fexFlowsR2[t] + myFex->fexFlowsR1[t] - myFex->R2Thresh),0); //How much water you can take without breaking R2 license (any further abstraction will have to be released from S1)
		double waterAvailableAtStoodleigh = myFex->fexFlowsR1[t]*myFex->r2r1ConversionFactor; 
		double waterAvailableForC1 = fmin(waterAvailableAtStoodleigh,waterAvailableAtR2);
		double maximumReleaseWithoutWaste = fmax(myFex->fexDemandD1Fixed - waterAvailableForC1,0);
		
		volumeS1R1[t] = fmax(targS1R1*myFex->linkCapS1R1,0);
		volumeS1R1[t] = fmin(volumeS1R1[t],maximumReleaseWithoutWaste);

		/// -------------------------------------------- ///
		
		
		/// ------------ RESOLVE CONSTRAINTS ----------- ///
		//Check that license allows abstraction
		volumeS1R1[t] = fmin(volumeS1R1[t],fmax(myFex->linkAnnualCapR1S1 - S1R1AnnualTotal,0));
		volumeR1S1[t] = fmin(volumeR1S1[t],fmax(myFex->linkAnnualCapR1S1 - R1S1AnnualTotal,0));
		//No simulataneous release and abstraction from R1
		if (volumeR1S1[t] > 0 && volumeS1R1[t] > 0)
		{
			if (volumeR1S1[t] > volumeS1R1[t])
			{
				volumeR1S1[t] = volumeR1S1[t] - volumeS1R1[t];
				volumeS1R1[t] = 0;
			}
			else
			{
				volumeS1R1[t] = volumeS1R1[t] - volumeR1S1[t];
				volumeR1S1[t] = 0;
			}
		}

		//Mass balance on s1leball from R1
		storageS1[t] = storageS1[t] + volumeR1S1[t];
		volumeS1R1[t] = fmin(volumeS1R1[t],storageS1[t]); // we currently draw s1R1 followed by s1D2 (on the basis that C1 own the dam)
		storageS1[t] = storageS1[t] - volumeS1R1[t]; 
		if (myFex->simulationType !=1)
		{
			volumeS1D2[t] = fmin(volumeS1D2[t],storageS1[t]);
		}

		//Check that license allows abstraction
		volumeS1D2[t] = fmin(volumeS1D2[t],fmax(myFex->linkAnnualCapS1D2 - S1D2AnnualTotal,0)); 
		volumeS2D2[t] = fmin(volumeS2D2[t],fmax(myFex->linkAnnualCapS2D2 - S2D2AnnualTotal,0)); 
		volumeS2D2[t] = fmin(volumeS2D2[t],fmax(myFex->linkMonthlyCapS2D2 - S2D2MonthlyTotal,0)); 

		//Find supply to D2down
		supplyToD2[t] = volumeS1D2[t] + volumeS2D2[t];
		if (supplyToD2[t] > myFex->fexDemandsD2[t] && myFex->simulationType != 2)
		{
			// Ensure no oversupply of D2
			// we currently reduce oversupply by reducing flows in proportion to the flow size
			
			// we also apply this after the minimum flow requirements (e.g. myFex->minS2D2)
			// so the flows may drop below. You can move this up, but then oversupply may happen
			
			
			
			double err = supplyToD2[t] - myFex->fexDemandsD2[t];
			volumeS2D2[t] = volumeS2D2[t] - err*volumeS2D2[t]/supplyToD2[t];
			volumeS1D2[t] = volumeS1D2[t] - err*volumeS1D2[t]/supplyToD2[t];
			supplyToD2[t] = supplyToD2[t] - err;
		}

		//Remove d2down supplies from reservoirs
		storageS2[t] = storageS2[t] - volumeS2D2[t];
		storageS1[t] = storageS1[t] - volumeS1D2[t];

		//Update annual counters
		S1D2AnnualTotal = S1D2AnnualTotal + volumeS1D2[t];
		S2D2AnnualTotal = S2D2AnnualTotal + volumeS2D2[t];
		S2D2MonthlyTotal = S2D2MonthlyTotal + volumeS2D2[t];
		R1S1AnnualTotal = R1S1AnnualTotal + volumeR1S1[t];
		S1R1AnnualTotal = S1R1AnnualTotal + volumeS1R1[t];
		
		//Calculate C1 supply
		supplyToC1[t] = fmin(myFex->fexDemandD1Fixed, waterAvailableForC1 + volumeS1R1[t] - volumeR1S1[t]); 
		volumeBorhole[t] = fmin(myFex->fexDemandD1Fixed - supplyToC1[t],myFex->linkCapBoreholeTimeVarying[t]);
		supplyToC1[t] = supplyToC1[t] + volumeBorhole[t];


		//Reset annual/monthly counters
		if (myFex->fexDoys[t] == DAYS_IN_YEAR)
		{
			R1S1AnnualTotal = 0;
			S1D2AnnualTotal = 0;
			S2D2AnnualTotal = 0;
			S1R1AnnualTotal = 0;
		}
		if ((t+1)%DAYS_IN_MONTH == 0)
		{
			S2D2MonthlyTotal = 0;
		}

		//apply spill
		spillS2[t] = fmax(0, storageS2[t] - myFex->storCapS2);
		storageS2[t] = storageS2[t] - spillS2[t];

		spillS1[t] = fmax(0, storageS1[t] - myFex->storCapS1);
		storageS1[t] = storageS1[t] - spillS1[t];

		/// -------------------------------------------- ///
		
		
		/// --------- CALCULATE TIME-STEP COSTS -------- ///
		if (t > myFex->warmup)
		{
			gsLinC1[t] = volumeBorhole[t];
			gsSqrC1[t] = volumeBorhole[t]*volumeBorhole[t];
			deficitSqrC1[t] = (myFex->fexDemandD1Fixed - supplyToC1[t])*(myFex->fexDemandD1Fixed    - supplyToC1[t]);
			deficitLinC1[t] = (myFex->fexDemandD1Fixed  - supplyToC1[t]);
			deficitSqrD2[t] = (myFex->fexDemandsD2[t] - supplyToD2[t])*(myFex->fexDemandsD2[t] - supplyToD2[t]);
			deficitLinD2[t] = (myFex->fexDemandsD2[t] - supplyToD2[t]);
			if (storageS2[t] < myFex->storageReliabilityThresholdS2)
			{
				storageReliabilityS2[t] = 1;
			}
			if (storageS1[t] < myFex->storageReliabilityThresholdS1)
			{
				storageReliabilityS1[t] = 1;
			}
			
			NPVc2[t] = (volumeS1D2[t]*myFex->linkCostS1D2)/pow((1+myFex->discountRate),t);
			NPVC1[t] = (volumeR1S1[t]*myFex->linkCostR1S1)/pow((1+myFex->discountRate),t);
			NPVtotal[t] = NPVc2[t] + NPVC1[t];
			pumpCostC2[t] = volumeS1D2[t]*myFex->linkCostS1D2;
			pumpCostC1[t] = volumeR1S1[t]*myFex->linkCostR1S1 + volumeBorhole[t]*myFex->linkCostGroundBorehole;
			pumpCostTotal[t] = pumpCostC2[t] + pumpCostC1[t];
		}
		else
		{
			gsLinC1[t] = 0;
			gsSqrC1[t] = 0;
			deficitSqrC1[t] = 0;
			deficitLinC1[t] = 0;
			deficitSqrD2[t] = 0;
			deficitLinD2[t] = 0;
			
				storageReliabilityS2[t] = 0;
			
				storageReliabilityS1[t] = 0;
			
			
			NPVc2[t] = 0;
			NPVC1[t] = 0;
			NPVtotal[t] = 0;
			pumpCostC2[t] = 0;
			pumpCostC1[t] = 0;
			pumpCostTotal[t] = 0;
		}
		/// -------------------------------------------- ///
		
		
		/// -------------- PRINT DYNAMICS -------------- ///
		if (myFex->printDynamics == 1)
		{
			if (t == 0)
			{
				fprintf(fid,"%f %f ", myFex->fexInitialStorS2,myFex->fexInitialStorS1);
			}
			else
			{
				fprintf(fid,"%f %f ", storageS2[t - 1], storageS1[t - 1]);
			}
			fprintf(fid,"%f %f %f %f ", myFex->fexFlowsS2[t], myFex->fexFlowsS1[t], myFex->fexFlowsR1[t],myFex->fexFlowsR2[t]);
			fprintf(fid,"%f %f ",myFex->fexDemandsD2[t],myFex->fexDemandD1Fixed);
			fprintf(fid,"%f %f %f %f ",volumeS2D2[t],volumeS1D2[t],volumeR1S1[t],volumeS1R1[t]);
			fprintf(fid,"%f %f %f %d ",volumeBorhole[t],spillS2[t],spillS1[t],myFex->fexDoys[t]);
			fprintf(fid,"%f %f %f %f\n",evapS2, evapS1,s2Env, s1Env);
		}
		
		/// -------------------------------------------- ///
		if (storageS2[t] < 0)
		{
			if ((storageS2[t] < -(DOUBLE_PRECISION*10000)) && myFex->simulationType != 2)
			{
				printf("warning: s2 storage below 0");
			}
			storageS2[t] = 0;
		}
		
		if (storageS1[t] < 0)
		{
			if ((storageS1[t] < -(DOUBLE_PRECISION*10000)) && myFex->simulationType != 1)
			{
				printf("warning: s1 storage below 0");
			}
			storageS1[t] = 0;
		}
	}
	
	

	free(spillS1);
	free(spillS2);
	free(supplyToD2);
	free(volumeS2D2);
	free(volumeS1D2);
	free(volumeR1S1);
	free(volumeS1R1);
	free(volumeBorhole);
	free(supplyToC1);
	
	if (myFex->simulationType == 3)
	{
		free(c2Weights);
		free(C1Weights);
	}
	
	
	/// ----------- CALCULATE OBJECTIVES ----------- ///
	
	//THIS MAKES ANGELS WEEP... YOU NEED AN OBJECTIVE SPECIFIER INDEX
	
	l = 0;
	if (myFex->objectiveSpecifier[0] == 1)
	{
		// timeObjectives[l] = storageReliabilityS2;
		objs[l] = sumDouble(storageReliabilityS2,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(storageReliabilityS2);
	// }
	
	if (myFex->objectiveSpecifier[1] == 1)
	{
		// timeObjectives[l] = storageReliabilityS1;
		objs[l] = sumDouble(storageReliabilityS1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(storageReliabilityS1);
	// }
	
	if (myFex->objectiveSpecifier[2] == 1)
	{
		// timeObjectives[l] = deficitLinD2;
		objs[l] = sumDouble(deficitLinD2,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitLinD2);
	// }
	
	if (myFex->objectiveSpecifier[3] == 1)
	{
		// timeObjectives[l] = deficitSqrD2;
		objs[l] = sumDouble(deficitSqrD2,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitSqrD2);
	// }
	
	if (myFex->objectiveSpecifier[4] == 1)
	{
		// timeObjectives[l] = deficitLinC1;
		objs[l] = sumDouble(deficitLinC1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitLinC1);
	// }
	
	if (myFex->objectiveSpecifier[5] == 1)
	{
		// timeObjectives[l] = deficitSqrC1;
		objs[l] = sumDouble(deficitSqrC1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(deficitSqrC1);
	// }
	
	if (myFex->objectiveSpecifier[6] == 1)
	{
		// timeObjectives[l] = NPVtotal;
		objs[l] = sumDouble(NPVtotal,myFex->fexT);
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVtotal);
	// }
	
	if (myFex->objectiveSpecifier[7] == 1)
	{
		// timeObjectives[l] = NPVc2;
		objs[l] = sumDouble(NPVc2,myFex->fexT);
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVc2);
	// }
	
	if (myFex->objectiveSpecifier[8] == 1)
	{
		// timeObjectives[l] = NPVC1;
		objs[l] = sumDouble(NPVC1,myFex->fexT);
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(NPVC1);
	// }
	
	if (myFex->objectiveSpecifier[9] == 1)
	{
		// timeObjectives[l] = gsLinC1;
		objs[l] = sumDouble(gsLinC1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(gsLinC1);
	// }
	
	if (myFex->objectiveSpecifier[10] == 1)
	{
		// timeObjectives[l] = gsSqrC1;
		objs[l] = sumDouble(gsSqrC1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(gsSqrC1);
	// }
	
	if (myFex->objectiveSpecifier[11] == 1)
	{
		// timeObjectives[l] = pumpCostTotal;
		objs[l] = sumDouble(pumpCostTotal,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostTotal);
	// }
	
	if (myFex->objectiveSpecifier[12] == 1)
	{
		// timeObjectives[l] = pumpCostC2;
		objs[l] = sumDouble(pumpCostC2,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostC1);
	// }
	
	if (myFex->objectiveSpecifier[13] == 1)
	{
		// timeObjectives[l] = pumpCostC1;
		objs[l] = sumDouble(pumpCostC1,myFex->fexT)/myFex->fexT;
		l=l+1;
	}
	// else if (includetimeObjectives == 1)
	// {
		// free(pumpCostC1);
	// }
	if (myFex->printDynamics == 1)
	{
		fprintf(fid,"#\n//Final Storages, Sc Sw:\n");
		fprintf(fid,"%f %f\n", storageS2[t - 1], storageS1[t - 1]);
		
		fprintf(fid,"#\n//Objectives:\n");
		for (i=0;i<myFex->numObjectives;i++)
		{
			fprintf(fid,"%f ",objs[i]);
		}
		fclose(fid);
	}
	
	free(storageS1);
	free(storageS2);

	// if (includetimeObjectives != 1)
	// {
	free(pumpCostTotal);
	free(pumpCostC2);
	free(pumpCostC1);
	free(NPVc2);
	free(NPVC1);
	free(NPVtotal);
	free(deficitSqrD2);
	free(deficitLinD2);
	free(deficitSqrC1);
	free(deficitLinC1);
	free(gsSqrC1);
	free(gsLinC1);
	free(storageReliabilityS2);
	free(storageReliabilityS1);
	// free(timeObjectives);
	// }
	
	
	/// -------------------------------------------- ///
}

void sortDescending(int N, double *vec)
{
	//Sort the vector (vec) of doubles (length = N) with largest first -- note that sorted values are stored in the original vector
	int i,j;
	double a;
	for (i = 0; i < N; ++i)
    {
        for (j = i + 1; j < N; ++j)
        {
            if (vec[i] < vec[j])
            {
                a =  vec[i];
                vec[i] = vec[j];
                vec[j] = a;
            }
        }
    }
}

void sortAscending(int N, double *vec)
{
	//Sort the vector (vec) of doubles (length = N) with smallest first -- note that sorted values are stored in the original vector
	int i,j;
	double a;
	for (i = 0; i < N; ++i)
    {
        for (j = i + 1; j < N; ++j)
        {
            if (vec[i] > vec[j])
            {
                a =  vec[i];
                vec[i] = vec[j];
                vec[j] = a;
            }
        }
    }
}

void calcPercentile(double *percentiles, double *ensembleObjectives, struct fexEnsemble *myEnsemble,int numObjectives)
{
	//Implementation of MATLAB prctile function
	double pctile = (myEnsemble->percentile/100)*myEnsemble->ensembleSize;
	int index = (int)(pctile+0.5);
	double remainder = pctile - index;
	int i, j;
	for (j=0;j< numObjectives;j++)
	{
		double *ensembleSingleObjective = calloc(myEnsemble->ensembleSize,sizeof(double));
		for (i=0;i<myEnsemble->ensembleSize;i++)
		{
			ensembleSingleObjective[i] = ensembleObjectives[j + i*numObjectives];
		}
		sortAscending(myEnsemble->ensembleSize,ensembleSingleObjective);
		if (index < 1)
		{
			percentiles[j] = ensembleSingleObjective[0];
		}
		else if (index >= myEnsemble->ensembleSize)
		{
			percentiles[j] = ensembleSingleObjective[myEnsemble->ensembleSize - 1];
		}
		else
		{
			percentiles[j] = (0.5 - remainder)*ensembleSingleObjective[index-1] + (0.5 + remainder)*ensembleSingleObjective[index];
		}
		free(ensembleSingleObjective);
	}
}

void calcMean(double *means, double *ensembleObjectives, struct fexEnsemble *myEnsemble,int numObjectives)
{
	double runsum;
	//Implementation of MATLAB prctile function
	int i, j;
	for (j=0;j< numObjectives;j++)
	{
		runsum = 0;
		for (i=0;i<myEnsemble->ensembleSize;i++)
		{
			runsum = runsum + ensembleObjectives[j + i*numObjectives];
		}
		means[j] = runsum/myEnsemble->ensembleSize;
	}
}

void evaluateEnsemble(double* weights, double* objs, double* consts, struct fexEnsemble *myEnsemble)
{
	//A wrapper around simFex that evaluates a fexEnsemble and returns the percentiles in objs.
	int i,j,numObjectives = myEnsemble->ensemble[0]->numObjectives;
	double *ensembleObjectives = calloc(myEnsemble->ensembleSize*numObjectives,sizeof(double));
	for (i=0;i< myEnsemble->ensembleSize;i++)
	{
		double *holderObjectives = calloc(numObjectives,sizeof(double));
		simFex(weights, holderObjectives, 0,myEnsemble->ensemble[i]);
		for (j = 0;j < numObjectives; j++)
		{
			ensembleObjectives[j + i*numObjectives]= holderObjectives[j];
		}
		free(holderObjectives);
	}
	// calcPercentile(objs, ensembleObjectives, myEnsemble, numObjectives);
	calcMean(objs, ensembleObjectives, myEnsemble, numObjectives);
	free(ensembleObjectives);
}


void repopulateFex(struct fexProblem *myFex, struct fexEnsemble *myEnsemble,int seed)
{
	destroyFex(myFex);
	createFex(myFex, myEnsemble, myEnsemble->ensembleT, seed, myEnsemble->simulationType);
	populateSeededVariables(myFex, myEnsemble, seed);
}
///	---------------- EXAMPLE SIMULATION ---------------- ///
// compiles (on cygwin) with command gcc -o fex.exe fex.c qr_solve.c r8lib.c -lm

// int main(){
	// struct fexEnsemble myEnsemble;
	// myEnsemble.ensembleSize = 100;
	// myEnsemble.ensembleT = 1000;
	// myEnsemble.baseSeed = 1;
	// myEnsemble.percentile = 99;
	// myEnsemble.simulationType = 0;
	
	// createEnsemble(&myEnsemble,6);
		
	// int i;
	// double *weights = calloc(myEnsemble.ensemble[0].cooperativeNet->numWeights,sizeof(double));
	// for (i = 0; i < myEnsemble.ensemble[0].cooperativeNet->numWeights;i++)
	// {
		// weights[i] = 0;
	// }
	// double *objective = calloc(myEnsemble.ensemble[0].numObjectives,sizeof(double));

	// evaluateEnsemble(weights, objective, 0, &myEnsemble);

	// for (i = 0; i < myEnsemble.ensemble[0].numObjectives;i++)
	// {
		// printf("%f ",objective[i]);
	// }
	// printf("\n");
	// destroyEnsemble(&myEnsemble);
	// free(objective);
	// free(weights);
// }

///ALT NET
// void evaluateNet(double *outputs, double *inputs, int numNetInputs, int numNetOutputs, int numNetHidden, double *weights)
// {
	// int i,j;
	
	// double *hiddenInputs = calloc(numNetHidden, sizeof(double));
	// for (i=0;i<numNetHidden;i++)
	// {
		// hiddenInputs[i] = weights[i];
		
		// for (j=0;j<numNetInputs;j++)
		// {
			// hiddenInputs[i] = hiddenInputs[i] + inputs[j]*weights[numNetHidden+numNetInputs*i+j];
		// }
		// hiddenInputs[i] = exp(-(hiddenInputs[i]*hiddenInputs[i]));
	// }
	// double *hiddenOutputs = calloc(numNetOutputs, sizeof(double));
	// for (i=0;i<numNetOutputs;i++)
	// {
		// hiddenOutputs[i] = weights[numNetHidden*(1+numNetInputs) + i];
		// for (j=0;j<numNetHidden;j++)
		// {
			// hiddenOutputs[i] = hiddenOutputs[i] + hiddenInputs[j]*weights[numNetHidden*(1+numNetInputs)+numNetOutputs+numNetHidden*i+j];
		// }
		// if (hiddenOutputs[i] < -1)
		// {
			// outputs[i] = 0;
		// }
		// else if (hiddenOutputs[i] < 1)
		// {
			// outputs[i] = (hiddenOutputs[i]+1)/2;
		// }
		// else
		// {
			// outputs[i] = 1;
		// }
	// }
	// free(hiddenInputs);
	// free(hiddenOutputs);
// }
