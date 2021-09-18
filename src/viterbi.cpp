#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <cassert>
#include <cmath>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

using boost::math::poisson;
using boost::math::pdf;
using boost::math::negative_binomial_distribution;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;

int MAX_CN=10;
float MISMAP_RATE=0.01;
using namespace std;
double lepsi=-10000;


static void printModel(vector<vector<double> > &transP,
                       int nStates)
{

  cerr << "\nTRANS: \n";
  for (int r=0;r<nStates;r++)
    {
      cerr << r <<": ";
      for (int c=0;c<nStates;c++)
        {
	  cerr << transP[r][c] << " ";
        }
      cerr << "\n";
    }
  cerr<< "\n";
}

static void printEmissions(vector<vector<double> > &emisP)

{
  cerr << "Emissions matrix: \n";
  for (size_t r=0;r<emisP.size(); r++) 
    {
      for (size_t c=0;c<emisP[r].size(); c++ )
	{
	  cerr << std::setw(4) << emisP[r][c];
	  if (c+1 < emisP[r].size()) 
	    {
	      cerr << " ";
	    }	  
	}
      cerr << endl;
    }
}


double max_over_row(vector<vector<double> > &v , size_t col ,size_t nStates ){

  double maxi=-1 * (std::numeric_limits<double>::max()) ;
  for(size_t i=0;i< nStates;i++){
    maxi=std::max(v[i][col], maxi);
  }
  return maxi;
}


double max_over_rows(vector<vector<double> > &v , size_t col ,vector<vector<double> > &v2 , size_t nextState,size_t nStates ){
  double maxi2=-1 * (std::numeric_limits<double>::max()) ;

  for(size_t i=0;i< nStates;i++){
    maxi2=std::max( v[i][col] + v2[i][nextState] , maxi2);

  }
  return maxi2;
}

double LgNegBinom(int cn, int cov, float Hmean, float Hvar) {
  double result=0;
  float r, p;

  if (Hmean == 0 or Hvar == 0) {
    return 0;
  }
  //
  // Since actual variance is unknown at high copy number states, assume linear increase.
  //
      
  if (Hmean==0) {//no alignment in contig
    if (cn!= 0)
      result=lepsi;
    else
      result=0;    
  }
  else if (cov >= MAX_CN*Hmean) {//max_obs filtered previously
    if(cn!= MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    //
    // Use poisson for 0-state assuming it's a random mismap process.
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else 
      {
	result=log(prob);
      }
  }
  else{
    Hmean*=cn;
    Hvar*=cn;
    p=Hmean/Hvar;
    r=Hmean*(Hmean/Hvar)/(1-Hmean/Hvar);
    
    negative_binomial_distribution<double> distribution(r,p);    
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else
      {
	result=log(prob);

      }
  }
  return result;  
}

double LgPrpoiss(int cn,  int cov, int Hmean) {
  double result=0;

  // Boundary case
  if (Hmean==0){//no alignment in contig
    if(cn!= 0)
      // Cannot emit 0 reads from non-zero states.
      result=lepsi;
    else
      // Can only emit 0
      result=0;
  }
  else if (cov >= MAX_CN*Hmean){//max_obs filtered previously
    //
    // Have state that catches ultra-high probability states.
    if(cn != MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) {
      result = lepsi;      
    }
    else {
      result = log(prob);
    }
  }
  else {
    poisson distribution(cn*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0)
      {
	result=lepsi;
      }
    else 
      {
	result=log(prob);	
      }
  }
  return result;
}




static void correctModel(vector<vector<double> > &transP,
                         int nStates)
{
  double sum;
  for (int i=0;i<nStates;i++)
    {
      sum = 0;
      for (int j=0;j<nStates;j++)
	sum+= std::exp(transP[i][j]);
      for (int j=0;j<nStates;j++)
	transP[i][j]= log(std::exp(transP[i][j])/sum);
    }
}//correctModel


void viterbi( vector<double> &startP,
	      vector<vector<double> > &transP,
	      vector<vector<double> > &emisP,
	      vector<size_t> &observations,
	      size_t nStates ,
	      size_t mean,
	      vector<size_t> &  viterbiPath ,size_t nObservations, size_t max_obs){

  //size_t  nObservations  = observations.size();
  vector<vector<double> > v(nStates, vector<double>(nObservations)  );
  vector<vector<double> > opt(nStates, vector<double>(nObservations)  );  
  // Init
  size_t obs = std::min(max_obs , observations[0]);
  for(size_t i=0;i<nStates;i++)
    {
      v[i][0] =  startP[i] + emisP[i][obs] ;
    }
  // Iteration
    
  for(size_t k=1 ; k<nObservations ; k++)
    {
      size_t obs = std::min(max_obs, observations[k]);
      for(size_t i=0;i<nStates;i++)
        {
	  double maxProb = v[0][k-1] + transP[0][i];
	  int maxState=0;
	  for(size_t j=0;j<nStates;j++)
            {
	      double rowProb = v[j][k-1] + transP[j][i];
	      if (rowProb > maxProb) {
		maxState=j;
		maxProb=rowProb;
	      }
            }
	  v[i][k] = emisP[i][obs] + maxProb;
	  opt[i][k] = maxState;
        }
    }
  /*
  for (size_t k=1; k <nObservations; k++) {
    cout << k << "\t" << observations[k] << "\t";
    for (size_t i=0; i < nStates; i++) {
      cout << std::setw(8) << v[i][k] << " ";
    }
    cout << endl;
  }

  for (size_t k=1; k <nObservations; k++) {
    cout << k << "\t" << observations[k] << "\t";
    for (size_t i=0; i < n States; i++) {
      cout << std::setw(4 ) << opt[i][k] << " ";
    }
    cout << endl;
  }
  */
      
  // Traceback
  for(size_t i=0;i<nStates;i++)
    {
      if( max_over_row(v,nObservations-1,nStates) == v[i][nObservations-1] )
        {
	  viterbiPath[nObservations-1] = i;
	  break;
        }
    }
  size_t lastObservation=nObservations-2;
  //
  // Find highest-scoring final entry.
  //
  size_t rowIndex=nObservations-2;
  double maxScore=v[0][rowIndex];
  int    maxState=0;
  for (size_t i=1;i < nStates; i++) {
    if (maxScore < v[i][rowIndex]) {
      maxState=i;
      maxScore=v[i][rowIndex];
    }
  }
  viterbiPath[lastObservation] = maxState;
  
  
  for( size_t f=lastObservation; f > 0; f--) {
      viterbiPath[f] = maxState;
      maxState=opt[maxState][f];
    }

}//viterbi




int main(int argc, const char * argv[]) {

  if (argc < 7) {
    std::cerr << "usage: " << argv[0] << " <coverage observation> <mean coverage> <output prefix> <scalar> <epsi> <Diploid(0)/Haploid(1)> " << endl
      //	      << " --model (string) poisson|negbinom . Use either Poisson or negative binomial for emission probability." << endl
	      << " --nbvar   (float)  var                Use a negative binomial emission model with the given data variance." << endl;    
    return EXIT_FAILURE;
  }

  const string prefix_file(argv[3]);
  //const string mean_file(argv[2]);


  //observations--------------------------------------------------------


  vector<size_t> observations;
  const string filename2(argv[1]);
  std::ifstream file;
  size_t coverage;
  file.open(filename2.c_str());
  if (file.good() == false) {
    cerr << "ERROR. Could not open " << filename2 << endl;
    exit(1);
  }
  int i=0;
  size_t maxx=0;
  size_t sum=0;
  int div=1;
  float var=-1;
  
  if (strcmp(argv[6], "0") == 0) {
    div=2;
  }
  int argi=7;
  typedef enum { POIS, NEG_BINOM  } MODEL_TYPE;
  MODEL_TYPE model=POIS;
    
  while (argi < argc ) 
    {
      if (strcmp(argv[argi], "--nbvar")==0) 
	{
	  argi++;
	  if (argi < argc) 
	    {
	      model=NEG_BINOM;
	      var=atof(argv[argi]);	      
	    }
	}
      ++argi;      
    }

  if (model == NEG_BINOM and var < 0) 
    {
      std::cerr << "ERROR. Using negative binomial model without specifying variance." << endl;
      exit(0);      
    }
        
  //int maxi=0;
  while(file >> coverage)
    {
      observations.push_back(coverage);
      sum=sum+coverage;
      if(coverage>maxx)
        {
	  maxx=coverage;
        }
      i++;
    }
    file.close();

    size_t nObservations=observations.size();


  size_t mean= std::floor( std::stoi(argv[2]) /div );
  const size_t calc_mean = std::floor( (sum/nObservations)  /div );

  //max cov value observed or upper cov bound -> max nState---------------
  size_t max_obs = std::min( ( (size_t) std::ceil(MAX_CN*mean)), maxx);
  int maxObsCN=std::ceil(max_obs/mean);
  size_t nStates= std::min( maxObsCN , MAX_CN  ) + 1; //+1 zeroth state
  MAX_CN=nStates+1;
  std::cerr << "nstates " <<nStates << endl;
  std::cerr << "max_cov "<< max_obs << endl;
  std::cerr << "mean " << mean << " calc mean " << calc_mean << endl;
  //----------------------------------------------------------------------


  if (mean == 0) {
    std::cerr << "mean is zero, Exiting" << endl;
    return EXIT_FAILURE;
  }

  vector<double> startP(nStates);

  for(size_t i=0;i<(nStates);i++){
    startP[i]=log(1./(nStates));
  }

  // trans prob

  poisson distribution1(3*mean);
  double result3=pdf(distribution1, 3*mean);

  poisson distribution2(2*mean);
  double result2=pdf(distribution2, 3*mean);

  double epsi23 = result2/result3;
  const double lepsi = stod(argv[5]);
  //*300
  const double scale = std::stod(argv[4]);

  //no epsi
  double beta =  lepsi + ( scale * log(epsi23))  ;
  //mean no. of bins for cn=3 call

  vector<vector<double> > transP(nStates, vector<double>(nStates));
  for (size_t i=0;i<nStates;i++)
    {
      for (size_t j=0;j<nStates;j++)
        {
	  if(i==j)
            {
	      transP[i][j]= log(1 - std::exp(beta) );
            }
            
            
	  else
            {
	      if ( j==3)
                {
		  transP[i][j]= beta - log(nStates-1) ;

                }
	      else
                {
		  transP[i][j]= beta - log(nStates-1);
                
                }
            }
        }
    }

  //  correctModel(transP,nStates);
  // printModel(transP,nStates);


  vector<vector<double> > emisP(nStates, vector<double>(max_obs+1));

  for (size_t i=0;i<nStates;i++){
    for (size_t j=0;j<=max_obs;j++){
      if (model == POIS) {	
	emisP[i][j]=LgPrpoiss( (int) i , j , (int) mean );
      }
      else 
	{
	  emisP[i][j]=LgNegBinom((int)i, (int) j, mean, var);
	}      
    }
  }
  
    
  printModel(transP,nStates);
  printEmissions(emisP); 

  vector<size_t> viterbiPath(nObservations);

  viterbi(startP, transP,emisP, observations, nStates,mean,viterbiPath,nObservations,max_obs);
  string v("viterout.txt");
  string filen = prefix_file + "." +v;

  std::ofstream output;
  output.open(filen.c_str());

  for (size_t i=0;i< nObservations; i++){
    output<<viterbiPath[i]<<endl;
  }
  output.close();

  return 0;
}//main
