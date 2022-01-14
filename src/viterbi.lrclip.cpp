

#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <istream>
#include <cassert>
#include <cmath>
#include <random>
#include <boost/math/distributions/poisson.hpp>
using boost::math::poisson;
using boost::math::pdf;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;


static double max_over_row(vector<vector<double>> &v , size_t col ,size_t nStates )
{

    double maxi=-1 * (std::numeric_limits<double>::max()) ;
    
    for(size_t i=0;i<= nStates;i++)
    {
            maxi=std::max(v[i][col], maxi);
    }
    return maxi;
}


static double max_over_rows(vector<vector<double>> &v , size_t col ,vector<vector<double>> &v2 , size_t nextState,size_t nStates )
{
    double maxi2=-1 * (std::numeric_limits<double>::max()) ;

    for(size_t i=0;i<= nStates;i++)
    {
            maxi2=std::max( v[i][col] + v2[i][nextState] , maxi2);

    }
    return maxi2;
}

/*
// try/add negB()
*/
//0left,1normal,2right
bool is_clip(int index, int nStates, int state)
{
    return (std::floor(   (double) (index/(nStates/3) )  ) ==state);

}

int return_cn(int i,int nStates)
{
    int result =std::floor(   (double) ( i/(nStates/3) )  );
    return (i-((nStates/3)*result));
}

double emissionPr( size_t index,
    size_t state,
    vector<size_t> & observations,
    vector<int> &cl_observations,
    vector<vector<double>> &emisP,
    size_t mean ,
    size_t max_obs,
    size_t nStates,
    int clip_thr )
{
    const double epsi=1e-99;
    size_t cov = observations[index];
    //cout<<cov;
    if (cov > max_obs)
    {
        //cov=max_obs;
        if(state == nStates){
            return log(1-epsi);
        }
        else{
            return log(epsi);
        }
    }//truncate possible coverage observations

    double clip_pr=0;
    double result;
    int clip = cl_observations[index+1]; //shift clip by -1 w.r.t. coverage
    // clip threshold 3

    if (is_clip((int) state, (int) nStates),2 ) //right clipping
    {
        if (clip <= (-1*(clip_thr)) )
            clip_pr=1-epsi;
        else
            clip_pr=epsi;
    }//right clip
    else if(is_clip((int) state, (int) nStates) ,0)
    {
        if (clip >= clip_thr)
            clip_pr=1-epsi;
        else
            clip_pr=epsi;
    }//left clip
    else if(is_clip((int) state, (int) nStates ), 1  ) 
    {
        if (clip > (-1*(clip_thr)) && clip < clip_thr)
            clip_pr=1-epsi;
        else
            clip_pr=epsi;
    }

    result = emisP[state][cov] + log(clip_pr) ;

    return result;
}




void viterbi( vector<double> &startP,
    vector<vector<double>> &transP,
    vector<vector<double>> &emisP,
    vector<size_t> &observations,
    vector<int> &cl_observations,
    size_t nStates ,
    size_t mean,
    vector<size_t> &viterbiPath ,
    size_t nObservations,
    size_t max_obs,
    vector<vector<double>> &v,
    int clip_thr)
{

    // Init
    for(size_t i=0;i<=nStates;i++)
    {
 //       double max_ini= -1 * (std::numeric_limits<double>::max());

        v[i][0] =  startP[i] + emissionPr( 0, i , observations,cl_observations,emisP, mean, max_obs,nStates ,clip_thr);
   
    }
    // Iteration
    for(size_t k=1 ; k<nObservations ; k++)
    {
        for(size_t i=0;i<=nStates;i++)
        {
            double maxi = -1 * (std::numeric_limits<double>::max());
            for(size_t j=0;j<=nStates;j++)
            {
                double temp = v[j][k-1] + transP[j][i];
                //maxi = std::max(maxi, temp);
                if(maxi<temp)
                {
                    maxi=temp;
          //          j_index=j;

                }

            }
            v[i][k] =  emissionPr( k,i,observations,cl_observations,emisP, mean,max_obs,nStates,clip_thr)  + maxi;
        }
 //       viterbiPath[k]=j_index;

    


    }
 //   viterbiPath[]



// Traceback
//can be improved to O(n)
    for(size_t i=0;i<=nStates;i++)
    {
        if( max_over_row(v,nObservations-1,nStates) == v[i][nObservations-1] )
        {
            viterbiPath[nObservations-1] = i;
            break;
        }
    }
    size_t k=nObservations-2;
    for( size_t f=0;f<nObservations-1;   f++ )
    {
        for(size_t i=0;i<=nStates;i++)
        {
            //compute max value in column
            double max_rows = max_over_rows(v,k-f,transP, viterbiPath[(k-f)+1], nStates );
            if( max_rows    == v[i][k-f]+ transP[i][viterbiPath[(k-f)+1] ]   )
            {
                viterbiPath[k-f] = i;
                break;
            }
        }
    }

}//viterbi

static void correctModel(vector<vector<double>> &transP,
                         int nStates)
{
    double sum;
    for (int i=0;i<=nStates;i++)
    {
        sum = 0;
        for (int j=0;j<=nStates;j++)
            sum+=exp(transP[i][j]);
        for (int j=0;j<=nStates;j++)
        {


            transP[i][j]=log(exp(transP[i][j])/sum);

        }
    }

/*

    cout<<"\n\n"<<"transP:after correction"<<endl;

    for(int i=0;i<=nStates;i++){
    cout<<i<<endl;
    for(int j=0;j<=nStates;j++){
    cout<<j<<":"<<transP[i][j]<<"  ";


    }cout<<endl;

    }

*/
// verbose



}//correct model



int main(int argc, const char * argv[]) 
{

    if (argc < 8) 
    {
        std::cerr << "usage: " << argv[0] << " <coverage observation>  <clip observation> <output prefix> <scaler> <epsi> <clip threshold> <provide mean[0/1]> <Hmean> " << endl;
        return EXIT_FAILURE;
    }

    const string prefix_file(argv[3]);
    //const string mean_file(argv[2]);

    //ouble input_epsi;
    //observations--------------------------------------------------------
    vector<size_t> observations;
    vector<int> cl_observations;

    const string filename2(argv[1]);
    std::ifstream file;
    size_t inputString;
    file.open(filename2);
//    size_t i=0;
    size_t max=0;
    size_t sum=0;
    //int maxi=0;

    while(file >> inputString)
    {
        observations.push_back(inputString);
        sum=sum+inputString;
        //observations[i]=inputString;
        if(inputString>=max)
        {
            max=inputString;
            //      maxi=i;
        }
       // i++;
    }
    file.close();


    const string filename(argv[2]);
    std::ifstream file2;

    int inputString2;
    file2.open(filename);
    while(file2 >> inputString2)
    {
        cl_observations.push_back(inputString2);
    }
    //pad last entry for offset of 1
    cl_observations.push_back(0);
    file2.close();


    size_t nObservations=observations.size();
    //const size_t mean = std::stoi(argv[4]);
    size_t mean = std::floor( (sum/nObservations)  /2 );
    if (argv[7]==1){
        mean=std::stoi(argv[8]);
    }

    //max cov value observed or upper cov bound -> max nState---------------
    size_t max_obs = std::min( ( (size_t) std::ceil(30.497*mean) ) , max);

    size_t nStates = (( ( (size_t) (std::ceil((double)max_obs/mean)) ) + 1 ) * 3);//+1 zeroth state

    int clip_thr= std::stoi(argv[6]);
    //size_t nStates_noclip = (((size_t) std::ceil(max_obs/mean))) +1;//+1 zeroth state

    cout<<std::ceil((double)max_obs/mean)<<" "<<max_obs/mean<<" "<<(double)sum/nObservations<<endl;

    cout<<"nstates "<<(nStates)<<endl;
    cout<<"cov_bound "<<max_obs<<" max_observed "<<max<<endl;
    cout<<"Haploid mean "<<mean<<endl;

    //----------------------------------------------------------------------

    const double epsi=1e-99;
    vector<double> startP(nStates+1);

    for(size_t i=0;i<=(nStates);i++)
    {
            startP[i]=log(1./(nStates+1));
            //nstates+1 => catch-all state
    }

    poisson distribution1(3*mean);
    double result3=pdf(distribution1, 3*mean);

    poisson distribution2(2*mean);
    double result2=pdf(distribution2, 3*mean);

    double epsi23 = result2/result3;

    //log(epsi);                                    
    const double input_epsi = std::stod(argv[5]);

    //*300
    const double scale = std::stod(argv[4]);

    double beta =  log(input_epsi) + ( scale *log(epsi23)) ;//log(nStates-1)
    int InStates=(int)nStates;

    // trans prob                   catch-all state
    vector<vector<double>> transP(nStates+1, vector<double>(nStates+1));
    for (int i=0;i<=InStates;i++)
    {
        for (int j=0;j<=InStates;j++)
        {
            //1 left clip
            if( is_lclip(i,InStates))// i< (InStates/3) )
            {
                if(is_lclip(j,InStates) ) //1-1
                {
                   // if(i-j>=0)
                        transP[i][j]=beta;
                }
                else if(is_nclip(j,InStates))///1-2
                {
                    int jj=return_cn(j,InStates);
                    if( (i - jj) > 0)
                        transP[i][j]=beta ;
                    else if( (i - jj) == 0)//diagonal
                        transP[i][j]=beta ;//+log(10);
                    else
                        transP[i][j]=log(1 - (( (2*(InStates/3))+ 1 + i ) * exp(beta) ) ) - log( (InStates/3)-i); //beta +log(100);
                }
                else if(is_rclip(j,InStates)) //1-3
                {
                   // if((i - ( j - ((InStates/3)*2)  )) == 0)
                        transP[i][j]=beta;
                }
                else if(j==InStates){//catch-all state
                    transP[i][j]=beta;//-log(1000);//log( 1-  (   ((2*i)+1+((InStates/3)-1 )  *exp(beta)  )   +  ((2+ ((InStates/3)-i) )* exp(log(10)+beta))  + (((InStates/3)-i)*exp(log(100)+beta) )   ));
                }
            }//1
            
            //2 - neutral
            else if(is_nclip(i,InStates))
            {
                if(is_lclip(j,InStates) ) //2-1
                    transP[i][j]=beta;
                else if(is_nclip(j,InStates) )
                { //2-2

                    if( (i -  j ) == 0)
                        transP[i][j]=log(1-  ((InStates) *exp(beta))  );
                    else
                        transP[i][j]=beta;
                }
                else if (is_rclip(j,InStates))//2-3
                    transP[i][j]=beta;
                else if(j==InStates)
                    transP[i][j]=beta;//-log(10000);
            }//2

            //3 right clip
            else if(is_rclip(i,InStates))//i rclip
            {
                if(is_lclip(j,InStates) )//3-1
                {
                    //if((j - ( i - ((InStates/3)*2)  )) == 0)
                    transP[i][j]=beta ;//+log(10);
                }
                else if(is_nclip(j,InStates))//3-2
                {
                
                    int ii=return_cn(i,InStates);
                    int jj=return_cn(j,InStates);

                    if(     (ii -  jj) < 0)  //upper diagonal
                        transP[i][j]=beta ;
                    else if( ii-jj == 0)//diagonal
                        transP[i][j]=beta ;//+log(10);
                    else //lower diagonal
                        transP[i][j]=log(1 - (((3*(InStates/3))+1 - ii) * exp(beta) ))-log(ii);//- ((InStates/3)*2));//beta +log(100);
                }
                else if (is_rclip(j,InStates))//3-3
                {

                    transP[i][j]=beta ;//+log(10);

                }
            //    else if (   (i -  j ) <= 0 && i!=InStates)//diagonal
                else if (j==InStates )
                {
                    transP[i][j]=beta;//-log(100000);//log( 1- (   (i*exp(log(100)+beta))  + ((2+i)* exp(log(10) +beta )  )  +  (  ((InStates-3)-(2*i)  )*exp(beta) )        )   );
                }

            }
            else if(i==InStates)
            {
                transP[i][j]=log(1./(nStates +1));

            }//cout<<endl;
        }
    }//transP

    poisson distribution3(0.01*mean);
    //poisson emmision
    
    vector<vector<double>> emisP(nStates, vector<double>(max_obs+1));
    for (size_t i=0;i<nStates;i++)
    {
        for (size_t j=0;j<=max_obs;j++)
        {
            //i%nStates ==2
            int ii=return_cn((int)i,InStates);
            if(ii==0)
            {
                double zero_result=pdf(distribution3, j);
                emisP[i][j]=log(zero_result);
            }
            else
            {
                poisson distribution4(ii*mean);
                double result=pdf(distribution4, j);
                emisP[i][j]=log(result);
            }
        }
    }

/*
cout<<"emisP"<<endl;
for(size_t i=0;i<nStates;i++){
    cout<<i<<endl;
    for(size_t j=0;j<=max_obs;j++){
        cout<<j<<":"<<emisP[i][j]<<"  ";
    }cout<<endl;
}


cout<<"transP:before"<<endl;

for(size_t i=0;i<=nStates;i++){
    cout<<i<<endl;
    for(size_t j=0;j<=nStates;j++){
        cout<<j<<":"<<transP[i][j]<<"  ";
    }cout<<endl;
}
*/
//if verbose

    correctModel(transP,InStates);

    vector<size_t> viterbiPath(nObservations);
    
    //posterior pr
    vector<vector<double>> v(nStates+1, vector<double>(nObservations) );

    viterbi(startP, transP, emisP, observations,cl_observations, nStates,mean,viterbiPath,nObservations,max_obs,v,clip_thr);

    string vv("viterout.txt");
    string filen = prefix_file + "." +vv;

    string p("posterior.txt");
    string pfilen = prefix_file + "." +p;


    //output posterior
    std::ofstream output;
    output.open(pfilen.c_str());

    for (size_t i=0;i< nObservations; i++)
    {
        output<<i+1<<"\t";
        for (size_t l=0; l<=nStates; l++)
        {
            output<<l<<":"<<v[l][i]<<"\t";

        }output<<endl;
    }
    output.close();

    std::ofstream output2;
    output2.open(filen.c_str());

    //decoder
    for (size_t i=0;i< nObservations; i++)
    {
        int path_i =return_cn( (int) viterbiPath[i], InStates);
        if(is_lclip(viterbiPath[i],InStates))
        {
            output2<<path_i<<"\tL\t"<<viterbiPath[i]<<endl;
        }
        else if(is_nclip(viterbiPath[i],InStates))
        {
            output2<<path_i<<"\tN\t"<<viterbiPath[i]<<endl;
        }
        else if(is_rclip(viterbiPath[i],InStates))
        {
            output2<<path_i<<"\tR\t"<<viterbiPath[i]<<endl;
        }
        else
            output2<<((nStates/3))<<"\tH\t"<<viterbiPath[i]<<endl;
    }
    output2.close();

return 0;
}//main
