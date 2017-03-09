#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector> /*I need dis so bad */
#include <iostream>
#include <algorithm>

#define MaxCycleLength 100
#define N 8
#define double_N 8.0

using namespace std;

int sigmaToInt(int16_t *Sigma);
/* source: http://c-faq.com/lib/gaussian.html */
double gaussrand(double variance){
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if(phase == 0) {
		do {
			double U1 = ((double)rand() / RAND_MAX);
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return sqrt(variance)*X;
}

int stepFunction(double x){
    if(x > 0) return 1;
    if(x < 0) return 0;
}

int sign(double x){
    if(x > 0) return 1;
    if(x < 0) return -1;
}
/*check if 2 arrays match element by element */
int checkIfEqual(int16_t *array1, int16_t *array2){
int ii;
    for(ii = 0; ii < N; ii++) {
        if (array1[ii] != array2[ii]) return 0;
    }
    return 1;
}
/* do the time evolution and look for a cycle (should I check a cycle by letting it run a little more?) */
int advanceTime(int16_t *Sigma, double J[N][N], vector<vector<uint64_t> > &attractors,vector<uint64_t> &basin){
    int i,j,k,k1,k2,k3;
    int16_t foundConfigs[N*100][N]; //It usually doesn't take that many steps to find a cycle I found 217 at most for N=16 so N*100 should be more than enough
    int configs = 1; //amount of time steps taken AKA configurations visited
    double h[N]; //the input equation (2)
    vector<uint64_t> newFoundAttractor;
    int attractor;
    for(i=0; i<N; i++) //save the initial configuration
        foundConfigs[0][i] = Sigma[i];

    while(configs < (int)pow(2,N) ){ //for safety while(true) seems too dangerous
        for(i=0; i<N; i++){
            h[i] = 0;
            for(j=0; j<N; j++){
                h[i] += J[i][j]*Sigma[j]; //equation (2)
            }
        }
        for(i=0; i<N; i++)
            foundConfigs[configs][i] = Sigma[i] = sign(h[i]); //equation (1)
        //printf("\t %d \n",sigmaToInt(Sigma)); //too useful to delete
        /*check if the newly found configuration Sigma matches a previous one, if so we have a cycle*/
        for(k=configs-1; k>=0; k--){
            if(checkIfEqual(foundConfigs[k],Sigma)){
                for(k1=k; k1<configs; k1++){ /*save all the sigmas in the attractor as numbers*/
                    newFoundAttractor.push_back(sigmaToInt(foundConfigs[k1]));
                }
                /* sort the array so that permutations of cycles are now regarded as equal (do we have to save the unsorted one?)*/
                sort(newFoundAttractor.begin(),newFoundAttractor.begin() + newFoundAttractor.size());
                for(k3=0; k3<attractors.size(); k3++){ /* check if we have found this attractor before */
                    if((attractors[k3] == newFoundAttractor)){
                        basin[k3]++; /* + one more initial configuration which flows to this attractor*/
                        return attractors[k3].size(); /* the length of the cycle */
                    }
                }
                /* if not, add it to the attractors*/
                basin.push_back(1); /* which means we have a new basin of attraction with 1 initial configuration flowing to it so far*/
                attractors.push_back(newFoundAttractor); /* add it to the array of attractors*/
                attractor = attractors.size() - 1;
                return attractors[attractor].size(); /* the length of the cycle */
            }
        }
        configs++;
    }
    return -1;
}

void setSigmaRandom(int16_t *Sigma,int x){
        int l;
        srand(222 + x);
        for(l=0; l<N; l++)
            Sigma[l] = sign( rand()-RAND_MAX*0.5 );
}
/* going from n=0 to n=2^N goes through every possible initial configuration*/
void setSigma(int16_t *Sigma,int n){
    int c,k;
    for (c = N-1; c >= 0; c--)
    {
    k = n >> c;
    if (k & 1)
      Sigma[c] = 1;
    else
      Sigma[c] = -1;
    }

}
/* goes through every single initial configuration and puts the cycle lengths in the res array*/
void getAllCycleLengths(double J[N][N],uint64_t *res,double &l_res, double &l_prime_res,double &fixed_points){
    int l;
    int cycle_length;
    double n,result;
    int16_t Sigma[N]; /*int8_t are single byte signed ints since it's only 1 or -1 anyway*/
    vector<vector<uint64_t> > attractors;
    vector<uint64_t> basin;
    for(l=0; l<pow(2,N); l++){
        setSigma(Sigma,l);
        cycle_length = advanceTime(Sigma,J,attractors,basin);
//        if(cycle_length<MaxCycleLength)
//            res[cycle_length]++;

    }
    n = attractors.size() ;
    for(l=0; l<n; l++){
        if(attractors[l].size()==1){
            fixed_points++;
        }
    }
    result = 0;
    for(l=0; l<n; l++)
        result += attractors[l].size(); /* cycle length of the attractor */
    l_res = (1.0/n)*result;

    result = 0;
    for(l=0; l<n; l++)
        result += attractors[l].size()*basin[l];
    l_prime_res = (1.0/pow(2,N))*result;
}
/* goes through n random initial configurations and puts the cycle lengths in the res array. x is a seed for the rand() function*/
void getRandomCycleLengths(double J[N][N],int configs,uint64_t *res,double &l_res, double &l_prime_res){
    int l;
    int l1;
    int cycle_length;
    int fixed_points = 0;
    double n,result;
    int16_t Sigma[N]; /*int8_t are single byte signed ints since it's only 1 or -1 anyway*/
    vector<vector<uint64_t> > attractors;
    vector<uint64_t> basin;

    for(l1=0; l1<configs; l1++){
        setSigmaRandom(Sigma,l1);
        cycle_length = advanceTime(Sigma,J,attractors,basin);
//        if(cycle_length<MaxCycleLength)
//            res[cycle_length]++;
        if(cycle_length==1)
            fixed_points++;
        //cout << "cycle length=" << cycle_length << endl;
    }
    n = attractors.size() ;
    result = 0;
    for(l=0; l<n; l++)
        result += attractors[l].size(); /* cycle length of the attractor */
    l_res = (1.0/n)*result;

    result = 0;
    for(l=0; l<n; l++)
        result += attractors[l].size()*basin[l];
    l_prime_res = (1.0/configs)*result; /* since we only have a fraction of the configurations this should work now*/
}

void saveCycles(char *fileName,uint64_t *cycle_length,int max){
    int l;
    FILE *file;
    file = fopen(fileName,"wb");
    if(file){
        for(l=0; l<max; l++){
            fprintf(file,"%d \t %d\n",l,cycle_length[l]);
        }
    }
    fclose(file);
}
/*generate a J with an symmetrical and antisymmetrical part, also k */
void setJSA(double J[N][N], double k,int x){
    int i,j;
    double Js[N][N]; //symmetric part
    double Ja[N][N]; //antisymmetric part
    srand(x + 111); /* seed for the J matrix*/
    /*first we fill the top right */
    for(i=0; i<N; i++){
        for(j=i; j<N; j++){
            Js[i][j] = gaussrand(1.0/(double_N*(1 + k*k)));
            Ja[i][j] = gaussrand(1.0/(double_N*(1 + k*k)));
        }
    }
    /*and then the bottom left */
    for(i=0; i<N; i++){
        for(j=0; j<i; j++){
            Ja[i][j] = Ja[j][i];
            Js[i][j] = -Js[j][i];
        }
    }
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            J[i][j] = Ja[i][j] + k*Js[i][j];
        }
    }
}
/* set J with just gaussian variables, x is the seed for srand() */
void setJSimple(double J[N][N], int x){
int i,j;
srand(x + 111); /* seed for the J matrix*/
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            J[i][j] = gaussrand(1/double_N);
            //printf("%lf  ",J[i][j]);
        }
    }
}
/* converts a given Sigma to a unique integer*/
int sigmaToInt(int16_t *Sigma){
    int i3;
    int res = 0 ;
    for(i3=0; i3<N; i3++){
        res += stepFunction(Sigma[i3])*pow(2,i3);
    }
    return res;
}

void saveNum(char *fileName,double l_avrg){
    int l;
    FILE *file;
    file = fopen(fileName,"wb");
    if(file){
            fprintf(file,"%d \t %lf\n",N,l_avrg );
    }
    fclose(file);
}

int main(){
int i,i2;
int total = 1000; /*total amount of different J's*/
int initial_configs = 1000; /*total amount of different configurations which are time evolved*/

double J[N][N];
char name[1024];
uint64_t lengths[MaxCycleLength]; /*cycles longer than MaxCycleLength steps are ignored*/
double l_avrg; /* average attractor size */
double l_prime_avrg; /* basin-size-weighted average attractor size*/
double l,l_prime;
double k = 0.5;
double fixed_points =0;
l_avrg = l_prime_avrg = 0;
for(i=0; i<MaxCycleLength; i++)
    lengths[i] = 0;
//srand(time(NULL)); /* seed for the J matrix*/

    for(i2=0; i2<total; i2++){
    /* progress indicator only works if total >= 100 */
        if (!(i2 % (total/100))){
            printf("\r");
            printf("%.3d%% Completed",100*i2/total);
        }
        //setJSimple(J,i2*2);
        setJSA(J,k,i2*2);
        //getRandomCycleLengths(J,initial_configs,lengths,l,l_prime);
        getAllCycleLengths(J,lengths,l,l_prime,fixed_points);
        l_avrg += l; l_prime_avrg += l_prime;
    }
    cout << endl;
    fixed_points /= (double)total;
    l_avrg /= (double)total; l_prime_avrg /= (double)total;
    cout << "log fixed_points=" <<log(fixed_points)<< endl;
//    cout << "l_avrg=" << l_avrg << "    l_prime_avrg=" << l_prime_avrg << endl;
//    sprintf(name,"res_N_cpp_%d_k=%lf.dat",N,k);
//    saveCycles(name,lengths,MaxCycleLength);
//    sprintf(name,"L_N_cpp_%d_k=%lf.dat",N,k);
//    saveNum(name,l_avrg);
//    sprintf(name,"L_prime_N_cpp_%d_k=%lf.dat",N,k);
//    saveNum(name,l_prime_avrg);
    sprintf(name,"log_fixed_points_%d_k=%lf.dat",N,k);
    saveNum(name,log(fixed_points));

}

