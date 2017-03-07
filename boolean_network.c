#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#define MaxCycleLength 1000
#define N 13
#define double_N 13.0

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

int sign(double x){
    if(x >= 0) return 1;
    if(x < 0) return -1;
}
/*check if 2 arrays match element by element */
int checkIfEqual(int8_t *array1, int8_t *array2){
int ii;
    for(ii = 0; ii < N; ii++) {
        if (array1[ii] != array2[ii]) return 0;
    }
    return 1;
}
/* do the time evolution and look for a cycle (should I check a cycle by letting it run a little more?) */
int advanceTime(int8_t *Sigma, double J[N][N] ){
    int i,j,k;
    int8_t foundConfigs[N*100][N]; //It usually doesn't take that many steps to find a cycle I found 217 at most for N=16 so N*100 should be more than enough
    int configs = 1; //total amount of found configurations
    double h[N]; //the input equation (2)

     for(i=0; i<N; i++) //save the initial configuration
        foundConfigs[0][i] = Sigma[i];

    while(configs < (int)pow(2,N) ){ //for safety while(true) seems dangerous
        for(i=0; i<N; i++){
            h[i] = 0;
            for(j=0; j<N; j++){
                h[i] += J[i][j]*Sigma[j]; //equation (2)
            }
        }
        for(i=0; i<N; i++){
            foundConfigs[configs][i] = Sigma[i] = sign(h[i]); //equation (1)
            //printf("%d  ",Sigma[i]);
        }
        //printf("\t %d \n",configs);
        /*check if the newly found configuration Sigma matches a previous one, if so we have a cycle (really? that's good enough?) */
        for(k=configs-1; k>=0; k--){
            if(checkIfEqual(foundConfigs[k],Sigma)){
                return configs - k; //cycle length
            }
        }
        configs++;
    }

    return -1;
}
/* going from n=0 to n=2^N goes through every possible initial configuration*/
void setSigma(int8_t *Sigma,int n){
    int i,j;
    for(i=0; i<N; i++)
        Sigma[i] = -1;
    for(j=0; j<n; j++){
        Sigma[0] += 2;
        for (i=0; i<N-1; i++){
            if(Sigma[i] > 1 ){
                Sigma[i] = -1;
                Sigma[i+1] += 2;
            }
        if(Sigma[N-1]>1)
            Sigma[N-1] = -1;
        }
    }
}
/* too many times have I had to add this manually so now it's a function */
void printSigma(int8_t *Sigma){
        int i1;
        printf("Sigma=");
        for(i1=0; i1<N; i1++)
            printf("%d ",Sigma[i1]);
        printf("\n");
}
/* goes through every single initial configuration and puts the cycle lengths in the res array*/
void getAllCycleLengths(char *fileName,int8_t *Sigma, double J[N][N],uint64_t *res){
    int l;
    int cycle_length;
    for(l=0; l<pow(2,N); l++){
        setSigma(Sigma,l);
        cycle_length = advanceTime(Sigma,J);
        if(cycle_length<MaxCycleLength)
            res[cycle_length]++;
    }
}
/* goes through n random initial configurations and puts the cycle lengths in the res array*/
void getRandomCycleLengths(char *fileName,int8_t *Sigma, double J[N][N],int n,int x,uint64_t *res){
    int l;
    int l1;
    int cycle_length;
    for(l1=0; l1<n; l1++){
        srand(x+l1);
        for(l=0; l<N; l++)
            Sigma[l] = sign( rand()-RAND_MAX*0.5 );
            cycle_length = advanceTime(Sigma,J);
            if(cycle_length<MaxCycleLength)
                res[cycle_length]++;
    }
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

int main(){
int i,j,i2;
int total = 10000; //total amount of different J's

double J[N][N];
int8_t Sigma[N]; //single byte signed int since it's only 1 or -1 anyway
char name[1024];

uint64_t lengths[MaxCycleLength]; //cycles longer than MaxCycleLength steps are ignored
    for(i=0; i<MaxCycleLength; i++)
        lengths[i] = 0;

    for(i2=0; i2<total; i2++){
        if (!(i2 % (total/100))){
            printf("\r");
            printf("%.3d%% Completed",100*i2/total);
        }
        srand(i2 + 111); /* seed for the J matrix*/
            for(i=0; i<N; i++){
                for(j=0; j<N; j++){
                    J[i][j] = gaussrand(1/double_N);
                    //printf("%lf  ",J[i][j]);
                }
            }
        sprintf(name,"dump%d.dat",i2);
        getRandomCycleLengths(name,Sigma,J,1000,i2,lengths);
        //getAllCycleLengths(name,Sigma,J,lengths);
    }
    sprintf(name,"res_N_%d.dat",N);
    saveCycles(name,lengths,100);
}

