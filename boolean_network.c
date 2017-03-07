#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#define N 10
#define double_N 10.0

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
/* calculate the next time step and look for cycles */
int advanceTime(int8_t *Sigma, double J[N][N] ){
    int i,j,k;
    int8_t foundConfigs[(int)pow(2,N)][N];
    int configs = 1; //total amount of found configurations
    double h[N];

     for(i=0; i<N; i++) //save the initial configuration
        foundConfigs[0][i] = Sigma[i];

    while(configs < (int)pow(2,N) ){
        for(i=0; i<N; i++){
            h[i] = 0;
            for(j=0; j<N; j++){
                h[i] += J[i][j]*Sigma[j]; //equation (2)
            }
        }
        for(i=0; i<N; i++){
            foundConfigs[configs][i] = Sigma[i] = sign(h[i]); //equation (1)
            printf("%d  ",Sigma[i]);
        }
        printf("\t %d \n",configs);
        /*check if the newly found configuration Sigma matches a previous one */
        for(k=configs-1; k>=0; k--){
            if(checkIfEqual(foundConfigs[k],Sigma))
                return configs - k; //cycle length
        }
        configs++;
    }
    return -1;
}


int main(){
int i,j;
double J[N][N];
int8_t Sigma[N]; //single byte signed int since it's only 1 or -1 anyway
int cycle_length;

srand(1); /* seed for the J matrix*/
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            J[i][j] = gaussrand(1/double_N);
            //printf("%lf  ",J[i][j]);
        }
    }

/* TODO: check every possible initial condition and then histogram the cycle lengths*/
srand(11); /* seed for the initial condition */
    printf("Sigma=");
    for(i=0; i<N; i++){
        Sigma[i] = sign( rand()-RAND_MAX*0.5 );
        printf("%d ",Sigma[i]);
    }
    printf("\t 0 \n");

cycle_length = advanceTime(Sigma,J);
printf("cycle length = %d \n",cycle_length);
}
