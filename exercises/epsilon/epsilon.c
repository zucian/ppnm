#include<stdio.h>
#include<math.h>
#include<limits.h>
#include<float.h>

int main(){
    printf("exercise epsilon \n");

    printf("1)\n");
    int i=1;
    while(i+1>i){
        i++;
    }
    printf("my INT_MAX = %i \n", INT_MAX);
    printf("my max int while = %i\n",i);


    i =1;
    do {
        i++;
    } while (i+1 > i);
    printf("my max int do = %i\n",i);

    i=1;
    while(i-1<i){
        i--;
    }
    printf("my INT_MIN = %i \n", INT_MIN);
    printf("my min int while = %i\n",i);

    i =1;
    do {
        i--;
    } while (i-1 < i);
    printf("my min int do = %i\n",i);

    double x=1;
    while(1+x!=1){
        x/=2;
    }
    x*=2;

    double e;
    for(e=1; 1+e!=1; e/=2){}
    e*=2;

    float fx=1;
    while(1+fx!=1){
        fx/=2;
    }
    fx*=2;

    float fe;
    for(fe=1; 1+fe!=1; fe/=2){}
    fe*=2;

    long double lx=1;
    while(1+fx!=1){
        fx/=2;
    }
    fx*=2;

    long double le;
    for(le=1; 1+le!=1; le/=2){}
    le*=2;

    printf("Calculated eps double while = %g\n",x);	
    printf("Calculated eps double for = %g\n",e);
    printf("Double eps = %g\n", DBL_EPSILON);
    printf("Calculated eps float while = %g\n",fx);	
    printf("Calculated eps float for = %g\n",fe);
    printf("Float eps = %g\n", FLT_EPSILON);
    printf("Calculated eps long double while = %Lg\n",lx);	
    printf("Calculated eps long double for = %Lg\n",le);
    printf("Long Double eps = %Lg\n", LDBL_EPSILON);

    int max =INT_MAX/2;
    float sum_up_float = 0;
    float sum_down_float = 0;
    for(int i=1; i<=max; i++){
        sum_up_float += 1.0f/i;
    }
    for(int i=0; i<max; i++){
        sum_down_float += 1.0f/(max-i);
    }

    printf("max = %i \n",max); 
	printf("sum_up_float=%f\n ",sum_up_float);
	printf("sum_down_float= %f \n", sum_down_float);

    double sum_up_double = 0;
    double sum_down_double = 0;
    for(int i=1; i<=max; i++){
        sum_up_double += 1.0f/i;
    }
    for(int i=0; i<max; i++){
        sum_down_double += 1.0f/(max-i);
    }

    printf("sum_up_double=%f\n ",sum_up_double);
    printf("sum_down_double= %f \n", sum_down_double);


    return 0;
}