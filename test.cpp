#include <stdio.h>

int main(int argc, char *argv[]){
    double a, b, c;

    while(scanf("%lf%lf", &a, &b) != EOF)
        printf("%lf, %lf, %lf\n", a, b, c);

    return 0;
}