#include <stdio.h>

int main(){
    int i = 3;
    int a[2] = {0, 1};

    printf("%d\n", i);
    printf("%d\n", a[i=1]);
    printf("%d\n", i);

    return 0;
}