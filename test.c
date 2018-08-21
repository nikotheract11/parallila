#include <stdio.h>

static inline int foo(){
	printf("Hello inline!\n");
}

int main(){
	foo();
	return 0;
}
