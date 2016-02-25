#include <stdio.h>

#define q 32
int mularray[]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,27,19,3,6,12,24,21,15,30,25,23,11,22,9,18}; 

void main(void)
{
	int i, j;

	for(i=1;i<q;i++)
	{
		for(j=0;j<q-1;j++)
			if( mularray[j] == i)
			{
				printf("\ntao_%d = %d\n", j, mularray[j]);
				break;
			}

	}



}