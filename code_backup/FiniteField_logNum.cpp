#include <stdio.h>

int mularray[]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,27,19,3,6,12,24,21,15,30,25,23,11,22,9,18}; 

void main(void)
{
	int i, j;

	for(i=1;i<n;i++)
	{
		for(j=0;j<n-1;j++)
			if( mularray[j] == i)
				printf("\n\ntao_%d = %d", j, mularray[j]);

	}



}