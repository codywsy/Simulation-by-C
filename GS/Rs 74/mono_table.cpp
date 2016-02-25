#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define k 5

int mono_order[20][400]; 



void main()
{
	int i,j;
	
	void mono_table();

	for(j=0;j<20;j++) 				
		for(i=0;i<400;i++)
			mono_order[j][i]=0;

	mono_table();
    

    for(i=0;i<3;i++)
    {
    	for(j=0;j<12;j++)
    	{
    		printf(" %d  ",mono_order[i][j]);
    	}

    	printf("\n");
    }
	
}



void mono_table(void)
{
	int i, j, v, l;

	j=0;	//represent row
	v=0;	//increasing counter for monomial order
	for(i=0;i<20;i++)
	{
		if(i<(k-1))	
		{
			mono_order[j][i]=v;
			v++;
		}
		else
		{
			l=floor((float)i/(float)(k-1));	
			for(j=0;j<=l;j++)
			{
				mono_order[j][i-(k-1)*j]=v;	//i-v*j
				v++;
			}
		}
	}
}