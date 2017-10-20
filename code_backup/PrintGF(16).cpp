#include <stdio.h>
#include <string.h>

#define q 16
#define p 4


int main(void)
{
	int i,j;
	int root[q], mularray[q-1], logNum[q];
	
	//initialization
	memset(root,-1,sizeof(root));
	memset(mularray,-1,sizeof(mularray));
	memset(logNum,-1,sizeof(logNum));

	//generate 
	//root[q]
	for(i=0;i<q;i++)
		root[i] = i;

	printf("\n\nroot[]=\n");
	for(j=0;j<q;j++)
	{
		printf("\t%d", root[j]);
	}
	printf("\n\n");

	//mularray[q-1]
	int poly[q];
	memset(poly,0,sizeof(poly));

	for(j=0;j<q-1;++j)
	{
		int temp = 0;
		unsigned int mask = 1;

		memset(poly,0,sizeof(poly));
		poly[j]=1;
		for(i=q-1;i>=p;i--)
			if(poly[i]!=0)
			{	//primitive polynomial x^4=x^2+1
				poly[i-3] = (poly[i-3]+1)%2;
				poly[i-4] = (poly[i-4]+1)%2;
				poly[i] = 0;
			}

		for(i=0;i<p;i++)
		{
			if(poly[i]!=0)
			{
				temp += mask;
			}

			mask = mask<<1;
		}

		mularray[j] = temp;
	}
		
	printf("\n\nmularray[]=\n");
	for(j=0;j<q-1;j++)
	{
		printf("\t%d", mularray[j]);
	}
	printf("\n\n");

	//logNum[q]
	for(i=0;i<q-1;i++)
	{
		for(j=0;j<q;j++)
			if(mularray[i]==j)
			{
				logNum[j]=i;
			}
	}

	printf("\n\nlogNum[]=\n");
	for(j=0;j<q;j++)
	{
		printf("\t%d", logNum[j]);
	}
	printf("\n\n");

	return 0;
}