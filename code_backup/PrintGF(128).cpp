#include <stdio.h>
#include <string.h>

#define q 128
#define p 7


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
		printf("%d,", root[j]);
	}
	printf("\n\n");

	FILE * fp;
	if((fp = fopen("root.txt","a")) == NULL)
	{
		printf("\nopen root.txt failed\n");
	}
	else
	{
		for(j=0;j<q;j++)
		{
			fprintf(fp, "%d,", root[j]);
		}

		fclose(fp);
	}

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
			{	//primitive polynomial x^7=x^3+1
				poly[i-4] = (poly[i-4]+1)%2;
				poly[i-7] = (poly[i-7]+1)%2;
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
		printf("%d,", mularray[j]);
	}
	printf("\n\n");

	if((fp = fopen("mularray.txt","a")) == NULL)
	{
		printf("\nopen mularray.txt failed\n");
	}
	else
	{
		for(j=0;j<q-1;j++)
		{
			fprintf(fp, "%d,", mularray[j]);
		}

		fclose(fp);
	}

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
		printf("%d,", logNum[j]);
	}
	printf("\n\n");

	if((fp = fopen("logNum.txt","a")) == NULL)
	{
		printf("\nopen logNum.txt failed\n");
	}
	else
	{
		for(j=0;j<q;j++)
		{
			fprintf(fp, "%d,", logNum[j]);
		}

		fclose(fp);
	}

	return 0;
}