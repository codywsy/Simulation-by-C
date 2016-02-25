/******************
function: print out the mono_table

var: (n,k)=(15,9)fg
********************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define n 15	//length of codeword
//#define k 9 //length of message
#define x_size 200
#define y_size 40

int n,k;
int m; // zero of multiplicity
int start,finish; 
long int mono_order[y_size][x_size]; // the size of monoTable is varible
int iterNum; // C
int tm;
int lm;



void mono_table(void);


void main(void)
{
	int i,j,temp;
	int para[3];  
	//para[0]--iterNum   para[1]--tm   para[2]--lm 

	printf("Enter n: ");
	scanf("%d",&n);
	printf("Enter k: ");
	scanf("%d",&k);

	mono_table();

	for ( i = 0; i < 60 ; ++i)
		printf("%d\t",i);
	printf("\n\n");
		
	

	printf("(%d,%d)'s mono_table:\n",n,k);
	for(j=0;j<y_size;j++)
	{
		for(i=0;i<x_size;i++)
			printf("%d\t",mono_order[j][i]);

		printf("\n\n");
	}

	for(i=0;i<3;i++)
		para[i]=0;
	temp=0;
	printf("Enter finish_m:");
	scanf("%d",&j);

	for(m=1;m<=j;m++)
	{	
		//cal iterNum
		para[0] = n*(m+1)*m/2;

		//cal tm
		i=0;
		while( mono_order[0][i] <= para[0] )
			i++;
		para[1]=n-1-(i-1)/m;

		//cal lm
		i=0;
		while( mono_order[i][0] <= para[0] )
			i++;
		para[2]=i-1;

		// if( para[1]>temp )
		// {
			// temp=para[1];
			printf("\n\nm=%d\titerNum=%d\ttm=%d\tlm=%d",m,para[0],para[1],para[2]);
		// }
		// else if(para[1]==-1)
		// {	
			// printf("\nerror");
			// break;
		// }
	
	}

	getchar();
	getchar();

}




void mono_table(void)
{
	long int i, j, v, l;

	for(j=0;j<y_size;j++)
		for(i=0;i<x_size;i++)
			mono_order[j][i]=-1;

	j=0;	//represent row
	v=0;	//increasing counter for monomial order
	for(i=0;i<x_size;i++)
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



