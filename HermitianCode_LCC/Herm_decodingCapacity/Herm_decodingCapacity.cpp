#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void main(void)
{
	unsigned long int i,j,u,v,w,gene;
	float error,n,k,r,w_temp,T;


	FILE *fp;
	if((fp=fopen("result.txt","a"))==NULL)
	{
		printf("Can't open the data.txt\n");
		exit(0);
	}
	fclose(fp);

	printf("\nPlease enter n: ");
	scanf("%f", &n);

	w_temp=1.0/3;
	w=(int)pow( (float)n,w_temp );

	gene=w*(w-1)/2;

	printf("gene=%d\n\n",gene);

	for(r=0.00;r<=1.02;r=r+0.05)
	{
		k=n*r;

		T=n-k-(float)gene;

		if(T<0.0)
			T=0.0;

		T=T/2;
		error=T/n;

		printf("r=%2.2f, k=%2.2f, n=%2.2f, T=%2.2f, error=%E\n",r,k,n,T,error);

		fp=fopen("result.txt","a");
		fprintf(fp,"r=%2.2f, k=%2.2f, n=%2.2f, T=%2.2f, error=%E\n",r,k,n,T,error);
		fclose(fp);
	}

	getchar();
	getchar();

}