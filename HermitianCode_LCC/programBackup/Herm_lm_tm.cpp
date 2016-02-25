#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define N 500
#define tg_size 2500



int tg_order_1[tg_size][2],tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
int n,k,w;

	
void tgorder(void);
int gama(int);


void main(void)
{
	unsigned long int i,j,u,temp,u_temp,m,m_temp;
	unsigned long int gene,iter_num,weiz;
	unsigned long int lm,tm,tao,tao_GS;
	unsigned long int last_lm,last_tm,last_tao;
	float w_temp;

	//scanf w,k
	printf("Enter n:");
	scanf("%d",&n);
	printf("Enter k:");
	scanf("%d",&k);
	printf("Enter m:");
	scanf("%d",&m_temp);
	//**********

	w_temp=1.0/3;
	w=(int)pow( (float)n,w_temp );

	tgorder();

	last_lm=-1;
	last_tm=-1;
	last_tao=-1;

	//calculate gene, C=iter_num, weiz
	for(m=1;m<=m_temp;m++)
	{	
		gene=w*(w-1)/2;

		iter_num=n*(m+1)*m/2;

		weiz=tg_order[k-1][0]*w + tg_order[k-1][1]*(w+1);

		//calculate lm
		u_temp=-1;
		for(u=0;u<N;u++)
		{
			temp=weiz*u*(u-1)/2-(u-1)*gene;
			if(temp<=iter_num)
				u_temp=u;
			else if(temp>iter_num)
				break;
		}

		lm=u_temp-1;

		//calculate tm
		u_temp=-1;
		for(u=0;u<N;u++)
		{
			temp=(lm+1)*u + weiz*(lm+1)*lm/2 - lm*gene - gama(u);
			if(temp<=iter_num)
				u_temp=u;
			if(temp>iter_num)
				break;
		}

		tm=u_temp;

		//calculate tao, tao_GS
		tao=n-((lm*weiz+tm)/m)-1;

		tao_GS=n-k-gene+1;	//distance
		tao_GS=(int)sqrt( (float)n*(n-tao_GS) );
		tao_GS=n-tao_GS-1;
	
		if(tao!=last_tao)
		{
			last_lm=lm;
			last_tm=tm;
			last_tao=tao;

//			printf("\nw=%d\tg=%d\tC=%d\tweiz=%d",w,gene,iter_num,weiz);

			if(tao_GS!=tao)
				//printf all the element
				printf("\nm=%d\tlm=%d\ttm=%d\tTm=%d\t",m,lm,tm,tao);
			else if(tao_GS==tao)
				printf("\nm=%d\tlm=%d\ttm=%d\tTgs=%d",m,lm,tm,tao_GS);
			printf("\n");
		}

	}
}



void tgorder(void)
{
	int i,j;

	for(i=0;i<tg_size;i++)
	{
		tg_order_1[i][0]=-1;
		tg_order_1[i][1]=-1;
		tg_order[i][0]=-1;
		tg_order[i][1]=-1;
	}

	//total graduate order
	tg_order_1[0][0]=0;
	tg_order_1[0][1]=0;

	for(i=1;i<tg_size;i++)
	{
		if(tg_order_1[i-1][0]==0)
		{
			tg_order_1[i][0]=tg_order_1[i-1][1]+1;
			tg_order_1[i][1]=0;
		}
		else
		{
			tg_order_1[i][0]=tg_order_1[i-1][0]-1;
			tg_order_1[i][1]=tg_order_1[i-1][1]+1;
		}
	}

	//eliminate terms containing x^3
	j=0;
	for(i=0;i<tg_size;i++)
	{
		if(tg_order_1[i][0]<w+1)	//w+1
		{
			tg_order[j][0]=tg_order_1[i][0];
			tg_order[j][1]=tg_order_1[i][1];
			j++;
		}
	}
	/*
	for(i=0;i<100;i++)
		printf("\n(%d, %d)", tg_order[i][0], tg_order[i][1]);
	*/
}

int gama(int u)
{
	int i,j,num,temp;

	j=0;
	num=0;
	for(i=0;i<=u;i++)
	{
		temp=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		if(temp!=i)
			num++;
		else if(temp==i)
			j++;
	}

	return num;
}