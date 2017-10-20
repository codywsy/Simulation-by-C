#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define N 500
#define tg_size 512



int tg_order[tg_size][2]={0};	//[promise the tg order beyond (w, max(deg_y))][2]
int n,k;
int w;

	
void tgorder(void);
int gama(int);


void main(void)
{
	unsigned long int i,j,u,u_temp,m,m_temp;
	unsigned long int gene,iter_num,weiz;
	unsigned long int lm,tm,tao,tao_GS;
	unsigned long int last_lm,last_tm,last_tao;
	float w_temp,temp;

	//scanf w,k
	printf("Enter n:");
	scanf("%d",&n);
	printf("Enter k:");
	scanf("%d",&k);
//	printf("Enter m:");
//	scanf("%d",&m_temp);
	m_temp=2;
	//**********
	//calculate w
	w=-1;
	w_temp=1.0/3;
	temp = pow( (float)n,w_temp );
	w = (int)temp;

	last_lm=-1;
	last_tm=-1;
	last_tao=-1;

	tgorder();

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
			temp=((float)(weiz*u)/2.0-gene)*(u-1);
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
			temp=(lm+1)*u + weiz*(lm+1)*lm/2.0 - lm*gene - gama(u);
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
	
		printf("\nw=%d\tg=%d\tC=%d\tweiz=%d\nm=%d\tlm=%d\ttm=%d\tTm=%d\t",w,gene,iter_num,weiz,m,lm,tm,tao);

/*
		if(tao!=last_tao)
		{
			last_lm=lm;
			last_tm=tm;
			last_tao=tao;

			printf("\nw=%d\tg=%d\tC=%d\tweiz=%d",w,gene,iter_num,weiz);

			if(tao_GS!=tao)
				//printf all the element
				printf("\nm=%d\tlm=%d\ttm=%d\tTm=%d\t",m,lm,tm,tao);
			else if(tao_GS==tao)
				printf("\nm=%d\tlm=%d\ttm=%d\tTgs=%d",m,lm,tm,tao_GS);
			printf("\n");
		}
*/
	}

	printf("\n");
	getchar();getchar();

}



void tgorder(void)
{
	int i, j, u, index_temp;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 59;

	j=0;
	for(i=0;i<=index_temp;i++)
		for(u=i;u>=0;u--)
			if(u<=w)
			{
				tg_order[j][0]=u;	// 0<deg_x<=w
				tg_order[j][1]=i-u;	// 0<deg_y
				j++;
			}
	
	//*****debug******
	for(i=0;i<k;i++)
		printf("\ntg[%d]=(%d, %d)",i,tg_order[i][0], tg_order[i][1]);
	//****************

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