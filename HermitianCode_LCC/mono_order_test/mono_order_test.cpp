#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define w 8
#define weiz 316
#define lm 1

//tgorder()
#define tg_size 4600	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 17
#define weight_XYsize 4600	//equals to the tg_size
#define mono_ordinSize (2*4600)	//choose the value with debugging
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 506	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 4600 //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size


//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 460	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize (lm+1)// same as the interpoly_Zsize
#define facpoly_Ysize 460		// same as the interpoly_Ysize
#define facpoly_Xsize (w+1)	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 505	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 17	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 36	//change with diff w and k, the max degY of probably used polebasis(polebasis_k-1)
#define faiMax_Xsize 8	//change with diff w and k, the max degX of probably used polebasis(polebasis_k-1)
//expoly()-->expanded polynomial
#define expoly_Ysize (faiMax_Ysize+1)	
#define expoly_Xsize (faiMax_Xsize+1)



void tgorder();
void mono_table();

int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
int tg_order[tg_size][2];
int k,n;


void main()
{
	int i;
/*
	printf("please enter n:");
	scanf("%d",&n);
	printf("please enter k:");
	scanf("%d",&k);
*/
	n=512;
	k=289;

	tgorder();
 
	mono_table();

	getchar();
	getchar();

}

void tgorder()
{
	int i, j, u, index_temp;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 515;

	j=0;
	for(i=0;i<=index_temp;i++)
		for(u=i;u>=0;u--)
			if(u<=w)
			{
				tg_order[j][0]=u;	// 0<deg_x<=w
				tg_order[j][1]=i-u;	// 0<deg_y
				j++;
			}

/*	//*****debug******
	for(i=0;i<tg_size;i++)
		printf("\ntg[%d]=(%d, %d)",i,tg_order[i][0], tg_order[i][1]);

	printf("\n\n");
*/	//****************

}

void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[weight_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]

	//Initialisation
	for(i=0;i<weight_Zsize;i++)	//(230/weight(z))+1
	{
		for(j=0;j<weight_XYsize;j++)	//pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][y_size][w+1].
		{
			weight[i][j]=-1;
			mono_order_1[i][j]=-1;
		}
	}

	for(i=0;i<monoTable_Zsize;i++)	// rs
		for(j=0;j<monoTable_Ysize;j++)	//max(deg_y)+1
			for(u=0;u<monoTable_Xsize;u++)	//w+1
				mono_order[i][j][u]=-1;




	//weight of monomial over function x^3-y^2-y=0
	for(i=0;i<weight_Zsize;i++)	
	{
		if(i==0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]=tg_order[j][0]*w + tg_order[j][1]*(w+1);
		else if(i>0)
			for(j=0;j<weight_XYsize;j++)
				weight[i][j]= weiz*i + weight[0][j];
	}


/*
	for(i=0;i<weight_Zsize;i++)	
	{
		weight[i][0]=weiz*i;
		for(j=1;j<weight_XYsize;j++)
			weight[i][j]=weiz*i+j+1;
	}
*/
	//construction 2-dimensional mono table
	z=0;
	for(v=0;v<mono_ordinSize;v++)	//for each possible weight until weight(term[rs-1][max(deg_y)+1][w-1])+1, note term[rs-1][max(deg_y)+1][w-1] is the term next to term[rs-1][max(deg_y)][w]
	{
		for(i=0;i<weight_Zsize;i++)	//230/weight(z)
		{
			for(j=0;j<weight_XYsize;j++)	//230, >the index of the largest term with deg(z)=0, and weight=weight(term[rs-1][max(deg_y)+1][w-1])
			{
				if(weight[i][j]==v)
				{
					mono_order_1[i][j]=z;
					z++;
					break;
				}
			}
		}
	}

	// 2-dimensional table transform into 3-dimensional table
	for(u=0;u<monoTable_Zsize;u++)	//rs
		for(z=0;z<monoTable_totalSize;z++)	//>index of term[max(deg_y)][w] in the pole basis + 1
			mono_order[u][tg_order[z][1]][tg_order[z][0]]=mono_order_1[u][z];

	//**********debug********
	for(j=0;j<weight_XYsize;j++)
		printf("\t%d",j);
	printf("\n");
	for(i=0;i<weight_Zsize;i++)
	{
		printf("\n");
		for(j=0;j<weight_XYsize;j++)
			printf("\t%d",weight[i][j]);
	}

	printf("\n\n\nmono_order_1\n");
	for(i=0;i<weight_Zsize;i++)
	{
		printf("\n");
		for(j=0;j<weight_XYsize;j++)
			printf("\t%d",mono_order_1[i][j]);
	}	
	
	printf("\nMonmial basis is:\n");
	for(i=0;i<monoTable_Zsize;i++)
	{
		printf("\n\nZ=%d",i);

		for(z=0;z<monoTable_Xsize;z++)
			printf("\t%d ",z);
		printf("\n\n");
		for(j=0;j<monoTable_Ysize;j++)
		{
			printf("\n%d",j);
			for(z=0;z<monoTable_Xsize;z++)
				printf("\t%d ", mono_order[i][j][z]);
		}
	}

/*	printf("\n\nweight degree for rcs():");
	for(j=0;j<interpoly_Ysize;j++)
	{
		printf("\n%d",j);
		for(i=0;i<interpoly_Xsize;i++)
			printf("\t%d",weiDeg_rcs[j][i]);
	}

	printf("\n\nweight degree index for rcs():");
	for(j=0;j<interpoly_Ysize;j++)
	{
		printf("\n%d",j);
		for(i=0;i<interpoly_Xsize;i++)
			printf("\t%d",weiDeg_rcs_index[j][i]);
	}
*/
	//******************

}