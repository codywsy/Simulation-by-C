/***********************************
program: 求解相关系数
statement: 核心模块是coefficientSearch()，输入输
出是coeffTable表，与zerobasis()，mono_table(),cal_max_a()存在耦合
关系.
*********************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//main
#define q 4	//GF(q)
#define p 2	//2^p
#define n 8
#define k 4
#define w 2
#define weiz 4

//*****need to be modify**************
//coefficientSearch()
#define max_m 3	//design large enough to cover 
#define max_alpha 8	//alpha < max_m, designing alpha large enough to cover all the m_ij of Mulplicity Matrix
#define tableSize_a 20	//more than max_a
#define tableSize_alpha (max_alpha+30) 	//more than max_alpha
//zerobasie()
#define zb_alpha 20	//size+1, larger than tableSize_alpha, according to the maxDeg_x and maxDeg_y of PbMax, we can judge the max num Zb we need
#define zb_Ysize (zb_alpha/(w+1)+1)	//this term must be changed when alpha greater than (w+1)
#define zb_Xsize (w+1)
//*****************************************


//tgorder()
#define tg_size 275	//tg_size represent the probably used pole basis num
//mono_table()
#define weight_Zsize 4
#define weight_XYsize 275
#define mono_ordinSize 100
#define monoTable_Zsize	2	//equals to the interpoly_Zsize, lm+1
#define monoTable_Ysize 67	//large than the interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize 200 //>index of term[max(deg_y)][w] in the pole basis + 1
//findpoint()

//cal_max_a()
#define N 500
//polebasis()
//#define pbNum 275	//num of pb
//#define pb_Ysize 52	//large than max(degY) + 1
//#define pb_Xsize (w+1)	//w+1




//funciton statement
void findpoints();
void tgorder();
void mono_table();
int cal_max_a();
int gama(int);
//void polebasis();
void zerobasis(int, int LM_zb[][zb_alpha]);
void coefficientSearch();
int mul(int,int);
int add(int,int);
int power(int,int);


//main
int mularray[]={1,2,3};
int root[]={0, 1, 2, 3};	//elements in GF(4)
//tgorder()
int tg_order[tg_size][2];	//[promise the tg order beyond (w, max(deg_y))][2]
//mono_table()
int mono_order[monoTable_Zsize][monoTable_Ysize][monoTable_Xsize];	//mono_order[z-deg+1][y-deg][x-deg+1], y-deg is greater than w-1+nw/(w+1)
//cal_max_a()

//polebasis()
//int pb[pbNum][pb_Ysize][pb_Xsize];	//pole basis[num_pb][y_size][x_size]
//findpoint()
int point[2][n];	//rational points[2][w^3]
//zerobasis
int zb[zb_alpha][zb_Ysize][zb_Xsize];
//coefficientSearch
int coeff_table[n][tableSize_alpha][tableSize_a];	//n points


void main(void)
{
	int i, u, v;

	//initilaization
	for(i=0;i<n;i++)
		for(u=0;u<tableSize_alpha;u++)
			for(v=0;v<tableSize_a;v++)
			{
				coeff_table[i][u][v] = -1;
			}

	findpoints();
	tgorder();
	mono_table();

	coefficientSearch();


}



void findpoints()
{
	int i, j, u, temp1, temp2, temp3, x, y;

	//Initialisation
	for(i=0;i<2;i++)
		for(j=0;j<n;j++)	//w^3
			point[i][j]=0;

	//find points over x^(w+1)-y^w-y=0
	u=0;
	for(i=0;i<q;i++)	//number of elements in GF(q)
	{
		x=root[i];
		for(j=0;j<q;j++)	//number of elements in GF(q)
		{
			y=root[j];
			temp1=power(x, w+1);
			temp2=power(y, w);
			temp3=y;
			if(add(temp3, add(temp1, temp2))==0)
			{
				point[0][u]=x;
				point[1][u]=y;
				u++;
			}
		}
	}

}

void tgorder()
{
	int i, j, u, index_temp;

	//judge the index's scale of coresponding tgsize
	/*pesudo code
	if(  0 <= index < w ) tgsize = (index+2) * (index+1) / 2
	else if( w <= index ) tgsize = w*(w+1)/2 + (w+1)*(index-w+1)
	*/
	index_temp = 92;	//according to the last formualtion to calculate out

	j=0;
	for(i=0;i<=index_temp;i++)
		for(u=i;u>=0;u--)
			if(u<=w)
			{
				tg_order[j][0]=u;	// 0<deg_x<=w
				tg_order[j][1]=i-u;	// 0<deg_y
				j++;
			}

}


void mono_table()
{
	int i, j, u, v, z, weight[weight_Zsize][weight_XYsize], mono_order_1[weight_Zsize][weight_XYsize];	//26=(100/weight(z))+1, 100=pre-set size, promise term[rs-1][max(deg_y)][w] is convered, so does every term in array[rs][max(deg_y)+1][w+1]
	int index_x, index_y;

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
		{
			index_x = tg_order[z][0];
			index_y = tg_order[z][1];
			mono_order[u][index_y][index_x]=mono_order_1[u][z];
		
		}
/*	//**********debug********
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
*/	//******************

}

int cal_max_a()
{
	unsigned long int u, gene, iter_num_temp, u_temp, deg_Q, a_temp;
	float temp;
	int lm_temp,tm;

	//calculate gene, C=iter_num, 
	gene=w*(w-1)/2;

	iter_num_temp=n*(max_m+1)*max_m/2;

	//calculate lm
	u_temp=-1;
	for(u=0;u<N;u++)
	{
		temp=((float)(weiz*u)/2.0-gene)*(u-1);
		if(temp<=iter_num_temp)
		{
			u_temp=u;
		}
		else if(temp>iter_num_temp)
		{
			break;
		}
	}

	lm_temp=u_temp-1;

	//calculate tm
	u_temp=-1;
	for(u=0;u<N;u++)
	{
		temp=(lm_temp+1)*u + weiz*(lm_temp+1)*lm_temp/2.0 - lm_temp*gene - gama(u);
		if(temp<=iter_num_temp)
		{	
			u_temp=u;
		}
		if(temp>iter_num_temp)
		{
			break;
		}
	}

	tm=u_temp;

	//calculate deg_Q
	deg_Q = lm_temp*weiz + tm;

	//search a
	a_temp=-1;
	for(u=0;u<tg_size;u++)
	{
		u_temp= tg_order[u][0]*w + tg_order[u][1]*(w+1);
		if( u_temp>=deg_Q )	//key part, determine can be "<" or "<="
		{
			a_temp=u;
			break;
		}
	}

	if(a_temp==-1)
	{
		printf("\n\ncalculate a is erro\n");
	}

	return a_temp;

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


//if(alpha<(w+1))
void zerobasis(int pointIndex, int LM_zb[][zb_alpha])
{
	int i,j,u,v;
	int temp,flag,xi,yi,index_flag,temp2;
	int zb_temp1[zb_Xsize], zb_temp2[zb_Ysize][zb_Xsize+1];
	int poly_temp1[zb_Ysize+w][zb_Xsize+w],poly_temp2[zb_Ysize+1][zb_Xsize];

	//Initialisation
	for(i=0;i<zb_alpha;i++)
		for(u=0;u<zb_Ysize;u++)		
			for(v=0;v<zb_Xsize;v++)
			{
				zb[i][u][v] = 0;
			}

	xi = point[0][pointIndex];
	yi = point[1][pointIndex];
	flag = 0;

	//calculate Zb
	for(i=0;i<zb_alpha;i++)
	{

		//set poly_temp1 to be 0
		for(u=0;u<zb_Ysize+w;u++)
			for(v=0;v<zb_Xsize+w;v++)
			{
				poly_temp1[u][v]=0;				
			}
		for(u=0;u<zb_Ysize+1;u++)
			for(v=0;v<zb_Xsize;v++)
			{
				poly_temp2[u][v]=0;
			}

		//zb_a
		temp=i%(w+1);
		if(temp==0)
		{
			//initialization
			for(v=0;v<zb_Xsize;v++)
			{
				zb_temp1[v]=0;
			}

			zb_temp1[0]=1;
			flag=1;	//change the zb_temp2
		}
		else if(temp==1)
		{
			zb_temp1[0]=xi;
			zb_temp1[1]=1;
		}
		else if(temp>1 && temp<=w)
		{
			//x*zb
			poly_temp1[0][0]=0;	//only the first line of poly_temp1 will be uesd
			for(v=0;v<zb_Xsize;v++)
			{
				poly_temp1[0][v+1]=zb_temp1[v];
			}

			//xi*zb+x*zb
			for(v=0;v<zb_Xsize;v++)
			{
				zb_temp1[v]=mul( xi,zb_temp1[v] );
				zb_temp1[v]=add( zb_temp1[v],poly_temp1[0][v] );				
			}

		}


		//zb_b;
		if(flag==1)
		{
			temp=i/(w+1);
			if(temp==0)
			{
				//initialization
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize+1;v++)
					{
						zb_temp2[u][v]=0;
					}

				//calculate zb
				zb_temp2[0][0] = 1;
			}
			else if(temp==1)
			{
				zb_temp2[1][0] = 1;
				zb_temp2[0][1] = power( xi,w );
				zb_temp2[0][0] = mul( zb_temp2[0][1],xi );
				zb_temp2[0][0] = add( zb_temp2[0][0],yi );
			}
			else if(temp>1)
			{
				//cal xi^w*x*zb_temp2
				temp2 = power( xi,w );
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp1[u][v+1]=mul( temp2,zb_temp2[u][v] );
						}
				//x^(w+1)=y^w+y
				for(u=0;u<zb_Ysize;u++)
					if( poly_temp1[u][w+1]!=0 )
					{
						poly_temp1[u+1][0] = add( poly_temp1[u+1][0],poly_temp1[u][w+1] );	//y^1
						poly_temp1[u+w][0] = add( poly_temp1[u+w][0],poly_temp1[u][w+1] );	//y^2
						poly_temp1[u][w+1] = 0;
					}

					
				//cal y+zb_temp2
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp2[u+1][v]=zb_temp2[u][v];
						}

				//real zb_temp2
				temp = power( xi,(w+1) );
				temp = add( temp,yi );
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						{
							zb_temp2[u][v] = mul( zb_temp2[u][v],temp );
							zb_temp2[u][v] = add( zb_temp2[u][v],poly_temp1[u][v] );
							zb_temp2[u][v] = add( zb_temp2[u][v],poly_temp2[u][v] );
						}
			}

			//close the changing of zb_temp2
			flag=0;
		}

		//zb = zb_temp1 * zb_temp2
		for(j=0;j<zb_Xsize;j++)
			if(	zb_temp1[j]!=0 )
			{

				//set poly_temp1 to be 0
				for(u=0;u<zb_Ysize+w;u++)
					for(v=0;v<zb_Xsize+w;v++)
					{
						poly_temp1[u][v]=0;
					}
				
				//cal zb_temp1[j]*zb_temp2
				temp = zb_temp1[j];
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb_temp2[u][v]!=0 )
						{
							poly_temp1[u][v+j] = mul( temp,zb_temp2[u][v] );
						}

				//x^(w+1)=y^w+y
				index_flag = zb_Xsize+j-1;
				while( index_flag>w )
				{
					temp = index_flag-(w+1);
					for(u=0;u<zb_Ysize;u++)
						if( poly_temp1[u][index_flag]!=0 )
						{
							poly_temp1[u+1][temp] = add( poly_temp1[u+1][temp],poly_temp1[u][index_flag] );	//y^1
							poly_temp1[u+w][temp] = add( poly_temp1[u+w][temp],poly_temp1[u][index_flag] );	//y^2
							poly_temp1[u][index_flag] = 0;
						}
					index_flag = index_flag-1;
				}

				
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb[i][u][v]!=0 || poly_temp1[u][v]!=0 )
						{
							zb[i][u][v] = add( zb[i][u][v],poly_temp1[u][v] );
						}				
			}

		//calculate LM_zb and LC_zb
		LM_zb[0][i] = i%(w+1);
		LM_zb[1][i] = i/(w+1);
	}
}


void coefficientSearch(void)
{
	int i, j, u, v, z, flag, temp, temp2;
	int poly[zb_Ysize][zb_Xsize];
	int max_a,u_temp, v_temp, index_x, index_y, choosen_index, lod_temp;
	int LM[2][zb_alpha];
	int effTable[n][tableSize_alpha][tableSize_a];

	//Initialisazion
	for(i=0;i<n;i++)
		for(u=0;u<tableSize_alpha;u++)
			for(v=0;v<tableSize_a;v++)
			{
				effTable[i][u][v]=0;
			}

	for(u=0;u<2;u++)
		for(v=0;v<zb_alpha;v++)
		{
			LM[u][v]=0;
		}

	max_a = cal_max_a();
	if( max_a >= tableSize_a )
	{
		printf("\n\nmax_a = %d, tableSize_a should be set larger\n\n", max_a);
	}


	//calculate coeff_table
	for(i=0;i<n;i++)
	{
		//initialization poly
		for(u=0;u<zb_Ysize;u++)
			for(v=0;v<zb_Xsize;v++)
			{
				poly[u][v]=0;
			}

		zerobasis(i,LM);	//calculate all the zb which index under tableSize_alpha and its leading monomial

		for(j=max_a;j>=0;j--)
		{
			//initialize coefficient
			for(z=0;z<(max_alpha+1);z++)
			{
				effTable[i][z][j]=0;
			}
			
			//search the choosen zb
			index_x = tg_order[j][0];
			index_y = tg_order[j][1];

			//cal the first coeffTable value
			for(z=zb_alpha-1;z>=0;z--)
				if( LM[0][z]==index_x && LM[1][z]==index_y )
				{
					choosen_index = z;
					if( z<tableSize_alpha )
					{
						effTable[i][z][j]=1;
					}
					else
					{
						printf("\n\n0.tableSize_alpha is not large enough\n\n");
					}
					
					break;
				}

			//initialize the first poly 
			for(u=0;u<zb_Ysize;u++)
				for(v=0;v<zb_Xsize;v++)
					if( zb[choosen_index][u][v]!=0 )
					{
						poly[u][v]=zb[choosen_index][u][v];
					}

			//set the item of poly as LM[choosen_index] to be 0
			poly[index_y][index_x]=0;

			//cal the coeffTable
			flag=0;	//judge if the poly is empty 
			for(u=0;u<zb_Ysize;u++)
			{
				for(v=0;v<zb_Xsize;v++)
					if( poly[u][v]!=0 )
					{
						flag=1;
						break;
					}

				if( flag==1 )
				{
					break;
				}
			}

			//start to search coefficient
			while( flag!=0 )
			{
				//find the poly's the second largest item
				lod_temp = -1;
				index_x = -1;
				index_y = -1;
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( poly[u][v]!=0 && lod_temp<mono_order[0][u][v] )
						{
							lod_temp=mono_order[0][u][v];
							index_x = v;
							index_y = u;
						}

				//check the correctess of the progress of search second largest item
				if( index_x==-1 || index_y==-1 )
				{
					printf("\n\nsearch progress has error\n\n");
				}
				
				//find the choose_index of zb
				for(z=zb_alpha-1;z>=0;z--)
					if( LM[0][z]==index_x && LM[1][z]==index_y )
					{
						choosen_index = z;
						if( z<tableSize_alpha )
						{
							effTable[i][z][j] = poly[index_y][index_x];
						}
						else
						{
							printf("\n\n1.tableSize_alpha is not large enough\n\n");
						}
						break;
					}

				//update the poly
				temp = poly[index_y][index_x];
				for(u=0;u<zb_Ysize;u++)
					for(v=0;v<zb_Xsize;v++)
						if( zb[choosen_index][u][v]!=0 )
						{
							temp2 =	mul( temp,zb[choosen_index][u][v] );
							poly[u][v] = add( temp2,poly[u][v] );
						}
				
				//check if the poly is empty
				flag=0;	//judge if the poly is empty 
				for(u=0;u<zb_Ysize;u++)
				{
					for(v=0;v<zb_Xsize;v++)
						if( poly[u][v]!=0 )
						{
							flag=1;
							break;
						}

					if( flag==1 )
					{
						break;
					}
				}
				//**************************

			}
			//************************
		}
		//**********************
	}

	for(i=0;i<n;i++)
		for(u=0;u<(max_alpha+1);u++)
			for(v=0;v<(max_a+1);v++)
			{
				coeff_table[i][u][v] = effTable[i][u][v];
			}


	/* debug: printf the coeffTable */
	printf("\ncoeffTable for Herm(8,4)");
	for(i=0;i<n;i++)
	{
		printf("\n**************************");
		printf("\npoint_%d(%d,%d):\n", i, point[i][0], point[i][1]);
		
		for(v=0;v<(max_a+1);v++)
		{
			printf("\t%d",v);
		}
		
		for(u=0;u<(max_alpha+1);u++)
		{
			printf("\n\n");
			for(v=0;v<(max_a+1);v++)
			{
				printf("\t%d", effTable[i][u][v]);
			}
		}
	}
	printf("\n");
	/**********************************/


}	

							
						
	


int mul(int fac1,int fac2)
{
	int mulresult=0;

	if(fac1==0||fac2==0)
		mulresult=0;
	else
	{
		for(int i=0;i<(q-1);i++)
			{
			for(int j=0;j<(q-1);j++)
				{
				if(fac1==mularray[i]&&fac2==mularray[j])
					mulresult=mularray[(i+j)%(q-1)];
				}
			}
	}

	return mulresult;
}

int add(int fac1,int fac2)
{
	int i,c[p],d[p],e[p],f;
	unsigned int mask=1;

	for(i=0;i<p;i++){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}

	mask=1;
	for(i=0;i<p;i++){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}

	for(i=0;i<p;i++)  //p=6
		e[i]=c[i]^d[i];
	
	int num=1;
	f=0;
	for(i=0;i<p;i++)
	{
		f+=(e[i]*num);
		num*=2;
	}

	return f;
}


int power(int a, int b)
{
	int i,temp,pow_result;

	if(b==0)
		pow_result=1;
	else if(a==0 && b!=0)
		pow_result=0;
	else if(a>0 && b!=0)
	{
		for(i=0;i<q-1;i++)
			if(a==mularray[i])
			{
				temp=(i*b)%(q-1);
				pow_result=mularray[temp];
			}
	}

	return pow_result;

}