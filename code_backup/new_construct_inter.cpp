#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 7	//length of codeword
#define k 2 //length of message
#define q 8 //GF(q)
#define p 3 //GF(2^p) = GF(q)
#define maxcost 300 //Designed maximal interpolation cost

//**************
#define m 2 // zero of multiplicity
#define iter_num 21     // iter_num = C
#define tm 3 
#define lm 5   //Desgined factorization output list size   (should be mended)
#define mono_size 27  // the size of inter_polynomial
#define fa_size 30
#define deg_size 20
//******************


int mul(int, int);
int add(int, int);
int power(int, int);
int inv(int);
void mono_table(void);
void encoder(void);
void modulation(void);
void channel(void);
void demodulation(void);
//******
void decoding(void);
int cal_poly_lod(int);	
int cal_minlod(void);	
int cal_delta(int,int,int);
void cal_Newinter_poly(int,int);
void cal_Newpoly_minlod(int,int);
void cal_Newpoly(int,int);
void refresh_poly_temp(int array[][mono_size+1]);
int fac(int);
//***********


unsigned long int seq_num;	//number of input binary sequences
int SNR;
double BER, FER;
int mularray[]={1,2,4,3,6,7,5};
int root[]={0,1,2,3,4,5,6,7};	//n+1	
int a[7]={1,2,4,3,6,7,5};

int bi_message[k*p];	//binary input message
int message[k];			//decimal input message
int codeword[n];		//encoded code word
int bi_codeword[n*p];   //encoded binary codeword
int dec_codeword[n*p];  //decoded binary codeword
float tx_symbol[n*p][2];	//modulated tranmitted signal
float rx_symbol[n*p][2];	//received destorted signal

int mono_order[40][200], channel_type, kv_signal;	//multiplexity matrix
float N0, sgm;

//*********
int rx_bicodeword[n*p];   //demodulated binary codeword
int rx_codeword[n];   //demodulated codeword to be decoded
int inter_point[n][2];  //interlolation points (pi,ri)
int inter_poly[(lm+1)][mono_size][mono_size];  //store the inter_poly of everygroup
int inter_poly_minlod[mono_size][mono_size];  // the poly with minimium used to factorization
int poly_lod[(lm+1)];			//store the lod of inter_poly of everygroup
int poly_lod_temptfg[(lm+1)];		// lod temp
int alpha_beta[iter_num/n][2];     //store alpha and beta
int delta[lm+1];					   //store the group of delta, totalling

int message_poly[fa_size][mono_size];   //the array of poly used in factorization
int time;     // factorization time;

//factorization's element
int last[fa_size],deg[fa_size];
int Q[fa_size][mono_size][mono_size];
int Coeff[mono_size];
int f[n][k];		//store the message polynomial after the factorization, note that may be overflow
int f_i=0;

//*********


void main()
{
	int i;

	mono_table();

	//initialize the rx_codeword
	rx_codeword[0]=1;
	rx_codeword[1]=4;
	rx_codeword[2]=2;
	rx_codeword[3]=1;
	rx_codeword[4]=5;
    rx_codeword[5]=2;
	rx_codeword[6]=2;

	printf("re_codeword=(");
	for(i=0;i<n-1;i++)
		printf("%d,",rx_codeword[i]);
	printf("%d)\n",rx_codeword[6]);

	decoding();


}

void mono_table(void)
{
	int i, j, v, l;

	for(j=0;j<40;j++)
		for(i=0;i<200;i++)
			mono_order[j][i]=-1;

	j=0;	//represent row
	v=0;	//increasing counter for monomial order
	for(i=0;i<200;i++)
	{
		if(i<(k-1))	
		{
			mono_order[j][i]=v;
			v++;
		}
		else
		{
			l=floor((float)i/(float)(k-1));	    // l=[i/(k-1)]
			for(j=0;j<=l;j++)
			{
				mono_order[j][i-(k-1)*j]=v;	//i-v*j
				v++;
			}
		}
	}
}





int power(int a, int b)
{
	int i, pow_result;

	pow_result=1;
	for(i=0;i<b;i++)
		pow_result=mul(pow_result,a);

	return pow_result;
}

int mul(int fac1,int fac2)
{
	int mulresult;

	if(fac1==0||fac2==0)
		mulresult=0;
	else
		for(int i=0;i<7;i++){
			for(int j=0;j<7;j++){
				if(fac1==mularray[i]&&fac2==mularray[j])
					mulresult=mularray[(i+j)%7];
			}
		}

	return mulresult;
}

int add(int fac1,int fac2)
{
	int i,c[3],d[3],e[3],f;
	unsigned int mask=1;

	for(i=2;i>=0;i--){
		if((fac1 & mask)>0)
			c[i]=1;
		else
			c[i]=0;

		mask=mask<<1; //shift 1 bit left
	}
	mask=1;
	for(i=2;i>=0;i--){
		if((fac2 & mask)>0)
			d[i]=1;
		else
			d[i]=0;

		mask=mask<<1;
	}
	for(i=0;i<3;i++)
		e[i]=c[i]^d[i];


	f=e[0]*4+e[1]*2+e[2]*1; //conver to decimal

	return f;
}



//*********
void decoding(void)
{
	int i,j,h,temp;
	int tk,t;   			// tk is the num of iteration, t is the number of inter_point
	int minlod_seq,minlod=0;	//minlod_seq is the j' with delta[j']!=0


	for(i=0;i<n;i++)					// interpolation points
	{
		inter_point[i][0]=a[i]; 		// the first word is pi--->xi
		inter_point[i][1]=rx_codeword[i];	// the second word is ri--->yi
	}

	//******检测项*******
		for( h=0;h<n;h++)
			printf("(%d,%d)\t",inter_point[h][0],inter_point[h][1]);
		printf("\n\n");
	//************


    //initialize group_inter_poly
    for(j=0;j<(lm+1);j++)
   	{
   		for(h=0;h<mono_size;h++)
   		{
   			for(i=0;i<mono_size;i++)
   			{
   				inter_poly[j][h][i]=0;
   			}
   		}
   	}

    
     //initialize alpha_beta points
    for(i=0,j=0;j<m;j++)     
    {
    	for(h=0;h<(m-j);h++)
    	{
    		if(i<(iter_num/n))
  			{
  				alpha_beta[i][0]=j;
   				alpha_beta[i][1]=h;
   				i++;
   			}
   		}
   	}

	//******检测项*******
		for( h=0;h<(iter_num/n);h++)
			printf("(a,b)_%d=(%d,%d)\t",h,alpha_beta[h][0],alpha_beta[h][1]);
		printf("\n\n");
	//************

    //initialize Group_0
    for(i=0;i<(lm+1);i++)
    	inter_poly[i][i][0]=1;

	//******检测项*******    
	for(i=0;i<(lm+1);i++)      
	{
		poly_lod[i]=cal_poly_lod(i); 
		printf("poly_lod[%d]=%d\t",i,poly_lod[i]);
	}
	printf("\n************begin********************\n\n");
	//****************

    //iteration of modifying the inter_polynomial   -------> optimize to a function
    for(tk=0;tk<n;tk++)
    {

    	for(t=0;t<(iter_num/n);t++)
    	{
    		//calculate delta_i
    		for(i=0;i<(lm+1);i++)
    		{
    			delta[i]=cal_delta(i,tk,t);
    		}

			//*******检查项*********
			for(i=0;i<(lm+1);i++)
				printf("delta_%d[%d]=%d\t",(tk+t),i,delta[i]);
			printf("\n");
			//*********************

       		//choose the delta_i with delta != 0 and j' = minlod_seq
  			minlod_seq=cal_minlod();

			//*******检查项*********
			printf("minlod_seq(%d)=%d\n\n",((tk+1)*t),minlod_seq);
			//****************

		 	//calculate the new inter polynomial -------> update
		 	cal_Newinter_poly(minlod_seq,tk);



		 	// work out the new poly_lod of inter_poly[i]
			for(i=0;i<(lm+1);i++)      
		 	poly_lod[i]=cal_poly_lod(i);  

			//******检测项*******
			for(i=0;i<(lm+1);i++)      
			{
				poly_lod[i]=cal_poly_lod(i); 
				printf("poly_lod[%d]=%d\t",i,poly_lod[i]);
			}
			printf("\n\n");
			//****************

 		 	//find the interpolation polynomial of Group_C with minimium lod-----------> inter_poly[minlod]
 		 	temp=poly_lod[0];
		 	for(i=0;i<(lm+1);i++)
		 		if( temp > poly_lod[i] )
		  		{
		  			temp=poly_lod[i];
		  			minlod=i;
		  		}
		}
    }

	// find out the inter_poly_minlod used to factorization
    for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)
			inter_poly_minlod[j][h]=inter_poly[minlod][j][h];

}



int cal_poly_lod(int i)   //input one inter_polynomial, output this poly's leading order
{
	int j,h;
	int max_ord=0;

	for(j=0;j<mono_size;j++)
		{
			for(h=0;h<mono_size;h++)
			{
				if( inter_poly[i][j][h] != 0 )
				{
					if(	max_ord < mono_order[j][h])
					{
						max_ord=mono_order[j][h];
					}
				}
			}
		}
	return max_ord;
}


int cal_minlod(void)				//calculate the minlod of this group
{
	int i,temp,flag=1,minlod_seqnum=0;

	//calculate the lod of poly and store it in poly_lod[i]
	for(i=0;i<(lm+1);i++)
  		poly_lod[i]=cal_poly_lod(i);

    //find out the minlod_seqnum of polynomial Group
  	for(i=0;i<(lm+1);i++)
  		if( delta[i]!=0 )
  		{
  			if(flag==1)
  			{
  				temp = poly_lod[i];
				minlod_seqnum = i;
  				flag=0;
  			}

  			if( temp > poly_lod[i] )
  			{
  				temp = poly_lod[i];
  				minlod_seqnum = i;
  			}
  		}

  	return minlod_seqnum;

}


int cal_delta(int i,int tk,int t)    //calculate the delta of inter_polynomial, c is inter_poly[i] and t is tk
{
	int delta=0,delta_temp,h,j,flag;
	int a,b,alpha,beta,x,y;

	alpha = alpha_beta[t][0];    //alpha
	beta = alpha_beta[t][1];		//beta
	x=inter_point[tk][0];			//xi=pi
	y=inter_point[tk][1];			//yi=ri

	/*
	//********检查项**********
	printf("alpha=%d,beta=%d,x_%d=%d,y_%d=%d\t",alpha,beta,tk,x,tk,y);
	printf("\n");
	//*********************** 
	*/

	for(j=0;j<mono_size; j++)
	{
		for(h=0;h<mono_size;h++)
		{
			if( (inter_poly[i][j][h] != 0) && ( h >= alpha) && ( j >= beta) )
			{
				a=h;
				b=j;
				flag= ( fac(a) / ( fac(a-alpha)*fac(alpha) ) ) * ( fac(b) / ( fac(b-beta)*fac(beta) ) );
				// if flag is even number , delta_temp is equal to 0
				// if flag is odd number, delta_temp can be calculated to the delta
				/*
				//***********
				 if( inter_poly[i][j][h] !=0 )
					 printf("flag[%d][%d]=%d\t",j,h,flag);
				//************
				*/
				if( (flag%2) !=0 )
				{
					delta_temp=mul(inter_poly[i][j][h],power(x,(a-alpha)));
					delta_temp=mul(delta_temp,power(y,(b-beta)));
					//*******
					printf("temp[%d][%d]=%d\t",a,b,delta_temp);
					//**********
					delta= add(delta,delta_temp);
				}
			}
		}	
	}
	printf("\n");


	return delta;
}



void cal_Newinter_poly(int minlod_seq,int tk)
{
	int i;

	// start to update
	// (delta != 0)-----> update  
    for( i=0 ; i<(lm+1) ; i++ )
        if( (delta[i]!=0) && (i!=minlod_seq) )
        	cal_Newpoly(i,minlod_seq);	

    // at last update the inter_poly[minlod_seq]
	cal_Newpoly_minlod(minlod_seq,tk);

	/*
	//************check*********
		printf("x%d=%d\n",tk,inter_point[tk][0]);
		printf("inter_poly[0]:\n");
		for(int j=0;j<mono_size;j++)
			{
				for(int h=0;h<mono_size;h++) 
					printf("%d\t",inter_poly[0][j][h]);
				printf("\n\n");
			}
		printf("inter_poly[1]:\n");
		for(int j=0;j<mono_size;j++)
			{
				for(int h=0;h<mono_size;h++) 
					printf("%d\t",inter_poly[1][j][h]);
				printf("\n\n");
			}
	//********************************
    */
}



void cal_Newpoly_minlod(int minlod_seq,int tk)					
{
	int j,h,temp;
	int poly_temp1[mono_size+1][mono_size+1];
	int poly_temp2[mono_size+1][mono_size+1];

	refresh_poly_temp(poly_temp1);
	refresh_poly_temp(poly_temp2);

	
	// delta[i]*x*inter_poly[i]
	for( j=0 ; j<mono_size ; j++ )
   		for( h=0 ; h<mono_size ; h++ )
       		if( inter_poly[minlod_seq][j][h] != 0)
       			poly_temp1[j][h+1]=mul(inter_poly[minlod_seq][j][h],delta[minlod_seq]);

    // delta[i] * xi * inter_poly[i]
  	temp=mul(delta[minlod_seq],inter_point[tk][0]);    //delta[i] * xi
  
	for( j=0 ; j<mono_size ; j++ )
   		for( h=0 ; h<mono_size ; h++ )
       		if( inter_poly[minlod_seq][j][h] != 0)
       			poly_temp2[j][h] = mul(inter_poly[minlod_seq][j][h],temp);

    // inter_poly + inter_poly_temp
    for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)
			if( (poly_temp1[j][h]==0) && (poly_temp2[j][h]==0) )  
				inter_poly[minlod_seq][j][h]=0;
			else
				inter_poly[minlod_seq][j][h]=add(poly_temp1[j][h],poly_temp2[j][h]);
}

void cal_Newpoly(int i,int minlod_seq)
{
	int h,j;
	int poly_temp1[mono_size+1][mono_size+1];
	int poly_temp2[mono_size+1][mono_size+1];

	refresh_poly_temp(poly_temp1);
	refresh_poly_temp(poly_temp2);

	// delta[j'] * inter_poly[i]
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  
			if( inter_poly[i][j][h] != 0 )
				poly_temp1[j][h]=mul(inter_poly[i][j][h],delta[minlod_seq]);
	// delta[i] * Q'(inter_poly_temp)
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  	
			if( inter_poly[minlod_seq][j][h] != 0 )
				poly_temp2[j][h]=mul(inter_poly[minlod_seq][j][h],delta[i]);
	//  new inter_poly[i+1] = new inter_poly[i] + new Q'
	for( j=0;j<mono_size;j++)
		for( h=0;h<mono_size;h++)  
			inter_poly[i][j][h]=add(poly_temp1[j][h],poly_temp2[j][h]);	


}


void refresh_poly_temp(int array[][mono_size+1])     //clean the poly_temp
{
	int h,j;

	for(j=0;j<(mono_size+1);j++)
   	{
   		for(h=0;h<(mono_size+1);h++)
   		{
   			array[j][h]=0;
   		}
   	}
}

int fac(int a)		//calculate the factorial of a
{
	int i,result=1;

	for(i=1;i<=a;i++)
		result=result*i;

	return result;
}
