//FiniteFieldBasisCal.cpp
#include <stdio.h>
#include "FiniteFieldBasisGF(16).h"

extern unsigned long int addNum, mulNum;
extern int flag_addNum, flag_mulNum;

extern int mularray[];
extern int root[];
extern int logNum[];

int mul(int fac1,int fac2)
{
	int mulresult=0,temp;

	//mulNum increasing one
	if(flag_mulNum==1)
	{
		mulNum++;
	}

	if(fac1==0||fac2==0)
	{
		mulresult = 0;
		return mulresult;
	}
	else 
	{
		temp = (logNum[fac1]+logNum[fac2]) % (q-1);
		mulresult = mularray[temp];
		return mulresult;
	}

}

int add(int fac1, int fac2)
{
	//addNum increasing one
	if( flag_addNum==1 )
	{
		addNum++;
	}
	return (fac1 ^ fac2);
}


int power(int a, int b)
{
	int temp,pow_result=-1;

	if(b==0)
	{
		pow_result=1;
		return pow_result;
	}
	else if(a==0 && b!=0)
	{
		pow_result=0;
		return pow_result;
	}
	else if(a>0 && b!=0)
	{
		temp = (logNum[a]*b) % (q-1);
		pow_result=mularray[temp];
		return pow_result;
	}

	if( a<0 )
	{	
		printf("\n\n power has error!!");
	}

}