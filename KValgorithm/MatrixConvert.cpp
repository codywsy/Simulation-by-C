//**********************************************
//作用：将reliability matrix 转化成 Multiplicity matrix
//输入：
//		1. reliability Matrix RM[q][n]
//		2. s(multiplicty sum)
//输出：Multiplicity matrix MM[q][n]
//***********************************************
void MatrixConvert()
{
	int i, j, index_i, index_j;
	float RM_temp[q][n], tempRM, max_temp;
	int MM[q][n], tempMM;

	//initialization
	for(i=0;i<n;i++)
		for(j=0;j<q;j++)
		{
			RM_temp[j][i] = RMj][i];
			MM[j][i] = 0; 
		}
		
	//start
	while(s>0)
	{
		max_temp = 0.0;
		index_j = -1;
		index_i = -1;

		//find the maximal entry factor in RM
		for(i=0;i<n;i++)
			for(j=0;j<q;j++)
			{
				if( max_temp < RM_temp[j][i] )
				{
					max_temp = RM_temp[j][i];
					index_j = j;
					index_i = i;

					//because the sum of a cloumn of RM is equal to 1. so 
					//for a cloumn, if there is a factor larger than 0.5, that must
					//be the largest one of the cloumn
					if( max_temp > 0.5 )
					{
						max_temp = RM_temp[j][i];
						index_j = j;
						index_i = i;
						break; //jump out the cloumn circle
					}
				}
			}

		//update
		RM_temp[index_j][index_i] = RM_temp[index_j][index_i] / (MM[index_j][index_i]+2);
		MM[index_j][index_i] += 1;
		s -= 1; 
	}
	


}