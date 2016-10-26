#ifndef _H_main_
#define _H_main_

//main.h
#define iterNum 7	//when m=1, C is equal to n
#define lm 1
#define able_correct ((n-k+1)/2)
#define pointNum 2
#define interval 1

//test_vec_construction
#define choose_num 2	//fixed

//mono_table()
#define monoTable_Ysize (lm+1+1)	//large than interpoly_Ysize
#define monoTable_Xsize (iterNum+1)	//large than the interpoly_Xsize
#define target_weight_degree ((lm+1)*(k-1)+(iterNum+1)-1)
#define weight_Ysize (target_weight_degree+1)	//larger than monoTable_Ysize, decided by the size of mono_Table
#define weight_Xsize (target_weight_degree+1)	//larger than monoTable_Ysize, decided by the size of mono_Table


//interpolation()
#define init_polyNum (lm+1)	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Ysize (lm+1)	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(iterNum+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Ysize interpoly_Ysize		// same as the interpoly_Ysize
#define facpoly_Xsize interpoly_Xsize	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 12	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 10	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 2	//change with diff w and k, the max degY of probably used polebasis
#define faiMax_Xsize 3	//change with diff w and k, the max degX of probably used polebasis
//expoly()-->expanded polynomial
#define expoly_Ysize (lm*faiMax_Ysize+1)	
#define expoly_Xsize (lm*faiMax_Xsize+1)


//conditional compile
//#define _GS_Normal_

#ifndef _GS_Normal_
#define eta 3
#define test_vec_num 8	//the num of test_vec is 2^eta
#else
#define eta 0
#define test_vec_num 1 //the num of test_vec is 2^eta
#endif


#define _Check_Com_Inter_
#define _Check_Uncom_Inter_
//#define _PrintGroupLod_

//#define _PrintLCCLod_

#endif