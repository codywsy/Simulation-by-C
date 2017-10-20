#ifndef _H_main_
#define _H_main_

//main.h
#define w 2
#define weiz 4
#define m_value_max 2
//#define iterNum 8	//when m=1, C is equal to n
#define lm 2
#define able_correct 1
#define pointNum 2
#define interval 1

//test_vec_construction
#define choose_num 2	//fixed
#define PROBA_W 0.5

//*****need to be modify_first**************
//coefficientSearch()
#define max_m 4	//design large enough to cover, dependent on the value of lm 
#define max_alpha 8	//alpha < max_m, designing alpha large enough to cover all the m_ij of Mulplicity Matrix
#define tableSize_a 360	//more than max_a, can be achieved by the hint of coefficientSearch()
#define tableSize_alpha 370	//more than max_alpha, can be achieved by the hint of coefficientSearch()
//zerobasie()
#define zb_alpha 400	//size+1, larger than tableSize_alpha, according to the maxDeg_x and maxDeg_y of PbMax, we can judge the max num Zb we need
#define zb_Ysize (zb_alpha/(w+1)+1)	//this term must be changed when alpha greater than (w+1)
#define zb_Xsize (w+1)
//cal_max_a()
#define N 800	//larger than tg_size
//*****************************************

//*****need to be modify_second**************
//tgorder()
#define tg_size 186	//more than k, tg_size represent the probably used pole basis num, (w+1)*w/2 + (w+1) * ( (w+interpoly_Ysize)-w+1 ) 
//mono_table()
#define weight_Zsize 8
#define weight_XYsize tg_size	//equals to the tg_size
#define mono_ordinSize 29	//choose the value with debugging, equals to max weigtdegree in weight[monoTable_Zsize][weight_XYsize]
#define monoTable_Zsize	(lm+1) 	//equals to the interpoly_Zsize
#define monoTable_Ysize 10	//large than interpoly_Ysize
#define monoTable_Xsize (w+1)	//large than the interpoly_Xsize
#define monoTable_totalSize tg_size //>index of term[max(deg_y)][w] in the pole basis + 1, nearly equals to tg_size
//*****************************************

//*****need to be modify_third**************
//interpolation()
#define init_polyNum (w*(lm+1))	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Zsize (lm+1)	//maxValue of z is lm=1, so the Zsize should add 1 more.
#define interpoly_Ysize 7	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(w+1)	//maxdeg of x is w, so the Xsize should add 1 more.
//factorization()
#define facpoly_Zsize interpoly_Zsize// same as the interpoly_Zsize
#define facpoly_Ysize interpoly_Ysize		// same as the interpoly_Ysize
#define facpoly_Xsize interpoly_Xsize	// same as the interpoly_Xsize
//rcs()
#define rcspoly_Ysize 14	//degY+1, degY = interpoly_Ysize + max(j_1) + w*( (max(i_1)+w)/(w+1) )
#define rcspoly_Xsize 12	//degX+1, degX = max(i_1) + w
#define faiMax_Ysize 2	//change with diff w and k, the max degY of probably used polebasis
#define faiMax_Xsize 3	//change with diff w and k, the max degX of probably used polebasis
//expoly()-->expanded polynomial
#define expoly_Ysize (lm*faiMax_Ysize+1)	
#define expoly_Xsize (lm*faiMax_Xsize+1)
//*******************************************




//conditional compile
//#define _GS_Normal_
#define _NoReduction_

#ifndef _GS_Normal_
#define eta 2
#define test_vec_num 4	//the num of test_vec is 2^eta
#else
#define eta 0
#define test_vec_num 1 //the num of test_vec is 2^eta
#endif

#define OpenFile fp=fopen("KV-LCC_Herm(8,4)_¦Ç=2.txt","a")
#define FrameError 309

//#define printDemodulation

#endif