#ifndef _H_BFLCC_Inter_
#define _H_BFLCC_Inter_

#define init_polyNum (lm+1)	//change with diff lm, the poly num of the init polyGroup
#define interpoly_Ysize (lm+1)	//maxdeg of y is (w-1) + w*(n/(w+1)), so the Ysize should add 1 more.
#define interpoly_Xsize	(iterNum+1)	//maxdeg of x is w, so the Xsize should add 1 more.

void BF_LCC_inter();
void InterOnce_2(int poly[][interpoly_Ysize][interpoly_Xsize], int interpoint[3], int poly_min_order[eta], int location, int poly_min_order_all[n], int location_all);
int ChooseTheMinPoly(int polySour[][interpoly_Ysize][interpoly_Xsize], int PolyNum);


#endif