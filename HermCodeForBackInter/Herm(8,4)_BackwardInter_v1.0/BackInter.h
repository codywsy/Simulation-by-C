#ifndef _H_BackInter_
#define _H_BackInter_

#define New_interpoly_Zsize (lm+1)
#define New_interpoly_Ysize w
#define New_interpoly_Xsize (iterNum+1)


void BackInterpolation(int Q_com_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int BackPoint[3]);
void BackInterpolation_2(int Q_com_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int backPoint[3], int poly_min_order[eta], int location, int lod[init_polyNum], int poly_min_order_all[n], int *location_all);
int ItsSize(int *A, int sizeA);
int *ComputeResult(const int Q[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int *A, int, int, int, int);
void BF_Update(const int alpha, int Q1[][New_interpoly_Ysize][New_interpoly_Xsize], const int beta, const int Q2[][New_interpoly_Ysize][New_interpoly_Xsize]);
void UpdatePolyWithDivision(int Q[][New_interpoly_Ysize][New_interpoly_Xsize], int alpha);
void DectectIfFactorInPoly(int Q_poly[][New_interpoly_Zsize][New_interpoly_Ysize][New_interpoly_Xsize], int polyNum, int xi);

#endif