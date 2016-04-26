#ifndef _H_BackInter_
#define _H_BackInter_

void BackInterpolation(int Q_com_poly[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int);
int ItsSize(int *A, int sizeA);
int *ComputeResult(const int Q[][interpoly_Zsize][interpoly_Ysize][interpoly_Xsize], int *A, int , int, int, int);
void Update2(const int alpha, int Q1[][interpoly_Ysize][interpoly_Xsize], const int beta, const int Q2[][interpoly_Ysize][interpoly_Xsize]);
void UpdatePolyWithDivision(int Q[][interpoly_Ysize][interpoly_Xsize], int alpha);


#endif