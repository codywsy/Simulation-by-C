//FiniteFieldBasisGF(16).cpp

//int mularray[] = {1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};  //this array is used to the infinite field element mul
//int root[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};	//n+1		be used to the factorization
//int logNum[] = {-1,0,1,4,2,8,5,10,3,14,9,7,6,13,11,12};	//used to locate the degree of finitefield element through the value of finitefield element

int mularray[] = { 1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20, 40,
19, 38, 15, 30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56, 51, 37,
9, 18, 36, 11, 22, 44, 27, 54, 47, 29, 58, 55, 45, 25, 50, 39,
13, 26, 52, 43, 21, 42, 23, 46, 31, 62, 63, 61, 57, 49, 33 };  //this array is used to the infinite field element mul

int root[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };	//n+1		be used to the factorization

int logNum[] = { -1, 0, 1, 6, 2, 12, 7, 26, 3, 32, 13, 35, 8, 48, 27, 18,
4, 24, 33, 16, 14, 52, 36, 54, 9, 45, 49, 38, 28, 41, 19,
56, 5, 62, 25, 11, 34, 31, 17, 47, 15, 23, 53, 51, 37, 44,
55, 40, 10, 61, 46, 30, 50, 22, 39, 43, 29, 60, 42, 21, 20,
59, 57, 58 };	//used to locate the degree of finitefield element through the value of finitefield element