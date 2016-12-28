#ifndef __TOOL_H__
#define __TOOL_H__

double Get_time();

void Insertion_ascend_sort_dvec(double *a, int left, int right);
void Quick_ascend_sort_ivec    (int    *a, int left, int right);

void Insertion_ascend_sort_ivec_dvec(int *a, double *b, int left, int right);
void Quick_ascend_sort_ivec_dvec    (int *a, double *b, int left, int right);

void Insertion_ascend_sort_dvec_dvecvec(double *a, double **b, int left, int right);

int Is_ascend_sorted_ivec(int *vec, int length);

int A_cap_B_sorted_ascend(int *A, int lengthA, int *B, int lengthB);
int A_cap_B(int *A, int lengthA, int *B, int lengthB);
int Get_ivec_cap_ivec(int *A, int lengthA, int *B, int lengthB, int *ncap, int **cap, int **indexA, int **indexB);

#endif
