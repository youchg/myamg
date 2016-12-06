#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "io.h"
#include "preprocess.h"
#include "amg_param.h"
#include "multigrid.h"

#if WITH_BMP
#include "bmp.h"
static void Write_dmatcsr_bmp1(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
static void Write_dmatcsr_bmp2(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
static void Write_dmatcsr_bmp3(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
static void Write_dmatcsr_bmp4(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double));
static void Pivot(double **B, int i, int j, double value);
static void Stretch(double **B, int height, int width);
#endif

dmatcsr *Read_dmatcsr(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
		exit(-1);
    }

    int i = 0;
    int nr, nc, nn;
    while(i<1 && EOF!=fscanf(file, "%d\n",   &nr)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n",   &nc)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n\n", &nn)) i++;
    i = 0;
   
    int    *ia = (int *)malloc((nr+1)*sizeof(int));
    int    *ja = (int *)malloc(nn*sizeof(int));
    double *va = (double *)malloc(nn*sizeof(double));
   
    while(i<nr+1 && EOF!=fscanf(file, "%d\n", ia+i)) i++;
    i = 0;
    while(i<1  && EOF!=fscanf(file, "\n")) i++;
    i = 0;
    while(i<nn && EOF!=fscanf(file, "%d\n", ja+i)) i++;
    i = 0;
    while(i<1  && EOF!=fscanf(file, "\n")) i++;
    i = 0;
    while(i<nn && EOF!=fscanf(file, "%lf\n", va+i)) i++;
    
    dmatcsr *A = (dmatcsr*)malloc(sizeof(dmatcsr));
    A->nr = nr;
    A->nc = nc;
    A->nn = nn;
    
    A->ia = ia;
    A->ja = ja;
    A->va = va;

    fclose(file);
    return A;
}

imatcsr *Read_imatcsr(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(0);
    }

    int i = 0;
    int nr, nc, nn;
    while(i<1 && EOF!=fscanf(file, "%d\n",   &nr)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n",   &nc)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "%d\n\n", &nn)) i++;
    i = 0;
    
   
    int *ia = (int *)malloc((nr+1)*sizeof(int));
    int *ja = (int *)malloc(nn*sizeof(int));
    int *va = (int *)malloc(nn*sizeof(int));
    
    while(i<nr+1 && EOF!=fscanf(file, "%d\n", ia+i)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "\n")) i++;
    i = 0;
    while(i<nn && EOF!=fscanf(file, "%d\n", ja+i)) i++;
    i = 0;
    while(i<1 && EOF!=fscanf(file, "\n")) i++;
    i = 0;
    while(i<nn && EOF!=fscanf(file, "%d\n", va+i)) i++;
    
    imatcsr *A = (imatcsr*)malloc(sizeof(imatcsr));
    A->nr = nr;
    A->nc = nc;
    A->nn = nn;
    
    A->ia = ia;
    A->ja = ja;
    A->va = va;

    fclose(file);
    return A;
}

void Print_dmatcsr(dmatcsr *A)
{
    printf("==================== Brief of the dmatcsr ====================\n");
    printf("nr = %d\n", A->nr);
    printf("nc = %d\n", A->nc);
    printf("nn = %d\n", A->nn);
    printf("============================= end ===========================\n");
}

void Print_imatcsr(imatcsr *A)
{
    printf("==================== Brief of the dmatcsr ====================\n");
    printf("nr = %d\n", A->nr);
    printf("nc = %d\n", A->nc);
    printf("nn = %d\n", A->nn);
    printf("============================= end ===========================\n");
}

void Print_ivec(int *vec, int length)
{
    printf("============================= ivec ===========================\n");
    printf(" index  value\n");
    int i;
    for(i=0; i<length; i++)
    {
	printf("%5d   %3d\n", i, *(vec+i));
    }
    printf("============================= end ============================\n");
}

void Write_dvec(double *x, int length, const char *filename)
{
    FILE *file = fopen(filename,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(-1);
    }
    
    int i;
    
    fprintf(file, "%d\n", length);
    for(i=0; i<length; i++) fprintf(file, "%15.12f\n", x[i]);
    
    fclose(file);
}

void Read_dvec(double *x, int length, const char *filename)
{
    FILE *file = fopen(filename,"r");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(-1);
    }
    
    int i = 0;
    int len = 0;
    
    while(i<1      && EOF!=fscanf(file, "%d\n", &len)) i++; 
    i = 0;
    if(length > 0) assert(length <= len);
    else           x = (double *)malloc(len * sizeof(double));
    while(i<length && EOF!=fscanf(file, "%lf\n", x+i)) i++;
    
    fclose(file);
}

void Write_dmatcsr_csr(dmatcsr *A, const char *filename)
{
    FILE *file = fopen(filename,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(0);
    }
    
    int i;
    
    fprintf(file, "%d\n", A->nr);
    fprintf(file, "%d\n", A->nc);
    fprintf(file, "%d\n", A->nn);
    fprintf(file, "\n");
    
    int nr = A->nr+1; // nr plus 1
    int nn = A->nn;
    
    for(i=0; i<nr; i++) fprintf(file, "%d\n",      A->ia[i]);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%d\n",      A->ja[i]);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%15.12f\n", A->va[i]);
    
    fclose(file);
}

void Write_imatcsr_csr(imatcsr *A, const char *filename)
{
    FILE *file = fopen(filename,"w");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(0);
    }
    
    int i;
    
    fprintf(file, "%d\n", A->nr);
    fprintf(file, "%d\n", A->nc);
    fprintf(file, "%d\n", A->nn);
    fprintf(file, "\n");
    
    int nr = A->nr+1; // nr plus 1
    int nn = A->nn;
    
    for(i=0; i<nr; i++) fprintf(file, "%d\n", A->ia[i]);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%d\n", A->ja[i]);
    fprintf(file, "\n");
    for(i=0; i<nn; i++) fprintf(file, "%d\n", A->va[i]);
    
    fclose(file);
}

void Write_dmatcsr_ps (dmatcsr *A, const char *filename)
{
    FILE *file;
    file = fopen(filename,"w");
    if(!file)
    {
        printf("\ncannot open %s!\n", filename);
	exit(0);
    }

    int size = 500;
    double dot1 = 20.0/(double)A->nr;
    double dot2 = (double)size/1000.0;
    double dot = dot2>dot1 ? dot2 : dot1;
    int lineWidth = 1;
    double step = (double)size/(double)(A->nr-1);
    int i, j;
    
    fprintf(file,"%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(file,"%%%%BoundingBox: 0 0 %d %d\n", size+20, size+20);
    fprintf(file,"%%%%EndComments\n");
    fprintf(file,"10 10 translate\n");
    fprintf(file,"newpath\n");
    fprintf(file,"%f %f moveto\n", 0-2*dot, 0-2*dot);
    fprintf(file,"0 %f rlineto\n", size+4*dot);
    fprintf(file,"%f 0 rlineto\n", size+4*dot );
    fprintf(file,"0 %f rlineto\n", -size-4*dot);
    fprintf(file,"closepath\n");
    fprintf(file,"%d setlinewidth\n",lineWidth);
    fprintf(file,"stroke\n\n");
    fprintf(file,"/box\n");
    fprintf(file,"{gsave\n");
    fprintf(file," translate\n");
    fprintf(file," newpath\n");
    fprintf(file," 0 0 %f 0 360 arc\n", dot);
    fprintf(file," closepath\n");
    fprintf(file," 0 setgray fill\n");
    fprintf(file," grestore} def\n\n");

    for ( i=0;i<A->nr;++i )
    {
	for ( j=A->ia[i];j<A->ia[i+1];++j )
	{
	    fprintf(file,"%f %f box\n", A->ja[j]*step, size-i*step);
	}
    }
    fprintf(file,"showpage\n");
    fclose(file);
}

#if WITH_BMP
void Write_dmatcsr_bmp(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double))
{
    const int nr = A->nr;
    const int nc = A->nc;
    if(width>10000 || height>10000 || width <=0 || height<=0)
    {
        printf("The output bmp file is too big or small! Return.\n");
        return;
    }
    if(nr>=height && nc>=width)
        Write_dmatcsr_bmp1(A, filename, width, height, ColorMap);
    else if(nr>=height && nc<width)
        Write_dmatcsr_bmp2(A, filename, width, height, ColorMap);
    else if(nr<height && nc>=width)
        Write_dmatcsr_bmp3(A, filename, width, height, ColorMap);
    else if(nr<height && nc<width)
        Write_dmatcsr_bmp4(A, filename, width, height, ColorMap);
}

/* 画图时二维数组 B 的最后一行在最上面，即所有行是倒序的，故下面每个函数都有类似下面
Pivot(B, ****height-1-****map_row_A2B[i], map_col_A2B[ja[j]], va[j]);
的措施保证画的矩阵是符合预期的 */
static void Write_dmatcsr_bmp1(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double))
{
    const int nr = A->nr;
    const int nc = A->nc;

    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;
    
    int i, j;
    int rq = nr/height;
    int rr = nr%height;
    int cq = nc/width;
    int cr = nc%width;
    
    double **B = (double**)malloc(height*sizeof(double*));
    for(i=0; i<height; i++)
        B[i] = (double*)calloc(width, sizeof(double));
    
    int *map_row_A2B = (int*)malloc(nr*sizeof(int));
    int *map_col_A2B = (int*)malloc(nc*sizeof(int));
    
    for(i=0; i<rr*(rq+1); i++)
        map_row_A2B[i] = i/(rq+1);
    for(i=rr*(rq+1); i<nr; i++)
        map_row_A2B[i] = (i-rr)/rq;
    for(i=0; i<cr*(cq+1); i++)
        map_col_A2B[i] = i/(cq+1);
    for(i=cr*(cq+1); i<nc; i++)
        map_col_A2B[i] = (i-cr)/cq;
    
    for(i=0; i<nr; i++)
    {
        for(j=ia[i]; j<ia[i+1]; j++)
        {
            Pivot(B, height-1-map_row_A2B[i], map_col_A2B[ja[j]], va[j]);
        }
    }
    Stretch(B, height, width);
    WriteBMPColorMap(filename, width, height, B, ColorMap);
    
    for(i=0; i<height; i++)
       free(B[i]);
    free(B);
    free(map_row_A2B);
    free(map_col_A2B);
}

static void Write_dmatcsr_bmp2(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double))
{
    const int nr = A->nr;
    const int nc = A->nc;

    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;
    
    int i, j, k;
    int rq = nr/height;
    int rr = nr%height;
    int cq = width/nc;
    int cr = width%nc;
    
    double **B = (double**)malloc(height*sizeof(double*));
    for(i=0; i<height; i++)
        B[i] = (double*)calloc(width, sizeof(double));
    
    int *map_row_A2B = (int*)malloc(nr*sizeof(int));
    int *map_col_A2B = (int*)malloc(nc*sizeof(int));
    
    for(i=0; i<rr*(rq+1); i++)
        map_row_A2B[i] = i/(rq+1);
    for(i=rr*(rq+1); i<nr; i++)
        map_row_A2B[i] = (i-rr)/rq;
        
    for(i=0; i<cr; i++)
        map_col_A2B[i] = cq+1;
    for(i=cr; i<nc; i++)
        map_col_A2B[i] = cq;  

    for(i=0; i<nr; i++)
    {
        for(j=ia[i]; j<ia[i+1]; j++)
        {
            if(ja[j] < cr)
            {
                for(k=0; k<cq+1; k++)
                    Pivot(B, height-1-map_row_A2B[i], ja[j]*(cq+1)+k, va[j]);
            }
            else
            {
                for(k=0; k<cq; k++)
                    Pivot(B, height-1-map_row_A2B[i], ja[j]*cq+cr+k, va[j]);
            }
        }
    }
    Stretch(B, height, width);
    WriteBMPColorMap(filename, width, height, B, ColorMap);
    
    for(i=0; i<height; i++)
       free(B[i]);
    free(B);
    free(map_row_A2B);
    free(map_col_A2B);
}

static void Write_dmatcsr_bmp3(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double))
{
    const int nr = A->nr;
    const int nc = A->nc;

    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;
    
    int i, j, k;
    int rq = height/nr;
    int rr = height%nr;
    int cq = nc/width;
    int cr = nc%width;
    
    double **B = (double**)malloc(height*sizeof(double*));
    for(i=0; i<height; i++)
        B[i] = (double*)calloc(width, sizeof(double));
    
    int *map_row_A2B = (int*)malloc(nr*sizeof(int));
    int *map_col_A2B = (int*)malloc(nc*sizeof(int));
    
    for(i=0; i<rr; i++)
        map_row_A2B[i] = rq+1;
    for(i=rr; i<nr; i++)
        map_row_A2B[i] = rq;
        
    for(i=0; i<cr*(cq+1); i++)
        map_col_A2B[i] = i/(cq+1);
    for(i=cr*(cq+1); i<nc; i++)
        map_col_A2B[i] = (i-cr)/cq;
    
    for(i=0; i<nr; i++)
    {
        for(j=ia[i]; j<ia[i+1]; j++)
        {
            if(i < rr)
            {
                for(k=0; k<rq+1; k++)
                    Pivot(B, height-1-(i*(rq+1)+k), map_col_A2B[ja[j]], va[j]);
            }
            else
            {
                for(k=0; k<rq; k++)
                    Pivot(B, height-1-(i*rq+rr+k), map_col_A2B[ja[j]], va[j]);
            }
        }
    }
    Stretch(B, height, width);
    WriteBMPColorMap(filename, width, height, B, ColorMap);
    
    for(i=0; i<height; i++)
       free(B[i]);
    free(B);
    free(map_row_A2B);
    free(map_col_A2B);
}

static void Write_dmatcsr_bmp4(dmatcsr *A, const char *filename, const int width, const int height, RGB_QUAD(*ColorMap)(double))
{
    const int nr = A->nr;
    const int nc = A->nc;

    int    *ia = A->ia;
    int    *ja = A->ja;
    double *va = A->va;
    
    int i, j, k, t;
    int rq = height/nr;
    int rr = height%nr;
    int cq = width/nc;
    int cr = width%nc;
    
    double **B = (double**)malloc(height*sizeof(double*));
    for(i=0; i<height; i++)
        B[i] = (double*)calloc(width, sizeof(double));
    
    int *map_row_A2B = (int*)malloc(nr*sizeof(int));
    int *map_col_A2B = (int*)malloc(nc*sizeof(int));
    
    for(i=0; i<rr; i++)
        map_row_A2B[i] = rq+1;
    for(i=rr; i<nr; i++)
        map_row_A2B[i] = rq;
        
    for(i=0; i<cr; i++)
        map_col_A2B[i] = cq+1;
    for(i=cr; i<nc; i++)
        map_col_A2B[i] = cq;
    
    for(i=0; i<nr; i++)
    {
        for(j=ia[i]; j<ia[i+1]; j++)
        {
            if(i<rr && ja[j]<cr)
            {
                for(k=0; k<rq+1; k++)
                    for(t=0; t<cq+1; t++)
                        Pivot(B, height-1-(i*(rq+1)+k), ja[j]*(cq+1)+t, va[j]);
            }
            else if(i>=rr && ja[j]<cr)
            {
                for(k=0; k<rq; k++)
                    for(t=0; t<cq+1; t++)
                        Pivot(B, height-1-(i*rq+rr+k), ja[j]*(cq+1)+t, va[j]);
            }
            else if(i<rr && ja[j]>=cr)
            {
                for(k=0; k<rq+1; k++)
                    for(t=0; t<cq; t++)
                        Pivot(B, height-1-(i*(rq+1)+k), ja[j]*cq+cr+t, va[j]);
            }
            else
            {
                for(k=0; k<rq; k++)
                    for(t=0; t<cq; t++)
                    {
                        Pivot(B, height-1-(i*rq+rr+k), ja[j]*cq+cr+t, va[j]);
                    }
            }
        }
    }
    Stretch(B, height, width);
    WriteBMPColorMap(filename, width, height, B, ColorMap);
    
    for(i=0; i<height; i++)
       free(B[i]);
    free(B);
    free(map_row_A2B);
    free(map_col_A2B);
}
static void Pivot(double **B, int i, int j, double value)
{
    B[i][j] += value;
}
static void Stretch(double **B, int height, int width)
{
    double max = SMALL;
    double min = LARGE;
    int i, j;
    for(i=0; i<height; i++)
    {
        for(j=0; j<width; j++)
        {
            max = (max<B[i][j])? B[i][j]: max;
            min = (min>B[i][j])? B[i][j]: min;
        }
    }
    assert(max >= min);
    /*For example, [-1, 2] -- [-0.5, 1]; [-2, 1] -- [-1, 0.5]*/
    double absmax = (fabs(max)>fabs(min))? fabs(max): fabs(min);
    if(absmax < EPS)
    {
        return;
    }
    for(i=0; i<height; i++)
    {
        for(j=0; j<width; j++)
        {
            B[i][j] /= absmax;
        }
    }
}
//void Write_imatcsr_bmp(imatcsr *A, const char *filename);

void Input_amg_param(const char *filename, amg_param *param, dmatcsr **A, dmatcsr **M, double **b)
{
    Init_amg_param(param);

    int read   = TRUE;
    int status = TRUE;
    char buffer[512];
    int    ivalue;
    double dvalue;

    FILE *file = fopen(filename,"r");
    if(!file)
    {
        printf("\nError: Cannot open %s!\n", filename);
	exit(-1);
    }

    while(read == TRUE)
    {
	status = fscanf(file, "%s", buffer);
	if(EOF == status)
	{
	    read = FALSE;
	    break;
	}
	if(('#'==buffer[0]) || ('%'==buffer[0]))
	{
	    if(!fgets(buffer, 512, file)) break;
	    continue;
	}

	if(0 == strcmp(buffer, "Afile"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%s", buffer);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    if(NULL != A)
	    {
		*A = Read_dmatcsr(buffer);
	    }
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "Mfile"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%s", buffer);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    if(NULL != M)
	    {
		*M = Read_dmatcsr(buffer);
	    }
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "bfile"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%s", buffer);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    if(NULL != b)
	    {
		Read_dvec(*b, -1, buffer);
	    }
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "strong_connection_threshold"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->strong_connection_threshold = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "strong_diagonally_dominant_threshold"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->strong_diagonally_dominant_threshold = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "truncation_threshold"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->truncation_threshold = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "positive_connected"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->positive_connected = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "interpolation_type"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->interpolation_type = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "max_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->max_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "max_coarsest_dof"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->max_coarsest_dof = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "setup_phase_print_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->setup_phase_print_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "linear_solver_base_type"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->linear_solver_base_type = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "cg_max_iter"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->cg_max_iter = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "cg_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->cg_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "gs_max_iter"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->gs_max_iter = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "gs_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->gs_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_coarsest_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_coarsest_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_mu"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_mu = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_pre_post_smooth"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_pre_post_smooth = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_coarsest_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_coarsest_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_max_coarsest_smooth"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_max_coarsest_smooth = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgcycle_print_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgcycle_print_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgsolver_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgsolver_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgsolver_max_cycle"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgsolver_max_cycle = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgsolver_max_convergence_factor"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgsolver_max_convergence_factor = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgsolver_nmax_convergence_factor"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgsolver_nmax_convergence_factor = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgsolver_print_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgsolver_print_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "pcg_amg_tol"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%lf", &dvalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->pcg_amg_tol = dvalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "pcg_amg_max_iter"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->pcg_amg_max_iter = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "pcg_amg_print_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->pcg_amg_print_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgeigen_coarsest_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgeigen_coarsest_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgeigen_nouter_iter"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgeigen_nouter_iter = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else if(0 == strcmp(buffer, "amgeigen_print_level"))
	{
	    status = fscanf(file, "%s", buffer);
	    if((EOF==status) || (0!=strcmp(buffer, "=")))
	    {
		read = FALSE;
		break;
	    }
	    status = fscanf(file, "%d", &ivalue);
	    if(EOF == status)
	    {
		read = FALSE;
		break;
	    }
	    param->amgeigen_print_level = ivalue;
	    if(!fgets(buffer, 512, file)) break;
	}
	else
	{
	    printf("Error in Input_amg_param: unknown input keyword \"%s\"!\n", buffer);
	    read = FALSE;
	    break;
	}
    }


    fclose(file);
}

void Print_amg_param(amg_param param)
{
    printf("********** amg parameters **********\n");
    printf("strong_connection_threshold          = %g\n", param.strong_connection_threshold);
    printf("strong_diagonally_dominant_threshold = %g\n", param.strong_diagonally_dominant_threshold);
    printf("truncation_threshold                 = %g\n", param.truncation_threshold);
    printf("positive_connected                   = %d\n", param.positive_connected);
    printf("interpolation_type                   = %d\n", param.interpolation_type);
    printf("max_level                            = %d\n", param.max_level);
    printf("max_coarsest_dof                     = %d\n", param.max_coarsest_dof);
    printf("setup_phase_print_level              = %d\n", param.setup_phase_print_level);
    printf("\n");

    printf("linear_solver_base_type = %d\n", param.linear_solver_base_type);
    printf("\n");
    
    printf("cg_max_iter = %d\n", param.cg_max_iter);
    printf("cg_tol      = %g\n", param.cg_tol);
    printf("\n");
    
    printf("gs_max_iter = %d\n", param.gs_max_iter);
    printf("gs_tol      = %g\n", param.gs_tol);
    printf("\n");

    printf("amgcycle_tol                 = %g\n", param.amgcycle_tol);
    printf("amgcycle_coarsest_tol        = %g\n", param.amgcycle_coarsest_tol);
    printf("amgcycle_mu                  = %d\n", param.amgcycle_mu);
    printf("amgcycle_pre_post_smooth     = %d\n", param.amgcycle_pre_post_smooth);
    printf("amgcycle_coarsest_level      = %d\n", param.amgcycle_coarsest_level);
    printf("amgcycle_max_coarsest_smooth = %d\n", param.amgcycle_max_coarsest_smooth);
    printf("amgcycle_print_level         = %d\n", param.amgcycle_print_level);
    printf("\n");

    printf("amgsolver_tol                     = %g\n", param.amgsolver_tol);
    printf("amgsolver_max_cycle               = %d\n", param.amgsolver_max_cycle);
    printf("amgsolver_max_convergence_factor  = %g\n", param.amgsolver_max_convergence_factor);
    printf("amgsolver_nmax_convergence_factor = %d\n", param.amgsolver_nmax_convergence_factor);
    printf("amgsolver_print_level             = %d\n", param.amgsolver_print_level);
    printf("\n");

    printf("pcg_amg_tol         = %g\n", param.pcg_amg_tol);
    printf("pcg_amg_max_iter    = %d\n", param.pcg_amg_max_iter);
    printf("pcg_amg_print_level = %d\n", param.pcg_amg_print_level);
    printf("\n");

    printf("amgeigen_coarsest_level = %d\n", param.amgeigen_coarsest_level);
    printf("amgeigen_nouter_iter    = %d\n", param.amgeigen_nouter_iter);
    printf("amgeigen_print_level    = %d\n", param.amgeigen_print_level);
    printf("************************************\n");
}

void Init_nev_argv(int argc, char* argv[], int *nev, int *nb, int *ne)
{
    int index  = 1;
    while(index < argc)
    {
	if(0 == strcmp(argv[index], "-nev"))
	    *nev = atoi(argv[index+1]);
	if(0 == strcmp(argv[index], "-nb"))
	    *nb = atoi(argv[index+1]);
	if(0 == strcmp(argv[index], "-ne"))
	    *ne = atoi(argv[index+1]);

	index++;
    }
}

void Init_amg_param_argv(int argc, char* argv[], amg_param *param, dmatcsr **A, dmatcsr **M, double **b)
{
    int index  = 1;
    while(index < argc)
    {
	if(0 == strcmp(argv[index], "-ini"))
	{
	    index++;
	    Input_amg_param(argv[index], param, A, M, b);
	    index++;
	}
	else
	{
	    break;
	}
    }
}

void Print_amg(multigrid *amg)
{
    double gc = 0.0;
    double oc = 0.0;
    printf("========================= multigrid =========================\n");
    int i;
    for(i=0; i<amg->actual_level; i++)
    {
        printf("level = %2d, nrow = %7d, nnz = %7d, sparse = %9.6f\n", 
               i, amg->A[i]->nr, amg->A[i]->nn, (double)amg->A[i]->nn/(double)amg->A[i]->nr);
        gc += (double)amg->A[i]->nr;
        oc += (double)amg->A[i]->nn;
    }
    printf("grid complexity = %f, operator complexity = %f\n", 
           gc/(double)amg->A[0]->nr, oc/(double)amg->A[0]->nn);
    printf("=============================================================\n");
}

#endif
