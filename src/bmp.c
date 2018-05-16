#include <stdlib.h>
#include <math.h>
#include "bmp.h"

#define PRINTLEVEL 0

BMP_IMAGE* ReadBMP(const char *bmpfile)
{
#if PRINTLEVEL > 3
    printf("begin read bmp...\n");
#endif
    FILE *fp = NULL;
    //fopen_s(&fp, bmpfile, "rb");
    fp = fopen(bmpfile, "rb");
    if (fp == NULL)
    {
        printf("fp == NULL\n");
        return NULL;
    }
    
#if PRINTLEVEL
    //===================== test data type size =================
    //size_t 在32位系统中是 unsigned int,
    //       在64位系统中是 long unsigned int
    printf("sizeof(BIT_MAP_FILE_HEADER) = %lu\n", sizeof(BIT_MAP_FILE_HEADER));
    printf("sizeof(BIT_MAP_INFO_HEADER) = %lu\n", sizeof(BIT_MAP_INFO_HEADER));
    printf("sizeof(long)                = %lu\n", sizeof(long));
    printf("sizeof(int)                 = %lu\n", sizeof(int));
    printf("sizeof(unsigned long)       = %lu\n", sizeof(unsigned long));
    printf("sizeof(unsigned int)        = %lu\n", sizeof(unsigned int));
    printf("sizeof(unsigned short)      = %lu\n", sizeof(unsigned short));
    //===========================================================
#endif
    int read_count = 0;
    BIT_MAP_FILE_HEADER bmpFileHeader;
    BIT_MAP_INFO_HEADER bmpInfoHeader;
    read_count += fread(&bmpFileHeader, sizeof(BIT_MAP_FILE_HEADER), 1, fp);
    read_count += fread(&bmpInfoHeader, sizeof(BIT_MAP_INFO_HEADER), 1, fp);

    if (bmpFileHeader.bfType != 0x4D42)
    {
        printf("error: bfType != 0x4D42\n");
        return NULL;
    }
#if PRINTLEVEL
    printf("=============== bmpFileHeader ===============\n");
    printf("bfType      = %X\n", bmpFileHeader.bfType);
    printf("bfSize      = %d\n", bmpFileHeader.bfSize);
    printf("bfReserved1 = %d\n", bmpFileHeader.bfReserved1);
    printf("bfReserved2 = %d\n", bmpFileHeader.bfReserved2);
    printf("bfOffBits   = %d\n", bmpFileHeader.bfOffBits);
    printf("=============================================\n\n");
    
    printf("=============== bmpInfoHeader ===============\n");
    printf("biSize          = %d\n", bmpInfoHeader.biSize);
    printf("biWidth         = %d\n", bmpInfoHeader.biWidth);
    printf("biHeight        = %d\n", bmpInfoHeader.biHeight);
    printf("biPlanes        = %d\n", bmpInfoHeader.biPlanes);
    printf("biBitCount      = %d\n", bmpInfoHeader.biBitCount);
    printf("biCompression   = %d\n", bmpInfoHeader.biCompression);
    printf("biSizeImage     = %d\n", bmpInfoHeader.biSizeImage);
    printf("biXPelsPerMeter = %d\n", bmpInfoHeader.biXPelsPerMeter);
    printf("biYPelsPerMeter = %d\n", bmpInfoHeader.biYPelsPerMeter);
    printf("biClrUsed       = %d\n", bmpInfoHeader.biClrUsed);
    printf("biClrImportant  = %d\n", bmpInfoHeader.biClrImportant);
    printf("=============================================\n\n");
#endif
	
    BMP_IMAGE *bmp = (BMP_IMAGE*)malloc(sizeof(BMP_IMAGE));
    bmp->rgb       = NULL;
    bmp->bmpData   = NULL;

    int width  = bmpInfoHeader.biWidth;
    int height = bmpInfoHeader.biHeight;
    int bit    = bmpInfoHeader.biBitCount;

    bmp->width  = width;
    bmp->height = height;
    bmp->bit    = bit;
	
    if (bit == 8)
    {
#if PRINTLEVEL > 3
        printf("grayscale image\n");
#endif
        //灰度图像有256项颜色表
        bmp->rgb = (RGB_QUAD*)malloc(sizeof(RGB_QUAD) * 256);
        read_count += fread(bmp->rgb, sizeof(RGB_QUAD), 256, fp);
    }
    else if (bit == 24)
    {
#if PRINTLEVEL > 3
        printf("truecolor image\n");
#endif
        //真彩色图像没有颜色表
        bmp->rgb = NULL;
    }
    else
    {
        printf("error: biBitCount is neither 8 nor 24!\n");
        free(bmp);
        return NULL;
    }
	
    //计算每行像素的字节数
    int lineByte = width * bit / 8;
    //由于bmp每行像素的实际字节数必须是4的整数倍
    //所以计算补齐需要添加的字节数
    int align = lineByte % 4;
    if (align != 0)
        align = 4 - align;
#if PRINTLEVEL > 3
    printf("lineByte = %d\n", lineByte);
    printf("align = %d\n", align);
#endif
    bmp->bmpData = (unsigned char*)malloc(sizeof(unsigned char)*lineByte*height);
    unsigned char *bmpData = bmp->bmpData;

#if PRINTLEVEL > 3
    printf("bmp->bmpData = %lu\n", sizeof(unsigned char)*lineByte*height);
#endif
    int i;
    for (i = 0; i<height; i++)
    {
        read_count += fread(bmpData, lineByte, 1, fp);
        bmpData += lineByte;
        fseek(fp, align, SEEK_CUR);
    }
    fclose(fp);
    
#if PRINTLEVEL > 3
    printf("read count = %d\n", read_count);
#endif
    return bmp;
}

int WriteBMP(const char *bmpfile, BMP_IMAGE *bmp)
{
#if PRINTLEVEL > 3
    printf("begin write bmp...\n");
#endif
    FILE *fp = NULL;
    //fopen_s(&fp, bmpfile, "wb");
    fp = fopen(bmpfile, "wb");
    if (fp == NULL)
    {
        printf("fp == NULL\n");
        return 1;
    }

    int width              = bmp->width;
    int height             = bmp->height;
    int bit                = bmp->bit;
    RGB_QUAD *rgb          = bmp->rgb;
    unsigned char *bmpData = bmp->bmpData;

    int lineByte = width * bit / 8;
    int align    = lineByte % 4;
    if (align != 0)
        align = 4 - align;
#if PRINTLEVEL > 3
    printf("lineByte = %d\n", lineByte);
    printf("align = %d\n", align);
#endif

    BIT_MAP_FILE_HEADER bmpFileHeader;
    bmpFileHeader.bfType      = 0x4D42;
    bmpFileHeader.bfReserved1 = 0;
    bmpFileHeader.bfReserved2 = 0;

    if (bit == 8)
    {
#if PRINTLEVEL > 3
        printf("grayscale image\n");
#endif
        bmpFileHeader.bfSize    = 54 + 256 * 4 + (lineByte + align)*height;
        bmpFileHeader.bfOffBits = 54 + 256 * 4;
    }
    else if (bit == 24)
    {
#if PRINTLEVEL > 3
        printf("truecolor image\n");
#endif
        bmpFileHeader.bfSize    = 54 + (lineByte + align)*height;
        bmpFileHeader.bfOffBits = 54;
    }
    else
    {
        printf("error: biBitCount is neither 8 nor 24!\n");
        return 2;
    }

    BIT_MAP_INFO_HEADER bmpInfoHeader;
    bmpInfoHeader.biSize          = 40;
    bmpInfoHeader.biWidth         = bmp->width;
    bmpInfoHeader.biHeight        = bmp->height;
    bmpInfoHeader.biPlanes        = 1;
    bmpInfoHeader.biBitCount      = bit;
    bmpInfoHeader.biCompression   = 0;
    bmpInfoHeader.biSizeImage     = (lineByte + align)*height;
    bmpInfoHeader.biXPelsPerMeter = 0;
    bmpInfoHeader.biYPelsPerMeter = 0;
    bmpInfoHeader.biClrUsed       = 0;
    bmpInfoHeader.biClrImportant  = 0;

    fwrite(&bmpFileHeader, sizeof(BIT_MAP_FILE_HEADER), 1, fp);
    fwrite(&bmpInfoHeader, sizeof(BIT_MAP_INFO_HEADER), 1, fp);

    int i;
    if (bit == 8)
    {
        if (rgb == NULL)
        {
            rgb = (RGB_QUAD*)malloc(sizeof(RGB_QUAD)*256);
            for (i = 0; i < 256; i++)
            {
                rgb[i].rgbRed      = i;
                rgb[i].rgbGreen    = i;
                rgb[i].rgbBlue     = i;
                rgb[i].rgbReserved = 0;
            }
        }
        fwrite(rgb, sizeof(RGB_QUAD), 256, fp);
        if (bmp->rgb == NULL)
            free(rgb);
    }

    unsigned char *p_align = (unsigned char*)calloc(align, sizeof(unsigned char));
    for (i = 0; i<height; i++)
    {
        fwrite(bmpData, lineByte, 1, fp);
        bmpData += lineByte;
        fwrite(p_align, align, 1, fp);
    }
    free(p_align);

    fclose(fp);
    return 0;
}

void FreeBMP_Image(BMP_IMAGE * bmp)
{
    if (bmp != NULL)
    {
        free(bmp->rgb);
        bmp->rgb = NULL;

        free(bmp->bmpData);
        bmp->bmpData = NULL;
    }
	
    free(bmp);
    bmp = NULL;
}

int WriteBMPColorMap(const char *bmpfile,
                     int width, int height,
                     double **value,
                     RGB_QUAD(*ColorMap)(double))
{
#if PRINTLEVEL > 3
    printf("begin write bmp...");
#endif
    FILE *fp = NULL;
    //fopen_s(&fp, bmpfile, "wb");
    fp = fopen(bmpfile, "wb");
    if (fp == NULL)
    {
        printf("fp == NULL\n");
        return 1;
    }

    int bit = 24;
    int lineByte = width * bit / 8;
    int align = lineByte % 4;
    if (align != 0)
        align = 4 - align;
#if 0
    printf("lineByte = %d\n", lineByte);
    printf("align = %d\n", align);
#endif
    BIT_MAP_FILE_HEADER bmpFileHeader;
    bmpFileHeader.bfType      = 0x4D42;
    bmpFileHeader.bfSize      = 54 + (lineByte + align)*height;
    bmpFileHeader.bfReserved1 = 0;
    bmpFileHeader.bfReserved2 = 0;
    bmpFileHeader.bfOffBits   = 54;

    BIT_MAP_INFO_HEADER bmpInfoHeader;
    bmpInfoHeader.biSize          = 40;
    bmpInfoHeader.biWidth         = width;
    bmpInfoHeader.biHeight        = height;
    bmpInfoHeader.biPlanes        = 1;
    bmpInfoHeader.biBitCount      = bit;
    bmpInfoHeader.biCompression   = 0;
    bmpInfoHeader.biSizeImage     = (lineByte + align)*height;
    bmpInfoHeader.biXPelsPerMeter = 0;
    bmpInfoHeader.biYPelsPerMeter = 0;
    bmpInfoHeader.biClrUsed       = 0;
    bmpInfoHeader.biClrImportant  = 0;

    fwrite(&bmpFileHeader, sizeof(BIT_MAP_FILE_HEADER), 1, fp);
    fwrite(&bmpInfoHeader, sizeof(BIT_MAP_INFO_HEADER), 1, fp);

    int i, j;
    RGB_QUAD color;
    //unsigned char *rgb = (unsigned char*)malloc(sizeof(unsigned char) * 3);
    unsigned char *p_align = (unsigned char*)calloc(align, sizeof(unsigned char));
    for (i = 0; i<height; i++)
    {
        for (j = 0; j < width; j++)
        {
            color = ColorMap(value[i][j]);
            fwrite(&color, 3, 1, fp);
        }
        fwrite(p_align, align, 1, fp);
    }
    free(p_align);

    fclose(fp);
#if PRINTLEVEL > 3
    printf("done.\n");
#endif
    return 0;
}

RGB_QUAD ColorMapSummer(double value)
{
    RGB_QUAD rgb;
    rgb.rgbRed      = (unsigned char)255.0*value;
    rgb.rgbGreen    = (unsigned char)128.0 + (255.0 - 128.0)*value;
    rgb.rgbBlue     = (unsigned char)102.0;
    rgb.rgbReserved = 0;
    return rgb;
}

RGB_QUAD ColorMapCustom(double value)
{
    RGB_QUAD rgb;
    rgb.rgbRed      = (unsigned char)255.0*value*value;
    rgb.rgbGreen    = (unsigned char)(255.0 - 128.0)*(value + value*value);
    rgb.rgbBlue     = (unsigned char)100.0*value;
    rgb.rgbReserved = 0;
    return rgb;
}

RGB_QUAD ColorMapJat(double value)
{
    RGB_QUAD rgb;
    if (value * 64 < 8)
    {
        rgb.rgbRed = 0;
        rgb.rgbGreen = 0;
        rgb.rgbBlue = 143 + (255 - 143)*value * 64 / 8;
    }
    else if (value * 64 < 24)
    {
        rgb.rgbRed = 0;
        rgb.rgbGreen = 255 * (value * 64 - 8) / (24 - 8);
        rgb.rgbBlue = 255;
    }
    else if (value * 64 < 40)
    {
        rgb.rgbRed = 255 * (value * 64 - 24) / (40 - 24);
        rgb.rgbGreen = 255;
        rgb.rgbBlue = 255 - 255 * (value * 64 - 24) / (40 - 24);
    }
    else if (value * 64 < 56)
    {
        rgb.rgbRed = 255;
        rgb.rgbGreen = 255 - 255 * (value * 64 - 40) / (56 - 40);
        rgb.rgbBlue = 0;
    }
    else
    {
        rgb.rgbRed = 255 + (128 - 255)*(value * 64 - 56) / (64 - 56);
        rgb.rgbGreen = 0;
        rgb.rgbBlue = 0;
    }
    rgb.rgbReserved = 0;
    return rgb;
}

RGB_QUAD ColorMapMat(double value)
{
    RGB_QUAD rgb;
    if(value > 0)
    {
        rgb.rgbRed      = (unsigned char)255.0;
        rgb.rgbGreen    = (unsigned char)255.0*(1-value*value);
        rgb.rgbBlue     = (unsigned char)255.0*cos(value*M_PI/2);
    }
    else if(value < 0)
    {
        //value
        rgb.rgbRed      = (unsigned char)255.0*cos(value*M_PI/2);
        rgb.rgbGreen    = (unsigned char)255.0*(1-value*value);
        rgb.rgbBlue     = (unsigned char)255.0;
    }
    else
    {
        rgb.rgbRed      = (unsigned char)0.0;
        rgb.rgbGreen    = (unsigned char)0.0;
        rgb.rgbBlue     = (unsigned char)0.0;
    }
    rgb.rgbReserved = 0;
    return rgb;
}

RGB_QUAD ColorMapMatStruct(double value)
{
    RGB_QUAD rgb;
    if(value > 0)
    {
        rgb.rgbRed      = (unsigned char)255.0;
        rgb.rgbGreen    = (unsigned char)0.0;
        rgb.rgbBlue     = (unsigned char)0.0;
    }
    else if(value < 0)
    {
        //value
        rgb.rgbRed      = (unsigned char)0.0;
        rgb.rgbGreen    = (unsigned char)255.0;
        rgb.rgbBlue     = (unsigned char)0.0;
    }
    else
    {
        rgb.rgbRed      = (unsigned char)0.0;
        rgb.rgbGreen    = (unsigned char)0.0;
        rgb.rgbBlue     = (unsigned char)0.0;
    }
    rgb.rgbReserved = 0;
    return rgb;
}
