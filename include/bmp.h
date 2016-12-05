#ifndef __SIMPLE_BMP__
#define __SIMPLE_BMP__

#include <stdint.h>
#include <stdio.h>

typedef  uint8_t  UINT8;
typedef uint16_t UINT16;
typedef uint32_t UINT32;
typedef  int16_t  INT16;
typedef  int32_t  INT32;

//文件信息头，14个字节
#pragma pack(push,2)
typedef struct BIT_MAP_FILE_HEADER_
{
    UINT16 bfType;
    UINT32 bfSize;
    UINT16 bfReserved1;
    UINT16 bfReserved2;
    UINT32 bfOffBits;
}BIT_MAP_FILE_HEADER;
#pragma pack(pop)
/*
bfType      说明文件的类型，必需为0x4D42，ASCII码为'BM'

bfSize      说明该位图文件的大小，以字节为单位

bfReserved1 保留，必需设置为0
bfReserved2 保留，必需设置为0

bfOffBits   说明从文件头开始到实际图像数据之间的字节偏移量。
            位图信息头和调色板的长度会根据不同情况而变化，
            因此可以用偏移值迅速地从文件中读取到位数据
*/

//位图信息头，40个字节
typedef struct BIT_MAP_INFO_HEADER_
{
    UINT32 biSize;
     INT32 biWidth;
     INT32 biHeight;
    UINT16 biPlanes;
    UINT16 biBitCount;
    UINT32 biCompression;
    UINT32 biSizeImage;
     INT32 biXPelsPerMeter;
     INT32 biYPelsPerMeter;
    UINT32 biClrUsed;
    UINT32 biClrImportant;
}BIT_MAP_INFO_HEADER;
/*
biSize 说明BitMapInfoHeader结构所需要的字节数

biWidth 说明图像的宽度，以像素为单位

biHeight        一个用处是说明图像的高度，以像素为单位。
                还有另外一个用处，指明该图像是
                倒向的位图还是正向的位图。
                如果该值是一个正数，说明图像是倒向的；
                如果该值是一个负数，说明图像是正向的。
                大多数的BMP文件都是倒向的位图，即高度为正。

biPlanes        为目标设备说明位面数，总是为1

biBitCount      说明比特数/像素，其值为1、4、8、16、24、32.
                绝大部分图像是24位和32位的。

biCompression   说明图像数据压缩的类型。不压缩设置为0(BI_RGB)

biSizeImage     说明图像的大小，以字节为单位。
                大小等于文件字节数减去文件头大小，再减去颜色表大小。
                当用BI_RGB格式时，可设置为0

biXPelsPerMeter 说明水平分辨率，用像素/米表示
biYPelsPerMeter 说明垂直分辨率，用像素/米表示

biClrUsed       说明位图实际使用的彩色表中的颜色索引数，
                0表示使用所有的调色板项

biClrImportant  说明对图像显示有重要影响的颜色索引数目，
                0表示都重要
*/

typedef struct RGB_QUAD_
{
    UINT8 rgbBlue;     //该颜色的蓝色分量
    UINT8 rgbGreen;    //该颜色的绿色分量  
    UINT8 rgbRed;      //该颜色的红色分量  
    UINT8 rgbReserved; //保留值  
}RGB_QUAD;

typedef struct BMP_IMAGE_
{
    int width;
    int height;
    int bit;
    RGB_QUAD *rgb;
    unsigned char *bmpData;
}BMP_IMAGE;

BMP_IMAGE* ReadBMP(const char *bmpfile);

int WriteBMP(const char *bmpfile, BMP_IMAGE *bmp);

int WriteBMPColorMap(const char *bmpfile,
                     int width, int height,
                     double **value,
                     RGB_QUAD(*ColorMap)(double));

void FreeBMP_Image(BMP_IMAGE * bmp);

/* 0 <= value <= 1 */
RGB_QUAD ColorMapSummer(double value);
RGB_QUAD ColorMapJat   (double value);
RGB_QUAD ColorMapCustom(double value);
RGB_QUAD ColorMapMat(double value);
RGB_QUAD ColorMapMatStruct(double value);
#endif
