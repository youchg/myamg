#ifndef __SIMPLE_BMP__
#define __SIMPLE_BMP__

#include <stdint.h>
#include <stdio.h>

typedef  uint8_t  UINT8;
typedef uint16_t UINT16;
typedef uint32_t UINT32;
typedef  int16_t  INT16;
typedef  int32_t  INT32;

//�ļ���Ϣͷ��14���ֽ�
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
bfType      ˵���ļ������ͣ�����Ϊ0x4D42��ASCII��Ϊ'BM'

bfSize      ˵����λͼ�ļ��Ĵ�С�����ֽ�Ϊ��λ

bfReserved1 ��������������Ϊ0
bfReserved2 ��������������Ϊ0

bfOffBits   ˵�����ļ�ͷ��ʼ��ʵ��ͼ������֮����ֽ�ƫ������
            λͼ��Ϣͷ�͵�ɫ��ĳ��Ȼ���ݲ�ͬ������仯��
            ��˿�����ƫ��ֵѸ�ٵش��ļ��ж�ȡ��λ����
*/

//λͼ��Ϣͷ��40���ֽ�
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
biSize ˵��BitMapInfoHeader�ṹ����Ҫ���ֽ���

biWidth ˵��ͼ��Ŀ�ȣ�������Ϊ��λ

biHeight        һ���ô���˵��ͼ��ĸ߶ȣ�������Ϊ��λ��
                ��������һ���ô���ָ����ͼ����
                �����λͼ���������λͼ��
                �����ֵ��һ��������˵��ͼ���ǵ���ģ�
                �����ֵ��һ��������˵��ͼ��������ġ�
                �������BMP�ļ����ǵ����λͼ�����߶�Ϊ����

biPlanes        ΪĿ���豸˵��λ����������Ϊ1

biBitCount      ˵��������/���أ���ֵΪ1��4��8��16��24��32.
                ���󲿷�ͼ����24λ��32λ�ġ�

biCompression   ˵��ͼ������ѹ�������͡���ѹ������Ϊ0(BI_RGB)

biSizeImage     ˵��ͼ��Ĵ�С�����ֽ�Ϊ��λ��
                ��С�����ļ��ֽ�����ȥ�ļ�ͷ��С���ټ�ȥ��ɫ���С��
                ����BI_RGB��ʽʱ��������Ϊ0

biXPelsPerMeter ˵��ˮƽ�ֱ��ʣ�������/�ױ�ʾ
biYPelsPerMeter ˵����ֱ�ֱ��ʣ�������/�ױ�ʾ

biClrUsed       ˵��λͼʵ��ʹ�õĲ�ɫ���е���ɫ��������
                0��ʾʹ�����еĵ�ɫ����

biClrImportant  ˵����ͼ����ʾ����ҪӰ�����ɫ������Ŀ��
                0��ʾ����Ҫ
*/

typedef struct RGB_QUAD_
{
    UINT8 rgbBlue;     //����ɫ����ɫ����
    UINT8 rgbGreen;    //����ɫ����ɫ����  
    UINT8 rgbRed;      //����ɫ�ĺ�ɫ����  
    UINT8 rgbReserved; //����ֵ  
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
