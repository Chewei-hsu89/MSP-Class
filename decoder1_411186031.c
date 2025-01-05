#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 8

// ---------------------------------------------------------
// clamp 將 double 轉為 [0,255] 之間的 unsigned char
// ---------------------------------------------------------
static inline unsigned char clamp(double val){
    if(val < 0.0)   return 0;
    if(val > 255.0) return 255;
    return (unsigned char)(val + 0.5); // 四捨五入
}

// ---------------------------------------------------------
// 讀取量化表 (ASCII 8×8)
// ---------------------------------------------------------
void load_quantization_table(const char* filename, float qt[N][N]){
    FILE* fp = fopen(filename,"r");
    if(!fp){
        fprintf(stderr,"Cannot open quantization table file: %s\n", filename);
        exit(1);
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(fscanf(fp, "%f", &qt[i][j]) != 1){
                fprintf(stderr,"Error reading quantization table: %s\n", filename);
                fclose(fp);
                exit(1);
            }
        }
    }
    fclose(fp);
}

// ---------------------------------------------------------
// 讀取 .raw (short / float)
// ---------------------------------------------------------
void read_raw_short(const char* fname, short* buffer, size_t length){
    FILE* fp = fopen(fname,"rb");
    if(!fp){
        fprintf(stderr,"Cannot open file: %s\n", fname);
        exit(1);
    }
    size_t rc = fread(buffer, sizeof(short), length, fp);
    if(rc != length){
        fprintf(stderr,"Error: read %zu short, expected %zu in %s\n", rc,length, fname);
        fclose(fp);
        exit(1);
    }
    fclose(fp);
}
void read_raw_float(const char* fname, float* buffer, size_t length){
    FILE* fp = fopen(fname,"rb");
    if(!fp){
        fprintf(stderr,"Cannot open file: %s\n", fname);
        exit(1);
    }
    size_t rc = fread(buffer, sizeof(float), length, fp);
    if(rc != length){
        fprintf(stderr,"Error: read %zu float, expected %zu in %s\n", rc,length, fname);
        fclose(fp);
        exit(1);
    }
    fclose(fp);
}

// ---------------------------------------------------------
// iDCT (double 精度計算，減少誤差)
// ---------------------------------------------------------
void idct_transform_double(const double in[8][8], double out[8][8])
{
    for(int x=0; x<8; x++){
        for(int y=0; y<8; y++){
            double sum = 0.0;
            for(int u=0; u<8; u++){
                for(int v=0; v<8; v++){
                    double alpha_u = (u==0)? sqrt(1.0/8.0) : sqrt(2.0/8.0);
                    double alpha_v = (v==0)? sqrt(1.0/8.0) : sqrt(2.0/8.0);
                    double cos_xu  = cos((M_PI/8.0)*((double)x + 0.5)*(double)u);
                    double cos_yv  = cos((M_PI/8.0)*((double)y + 0.5)*(double)v);

                    sum += alpha_u * alpha_v * in[u][v] * cos_xu * cos_yv;
                }
            }
            out[x][y] = sum;
        }
    }
}

// ---------------------------------------------------------
// YCbCr -> RGB (BT.601 / JPEG) - double 版
// 對應到 encoder 的 rgb_to_ycbcr
//   R = Y + 1.402   *(Cr - 128)
//   G = Y - 0.34414 *(Cb - 128) - 0.71414*(Cr -128)
//   B = Y + 1.772   *(Cb - 128)
// ---------------------------------------------------------
void ycbcr_to_rgb_double(double Y, double Cb, double Cr,
                         unsigned char* R, unsigned char* G, unsigned char* B)
{
    double cB = Cb - 128.0;
    double cR = Cr - 128.0;

    // 對應 encoder 端:
    // Y  = 0.299R + 0.587G + 0.114B
    // Cb = -0.168736R -0.331264G + 0.5B +128
    // Cr =  0.5R -0.418688G -0.081312B +128
    // => 反變換:
    double r = Y + 1.402   * cR;
    double g = Y - 0.34414 * cB - 0.71414 * cR;
    double b = Y + 1.772   * cB;

    *R = clamp(r);
    *G = clamp(g);
    *B = clamp(b);
}

// ---------------------------------------------------------
// 寫出 24 bits BMP (簡易版)
// ---------------------------------------------------------
#pragma pack(push,1)
typedef struct {
    unsigned short bfType; 
    unsigned int   bfSize;
    unsigned short bfReserved1;
    unsigned short bfReserved2;
    unsigned int   bfOffBits;
} BMPFileHeader;

typedef struct {
    unsigned int   biSize;
    int            biWidth;
    int            biHeight;
    unsigned short biPlanes;
    unsigned short biBitCount;
    unsigned int   biCompression;
    unsigned int   biSizeImage;
    int            biXPelsPerMeter;
    int            biYPelsPerMeter;
    unsigned int   biClrUsed;
    unsigned int   biClrImportant;
} BMPInfoHeader;
#pragma pack(pop)

void save_bmp(const char* filename, unsigned char* rgbData,
              int width, int height)
{
    FILE* fp = fopen(filename, "wb");
    if(!fp){
        fprintf(stderr,"Cannot open BMP file: %s\n", filename);
        exit(1);
    }
    int rowSize = (width*3 +3)&~3;
    int dataSize= rowSize*height;
    int fileSize= 54 + dataSize;

    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;

    fileHeader.bfType    = 0x4D42; // 'BM'
    fileHeader.bfSize    = fileSize;
    fileHeader.bfReserved1=0;
    fileHeader.bfReserved2=0;
    fileHeader.bfOffBits = 54;

    infoHeader.biSize          = 40;
    infoHeader.biWidth         = width;
    infoHeader.biHeight        = height;
    infoHeader.biPlanes        = 1;
    infoHeader.biBitCount      = 24;
    infoHeader.biCompression   = 0;
    infoHeader.biSizeImage     = dataSize;
    infoHeader.biXPelsPerMeter = 2835;
    infoHeader.biYPelsPerMeter = 2835;
    infoHeader.biClrUsed       = 0;
    infoHeader.biClrImportant  = 0;

    // 寫 header
    fwrite(&fileHeader, sizeof(fileHeader), 1, fp);
    fwrite(&infoHeader, sizeof(infoHeader), 1, fp);

    // 寫像素 (bottom->top)
    for(int y=0; y<height; y++){
        int rowIndex = (height-1-y)*width*3;
        fwrite(&rgbData[rowIndex],1, width*3, fp);
        // padding
        for(int p=width*3; p<rowSize; p++){
            fputc(0,fp);
        }
    }

    fclose(fp);
    printf("BMP saved to %s\n", filename);
}

// ---------------------------------------------------------
// main
//   Usage:
//     (no eF) => 9 參數
//       ./decoder1 <outBMP> <QtY> <QtCb> <QtCr> <dim.txt> <qF_Y> <qF_Cb> <qF_Cr>
//     (with eF) => 12 參數
//       ./decoder1 <outBMP> <QtY> <QtCb> <QtCr> <dim.txt> \
//                  <qF_Y> <qF_Cb> <qF_Cr> <eF_Y> <eF_Cb> <eF_Cr>
// ---------------------------------------------------------
int main(int argc, char** argv)
{
    if(argc!=9 && argc!=12){
        fprintf(stderr,
          "Usage:\n"
          "  (no eF ): %s <outBMP> <QtY.txt> <QtCb.txt> <QtCr.txt> <dim.txt> <qF_Y> <qF_Cb> <qF_Cr>\n"
          "  (with eF):%s <outBMP> <QtY.txt> <QtCb.txt> <QtCr.txt> <dim.txt> <qF_Y> <qF_Cb> <qF_Cr> <eF_Y> <eF_Cb> <eF_Cr>\n",
          argv[0], argv[0]);
        return 1;
    }

    // 是否含誤差檔
    int withErrorFile = (argc==12)? 1 : 0;

    const char* outBMP   = argv[1];
    const char* qtYFile  = argv[2];
    const char* qtCbFile = argv[3];
    const char* qtCrFile = argv[4];
    const char* dimFile  = argv[5];

    const char* qFYFile  = argv[6];
    const char* qFCbFile = argv[7];
    const char* qFCrFile = argv[8];

    const char* eFYFile  = NULL;
    const char* eFCbFile = NULL;
    const char* eFCrFile = NULL;

    if(withErrorFile){
        eFYFile  = argv[9];
        eFCbFile = argv[10];
        eFCrFile = argv[11];
    }

    // (1) 讀取 dim.txt => (width,height)
    FILE* fDim = fopen(dimFile,"r");
    if(!fDim){
        fprintf(stderr,"Cannot open dim file: %s\n", dimFile);
        return 1;
    }
    int realWidth, realHeight;
    fscanf(fDim, "%d %d", &realWidth, &realHeight);
    fclose(fDim);

    // (2) padding => paddedWidth, paddedHeight
    int paddedWidth  = ((realWidth + (N-1))/N)*N;
    int paddedHeight = ((realHeight+(N-1))/N)*N;
    int blockCountX  = paddedWidth / N;
    int blockCountY  = paddedHeight/ N;
    int totalBlocks  = blockCountX* blockCountY;
    size_t totalSize = (size_t)totalBlocks * (N*N);

    // (3) 讀取量化表 => qtY,qtCb,qtCr
    float qtY [N][N], qtCb[N][N], qtCr[N][N];
    load_quantization_table(qtYFile,  qtY);
    load_quantization_table(qtCbFile, qtCb);
    load_quantization_table(qtCrFile, qtCr);

    // (4) 分配 + 讀 qF
    short* qF_Y  = (short*)malloc(sizeof(short)* totalSize);
    short* qF_Cb = (short*)malloc(sizeof(short)* totalSize);
    short* qF_Cr = (short*)malloc(sizeof(short)* totalSize);
    if(!qF_Y || !qF_Cb || !qF_Cr){
        fprintf(stderr,"Memory allocation fail for qF.\n");
        return 1;
    }
    read_raw_short(qFYFile,  qF_Y,  totalSize);
    read_raw_short(qFCbFile, qF_Cb, totalSize);
    read_raw_short(qFCrFile, qF_Cr, totalSize);

    float* eF_Y=NULL;
    float* eF_Cb=NULL;
    float* eF_Cr=NULL;
    if(withErrorFile){
        eF_Y  = (float*)malloc(sizeof(float)* totalSize);
        eF_Cb = (float*)malloc(sizeof(float)* totalSize);
        eF_Cr = (float*)malloc(sizeof(float)* totalSize);
        if(!eF_Y || !eF_Cb || !eF_Cr){
            fprintf(stderr,"Memory alloc fail for eF.\n");
            return 1;
        }
        read_raw_float(eFYFile,  eF_Y,  totalSize);
        read_raw_float(eFCbFile, eF_Cb, totalSize);
        read_raw_float(eFCrFile, eF_Cr, totalSize);
    }

    // (5) 建立 (double) Y,Cb,Cr 陣列
    double** Y_plane  =(double**)malloc(paddedHeight*sizeof(double*));
    double** Cb_plane =(double**)malloc(paddedHeight*sizeof(double*));
    double** Cr_plane =(double**)malloc(paddedHeight*sizeof(double*));

    for(int i=0;i<paddedHeight;i++){
        Y_plane[i]  =(double*)malloc(paddedWidth*sizeof(double));
        Cb_plane[i] =(double*)malloc(paddedWidth*sizeof(double));
        Cr_plane[i] =(double*)malloc(paddedWidth*sizeof(double));
        for(int j=0;j<paddedWidth;j++){
            Y_plane[i][j]  = 0.0;
            Cb_plane[i][j] =128.0;
            Cr_plane[i][j] =128.0;
        }
    }

    // (6) 反量化 + iDCT (double 精度)
    int blockIndex=0;
    for(int by=0; by<blockCountY; by++){
        for(int bx=0; bx<blockCountX; bx++){
            // 重建 F(u,v)
            double F_Y[8][8], F_Cb_blk[8][8], F_Cr_blk[8][8];
            for(int k=0;k<64;k++){
                int u = k/8;
                int v = k%8;
                double qValY  =(double)qF_Y [blockIndex*64 + k];
                double qValCb =(double)qF_Cb[blockIndex*64 + k];
                double qValCr =(double)qF_Cr[blockIndex*64 + k];

                double eValY=0.0, eValCb=0.0, eValCr=0.0;
                if(withErrorFile){
                    eValY  =(double)eF_Y [ blockIndex*64 + k];
                    eValCb =(double)eF_Cb[ blockIndex*64 + k];
                    eValCr =(double)eF_Cr[ blockIndex*64 + k];
                }

                // F(u,v)= qF(u,v)* Qt[u][v] + eF(u,v)
                F_Y[u][v]     = qValY  * qtY [u][v] + eValY;
                F_Cb_blk[u][v]= qValCb * qtCb[u][v] + eValCb;
                F_Cr_blk[u][v]= qValCr * qtCr[u][v] + eValCr;
            }

            // iDCT => blockY,blockCb,blockCr
            double blockY[8][8], blockCb[8][8], blockCr[8][8];
            idct_transform_double(F_Y,      blockY);
            idct_transform_double(F_Cb_blk, blockCb);
            idct_transform_double(F_Cr_blk, blockCr);

            // 放回 plane
            int startY= by*8;
            int startX= bx*8;
            for(int x=0;x<8;x++){
                for(int y=0;y<8;y++){
                    int dstY= startY + x;
                    int dstX= startX + y;
                    if(dstY<paddedHeight && dstX<paddedWidth){
                        Y_plane[dstY][dstX]  = blockY [x][y];
                        Cb_plane[dstY][dstX] = blockCb[x][y];
                        Cr_plane[dstY][dstX] = blockCr[x][y];
                    }
                }
            }
            blockIndex++;
        }
    }

    // (7) YCbCr->RGB => realWidth×realHeight
    unsigned char* rgbData=(unsigned char*)malloc(realWidth*realHeight*3);
    if(!rgbData){
        fprintf(stderr,"Memory error (rgbData).\n");
        return 1;
    }
    for(int r=0;r<realHeight;r++){
        for(int c=0;c<realWidth;c++){
            double Yv = Y_plane[r][c];
            double Cbv= Cb_plane[r][c];
            double Crv= Cr_plane[r][c];

            // 以 "ycbcr_to_rgb_double" 做 color transform
            unsigned char R,G,B;
            ycbcr_to_rgb_double(Yv, Cbv, Crv, &R, &G, &B);

            int idx= (r*realWidth + c)*3;
            rgbData[idx+0]= R;
            rgbData[idx+1]= G;
            rgbData[idx+2]= B;
        }
    }

    // (8) 寫 BMP
    save_bmp(outBMP, rgbData, realWidth, realHeight);
    free(rgbData);

    // 釋放
    free(qF_Y);  free(qF_Cb);  free(qF_Cr);
    if(withErrorFile){
        free(eF_Y); free(eF_Cb); free(eF_Cr);
    }
    for(int i=0;i<paddedHeight;i++){
        free(Y_plane[i]);
        free(Cb_plane[i]);
        free(Cr_plane[i]);
    }
    free(Y_plane);
    free(Cb_plane);
    free(Cr_plane);

    printf("Decoder finished. Output saved to %s\n", outBMP);
    return 0;
}
