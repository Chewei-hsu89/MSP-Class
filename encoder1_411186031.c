#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 8
#define CHANNELS 3  // 0->Y, 1->Cb, 2->Cr

#pragma pack(push,1)
/* --- BMP File Header (14 bytes) --- */
typedef struct {
    unsigned short bfType;      // 'BM' = 0x4D42
    unsigned int   bfSize;      // 檔案大小 (含 Header)
    unsigned short bfReserved1;
    unsigned short bfReserved2;
    unsigned int   bfOffBits;   // 圖像資料在檔案的起始位址 (通常 54)
} BMPFileHeader;

/* --- BMP Info Header (40 bytes, Windows 3.x) --- */
typedef struct {
    unsigned int   biSize;          
    int            biWidth;
    int            biHeight;
    unsigned short biPlanes;        
    unsigned short biBitCount;      // 24 => 24bits
    unsigned int   biCompression;   // 0 => BI_RGB (無壓縮)
    unsigned int   biSizeImage;     
    int            biXPelsPerMeter; 
    int            biYPelsPerMeter;
    unsigned int   biClrUsed;
    unsigned int   biClrImportant;
} BMPInfoHeader;
#pragma pack(pop)

/* 
   ====================================================================
   1) 讀取 24-bit BMP (無壓縮)
   ====================================================================
*/
unsigned char* read_bmp24(const char *filename, int *width, int *height)
{
    FILE *fp = fopen(filename, "rb");
    if(!fp) {
        fprintf(stderr, "Cannot open BMP: %s\n", filename);
        return NULL;
    }

    BMPFileHeader fileHeader;
    if(fread(&fileHeader, sizeof(fileHeader), 1, fp)!=1){
        fprintf(stderr,"Error read BMPFileHeader\n");
        fclose(fp);
        return NULL;
    }

    BMPInfoHeader infoHeader;
    if(fread(&infoHeader, sizeof(infoHeader), 1, fp)!=1){
        fprintf(stderr,"Error read BMPInfoHeader\n");
        fclose(fp);
        return NULL;
    }

    if(fileHeader.bfType != 0x4D42){  // 'BM'
        fprintf(stderr,"Not a valid BMP (bfType=0x%X)\n", fileHeader.bfType);
        fclose(fp);
        return NULL;
    }
    if(infoHeader.biBitCount!=24 || infoHeader.biCompression!=0){
        fprintf(stderr,
          "Only support 24-bit uncompressed BMP (biBitCount=%d, biCompression=%u)\n",
          infoHeader.biBitCount, infoHeader.biCompression);
        fclose(fp);
        return NULL;
    }

    int w = infoHeader.biWidth;
    int h = infoHeader.biHeight;
    int topDown = 0;
    if(h<0){
        h = -h;
        topDown=1;
    }
    *width  = w;
    *height = h;

    int rowSize = (w*3 +3)&~3;
    unsigned char* bmpData=(unsigned char*)malloc(rowSize*h);
    if(!bmpData){
        fprintf(stderr,"malloc fail for bmpData\n");
        fclose(fp);
        return NULL;
    }

    fseek(fp, fileHeader.bfOffBits, SEEK_SET);
    if(fread(bmpData,1,rowSize*h, fp) != (size_t)(rowSize*h)){
        fprintf(stderr,"fail read BMP pixel data\n");
        free(bmpData);
        fclose(fp);
        return NULL;
    }
    fclose(fp);

    if(topDown){
        // flip
        unsigned char* tempRow=(unsigned char*)malloc(rowSize);
        if(!tempRow){
            fprintf(stderr,"malloc fail for topDown flip\n");
            free(bmpData);
            return NULL;
        }
        for(int row=0; row<h/2; row++){
            unsigned char* rowA = bmpData + rowSize* row;
            unsigned char* rowB = bmpData + rowSize*(h-1-row);
            memcpy(tempRow, rowA, rowSize);
            memcpy(rowA, rowB, rowSize);
            memcpy(rowB, tempRow, rowSize);
        }
        free(tempRow);
    }

    return bmpData;
}

/*
   ====================================================================
   2) RGB -> YCbCr (JPEG公式)
   ====================================================================
*/
void rgb_to_ycbcr(unsigned char R,unsigned char G,unsigned char B,
                  float* Y, float* Cb, float* Cr)
{
    *Y  =  0.299f*R + 0.587f*G + 0.114f*B;
    *Cb = -0.168736f*R - 0.331264f*G + 0.5f*B +128.f;
    *Cr =  0.5f*R - 0.418688f*G - 0.081312f*B +128.f;
}

/*
   ====================================================================
   3) DCT
   ====================================================================
*/
void dct_transform(float input[N][N], float output[N][N]){
    for(int u=0; u<N; u++){
        for(int v=0; v<N; v++){
            float sum=0.f;
            for(int x=0;x<N;x++){
                for(int y=0;y<N;y++){
                    sum += input[x][y]
                         * cosf((M_PI/N)*(x+0.5f)*u)
                         * cosf((M_PI/N)*(y+0.5f)*v);
                }
            }
            float alpha_u=(u==0)? sqrtf(1.f/N): sqrtf(2.f/N);
            float alpha_v=(v==0)? sqrtf(1.f/N): sqrtf(2.f/N);
            output[u][v]= alpha_u*alpha_v* sum;
        }
    }
}

/*
   ====================================================================
   4) 處理 8×8 區塊 (DCT + 量化 + 誤差)
   ====================================================================
*/
void process_block(float block[N][N], float qt[N][N],
                   float F_out[N][N],
                   short qF_out[N*N],
                   float eF_out[N*N])
{
    float dct_blk[N][N];
    dct_transform(block, dct_blk);
    for(int u=0;u<N;u++){
        for(int v=0;v<N;v++){
            float F_val= dct_blk[u][v];
            float Q_val= qt[u][v];
            short q_val= (short)roundf(F_val/ Q_val);
            float e_val= F_val - q_val*Q_val;

            F_out[u][v]      = F_val;
            qF_out[u*N + v]  = q_val;
            eF_out[u*N + v]  = e_val;
        }
    }
}

/*
   ====================================================================
   5) 寫檔 (qF, eF)
   ====================================================================
*/
void save_quantized_to_file(const short* arr,int totalBlocks,const char* fname)
{
    FILE* fp=fopen(fname,"wb");
    if(!fp){
        fprintf(stderr,"Cannot open %s\n", fname);
        exit(1);
    }
    fwrite(arr, sizeof(short), totalBlocks*N*N, fp);
    fclose(fp);
    printf("Saved quantized data to %s\n", fname);
}
void save_error_to_file(const float* arr,int totalBlocks,const char* fname)
{
    FILE* fp=fopen(fname,"wb");
    if(!fp){
        fprintf(stderr,"Cannot open %s\n", fname);
        exit(1);
    }
    fwrite(arr, sizeof(float), totalBlocks*N*N, fp);
    fclose(fp);
    printf("Saved error data to %s\n", fname);
}

/*
   ====================================================================
   6) 產生量化表 & 輸出
   ====================================================================
*/
void generate_quantization_table(float qt[N][N], int offset)
{
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            qt[i][j] = 1.f + (float)(i+j+offset);
        }
    }
}
void save_quantization_table_to_txt(const float qt[N][N],const char* fname)
{
    FILE* fp=fopen(fname,"w");
    if(!fp){
        fprintf(stderr,"Cannot open %s\n", fname);
        exit(1);
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            fprintf(fp,"%.3f", qt[i][j]);
            if(j<N-1) fprintf(fp," ");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    printf("Saved quant table to %s\n", fname);
}

/*
   ====================================================================
   7) SQNR 計算：3×64
       sumF[ch][k], sumE[ch][k] => SQNR[k] = 10log10(sumF/sumE)
   ====================================================================
*/
static double sumF[CHANNELS][64];
static double sumE[CHANNELS][64];

void init_SQNR_arrays()
{
    for(int ch=0; ch<CHANNELS; ch++){
        for(int k=0; k<64; k++){
            sumF[ch][k]=0.0;
            sumE[ch][k]=0.0;
        }
    }
}
void accumulate_SQNR(int ch, const float F_out[8][8], const float eF_out[64])
{
    for(int k=0;k<64;k++){
        int u = k/8;
        int v = k%8;
        double valF = (double)F_out[u][v];
        double valE = (double)eF_out[k];
        sumF[ch][k] += valF*valF;
        sumE[ch][k] += valE*valE;
    }
}
void print_SQNR()
{
    printf("\n===== SQNR(dB) for each frequency bin (u,v) =====\n");
    const char* chName[3]={"Y","Cb","Cr"};
    for(int ch=0; ch<CHANNELS; ch++){
        printf("Channel: %s\n", chName[ch]);
        for(int k=0;k<64;k++){
            double sF= sumF[ch][k];
            double sE= sumE[ch][k];
            double sqnr = (sE==0.0)? 999.99 : 10.0*log10(sF/sE);

            int u = k/8;
            int v = k%8;
            printf(" (%d,%d): %7.2f dB", u,v, sqnr);
            if((k+1)%8==0) printf("\n");
        }
        printf("\n");
    }
}

/*
   ====================================================================
   8) main():
   假設呼叫方式：
   ./encoder1 1 <BMPfile> <QtY> <QtCb> <QtCr> <dim.txt>
              <qF_Y.raw> <qF_Cb.raw> <qF_Cr.raw> <eF_Y.raw> <eF_Cb.raw> <eF_Cr.raw>
   總共有 13 個 argv
   => argv[1] = "1" 保留無用
   => argv[2] = BMPfile
   => argv[3..5] = QtY/Cb/Cr
   => argv[6] = dim.txt
   => argv[7..9] = qF_*.raw
   => argv[10..12] = eF_*.raw
   ====================================================================
*/
int main(int argc, char** argv)
{
    if(argc!=13){
        fprintf(stderr,"Usage: %s 1 <BMP> <QtY> <QtCb> <QtCr> <dim.txt> <qF_Y> <qF_Cb> <qF_Cr> <eF_Y> <eF_Cb> <eF_Cr>\n", argv[0]);
        return 1;
    }
    // 跳過 argv[1]="1"
    char** realargv = argv+2;
    // realargv[0] => BMPfile
    // realargv[1] => QtY.txt
    // realargv[2] => QtCb.txt
    // realargv[3] => QtCr.txt
    // realargv[4] => dim.txt
    // realargv[5] => qF_Y.raw
    // realargv[6] => qF_Cb.raw
    // realargv[7] => qF_Cr.raw
    // realargv[8] => eF_Y.raw
    // realargv[9] => eF_Cb.raw
    // realargv[10]=> eF_Cr.raw

    const char* bmpFile   = realargv[0];
    const char* qtYFile   = realargv[1];
    const char* qtCbFile  = realargv[2];
    const char* qtCrFile  = realargv[3];
    const char* dimFile   = realargv[4];
    const char* qFYFile   = realargv[5];
    const char* qFCbFile  = realargv[6];
    const char* qFCrFile  = realargv[7];
    const char* eFYFile   = realargv[8];
    const char* eFCbFile  = realargv[9];
    const char* eFCrFile  = realargv[10];

    // (1) 讀 dim.txt => (width, height)
    FILE* fDim=fopen(dimFile,"r");
    if(!fDim){
        fprintf(stderr,"Cannot open dim file: %s\n", dimFile);
        return 1;
    }
    int width, height;
    fscanf(fDim,"%d %d", &width, &height);
    fclose(fDim);
    printf("From %s => width=%d, height=%d\n", dimFile,width,height);

    // (2) 讀 BMP => "Kimberly.bmp"
    int bmpW,bmpH;
    unsigned char* bmpData = read_bmp24(bmpFile, &bmpW, &bmpH);
    if(!bmpData){
        fprintf(stderr,"Fail read BMP: %s\n", bmpFile);
        return 1;
    }
    printf("BMP => w=%d, h=%d\n", bmpW,bmpH);

    // 可選: 以 BMP => override dim.txt
    width  = bmpW;
    height = bmpH;

    // (3) padding
    int paddedWidth  = ((width +(N-1))/N)*N;
    int paddedHeight = ((height+(N-1))/N)*N;
    int blockCountX= paddedWidth/N;
    int blockCountY= paddedHeight/N;
    int totalBlocks= blockCountX* blockCountY;

    float** Y_plane =(float**)malloc(paddedHeight*sizeof(float*));
    float** Cb_plane=(float**)malloc(paddedHeight*sizeof(float*));
    float** Cr_plane=(float**)malloc(paddedHeight*sizeof(float*));

    for(int i=0;i<paddedHeight;i++){
        Y_plane[i]  =(float*)malloc(paddedWidth*sizeof(float));
        Cb_plane[i] =(float*)malloc(paddedWidth*sizeof(float));
        Cr_plane[i] =(float*)malloc(paddedWidth*sizeof(float));
        for(int j=0;j<paddedWidth;j++){
            Y_plane[i][j]  = 0.f;
            Cb_plane[i][j] =128.f;
            Cr_plane[i][j] =128.f;
        }
    }

    // (4) BMP -> YCbCr
    int rowSize = (bmpW*3 +3)&~3;
    for(int r=0; r<bmpH; r++){
        for(int c=0; c<bmpW; c++){
            int bmp_y= (bmpH-1 - r);
            int idx= bmp_y*rowSize + c*3;
            unsigned char B= bmpData[idx+0];
            unsigned char G= bmpData[idx+1];
            unsigned char R= bmpData[idx+2];

            float Yv,Cbv,Crv;
            rgb_to_ycbcr(R,G,B, &Yv, &Cbv, &Crv);
            Y_plane[r][c]  = Yv;
            Cb_plane[r][c] = Cbv;
            Cr_plane[r][c] = Crv;
        }
    }
    free(bmpData);

    // (5) 建立量化表
    float qtY [N][N], qtCb[N][N], qtCr[N][N];
    generate_quantization_table(qtY,  0);
    generate_quantization_table(qtCb, 5);
    generate_quantization_table(qtCr,10);

    // 輸出量化表 => QtY, QtCb, QtCr (若你要保留)
    save_quantization_table_to_txt(qtY ,  qtYFile);
    save_quantization_table_to_txt(qtCb,  qtCbFile);
    save_quantization_table_to_txt(qtCr,  qtCrFile);

    // 分配 qF,eF
    short *qF_Y  =(short*)malloc(sizeof(short)* totalBlocks*N*N);
    short *qF_Cb =(short*)malloc(sizeof(short)* totalBlocks*N*N);
    short *qF_Cr =(short*)malloc(sizeof(short)* totalBlocks*N*N);

    float *eF_Y  =(float*)malloc(sizeof(float)* totalBlocks*N*N);
    float *eF_Cb =(float*)malloc(sizeof(float)* totalBlocks*N*N);
    float *eF_Cr =(float*)malloc(sizeof(float)* totalBlocks*N*N);

    // (6) 為計算 SQNR 累加, 建立 sumF, sumE
    static double sumF[CHANNELS][64];
    static double sumE[CHANNELS][64];
    for(int ch=0; ch<CHANNELS; ch++){
        for(int k=0; k<64; k++){
            sumF[ch][k] = 0.0;
            sumE[ch][k] = 0.0;
        }
    }

    // 逐區塊 => DCT + 量化 => 累加 SQNR
    int idxY=0, idxCb=0, idxCr=0;
    for(int by=0; by<blockCountY; by++){
        for(int bx=0; bx<blockCountX; bx++){
            // ---------------- Y block ----------------
            float blockY[N][N];
            for(int x=0;x<N;x++){
                for(int y=0;y<N;y++){
                    blockY[x][y]= Y_plane[ by*N + x ][ bx*N + y ];
                }
            }
            float F_Y[N][N];
            short qF_Y_blk[64];
            float eF_Y_blk[64];
            process_block(blockY, qtY, F_Y, qF_Y_blk, eF_Y_blk);
            for(int k=0;k<64;k++){
                qF_Y[idxY*64 +k] = qF_Y_blk[k];
                eF_Y[idxY*64 +k] = eF_Y_blk[k];
            }
            // 累加 SQNR => ch=0 (Y)
            for(int k=0;k<64;k++){
                int u=k/8, v=k%8;
                double valF = (double)F_Y[u][v];
                double valE = (double)eF_Y_blk[k];
                sumF[0][k]+= valF*valF;
                sumE[0][k]+= valE*valE;
            }
            idxY++;

            // ---------------- Cb block ----------------
            float blockCb[N][N];
            for(int x=0;x<N;x++){
                for(int y=0;y<N;y++){
                    blockCb[x][y]= Cb_plane[ by*N + x ][ bx*N + y ];
                }
            }
            float F_Cb[N][N];
            short qF_Cb_blk[64];
            float eF_Cb_blk[64];
            process_block(blockCb, qtCb, F_Cb, qF_Cb_blk, eF_Cb_blk);
            for(int k=0;k<64;k++){
                qF_Cb[idxCb*64 +k] = qF_Cb_blk[k];
                eF_Cb[idxCb*64 +k] = eF_Cb_blk[k];
            }
            // 累加 SQNR => ch=1 (Cb)
            for(int k=0;k<64;k++){
                int u=k/8, v=k%8;
                double valF = (double)F_Cb[u][v];
                double valE = (double)eF_Cb_blk[k];
                sumF[1][k]+= valF*valF;
                sumE[1][k]+= valE*valE;
            }
            idxCb++;

            // ---------------- Cr block ----------------
            float blockCr[N][N];
            for(int x=0;x<N;x++){
                for(int y=0;y<N;y++){
                    blockCr[x][y]= Cr_plane[ by*N + x ][ bx*N + y ];
                }
            }
            float F_Cr[N][N];
            short qF_Cr_blk[64];
            float eF_Cr_blk[64];
            process_block(blockCr, qtCr, F_Cr, qF_Cr_blk, eF_Cr_blk);
            for(int k=0;k<64;k++){
                qF_Cr[idxCr*64 +k] = qF_Cr_blk[k];
                eF_Cr[idxCr*64 +k] = eF_Cr_blk[k];
            }
            // 累加 SQNR => ch=2 (Cr)
            for(int k=0;k<64;k++){
                int u=k/8, v=k%8;
                double valF = (double)F_Cr[u][v];
                double valE = (double)eF_Cr_blk[k];
                sumF[2][k]+= valF*valF;
                sumE[2][k]+= valE*valE;
            }
            idxCr++;
        }
    }

    // (7) 寫出 qF_*.raw & eF_*.raw
    save_quantized_to_file(qF_Y,  idxY,  qFYFile);
    save_quantized_to_file(qF_Cb, idxCb, qFCbFile);
    save_quantized_to_file(qF_Cr, idxCr, qFCrFile);

    save_error_to_file(eF_Y,  idxY,  eFYFile);
    save_error_to_file(eF_Cb, idxCb, eFCbFile);
    save_error_to_file(eF_Cr, idxCr, eFCrFile);

    // (8) 計算並印出 3×64 SQNR
    printf("\n===== SQNR(dB) for each frequency bin (u,v) =====\n");
    const char* chName[3]={"Y","Cb","Cr"};
    for(int ch=0; ch<CHANNELS; ch++){
        printf("Channel: %s\n", chName[ch]);
        for(int k=0;k<64;k++){
            double sF= sumF[ch][k];
            double sE= sumE[ch][k];
            double sqnr= (sE==0.0)? 999.99 : 10.0*log10(sF/sE);

            int u=k/8, v=k%8;
            printf("  (%d,%d): %7.2f dB", u,v, sqnr);
            if((k+1)%8==0) printf("\n");
        }
        printf("\n");
    }

    // 釋放
    free(qF_Y);  free(qF_Cb);  free(qF_Cr);
    free(eF_Y);  free(eF_Cb);  free(eF_Cr);

    for(int i=0;i<paddedHeight;i++){
        free(Y_plane[i]);
        free(Cb_plane[i]);
        free(Cr_plane[i]);
    }
    free(Y_plane);
    free(Cb_plane);
    free(Cr_plane);

    printf("\nEncoder finished.\n");
    return 0;
}
