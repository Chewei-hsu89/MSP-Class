#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 定義 BMP 結構
typedef struct _bmp {
    int Hpixels;
    int Vpixels;
    unsigned char HeaderInfo[54];
    unsigned long int Hbytes;
    unsigned char **data;
    unsigned char **R;
    unsigned char **G;
    unsigned char **B;
} bmp;

// 函數宣告
int bmp_free(bmp *p_bmp);

void decoder_method_0(char *output_bmp, char *r_file, char *g_file, char *b_file, char *dim_file) {
    FILE *dim_fp = fopen(dim_file, "r");
    if (!dim_fp) {
        fprintf(stderr, "無法開啟尺寸檔案: %s\n", dim_file);
        exit(1);
    }

    bmp image;
    fscanf(dim_fp, "%d %d", &image.Hpixels, &image.Vpixels);
    fclose(dim_fp);

    image.R = (unsigned char **)malloc(image.Vpixels * sizeof(unsigned char *));
    image.G = (unsigned char **)malloc(image.Vpixels * sizeof(unsigned char *));
    image.B = (unsigned char **)malloc(image.Vpixels * sizeof(unsigned char *));
    for (int i = 0; i < image.Vpixels; i++) {
        image.R[i] = (unsigned char *)malloc(image.Hpixels);
        image.G[i] = (unsigned char *)malloc(image.Hpixels);
        image.B[i] = (unsigned char *)malloc(image.Hpixels);
    }

    FILE *r_fp = fopen(r_file, "r");
    FILE *g_fp = fopen(g_file, "r");
    FILE *b_fp = fopen(b_file, "r");

    if (!r_fp || !g_fp || !b_fp) {
        fprintf(stderr, "無法開啟一個或多個輸入檔案。\n");
        exit(1);
    }

    for (int i = 0; i < image.Vpixels; i++) {
        for (int j = 0; j < image.Hpixels; j++) {
            fscanf(r_fp, "%hhu", &image.R[i][j]);
            fscanf(g_fp, "%hhu", &image.G[i][j]);
            fscanf(b_fp, "%hhu", &image.B[i][j]);
        }
    }

    fclose(r_fp);
    fclose(g_fp);
    fclose(b_fp);

    FILE *output_fp = fopen(output_bmp, "wb");
    if (!output_fp) {
        fprintf(stderr, "無法建立輸出 BMP 檔案: %s\n", output_bmp);
        exit(1);
    }

    unsigned char bmp_header[54] = {
        0x42, 0x4D,                // Signature "BM"
        0, 0, 0, 0,               // File size
        0, 0, 0, 0,               // Reserved
        54, 0, 0, 0,              // Offset to image data
        40, 0, 0, 0,              // Header size
        0, 0, 0, 0,               // Image width
        0, 0, 0, 0,               // Image height
        1, 0,                     // Planes
        24, 0,                    // Bits per pixel
        0, 0, 0, 0,               // Compression
        0, 0, 0, 0,               // Image size
        0, 0, 0, 0,               // X pixels per meter
        0, 0, 0, 0,               // Y pixels per meter
        0, 0, 0, 0,               // Total colors
        0, 0, 0, 0                // Important colors
    };

    // 設定寬高
    *((int *)&bmp_header[18]) = image.Hpixels;
    *((int *)&bmp_header[22]) = image.Vpixels;

    fwrite(bmp_header, sizeof(unsigned char), 54, output_fp);

    int padding = (4 - (image.Hpixels * 3) % 4) % 4;

    for (int i = 0; i < image.Vpixels; i++) {
        for (int j = 0; j < image.Hpixels; j++) {
            fputc(image.B[i][j], output_fp);
            fputc(image.G[i][j], output_fp);
            fputc(image.R[i][j], output_fp);
        }
        for (int p = 0; p < padding; p++) {
            fputc(0, output_fp);
        }
    }

    fclose(output_fp);
    bmp_free(&image);
    printf("Decoder 方法 0 執行完成，輸出檔案為 %s。\n", output_bmp);
}

int bmp_free(bmp *p_bmp) {
    for (int i = 0; i < p_bmp->Vpixels; i++) {
        free(p_bmp->R[i]);
        free(p_bmp->G[i]);
        free(p_bmp->B[i]);
    }
    free(p_bmp->R);
    free(p_bmp->G);
    free(p_bmp->B);
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr, "用法: decoder 0 <ResKimberly.bmp> <R.txt> <G.txt> <B.txt> <dim.txt>\n");
        return 1;
    }

    int method = atoi(argv[1]);
    if (method == 0) {
        decoder_method_0(argv[2], argv[3], argv[4], argv[5], argv[6]);
    } else {
        fprintf(stderr, "無效的方法: %d\n", method);
        return 1;
    }

    return 0;
}
