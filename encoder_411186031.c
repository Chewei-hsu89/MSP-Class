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
int bmp_load_fn(char *filename, bmp *p_bmp);
int bmp_save_ascii_fn(char *fn, bmp *p_bmp);
int bmp_free(bmp *p_bmp);

void encoder_method_0(char *input_bmp, char *r_file, char *g_file, char *b_file, char *dim_file) {
    bmp image;
    bmp_load_fn(input_bmp, &image);

    // 儲存 R, G, B 通道到文件
    FILE *r_fp = fopen(r_file, "w");
    FILE *g_fp = fopen(g_file, "w");
    FILE *b_fp = fopen(b_file, "w");
    FILE *dim_fp = fopen(dim_file, "w");
    fprintf(dim_fp, "%d %d\n", image.Hpixels, image.Vpixels);

    for (int i = 0; i < image.Vpixels; i++) {
        for (int j = 0; j < image.Hpixels; j++) {
            fprintf(r_fp, "%d ", image.R[i][j]);
            fprintf(g_fp, "%d ", image.G[i][j]);
            fprintf(b_fp, "%d ", image.B[i][j]);
        }
        fprintf(r_fp, "\n");
        fprintf(g_fp, "\n");
        fprintf(b_fp, "\n");
    }

    fclose(r_fp);
    fclose(g_fp);
    fclose(b_fp);
    fclose(dim_fp);

    bmp_free(&image);
    printf("Method 0 completed.\n");
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: encoder <method> <args...>\n");
        return 1;
    }

    int method = atoi(argv[1]);
    if (method == 0) {
        if (argc != 7) {
            fprintf(stderr, "Usage: encoder 0 <input.bmp> <R.txt> <G.txt> <B.txt> <dim.txt>\n");
            return 1;
        }
        encoder_method_0(argv[2], argv[3], argv[4], argv[5], argv[6]);
    } else {
        fprintf(stderr, "Invalid method: %d\n", method);
        return 1;
    }

    return 0;
}

// BMP 加載函數
int bmp_load_fn(char *filename, bmp *p_bmp) {
    FILE *f = fopen(filename, "rb");
    if (!f) {
        fprintf(stderr, "File %s not found.\n", filename);
        return 0;
    }

    fread(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f);
    p_bmp->Hpixels = *((int *)&p_bmp->HeaderInfo[18]);
    p_bmp->Vpixels = *((int *)&p_bmp->HeaderInfo[22]);

    int RowBytes = (p_bmp->Hpixels * 3 + 3) & (~3);
    p_bmp->Hbytes = RowBytes;

    p_bmp->data = (unsigned char **)malloc(p_bmp->Vpixels * sizeof(unsigned char *));
    for (int i = 0; i < p_bmp->Vpixels; i++) {
        p_bmp->data[i] = (unsigned char *)malloc(RowBytes);
        fread(p_bmp->data[i], sizeof(unsigned char), RowBytes, f);
    }

    p_bmp->R = (unsigned char **)malloc(p_bmp->Vpixels * sizeof(unsigned char *));
    p_bmp->G = (unsigned char **)malloc(p_bmp->Vpixels * sizeof(unsigned char *));
    p_bmp->B = (unsigned char **)malloc(p_bmp->Vpixels * sizeof(unsigned char *));
    for (int i = 0; i < p_bmp->Vpixels; i++) {
        p_bmp->R[i] = (unsigned char *)malloc(p_bmp->Hpixels);
        p_bmp->G[i] = (unsigned char *)malloc(p_bmp->Hpixels);
        p_bmp->B[i] = (unsigned char *)malloc(p_bmp->Hpixels);
        for (int j = 0; j < p_bmp->Hpixels; j++) {
            p_bmp->B[i][j] = p_bmp->data[i][j * 3];
            p_bmp->G[i][j] = p_bmp->data[i][j * 3 + 1];
            p_bmp->R[i][j] = p_bmp->data[i][j * 3 + 2];
        }
    }

    fclose(f);
    return 1;
}

int bmp_free(bmp *p_bmp) {
    for (int i = 0; i < p_bmp->Vpixels; i++) {
        free(p_bmp->data[i]);
        free(p_bmp->R[i]);
        free(p_bmp->G[i]);
        free(p_bmp->B[i]);
    }
    free(p_bmp->data);
    free(p_bmp->R);
    free(p_bmp->G);
    free(p_bmp->B);

    return 1;
}


