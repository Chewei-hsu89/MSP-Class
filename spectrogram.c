#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define WAV_HEADER_SIZE 44
#define PI 3.14159265358979323846

int read_wav_header(FILE *f, int *sample_rate, int *total_samples, int *bits_per_sample) {
    uint8_t header[WAV_HEADER_SIZE];
    if (fread(header, 1, WAV_HEADER_SIZE, f) < WAV_HEADER_SIZE) return 0;
    if (memcmp(header, "RIFF", 4) != 0 || memcmp(header+8, "WAVE", 4) != 0) {
        fprintf(stderr, "Invalid WAV file or unsupported format.\n");
        return 0;
    }

    memcpy(sample_rate, header+24, 4);
    memcpy(bits_per_sample, header+34, 2);
    uint32_t data_size;
    memcpy(&data_size, header+40, 4);

    if (*bits_per_sample != 16) {
        fprintf(stderr, "Warning: This program is designed for 16-bit PCM. Detected %d-bit. Result may be incorrect.\n", *bits_per_sample);
    }

    int bytes_per_sample = *bits_per_sample / 8;
    if (bytes_per_sample == 0) {
        fprintf(stderr, "Invalid bits_per_sample: %d\n", *bits_per_sample);
        return 0;
    }

    *total_samples = data_size / bytes_per_sample;
    return 1;
}

void dft(const double *in, double *out_re, double *out_im, int N) {
    for (int k=0; k<N; k++) {
        double re = 0.0;
        double im = 0.0;
        for (int n=0; n<N; n++) {
            double theta = -2.0 * PI * k * n / N;
            re += in[n]*cos(theta);
            im += in[n]*sin(theta);
        }
        out_re[k] = re;
        out_im[k] = im;
    }
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s w_size w_type dft_size f_itv wav_in spec_out\n", argv[0]);
        return 1;
    }

    int w_size_ms = atoi(argv[1]);
    char *w_type_str = argv[2];
    int dft_size_ms = atoi(argv[3]);
    int f_itv_ms = atoi(argv[4]);
    char *wav_in = argv[5];
    char *spec_out = argv[6];

    FILE *fin = fopen(wav_in, "rb");
    if (!fin) {
        perror("Failed to open input wav");
        return 1;
    }

    int sample_rate, total_samples, bits_per_sample;
    if (!read_wav_header(fin, &sample_rate, &total_samples, &bits_per_sample)) {
        fclose(fin);
        return 1;
    }

    int bytes_per_sample = bits_per_sample / 8;
    if (bytes_per_sample == 0) {
        fprintf(stderr, "Invalid bits_per_sample.\n");
        fclose(fin);
        return 1;
    }

    int w_size_samps = w_size_ms * sample_rate / 1000;
    int dft_size_samps = dft_size_ms * sample_rate / 1000;
    int f_itv_samps = f_itv_ms * sample_rate / 1000;

    uint8_t *raw_data = malloc(total_samples * bytes_per_sample);
    if (!raw_data) {
        fprintf(stderr, "Memory allocation error.\n");
        fclose(fin);
        return 1;
    }

    if (fread(raw_data, bytes_per_sample, total_samples, fin) < (size_t)total_samples) {
        fprintf(stderr, "Error reading wav data.\n");
        fclose(fin);
        free(raw_data);
        return 1;
    }
    fclose(fin);

    // 僅正確處理16-bit PCM，其他格式雖不退出，但結果可能不正確
    double *pcm = malloc(sizeof(double)*total_samples);
    if (!pcm) {
        fprintf(stderr, "Memory allocation error.\n");
        free(raw_data);
        return 1;
    }

    if (bits_per_sample == 16) {
        int16_t *pcm_data = (int16_t*)raw_data;
        for (int i=0; i<total_samples; i++) {
            pcm[i] = pcm_data[i];
        }
    } else {
        // 非16-bit暫不精確處理，只是示範
        uint32_t max_val = (1U << bits_per_sample) - 1;
        for (int i=0; i<total_samples; i++) {
            uint32_t val = 0;
            for (int b=0; b<bytes_per_sample; b++) {
                val |= (uint32_t)raw_data[i*bytes_per_sample+b] << (8*b);
            }
            double normalized = (double)val / (double)max_val;
            pcm[i] = (normalized - 0.5)*32767.0;
        }
    }

    free(raw_data);

    double *window = malloc(sizeof(double)*w_size_samps);
    if (!window) {
        fprintf(stderr, "Memory allocation error.\n");
        free(pcm);
        return 1;
    }

    if (strcmp(w_type_str, "hamming") == 0) {
        for (int i=0; i<w_size_samps; i++) {
            window[i] = 0.54 - 0.46*cos(2.0*PI*i/(w_size_samps-1));
        }
    } else {
        for (int i=0; i<w_size_samps; i++) {
            window[i] = 1.0;
        }
    }

    FILE *fout = fopen(spec_out, "w");
    if (!fout) {
        perror("Failed to open spec_out");
        free(pcm);
        free(window);
        return 1;
    }

    double *frame = malloc(sizeof(double)*dft_size_samps);
    double *out_re = malloc(sizeof(double)*dft_size_samps);
    double *out_im = malloc(sizeof(double)*dft_size_samps);

    if (!frame || !out_re || !out_im) {
        fprintf(stderr, "Memory allocation error.\n");
        fclose(fout);
        free(pcm);
        free(window);
        if (frame) free(frame);
        if (out_re) free(out_re);
        if (out_im) free(out_im);
        return 1;
    }

    for (int start=0; start+w_size_samps<=total_samples; start+=f_itv_samps) {
        for (int i=0; i<dft_size_samps; i++) {
            if (i < w_size_samps) {
                frame[i] = pcm[start+i]*window[i];
            } else {
                frame[i] = 0.0;
            }
        }

        dft(frame, out_re, out_im, dft_size_samps);

        for (int k=0; k<dft_size_samps; k++) {
            double mag = sqrt(out_re[k]*out_re[k] + out_im[k]*out_im[k]);
            double mag_dB = 20.0*log10(mag+1e-12);
            fprintf(fout, "%f ", mag_dB);
        }
        fprintf(fout, "\n");
    }

    fclose(fout);
    free(pcm);
    free(window);
    free(frame);
    free(out_re);
    free(out_im);

    return 0;
}
