#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define NUM_COMPONENTS 10

// Waveform generator function
static double generate_waveform(int j, double freq, double time) {
    if (j == 0) // sine wave
        return sin(2 * PI * freq * time);
    else if (j == 1) // sawtooth wave
        return 2 * (freq * time - floor(freq * time + 0.5));
    else if (j == 2) // square wave
        return sin(2 * PI * freq * time) >= 0 ? 1 : -1;
    else if (j == 3) // triangle wave
        return 2 * fabs(2 * (freq * time - floor(freq * time + 0.5))) - 1;
    return 0.0;
}

// Unit step function u(t)
static int unit_step(double t) {
    return t >= 0 ? 1 : 0;
}

// Generate x(t)
static double generate_xt(double time, int wave_type) {
    double a[NUM_COMPONENTS] = {100, 2000, 1000, 500, 250, 100, 2000, 1000, 500, 250};
    double f[NUM_COMPONENTS] = {0, 31.25, 500, 2000, 4000, 44, 220, 440, 1760, 3960};
    double xt = 0.0;

    for (int j = 0; j <= 3; ++j) {
        for (int i = 0; i < NUM_COMPONENTS; ++i) {
            double t_offset = time - 0.1 * i - j;
            xt += a[i] * (unit_step(t_offset) - unit_step(t_offset - 0.1)) * generate_waveform(wave_type, f[i], t_offset);
        }
    }
    return xt;
}

// Write WAV header
void write_wav_header(FILE *file, int sample_rate, int bits_per_sample, int channels, int samples) {
    int byte_rate = sample_rate * channels * (bits_per_sample / 8);
    int block_align = channels * (bits_per_sample / 8);
    int data_size = samples * channels * (bits_per_sample / 8);
    int chunk_size = 36 + data_size;

    fwrite("RIFF", 1, 4, file);
    fwrite(&chunk_size, 4, 1, file);
    fwrite("WAVE", 1, 4, file);
    fwrite("fmt ", 1, 4, file);

    int fmt_chunk_size = 16;
    short audio_format = 1;
    fwrite(&fmt_chunk_size, 4, 1, file);
    fwrite(&audio_format, 2, 1, file);
    fwrite(&channels, 2, 1, file);
    fwrite(&sample_rate, 4, 1, file);
    fwrite(&byte_rate, 4, 1, file);
    fwrite(&block_align, 2, 1, file);
    fwrite(&bits_per_sample, 2, 1, file);

    fwrite("data", 1, 4, file);
    fwrite(&data_size, 4, 1, file);
}

void generate_waveform_file(const char *filename, int sample_rate, int bit_depth, int wave_type, double start_time, double duration) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file %s for writing.\n", filename);
        return;
    }

    int total_samples = (int)(sample_rate * duration);
    int max_amplitude = (1 << (bit_depth - 1)) - 1;

    write_wav_header(file, sample_rate, bit_depth, 1, total_samples);

    for (int i = 0; i < total_samples; ++i) {
        double t = start_time + (double)i / sample_rate;
        double sample = generate_xt(t, wave_type);
        short quantized_sample = (short)(sample * max_amplitude / 1000.0);
        fwrite(&quantized_sample, sizeof(short), 1, file);
    }

    fclose(file);
}

int main() {
    const double segment_duration = 0.1; // 每段0.1秒，共40段(0~39)
    const double total_duration = 4.0;
    int sample_rates[] = {8000, 16000};
    int bit_depth = 16;

    for (int rate_idx = 0; rate_idx < 2; ++rate_idx) {
        int sample_rate = sample_rates[rate_idx];
        int total_segments = (int)(total_duration / segment_duration);

        for (int i = 0; i < total_segments; ++i) {
            // 使用整數序號命名，而不是小數點，避免3.9檔案顯示問題
            // 最後一段 i=39 -> seg_39_8k.wav or seg_39_16k.wav
            char filename[50];
            if (sample_rate == 8000) {
                sprintf(filename, "seg_%02d_8k.wav", i);
            } else {
                sprintf(filename, "seg_%02d_16k.wav", i);
            }

            double start_time = i * segment_duration;
            generate_waveform_file(filename, sample_rate, bit_depth, 0, start_time, segment_duration);
            printf("Generated: %s (start=%.1f)\n", filename, start_time);
        }
    }

    return 0;
}
