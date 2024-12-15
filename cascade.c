#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define WAV_HEADER_SIZE 44

int is_wav_file(FILE *f) {
    char riff[4], wave[4];
    if (fread(riff, 1, 4, f) < 4) return 0;
    fseek(f, 4, SEEK_CUR);
    if (fread(wave, 1, 4, f) < 4) return 0;
    fseek(f, 0, SEEK_SET);
    return (memcmp(riff, "RIFF", 4) == 0 && memcmp(wave, "WAVE", 4) == 0);
}

void concatenate_wav(const char *list_file, const char *output_file) {
    FILE *output = fopen(output_file, "wb");
    if (!output) {
        perror("Failed to open output file");
        return;
    }

    FILE *list = fopen(list_file, "r");
    if (!list) {
        perror("Failed to open list file");
        fclose(output);
        return;
    }

    uint8_t buffer[1024];
    size_t read_size;
    uint32_t total_data_size = 0;
    int first_file = 1;

    uint8_t header[WAV_HEADER_SIZE] = {0};
    fwrite(header, 1, WAV_HEADER_SIZE, output);

    char line[256];
    while (fgets(line, sizeof(line), list)) {
        line[strcspn(line, "\n")] = 0;
        if (line[0] == '\0') continue;

        FILE *input = fopen(line, "rb");
        if (!input) {
            fprintf(stderr, "Warning: Could not open %s\n", line);
            continue;
        }

        if (!is_wav_file(input)) {
            fprintf(stderr, "Warning: File %s is not a valid WAV file.\n", line);
            fclose(input);
            continue;
        }

        if (first_file) {
            if (fread(header, 1, WAV_HEADER_SIZE, input) < WAV_HEADER_SIZE) {
                fprintf(stderr, "Warning: %s is too short to be a valid wav.\n", line);
                fclose(input);
                continue;
            }
            fseek(output, 0, SEEK_SET);
            fwrite(header, 1, WAV_HEADER_SIZE, output);
            fseek(output, 0, SEEK_END);
            first_file = 0;
        } else {
            fseek(input, WAV_HEADER_SIZE, SEEK_SET);
        }

        while ((read_size = fread(buffer, 1, sizeof(buffer), input)) > 0) {
            fwrite(buffer, 1, read_size, output);
            total_data_size += (uint32_t)read_size;
        }

        fclose(input);
    }

    fclose(list);

    if (first_file) {
        fprintf(stderr, "No valid WAV files listed in %s.\n", list_file);
        freopen(output_file, "wb", output);
        fclose(output);
        return;
    }

    fseek(output, 4, SEEK_SET);
    uint32_t chunk_size = 36 + total_data_size;
    fwrite(&chunk_size, 4, 1, output);

    fseek(output, 40, SEEK_SET);
    fwrite(&total_data_size, 4, 1, output);

    fclose(output);
    printf("Output file %s generated successfully.\n", output_file);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <list_file> <output_file>\n", argv[0]);
        return 1;
    }

    concatenate_wav(argv[1], argv[2]);
    return 0;
}
