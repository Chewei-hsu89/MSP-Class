#!/bin/bash

# 編譯 C 程式
gcc -o sinegen.exe sinegen.c -lm
gcc -o cascade.exe cascade.c -lm
gcc -o spectrogram.exe spectrogram.c -lm

# 1. 使用 sinegen.exe 產生 0~39 段的 8kHz 與16kHz檔案 (seg_00_8k.wav ... seg_39_8k.wav)
./sinegen.exe

# 建立 list-8k.txt 與 list-16k.txt
> list-8k.txt
for ((i=0;i<40;i++)); do
    printf "seg_%02d_8k.wav\n" $i >> list-8k.txt
done

> list-16k.txt
for ((i=0;i<40;i++)); do
    printf "seg_%02d_16k.wav\n" $i >> list-16k.txt
done

# 串接成 s-8k.wav 與 s-16k.wav
./cascade.exe list-8k.txt s-8k.wav
./cascade.exe list-16k.txt s-16k.wav

# 假設 aeueo-8kHz.wav、aeueo-16kHz.wav 已存在
wav_files=("s-8k.wav" "s-16k.wav" "aeueo-8kHz.wav" "aeueo-16kHz.wav")

# 生成 ASCII 頻譜檔案
for wav_in in "${wav_files[@]}"; do
    base_name=$(basename "$wav_in" .wav)
    ./spectrogram.exe 32 rectangular 32 10 "$wav_in" "${base_name}.Set1.txt"
    ./spectrogram.exe 32 hamming    32 10 "$wav_in" "${base_name}.Set2.txt"
    ./spectrogram.exe 30 rectangular 32 10 "$wav_in" "${base_name}.Set3.txt"
    ./spectrogram.exe 30 hamming    32 10 "$wav_in" "${base_name}.Set4.txt"
done

# 生成 PDF 圖形
for wav_in in "${wav_files[@]}"; do
    base_name=$(basename "$wav_in" .wav)
    for set_num in 1 2 3 4; do
        in_txt="${base_name}.Set${set_num}.txt"
        out_pdf="${base_name}.Set${set_num}.pdf"
        python3 spectshow.py "$wav_in" "$in_txt" "$out_pdf"
    done
done

echo "All processes completed successfully."
