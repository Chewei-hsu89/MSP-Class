import sys
import matplotlib.pyplot as plt
import numpy as np
import wave

if len(sys.argv) != 4:
    print("Usage: python3 spectshow.py in_wav in_txt out_pdf")
    sys.exit(1)

in_wav, in_txt, out_pdf = sys.argv[1], sys.argv[2], sys.argv[3]

with wave.open(in_wav, 'rb') as w:
    frames = w.readframes(-1)
    samples = np.frombuffer(frames, dtype=np.int16)
    time = np.linspace(0, len(samples) / w.getframerate(), len(samples))

spec = np.loadtxt(in_txt)

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(time, samples, color='black')
plt.title("Waveform")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")

plt.subplot(2, 1, 2)
plt.imshow(spec.T, origin='lower', aspect='auto', cmap='jet')
plt.title("Spectrogram")
plt.xlabel("Time")
plt.ylabel("Frequency")
plt.colorbar

