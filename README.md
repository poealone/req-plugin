# ReQ — Parametric EQ for PocketDAW

4-band stereo parametric EQ plugin for [PocketDAW](https://pocketdaw.net). Drag the nodes directly on the curve to set frequency and gain. Built with the PocketDAW SDK v3.1.

![ReQ Screenshot](https://pocketdaw.net/assets/screenshots/pdaw-param-eq.png)

## Bands

| Band | Type | Frequency Range | Gain | Q |
|------|------|----------------|------|---|
| 1 | Low Shelf | 20–500 Hz | ±12 dB | fixed |
| 2 | Peak | 80–2000 Hz | ±12 dB | 0.3–4.0 |
| 3 | Peak | 800–12000 Hz | ±12 dB | 0.3–4.0 |
| 4 | High Shelf | 2000–20000 Hz | ±12 dB | fixed |

Plus a **Master Gain** control (±6 dB).

## Presets

- **Flat** — all bands at 0 dB
- **Low Cut** — attenuate lows, presence boost
- **Vocal Presence** — dip 300 Hz, boost 3 kHz
- **Bass Boost** — low shelf up
- **Air** — high shelf boost

## Custom UI

The full-screen editor renders a **live frequency response curve** with draggable nodes per band. On desktop (v3.1), click and drag any node to sweep frequency (left/right) and adjust gain (up/down).

Band colors:
- 🟠 Band 1 (Low Shelf)
- 🔵 Band 2 (Peak)
- 🟣 Band 3 (Peak)
- 🔴 Band 4 (High Shelf)

## Build

```bash
# Linux x64
gcc -shared -fPIC -O2 -o req.so req.c -lm

# Anbernic (aarch64)
aarch64-none-linux-gnu-gcc -shared -fPIC -O2 -o req.so req.c -lm

# Windows x64
x86_64-w64-mingw32-gcc -shared -o req.dll req.c -lm
```

## Install

```bash
# Desktop
cp req.so path/to/PocketDAW/plugins/fx/

# Anbernic
scp req.so root@<device-ip>:/mnt/mmc/MUOS/application/PocketDAW/plugins/fx/
```

ReQ appears in the FX chain under `filter`. Works in both the mixer strip FX chain and per-track insert FX slots.

## License

Proprietary — free to use. Built for PocketDAW.
