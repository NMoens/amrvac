
#!/usr/bin/env bash
# blk2gif.sh — Generate an animated GIF from cak1d####.blk files using gnuplot + ffmpeg.
# Usage:
#   ./blk2gif.sh /path/to/data out.gif [FPS] [WIDTH] [HEIGHT] [FRAMES_DIR]
# Defaults:
#   FPS=20, WIDTH=1800, HEIGHT=600, FRAMES_DIR=./frames_out
#
# Behavior:
# - Renders each frame with gnuplot (3 panels: 3 vs 1, 5 vs 1, 6 vs 1).
# - Overlays the first snapshot (cak1d0000.blk) as a dotted reference in every frame.
# - If gnuplot fails for any frame, writes a placeholder PNG and continues.
# - Keeps frames in FRAMES_DIR; logs per-frame gnuplot stderr (unless cleaned).
# - Two-pass palette for high-quality GIF.
# - Cleans up per-frame .gp and .log files at the end unless KEEP_ARTIFACTS=1.

# Robust mode: do not exit on non-zero; handle errors explicitly.
set -u
set -o pipefail

DATA_DIR="${1:-.}"
OUT_GIF="${2:-plot.gif}"
FPS="${3:-20}"
WIDTH="${4:-1800}"
HEIGHT="${5:-600}"
FRAMES_DIR="${6:-./frames_out}"

# Set KEEP_ARTIFACTS=1 to keep .gp and .log files.
KEEP_ARTIFACTS="${KEEP_ARTIFACTS:-0}"

echo "==> Data dir: $DATA_DIR"
echo "==> Output GIF: $OUT_GIF"
echo "==> FPS: $FPS  Size: ${WIDTH}x${HEIGHT}"
echo "==> Frames dir: $FRAMES_DIR"
echo "==> Keep artifacts (.gp/.log): ${KEEP_ARTIFACTS}"

# --- Dependencies ---
for cmd in gnuplot ffmpeg; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: '$cmd' not found. Please install it and re-run." >&2
    exit 1
  fi
done

# --- Inputs ---
mapfile -t files < <(LC_ALL=C find "$DATA_DIR" -maxdepth 1 -type f -name 'Ostar_1d[0-9][0-9][0-9][0-9].blk' | sort -V)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No files matching 'Ostar_1d####.blk' found in '$DATA_DIR'." >&2
  exit 1
fi
echo "==> Found ${#files[@]} input files."

# Reference file (first snapshot)
REF_FILE="${files[0]}"
REF_BASE="$(basename "$REF_FILE")"
echo "==> Using reference snapshot: ${REF_BASE}"

# --- Prepare frames directory ---
mkdir -p "$FRAMES_DIR"
# Uncomment to start clean:
# rm -f "$FRAMES_DIR"/frame_*.png "$FRAMES_DIR"/gnuplot_err_*.log "$FRAMES_DIR"/palette.png "$FRAMES_DIR"/*.gp

# --- Detect data separator (comma vs whitespace) using the reference file ---
if grep -q ',' "$REF_FILE"; then
  DATA_SEP=","
else
  DATA_SEP="whitespace"
fi
echo "==> Data separator detected: '${DATA_SEP}' (from ${REF_BASE})"

# --- Determine a working PNG terminal for gnuplot: try pngcairo -> png ---
GNUPLOT_TERM="pngcairo"
if ! gnuplot >/dev/null 2>&1 <<'GNUTERM'
set terminal pngcairo size 100,50
set output ".__gnuterm_test.png"
plot sin(x)
unset output
GNUTERM
then
  echo "!! gnuplot 'pngcairo' terminal not available; falling back to 'png'."
  GNUPLOT_TERM="png"
  if ! gnuplot >/dev/null 2>&1 <<'GNUTERM'
set terminal png size 100,50
set output ".__gnuterm_test.png"
plot sin(x)
unset output
GNUTERM
  then
    echo "ERROR: Neither 'pngcairo' nor 'png' gnuplot terminals are available." >&2
    exit 1
  fi
fi
rm -f .__gnuterm_test.png
echo "==> Using gnuplot terminal: $GNUPLOT_TERM"

# --- Optional fixed axes (leave empty for auto) ---
X_RANGE=""
Y3_RANGE=""
Y5_RANGE=""
Y6_RANGE=""

# --- Helper: placeholder frame on failure (keeps numbering continuous) ---
write_placeholder() {
  local out_png="$1"
  local title="$2"
  local msg="$3"
  local gp="$FRAMES_DIR/__placeholder.gp"

  cat >"$gp" <<GNUPLOT
set terminal ${GNUPLOT_TERM} size ${WIDTH},${HEIGHT} font ",14"
set output "${out_png}"
unset key
set border lw 1.5
set tics out nomirror
set multiplot layout 2,3 title "${title}"
# Left panel: message
set xrange [0:1]; set yrange [0:1]
set label 1 "${msg}" at 0.5,0.6 center tc rgb "#aa0000"
set label 2 "placeholder frame" at 0.5,0.45 center tc rgb "#555555"
plot NaN notitle
# Middle panel
set xrange [0:1]; set yrange [0:1]
plot NaN notitle
# Right panel
set xrange [0:1]; set yrange [0:1]
plot NaN notitle
unset multiplot
unset output
GNUPLOT

  gnuplot "$gp" >/dev/null 2>&1 || true
}

# --- Render frames (never abort; log errors, continue) ---
i=0
for f in "${files[@]}"; do
  printf -v frame "frame_%04d.png" "$i"
  out="$FRAMES_DIR/$frame"
  base="$(basename "$f")"
  err_log="$FRAMES_DIR/gnuplot_err_${frame}.log"
  gp_file="$FRAMES_DIR/plot_${frame%.png}.gp"

  # Extract frame number (#### from filename)
  frame_num=""
  if [[ "$base" =~ ^Ostar_1d([0-9]{4})\.blk$ ]]; then
    frame_num="${BASH_REMATCH[1]}"
  fi

  echo "   -> Rendering: $base  -> $frame"

  # Build gnuplot script for this frame (with dotted reference overlay)
  cat >"$gp_file" <<GNUPLOT
# Auto-generated gnuplot script for $base -> $frame (with reference overlay: ${REF_BASE})
set terminal ${GNUPLOT_TERM} size ${WIDTH},${HEIGHT} font ",14"
set output "${out}"

# Data parsing
set datafile separator ${DATA_SEP}
set datafile commentschars "#"

# Styling
set grid
set key top left
set border lw 1.5
set tics out nomirror
# Current frame styles
set style line 1 lw 2 lc rgb "#1f77b4"  # blue solid
set style line 2 lw 2 lc rgb "#ff7f0e"  # orange solid
set style line 3 lw 2 lc rgb "#2ca02c"  # green solid
# Reference (first snapshot) dotted styles
set style line 11 lw 2 lc rgb "#1f77b4" dt 2  # blue dotted
set style line 12 lw 2 lc rgb "#ff7f0e" dt 2  # orange dotted
set style line 13 lw 2 lc rgb "#2ca02c" dt 2  # green dotted

# Fixed ranges if supplied
if (strlen("${X_RANGE}") > 0) set xrange ${X_RANGE}

set multiplot layout 2,3 title sprintf("File: %s%s", "${base}", (strlen("${frame_num}") > 0 ? sprintf("  |  frame %s", "${frame_num}") : ""))

# Panel 1: column 3 vs 1
set xlabel "radius"
set ylabel "velocity"
if (strlen("${Y3_RANGE}") > 0) set yrange ${Y3_RANGE}
plot "${f}"       using 1:3 with lines ls 1 notitle, \
     "${REF_FILE}" using 1:3 with lines ls 11 notitle
unset yrange

# Panel 2: column 5 vs 1
set ylabel "Radiation Energy"
if (strlen("${Y5_RANGE}") > 0) set yrange ${Y5_RANGE}
set logscale y
plot "${f}"       using 1:5 with lines ls 2 notitle, \
     "${REF_FILE}" using 1:5 with lines ls 12 notitle
unset logscale y

# Panel 3: column 6 vs 1
set ylabel "Radiation Flux"
if (strlen("${Y6_RANGE}") > 0) set yrange ${Y6_RANGE}
plot "${f}"       using 1:6 with lines ls 3 notitle, \
     "${REF_FILE}" using 1:6 with lines ls 13 notitle
unset yrange

# Panel 4:
set ylabel "Density"
set logscale y
plot "${f}"       using 1:2 with lines ls 1 notitle , \
     "${REF_FILE}" using 1:2 with lines ls 11 notitle
unset logscale y

# Panel 5:
set ylabel "Gamma"
if (strlen("${Y6_RANGE}") > 0) set yrange ${Y6_RANGE}
plot "${f}"       using 1:10 with lines ls 2 notitle
unset yrange

# Panel 5:
set ylabel "Temperature"
if (strlen("${Y6_RANGE}") > 0) set yrange ${Y6_RANGE}
plot "${f}"       using 1:15 with lines ls 3 title 'gas', \
     "${f}"       using 1:16 with lines ls 13 title 'rad'
unset yrange

unset multiplot
unset output
GNUPLOT

  # Run gnuplot and capture errors (do NOT exit on failure)
  gnuplot "$gp_file" > /dev/null 2> "$err_log"
  gp_status=$?

  # If plotting failed or output is missing, write a placeholder and continue
  if [[ $gp_status -ne 0 || ! -s "$out" ]]; then
    echo "      WARN: gnuplot failed for $base (status=$gp_status). -> placeholder"
    tail -n 5 "$err_log" 2>/dev/null | sed 's/^/            gnuplot: /' || true
    write_placeholder "$out" "$base | frame ${frame_num}" "Plot error; see log: $(basename "$err_log")"
  fi

  echo "      OK: ${frame} written"
  ((i++))
  if (( i % 50 == 0 )); then echo "      ...rendered $i frames"; fi
done

echo "==> Rendered $i frames to $FRAMES_DIR."

# --- Verify frames exist ---
num_png=$(ls "$FRAMES_DIR"/frame_*.png 2>/dev/null | wc -l || true)
if [[ "$num_png" -eq 0 ]]; then
  echo "ERROR: No PNG frames found in $FRAMES_DIR. GIF cannot be created." >&2
  exit 1
fi
echo "==> Found $num_png PNG frames."

# --- Generate palette ---
PALETTE="$FRAMES_DIR/palette.png"
echo "==> Generating GIF palette..."
ffmpeg -hide_banner -loglevel info -y -framerate "$FPS" -i "$FRAMES_DIR/frame_%04d.png" \
  -vf "palettegen=stats_mode=single" "$PALETTE"

if [[ ! -s "$PALETTE" ]]; then
  echo "ERROR: Palette generation failed (palette file is missing or empty)." >&2
  exit 1
fi

# --- Create GIF using the palette ---
echo "==> Creating GIF '$OUT_GIF'..."
ffmpeg -hide_banner -loglevel info -y -framerate "$FPS" -i "$FRAMES_DIR/frame_%04d.png" -i "$PALETTE" \
  -filter_complex "paletteuse=dither=bayer:bayer_scale=5" \
  "$OUT_GIF"

if [[ ! -s "$OUT_GIF" ]]; then
  echo "ERROR: GIF was not created or is empty: $OUT_GIF" >&2
  exit 1
fi

echo "✅ Done! Animated GIF created: $OUT_GIF"
echo "Frames kept in: $FRAMES_DIR"

# --- Cleanup artifacts unless KEEP_ARTIFACTS=1 ---
if [[ "${KEEP_ARTIFACTS}" != "1" ]]; then
  echo "==> Cleaning up per-frame artifacts (.gp and .log files)..."
  rm -f "$FRAMES_DIR"/plot_frame_*.gp \
        "$FRAMES_DIR"/gnuplot_err_frame_*.log \
        "$FRAMES_DIR"/__placeholder.gp 2>/dev/null || true
  echo "==> Cleanup complete."
else
  echo "==> Artifacts retained (.gp/.log). Set KEEP_ARTIFACTS=0 to remove next run."
fi
