parse_lat() {
    local s="$1"
    local deg="${s:0:2}"
    local hemi="${s:2:1}"
    if [[ "$hemi" == "S" ]]; then
        echo $((-10#$deg))
    else
        echo $((10#$deg))
    fi
}

parse_lon() {
    local s="$1"
    local deg="${s:0:3}"
    local hemi="${s:3:1}"
    if [[ "$hemi" == "W" ]]; then
        echo $((-10#$deg))
    else
        echo $((10#$deg))
    fi
}

fmt_lat() {
    local x=$1
    if (( x < 0 )); then
        printf '%02dS' $((-x))
    else
        printf '%02dN' "$x"
    fi
}

fmt_lon() {
    local x=$1
    if (( x < 0 )); then
        printf '%03dW' $((-x))
    else
        printf '%03dE' "$x"
    fi
}

mkdir -p split

for infile in *.tif; do
    base="${infile%.tif}"
    lat_tag="${base%_*}"
    lon_tag="${base#*_}"

    lat=$(parse_lat "$lat_tag")
    lon=$(parse_lon "$lon_tag")

    lat_s=$((lat - 5))
    lon_e=$((lon + 5))

    gdal_translate -srcwin 0 0 20000 20000 \
		   -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=2 \
      "$infile" "split/$(fmt_lat "$lat")_$(fmt_lon "$lon").tif"

    gdal_translate -srcwin 20000 0 20000 20000 \
		   -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=2 \
      "$infile" "split/$(fmt_lat "$lat")_$(fmt_lon "$lon_e").tif"

    gdal_translate -srcwin 0 20000 20000 20000 \
		   -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=2 \
      "$infile" "split/$(fmt_lat "$lat_s")_$(fmt_lon "$lon").tif"

    gdal_translate -srcwin 20000 20000 20000 20000 \
		   -co COMPRESS=LZW -co TILED=YES -co PREDICTOR=2 \
      "$infile" "split/$(fmt_lat "$lat_s")_$(fmt_lon "$lon_e").tif"
done
