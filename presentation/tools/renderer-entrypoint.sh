#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: render-slides DECK.pptx" >&2
    exit 2
fi

source_path="/source/$1"
source_stem=${1%.pptx}
profile_path="/tmp/libreoffice-profile-$$"
mkdir -p "$profile_path"

soffice \
    --headless \
    --nologo \
    --nodefault \
    --nofirststartwizard \
    "-env:UserInstallation=file://$profile_path" \
    --convert-to pdf \
    --outdir /output \
    "$source_path"

pdf_path="/output/$source_stem.pdf"
if [ ! -f "$pdf_path" ]; then
    echo "LibreOffice did not create $pdf_path" >&2
    exit 1
fi

pdftoppm \
    -png \
    -scale-to-x 3840 \
    -scale-to-y 2160 \
    "$pdf_path" \
    /output/slide
