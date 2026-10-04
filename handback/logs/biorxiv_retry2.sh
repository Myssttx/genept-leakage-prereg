#!/bin/bash
# Retry bioRxiv pages, honouring Retry-After; full browser headers. Logs one line per attempt.
cd ~/Documents/genept-leakage
UA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0 Safari/537.36'
declare -a NAMES=(zhong_full zhong_pdf supp)
declare -a URLS=('https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2.full' 'https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2.full.pdf' 'https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.supplementary-material')
declare -a OUTS=(data/papers/zhong2025_biorxiv_v2.full.html data/papers/zhong2025_biorxiv_v2.full.pdf data/papers/genept_biorxiv_v2_supplementary-material.html)
declare -a DONE=(0 0 0)
end=$(( $(date +%s) + 1500 ))
while [ $(date +%s) -lt $end ]; do
  for i in 0 1 2; do
    [ ${DONE[$i]} = 1 ] && continue
    # zhong_pdf only needed if zhong_full not done
    [ $i = 1 ] && [ ${DONE[0]} = 1 ] && continue
    t=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    c=$(curl -sSL -m 120 -A "$UA" -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,application/pdf,*/*;q=0.8' -H 'Accept-Language: en-US,en;q=0.9' -D /tmp/biorxiv_hdr.txt -o "${OUTS[$i]}.tmp" -w "%{http_code}" "${URLS[$i]}")
    ra=$(grep -i '^retry-after:' /tmp/biorxiv_hdr.txt | tail -1 | tr -dc '0-9')
    echo "$t ${NAMES[$i]} HTTP $c size=$(stat -f %z "${OUTS[$i]}.tmp") retry-after=${ra:-none}"
    if [ "$c" = 200 ]; then mv "${OUTS[$i]}.tmp" "${OUTS[$i]}"; DONE[$i]=1; else rm -f "${OUTS[$i]}.tmp"; fi
    sleep $(( ${ra:-30} + 10 ))
  done
  if [ ${DONE[2]} = 1 ] && { [ ${DONE[0]} = 1 ] || [ ${DONE[1]} = 1 ]; }; then echo "ALL FETCHED"; exit 0; fi
  sleep 60
done
echo "GAVE UP after 25 min"
