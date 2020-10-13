#!/bin/bash
tempdir=$(mktemp --directory)
trap "rm -r $tempdir" EXIT SIGTERM

for infile in "$@"; do
  sort "$infile" > "${tempdir}/${infile}.sorted"
  if [ -e "${tempdir}/final.results" ]
  then
    join -a1 -a2 -e "NULL" -o auto \
      "${tempdir}/final.results" "${tempdir}/${infile}.sorted" \
      > "${tempdir}/res"
    mv "${tempdir}/res" "${tempdir}/final.results"
  else
    cp "${tempdir}/${infile}.sorted" "${tempdir}/final.results"
  fi
done
cat "${tempdir}/final.results"