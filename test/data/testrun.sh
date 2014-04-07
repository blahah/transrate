mkdir tmp
cd tmp
bundle exec ../../../bin/transrate \
  --assembly ../assembly.fa \
  --reference ../reference.fa \
  --left ../72pc.l.fq \
  --right ../72pc.r.fq \
  --threads 4
rm -rf tmp
