projDir="../"
nextflow run \
  ${projDir} \
  -with-dag gbrs_nxf.png \
  -profile sumner \
  -ansi-log true \
  -params-file gbrs_test01.yaml \
  $*

  #-dump-channels true \
  #-dump-hashes true \

