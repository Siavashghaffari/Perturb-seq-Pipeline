#source /gstore/home/yangy197/miniconda3/bin/activate cumulus
alto cromwell run -s cumulus.gcp.science.roche.com \
  -m cumulus/cellranger:4.0 \
  -i inputs_SAM24456643.json \
  -o inputs_updated.json \
  -b gs://gred-melo-carlos-lab/NGS5888/SAM24456643_1/uploads/ \
  --no-ssl-verify \
  --no-cache \
  --profile dssc-cumulus | tee NGS5888.txt
  