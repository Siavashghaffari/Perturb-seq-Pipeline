#source /gstore/home/yangy197/miniconda3/bin/activate cumulus
alto cromwell run -s cumulus.gcp.science.roche.com \
  -m cumulus/demultiplexing \
  -i demux_inputs_SAM24456643.json \
  -o demux_inputs_updated.json \
  -b gs://gred-melo-carlos-lab/NGS5888/SAM24456643_1/uploads/ \
  --profile dssc-cumulus 
  