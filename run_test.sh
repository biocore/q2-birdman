#!/bin/bash

#qiime birdman run --i-table q2_birdman/tests/data/94270_filtered.qza --m-metadata-file q2_birdman/tests/data/11913_filtered.tsv --p-formula "age" --o-output-dir q2_birdman/tests/out --p-threads 32
qiime birdman plot --i-results-artifact q2_birdman/tests/out.qza --p-plot-var "age" --o-output-dir q2_birdman/tests/viz

#qiime birdman run --i-table q2_birdman/tests/data/metag.qza --m-metadata-file q2_birdman/tests/data/metadata.tsv --p-formula "host_age+sex+bmi_score+mind_score+assignment" --o-output-dir q2_birdman/tests/out --p-longitudinal True --p-subject-column "host_subject_id"

#qiime birdman run --i-table q2_birdman/tests/data/metag.qza --m-metadata-file q2_birdman/tests/data/metadata.tsv --p-formula "host_age+sex+bmi_score+mind_score+assignment+C(timepoint, Treatment('0'))" --o-output-dir q2_birdman/tests/out --p-longitudinal True --p-subject-column "host_subject_id"
