RSCRIPT ?= Rscript

.PHONY: all reproduce refresh-data analyze test check clean

all: analyze test

reproduce: refresh-data analyze test

refresh-data:
	$(RSCRIPT) --vanilla scripts/prepare_data.R

analyze: data/blueberry_2013_sensory.csv R/analysis_helpers.R scripts/run_analysis.R
	$(RSCRIPT) --vanilla scripts/run_analysis.R

test: analyze
	$(RSCRIPT) --vanilla tests/test_helpers.R
	$(RSCRIPT) --vanilla tests/test_outputs.R

check:
	$(RSCRIPT) --vanilla scripts/check_environment.R

clean:
	$(RM) results/*.csv results/*.json results/*.txt figures/*.png REPORT.md
