.PHONY : casestudy

rawdata/IE/IE-convert-15050.csv : convert_IE_data.py IEDATA
	@mkdir -p $(@D) && \
	python $< M=50 P=50 exinh=M3=1 exinh=M5=1 exinh=M12=1 exinh=M14=1exinh=M15=1 exinh=P3=-1 exc=M11

outputs/AA/MaxSyn/%.txt : fit_data.py rawdata/AA/%.csv
	@mkdir -p $(@D) && \
	python $< $(word 2,$^)

outputs/AA/comp/%.txt : run_paup.py rawdata/AA/%.csv
	@mkdir -p outputs/AA/comp && \
	mkdir -p outputs/AA/wmp && \
	mkdir -p nexfiles/AA && \
	python $< $(word 2,$^)

outputs/IE/MaxSyn/%.txt : fit_data.py rawdata/IE/%.csv
	@mkdir -p $(@D) && \
	python $< $(word 2,$^)

outputs/IE/comp/%.txt : run_paup.py rawdata/IE/%.csv
	@mkdir -p outputs/IE/comp && \
	mkdir -p outputs/IE/wmp && \
	mkdir -p nexfiles/IE && \
	python $< $(word 2,$^)

casestudy : outputs/AA/MaxSyn/AA-Attested-15050.txt outputs/AA/comp/AA-Attested-15050.txt outputs/AA/MaxSyn/AA-Proto-15050.txt outputs/AA/comp/AA-Proto-15050.txt outputs/IE/MaxSyn/IE-convert-15050.txt outputs/IE/comp/IE-convert-15050.txt
