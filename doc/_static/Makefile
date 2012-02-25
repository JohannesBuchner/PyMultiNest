
# you need to run pyapemost_demo.py to obtain the three pdf files.

all: marginals.png chain0.png data.png
clean: 
	rm -f marginals.png chain0.png data.png

%.png: %.pdf
	convert -trim $< $@
.PHONY: all clean
