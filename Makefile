build:
	+make -C src/ build ; \
	mkdir -p bin/ ; \
	cp src/complete-contigs src/maximality bin/	

clean:
	+make -C src/ clean ; \
	rm -rf bin/
