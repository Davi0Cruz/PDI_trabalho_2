
run:
	gcc src/main.c -ltiff -o app
	./app hist 333 
	mv *.dat dat/

test:
	gcc src/main.c -ltiff -o app -lm
	for img in imags/333 imags/338 imags/praia imags/mulher; do \
		convert $$img.tif $$img.png; \
		./app hist $$img; \
		./app high-boost $$img 4.5; \
		./app hist $${img}_high-boost; \
		./app laplacian $$img -1; \
		./app hist $${img}_laplacian; \
	done
	mv imags/*.dat dat/
	for img in imags/*_laplacian.tif imags/*_high-boost.tif; do \
		convert $$img $${img%.*}.png; \
		convert $${img%.*}.png -flop $${img%.*}.png; \
	done

rmpng:
	rm imags/*.png