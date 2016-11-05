all:
	(cd DUDE-Seq-1 && make)
	(cd DUDE-Seq-2 && make)

clean:
	(cd DUDE-Seq-1 && make clean)
	(cd DUDE-Seq-2 && make clean)
	(cd bin && rm -rf *)

install:
	cp DUDE-Seq-1/DUDE-Seq-1 bin/
	cp DUDE-Seq-2/DUDE-Seq-2 bin/
